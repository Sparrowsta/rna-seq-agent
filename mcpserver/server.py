from fastapi import FastAPI, HTTPException, BackgroundTasks, File, UploadFile, Request
from fastapi.staticfiles import StaticFiles
from sse_starlette.sse import EventSourceResponse
import os
import pandas as pd
import io
import json
import asyncio
import sys
from pydantic import BaseModel, Field
from typing import Dict, Any, List, Optional
import uuid
import datetime

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
import launch

# --- Pydantic Models for API Data Structure ---

class PipelineRequest(BaseModel):
    srr_list: Optional[str] = Field(None, description="Path to the SRR list file.")
    fastq_dir: Optional[str] = Field(None, description="Path to the directory with local fastq files.")
    seq_mode: bool = Field(False, description="True for single-end, False for paired-end.")
    species: str = Field(..., description="Species to process (e.g., 'human', 'mouse').")
    genome_version: str = Field(..., description="Genome version (e.g., 'hg38').")
    seq_type: str = Field("rna-seq", description="Type of sequencing analysis.")

class Task(BaseModel):
    task_id: str
    status: str
    params: PipelineRequest
    created_at: datetime.datetime
    started_at: Optional[datetime.datetime] = None
    completed_at: Optional[datetime.datetime] = None
    results: Dict[str, Any] = {}

class TaskInfo(BaseModel):
    task_id: str
    status: str
    message: str

class UpdateResultRequest(BaseModel):
    step_name: str = Field(..., description="The name of the workflow step (e.g., 'fastp', 'star').")
    data: Dict[str, Any] = Field(..., description="The result data from the step.")

class GroupInfo(BaseModel):
    control_group: List[str] = Field(..., description="List of sample IDs in the control group.")
    experiment_group: List[str] = Field(..., description="List of sample IDs in the experiment group.")

# --- FastAPI App Initialization ---

app = FastAPI(
    title="MCP Server",
    description="A server to manage and orchestrate bioinformatics workflows with LLM integration.",
    version="0.2.0",
)

# Mount the '/data' directory to serve all result files statically
# This allows the client to access any file in the data directory via the /results URL path
# e.g., /results/de_analysis/{task_id}/volcano_plot.png
data_dir = "/data"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)
app.mount("/results", StaticFiles(directory=data_dir), name="results")


tasks_db: Dict[str, Task] = {}
# Global dictionary to hold event queues for each task
task_event_queues: Dict[str, asyncio.Queue] = {}


# --- SSE and Background Task Logic ---

async def send_event(queue: asyncio.Queue, event_name: str, data: dict):
    """Helper to put a formatted event into the queue."""
    await queue.put({"event": event_name, "data": json.dumps(data)})

async def stream_subprocess(command: List[str], task_id: str, queue: asyncio.Queue):
    """
    Runs a subprocess and streams its stdout/stderr to an asyncio queue.
    """
    process = await asyncio.create_subprocess_exec(
        *command,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd='/data'
    )

    # Concurrently read stdout and stderr
    async def stream_to_queue(stream, event_name):
        while True:
            line = await stream.readline()
            if not line:
                break
            await send_event(queue, event_name, {"line": line.decode().strip()})

    await asyncio.gather(
        stream_to_queue(process.stdout, "log"),
        stream_to_queue(process.stderr, "error_log")
    )

    return_code = await process.wait()
    return return_code

async def run_pipeline_in_background(task_id: str, params: PipelineRequest):
    """Wrapper function to run the pipeline and push status updates via SSE."""
    task = tasks_db[task_id]
    queue = task_event_queues[task_id]
    
    try:
        task.status = "running"
        task.started_at = datetime.datetime.now()
        await send_event(queue, "status_update", {"status": "running", "message": f"Task {task_id} started."})

        # Convert Pydantic model to a dictionary for the launch script
        launch_params = params.dict()
        launch_params['task_id'] = task_id
        
        # This function in launch.py now needs to return the command list
        command = launch.prepare_pipeline_command(launch_params)
        
        return_code = await stream_subprocess(command, task_id, queue)

        if return_code == 0:
            task.status = "completed"
            await send_event(queue, "status_update", {"status": "completed", "message": f"Task {task_id} completed successfully."})
        else:
            raise Exception(f"Pipeline process exited with non-zero code: {return_code}")

    except Exception as e:
        task.status = "failed"
        task.results["error"] = str(e)
        await send_event(queue, "status_update", {"status": "failed", "message": f"Task {task_id} failed: {e}"})
    finally:
        task.completed_at = datetime.datetime.now()
        await send_event(queue, "task_end", {"status": task.status})
        # Signal the SSE stream to close
        await queue.put(None)

# --- API Endpoints ---

@app.get("/", tags=["General"])
async def read_root():
    return {"message": "Welcome to the MCP Server!"}


@app.post("/pipeline/run", response_model=TaskInfo, status_code=202, tags=["Pipeline"])
async def run_pipeline(params: PipelineRequest, background_tasks: BackgroundTasks):
    """
    Creates and runs a new analysis pipeline as a background task.
    A queue is created for this task to stream live updates via SSE.
    """
    task_id = str(uuid.uuid4())
    new_task = Task(
        task_id=task_id,
        status="pending",
        params=params,
        created_at=datetime.datetime.now(),
    )
    tasks_db[task_id] = new_task
    task_event_queues[task_id] = asyncio.Queue()
    
    # Add the long-running job to the background
    background_tasks.add_task(run_pipeline_in_background, task_id, params)
    
    return TaskInfo(
        task_id=task_id,
        status="pending",
        message="Pipeline task has been accepted. Connect to /tasks/{task_id}/stream for live updates."
    )

@app.get("/tasks/{task_id}", response_model=Task, tags=["Tasks"])
async def get_task_status(task_id: str):
    task = tasks_db.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")
    return task

@app.post("/tasks/{task_id}/results", response_model=Task, tags=["Tasks"])
async def update_task_results(task_id: str, result_update: UpdateResultRequest):
    """
    Endpoint for Nextflow to post back results of individual steps.
    This now also pushes an update to the SSE queue.
    """
    task = tasks_db.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    new_status = f"completed_{result_update.step_name}"
    task.results[result_update.step_name] = result_update.data
    task.status = new_status
    
    if queue := task_event_queues.get(task_id):
        await send_event(queue, "status_update", {"status": new_status, "message": f"Step '{result_update.step_name}' completed."})

    print(f"Task {task_id} updated with results from step: {result_update.step_name}")
    return task

@app.post("/tasks/{task_id}/meta", response_model=Task, tags=["Tasks"])
async def upload_meta_file(task_id: str, file: UploadFile = File(...)):
    """
    Uploads and processes a metadata file (CSV or TSV) for a specific task.
    """
    task = tasks_db.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    try:
        contents = await file.read()
        # Use pandas to read the file content from memory
        # This handles both CSV and TSV by sniffing the separator
        df = pd.read_csv(io.StringIO(contents.decode('utf-8')), sep=None, engine='python')
        
        # Convert dataframe to a list of records (JSON)
        meta_data = df.to_dict(orient='records')
        
        task.results["meta_data"] = meta_data
        task.status = "meta_uploaded"
        print(f"Metadata for task {task_id} uploaded and processed.")
        
        return task
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to process meta file: {e}")

@app.post("/tasks/{task_id}/groups", response_model=Task, tags=["Tasks"])
async def define_groups(task_id: str, groups: GroupInfo):
    """
    Defines the control and experiment groups for differential analysis.
    """
    task = tasks_db.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    task.results["groups"] = groups.dict()
    task.status = "groups_defined"
    print(f"Groups for task {task_id} have been defined.")
    
    return task

async def run_de_analysis_in_background(task_id: str):
    """
    Resumes a Nextflow pipeline to run the downstream DE analysis, with SSE.
    """
    task = tasks_db[task_id]
    queue = task_event_queues[task_id]

    try:
        task.status = "running_de_analysis"
        await send_event(queue, "status_update", {"status": "running_de_analysis", "message": "Starting DE analysis."})

        # --- 1. Prepare parameters for the new Nextflow run ---
        resume_params = task.params.dict()
        resume_params['run_de_analysis'] = True
        
        meta_data = task.results.get("meta_data")
        if not meta_data:
            raise ValueError("Metadata not found.")
        
        meta_dir = f"/data/meta/{task_id}"
        os.makedirs(meta_dir, exist_ok=True)
        meta_file_path = os.path.join(meta_dir, "meta.csv")
        pd.DataFrame(meta_data).to_csv(meta_file_path, index=False)
        resume_params['meta_file'] = meta_file_path

        groups = task.results.get("groups")
        if not groups:
            raise ValueError("Group definitions not found.")
        resume_params['control_group'] = ','.join(groups['control_group'])
        resume_params['experiment_group'] = ','.join(groups['experiment_group'])

        # --- 2. Construct the Nextflow command ---
        nextflow_params = []
        for key, value in resume_params.items():
            if isinstance(value, bool):
                if value: nextflow_params.append(f"--{key}")
            elif value is not None:
                nextflow_params.extend([f"--{key}", str(value)])

        command = [
            "nextflow", "run", "/app/nextflow/main.nf",
            "-profile", "inside_container,standard",
            "-resume",
            "-work-dir", "/data/work",
            "--task_id", task_id
        ] + nextflow_params

        # --- 3. Execute the command ---
        return_code = await stream_subprocess(command, task_id, queue)
        if return_code != 0:
            raise Exception(f"DE analysis process exited with non-zero code: {return_code}")

        task.status = "completed_de_analysis"
        await send_event(queue, "status_update", {"status": "completed_de_analysis", "message": "DE analysis finished successfully."})

    except Exception as e:
        task.status = "failed_de_analysis"
        task.results["error"] = str(e)
        await send_event(queue, "status_update", {"status": "failed_de_analysis", "message": f"DE analysis failed: {e}"})
    finally:
        await send_event(queue, "task_end", {"status": task.status})
        await queue.put(None) # Signal end

@app.post("/tasks/{task_id}/analyze", response_model=TaskInfo, status_code=202, tags=["Analysis"])
async def start_de_analysis(task_id: str, background_tasks: BackgroundTasks):
    """
    Triggers the downstream differential expression analysis for a task
    by resuming the original Nextflow pipeline.
    """
    task = tasks_db.get(task_id)
    if not task:
        raise HTTPException(status_code=404, detail="Task not found")

    if "groups" not in task.results or "meta_data" not in task.results:
        raise HTTPException(status_code=400, detail="Metadata and group definitions must be provided before starting analysis.")

    background_tasks.add_task(run_de_analysis_in_background, task_id)

    return TaskInfo(
        task_id=task_id,
        status="de_analysis_started",
        message="DE analysis started. Connect to the stream for live updates."
    )

@app.get("/tasks/{task_id}/stream", tags=["Tasks"])
async def stream_task_events(request: Request, task_id: str):
    """
    SSE endpoint to stream live updates for a given task.
    Handles invalid task_id gracefully by sending an SSE error event.
    """
    queue = task_event_queues.get(task_id)

    async def event_generator():
        # If the queue doesn't exist for the task, send an error and close the stream.
        if not queue:
            yield {
                "event": "error",
                "data": json.dumps({"message": f"Task ID '{task_id}' not found or has already completed."})
            }
            return

        # If the queue exists, proceed with streaming events.
        try:
            while True:
                if await request.is_disconnected():
                    print(f"Client disconnected from task {task_id} stream.")
                    break
                
                event = await queue.get()
                if event is None: # A None event is the signal to close the connection.
                    break
                
                yield event
                
        except asyncio.CancelledError:
            print(f"SSE stream for task {task_id} was cancelled by the client.")
        finally:
            print(f"Closing SSE stream for task {task_id}.")
            # Clean up the queue when the stream is closed.
            # This part is important to prevent memory leaks.
            if task_id in task_event_queues:
                # Ensure the queue is empty.
                while not task_event_queues[task_id].empty():
                    task_event_queues[task_id].get_nowait()
                # Remove the queue from the global dictionary.
                del task_event_queues[task_id]
                print(f"Cleaned up and removed queue for task {task_id}.")


    return EventSourceResponse(event_generator())


# --- Main Execution ---

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=48001)