import os
import requests
from langchain_openai import ChatOpenAI
from langchain.agents import AgentExecutor, create_react_agent
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.tools import tool
from pydantic import BaseModel, Field
from typing import Optional, List
from dotenv import load_dotenv

# --- Environment and API Configuration ---
load_dotenv()
MCPSERVER_URL = os.environ.get("MCPSERVER_URL", "http://localhost:48001")

# --- Pydantic Models for Tool Inputs ---
class RNASeqPipelineInput(BaseModel):
    """Input model for the run_rnaseq_pipeline tool."""
    srr_list: Optional[str] = Field(None, description="Path to the SRR list file.")
    fastq_dir: Optional[str] = Field(None, description="Path to the directory with local fastq files.")
    seq_mode: bool = Field(False, description="True for single-end, False for paired-end.")
    species: str = Field("human", description="Species to process (e.g., 'human', 'mouse').")
    genome_version: str = Field("hg38", description="Genome version (e.g., 'hg38').")
    seq_type: str = Field("rna-seq", description="Type of sequencing analysis.")

class TaskStatusInput(BaseModel):
    """Input model for the get_task_status tool."""
    task_id: str = Field(..., description="The ID of the task to check.")

class UploadMetaInput(BaseModel):
    """Input model for the upload_meta_file tool."""
    task_id: str = Field(..., description="The ID of the task.")
    file_path: str = Field(..., description="The local path to the metadata file.")

class DefineGroupsInput(BaseModel):
    """Input model for the define_groups tool."""
    task_id: str = Field(..., description="The ID of the task.")
    control_group: List[str] = Field(..., description="List of sample IDs in the control group.")
    experiment_group: List[str] = Field(..., description="List of sample IDs in the experiment group.")

class StartDEAnalysisInput(BaseModel):
    """Input model for the start_de_analysis tool."""
    task_id: str = Field(..., description="The ID of the task.")

class GetAnalysisResultsInput(BaseModel):
    """Input model for the get_analysis_results tool."""
    task_id: str = Field(..., description="The ID of the task.")

# --- Agent Tools ---
@tool(args_schema=RNASeqPipelineInput)
def run_rnaseq_pipeline(srr_list: Optional[str] = None, fastq_dir: Optional[str] = None, seq_mode: bool = False, species: str = "human", genome_version: str = "hg38", seq_type: str = "rna-seq") -> str:
    """
    Calls the mcpserver to run a new RNA-seq analysis pipeline.
    """
    endpoint = f"{MCPSERVER_URL}/pipeline/run"
    params = {
        "srr_list": srr_list,
        "fastq_dir": fastq_dir,
        "seq_mode": seq_mode,
        "species": species,
        "genome_version": genome_version,
        "seq_type": seq_type,
    }
    try:
        response = requests.post(endpoint, json=params)
        response.raise_for_status()
        return f"Successfully started pipeline. Response: {response.json()}"
    except requests.exceptions.RequestException as e:
        return f"Error calling mcpserver: {e}"

@tool(args_schema=TaskStatusInput)
def get_task_status(task_id: str) -> str:
    """
    Retrieves the status and results of a specific task from the mcpserver.
    """
    endpoint = f"{MCPSERVER_URL}/tasks/{task_id}"
    try:
        response = requests.get(endpoint)
        response.raise_for_status()
        return f"Task status: {response.json()}"
    except requests.exceptions.RequestException as e:
        return f"Error retrieving task status: {e}"

@tool(args_schema=UploadMetaInput)
def upload_meta_file(task_id: str, file_path: str) -> str:
    """
    Uploads a metadata file for a specific task to the mcpserver.
    """
    endpoint = f"{MCPSERVER_URL}/tasks/{task_id}/meta"
    if not os.path.exists(file_path):
        return f"Error: File not found at {file_path}"
    
    with open(file_path, 'rb') as f:
        files = {'file': (os.path.basename(file_path), f)}
        try:
            response = requests.post(endpoint, files=files)
            response.raise_for_status()
            return f"Successfully uploaded metadata file. Response: {response.json()}"
        except requests.exceptions.RequestException as e:
            return f"Error uploading metadata file: {e}"

@tool(args_schema=DefineGroupsInput)
def define_groups(task_id: str, control_group: List[str], experiment_group: List[str]) -> str:
    """
    Defines the control and experiment groups for a specific task on the mcpserver.
    """
    endpoint = f"{MCPSERVER_URL}/tasks/{task_id}/groups"
    params = {
        "control_group": control_group,
        "experiment_group": experiment_group,
    }
    try:
        response = requests.post(endpoint, json=params)
        response.raise_for_status()
        return f"Successfully defined groups. Response: {response.json()}"
    except requests.exceptions.RequestException as e:
        return f"Error defining groups: {e}"

@tool(args_schema=StartDEAnalysisInput)
def start_de_analysis(task_id: str) -> str:
    """
    Starts the downstream differential expression analysis for a specific task.
    """
    endpoint = f"{MCPSERVER_URL}/tasks/{task_id}/analyze"
    try:
        response = requests.post(endpoint)
        response.raise_for_status()
        return f"Successfully started DE analysis. Response: {response.json()}"
    except requests.exceptions.RequestException as e:
        return f"Error starting DE analysis: {e}"

@tool(args_schema=GetAnalysisResultsInput)
def get_analysis_results(task_id: str) -> str:
    """
    Retrieves the downstream analysis results for a specific task.
    This tool should be used after the DE analysis is completed.
    The results will contain URLs to the generated plots.
    """
    endpoint = f"{MCPSERVER_URL}/tasks/{task_id}"
    try:
        response = requests.get(endpoint)
        response.raise_for_status()
        task_data = response.json()
        # Assuming the results are stored in the 'results' field of the task
        analysis_results = task_data.get("results", {}).get("de_analysis", {})
        if not analysis_results:
            return "DE analysis results not found. Please ensure the analysis has completed successfully."
        return f"DE analysis results: {analysis_results}"
    except requests.exceptions.RequestException as e:
        return f"Error retrieving analysis results: {e}"

# Define a safe base path for file operations
SAFE_BASE_PATH = os.path.abspath("data")

@tool
def list_files(path: str = ".") -> str:
    """
    Lists the files and directories in a specified path relative to the data directory.
    Use '.' to list the contents of the base data directory.
    Example: list_files(path="bam") will list contents of the 'data/bam' directory.
    This tool is for exploring the available data files.
    """
    try:
        # Prevent directory traversal attacks
        full_path = os.path.abspath(os.path.join(SAFE_BASE_PATH, path))
        if not full_path.startswith(SAFE_BASE_PATH):
            return "Error: Access denied. You can only list files within the data directory."

        if not os.path.exists(full_path) or not os.path.isdir(full_path):
            return f"Error: Path '{path}' does not exist or is not a directory."

        files = os.listdir(full_path)
        if not files:
            return f"The directory '{path}' is empty."
        
        return f"Contents of '{path}':\n" + "\n".join(files)
    except Exception as e:
        return f"An error occurred: {e}"

# --- Main Agent Logic ---
def main():
    """
    Main function to initialize and run the agent.
    """
    # 1. Initialize LLM
    llm = ChatOpenAI(
        temperature=0,
        model_name=os.environ.get("OPENAI_MODEL_NAME"),
        openai_api_key=os.environ.get("OPENAI_API_KEY"),
        openai_api_base=os.environ.get("OPENAI_API_BASE"),
    )

    # 2. Define tools
    tools = [run_rnaseq_pipeline, get_task_status, upload_meta_file, define_groups, start_de_analysis, get_analysis_results, list_files]

    # 3. Create prompt
    prompt_template = """
    You are an expert bioinformatics assistant.
    Your job is to help users run RNA-seq analysis pipelines and analyze the results.
    You have access to the following tools:
    {tools}

    **Workflow:**

    1.  **Run Pipeline:** When a user asks to run a pipeline, use the `run_rnaseq_pipeline` tool.
    2.  **Monitor Status:** After starting the pipeline, you will get a `task_id`. Use the `get_task_status` tool with this `task_id` to monitor the progress.
    3.  **Downstream Analysis:** Once the pipeline is complete, guide the user through the downstream analysis:
        *   Ask the user for the path to their metadata file and use the `upload_meta_file` tool.
        *   Ask the user to define the control and experiment groups and use the `define_groups` tool.
        *   Use the `start_de_analysis` tool to begin the analysis.
    4.  **Get Results:** After the DE analysis is complete, use the `get_analysis_results` tool to retrieve the results, such as plots and tables.
    5.  **Interpret and Report:** Once you have the analysis results (which will include URLs for plots), your final step is to interpret them. Describe what the plots (e.g., volcano plot, heatmap) show and provide a summary of the key findings. Present this as a final, comprehensive report to the user.

    **Important:** Always keep the user informed about the progress and the next steps.

    Use the following format:
    
    Question: the input question you must answer
    Thought: you should always think about what to do
    Action: the action to take, should be one of [{tool_names}]
    Action Input: the input to the action
    Observation: the result of the action
    ... (this Thought/Action/Action Input/Observation can repeat N times)
    Thought: I now know the final answer
    Final Answer: the final answer to the original input question
    
    Begin!
    
    Question: {input}
    Thought: {agent_scratchpad}
    """
    prompt = ChatPromptTemplate.from_template(prompt_template)

    # 4. Create agent
    agent = create_react_agent(llm, tools, prompt)

    # 5. Create agent executor
    agent_executor = AgentExecutor(agent=agent, tools=tools, verbose=True)

    # 6. Run agent with a sample query
    # print("Starting agent...")
    # query = "Please run an RNA-seq analysis for human samples using the hg38 genome."
    # result = agent_executor.invoke({"input": query})
    # print(f"Agent finished with result: {result}")

tools = [run_rnaseq_pipeline, get_task_status, upload_meta_file, define_groups, start_de_analysis, get_analysis_results]

prompt_template = """
You are an expert bioinformatics assistant.
Your job is to help users run RNA-seq analysis pipelines and analyze the results.
You have access to the following tools:
{tools}

**Workflow:**

1.  **Explore Data:** If the user asks what data is available, how many files are in a directory, or wants to see the contents of a directory, use the `list_files` tool. This helps understand the available data before starting a pipeline.
2.  **Run Pipeline:** When a user asks to run a pipeline, use the `run_rnaseq_pipeline` tool.
3.  **Monitor Status:** After starting the pipeline, you will get a `task_id`. Use the `get_task_status` tool with this `task_id` to monitor the progress.
4.  **Downstream Analysis:** Once the pipeline is complete, guide the user through the downstream analysis:
    *   Ask the user for the path to their metadata file (you can use `list_files` to help them find it) and use the `upload_meta_file` tool.
    *   Ask the user to define the control and experiment groups and use the `define_groups` tool.
    *   Use the `start_de_analysis` tool to begin the analysis.
5.  **Get Results:** After the DE analysis is complete, use the `get_analysis_results` tool to retrieve the results, such as plots and tables.
6.  **Interpret and Report:** Once you have the analysis results (which will include URLs for plots), your final step is to interpret them. Describe what the plots (e.g., volcano plot, heatmap) show and provide a summary of the key findings. Present this as a final, comprehensive report to the user.

**Important:** Always keep the user informed about the progress and the next steps.

Use the following format:

Question: the input question you must answer
Thought: you should always think about what to do
Action: the action to take, should be one of [{tool_names}]
Action Input: the input to the action
Observation: the result of the action
... (this Thought/Action/Action Input/Observation can repeat N times)
Thought: I now know the final answer
Final Answer: the final answer to the original input question

Begin!

Question: {input}
Thought: {agent_scratchpad}
"""

if __name__ == "__main__":
    main()