import asyncio
import subprocess
import uuid
from typing import Dict, Any

class NextflowManager:
    """
    Manages Nextflow workflow processes.
    """
    def __init__(self):
        self.tasks: Dict[str, Dict[str, Any]] = {}

    async def run_workflow(self, command: list[str]) -> str:
        """
        Runs a Nextflow workflow in the background.

        Args:
            command: The command to execute.

        Returns:
            The task ID.
        """
        task_id = str(uuid.uuid4())
        process = await asyncio.create_subprocess_exec(
            *command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        self.tasks[task_id] = {
            "process": process,
            "status": "running",
            "log": [],
            "task_id": task_id
        }
        asyncio.create_task(self._log_monitor(task_id))
        return task_id

    async def _log_monitor(self, task_id: str):
        """
        Monitors the stdout and stderr of a running process.
        """
        process = self.tasks[task_id]["process"]
        while process.returncode is None:
            if process.stdout:
                line = await process.stdout.readline()
                if line:
                    self.tasks[task_id]["log"].append(line.decode().strip())
            await asyncio.sleep(0.1)
        
        # After process finishes, update status
        if process.returncode == 0:
            self.tasks[task_id]["status"] = "completed"
        else:
            self.tasks[task_id]["status"] = "failed"
            if process.stderr:
                stderr_output = await process.stderr.read()
                self.tasks[task_id]["log"].append(stderr_output.decode().strip())


    def get_task_status(self, task_id: str) -> dict:
        """
        Gets the status of a specific task.

        Args:
            task_id: The ID of the task.

        Returns:
            A dictionary with the task's status and log.
        """
        if task_id not in self.tasks:
            return {"error": "Task not found"}
        
        task_info = self.tasks[task_id]
        return {
            "task_id": task_id,
            "status": task_info["status"],
            "log": task_info["log"][-100:]  # Return last 100 lines of log
        }

    def cancel_task(self, task_id: str) -> dict:
        """
        Cancels a running task.

        Args:
            task_id: The ID of the task.

        Returns:
            A dictionary with the result of the cancellation.
        """
        if task_id not in self.tasks:
            return {"error": "Task not found"}

        task_info = self.tasks[task_id]
        if task_info["status"] == "running":
            task_info["process"].terminate()
            task_info["status"] = "cancelled"
            return {"status": "cancellation request sent"}
        else:
            return {"status": f"task is already {task_info['status']}"}
