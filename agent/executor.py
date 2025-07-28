import subprocess
import threading
from pathlib import Path

def _stream_reader(stream, log_file_path):
    """Reads a stream line by line and appends to a log file."""
    try:
        with open(log_file_path, 'ab') as f: # Open in append-binary mode
            for line in iter(stream.readline, b''):
                f.write(line)
                f.flush()
    finally:
        stream.close()

def execute_with_logging(command: list[str], log_file_path: Path, cwd: Path, env: dict | None = None) -> subprocess.Popen:
    """
    Executes a command, streams its stdout and stderr to a log file,
    and returns the process handle.

    Args:
        command: The command to execute, as a list of strings.
        log_file_path: The path to the file where logs should be written.
        cwd: The working directory for the command.
        env: A dictionary of environment variables to set for the process.

    Returns:
        A subprocess.Popen object representing the running process.
    """
    log_file_path.parent.mkdir(parents=True, exist_ok=True)
    
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, # Redirect stderr to stdout
        cwd=str(cwd),
        env=env,
        bufsize=1, # Line-buffered
    )

    # Start a thread to read the process's output stream
    # and write it to the log file.
    threading.Thread(
        target=_stream_reader,
        args=(process.stdout, log_file_path),
        daemon=True
    ).start()

    return process