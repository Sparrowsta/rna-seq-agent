import os
from langchain_core.tools import tool
from pydantic.v1 import BaseModel, Field # LangChain's tool decorator requires Pydantic v1 for argument validation.

class DirectoryPath(BaseModel):
    """Input model for the list_directory tool."""
    path: str = Field(description="The path to the directory you want to inspect.")

@tool(args_schema=DirectoryPath)
def list_directory(path: str) -> str:
    """
    Lists the contents (files and subdirectories) of a specified directory.
    """
    try:
        if not os.path.exists(path):
            return f"Error: Directory not found at '{path}'."
        if not os.path.isdir(path):
            return f"Error: The path '{path}' is not a directory."
        contents = os.listdir(path)
        if not contents:
            return f"The directory '{path}' is empty."
        return "Directory contents:\n" + "\n".join(contents)
    except PermissionError:
        return f"Error: Permission denied to access the directory '{path}'."
    except Exception as e:
        return f"An unexpected error occurred: {str(e)}"