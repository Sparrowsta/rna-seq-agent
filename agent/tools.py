import os
import json
from typing import Dict, Any
from langchain_core.tools import tool
from pydantic.v1 import BaseModel, Field

class DirectoryPath(BaseModel):
    """Input model for the list_directory tool."""
    path: str = Field(description="The path to the directory you want to inspect.")

class ReadGenomeConfigArgs(BaseModel):
    """Input model for the read_genomes_config tool."""
    config_path: str = Field(description="The path to the genome configuration file.")

class UpdateGenomeConfigArgs(BaseModel):
    """Input model for the update_genomes_config tool."""
    config_path: str = Field(description="The path to the genome configuration file.")
    new_config: Dict[str, Any] = Field(description="A dictionary containing the new genomes configuration to overwrite the existing one.")

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

@tool(args_schema=ReadGenomeConfigArgs)
def read_genomes_config(config_path: str) -> str:
    """
    Reads the contents of the specified genomes configuration file.
    """
    try:
        with open(config_path, 'r',encoding="utf-8") as f:
            return json.dumps(json.load(f), indent=4)
    except FileNotFoundError:
        return f"Error: The configuration file was not found at '{config_path}'."
    except json.JSONDecodeError:
        return f"Error: The configuration file at '{config_path}' is not a valid JSON."
    except Exception as e:
        return f"An unexpected error occurred: {str(e)}"

@tool(args_schema=UpdateGenomeConfigArgs)
def update_genomes_config(config_path: str, new_config: Dict[str, Any]) -> str:
    """
    Overwrites the entire specified genome configuration file with the provided configuration.
    """
    try:
        # Write the new config directly to the file, overwriting existing content
        with open(config_path, 'w', encoding="utf-8", errors='surrogateescape') as f:
            json.dump(new_config, f, indent=4)
        return f"Successfully overwrote genomes configuration in '{config_path}'."
    except TypeError as e:
        return f"Error: The provided configuration is not serializable to JSON. Details: {str(e)}"
    except Exception as e:
        return f"An unexpected error occurred while updating the configuration: {str(e)}"