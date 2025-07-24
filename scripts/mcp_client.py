import requests
import argparse
import json

def get_args():
    """Parses command-line arguments for the MCP client."""
    parser = argparse.ArgumentParser(description="MCP Client to launch NGS pipelines.")
    parser.add_argument("--srr_list", type=str, help="Path to the SRR list file.")
    parser.add_argument("--fastq_dir", type=str, help="Path to the directory with local fastq files.")
    parser.add_argument("--seq_mode", action="store_true", help="Flag for single-end sequencing (default is paired-end).")
    parser.add_argument("--species", type=str, required=True, help="Species to process (e.g., 'human', 'mouse').")
    parser.add_argument("--genome_version", type=str, required=True, help="Genome version (e.g., 'hg38').")
    return parser.parse_args()

def main():
    """Main function to send a pipeline request to the MCP server."""
    args = get_args()

    # Construct the request payload from command-line arguments
    payload = {
        "srr_list": args.srr_list,
        "fastq_dir": args.fastq_dir,
        "seq_mode": args.seq_mode,
        "species": args.species,
        "genome_version": args.genome_version,
    }

    # The server URL
    url = "http://localhost:8001/pipeline/run"

    print("--- Sending Pipeline Request to MCP Server ---")
    print(f"URL: {url}")
    print("Payload:")
    print(json.dumps(payload, indent=4))
    print("-" * 45)

    try:
        response = requests.post(url, json=payload)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)

        response_data = response.json()
        print("--- Server Response ---")
        print(f"Status Code: {response.status_code}")
        print("Response Body:")
        print(json.dumps(response_data, indent=4))
        print("-" * 25)
        print("Pipeline task successfully submitted.")

    except requests.exceptions.RequestException as e:
        print(f"\n--- Error ---")
        print(f"Failed to connect to MCP server: {e}")
        print("Please ensure the MCP server is running and accessible.")
    except Exception as e:
        print(f"\n--- An Unexpected Error Occurred ---")
        print(e)

if __name__ == "__main__":
    main()