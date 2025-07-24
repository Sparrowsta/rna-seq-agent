process PUBLISH_RESULTS {
    // This process does not generate files, so it doesn't need a publishDir
    // It's a utility process to send data out
    
    input:
    val task_id
    val step_name
    path result_file

    script:
    """
    #!/bin/bash
    # This script sends step results to the MCP server.
    # It handles both JSON and plain text result files.

    # Exit if essential commands are missing
    if ! command -v curl &> /dev/null || ! command -v jq &> /dev/null; then
        echo "Error: curl and jq are required for publishing results." >&2
        exit 1
    fi

    # If no task_id is provided, we can't publish, so exit gracefully.
    # This allows the pipeline to run without server integration.
    if [ -z "${task_id}" ]; then
        echo "Warning: task_id is not set. Skipping result publication for step '${step_name}'."
        exit 0
    fi

    # Prepare the JSON payload.
    # If the result file is valid JSON, its content is used as the 'data' object.
    # Otherwise, the raw text content is wrapped inside a 'data.raw_content' field.
    if jq -e . >/dev/null 2>&1 < "${result_file}"; then
        # Valid JSON file
        echo "Result file '${result_file}' is valid JSON. Embedding directly."
        json_payload=\$(jq -n --arg step "${step_name}" --argjson data "\$(cat ${result_file})" \\
          '{step_name: \$step, data: \$data}')
    else
        # Not a JSON file, wrap content as a string
        echo "Result file '${result_file}' is not JSON. Wrapping content as raw text."
        json_payload=\$(jq -n --arg step "${step_name}" --arg data_raw "\$(cat ${result_file})" \\
          '{step_name: \$step, data: {raw_content: \$data_raw}}')
    fi

    # Send the data to the MCP server.
    # curl options:
    # -f: Fail silently (no output) on HTTP errors (4xx, 5xx)
    # -s: Silent mode
    # -S: Show error message on failure
    echo "Sending payload for step '${step_name}' to MCP server for task '${task_id}'..."
    curl -f -s -S -X POST "${params.mcp_server_url}/tasks/${task_id}/results" \\
         -H "Content-Type: application/json" \\
         --data "\$json_payload"
    
    echo "Publication for step '${step_name}' finished with exit code \$?."
    """
}

process PUBLISH_DE_RESULTS {
    tag "Publish DE results for ${task_id}"
    
    publishDir "${params.outdir}/de_analysis/${task_id}", mode: 'copy', overwrite: true

    input:
    val task_id
    path de_results // This will be the directory from DE_ANALYSIS

    output:
    path "."

    script:
    """
    #!/bin/bash
    # This process simply ensures that the results from the DE_ANALYSIS step
    # are published to the correct, task-specific directory.
    # The publishDir directive handles the actual file copying.
    echo "Publishing DE analysis results for task ${task_id} to ${params.outdir}/de_analysis/${task_id}"
    # Create a dummy file to satisfy the output requirement if needed,
    # though publishing the directory itself is the main goal.
    touch publish_de_complete.txt
    """
}