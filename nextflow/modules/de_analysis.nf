process DE_ANALYSIS {
    // Define a new label for processes that use the 'de_env' conda environment
    label 'de_analysis_process'

    // The output will be a directory, so we can publish it easily
    publishDir "${params.outdir}/de_analysis/${task_id}", mode: 'copy'

    input:
    val task_id
    path counts_file
    path meta_file
    val control_group
    val experiment_group

    output:
    path "results"

    script:
    """
    #!/bin/bash
    # This process orchestrates the downstream differential expression analysis.

    # Activate the correct conda environment
    source activate ngs_env

    # Run the Python script that in turn calls the R script
    python ${baseDir}/differential_analysis.py \\
        --counts_file ${counts_file} \\
        --meta_data '${meta_file.text}' \\
        --groups '{"control_group": ["${control_group.join('","')}"], "experiment_group": ["${experiment_group.join('","')}"]}' \\
        --output_dir results \\
        --task_id ${task_id}
    """
}