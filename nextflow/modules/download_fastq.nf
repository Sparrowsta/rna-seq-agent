//
// Module: Download FASTQ from SRA
// Responsibility: Unconditionally download FASTQ files for a given list of SRR IDs.
//
nextflow.enable.dsl=2

workflow DOWNLOAD_FASTQ {
    take:
        srr_ids_ch // channel: ['SRR123', 'SRR456']

    main:
        srr_ids_ch
            .splitText()
            .map { it.trim() }
            .set { srr_id_ch_for_download }

        FASTERQ_DUMP(srr_id_ch_for_download)

    emit:
        reads = FASTERQ_DUMP.out.reads // channel: [ [id, [/path/to/r1, /path/to/r2]], ... ]
}

process FASTERQ_DUMP {
    // Publish to a central 'fastq' directory, organized by SRR ID.
    publishDir "${params.outdir}/fastq/${srr_id}", mode: 'copy'

    label 'default_process'

    input:
    val srr_id

    output:
    tuple val(srr_id), path("*.fastq.gz"), emit: reads

    script:
    """
    source activate qc_env
    fasterq-dump ${srr_id} --split-files --gzip -O . -p
    """
}