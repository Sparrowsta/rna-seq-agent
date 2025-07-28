#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- 1. Parameter Definitions ---
// These parameters are the public API of the pipeline, controlled by the agent.

// Input data
params.reads = null           // Glob pattern for existing FASTQ files, e.g., 'data/fastq/SRR*/*_{1,2}.fastq.gz'
params.srr_ids = null         // Comma-separated list of SRR IDs to download, e.g., "SRR123,SRR456"

// Genome resources
params.fasta = null           // Path to genome FASTA file
params.gtf = null             // Path to gene annotation GTF file
params.species = "unknown"
params.genome_version = "unknown"

// Output directory
params.outdir = './results'

// Tool-specific parameters
params.fastp_q_val = 20
params.fastp_unqual_pct = 40
params.fastp_len_req = 50
params.fc_is_paired_end = true
params.fc_strand_spec = 2

// Orchestration parameters (controlled by the agent's plan)
params.run_genome_prep = false
params.run_download_fastq = false
params.run_analysis = false


// --- 2. Module Imports ---
// Include the building blocks of our pipeline.
include { GENOME_PREP } from './modules/genome_prep'
include { DOWNLOAD_FASTQ } from './modules/download_fastq'
include { ANALYSIS } from './modules/analysis'


// --- 3. Workflow Orchestration ---
workflow {
    log.info """
    --- RNA-Seq Pipeline ---
    Run Flags:
      - Genome Prep:    ${params.run_genome_prep}
      - FASTQ Download: ${params.run_download_fastq}
      - Analysis:       ${params.run_analysis}
    Input FASTQ:      ${params.reads ?: 'None'}
    Input SRR IDs:    ${params.srr_ids ?: 'None'}
    Genome FASTA:     ${params.fasta ?: 'None'}
    ---------------------------
    """

    // Channel for the STAR index. It will be populated either by the
    // GENOME_PREP module or by an existing path provided by the agent.
    ch_star_index = Channel.empty()

    // Channel for the input reads. It will be populated either by the
    // DOWNLOAD_FASTQ module or by a glob pattern of existing files.
    ch_reads_for_analysis = Channel.empty()

    // --- Module Execution Logic ---

    // Step 1: Genome Preparation (Optional)
    if (params.run_genome_prep) {
        if (!params.fasta || !params.gtf) {
            error "To run genome prep, --fasta and --gtf must be provided."
        }
        fasta_ch = Channel.fromPath(params.fasta)
        gtf_ch = Channel.fromPath(params.gtf)
        GENOME_PREP(fasta_ch, gtf_ch)
        ch_star_index = GENOME_PREP.out.star_index
    } else if (params.run_analysis) {
        // If not running prep but running analysis, the index must already exist.
        // The agent is responsible for providing the correct path via params.fasta.
        def existing_index_path = file(params.fasta).parent.resolve('star_index')
        if (!existing_index_path.exists() || !existing_index_path.isDirectory()) {
            error "Analysis requires a STAR index. Path not found: ${existing_index_path}"
        }
        ch_star_index = Channel.fromPath(existing_index_path)
    }

    // Step 2: FASTQ Download (Optional)
    if (params.run_download_fastq) {
        if (!params.srr_ids) {
            error "To run FASTQ download, --srr_ids must be provided."
        }
        srr_ids_ch = Channel.from(params.srr_ids)
        DOWNLOAD_FASTQ(srr_ids_ch)
        ch_reads_for_analysis = DOWNLOAD_FASTQ.out.reads
    } else if (params.run_analysis) {
        // If not running download but running analysis, the reads must already exist.
        // The agent is responsible for providing the correct glob pattern via params.reads.
        if (!params.reads) {
            error "Analysis requires input reads. Please provide a glob pattern via --reads."
        }
        ch_reads_for_analysis = Channel
            .fromFilePairs(params.reads, size: -1) { file ->
                file.name.replaceAll(/_?[12]\.fastq\.gz$|_\.fastq\.gz$/, "")
            }
            .map { sample_id, files ->
                tuple(sample_id, files.sort())
            }
    }

    // Step 3: Core Analysis (Optional)
    if (params.run_analysis) {
        if (!params.gtf) {
            error "Analysis requires a GTF file. Please provide it via --gtf."
        }
        gtf_ch_for_analysis = Channel.fromPath(params.gtf)
        ANALYSIS(ch_reads_for_analysis, ch_star_index, gtf_ch_for_analysis)
    }
}
