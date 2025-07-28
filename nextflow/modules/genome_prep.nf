//
// Module: Genome Preparation
// Responsibility: Unconditionally generate a STAR index from a given FASTA and GTF file.
//
nextflow.enable.dsl=2

workflow GENOME_PREP {
    take:
        fasta_ch // channel: /path/to/genome.fa
        gtf_ch   // channel: /path/to/genes.gtf

    main:
        STAR_GENOME_GENERATE(fasta_ch, gtf_ch)

    emit:
        star_index = STAR_GENOME_GENERATE.out.star_index // channel: /path/to/star_index
}

process STAR_GENOME_GENERATE {
    // The publishDir now points to the genomes directory.
    publishDir "data/genomes/${params.species}/${params.genome_version}/", mode: 'copy'
    
    label 'large_mem_process'

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: star_index

    script:
    """
    source activate align_env
    mkdir star_index
    STAR --runMode genomeGenerate \\
         --genomeDir star_index \\
         --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile ${gtf} \\
         --runThreadN ${task.cpus}
    """
}