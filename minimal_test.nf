#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads = "data/fastq/*.fastq.gz"

workflow {
    Channel.fromPath(params.reads)
        .map { file -> tuple(file.baseName, file) }
        .set { fastq_ch }
    
    FASTQC(fastq_ch)
    MULTIQC(FASTQC.out.zip)
}

process FASTQC {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.zip", emit: zip
    
    script:
    """
    fastqc --quiet --outdir . "$reads"
    """
}

process MULTIQC {
    tag "MultiQC"
    
    input:
    path fastqc_results
    
    output:
    path "multiqc_report.html"
    
    script:
    """
    multiqc . -f
    """
}
