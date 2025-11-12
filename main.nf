#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * ============================================================
 * RNA-seq Realtime Pipeline (FastQC → MultiQC)
 * Author: MD Tariqul Islam (Tariq)
 * ============================================================
 */

// ---------------- PARAMETERS ----------------
params.reads  = "data/test.fastq.gz"
params.outdir = "results"

// ---------------- CHANNELS ----------------
Channel
    .fromPath(params.reads)
    .ifEmpty { error "❌ No FASTQ file found in ${params.reads}" }
    .set { reads_ch }

// ---------------- PROCESS: FASTQC ----------------
process FASTQC {
    tag "$reads.baseName"

    input:
    path reads

    output:
    tuple val(reads.baseName), path("*.zip"), path("*.html")

    script:
    """
    echo "🔍 Running FastQC on ${reads}"
    fastqc --threads 2 --outdir . ${reads}
    """
}

// ---------------- PROCESS: MULTIQC ----------------
process MULTIQC {
    tag "multiqc"

    input:
    tuple val(sample), path(zips), path(htmls)

    output:
    path "multiqc_report.html"

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    echo "📊 Aggregating FastQC results with MultiQC"
    multiqc . -o .
    """
}

// ---------------- WORKFLOW ----------------
workflow {
    // Step 1: run FASTQC
    fastqc_results = FASTQC(reads_ch)

    // Step 2: run MULTIQC on FASTQC results
    MULTIQC(fastqc_results)
}
