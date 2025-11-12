#!/bin/bash -ue
echo "🔍 Running FastQC on test.fastq.gz"
fastqc --quiet --threads 2 --outdir . test.fastq.gz
