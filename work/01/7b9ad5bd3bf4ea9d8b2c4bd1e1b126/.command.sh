#!/bin/bash -ue
fastqc --quiet --threads 2 --outdir . "sample1.fastq.gz"
