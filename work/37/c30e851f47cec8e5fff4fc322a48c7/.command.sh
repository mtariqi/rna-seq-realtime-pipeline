#!/bin/bash -ue
fastqc --quiet --threads 2 --outdir . "sample2.fastq.gz"
