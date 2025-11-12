#!/bin/bash -ue
fastqc --quiet --threads 2 --outdir . "test.fastq.gz"
