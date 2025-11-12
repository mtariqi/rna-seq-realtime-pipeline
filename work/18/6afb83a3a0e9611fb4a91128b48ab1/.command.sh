#!/bin/bash -ue
echo "📊 Aggregating FastQC results with MultiQC"
multiqc . -o .
