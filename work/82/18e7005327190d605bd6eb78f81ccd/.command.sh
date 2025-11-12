#!/bin/bash -ue
echo "Quantifying sample1..."
featureCounts -a "annotation.gtf" -o "sample1_gene_counts.txt" -T 2 "sample1_aligned.bam"
echo "Quantification complete for sample1"
