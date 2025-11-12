#!/bin/bash -ue
echo "Quantifying sample2..."
featureCounts -a "annotation.gtf" -o "sample2_gene_counts.txt" -T 2 "sample2_aligned.bam"
echo "Quantification complete for sample2"
