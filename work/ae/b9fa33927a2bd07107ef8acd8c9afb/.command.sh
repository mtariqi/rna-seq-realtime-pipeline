#!/bin/bash -ue
echo "Quantifying test..."
featureCounts -a "annotation.gtf" -o "test_gene_counts.txt" -T 2 "test_aligned.bam"
echo "Quantification complete for test"
