#!/bin/bash -ue
# Create HISAT2 index if not exists
if [ ! -f "genome.fa.1.ht2" ]; then
    echo "Building HISAT2 index..."
    hisat2-build "genome.fa" "genome.fa"
fi

echo "Aligning test..."
hisat2 -x "genome.fa" -U "test.fastq.gz" --threads 2         | samtools view -@ 2 -Sb -         | samtools sort -@ 2 -o "test_aligned.bam" -

samtools index "test_aligned.bam"
echo "Alignment complete for test"
