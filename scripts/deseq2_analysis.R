#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2) {
    stop("Usage: Rscript deseq2_analysis.R <counts_file> <output_file>")
}

counts_file <- args[1]
output_file <- args[2]

cat("Running DESeq2 analysis...\n")
cat("Input file:", counts_file, "\n")
cat("Output file:", output_file, "\n")

# Load required libraries
if(!require(DESeq2)) {
    stop("DESeq2 package not installed")
}

# Here you would normally:
# 1. Read the count matrix from featureCounts output
# 2. Create sample metadata
# 3. Run DESeq2 analysis
# 4. Save results

# For now, create a dummy result
dummy_results <- data.frame(
    gene = paste0("gene", 1:100),
    log2FoldChange = rnorm(100),
    pvalue = runif(100),
    padj = runif(100)
)

write.csv(dummy_results, output_file, row.names = FALSE)
cat("Analysis complete. Results saved to:", output_file, "\n")

