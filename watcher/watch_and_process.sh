#!/usr/bin/env bash
# ==========================================================
#  Real-time RNA-seq Watcher
#  Author: MD Tariqul Islam (Tariq)
#  GitHub: https://github.com/mtariqi/rna-seq-realtime-pipeline
# ==========================================================
#  Description:
#   Continuously watches a directory for new FASTQ/FASTQ.GZ files
#   and triggers the Nextflow pipeline automatically when
#   new data appear (ideal for Nanopore or Illumina live runs).
# ==========================================================

# -------- CONFIGURATION -----------------------------------
WATCH_DIR="data/live"              # Folder to monitor
LOG_DIR="logs"                     # Folder for watcher logs
PIPELINE="main.nf"                 # Nextflow pipeline file
CONFIG="nextflow.config"           # Nextflow config file
ENV_NAME="rna_realtime_env"        # Conda environment
SLEEP_INTERVAL=10                  # Seconds between checks
# -----------------------------------------------------------

mkdir -p "$WATCH_DIR" "$LOG_DIR"

echo "👀 Starting watcher on: $WATCH_DIR"
echo "🕐 Checking every $SLEEP_INTERVAL seconds for new FASTQ files..."
echo "--------------------------------------------------------------"

# Track previously seen files
touch "$LOG_DIR/.seen_files.txt"

# Activate conda environment if available
if command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME" || echo "⚠️  Conda environment not activated; using system env"
fi

while true; do
    # Find new FASTQ files not yet processed
    NEW_FILES=$(find "$WATCH_DIR" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | grep -v -F -f "$LOG_DIR/.seen_files.txt")

    if [ -n "$NEW_FILES" ]; then
        echo "🧬 New FASTQ files detected:"
        echo "$NEW_FILES"
        echo "$NEW_FILES" >> "$LOG_DIR/.seen_files.txt"

        echo "🚀 Launching Nextflow pipeline..."
        TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
        nextflow run "$PIPELINE" -c "$CONFIG" -resume \
            > "$LOG_DIR/run_${TIMESTAMP}.log" 2>&1 &

        echo "✅ Pipeline triggered at $TIMESTAMP"
        echo "--------------------------------------------------------------"
    fi

    sleep "$SLEEP_INTERVAL"
done
