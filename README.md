[![RNA-seq Pipeline CI](https://github.com/mtariqi/rna-seq-realtime-pipeline/actions/workflows/test.yml/badge.svg)](https://github.com/mtariqi/rna-seq-realtime-pipeline/actions/workflows/test.yml)
![Conda](https://img.shields.io/badge/Conda-ready-blue)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-success)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)



# ⚡ Real-time RNA-seq Pipeline

**Author:** MD Tariqul Islam (Tariq)  
**GitHub:** [@mtariqi](https://github.com/mtariqi)  
**Version:** 1.0 • **License:** MIT • **Last Updated:** November 2025  

---

## 📘 Overview
This repository implements a **modular, streaming-aware workflow** for **real-time Nanopore RNA-seq analysis**.  
It continuously monitors an input directory, basecalls new reads with **Dorado**, performs **on-the-fly alignment** using `minimap2`, and updates **gene-level counts** and **fusion detections** as data arrive.  
The pipeline is written in **Nextflow DSL2** with optional Bash + Python watchers for low-latency updates.

### 🔑 Key Features
- 🧬 **Dorado basecalling** — GPU-optimized, accurate ONT basecaller  
- 🧭 **Streaming alignment** — fast spliced alignment via `minimap2`  
- 📊 **Incremental gene counting** — real-time quantification with `featureCounts`  
- 🔍 **Fusion detection** — continuous structural transcript analysis using `JAFFAL`  
- ⏱️ **Near real-time feedback** — process reads as soon as they are generated  
- ☁️ **Portable** — run locally, on HPC (SLURM), or in the cloud (AWS Batch)

---

## 🧩 Workflow Diagram

```mermaid
graph TD
    A[Incoming FAST5/FASTQ (Nanopore run)] --> B[Dorado Basecalling]
    B --> C[minimap2 Alignment]
    C --> D[featureCounts Quantification]
    C --> E[JAFFAL Fusion Detection]
    D --> F[Real-time Expression Dashboard]
    E --> F

## ⚙️ Installation & Setup
1️⃣ Create the Environment
conda env create -f environment.yml
conda activate rna_realtime_env
2️⃣ Prepare Input
Place your live Nanopore output (FAST5 or FASTQ) under:
data/live/
3️⃣ Launch Real-time Processing
bash watcher/watch_and_process.sh

🧠 Configuration
```
| Parameter        | Description              | Default         |
| ---------------- | ------------------------ | --------------- |
| `params.in_dir`  | Input directory to watch | `data/live`     |
| `params.out_dir` | Output directory         | `results`       |
| `params.genome`  | Reference genome FASTA   | `ref/genome.fa` |
| `params.gtf`     | Annotation file          | `ref/genes.gtf` |
| `params.device`  | GPU ID for Dorado        | `cuda:0`        |
```
## 📊 Output Overview
```
| Folder                 | Contents                                 |
| ---------------------- | ---------------------------------------- |
| `results/basecalling/` | FASTQ files and Dorado logs              |
| `results/alignment/`   | BAM + index files from minimap2          |
| `results/counts/`      | Incremental count tables (featureCounts) |
| `results/fusions/`     | JAFFAL fusion calls (JSON / TSV)         |
| `results/reports/`     | MultiQC summaries and runtime stats      |
```
## 🧬 Core Tools & Citations
```
| Tool              | Purpose                             | Reference                                                                                                                                                                                                                                                               |
| ----------------- | ----------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Dorado (ONT)**  | GPU basecaller for Nanopore signals | Oxford Nanopore Technologies (2024). *Dorado basecaller*. [GitHub: nanoporetech/dorado](https://github.com/nanoporetech/dorado)                                                                                                                                         |
| **minimap2**      | Fast alignment of long reads        | Li H. (2018). *Minimap2: pairwise alignment for nucleotide sequences.* **Bioinformatics 34(18)**: 3094–3100. [https://doi.org/10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)                                                             |
| **featureCounts** | Read quantification                 | Liao Y., Smyth G.K., Shi W. (2014). *featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.* **Bioinformatics 30(7)**: 923–930.                                                                                          |
| **JAFFAL**        | Fusion detection from long reads    | Davidson N.M., Schroder J., Robinson M.D. (2022). *JAFFAL: detecting fusion genes from long-read transcriptome sequencing data.* **Bioinformatics 38(12)**: 3312–3318. [https://doi.org/10.1093/bioinformatics/btac321](https://doi.org/10.1093/bioinformatics/btac321) |
```

## 🖥️ HPC & Cloud Support
```
process {
  executor = 'slurm'
  cpus     = 8
  memory   = '32 GB'
  time     = '4h'
  clusterOptions = '--account=your_lab_account'
}
```
Then execute:
nextflow run main.nf -c nextflow_slurm.config
## 🧪 Reproducibility

Deterministic outputs via Nextflow DSL2 caching

Environment captured in environment.yml

Configurable resources per process (CPUs, RAM, GPU)

Re-run using:

nextflow run main.nf -resume

📜 Citation

If you use or adapt this repository, please cite:

Islam, M.T. (2025). Real-time RNA-seq Pipeline for Nanopore Data Processing and Incremental Analysis using Nextflow DSL2. GitHub: https://github.com/mtariqi/rna-seq-realtime-pipeline

🧩 Contact

📧 tariqul@scired.com

🔗 LinkedIn
 | GitHub

“Turning live Nanopore signals into biological insight — in real time.”















