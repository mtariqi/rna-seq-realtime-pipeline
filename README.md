[![RNA-seq Pipeline CI](https://github.com/mtariqi/rna-seq-realtime-pipeline/actions/workflows/test.yml/badge.svg)](https://github.com/mtariqi/rna-seq-realtime-pipeline/actions/workflows/test.yml)
![Conda](https://img.shields.io/badge/Conda-ready-blue)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-green)
![License](https://img.shields.io/badge/License-MIT-yellow)
![Version](https://img.shields.io/badge/version-1.0-brightgreen)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)




# ⚡ Real-time RNA-seq Pipeline
**Author:** MD Tariqul Islam (Tariq)  
**GitHub:** [@mtariqi](https://github.com/mtariqi)  
**LinkedIn:** [www.linkedin.com/in/mdtariqulscired](https://www.linkedin.com/in/mdtariqulscired) 
**License:** MIT • **Version:** 1.0 • **Last Updated:** November 2025  



---

## 📘 Overview
This repository implements a **modular, streaming-aware RNA-seq workflow** for *real-time Nanopore* or *Illumina* sequencing analysis.  
The pipeline continuously monitors an input directory for new reads, performs **on-the-fly basecalling**, **alignment**, and **gene-level quantification**, and reports **fusion events** as data arrive.

It is designed to operate both in:
- 🧑‍🔬 **Interactive mode** (local/academic HPC)  
- ☁️ **Cloud/HPC environments** (AWS Batch, SLURM, Google Cloud)  

Built entirely in **Nextflow DSL2** with Conda-based reproducibility, this pipeline supports **streaming bioinformatics** and **low-latency diagnostics** applications.

---

## 🔑 Key Features

| Category | Description |
|-----------|--------------|
| 🧬 **Dorado Basecalling** | GPU-optimized Oxford Nanopore basecaller |
| 🧭 **Streaming Alignment** | On-the-fly spliced alignment using [`minimap2`](https://github.com/lh3/minimap2) |
| 📊 **Incremental Quantification** | Gene-level counting via [`featureCounts`](https://subread.sourceforge.net) |
| 🔍 **Fusion Detection** | Continuous fusion transcript detection using [`JAFFAL`](https://github.com/Oshlack/JAFFAL) |
| 🧠 **Automation** | Real-time file watcher triggers the pipeline automatically as FASTQ files appear |
| ⚙️ **Nextflow DSL2 Modularity** | Scalable, maintainable processes and channels |
| 🧪 **CI Integration** | Automated GitHub Actions test with environment validation and dry-run simulation |
| 🔐 **Reproducibility** | Environment-locked `environment.yml` for fully deterministic runs |

---

## 🧩 Project Structure
```
rna-seq-realtime-pipeline/
├── .github/workflows/ # CI automation (Nextflow validation + dry-run)
├── data/ # Example FASTQ input data
├── results/ # Pipeline output
├── scripts/ # Helper scripts & utilities
├── watcher/ # File watcher for real-time streaming mode
├── nextflow.config # Runtime configuration
├── main.nf # Core Nextflow workflow (DSL2)
├── environment.yml # Conda environment definition
└── README.md # Project documentation
```
---

## ⚙️ Installation

### 1️⃣ Clone the repository
```bash
git clone https://github.com/mtariqi/rna-seq-realtime-pipeline.git
cd rna-seq-realtime-pipeline

### 2️⃣ Create Conda environment
conda env create -f environment.yml
conda activate rna_realtime_env

### 3️⃣ Run a quick test (dry-run)

nextflow run main.nf -stub-run

### 4️⃣ Enable real-time watcher
bash watcher/watch_and_process.sh


This continuously monitors your FASTQ directory and launches analysis as new data arrive.

🧪 Continuous Integration (CI)

The repository includes an automated workflow using GitHub Actions to:

Build and validate the rna_realtime_env Conda environment

Install and verify Nextflow, FastQC, MultiQC, and Minimap2

Perform a Nextflow dry-run simulation (-stub-run)

Upload diagnostic logs as CI artifacts

You can view the live CI status under the Actions tab
or by following this badge:
👉

📈 Example Output
Module	Output	Description
FastQC	results/fastqc/	Quality metrics per read
Minimap2	results/alignment/	Spliced alignments (BAM/SAM)
FeatureCounts	results/counts.txt	Gene-level counts
JAFFAL	results/fusions/	Fusion gene candidates
MultiQC	results/multiqc_report.html	Aggregated QC report
🧮 Reproducibility

Workflow Language: Nextflow DSL2

Environment Manager: Conda

Container Support: Docker / Singularity (optional)

Validation: Continuous Integration via GitHub Actions

For full reproducibility, freeze all package versions before deployment:

conda env export > environment.lock.yml

🧠 Citation & References

Islam, M. T. (2025). *Real-time RNA-seq Pipeline (v1.0.0)* [Software]. Zenodo.  
https://doi.org/10.5281/zenodo.17603512

Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18): 3094–3100.

Davidson N. et al. (2022). JAFFAL: Fusion gene detection from long-read transcriptome data. Bioinformatics, 38(6): 1577–1583.

Oxford Nanopore Technologies (2024). Dorado Basecaller.

Subread Team (2014). featureCounts: efficient read summarization program.

🤝 Contributing

Contributions and pull requests are welcome!
Please fork the repo, create a feature branch, and submit a pull request.

git checkout -b feature/new-module
git commit -m "Add new module"
git push origin feature/new-module

🧾 How to Cite This Repository

If you use this pipeline in your research, please cite it as follows:

APA (Recommended):

Islam, M. T. (2025). Real-time RNA-seq Pipeline: A Nextflow DSL2 framework for live transcriptomic analysis using Dorado, Minimap2, FeatureCounts, and JAFFAL. Zenodo. https://doi.org/10.5281/zenodo.xxxxxxx

BibTeX:
```
@software{islam_rnaseq_realtime_2025,
  author       = {Islam, MD Tariqul},
  title        = {Real-time RNA-seq Pipeline: A Nextflow DSL2 framework for live transcriptomic analysis using Dorado, Minimap2, FeatureCounts, and JAFFAL},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.xxxxxxx},
  url          = {https://doi.org/10.5281/zenodo.xxxxxxx}
}
```
---------------------------------------------------------------------
```
📧 Contact

For technical inquiries or collaborations:
📩 tariqul@scired.com

🌐 LinkedIn: www.linkedin.com/in/mdtariqulscired

💻 GitHub: https://github.com/mtariqi
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17603512.svg)](https://doi.org/10.5281/zenodo.17603512)
📧 Contact

🧾 License

This project is released under the MIT License
.

“Real-time RNA-seq analysis is not just computation — it’s precision medicine in motion.”
— MD Tariqul Islam (Tariq)














