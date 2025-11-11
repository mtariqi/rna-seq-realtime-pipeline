# Real-time RNA-seq Pipeline

**Author:** MD Tariqul Islam (Tariq)  
**GitHub:** [@mtariqi](https://github.com/mtariqi)

This repository provides a modular workflow for *real-time* Nanopore RNA-seq analysis, including:

- 🧬 **Dorado basecalling**
- 🧭 **Streaming alignment** with `minimap2`
- 📊 **Incremental gene counting** with `featureCounts`
- 🔍 **Fusion detection** with `JAFFAL`
- ⏱️ Designed for near-real-time analysis as data arrive

### Quick start
```bash
conda env create -f environment.yml
conda activate rna_realtime_env
bash watcher/watch_and_process.sh

Citation / References

Dorado (ONT) 2024

Li H. minimap2 Bioinformatics 2018

Davidson N. et al. JAFFAL Bioinformatics 2022



