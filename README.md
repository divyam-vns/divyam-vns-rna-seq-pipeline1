# divyam-vns-rna-seq-pipeline1
# RNA-seq Analysis Pipeline (FASTQ → Differential Expression)

This repository contains a full RNA-seq pipeline for biologists: from raw FASTQ files to differential gene expression results.  
It uses `FastQC`, `fastp`, `STAR`, `featureCounts`, and `DESeq2`.

---

## 🧬 Overview
1. **Quality Control (FastQC)**
2. **Trimming (fastp)**
3. **Alignment (STAR 2-pass)**
4. **Counting (featureCounts)**
5. **Differential Expression (DESeq2)**

---

## ⚙️ Requirements
- Linux / macOS
- Conda or Mamba
- 8+ threads recommended

Install dependencies:
```bash
mamba env create -f environment.yml
mamba activate rnaseq
