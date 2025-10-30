# Reference Genome Folder

This folder should contain all genome reference and annotation files required for RNA-seq alignment and quantification.

## Expected files
- `GRCh38.fa` — Human genome FASTA file (primary assembly)
- `GRCh38.fa.fai` — FASTA index file
- `GRCh38.gtf` — Gene annotation file
- `HISAT2_index/` — Folder with HISAT2 genome index files

## Download Instructions
You can download these files using the commands below:

```bash
# Create genome directory
mkdir -p config/genomes && cd config/genomes

# Download GRCh38 FASTA and GTF from Ensembl
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz

# Unzip them
gunzip *.gz

# Build HISAT2 index
hisat2-build GRCh38.dna.primary_assembly.fa HISAT2_index/GRCh38
