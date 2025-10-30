#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# RNA-seq Pipeline: QC → Alignment → Counting → DE Analysis
# ============================================================

CONFIG="config/params.yml"

# --- Helper function to read YAML values ---
read_yaml() {
    yq e "$1" "$CONFIG"
}

# --- Load parameters from config/params.yml ---
PROJECT=$(read_yaml '.project_name')
THREADS=$(read_yaml '.threads')
OUTDIR=$(read_yaml '.outdir')
SAMPLES=$(read_yaml '.samples')

GENOME_FASTA=$(read_yaml '.genome.fasta')
GTF=$(read_yaml '.genome.gtf')
INDEX=$(read_yaml '.genome.hisat2_index')

# --- Create results directories ---
mkdir -p "$OUTDIR" "$OUTDIR"/{01_fastqc,02_multiqc,03_alignment,04_counts,05_deseq2,logs}

echo "==============================================="
echo "Starting RNA-seq Pipeline: $PROJECT"
echo "Threads: $THREADS"
echo "Output Directory: $OUTDIR"
echo "==============================================="

# ============================================================
# STEP 1. Quality Control
# ============================================================
if [ "$(read_yaml '.fastqc.enabled')" = "true" ]; then
    echo "[1/5] Running FastQC..."
    while IFS=, read -r SAMPLE CONDITION R1 R2; do
        [ "$SAMPLE" = "sample" ] && continue
        fastqc -t "$THREADS" -o "$OUTDIR/01_fastqc" "$R1" "$R2"
    done < "$SAMPLES"
fi

# MultiQC
if [ "$(read_yaml '.multiqc.enabled')" = "true" ]; then
    echo "[2/5] Running MultiQC..."
    multiqc "$OUTDIR/01_fastqc" -o "$OUTDIR/02_multiqc"
fi

# ============================================================
# STEP 2. Alignment with HISAT2
# ============================================================
if [ "$(read_yaml '.hisat2.enabled')" = "true" ]; then
    echo "[3/5] Running HISAT2 alignments..."
    while IFS=, read -r SAMPLE CONDITION R1 R2; do
        [ "$SAMPLE" = "sample" ] && continue
        hisat2 -p "$THREADS" --dta -x "$INDEX" -1 "$R1" -2 "$R2" \
            -S "$OUTDIR/03_alignment/${SAMPLE}.sam" \
            2> "$OUTDIR/logs/${SAMPLE}_hisat2.log"
        samtools sort -@ "$THREADS" -o "$OUTDIR/03_alignment/${SAMPLE}.bam" "$OUTDIR/03_alignment/${SAMPLE}.sam"
        samtools index "$OUTDIR/03_alignment/${SAMPLE}.bam"
        rm "$OUTDIR/03_alignment/${SAMPLE}.sam"
    done < "$SAMPLES"
fi

# ============================================================
# STEP 3. Counting with featureCounts
# ============================================================
if [ "$(read_yaml '.featureCounts.enabled')" = "true" ]; then
    echo "[4/5] Counting reads with featureCounts..."
    BAM_LIST=$(find "$OUTDIR/03_alignment" -name "*.bam" | tr '\n' ' ')
    featureCounts -T "$THREADS" -a "$GTF" -o "$OUTDIR/04_counts/gene_counts.txt" $BAM_LIST \
        -p -s 0 -g gene_id 2> "$OUTDIR/logs/featureCounts.log"
fi

# ============================================================
# STEP 4. Differential Expression (DESeq2)
# ============================================================
if [ "$(read_yaml '.deseq2.enabled')" = "true" ]; then
    echo "[5/5] Running DESeq2 (R-based differential analysis)..."
    Rscript - <<'EOF'
library(DESeq2)
library(tidyverse)

# Read featureCounts output
counts <- read.delim("results/04_counts/gene_counts.txt", comment.char="#")
rownames(counts) <- counts$Geneid
counts <- counts[, -(1:6)]  # remove extra columns

# Read sample metadata
samples <- read.csv("samples.csv")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=samples,
                              design=~condition)

dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), "results/05_deseq2/deseq2_results.csv")
EOF
fi

# ============================================================
# Done
# ============================================================
echo "✅ RNA-seq Pipeline Completed Successfully!"
echo "Results saved in: $OUTDIR"
