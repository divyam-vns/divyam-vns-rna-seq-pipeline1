#!/usr/bin/env bash
# RNA-seq pipeline: FASTQ -> counts -> DESeq2
# Usage: bash main.sh

set -e
THREADS=8
GENOME=config/genomes/GRCh38.fa
GTF=config/genomes/gencode.v44.gtf

mkdir -p results/qc results/star results/counts

echo "=== Step 1: Quality Control ==="
fastqc -t $THREADS data/*.fastq.gz -o results/qc
multiqc results/qc -o results/qc

echo "=== Step 2: Trimming ==="
for s in data/*_R1.fastq.gz; do
  base=$(basename $s _R1.fastq.gz)
  fastp -i data/${base}_R1.fastq.gz -I data/${base}_R2.fastq.gz \
    -o results/${base}_R1.trim.fastq.gz -O results/${base}_R2.trim.fastq.gz \
    -h results/${base}.fastp.html -j results/${base}.fastp.json -w $THREADS
done

echo "=== Step 3: STAR alignment ==="
mkdir -p star_index
STAR --runThreadN $THREADS --runMode genomeGenerate \
  --genomeDir star_index --genomeFastaFiles $GENOME --sjdbGTFfile $GTF --sjdbOverhang 100

for s in results/*_R1.trim.fastq.gz; do
  base=$(basename $s _R1.trim.fastq.gz)
  STAR --runThreadN $THREADS --genomeDir star_index \
    --readFilesIn results/${base}_R1.trim.fastq.gz results/${base}_R2.trim.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix results/star/${base}. \
    --outSAMtype BAM SortedByCoordinate --twopassMode Basic
  samtools index results/star/${base}.Aligned.sortedByCoord.out.bam
done

echo "=== Step 4: Counting (featureCounts) ==="
featureCounts -T $THREADS -a $GTF -o results/counts.txt results/star/*.bam

echo "=== Pipeline complete! ==="
echo "Check results/counts.txt and run DESeq2 in R."
