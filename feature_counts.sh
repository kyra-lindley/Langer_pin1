#!/bin/bash
#SBATCH --job-name=featureCounts_bulk
#SBATCH --output=logs/featureCounts_bulk.out
#SBATCH --error=logs/featureCounts_bulk.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8



# Paths
BAM_DIR=/home/groups/CEDAR/archive/seq/Langer/bowman_bulk/aligned
OUT_DIR=/home/groups/CEDAR/archive/seq/Langer/bowman_bulk/counts
GTF=/home/groups/CEDAR/lindley/genome/GRCh38/Homo_sapiens.GRCh38.109.gtf

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Run featureCounts
featureCounts \
  -T 8 \
  -p \
  -t exon \
  -g gene_id \
  -a "$GTF" \
  -o "$OUT_DIR"/counts.txt \
  "$BAM_DIR"/*Aligned.sortedByCoord.out.bam
