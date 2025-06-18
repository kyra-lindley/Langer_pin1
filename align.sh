#!/bin/bash
#SBATCH --job-name=star_human
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err
#SBATCH --array=1-10        # Adjust to number of samples
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=12:00:00

module load STAR  # or use your STAR path if not via modules

# Paths
FASTQ_DIR=/home/groups/CEDAR/lindley/bulk/Bulk-RNA-seq-pipeline-PE/samples/trimmed
GENOME_DIR=/path/to/genome_index
OUT_DIR=/home/groups/CEDAR/archive/seq/Langer/bowman_bulk/aligned
SAMPLES_FILE=samples.txt

# Get sample name from file
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLES_FILE)

# Input files
R1=${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz
R2=${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz

# Create sample-specific output directory
mkdir -p ${OUT_DIR}/${SAMPLE}

# Run STAR
STAR \
  --genomeDir $GENOME_DIR \
  --readFilesIn $R1 $R2 \
  --readFilesCommand zcat \
  --runThreadN 8 \
  --outFileNamePrefix ${OUT_DIR}/${SAMPLE}/${SAMPLE}_ \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts
