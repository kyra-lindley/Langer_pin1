#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --output=results/trimming/trimmomatic.out
#SBATCH --error=results/trimming/trimmomatic.err



# Set paths
ADAPTERS="/home/groups/CEDAR/lindley/bulk/Bulk-RNA-seq-pipeline-PE/data/adapters.fa"
RAW_DIR="/home/groups/CEDAR/lindley/bulk/Bulk-RNA-seq-pipeline-PE/samples/raw"
TRIMMED_DIR="samples/trimmed"

mkdir -p $TRIMMED_DIR

for fq1 in $RAW_DIR/*_R1.fastq.gz; do
    fq2="${fq1/_R1/_R2}"
    sample=$(basename $fq1 | cut -d_ -f1)

    java -jar $TRIMMOMATIC PE \
        -threads 4 \
        $fq1 $fq2 \
        $TRIMMED_DIR/${sample}_R1_paired.fq.gz $TRIMMED_DIR/${sample}_R1_unpaired.fq.gz \
        $TRIMMED_DIR/${sample}_R2_paired.fq.gz $TRIMMED_DIR/${sample}_R2_unpaired.fq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
