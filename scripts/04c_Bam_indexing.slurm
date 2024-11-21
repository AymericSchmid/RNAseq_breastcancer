#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_bam_indexing
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/bam_indexing_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/bam_indexing_error_%J.o

BAM_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/BAM_processing"
SAMPLE_LIST="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/metadata/reads_list.tsv"
HISAT2_IMG="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Get sample name
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLE_LIST)

echo "Indexing $SAMPLE.bam">&2

apptainer exec --bind $BAM_DIR $HISAT2_IMG samtools index $BAM_DIR/$SAMPLE.bam