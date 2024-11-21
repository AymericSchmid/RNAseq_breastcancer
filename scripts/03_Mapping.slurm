#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=15:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_mapping
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/mapping_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/mapping_error_%J.o

# Define variables
READS_LIST="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/metadata/reads_list.tsv"
INDEX_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Reference_genome"
INDEX_BASENAME="$INDEX_DIR/genome_index"
RESULTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Mapping"
READS_DIR="/data/courses/rnaseq_course/breastcancer_de/reads"
HISAT2_IMG="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p $RESULTS_DIR

# Get sample name and reads paths related to the running job in the reads list
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $READS_LIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $READS_LIST)
READ2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $READS_LIST)

echo "Mapping $SAMPLE">&2

# Mapping
apptainer exec --bind $READS_DIR,$RESULTS_DIR,$INDEX_DIR $HISAT2_IMG hisat2 -x $INDEX_BASENAME -1 $READ1 -2 $READ2 -S $RESULTS_DIR/$SAMPLE.sam -p $SLURM_CPUS_PER_TASK