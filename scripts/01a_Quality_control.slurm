#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_quality_control
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/quality_control_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/quality_control_error_%J.o

# Define variables
READS_LIST="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/metadata/reads_list.tsv"
RESULTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Quality_control"
READS_DIR="/data/courses/rnaseq_course/breastcancer_de/reads"
FASTQC_IMG="/containers/apptainer/fastqc-0.12.1.sif"

mkdir -p $RESULTS_DIR

# Get sample name and reads paths related to the running job in the reads list
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $READS_LIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $READS_LIST)
READ2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $READS_LIST)

# Running fastqc from container
# Each job runs sample's R1 and R2
apptainer run --bind $READS_DIR $FASTQC_IMG fastqc -t $SLURM_CPUS_PER_TASK -o $RESULTS_DIR $READ1 $READ2