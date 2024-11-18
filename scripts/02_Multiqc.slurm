#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10M
#SBATCH --time=00:30:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_MultiQC
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/multiqc_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/multiqc_control_error_%J.o

# Define variables
QC_REPORTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Quality_control"
RESULTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Quality_control"        # Same as QC reports for now
MULTIQC_IMG="/containers/apptainer/multiqc-1.19.sif"

# Running MultiQC
apptainer run --bind $QC_REPORTS_DIR $MULTIQC_IMG multiqc $QC_REPORTS_DIR -o $RESULTS_DIR