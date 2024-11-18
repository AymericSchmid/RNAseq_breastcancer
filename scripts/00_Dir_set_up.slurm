#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:05:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_dir_set_up
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/dir_set_up_output_%j.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/dir_set_up_error_%j.o

# Define variables
RESULTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results"
TOOLS_DIR="/data/courses/rnaseq_course/tools"
READS_DIR="/data/courses/rnaseq_course/breastcancer_de/reads"
METADATA_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/metadata"

# Creating a directory for the outpout files of the worklfow
mkdir $RESULTS_DIR
mkdir $METADATA_DIR

# Generating a file that will contain a list of path to the fastq files
$TOOLS_DIR/get_samplelist.sh $READS_DIR > $METADATA_DIR/reads_list.tsv