#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=15:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_reads_counting
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/reads_counting_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/reads_counting_error_%J.o

BAM_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/BAM_processing"
SAMPLE_LIST="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/metadata/reads_list.tsv"
ANNOTATION="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Reference_genome/Homo_sapiens.GRCh38.113.gtf.gz"
RESULTS_FILE="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Quantification/gene_counts_table.txt"
FEATURECOUNTS_IMG="/containers/apptainer/subread_2.0.1--hed695b0_0.sif"

mkdir -p $(dirname $RESULTS_FILE)

# Get the list of BAM files by concatenation of the BAM directory + / + the sample name + .bam
BAM_FILES=$(awk -v dir="$BAM_DIR" '{print dir "/" $1 ".bam"}' $SAMPLE_LIST)

# Calling featureCounts withe the list of all related bam files. This will create a gene_counts_table.txt file that with all the samples in it
apptainer exec --bind $BAM_DIR,$(dirname $RESULTS_FILE),$(dirname $ANNOTATION) $FEATURECOUNTS_IMG featureCounts -T $SLURM_CPUS_PER_TASK -a $ANNOTATION -o $RESULTS_FILE $BAM_FILES