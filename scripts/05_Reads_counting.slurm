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
RESULTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Quantification"
FEATURECOUNTS_IMG="/containers/apptainer/subread_2.0.1--hed695b0_0.sif"

mkdir -p $RESULTS_DIR

# Get the list of BAM files by concatenation of the BAM directory + / + the sample name + .bam
BAM_FILES=$(awk -v dir="$BAM_DIR" '{print dir "/" $1 ".bam"}' $SAMPLE_LIST)

# Calling featureCounts withe the list of all related bam files. This will create a gene_counts_table.txt file with all the samples in it
# Using -p to specify that input data contain paired-end reads.
# Using -Q to set a mapping quality threshold
# If the parameter of the script is "multi" => multi-mapped reads will be counted
if [ "$1" == "multi" ]
then
    echo "Taking multi-mapped reads in account">&2
    results_file=$RESULTS_DIR/gene_counts_multi_table.txt
    # Using -M each alignment reported for a multi-mapping read will be counted. Each alignment will carry a fractional count (--fraction)
    apptainer exec --bind $BAM_DIR,$RESULTS_DIR,$(dirname $ANNOTATION) $FEATURECOUNTS_IMG featureCounts -M --fraction -p -Q 10 -T $SLURM_CPUS_PER_TASK -a $ANNOTATION -o $results_file $BAM_FILES
else
    results_file=$RESULTS_DIR/gene_counts_table.txt
    apptainer exec --bind $BAM_DIR,$RESULTS_DIR,$(dirname $ANNOTATION) $FEATURECOUNTS_IMG featureCounts -p -Q 10 -T $SLURM_CPUS_PER_TASK -a $ANNOTATION -o $results_file $BAM_FILES
fi