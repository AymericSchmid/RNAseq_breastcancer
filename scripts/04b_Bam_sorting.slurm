#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25G
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_bam_sorting
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/bam_sorting_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/bam_sorting_error_%J.o

UNSORTED_BAM_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/BAM_processing"
SAMPLE_LIST="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/metadata/reads_list.tsv"
RESULTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/BAM_processing"         # Same as unsorted bam folder for now
HISAT2_IMG="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p $RESULTS_DIR

# Get sample name
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLE_LIST)

echo "Sorting bam file for $SAMPLE">&2

apptainer exec --bind $UNSORTED_BAM_DIR,$RESULTS_DIR $HISAT2_IMG samtools sort -m 25G -@ $SLURM_CPUS_PER_TASK -o $RESULTS_DIR/$SAMPLE.bam -T temp $RESULTS_DIR/$SAMPLE"_unsorted.bam"

# If sorted file exists, remove the unsorted one.
# The test avoid to delete the unsorted one if something went wrong with the sorting command
if [[ -f "$RESULTS_DIR/$SAMPLE.bam" ]]
then
    rm $RESULTS_DIR/$SAMPLE"_unsorted.bam"
else
    echo "Error: Something went wrong with the sorting..." >&2
fi