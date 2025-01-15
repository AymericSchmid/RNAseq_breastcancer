#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_sam_to_bam
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/sam_to_bam_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/sam_to_bam_error_%J.o

SAM_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Mapping"
SAMPLE_LIST="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/metadata/reads_list.tsv"
RESULTS_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/BAM_processing"
HISAT2_IMG="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p $RESULTS_DIR

# Get sample name
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLE_LIST)

echo "Converting SAM to BAM for $SAMPLE">&2

apptainer exec --bind $SAM_DIR,$RESULTS_DIR $HISAT2_IMG samtools view -@ $SLURM_CPUS_PER_TASK -hbS $SAM_DIR/$SAMPLE.sam > $RESULTS_DIR/$SAMPLE"_unsorted.bam"