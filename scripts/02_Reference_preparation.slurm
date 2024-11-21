#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=05:00:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=RNAseq_reference_preparation
#SBATCH --output=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/reference_preparation_output_%J.o
#SBATCH --error=/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/output/reference_preparation_error_%J.o

# Define variables
REFERENCE_DIR="/data/users/aschmid/RNAseq_project/RNAseq_breastcancer/results/Reference_genome"     
REFERENCE_GENOME_LINK="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REFERENCE_ANNOTATION_LINK="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz"
GENOME_CHECKSUM_LINK="https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/CHECKSUMS"
ANNOTATION_CHECKSUM_LINK="https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/CHECKSUMS"
HISAT2_IMG="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

reference_genome=$(basename $REFERENCE_GENOME_LINK)
reference_annotation=$(basename $REFERENCE_ANNOTATION_LINK)

mkdir -p $REFERENCE_DIR
cd $REFERENCE_DIR

# Downloading reference genome and genome annotation
wget $REFERENCE_GENOME_LINK
wget $REFERENCE_ANNOTATION_LINK

# Check file integrity

# download checksum file without saving it -> get the line corresponding to the files we use -> keep the 2 first columns
genome_expected_checksum=$(wget -O - $GENOME_CHECKSUM_LINK | grep $reference_genome | awk '{print $1,$2}')
annotation_expected_checksum=$(wget -O - $ANNOTATION_CHECKSUM_LINK | grep $reference_annotation | awk '{print $1,$2}')

# Testing if checksums calculated on downloaded files are as expected
if [ "$genome_expected_checksum" != "$(sum $reference_genome)" ] || \
    [ "$annotation_expected_checksum" != "$(sum $reference_annotation)" ]
then 
    # throw an error and exit
    echo "Checksums are not as expected !">&2
    rm $reference_genome
    rm $reference_annotation
    exit 1
fi

# Unzipping the reference genome if needed
if [[ "$reference_genome" == *.gz ]]
then
    echo "Unzipping the reference genome..."
    gunzip $reference_genome
    reference_genome=${reference_genome%.*}         # After unzipping, remove the .gz extension from the file name
fi

# Indexing the reference genome
# the last parameter is the prefix of output files
apptainer run --bind $REFERENCE_DIR $HISAT2_IMG hisat2-build -p $SLURM_CPUS_PER_TASK $REFERENCE_DIR/$reference_genome $REFERENCE_DIR/genome_index