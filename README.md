# **Breast Cancer Transcriptomic Analysis**

## **Project Overview**
This repository contains the code and instructions for analyzing transcriptomic differences between TNBC, Non-TNBC and HER2 cancer subtypes using RNA-seq data. The study focuses on differential expression and gene ontology enrichment analyses to identify key molecular differences between these subtypes.

## **Dataset Information**
The dataset is a subset of the **Eswaran et al. (2012)** study, containing RNA-seq data from **TNBC, Non-TNBC, HER2, and Normal** breast tissue samples. The RNA sequencing was performed using **paired-end sequencing** on an **Illumina HiSeq 2000** platform.

- **12 total samples:**
  - **TNBC (n=3)**
  - **Non-TNBC (n=3)**
  - **HER2 (n=3)**
  - **Normal (n=3)**
- **Read lengths:**
  - **Cancer samples:** 60bp paired-end reads
  - **Normal samples:** 75bp paired-end reads (derived from organoid tissues)
- **Reference Genome:** *Homo sapiens* (GRCh38, Ensembl Release 113)
- **Raw Data Access:** Available in the Gene Expression Omnibus (GEO) under **accession GSE52194**.

---

## **Repository Structure**
```
├── DE_analysis/
│   ├── scripts/                 # R scripts for data analysis
│   │   ├── 01_preprocessing.R   # Preprocessing
│   │   ├── 02_PCA_Analysis.R    # Principal Component Analysis
│   │   ├── 03_DE_analysis.R     # Differential Expression analysis with DESeq2
│   │   ├── 04_GO_analysis.R     # Gene Ontology  enrichment analysis
│   │   ├── utils.R              # Utility functions and constants
├── metadata/
│   ├── reads_list.tsv           # List of sequencing reads and file locations
├── scripts/                     # SLURM scripts for execution on cluster
│   ├── 01a_Quality_control.slurm
|   ├── 01b_Multiqc.slurm
│   ├── 02_Reference_preparation.slurm
│   ├── 03_Mapping.slurm
│   ├── 04a_Sam_to_bam.slurm
│   ├── 04b_Bam_sorting.slurm
│   ├── 04c_Bam_indexing.slurm
│   ├── 05_Reads_counting.slurm
```

---

## **Pipeline Execution**
### **On a HPC cluster**
The pipeline is designed to be executed on an **HPC cluster** (using **SLURM** and **Apptainer containers**) in the following steps:
- Run `01a_Quality_control.slurm` and `01b_Multiqc.slurm`
- Prepare the reference genome with `02_Reference_preparation.slurm`
- Align reads using `03_Mapping.slurm`
- Process BAM files:
  - `04a_Sam_to_bam.slurm`
  - `04b_Bam_sorting.slurm`
  - `04c_Bam_indexing.slurm`
- Generate raw count matrices with `05_Reads_counting.slurm`

### **Locally**
- Preprocess quantification data with `01_preprocessing.R`.
- Perform PCA: `02_PCA_Analysis.R`
- Run differential expression analysis: `03_DE_analysis.R`
- Conduct GO enrichment analysis: `04_GO_analysis.R`

---

## **Software & Dependencies**
This pipeline requires the following tools:

### **System Requirements**
- **HPC cluster** with **SLURM job scheduler**
- **Apptainer**

### **RNA-seq Processing**
- FastQC (v0.12.1)
- MultiQC (v1.19)
- HISAT2 (v2.2.1)
- SAMtools (v1.20)
- featureCounts (Subread v2.0.1)

### **Differential Expression & GO Analysis (R v4.4.1)**
- DESeq2 (v1.46.0, Bioconductor v3.20.0)
- org.Hs.eg.db (v3.20.0)
- pheatmap (v1.0.12)
- clusterProfiler (v4.14.4)
- enrichplot (v1.26.3)
