setwd("C:/Users/aymer/Documents/Master/1st semester/467713 RNA-sequencing/RNAseq_breastcancer/DE_analysis")
source("scripts/utils.R")

# Load both tables
counts_simple <- read_featurecounts_table(COUNTS_FILE_SIMPLE)
counts_multi <- read_featurecounts_table(COUNTS_FILE_MULTI)


# Create DESeq2 datasets for both cases
dds_simple <- DESeqDataSetFromMatrix(countData = counts_simple,
                                     colData = generate_sample_metadata(counts_simple),
                                     design = ~ type)
dds_multi <- DESeqDataSetFromMatrix(countData = counts_multi,
                                    colData = generate_sample_metadata(counts_multi),
                                    design = ~ type)

# Run DESeq for both datasets
dds_simple <- DESeq(dds_simple)
dds_multi <- DESeq(dds_multi)

# Saving the dds objects
saveRDS(dds_simple, DDS_FILE_SIMPLE)
saveRDS(dds_multi, DDS_FILE_MULTI)