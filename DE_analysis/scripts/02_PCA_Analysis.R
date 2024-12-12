setwd("C:/Users/aymer/Documents/Master/1st semester/467713 RNA-sequencing/RNAseq_breastcancer/DE_analysis")
source("scripts/utils.R")

dds_simple <- readRDS(DDS_FILE_SIMPLE)
dds_multi <- readRDS(DDS_FILE_MULTI)

# Variance stabilizing transformation (VST) for both datasets
vsd_simple <- vst(dds_simple, blind = TRUE)
vsd_multi <- vst(dds_multi, blind = TRUE)

# PCA plots for both cases
# Extract PCA data for both datasets
pca_simple <- plotPCA(vsd_simple, intgroup = "type", returnData = TRUE)
pca_multi <- plotPCA(vsd_multi, intgroup = "type", returnData = TRUE)

# Create ggplot objects for both PCA plots
plot_simple <- custom_pca_plot(pca_simple, "PCA - Without Multimapped Reads")
plot_multi <- custom_pca_plot(pca_multi, "PCA - With Multimapped Reads")

# Arrange the plots side by side
grid.arrange(plot_simple, plot_multi, ncol = 2)

# Compare mapping metrics
assigned_simple <- sum(counts(dds_simple))
assigned_multi <- sum(counts(dds_multi))

cat("Number of genes assigned in samples (without multimapped reads):", assigned_simple)
cat("Number of genes assigned in samples (with multimapped reads):", assigned_multi)
