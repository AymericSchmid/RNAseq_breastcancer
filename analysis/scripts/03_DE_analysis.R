setwd("C:/Users/aymer/Documents/Master/1st semester/467713 RNA-sequencing/Project/Analysis")
source("scripts/utils.R")

# Only load the "simple" DESeqDataSet object, as we will use it for the analysis
dds <- readRDS(DDS_FILE_SIMPLE)

# Perform differential expression analysis
res <- results(dds, contrast = c("type", "TNBC", "NonTNBC"))

# Filter significant genes
res_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]

# Number of differentially expressed genes
num_de_genes <- nrow(res_sig)
cat("Number of DE genes (padj < 0.05):", num_de_genes, "\n")

# Number of upregulated and downregulated genes
num_upregulated <- sum(res_sig$log2FoldChange > 0)
num_downregulated <- sum(res_sig$log2FoldChange < 0)
cat("Upregulated genes:", num_upregulated, "\n")
cat("Downregulated genes:", num_downregulated, "\n")

# TODO: For each gene of interest
plotCounts(dds, gene=which.min(res$padj), intgroup="type")


sum(counts(dds)[,"HER21"]) + sum(counts(dds)[,"HER22"]) + sum(counts(dds)[,"HER23"])
sum(counts(dds)[,"TNBC1"]) + sum(counts(dds)[,"TNBC2"]) + sum(counts(dds)[,"TNBC3"])
sum(counts(dds)[,"NonTNBC1"]) + sum(counts(dds)[,"NonTNBC2"]) + sum(counts(dds)[,"NonTNBC3"])
sum(counts(dds)[,"Normal1"]) + sum(counts(dds)[,"Normal2"]) + sum(counts(dds)[,"Normal3"])


# TODO: Based on the original publication, select 2-3 genes that are of particular interest and investigate their expression level. You could use, for example, the normalised counts (see DESeq2::counts) where the effect of between-sample differences in sequencing depth has been removed.