setwd("C:/Users/aymer/Documents/Master/1st semester/467713 RNA-sequencing/RNAseq_breastcancer/DE_analysis")
source("scripts/utils.R")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# Only load the "simple" DESeqDataSet object, as we will use it for the analysis
dds <- readRDS(DDS_FILE_SIMPLE)

group_interest <- c("TNBC", "NonTNBC")

# Perform differential expression analysis
res <- results(dds, contrast = c("type", group_interest))

# Annotating the results with gene names
res$geneName <- mapIds(org.Hs.eg.db,
                           keys = rownames(res),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "list")

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

# Defining genes of interest
genes_interest <- head(res[order(res$padj), "geneName"], 3)
genes_interest <- c("SPARC", "RACK1", "APOE")

# For each gene of interest, plot the expression levels across the different sample groups (TNBC, NonTNBC)
counts_plots <- lapply(genes_interest, generate_gene_counts_plot, dds = dds, res = res, group_interest = group_interest)

grid.arrange(grobs = counts_plots, ncol = length(counts_plots))
