setwd("C:/Users/aymer/Documents/Master/1st semester/467713 RNA-sequencing/RNAseq_breastcancer/DE_analysis")
source("scripts/utils.R")

library(org.Hs.eg.db)
library(pheatmap)

# Only load the "simple" DESeqDataSet object, as we will use it for the analysis
dds <- readRDS(DDS_FILE_SIMPLE)

# Perform differential expression analysis
res <- results(dds, contrast = c("type", PAIR_INTEREST))

# Annotating the results with gene names
res$geneName <- mapIds(org.Hs.eg.db,
                           keys = rownames(res),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
# Saving the significant genes to a file
write.csv(res, get_DE_genes_file(PAIR_INTEREST))

# Filter significant genes
res_sig <- filter_DE_genes(res)

# Number of differentially expressed genes
num_de_genes <- nrow(res_sig)
cat("Number of DE genes (padj < 0.05):", num_de_genes, "\n")

# Number of upregulated and downregulated genes
num_upregulated <- sum(res_sig$log2FoldChange > 0)
num_downregulated <- sum(res_sig$log2FoldChange < 0)
cat("Upregulated genes:", num_upregulated, "\n")
cat("Downregulated genes:", num_downregulated, "\n")

# Defining genes of interest
genes_interest_ids <- head(rownames(res[order(res$padj),]), 10) # Getting the top genes with the lowest p-value
genes_interest_ids <- c(genes_interest_ids, 
                    rownames(res[res$geneName %in% c("SPARC", "RACK1", "APOE"),]))             # Manually defining genes of interest

# For each gene of interest, plot the expression levels across the different sample groups (TNBC, NonTNBC)
counts_plots <- lapply(genes_interest_ids, generate_gene_counts_plot, dds = dds, res = res, pair_interest = PAIR_INTEREST)

grid.arrange(grobs = counts_plots)

# Subset the vst-transformed matrix
vsd_subset <- assay(vst(dds, blind = TRUE))[genes_interest_ids,]
rownames(vsd_subset) <- res[genes_interest_ids, "geneName"]

# Annotation for heatmap columns
pheatmap_annotation_col <- data.frame(type=colData(dds)[,"type"])
rownames(pheatmap_annotation_col) <- rownames(colData(dds))


# Generate heatmap
pheatmap(vsd_subset, 
         cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=FALSE,
         annotation_col=pheatmap_annotation_col,
         color = colorRampPalette(c("royalblue1", "ivory", "violetred"))(50))

# Do a volcano or heatmap (specific to some genes like ESR1 and ESR2 or top5 most expressed
#                           The heatmap take the counts as input)
# ESR1 and ESR2 are the ER genes
# ERBB2 is the HER2 gene

# Next step is do a dotplot ?
