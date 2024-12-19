setwd("C:/Users/aymer/Documents/Master/1st semester/467713 RNA-sequencing/RNAseq_breastcancer/DE_analysis")
source("scripts/utils.R")

library(org.Hs.eg.db) 
library(clusterProfiler)
library(DOSE)
library(enrichplot)

# Load the DE genes file
DE_genes <- read.csv(get_DE_genes_file(PAIR_INTEREST), row.names = 1)
# Filter significant DE genes
DE_genes_sig <- filter_DE_genes(DE_genes)

# Perform over representation analysis
ego <- enrichGO(gene          = rownames(DE_genes_sig),
                universe      = rownames(DE_genes),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP", # I think we are only interested in biological process
                keyType       = "ENSEMBL")

# Plot the top 20 GO terms
barplot(ego, showCategory=20) 

ego_df <- as.data.frame(ego)
# Get the top 20 GO terms
top_go_terms <- ego_df[order(ego_df$p.adjust), ][1:20, c("ID", "Description", "p.adjust")]
