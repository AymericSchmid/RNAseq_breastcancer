# ------ Libraries ------

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2") 
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra") 
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("DOSE", quietly = TRUE)) BiocManager::install("DOSE")
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

# Load libraries used in multiple scripts
library(stringr)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggrepel)

# ------ Constants ------

# FeatureCounts table without multimapped reads
COUNTS_FILE_SIMPLE <- "data/gene_counts_table.txt"
# FeatureCounts table with multimapped reads
COUNTS_FILE_MULTI <- "data/gene_counts_multi_table.txt"
# DDS object without multimapped reads file
DDS_FILE_SIMPLE <- "results/dds_simple.rds"
# DDS object with multimapped reads file
DDS_FILE_MULTI <- "results/dds_multi.rds"
# Cancer types of interest, used in the DE analysis
PAIR_INTEREST <- c("TNBC", "NonTNBC")

# ------ Functions ------

read_featurecounts_table <- function(path){
  # Loading quantification table
  counts <- read.table(path, header = T)
  # Setting the row name with gene id
  rownames(counts) <- counts$Geneid
  # Select only columns that are sample-specific (.bam files)
  counts <- counts[, grepl("\\.bam$", colnames(counts))]
  # Extract sample names from column names (removing path and .bam suffix)
  names(counts) <- str_extract(names(counts), "[^\\.]+(?=\\.bam$)")
  # Rounding counts to integers (in case of fractional counts for multimapped reads)
  counts <- round(counts)
  
  return(counts)
}

generate_sample_metadata <- function(counts){
  # Extract sample types by removing trailing digits
  sample_types <- sub("\\d$", "", names(counts))
  
  return(data.frame(
    row.names = names(counts),
    type = factor(sample_types)
  ))
}

custom_pca_plot <- function(pca_data, title){
  percents <- round(100 * attr(pca_data, "percentVar"))
  plot <- ggplot(pca_data, aes(PC1, PC2, color=type, label=name)) +
    geom_point(size=3) +
    labs(title = title, color = "Sample Type") +
    xlab(paste0("PC1: ",percents[1],"% variance")) +
    ylab(paste0("PC2: ",percents[2],"% variance")) + 
    coord_fixed() +
    geom_label_repel(aes(label = name),
                     size=3,
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     show.legend = F) +
    theme_minimal()
  return(plot)
}

generate_gene_counts_plot <- function(gene_id, dds, res, pair_interest=c()) {
  # Get the gene name
  gene_name <- res[gene_id,"geneName"]
  # Extract the counts for the gene of interest
  plot_data <- plotCounts(dds, gene=gene_id, intgroup="type", returnData=T, normalized = T)
  # Filter the data for the groups of interest
  if(length(pair_interest) > 0)
    plot_data <- plot_data[plot_data$type %in% pair_interest,]
  
  # Create the ggplot object
  return(ggplot(plot_data, aes(x=type, y=count)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    #scale_y_log10(breaks = scales::breaks_log(10)(range(plot_data$count))) + # Dynamic breaks
    labs(title = gene_name))
}

# Function to get the file name where the DE genes will be saved for a given pair of interest
get_DE_genes_file <- function(pair_interest){
  return(paste0("results/DE_genes_", pair_interest[1], "_vs_", pair_interest[2], ".txt"))
}

get_interesting_DE_genes_file <- function(pair_interest){
  return(paste0("results/DE_genes_interesting_", pair_interest[1], "_vs_", pair_interest[2], ".txt"))
}

# Filter significant DE genes based on adjusted p-value threshold
filter_DE_genes <- function(res, padj_threshold = 0.05){
  return(res[!is.na(res$padj) & res$padj < padj_threshold,])
}