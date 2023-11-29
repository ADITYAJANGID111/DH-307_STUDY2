# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

# Load count matrices
gene_counts <- read.csv("gene_count_matrix.csv", row.names = 1)

# Adding small count to address zero gene-counts issue
gene_counts <- gene_counts + 1

# Filter out low-count genes
gene_counts_filtered <- gene_counts[rowSums(gene_counts) > 1, ]

# Sample information (replace with your actual sample information)
col_data <- data.frame(
  SampleID = colnames(gene_counts),  # Assuming column names in gene_counts are the sample IDs
  cell_line = rep(c("U251", "U251", "U343", "U343"), length.out = 17),
  time = rep(c("T1", "T2", "T1", "T2"), length.out = 17)
)

# Print matirx dimensions
print(dim(gene_counts))
print(dim(col_data))

# Create DESeqDataSet objects
dds_gene <- DESeqDataSetFromMatrix(countData = gene_counts, colData = col_data, design = ~ cell_line + time)

# Run DESeq analysis
dds_gene <- DESeq(dds_gene)

# Perform time-specific contrasts
contrast_gene <- c("time", "T2", "T1")
contrast_transcriptome <- c("time", "T2", "T1")

dds_gene_results <- results(dds_gene, contrast = contrast_gene)

# Adjust p-values
dds_gene_results <- results(dds_gene, contrast = contrast_gene, alpha = 0.05)

# Filter differentially expressed genes
DE_genes <- subset(dds_gene_results, padj < 0.05)

# Save results
#write.csv(DE_genes, "differential_genes.csv")

# Top 10 genes arranged according to the ascending order of adjusted p-values(padj)
top10_genes <- head(dds_gene_results[order(dds_gene_results$padj), ], 10)
write.csv(top10_genes, "top10_genes.csv")
