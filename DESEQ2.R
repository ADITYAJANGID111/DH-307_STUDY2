# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

# Load count matrices
gene_counts <- read.csv("gene_count_matrix.csv", row.names = 1)
transcriptome_counts <- read.csv("transcript_count_matrix.csv", row.names = 1)

# Sample information (replace with your actual sample information)
col_data <- data.frame(
  SampleID = colnames(gene_counts),  # Assuming column names in gene_counts are the sample IDs
  cell_line = rep(c("U251", "U251", "U343", "U343"), length.out = 17),
  time = rep(c("T1", "T2", "T1", "T2"), length.out = 17)
)

# print the deimensions of matrix
print(dim(gene_counts))
print(dim(transcriptome_counts))
print(dim(col_data))

# Create DESeqDataSet objects
dds_gene <- DESeqDataSetFromMatrix(countData = gene_counts, colData = col_data, design = ~ cell_line + time)
dds_transcriptome <- DESeqDataSetFromMatrix(countData = transcriptome_counts, colData = col_data, design = ~ cell_line + time)

# Run DESeq analysis
dds_gene <- DESeq(dds_gene)
dds_transcriptome <- DESeq(dds_transcriptome)

# Perform time-specific contrasts
contrast_gene <- c("time", "T2", "T1")
contrast_transcriptome <- c("time", "T2", "T1")

dds_gene_results <- results(dds_gene, contrast = contrast_gene)
dds_transcriptome_results <- results(dds_transcriptome, contrast = contrast_transcriptome)

# Adjust p-values
dds_gene_results <- results(dds_gene, contrast = contrast_gene, alpha = 0.05)
dds_transcriptome_results <- results(dds_transcriptome, contrast = contrast_transcriptome, alpha = 0.05)

# Filter differentially expressed genes
DE_genes <- subset(dds_gene_results, padj < 0.05)
DE_transcripts <- subset(dds_transcriptome_results, padj < 0.05)

# Explore and visualize results (example: MA plot)
plotMA(dds_gene_results, main="DESeq2", ylim=c(-2,2))

# Save results
write.csv(DE_genes, "differential_genes.csv")
write.csv(DE_transcripts, "differential_transcripts.csv")