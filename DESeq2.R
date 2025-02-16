# Load necessary packages
required_packages <- c("DESeq2", "ggplot2", "pheatmap", "BiocManager")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg == "BiocManager") {
      install.packages("BiocManager")
    } else {
      BiocManager::install(pkg)
    }
  }
}

# Load installed packages
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Set working directory
setwd("/Users/shaonihui/desktop")

# Read count matrix
counts_data <- read.table("reformatted_counts.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Extract sample names from file paths
colnames(counts_data) <- gsub("_mappedReads_sorted\\.bam", "", basename(colnames(counts_data)))
colnames(counts_data) <- gsub("_counts\\.txt", "", colnames(counts_data)) # If necessary

# Define sample information table
sample_info <- data.frame(
  sample = colnames(counts_data),
  condition = c(rep("HER2", 3), rep("NonTNBC", 3), rep("Normal", 3), rep("TNBC", 3)),
  replicate = rep(1:3, 4)  # 3 replicates per condition
)

rownames(sample_info) <- sample_info$sample
print(sample_info)

# Ensure column names match sample info row names
if (!all(colnames(counts_data) %in% rownames(sample_info))) {
  stop("Mismatch between count matrix columns and sample info rows!")
}

# Reorder count matrix columns to match sample info
counts_data <- counts_data[, rownames(sample_info)]
if (!all(colnames(counts_data) == rownames(sample_info))) {
  stop("Sample names in count matrix and sample info are not in the same order!")
}

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = sample_info,
  design = ~ condition
)

# Filter low-expression genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract differential expression results
res <- results(dds)
res <- lfcShrink(dds, coef = 2, type = "apeglm")

# Save results
write.csv(as.data.frame(res), "differential_expression_results.csv")

# **Visualization**

# Volcano plot
res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "yes", "no")
volcano_plot <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value")
ggsave("Volcano_plot.jpeg", plot = volcano_plot, width = 6, height = 4)

# PCA plot
vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()
ggsave("PCA_plot_condition.jpeg", plot = pca_plot, width = 6, height = 4)

# Heatmap (Top 20 most variable genes)
top_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
pheatmap(assay(vsd)[top_genes, ], cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, annotation_col = as.data.frame(colData(dds)["condition"]),
         filename = "Heatmap_top20_genes.jpeg")

# Save normalized counts
normalized_counts <- as.data.frame(assay(vsd))
write.csv(normalized_counts, "normalized_counts_matrix.csv")

# Sample-to-sample distance heatmap
sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)
pheatmap(sample_dist_matrix, clustering_distance_rows = sample_dist, clustering_distance_cols = sample_dist,
         annotation_col = as.data.frame(colData(dds)["condition"]), main = "Sample-to-Sample Distance Heatmap",
         filename = "Sample_to_Sample_Heatmap.jpeg")