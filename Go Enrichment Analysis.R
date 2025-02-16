# Set working directory
setwd("/Users/shaonihui/desktop")

# Ensure Bioconductor and necessary R packages are installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}

# Automatically install required Bioconductor packages
required_packages <- c("clusterProfiler", "org.Hs.eg.db", "ggplot2", "enrichplot", "cowplot", "patchwork")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Load installed R packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(cowplot)
library(patchwork)

# Function to perform GO enrichment analysis
run_go_enrichment <- function(de_genes, all_genes, comparison_name, gene_list) {
  if (length(de_genes) == 0) {
    warning(paste0("No significant genes found for ", comparison_name))
    return(NULL)
  }
  
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(
    gene = de_genes,
    universe = all_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",  # Biological Process
    keyType = "ENSEMBL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # Save GO enrichment results
  write.csv(go_enrichment@result, file = paste0("GO_Enrichment_", comparison_name, ".csv"))
  
  # Generate Dot Plot
  dot_p <- dotplot(go_enrichment, showCategory = 30, font.size = 14) + ggtitle(paste0("GO Enrichment: ", comparison_name))
  print(dot_p)  # Display in RStudio Plots window
  ggsave(paste0("GO_Dotplot_", comparison_name, ".jpeg"), width = 16, height = 12, dpi = 300)
  
  # Generate Bar Plot
  bar_p <- barplot(go_enrichment, showCategory = 20, title = paste0("GO Enrichment: ", comparison_name))
  print(bar_p)
  ggsave(paste0("GO_Barplot_", comparison_name, ".jpeg"), width = 16, height = 12, dpi = 300)
  
  # Generate Gene-Concept Network Plot
  edox <- setReadable(go_enrichment, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")
  cnet_p <- cnetplot(edox, foldChange = gene_list)
  print(cnet_p)
  ggsave(paste0("GO_Cnetplot_", comparison_name, ".jpeg"), width = 12, height = 12, dpi = 300)
  
  # Generate Heatmap
  heat_p <- heatplot(edox, showCategory = 10, foldChange = gene_list)
  print(heat_p)
  ggsave(paste0("GO_Heatplot_", comparison_name, ".jpeg"), width = 16, height = 10, dpi = 300)
  
  # Generate Tree Plot
  edox2 <- pairwise_termsim(edox)
  tree_p <- treeplot(edox2)
  print(tree_p)
  ggsave(paste0("GO_Treeplot_", comparison_name, ".jpeg"), width = 12, height = 12, dpi = 300)
  
  print(paste0("Top 10 GO terms for ", comparison_name, ":"))
  print(head(go_enrichment@result, 10))
  
  return(go_enrichment)
}

# Define background gene set
background_genes <- rownames(dds)

# TNBC vs Normal
res_tnbc_vs_normal <- results(dds, contrast = c("condition", "TNBC", "Normal"))
de_genes_tnbc <- rownames(subset(res_tnbc_vs_normal, padj < 0.05 & abs(log2FoldChange) > 1))
if (length(de_genes_tnbc) > 0) {
  gene_list_tnbc <- res_tnbc_vs_normal$log2FoldChange
  names(gene_list_tnbc) <- rownames(res_tnbc_vs_normal)
  go_tnbc <- run_go_enrichment(de_genes_tnbc, background_genes, "TNBC_vs_Normal", gene_list_tnbc)
}

# HER2 vs Normal
res_her2_vs_normal <- results(dds, contrast = c("condition", "HER2", "Normal"))
de_genes_her2 <- rownames(subset(res_her2_vs_normal, padj < 0.05 & abs(log2FoldChange) > 1))
if (length(de_genes_her2) > 0) {
  gene_list_her2 <- res_her2_vs_normal$log2FoldChange
  names(gene_list_her2) <- rownames(res_her2_vs_normal)
  go_her2 <- run_go_enrichment(de_genes_her2, background_genes, "HER2_vs_Normal", gene_list_her2)
}

# NonTNBC vs Normal
res_nontnbc_vs_normal <- results(dds, contrast = c("condition", "NonTNBC", "Normal"))
de_genes_nontnbc <- rownames(subset(res_nontnbc_vs_normal, padj < 0.05 & abs(log2FoldChange) > 1))
if (length(de_genes_nontnbc) > 0) {
  gene_list_nontnbc <- res_nontnbc_vs_normal$log2FoldChange
  names(gene_list_nontnbc) <- rownames(res_nontnbc_vs_normal)
  go_nontnbc <- run_go_enrichment(de_genes_nontnbc, background_genes, "NonTNBC_vs_Normal", gene_list_nontnbc)
}
