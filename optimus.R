# Set relative directories for saving data and images

use_my_local_paths <- FALSE  # <- Dont change this please

## For ME :(

if (use_my_local_paths) {
  image_dir <- "D:/BIN_BCDR_HER2/Images"
  data_dir <- "D:/BIN_BCDR_HER2/Data/R_Data"
  ML_data_dir <- "D:/BIN_BCDR_HER2/Data/ML_Data"
} else {
  image_dir <- "Images"
  data_dir <- "Data/R_Data"
  ML_data_dir <- "Data/ML_Data"
}

## For You :)

# Create directories 
dir.create(image_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ML_data_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1: Data Inspection and Preprocessing

# Load libraries
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# Load dataset
gse <- getGEO("GSE65194", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
expression_data <- exprs(gse[[1]])

# Check for missing values
cat("Missing values in expression data:", sum(is.na(expression_data)), "\n")

# Boxplot before and after normalization
png(file = file.path(image_dir, "boxplot_before_normalization.png"))
boxplot(expression_data, main = "Before Normalization", las = 2, outline = FALSE)
dev.off()

expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")

png(file = file.path(image_dir, "boxplot_after_normalization.png"))
boxplot(expression_data, main = "After Normalization", las = 2, outline = FALSE)
dev.off()

# Filter low-expressed genes
gene_means <- rowMeans(expression_data)
threshold <- quantile(gene_means, probs = 0.25)
expression_data <- expression_data[gene_means > threshold, ]
cat("Dimensions after filtering:", dim(expression_data), "\n")

# Extract subtype
subtype <- metadata$`sample_group:ch1`
subtype <- as.factor(subtype)
cat("Subtype distribution:\n")
table(subtype)

# PCA plot
pca_result <- prcomp(t(expression_data), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Subtype = subtype)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point() +
  ggtitle("PCA Plot: Subtype Clustering") +
  theme_minimal()

ggsave(file.path(image_dir, "pca_plot.png"), plot = pca_plot)

# Save processed data
save(expression_data, metadata, subtype, file = file.path(data_dir, "processed_data_step1.RData"))

# Step 2: Exploratory Data Analysis (EDA)

# Heatmap of top 100 most variable genes
library(pheatmap)
top_var_genes <- order(apply(expression_data, 1, var), decreasing = TRUE)[1:50]
heatmap_data <- expression_data[top_var_genes, ]
annotation_col <- data.frame(Subtype = subtype)
rownames(annotation_col) <- colnames(expression_data)

pheatmap_plot <- pheatmap(heatmap_data,
                          annotation_col = annotation_col,
                          show_rownames = FALSE,
                          show_colnames = FALSE,
                          main = "Heatmap of Top 50 Most Variable Genes")

ggsave(file.path(image_dir, "heatmap_top_50_genes.png"))

# Hierarchical clustering
dist_matrix <- dist(t(expression_data))
hclust_result <- hclust(dist_matrix, method = "complete")
png(file = file.path(image_dir, "hierarchical_clustering_dendrogram.png"))
plot(hclust_result, labels = subtype, main = "Hierarchical Clustering Dendrogram")
dev.off()

# Step 3: Differential Expression Analysis (DEA)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Prepare DESeq2 data
dds <- DESeqDataSetFromMatrix(countData = round(expression_data), 
                              colData = data.frame(condition = subtype), 
                              design = ~ condition)
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "Her2", "TNBC"))

# Filter significant DEGs
significant_genes <- results[!is.na(results$padj) & results$padj < 0.05, ]
cat("Number of significant DEGs:", nrow(significant_genes), "\n")
write.csv(as.data.frame(significant_genes), file = file.path(data_dir, "significant_genes_her2_vs_tnbc.csv"))

# Volcano plot
library(ggplot2)
library(ggrepel)

volcano_df <- as.data.frame(results)
volcano_df$gene <- rownames(volcano_df)
volcano_df$significant <- ifelse(volcano_df$padj < 0.05, "Significant", "Not significant")

# Volcano plot with top 10 DEGs labeled
top_genes <- volcano_df %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(10)

volcano_with_labels_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_label_repel(data = top_genes, aes(label = gene), box.padding = 0.5, max.overlaps = Inf) +
  labs(title = "Volcano Plot with Top 10 DEGs",
       x = "log2 Fold Change (Her2/TNBC)",
       y = "-log10(Adjusted p-value)") +
  theme_minimal()

ggsave(file.path(image_dir, "volcano_with_labels_plot.png"), plot = volcano_with_labels_plot)

# MA Plot (Log2 fold change vs mean expression)
ma_plot <- ggplot(volcano_df, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(alpha = 0.7, color = "gray") +
  geom_point(data = subset(volcano_df, padj < 0.05), aes(color = significant), size = 2) +
  scale_color_manual(values = c("red")) +
  labs(title = "MA Plot", x = "Log2 Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

ggsave(file.path(image_dir, "ma_plot.png"), plot = ma_plot)

# Step 4: Pathway Enrichment Analysis

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

# Map probes to Entrez IDs
gpl <- getGEO("GPL570", destdir = ".")
probe_to_entrez <- Table(gpl)[, c("ID", "ENTREZ_GENE_ID")]
probe_to_entrez <- probe_to_entrez[probe_to_entrez$ENTREZ_GENE_ID != "", ]
significant_genes_df <- as.data.frame(significant_genes)
significant_genes_df$probe_id <- rownames(significant_genes_df)
significant_genes_df <- merge(significant_genes_df, probe_to_entrez, by.x = "probe_id", by.y = "ID")
entrez_ids <- significant_genes_df$ENTREZ_GENE_ID

# GO enrichment analysis
go_enrichment <- enrichGO(gene = entrez_ids,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
summary(go_enrichment)
dotplot(go_enrichment, showCategory = 8, title = "GO Enrichment Analysis (Biological Process)") +
  theme(axis.text.y = element_text(size = 10))

# Save GO results
write.csv(as.data.frame(go_enrichment), file = file.path(data_dir, "go_enrichment_results.csv"))


##################  WILL FIX THIS LATER - NOT WORKING ENTRZ IDS NOT LABELLING #################

# KEGG enrichment analysis
# kegg_mapping <- bitr(entrez_ids, 
#                      fromType = "ENTREZID", 
#                      toType = "PATH", 
#                      OrgDb = org.Hs.eg.db)
# valid_kegg_ids <- kegg_mapping$PATH[!is.na(kegg_mapping$PATH)]
# cat("Number of valid KEGG IDs:", length(valid_kegg_ids), "\n")
# if (length(valid_kegg_ids) == 0) {
#   stop("No valid KEGG IDs found.")
# }
# kegg_enrichment <- enrichKEGG(
#   gene = valid_kegg_ids,
#   organism = "hsa",
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.5,
#   qvalueCutoff = 0.5,
# )
# if (nrow(as.data.frame(kegg_enrichment)) > 0) {
#   print(summary(kegg_enrichment))
#   dotplot(kegg_enrichment, 
#           showCategory = 15, 
#           title = "KEGG Pathway Enrichment Analysis") +
#     theme(axis.text.y = element_text(size = 10))
#   write.csv(as.data.frame(kegg_enrichment), file = "kegg_enrichment_results.csv")
# } else {
#   cat("No significant KEGG pathways found.\n")
# }

###################### END OF BUG (FIX IT LATER) ####################################


# Step 5: Machine Learning Data Export

# Save expression data (transposed for Python: samples x genes)
her2_tnbc_samples <- subtype %in% c("Her2", "TNBC")
expression_subset <- t(expression_data[, her2_tnbc_samples])  # Transpose for Python
write.csv(expression_subset, file.path(ML_data_dir, "expression_data_her2_tnbc.csv"))

# Save subtype labels
subtype_subset <- subtype[her2_tnbc_samples]

##############  ERROR IN DATA FRAM AS DIFFERING NUMBER OF ROWS: 41004,94  - WILL FIX LATER ################


#write.csv(data.frame(Sample = colnames(expression_subset), Subtype = subtype_subset), 
#          file.path(ML_data_dir, "subtype_labels.csv"), row.names = FALSE)


#################################################################################################

# Save significant DEGs

significant_genes <- rownames(significant_genes)  # From DESeq2 results
write.csv(data.frame(Gene = significant_genes), file.path(ML_data_dir, "significant_genes.csv"), row.names = FALSE)