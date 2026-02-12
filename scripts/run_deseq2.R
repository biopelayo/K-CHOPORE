#!/usr/bin/env Rscript
# =============================================================
# K-CHOPORE DESeq2 Differential Expression Analysis
# =============================================================
# Performs differential gene/isoform expression analysis from
# FLAIR counts matrix output.
#
# Usage:
#   Rscript run_deseq2.R <counts_matrix> <output_dir> <padj_threshold> <lfc_threshold>
#
# Input:  FLAIR counts_matrix.tsv (genes/isoforms x samples)
# Output: DESeq2 results CSV, MA plot, volcano plot, PCA plot
# =============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript run_deseq2.R <counts_matrix> <output_dir> <padj_threshold> <lfc_threshold>")
}

counts_file  <- args[1]
output_dir   <- args[2]
padj_cutoff  <- as.numeric(args[3])
lfc_cutoff   <- as.numeric(args[4])

cat("[K-CHOPORE] DESeq2 Analysis\n")
cat(sprintf("[K-CHOPORE] Counts file: %s\n", counts_file))
cat(sprintf("[K-CHOPORE] Output dir: %s\n", output_dir))
cat(sprintf("[K-CHOPORE] padj threshold: %s\n", padj_cutoff))
cat(sprintf("[K-CHOPORE] LFC threshold: %s\n", lfc_cutoff))

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read counts matrix from FLAIR output
counts_raw <- read.table(counts_file, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)

# Round counts to integers (FLAIR may output non-integer TPMs)
counts_matrix <- round(as.matrix(counts_raw))

# Remove rows with zero counts across all samples
counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, ]

cat(sprintf("[K-CHOPORE] Loaded %d features across %d samples\n",
            nrow(counts_matrix), ncol(counts_matrix)))

# Build sample metadata from column names
# Expected convention: condition_replicate (e.g., WT_C_R1, MUT_T_R1)
sample_names <- colnames(counts_matrix)

# Attempt to extract condition from sample names
# Strategy: everything except last _RN part is the condition
extract_condition <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 3) {
    # Join all parts except the last one (replicate identifier)
    paste(parts[1:(length(parts)-1)], collapse = "_")
  } else {
    name
  }
}

conditions <- sapply(sample_names, extract_condition)

coldata <- data.frame(
  condition = factor(conditions),
  row.names = sample_names
)

cat("[K-CHOPORE] Sample conditions:\n")
print(coldata)

# Check we have at least 2 conditions for differential analysis
if (length(unique(coldata$condition)) < 2) {
  cat("[K-CHOPORE] WARNING: Only one condition detected. DESeq2 requires at least 2 conditions.\n")
  cat("[K-CHOPORE] Generating descriptive statistics only.\n")

  # Write basic stats
  stats_df <- data.frame(
    gene = rownames(counts_matrix),
    mean_counts = rowMeans(counts_matrix),
    sd_counts = apply(counts_matrix, 1, sd)
  )
  write.csv(stats_df, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)

  # Create placeholder plots
  pdf(file.path(output_dir, "MA_plot.pdf"))
  plot(1, type = "n", main = "MA Plot - Requires 2+ conditions")
  dev.off()

  pdf(file.path(output_dir, "volcano_plot.pdf"))
  plot(1, type = "n", main = "Volcano Plot - Requires 2+ conditions")
  dev.off()

  cat("[K-CHOPORE] DESeq2 analysis completed (single condition mode).\n")
  quit(save = "no", status = 0)
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~ condition
)

# Run DESeq2
cat("[K-CHOPORE] Running DESeq2 differential expression...\n")
dds <- DESeq(dds)

# Extract results
res <- results(dds, alpha = padj_cutoff)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Add significance classification
res_df$significant <- ifelse(
  !is.na(res_df$padj) & res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) > lfc_cutoff,
  "significant", "not_significant"
)

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

# Write results
write.csv(res_df, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)

n_sig <- sum(res_df$significant == "significant", na.rm = TRUE)
cat(sprintf("[K-CHOPORE] Found %d significant genes (padj < %s, |LFC| > %s)\n",
            n_sig, padj_cutoff, lfc_cutoff))

# MA Plot
pdf(file.path(output_dir, "MA_plot.pdf"), width = 8, height = 6)
plotMA(res, main = "K-CHOPORE DESeq2 MA Plot", ylim = c(-5, 5))
dev.off()
cat("[K-CHOPORE] MA plot saved.\n")

# Volcano Plot
pdf(file.path(output_dir, "volcano_plot.pdf"), width = 8, height = 6)
volcano_df <- res_df[!is.na(res_df$padj), ]
p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("not_significant" = "grey", "significant" = "red")) +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "K-CHOPORE Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "bottom")
print(p)
dev.off()
cat("[K-CHOPORE] Volcano plot saved.\n")

# PCA Plot
vsd <- vst(dds, blind = FALSE)
pdf(file.path(output_dir, "PCA_plot.pdf"), width = 8, height = 6)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pca_var[1], "% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "% variance")) +
  theme_minimal() +
  labs(title = "K-CHOPORE PCA Plot")
print(p_pca)
dev.off()
cat("[K-CHOPORE] PCA plot saved.\n")

# Write normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(output_dir, "normalized_counts.csv"))

cat("[K-CHOPORE] DESeq2 analysis completed successfully.\n")
