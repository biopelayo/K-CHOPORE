#!/usr/bin/env Rscript
# =============================================================
# K-CHOPORE DESeq2 Differential Expression Analysis
# 2x2 Factorial Design: Genotype x Treatment
# =============================================================
# Performs differential gene/isoform expression analysis using
# a factorial design with explicit sample sheet.
#
# Usage:
#   Rscript run_deseq2.R <counts_matrix> <sample_sheet> <output_dir> <padj> <lfc>
#
# Input:
#   counts_matrix: FLAIR counts_matrix.tsv (genes/isoforms x samples)
#   sample_sheet:  TSV with columns: sample, genotype, treatment
#   output_dir:    Directory for results
#   padj:          Adjusted p-value threshold (e.g., 0.05)
#   lfc:           Log2 fold change threshold (e.g., 1.0)
#
# Output:
#   deseq2_genotype_results.csv   - Genotype effect (anac017_1 vs WT)
#   deseq2_treatment_results.csv  - Treatment effect (AA vs C)
#   deseq2_interaction_results.csv - Genotype x Treatment interaction
#   PCA_plot.pdf, MA_plot_*.pdf, volcano_*.pdf
#   normalized_counts.csv
# =============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_deseq2.R <counts_matrix> <sample_sheet> <output_dir> <padj> <lfc>")
}

counts_file    <- args[1]
sample_sheet   <- args[2]
output_dir     <- args[3]
padj_cutoff    <- as.numeric(args[4])
lfc_cutoff     <- as.numeric(args[5])

cat("[K-CHOPORE] DESeq2 Factorial Analysis (Genotype x Treatment)\n")
cat(sprintf("[K-CHOPORE] Counts file: %s\n", counts_file))
cat(sprintf("[K-CHOPORE] Sample sheet: %s\n", sample_sheet))
cat(sprintf("[K-CHOPORE] Output dir: %s\n", output_dir))
cat(sprintf("[K-CHOPORE] padj threshold: %s\n", padj_cutoff))
cat(sprintf("[K-CHOPORE] LFC threshold: %s\n", lfc_cutoff))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------
# Load data
# -------------------------------------------------------------

# Read counts matrix
counts_raw <- read.table(counts_file, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)
counts_matrix <- round(as.matrix(counts_raw))
counts_matrix <- counts_matrix[rowSums(counts_matrix) > 0, ]

cat(sprintf("[K-CHOPORE] Loaded %d features across %d samples\n",
            nrow(counts_matrix), ncol(counts_matrix)))

# Read sample sheet
coldata <- read.table(sample_sheet, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample

# Verify sample names match between counts and sample sheet
counts_samples <- colnames(counts_matrix)
sheet_samples <- coldata$sample

missing_in_sheet <- setdiff(counts_samples, sheet_samples)
missing_in_counts <- setdiff(sheet_samples, counts_samples)

if (length(missing_in_sheet) > 0) {
  cat(sprintf("[K-CHOPORE] WARNING: Samples in counts but not in sample sheet: %s\n",
              paste(missing_in_sheet, collapse = ", ")))
}
if (length(missing_in_counts) > 0) {
  cat(sprintf("[K-CHOPORE] WARNING: Samples in sheet but not in counts: %s\n",
              paste(missing_in_counts, collapse = ", ")))
}

# Use only samples present in both
common_samples <- intersect(counts_samples, sheet_samples)
cat(sprintf("[K-CHOPORE] Using %d samples present in both counts and sample sheet\n",
            length(common_samples)))

counts_matrix <- counts_matrix[, common_samples]
coldata <- coldata[common_samples, ]

# Set up factors with explicit reference levels
coldata$genotype <- factor(coldata$genotype, levels = c("WT", sort(setdiff(unique(coldata$genotype), "WT"))))
coldata$treatment <- factor(coldata$treatment, levels = c("C", sort(setdiff(unique(coldata$treatment), "C"))))

cat("[K-CHOPORE] Experimental design:\n")
cat(sprintf("  Genotypes: %s (reference: %s)\n",
            paste(levels(coldata$genotype), collapse = ", "),
            levels(coldata$genotype)[1]))
cat(sprintf("  Treatments: %s (reference: %s)\n",
            paste(levels(coldata$treatment), collapse = ", "),
            levels(coldata$treatment)[1]))
cat("[K-CHOPORE] Sample sheet:\n")
print(coldata[, c("sample", "genotype", "treatment")])

# Verify replication
design_table <- table(coldata$genotype, coldata$treatment)
cat("\n[K-CHOPORE] Replicates per group:\n")
print(design_table)

if (any(design_table < 2)) {
  stop("[K-CHOPORE] ERROR: At least 2 replicates per group required for DESeq2.")
}

# -------------------------------------------------------------
# DESeq2 factorial analysis
# -------------------------------------------------------------

# Create DESeq2 dataset with factorial design
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~ genotype + treatment + genotype:treatment
)

cat("[K-CHOPORE] Running DESeq2 with design: ~ genotype + treatment + genotype:treatment\n")
dds <- DESeq(dds)

cat("[K-CHOPORE] Available result names:\n")
print(resultsNames(dds))

# -------------------------------------------------------------
# Helper: extract and save results for a contrast
# -------------------------------------------------------------
save_contrast_results <- function(dds, name, contrast_name, output_dir, padj_cutoff, lfc_cutoff) {
  res <- results(dds, name = name, alpha = padj_cutoff)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)

  res_df$significant <- ifelse(
    !is.na(res_df$padj) & res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) > lfc_cutoff,
    "significant", "not_significant"
  )

  res_df <- res_df[order(res_df$padj), ]
  n_sig <- sum(res_df$significant == "significant", na.rm = TRUE)
  cat(sprintf("[K-CHOPORE] %s: %d significant genes (padj < %s, |LFC| > %s)\n",
              contrast_name, n_sig, padj_cutoff, lfc_cutoff))

  # Write results CSV
  write.csv(res_df,
            file.path(output_dir, paste0("deseq2_", contrast_name, "_results.csv")),
            row.names = FALSE)

  # MA plot
  pdf(file.path(output_dir, paste0("MA_plot_", contrast_name, ".pdf")), width = 8, height = 6)
  plotMA(res, main = paste("MA Plot -", contrast_name), ylim = c(-5, 5))
  dev.off()

  # Volcano plot
  volcano_df <- res_df[!is.na(res_df$padj), ]
  if (nrow(volcano_df) > 0) {
    pdf(file.path(output_dir, paste0("volcano_", contrast_name, ".pdf")), width = 8, height = 6)
    p <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(values = c("not_significant" = "grey60", "significant" = "red")) +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "blue") +
      theme_minimal() +
      labs(title = paste("Volcano Plot -", contrast_name),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      theme(legend.position = "bottom")
    print(p)
    dev.off()
  }

  return(n_sig)
}

# -------------------------------------------------------------
# Extract contrasts
# -------------------------------------------------------------

# Get result names to determine available contrasts
result_names <- resultsNames(dds)

# 1. Genotype effect (e.g., genotype_anac017_1_vs_WT)
genotype_name <- grep("^genotype_", result_names, value = TRUE)
if (length(genotype_name) > 0) {
  cat(sprintf("\n[K-CHOPORE] --- Genotype contrast: %s ---\n", genotype_name[1]))
  save_contrast_results(dds, genotype_name[1], "genotype", output_dir, padj_cutoff, lfc_cutoff)
} else {
  cat("[K-CHOPORE] WARNING: No genotype contrast found in results.\n")
}

# 2. Treatment effect (e.g., treatment_AA_vs_C)
treatment_name <- grep("^treatment_", result_names, value = TRUE)
if (length(treatment_name) > 0) {
  cat(sprintf("\n[K-CHOPORE] --- Treatment contrast: %s ---\n", treatment_name[1]))
  save_contrast_results(dds, treatment_name[1], "treatment", output_dir, padj_cutoff, lfc_cutoff)
} else {
  cat("[K-CHOPORE] WARNING: No treatment contrast found in results.\n")
}

# 3. Interaction effect (genotype:treatment)
interaction_name <- grep("genotype.*treatment|treatment.*genotype", result_names, value = TRUE)
if (length(interaction_name) > 0) {
  cat(sprintf("\n[K-CHOPORE] --- Interaction contrast: %s ---\n", interaction_name[1]))
  save_contrast_results(dds, interaction_name[1], "interaction", output_dir, padj_cutoff, lfc_cutoff)
} else {
  cat("[K-CHOPORE] WARNING: No interaction contrast found in results.\n")
}

# -------------------------------------------------------------
# PCA Plot (colored by genotype, shaped by treatment)
# -------------------------------------------------------------
cat("\n[K-CHOPORE] Generating PCA plot...\n")
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("genotype", "treatment"), returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"))

pdf(file.path(output_dir, "PCA_plot.pdf"), width = 10, height = 7)
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = genotype, shape = treatment)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", pca_var[1], "% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "% variance")) +
  theme_minimal(base_size = 14) +
  labs(title = "K-CHOPORE PCA: Genotype x Treatment",
       color = "Genotype", shape = "Treatment") +
  theme(legend.position = "right")
print(p_pca)
dev.off()
cat("[K-CHOPORE] PCA plot saved.\n")

# -------------------------------------------------------------
# Normalized counts
# -------------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(output_dir, "normalized_counts.csv"))
cat("[K-CHOPORE] Normalized counts saved.\n")

cat("\n[K-CHOPORE] DESeq2 factorial analysis completed successfully.\n")
