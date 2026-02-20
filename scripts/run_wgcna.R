#!/usr/bin/env Rscript
# =============================================================
# K-CHOPORE v3.0 — WGCNA Co-expression Network Analysis
# =============================================================
# Builds a weighted gene co-expression network from combined
# mRNA + lncRNA + miRNA count matrices. Identifies modules and
# correlates them with experimental conditions.
#
# Usage:
#   Rscript run_wgcna.R <mrna_counts> <sample_sheet> <lncrna_counts> <mirna_counts> <outdir> <soft_power> <min_module_size>
# =============================================================

suppressPackageStartupMessages({
  library(WGCNA)
})

# Allow multi-threading
allowWGCNAThreads()

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript run_wgcna.R <mrna_counts> <sample_sheet> <lncrna_counts> <mirna_counts> <outdir> <soft_power> <min_module_size>")
}

mrna_file     <- args[1]
sample_file   <- args[2]
lncrna_file   <- args[3]  # May be empty string
mirna_file    <- args[4]  # May be empty string
output_dir    <- args[5]
soft_power    <- args[6]  # "auto" or integer
min_mod_size  <- as.integer(args[7])

cat("[K-CHOPORE] WGCNA Co-expression Network Analysis\n")
cat(sprintf("[K-CHOPORE] mRNA counts: %s\n", mrna_file))
cat(sprintf("[K-CHOPORE] Output dir: %s\n", output_dir))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "wgcna_plots"), recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------
# Load and merge count matrices
# -------------------------------------------------------------
cat("[K-CHOPORE] Loading count matrices...\n")

# mRNA counts (required)
mrna_counts <- read.table(mrna_file, header = TRUE, sep = "\t",
                          row.names = 1, check.names = FALSE)
cat(sprintf("[K-CHOPORE] mRNA features: %d\n", nrow(mrna_counts)))

# lncRNA counts (optional)
if (nchar(lncrna_file) > 0 && file.exists(lncrna_file)) {
  lncrna_counts <- read.table(lncrna_file, header = TRUE, sep = "\t",
                              row.names = 1, check.names = FALSE)
  cat(sprintf("[K-CHOPORE] lncRNA features: %d\n", nrow(lncrna_counts)))
  # Ensure same samples
  common <- intersect(colnames(mrna_counts), colnames(lncrna_counts))
  if (length(common) > 0) {
    mrna_counts <- rbind(mrna_counts[, common], lncrna_counts[, common])
  }
}

# miRNA counts (optional)
if (nchar(mirna_file) > 0 && file.exists(mirna_file)) {
  mirna_counts <- read.table(mirna_file, header = TRUE, sep = "\t",
                             row.names = 1, check.names = FALSE)
  cat(sprintf("[K-CHOPORE] miRNA features: %d\n", nrow(mirna_counts)))
  common <- intersect(colnames(mrna_counts), colnames(mirna_counts))
  if (length(common) > 0) {
    mrna_counts <- rbind(mrna_counts[, common], mirna_counts[, common])
  }
}

cat(sprintf("[K-CHOPORE] Combined features: %d across %d samples\n",
            nrow(mrna_counts), ncol(mrna_counts)))

# Filter low-expression genes (keep if >10 counts in at least 3 samples)
keep <- rowSums(mrna_counts > 10) >= 3
mrna_counts <- mrna_counts[keep, ]
cat(sprintf("[K-CHOPORE] After filtering: %d features\n", nrow(mrna_counts)))

# Load sample metadata
coldata <- read.table(sample_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample

# Match samples
common_samples <- intersect(colnames(mrna_counts), coldata$sample)
mrna_counts <- mrna_counts[, common_samples]
coldata <- coldata[common_samples, ]

# Log-transform
datExpr <- t(log2(as.matrix(mrna_counts) + 1))

cat(sprintf("[K-CHOPORE] Expression matrix: %d samples x %d features\n",
            nrow(datExpr), ncol(datExpr)))

# -------------------------------------------------------------
# Choose soft-thresholding power
# -------------------------------------------------------------
if (soft_power == "auto") {
  cat("[K-CHOPORE] Auto-detecting soft-thresholding power...\n")
  powers <- c(1:20)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

  # Pick power where R^2 > 0.8
  best_power <- sft$powerEstimate
  if (is.na(best_power)) best_power <- 6  # Default fallback

  cat(sprintf("[K-CHOPORE] Selected soft power: %d\n", best_power))

  # Save soft threshold plot
  pdf(file.path(output_dir, "wgcna_plots", "soft_threshold.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
       main = "Scale independence", type = "n")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, col = "red")
  abline(h = 0.8, col = "blue")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       main = "Mean connectivity", type = "n")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
  dev.off()
} else {
  best_power <- as.integer(soft_power)
  cat(sprintf("[K-CHOPORE] Using specified soft power: %d\n", best_power))
}

# -------------------------------------------------------------
# Build network and detect modules
# -------------------------------------------------------------
cat("[K-CHOPORE] Building network (this may take a while)...\n")

net <- blockwiseModules(
  datExpr,
  power = best_power,
  TOMType = "unsigned",
  minModuleSize = min_mod_size,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 0
)

cat(sprintf("[K-CHOPORE] Modules detected: %d\n", max(net$colors)))
module_table <- table(net$colors)
cat("[K-CHOPORE] Module sizes:\n")
print(module_table)

# Assign module colors
moduleColors <- labels2colors(net$colors)

# Save dendrogram
pdf(file.path(output_dir, "wgcna_plots", "dendrogram_modules.pdf"), width = 12, height = 8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "K-CHOPORE: Gene dendrogram and module colors")
dev.off()

# -------------------------------------------------------------
# Module-trait correlations
# -------------------------------------------------------------
cat("[K-CHOPORE] Correlating modules with traits...\n")

# Create binary trait matrix
traits <- data.frame(
  genotype = as.numeric(coldata$genotype != "WT"),
  treatment = as.numeric(coldata$treatment != "C"),
  row.names = rownames(coldata)
)

# Module eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

# Correlations
moduleTraitCor <- cor(MEs, traits[rownames(datExpr), ], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Save heatmap
pdf(file.path(output_dir, "wgcna_plots", "module_trait_heatmap.pdf"), width = 8, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = "K-CHOPORE: Module-Trait Relationships")
dev.off()

# -------------------------------------------------------------
# Save results
# -------------------------------------------------------------

# Module assignments
module_df <- data.frame(
  gene = colnames(datExpr),
  module_number = net$colors,
  module_color = moduleColors,
  stringsAsFactors = FALSE
)

write.table(module_df,
            file.path(output_dir, "wgcna_modules.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Module-trait correlations
trait_df <- data.frame(
  module = rownames(moduleTraitCor),
  cor_genotype = moduleTraitCor[, "genotype"],
  pval_genotype = moduleTraitPvalue[, "genotype"],
  cor_treatment = moduleTraitCor[, "treatment"],
  pval_treatment = moduleTraitPvalue[, "treatment"],
  stringsAsFactors = FALSE
)

write.table(trait_df,
            file.path(output_dir, "wgcna_module_traits.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Module eigengenes
write.table(MEs,
            file.path(output_dir, "wgcna_eigengenes.tsv"),
            sep = "\t", quote = FALSE)

cat(sprintf("\n[K-CHOPORE] WGCNA completed. %d modules across %d features.\n",
            max(net$colors), ncol(datExpr)))
cat("[K-CHOPORE] Output files:\n")
cat(sprintf("[K-CHOPORE]   Modules: %s\n", file.path(output_dir, "wgcna_modules.tsv")))
cat(sprintf("[K-CHOPORE]   Traits: %s\n", file.path(output_dir, "wgcna_module_traits.tsv")))
cat(sprintf("[K-CHOPORE]   Plots: %s\n", file.path(output_dir, "wgcna_plots")))
