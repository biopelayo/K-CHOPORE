# Direct RNA Sequencing of *anac017-1* Mutant Under Antimycin A Treatment in *Arabidopsis thaliana*: A K-CHOPORE Pipeline Analysis

## SHORT VERSION

---

## 1. Introduction

Mitochondrial retrograde regulation (MRR) enables plant cells to reprogram nuclear gene expression in response to organellar stress. In *Arabidopsis thaliana*, the ER-anchored NAC transcription factor ANAC017 is a central mediator of this pathway: upon mitochondrial dysfunction, ANAC017 is proteolytically released and translocates to the nucleus to activate stress-responsive genes, including Alternative Oxidase 1A (AOX1A). Antimycin A (AA), a Complex III inhibitor, is the standard pharmacological inducer of MRR.

Oxford Nanopore Technologies (ONT) direct RNA sequencing (DRS) enables native, full-length RNA analysis without reverse transcription, offering simultaneous access to gene expression, isoform diversity, and RNA modifications. Here, we applied the K-CHOPORE pipeline — a Dockerized Snakemake workflow for ONT DRS — to characterize the transcriptome of wild-type and *anac017-1* plants under control and AA conditions using a 2×2 factorial design.

---

## 2. Materials and Methods

### 2.1 Experimental Design

*A. thaliana* Col-0 (WT) and *anac017-1* rosette leaves (3 weeks old) were treated with 50 μM Antimycin A or DMSO control for 3 h. The final dataset comprised **10 samples**: WT×C (n=3), WT×AA (n=3), *anac017-1*×C (n=3), *anac017-1*×AA (n=1). Two additional *anac017-1*×AA replicates were excluded due to GPU basecalling requirements.

### 2.2 Sequencing and Analysis

Direct RNA libraries (SQK-RNA002) were sequenced on MinION R9.4.1 flow cells. The K-CHOPORE pipeline (v1.0) processed all data through: quality filtering (**NanoFilt**, Q≥7, ≥200 nt), splice-aware alignment to TAIR10 (**minimap2**, `-ax splice -uf -k14`), isoform analysis (**FLAIR** v2.0.0: align → correct → collapse → quantify), and differential expression (**DESeq2**, additive model `~ genotype + treatment`). QC was aggregated via **NanoPlot**, **NanoComp**, and **MultiQC**. The pipeline ran on a Dell PowerEdge T440 (40 threads) within a Docker container.

---

## 3. Results

### 3.1 Sequencing and Alignment

A total of **12.8 million raw reads** were generated (mean 1.28M/sample), with 98.8% passing quality filters. Mean read length was 997 nt (N50: 1,148 nt) and mean quality 9.1. Alignment to TAIR10 achieved a **90.5% mean mapping rate** (range: 77.1–94.7%). FLAIR quantified **20,958 isoforms** across all samples.

### 3.2 Principal Component Analysis

PCA on variance-stabilized expression revealed genotype as the dominant source of variation (PC1: 44%), with clear WT vs *anac017-1* separation. PC2 (28%) captured treatment-related variability. Together, PC1+PC2 explained 72% of total variance.

### 3.3 Genotype Effect (*anac017-1* vs WT)

DESeq2 identified **435 DE isoforms** (padj<0.05, |log₂FC|>1): 303 upregulated and 132 downregulated in *anac017-1*. Top hits included AT5G64120 (peroxidase, +1.46), AT4G09040 (RNA-binding protein, +1.96), and AT1G53480 (mitochondrial Mg²⁺ transporter, −4.16 — the strongest downregulated gene, consistent with disrupted mitochondrial function).

### 3.4 Treatment Effect (Antimycin A vs Control)

AA treatment induced **266 DE isoforms** with a striking directional bias: **245 upregulated (92.1%)** vs only 21 downregulated. Top genes included OPR1 (AT1G76680, +1.96), DTX1 (AT2G04040, +3.43), and critically, **AOX1A** (AT3G22370, +2.24, padj=1.3×10⁻¹⁷) — the canonical MRR marker, validating both the experimental system and the pipeline.

---

## 4. Discussion

This study demonstrates that ONT direct RNA-seq, analyzed through the K-CHOPORE pipeline, effectively resolves transcriptomic differences in the *Arabidopsis* MRR system. Three key findings emerge:

**First**, ANAC017 loss causes substantial basal transcriptomic reprogramming (435 DE isoforms), with the strong downregulation of a mitochondrial Mg²⁺ transporter (AT1G53480, log₂FC=−4.16) suggesting disrupted mitochondrial homeostasis even without stress. The 2:1 upregulation bias indicates ANAC017 may function as both activator and repressor.

**Second**, the AA response is overwhelmingly unidirectional (92% upregulated), reflecting coordinated activation of defense programs. Recovery of AOX1A as a top hit confirms the biological validity of the dataset and analysis.

**Third**, FLAIR-based isoform quantification (20,958 isoforms) enables resolution beyond gene-level analysis, exploiting the full-length read advantage of DRS.

**Limitations** include the unbalanced design (n=1 for *anac017-1*×AA, precluding interaction estimation), ELIGOS2 failure (rpy2/R compatibility issue preventing epitranscriptomic analysis), and the inability to perform m6Anet modification detection due to FAST5 storage constraints.

Future work should complete the factorial design with all 12 samples, resolve modification detection tools, and integrate m6A analysis to fully exploit the epitranscriptomic potential of direct RNA sequencing.

---

*Pipeline: K-CHOPORE v1.0 (https://github.com/biopelayo/K-CHOPORE) — Analysis: February 2026*
