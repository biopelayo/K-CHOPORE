# K-CHOPORE: A Comprehensive Non-Coding RNA and Epitranscriptomic Analysis Pipeline for ONT Direct RNA Sequencing

**Pelayo Gonzalez de Lena Rodriguez**
Cancer Epigenetics and Nanomedicine Lab | FINBA / Systems Biology Lab | University of Oviedo

---

## Abstract

K-CHOPORE (Keen Comprehensive High-throughput Omics Pipeline Organizer) is an open-source, containerized bioinformatics pipeline for the integrated analysis of Oxford Nanopore Technologies (ONT) direct RNA sequencing data. Version 3.0 extends the original basecalling-to-differential-expression workflow with five new modules for non-coding RNA biology: (1) long non-coding RNA discovery via multi-tool consensus classification, (2) small RNA profiling from complementary Illumina sRNA-seq data, (3) miRNA target prediction with degradome validation, (4) enhanced epitranscriptomic modification detection through multi-tool consensus, and (5) multi-omics data integration via co-expression network analysis and cis-regulatory element detection. The pipeline processes a 2x2 factorial experiment in *Arabidopsis thaliana* (WT vs *anac017-1* mutant x Control vs Antimycin A) and is designed for reproducibility through Snakemake workflow management and Docker containerization with 46 analysis rules and 136 automated tests.

---

## 1. Introduction

Direct RNA sequencing (DRS) on the Oxford Nanopore Technologies (ONT) platform enables the native analysis of full-length RNA molecules without reverse transcription or amplification, preserving information about RNA modifications, polyadenylation, and full-length isoform structure. This technology is particularly suited to the study of epitranscriptomic modifications such as N6-methyladenosine (m6A), which play critical roles in plant stress responses and development.

The K-CHOPORE pipeline was developed to provide a comprehensive, reproducible analysis framework that maximizes the biological information extractable from DRS experiments, while integrating complementary sequencing modalities (Illumina small RNA-seq, degradome-seq) for a systems-level view of post-transcriptional gene regulation.

### 1.1 Experimental Design

- **Organism**: *Arabidopsis thaliana* (Columbia-0 ecotype)
- **Reference genome**: TAIR10 with AtRTDv2 annotation
- **Genotypes**: Wild-type (WT) vs *anac017-1* mitochondrial retrograde signaling mutant
- **Treatment**: Control (C) vs Antimycin A (AA), a mitochondrial complex III inhibitor
- **Design**: 2x2 factorial with 10 samples (WT-C: 3, WT-AA: 3, anac017-1-C: 3, anac017-1-AA: 1)
- **Sequencing**: ONT MinION, R9.4.1 FLO-MIN106 flowcell, Direct RNA Sequencing Kit (SQK-RNA002)

---

## 2. Materials and Methods

### 2.1 Core Pipeline (v2.0)

#### 2.1.1 Basecalling
Raw FAST5/POD5 signal files were basecalled using Dorado (ONT) with the `rna004_130bps_hac` high-accuracy model. Legacy support for Guppy CPU basecalling is maintained for backward compatibility.

#### 2.1.2 Read Quality Control and Filtering
Basecalled reads were filtered using NanoFilt (Q >= 7, length >= 200 nt). Per-sample quality metrics were generated with NanoPlot, and cross-sample comparisons were performed with NanoComp.

#### 2.1.3 Splice-Aware Alignment
Filtered reads were aligned to the TAIR10 reference genome using minimap2 v2.24 with direct RNA-optimized parameters:
```
minimap2 -ax splice -uf -k14 --secondary=no --MD --junc-bed <junction_bed>
```
The `-uf` flag is critical for ONT direct RNA reads, which are forward-stranded. The `--MD` tag is required by downstream modification detection tools. Known splice junctions from AtRTDv2 were provided via `--junc-bed` to guide alignment. Alignments were sorted and indexed with samtools.

#### 2.1.4 Isoform Analysis
Full-length isoform analysis was performed with FLAIR v1.7:
1. **flair align**: Initial alignment to reference
2. **flair correct**: Splice site correction using AtRTDv2 annotation
3. **flair collapse**: Collapsing reads into isoform models (minimum support = 3)
4. **flair quantify**: Isoform-level quantification across all 10 samples
5. **flair diffExp/diffSplice**: Differential isoform expression and splicing

Alternative isoform assembly was performed with StringTie2 in long-read mode (`-L` flag).

#### 2.1.5 Epitranscriptomic Modification Detection (v2.0)
Three complementary approaches were employed:
- **ELIGOS2**: Identifies RNA modifications from base-level error patterns by comparing observed error rates against an IVT reference model (p < 0.05, oddR > 2.5, ESB > 0.1)
- **m6Anet**: Deep learning-based m6A detection from Nanopolish eventalign signal data. Identifies modification probability at DRACH motif sites
- **xPore**: Differential modification analysis between conditions using signal-level data

#### 2.1.6 Differential Expression Analysis
Differential expression was analyzed using DESeq2 with an additive model (`~ genotype + treatment`) due to the unbalanced design (n=1 in anac017-1 x AA group). Three contrasts were extracted: genotype effect (anac017-1 vs WT), treatment effect (AA vs C), and interaction (when estimable). Significance thresholds: padj < 0.05, |log2FoldChange| > 1.

### 2.2 Module 1: lncRNA Discovery (v3.0)

Long non-coding RNAs were identified from FLAIR-collapsed novel isoforms using a multi-tool consensus approach:

1. **Candidate identification**: FLAIR isoforms were classified against Araport11 using gffcompare. Transcripts with class codes u (intergenic), x (antisense), i (intronic), p (polymerase run-on), and y (containing a reference within intron) were selected as lncRNA candidates (>200 nt, >= 1 exon).

2. **ORF filtering**: TransDecoder was used to predict open reading frames. Candidates with ORFs >100 amino acids were excluded as likely protein-coding.

3. **Consensus coding potential assessment**: Three independent tools were applied:
   - **FEELnc**: Three-step pipeline (filter, codpot, classifier) for lncRNA classification into lincRNA, antisense, and intronic subtypes
   - **CPC2**: Coding Potential Calculator based on sequence features and ORF integrity
   - **CPAT**: Coding Potential Assessment Tool using Arabidopsis-trained hexamer and Fickett score models (non-coding threshold: coding probability < 0.39)

4. **Consensus rule**: A transcript was classified as lncRNA if >= 2 of 3 tools (FEELnc, CPC2, CPAT) agreed it was non-coding AND TransDecoder found no ORF > 100 aa.

5. **Annotation enrichment**: Consensus lncRNAs were cross-referenced against CANTATAdb v3 (known plant lncRNA database) and classified as known or novel. TE-derived lncRNAs were identified by intersection with TAIR10 transposable element annotations using bedtools.

6. **Differential expression**: lncRNA expression was analyzed with the same DESeq2 framework as coding genes, using extracted count rows from the FLAIR quantification matrix.

### 2.3 Module 2: Small RNA Analysis (v3.0)

Small RNA profiling was designed for complementary Illumina sRNA-seq data:

1. **Preprocessing**: Adapter trimming and size selection (18-30 nt) using Trim Galore
2. **sRNA locus discovery**: ShortStack v4 for joint analysis of all samples, identifying sRNA-producing loci and annotating MIRNA loci using strict plant criteria (miRNA/miRNA* duplex, stem-loop structure, precise Dicer processing)
3. **Novel miRNA prediction**: miRDeep-P2 for plant-specific novel miRNA discovery with secondary structure validation
4. **Annotation**: Cross-reference with miRBase v22 (Arabidopsis) and PmiREN to classify known vs novel miRNAs
5. **Quantification**: miRNA count matrix generation and sRNA size distribution analysis

### 2.4 Module 3: miRNA Target Prediction and Degradome Validation (v3.0)

1. **Computational prediction**: Two complementary approaches:
   - **TargetFinder**: Complementarity-based target prediction (score cutoff <= 4.0)
   - **psRNATarget**: Energy-based prediction with target-site accessibility model (expectation <= 3.0)

2. **Experimental validation**: CleaveLand4 analysis of PARE/degradome-seq data to identify miRNA-guided cleavage sites (category <= 2)

3. **Evidence integration**: Predictions were merged, cross-referenced with DESeq2 results to identify anti-correlated miRNA-target pairs, and lncRNA targets were specifically flagged for the regulatory network.

### 2.5 Module 4: Enhanced Epitranscriptomics (v3.0)

Multi-tool modification consensus was built by:

1. **Nanocompore analysis**: Signal-level differential modification detection for configured pairwise comparisons (min coverage = 30, GMM logit p-value < 0.05)

2. **Output standardization**: ELIGOS2, m6Anet, and Nanocompore results were parsed into a common format (chr, position, strand, modification type, p-value/probability)

3. **Consensus calling**: Sites called significant by >= 2 of 3 tools were designated high-confidence modification sites

4. **Biotype stratification**: Modification sites were classified by transcript biotype (mRNA vs lncRNA) using the Module 1 lncRNA annotation, enabling analysis of the epitranscriptomic landscape of non-coding RNAs

### 2.6 Module 5: Multi-Omics Data Integration (v3.0)

1. **Co-expression network analysis**: WGCNA (Weighted Gene Co-expression Network Analysis) was applied to combined mRNA + lncRNA + miRNA expression matrices. Soft-thresholding power was auto-detected, and network modules were correlated with genotype and treatment traits.

2. **Cis-regulatory lncRNA detection**: lncRNAs within 10 kb of differentially expressed protein-coding genes were identified and classified by orientation (sense/antisense) and relative position (upstream/downstream/overlapping).

3. **Population variation context**: ncRNA loci were intersected with 1001 Genomes Project variant data to identify loci with high natural variation, suggesting relaxed selective constraint or adaptive evolution.

4. **Integration report**: A comprehensive HTML report combining all module outputs for biological interpretation.

---

## 3. Pipeline Implementation

### 3.1 Architecture

K-CHOPORE v3.0 comprises 46 Snakemake rules organized in a directed acyclic graph (DAG). Each module is independently toggleable via boolean parameters in `config/config.yml`, enabling phased execution:

| Phase | Modules | Data Required |
|-------|---------|---------------|
| A | M1 (lncRNA) + M4 (Epitx) | ONT DRS data (existing) |
| B | M2 (Small RNA) + M3 (Targets) | Illumina sRNA-seq + degradome FASTQ |
| C | M5 (Integration) | All module outputs |

### 3.2 Containerization

All 35+ bioinformatics tools are pre-installed in a single Docker image. The Dockerfile follows a strict build order with tool-specific compilation strategies (e.g., Nanopolish serial compilation to avoid HDF5 race conditions). Critical NTFS/Windows compatibility patterns include latency-wait for mount propagation and `mkdir -p` guards in every rule.

### 3.3 Quality Assurance

The test suite comprises **136 automated tests** across 8 test files:
- Configuration validation (31 tests)
- Snakefile rule consistency and DAG structure (21 tests)
- Module-specific unit tests (46 tests)
- End-to-end integration tests with toy data (38 tests)

Integration tests exercise core Python scripts (lncrna_consensus.py, epitx_consensus.py, annotate_mirna.py) with synthetic mock data representing realistic pipeline outputs.

### 3.4 Reproducibility

- **FAIR-compliant**: Findable (GitHub), Accessible (open-source MIT license), Interoperable (standard file formats), Reusable (Docker + config-driven)
- **Version control**: All parameters, tool versions, and reference databases captured in config.yml
- **Logging**: Every rule logs to dedicated log files with `[K-CHOPORE]` prefixed audit messages
- **Annotation provenance**: `scripts/download_annotations.sh` records checksums and download dates

---

## 4. Results Summary (v2.0 Production Run)

### 4.1 Sequencing and Alignment
- 10 samples processed through the core pipeline
- Splice-aware alignment with minimap2 (average mapping rate: >95%)

### 4.2 Isoform Analysis
- **20,958 isoforms** quantified across samples by FLAIR collapse/quantify

### 4.3 Differential Expression
- **435 DE genes** by genotype (anac017-1 vs WT, padj < 0.05, |LFC| > 1)
- **266 DE genes** by treatment (AA vs C, padj < 0.05, |LFC| > 1)

---

## 5. Software Availability

K-CHOPORE is available at: https://github.com/biopelayo/K-CHOPORE

**Requirements**: Docker, Snakemake, Python 3.8+, Git

---

## 6. Key References

1. De Coster W, et al. (2018) NanoPack. *Bioinformatics*, 34(15):2666-2669.
2. Li H (2018) Minimap2. *Bioinformatics*, 34(18):3094-3100.
3. Tang AD, et al. (2020) Full-length transcript characterization of SF3B1 mutation in CLL. *Nature Communications*, 11:1438.
4. Jenjaroenpun P, et al. (2021) Decoding the epitranscriptional landscape from native RNA sequences. *Nature Communications*, 12:7091.
5. Hendra C, et al. (2022) Detection of m6A from direct RNA sequencing using a multiple instance learning framework. *Nature Methods*, 19:1590-1598.
6. Leger A, et al. (2021) RNA modifications detection by comparative Nanopore direct RNA sequencing. *Nature Communications*, 12:7025.
7. Love MI, et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550.
8. Wucher V, et al. (2017) FEELnc: a tool for long non-coding RNA annotation. *Nucleic Acids Research*, 45(5):e21.
9. Kang YJ, et al. (2017) CPC2: a fast and accurate coding potential calculator. *Nucleic Acids Research*, 45(W1):W12-W16.
10. Wang L, et al. (2013) CPAT: Coding-Potential Assessment Tool. *Nucleic Acids Research*, 41(6):e74.
11. Axtell MJ (2013) ShortStack: comprehensive annotation and quantification of small RNA genes. *RNA*, 19(6):740-751.
12. Fahlgren N, Carrington JC (2010) miRNA target prediction in plants. *Methods in Molecular Biology*, 592:51-57.
13. Dai X, et al. (2018) psRNATarget: a plant small RNA target analysis server. *Nucleic Acids Research*, 46(W1):W49-W54.
14. Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*, 9:559.
