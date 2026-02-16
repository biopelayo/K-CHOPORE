# Bioinformatics Methods — K-CHOPORE Pipeline

## 2.X Bioinformatics Analysis of Direct RNA Sequencing Data

### 2.X.1 Pipeline Overview

Direct RNA sequencing (dRNA-seq) data from Oxford Nanopore Technologies (ONT) MinION was analyzed using the K-CHOPORE pipeline, a reproducible Snakemake-based workflow (v7.32.4) encapsulated in a Docker container (Ubuntu 22.04). The pipeline integrates quality control, read alignment, isoform analysis, and epitranscriptomic modification detection for *Arabidopsis thaliana* (TAIR10 reference genome). All analyses were performed using the TAIR10 chromosome assembly and AtRTDv2 gene annotation.

### 2.X.2 Samples

Two ONT direct RNA sequencing libraries were analyzed:
- **WT_C_R1** (control): 1,343,980 reads, 1.35 Gb total bases
- **WT_C_R2** (treatment): 1,649,483 reads, 1.57 Gb total bases

FAST5 signal-level data was available for WT_C_R2, enabling nanopolish-based modification detection (m6Anet) for this sample.

### 2.X.3 Quality Control and Read Filtering

Raw reads were quality-filtered using NanoFilt (v2.8.0; minimum quality score Q7, minimum length 200 bp). Read quality statistics were assessed using NanoPlot (v1.43.0), which reported median read qualities of Q11.5 (R1) and Q11.1 (R2) with median read lengths of 881 bp and 836 bp, respectively. Read length N50 values were 1,152 bp (R1) and 1,070 bp (R2). Comparative analysis between samples was performed using NanoComp. Sequencing run quality was additionally evaluated using pycoQC using the MinION sequencing summary files.

### 2.X.4 Read Alignment

Filtered reads were aligned to the TAIR10 reference genome using minimap2 (v2.28) with splice-aware parameters optimized for ONT direct RNA sequencing (`-ax splice -uf -k14 --MD`). The `-uf` flag was used to indicate forward-stranded reads (characteristic of ONT dRNA-seq), and the `--MD` tag was included for downstream modification detection tools. Alignment files were sorted and indexed using samtools (v1.21).

Mapping rates were 92.4% (1,241,286 primary mapped reads) for R1 and 93.4% (1,540,497 primary mapped reads) for R2. Supplementary alignments accounted for 47,973 (R1) and 38,855 (R2) reads. The overall base-level error rate from CIGAR-based alignment was approximately 9.6% (R1) and 10.2% (R2), consistent with expected ONT direct RNA sequencing accuracy.

### 2.X.5 Isoform Analysis

Transcript isoform identification was performed using FLAIR (Full-Length Alternative Isoform analysis of RNA; v2.0). The FLAIR workflow consisted of four stages:

1. **Alignment**: Reads were aligned using FLAIR's splice-aware aligner against a chromosome-renamed genome (matching the AtRTDv2 GTF chromosome nomenclature: Chr1–Chr5, ChrM, ChrC).

2. **Correction**: Splice junctions were corrected using the AtRTDv2 annotation, yielding 1,214,468 corrected reads (R1) and 1,506,810 corrected reads (R2).

3. **Collapse**: Corrected reads were collapsed into unique transcript isoforms, identifying 21,286 isoforms (R1) and 23,294 isoforms (R2).

4. **Quantification**: Isoform expression was quantified across both samples using the FLAIR quantify module, producing a count matrix of 21,286 isoforms. Of these, 20,989 isoforms (98.6%) were detected in both samples. A total of 1,146,310 reads (R1) and 1,398,942 reads (R2) were assigned to isoforms.

The most highly expressed isoforms included AT1G67090 (rubisco small subunit family; ~25,000–28,000 reads per sample), consistent with the photosynthetic tissue origin of the RNA.

### 2.X.6 Epitranscriptomic Modification Detection

#### m6Anet (N6-methyladenosine detection)

Signal-level m6A modification detection was performed on WT_C_R2 (which had available FAST5 data) using the following workflow:

1. **Signal indexing**: FAST5 files were indexed and linked to basecalled reads using nanopolish index.

2. **Event alignment**: Nanopolish eventalign was used to align ionic current signal events to the reference sequence, producing a 70 GB event-level alignment file.

3. **m6Anet dataprep**: Event-level data was processed into the m6Anet input format.

4. **m6Anet inference**: m6A modification probabilities were computed using m6Anet (v2.1.0) with 1,000 iterations.

A total of 228 DRACH-motif sites were tested across all chromosomes. The per-chromosome distribution of tested sites was: Chr4 (100), Chr3 (40), chloroplast (33), Chr2 (18), Chr1 (18), Chr5 (15), and mitochondria (4). The relatively low number of sites reflects the stringent read depth requirement (minimum 20 reads at a given position).

Among the tested sites, 18 sites exhibited modification probabilities above 0.5, with 2 sites exceeding the 0.9 confidence threshold. The highest-confidence sites were located on chromosome 4 (position 130,288; probability 0.94, AGACT motif) and chromosome 3 (position 87,790; probability 0.91, AAACT motif). Chloroplast sites also showed elevated modification signals (probability 0.87 at position 104,865). The predominant kmer context for modified sites was GAACT (7/18 sites), followed by GAACA (3/18), consistent with the canonical DRACH m6A motif.

#### ELIGOS2 (Error-linked identification of genomic RNA modifications)

ELIGOS2 was run in `rna_mod` mode, comparing per-position error signatures in aligned reads against the built-in rBEM5+2 error model. FLAIR-collapsed isoform BED coordinates were used as the search regions (16,017 genomic intervals per sample). After processing all regions across both strands, ELIGOS2 did not identify statistically significant modification sites above the configured thresholds (p-value < 0.05, odds ratio > 5, excess error rate > 0.2). This may reflect the stringent filtering criteria or the comparison against a generic error model rather than a matched control condition.

### 2.X.7 Differential Expression Analysis

Differential isoform expression analysis using DESeq2 (v1.34.0) was not performed due to insufficient biological replicates (n=1 per condition). DESeq2 requires a minimum of two biological replicates per condition for proper dispersion estimation. The FLAIR count matrix is available for future analysis when additional replicates are generated.

### 2.X.8 Aggregated Quality Report

All quality metrics were aggregated into a unified report using MultiQC (v1.33), incorporating NanoPlot, samtools flagstat, and samtools stats outputs.

### 2.X.9 Software and Reproducibility

The complete analysis pipeline is available as a Docker image and Snakemake workflow at https://github.com/biopelayo/K-CHOPORE. Key software versions: Snakemake 7.32.4, minimap2 2.28, samtools 1.21, FLAIR 2.0, nanopolish 0.14.0, m6Anet 2.1.0, ELIGOS2 (GitLab: piroonj/eligos2), NanoPlot 1.43.0, NanoFilt 2.8.0, DESeq2 1.34.0, MultiQC 1.33. All tools were containerized in a single Docker image (Ubuntu 22.04, Python 3.10, R 4.1.2) for reproducibility.

---

### Summary Table

| Metric | WT_C_R1 (Control) | WT_C_R2 (Treatment) |
|--------|-------------------|---------------------|
| Total reads | 1,343,980 | 1,649,483 |
| Total bases | 1.35 Gb | 1.57 Gb |
| Mean read length | 1,007 bp | 953 bp |
| Median read quality | Q11.5 | Q11.1 |
| Read length N50 | 1,152 bp | 1,070 bp |
| Reads >Q10 | 1,113,095 (82.8%) | 1,259,301 (76.3%) |
| Mapped reads | 1,241,286 (92.4%) | 1,540,497 (93.4%) |
| Error rate (CIGAR) | 9.6% | 10.2% |
| FLAIR corrected reads | 1,214,468 | 1,506,810 |
| FLAIR isoforms | 21,286 | 23,294 |
| Reads assigned to isoforms | 1,146,310 | 1,398,942 |
| FAST5 available | No | Yes |
| m6Anet sites tested | — | 228 |
| m6Anet sites (prob > 0.5) | — | 18 |
| m6Anet sites (prob > 0.9) | — | 2 |
