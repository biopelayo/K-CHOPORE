# Direct RNA Sequencing Reveals Transcriptomic Reprogramming in the *anac017-1* Mutant Under Mitochondrial Stress in *Arabidopsis thaliana*

## LONG VERSION

---

## 1. Introduction

Mitochondria are essential organelles in plant cells, responsible for energy production through oxidative phosphorylation, as well as playing key roles in stress signaling, programmed cell death, and metabolic regulation. In *Arabidopsis thaliana*, mitochondrial dysfunction triggers a retrograde signaling pathway that communicates organellar stress to the nucleus, leading to widespread transcriptional reprogramming. This mitochondrial retrograde regulation (MRR) is critical for plant acclimation to both biotic and abiotic stresses.

ANAC017 is a NAC-domain transcription factor anchored to the endoplasmic reticulum membrane that functions as a key mediator of mitochondrial retrograde signaling in *Arabidopsis*. Under conditions of mitochondrial dysfunction, ANAC017 is proteolytically cleaved and translocates to the nucleus, where it activates the expression of genes involved in alternative respiration, antioxidant defense, and stress tolerance. Antimycin A (AA), a specific inhibitor of Complex III of the mitochondrial electron transport chain, is widely used as a pharmacological tool to induce mitochondrial dysfunction and activate MRR pathways.

While previous studies have employed short-read RNA-seq to characterize transcriptional responses to mitochondrial stress, these approaches are limited in their ability to detect full-length transcript isoforms and RNA modifications. Oxford Nanopore Technologies (ONT) direct RNA sequencing (DRS) offers the unique advantage of sequencing native, full-length RNA molecules without reverse transcription or amplification, enabling simultaneous analysis of gene expression, isoform usage, and epitranscriptomic modifications in a single experiment.

In this study, we applied the K-CHOPORE pipeline — a comprehensive Snakemake-based workflow for ONT direct RNA-seq analysis — to investigate the transcriptomic landscape of wild-type (WT) and *anac017-1* mutant *Arabidopsis* plants under control conditions and Antimycin A treatment. Using a 2×2 factorial experimental design, we aimed to: (1) identify genes differentially expressed between genotypes and treatments, (2) characterize isoform-level expression patterns using FLAIR, and (3) assess the contribution of ANAC017 to the mitochondrial stress response.

---

## 2. Materials and Methods

### 2.1 Plant Material and Experimental Design

*Arabidopsis thaliana* wild-type (Col-0) and *anac017-1* T-DNA insertion mutant plants were grown under controlled conditions (22°C, 16h light/8h dark, 120 μmol·m⁻²·s⁻¹). Three-week-old rosette leaves were treated with either 50 μM Antimycin A (AA) dissolved in DMSO or DMSO alone (control, C) for 3 hours before harvest.

The experiment followed a **2×2 factorial design** with two genotypes (WT, *anac017-1*) and two treatments (Control, Antimycin A), with three biological replicates per group. Due to data availability constraints (two *anac017-1* × AA replicates required GPU-dependent basecalling and were excluded), the final dataset comprised **10 samples**: WT×C (n=3), WT×AA (n=3), *anac017-1*×C (n=3), and *anac017-1*×AA (n=1).

### 2.2 RNA Extraction and Library Preparation

Total RNA was extracted using TRIzol reagent followed by poly(A) selection. Direct RNA sequencing libraries were prepared using the Oxford Nanopore Technologies SQK-RNA002 kit following the manufacturer's protocol. Libraries were sequenced on an ONT MinION device using R9.4.1 flow cells (FLO-MIN106D).

### 2.3 K-CHOPORE Bioinformatic Pipeline

All analyses were performed using the **K-CHOPORE pipeline** (v1.0), a Docker-containerized Snakemake workflow integrating multiple bioinformatics tools for ONT direct RNA-seq analysis. The pipeline was executed on a Dell PowerEdge T440 server (40 CPU threads, 480 GB SSD).

#### 2.3.1 Read Processing and Quality Control
Pre-basecalled FASTQ files from MinKNOW (Guppy HAC, RNA model `rna_r9.4.1_70bps_hac`) were merged per sample (multiple flow cell runs concatenated). Quality filtering was performed with **NanoFilt** v2.8.0 (minimum quality score Q ≥ 7, minimum read length ≥ 200 nt). Per-sample quality metrics were generated with **NanoPlot** v1.42.0 and cross-sample comparisons with **NanoComp** v1.24.0.

#### 2.3.2 Alignment
Filtered reads were aligned to the *Arabidopsis thaliana* TAIR10 reference genome using **minimap2** v2.28 with splice-aware settings optimized for ONT direct RNA sequencing: `-ax splice -uf -k14 --junc-bed <GTF> --secondary=no --MD`. The `-uf` flag ensures forward-strand alignment, critical for direct RNA reads. BAM files were sorted and indexed with **samtools** v1.21. Alignment QC was assessed with samtools flagstat/stats and **pycoQC** v2.5.2.

#### 2.3.3 Isoform Analysis
Isoform detection and quantification were performed with **FLAIR** v2.0.0:
- **flair align**: Initial read alignment to the TAIR10 genome
- **flair correct**: Splice junction correction using the AtRTDv2 annotation (April 2016)
- **flair collapse**: Isoform model generation with a minimum support of 3 reads
- **flair quantify**: Isoform-level expression quantification across all samples using a reads manifest

Chromosome nomenclature was harmonized between TAIR10 (numeric: 1–5) and the AtRTDv2 GTF annotation (prefixed: Chr1–Chr5) at each pipeline stage.

#### 2.3.4 Epitranscriptomic Modification Detection
RNA modification analysis was attempted using **ELIGOS2** v2.1.0 for base-level error signal analysis. ELIGOS2 parameters: p-value threshold = 0.05, odds ratio ≥ 2.5, error signal bias ≥ 0.1, minimum depth = 50. **Note:** ELIGOS2 analyses failed across all samples due to a known Cochran-Mantel-Haenszel (CMH) test failure in the rpy2/R bridge within the Docker environment. This represents a known limitation of the current pipeline version.

Signal-level modification detection with **m6Anet** (via **Nanopolish** eventalign) was not performed as it requires raw FAST5 signal files on local disk (~644 GB total), exceeding available storage.

#### 2.3.5 Differential Expression Analysis
Differential expression analysis was performed with **DESeq2** v1.42.0 in R using the FLAIR isoform counts matrix (20,958 isoforms × 10 samples). Given the unbalanced design (1 replicate in the *anac017-1* × AA group), an **additive model** (`~ genotype + treatment`) was used instead of a full factorial model with interaction term. Reference levels were set to WT (genotype) and Control (treatment). Genes were considered significant at adjusted p-value (Benjamini-Hochberg) < 0.05 and |log₂FC| > 1. Variance-stabilized counts (VST) were used for PCA visualization.

#### 2.3.6 Quality Aggregation
An aggregate QC report was generated with **MultiQC** v1.25.1 combining outputs from NanoPlot, samtools, and pycoQC.

### 2.4 Software and Reproducibility

The complete analysis environment was containerized in a Docker image (k-chopore:latest, 22.7 GB) containing all tools and dependencies. The pipeline source code, configuration files, and Snakefile are available at https://github.com/biopelayo/K-CHOPORE.

---

## 3. Results

### 3.1 Sequencing Yield and Quality

Direct RNA sequencing of 10 *Arabidopsis* samples generated a total of **12,763,204 raw reads** (mean: 1,276,320 reads/sample; range: 364,047–1,671,811). After quality filtering (Q ≥ 7, length ≥ 200 nt), **12,603,968 reads** were retained (98.8% retention rate).

| Sample | Genotype | Treatment | Raw Reads | Filtered Reads | Mean Length (nt) | Mean Quality | N50 (nt) |
|--------|----------|-----------|-----------|----------------|------------------|--------------|----------|
| WT_C_R1 | WT | C | 1,271,877 | 1,258,927 | 1,012 | 9.3 | 1,162 |
| WT_C_R2 | WT | C | 1,533,146 | 1,511,111 | 960 | 9.2 | 1,076 |
| WT_C_R3 | WT | C | 1,195,093 | 1,183,990 | 1,011 | 9.2 | 1,155 |
| WT_AA_R1 | WT | AA | 1,290,274 | 1,281,147 | 1,032 | 9.4 | 1,293 |
| WT_AA_R2 | WT | AA | 1,352,754 | 1,332,072 | 961 | 8.9 | 1,068 |
| WT_AA_R3 | WT | AA | 1,537,278 | 1,515,132 | 900 | 8.9 | 1,004 |
| anac017_1_C_R1 | anac017-1 | C | 1,671,811 | 1,651,309 | 1,046 | 9.0 | 1,252 |
| anac017_1_C_R2 | anac017-1 | C | 1,126,095 | 1,108,328 | 998 | 8.9 | 1,156 |
| anac017_1_C_R3 | anac017-1 | C | 364,047 | 357,862 | 958 | 8.8 | 1,083 |
| anac017_1_AA_R1 | anac017-1 | AA | 1,420,829 | 1,404,090 | 1,047 | 8.9 | 1,235 |

**Table 1.** Sequencing statistics for 10 direct RNA-seq samples. Mean read length: 997 nt (range: 900–1,047). Mean quality: 9.1 (range: 8.8–9.4). Mean N50: 1,148 nt.

One sample (anac017_1_C_R3) yielded notably fewer reads (357,862) compared to the others, likely due to lower input RNA quantity or flow cell performance.

### 3.2 Alignment

Filtered reads were aligned to the TAIR10 reference genome with an overall mapping rate of **90.5%** (range: 77.1–94.7%).

| Sample | Total Alignments | Mapped | Mapping Rate |
|--------|-----------------|--------|-------------|
| WT_C_R1 | 1,303,295 | 1,203,856 | 92.4% |
| WT_C_R2 | 1,546,174 | 1,441,919 | 93.3% |
| WT_C_R3 | 1,220,299 | 1,132,427 | 92.8% |
| WT_AA_R1 | 1,296,023 | 999,132 | 77.1% |
| WT_AA_R2 | 1,378,672 | 1,306,186 | 94.7% |
| WT_AA_R3 | 1,559,371 | 1,468,102 | 94.2% |
| anac017_1_C_R1 | 1,677,460 | 1,511,222 | 90.1% |
| anac017_1_C_R2 | 1,128,366 | 1,047,368 | 92.8% |
| anac017_1_C_R3 | 360,482 | 315,406 | 87.5% |
| anac017_1_AA_R1 | 1,423,044 | 1,288,809 | 90.6% |

**Table 2.** Alignment statistics. Sample WT_AA_R1 showed a lower mapping rate (77.1%), possibly due to adapter contamination or degraded RNA. All other samples exceeded 87%.

### 3.3 Isoform Detection and Quantification

FLAIR identified between **11,548 and 25,611 isoform models per sample** after collapse (minimum 3 supporting reads). The combined quantification across all samples yielded **20,958 isoforms** in the counts matrix, representing expression across the *Arabidopsis* transcriptome.

### 3.4 Principal Component Analysis

PCA of variance-stabilized expression data revealed clear separation by genotype along PC1, which accounted for **44% of total variance** (Figure 1). PC2 (28% variance) captured additional variability including treatment effects. Key observations:

- WT samples (control and AA-treated) formed a cluster on the positive PC1 axis, with AA-treated samples shifted slightly relative to controls.
- *anac017-1* control samples clustered tightly on the negative PC1 axis.
- The single *anac017-1* × AA sample occupied an outlier position at the bottom of the PCA space, possibly reflecting a unique transcriptomic state or the limitation of n=1 representation.
- Together, PC1 and PC2 explained **72% of total variance**, indicating strong biological signals.

### 3.5 Differential Expression: Genotype Effect (*anac017-1* vs WT)

DESeq2 analysis identified **435 differentially expressed isoforms** (padj < 0.05, |log₂FC| > 1) between *anac017-1* and WT, with a notable bias toward upregulation in the mutant:

- **303 upregulated** in *anac017-1* (69.7%)
- **132 downregulated** in *anac017-1* (30.3%)

The top significant genes included:

| Gene | log₂FC | padj | Description |
|------|--------|------|-------------|
| AT5G64120 | +1.46 | 4.7×10⁻²¹ | Peroxidase superfamily |
| AT1G70780 | +1.33 | 1.6×10⁻¹⁹ | Unknown protein |
| AT4G09040 | +1.96 | 1.5×10⁻¹⁸ | RNA-binding protein |
| AT1G53480 | −4.16 | 1.5×10⁻¹⁸ | MRS2-type magnesium transporter |
| AT1G34190 | +1.79 | 6.7×10⁻¹⁷ | NAC-domain protein |
| AT1G21525 | +2.61 | 1.2×10⁻¹⁶ | Unknown protein |

**Table 3.** Top differentially expressed genes by genotype (anac017-1 vs WT). The strongest downregulated gene, AT1G53480 (log₂FC = −4.16), encodes a mitochondrial magnesium transporter, consistent with disrupted mitochondrial function in the *anac017-1* mutant.

### 3.6 Differential Expression: Treatment Effect (Antimycin A vs Control)

The Antimycin A treatment induced **266 differentially expressed isoforms**, with a striking directional bias toward upregulation:

- **245 upregulated** by AA (92.1%)
- **21 downregulated** by AA (7.9%)

This strongly asymmetric response reflects the known activation of mitochondrial retrograde signaling by Complex III inhibition.

Top AA-responsive genes:

| Gene | log₂FC | padj | Description |
|------|--------|------|-------------|
| AT1G76680 | +1.96 | 5.3×10⁻³⁷ | 12-oxophytodienoate reductase (OPR1) |
| AT2G26560 | +1.70 | 1.6×10⁻²⁶ | PLP-dependent transferase |
| AT1G76600 | +2.35 | 8.8×10⁻²³ | Unknown protein |
| AT4G15760 | +2.91 | 1.7×10⁻²² | Monooxygenase (MO1) |
| AT1G05675 | +3.23 | 2.2×10⁻²⁰ | UDP-glucuronosyltransferase |
| AT4G30280 | +2.52 | 2.2×10⁻¹⁹ | Xyloglucan endotransglucosylase (XTH18) |
| AT2G04040 | +3.43 | 2.5×10⁻¹⁸ | MATE efflux family (DTX1) |
| AT3G02040 | +2.56 | 1.1×10⁻¹⁷ | Unknown protein (SRG3-like) |
| AT3G22370 | +2.24 | 1.3×10⁻¹⁷ | Alternative oxidase 1A (AOX1A) |

**Table 4.** Top differentially expressed genes by treatment (AA vs Control). The identification of AT3G22370 (AOX1A, log₂FC = +2.24) among the top hits validates the experimental system, as Alternative Oxidase 1A is the canonical marker of mitochondrial retrograde signaling activated by Complex III inhibition.

### 3.7 Epitranscriptomic Modification Detection

ELIGOS2 analysis for RNA modification detection failed across all 10 samples due to a CMH statistical test error within the rpy2/R interface in the Docker environment. Placeholder output files were generated by the pipeline's fault-tolerant error handling. Signal-level m6A detection with m6Anet was not performed due to FAST5 storage constraints. Epitranscriptomic analysis remains a priority for future work.

---

## 4. Discussion

### 4.1 Genotype-Dependent Transcriptomic Reprogramming

Our direct RNA-seq analysis reveals that the loss of ANAC017 function causes substantial transcriptomic alterations even under basal (non-stress) conditions, with 435 differentially expressed isoforms. The dominance of upregulated genes (303 of 435, 69.7%) in *anac017-1* suggests that ANAC017 may act as both an activator and repressor of gene expression, or that its absence triggers compensatory transcriptional programs. The strong downregulation of AT1G53480 (mitochondrial magnesium transporter, log₂FC = −4.16) is particularly noteworthy, as magnesium homeostasis is critical for mitochondrial function and this gene's suppression may reflect or contribute to altered mitochondrial physiology in the mutant.

The clear separation of genotypes along PC1 (44% variance) confirms that ANAC017 is a major determinant of the *Arabidopsis* transcriptome, consistent with its established role as a master regulator of mitochondrial retrograde signaling.

### 4.2 Antimycin A Activates a Robust Unidirectional Stress Response

The Antimycin A treatment elicited a highly directional transcriptional response, with 92.1% of differentially expressed genes being upregulated. This pattern is consistent with the activation of mitochondrial retrograde signaling by Complex III inhibition and reflects a coordinated stress defense program.

The identification of **AOX1A** (AT3G22370, log₂FC = +2.24, padj = 1.3×10⁻¹⁷) among the most significantly upregulated genes provides strong biological validation of our experimental system and analytical pipeline. AOX1A encodes Alternative Oxidase 1A, the best-characterized target of ANAC017-mediated mitochondrial retrograde signaling, whose induction allows electrons to bypass the inhibited cytochrome pathway. Additional highly upregulated genes included **OPR1** (AT1G76680), involved in jasmonate biosynthesis and oxidative stress, and **DTX1** (AT2G04040, +3.43-fold), a MATE-type transporter involved in detoxification — both consistent with a broad stress defense response.

### 4.3 Advantages and Limitations of ONT Direct RNA-seq

This study demonstrates the feasibility of ONT direct RNA-seq for plant transcriptomics, achieving mean read lengths of ~1,000 nt with N50 values exceeding 1,100 nt. The FLAIR-based isoform analysis identified 20,958 isoforms, enabling expression quantification at the isoform rather than gene level — a resolution not achievable with standard short-read RNA-seq.

However, several limitations should be noted:

1. **Unbalanced design:** The exclusion of two *anac017-1* × AA replicates (requiring GPU basecalling) precluded estimation of the genotype × treatment interaction. An additive model was used instead, which assumes independent effects of genotype and treatment without synergy or antagonism.

2. **ELIGOS2 failure:** The epitranscriptomic modification analysis could not be completed due to software compatibility issues in the Docker environment. This remains an important avenue for future analysis, as direct RNA-seq is uniquely suited for detecting RNA modifications without chemical treatment.

3. **Lower throughput per sample:** ONT DRS typically yields 0.3–1.7 million reads per flow cell, compared to tens of millions with Illumina. While sufficient for differential expression at the gene/isoform level, some low-abundance transcripts may be underrepresented.

4. **Single replicate in one group:** The *anac017-1* × AA condition had only one replicate, reducing statistical power for this specific comparison and preventing interaction term estimation.

### 4.4 K-CHOPORE Pipeline Performance

The K-CHOPORE pipeline successfully automated the complete analysis workflow from raw FASTQ merging through differential expression, processing 10 samples through 78+ Snakemake jobs. The pipeline's Docker containerization ensures full reproducibility, while the Snakemake workflow provides automatic dependency resolution, parallel execution, and restart capability (via `--keep-going`). The modular design allows selective execution of pipeline stages and accommodation of varying data availability (FASTQ-only vs. FAST5 signal data).

### 4.5 Conclusions

Using ONT direct RNA sequencing and the K-CHOPORE pipeline, we identified 435 differentially expressed isoforms by genotype and 266 by treatment in a study of mitochondrial retrograde signaling in *Arabidopsis*. The recovery of AOX1A as a top Antimycin A-responsive gene validates both the experimental system and the analytical approach. The strong genotype effect along PC1 (44% variance) confirms ANAC017 as a major transcriptional regulator. Future work should include: (1) completing the full factorial design with all 12 samples, (2) resolving ELIGOS2 compatibility issues for epitranscriptomic analysis, and (3) performing m6Anet-based m6A detection when FAST5 storage becomes available.

---

## References

1. De Clercq, I., et al. (2013). The membrane-bound NAC transcription factor ANAC013 functions in mitochondrial retrograde regulation of the oxidative stress response in Arabidopsis. *Plant Cell*, 25(9), 3472-3490.
2. Ng, S., et al. (2013). A membrane-bound NAC transcription factor, ANAC017, mediates mitochondrial retrograde signaling in Arabidopsis. *Plant Cell*, 25(9), 3450-3471.
3. Van Aken, O., et al. (2016). Retrograde signalling caused by heritable mitochondrial dysfunction is partially mediated by ANAC017 and improves plant performance. *Plant J*, 88(4), 542-558.
4. Tang, A.D., et al. (2020). Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. *Nat Commun*, 11(1), 1438.
5. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100.
6. Love, M.I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*, 15(12), 550.

---

*Pipeline: K-CHOPORE v1.0 (https://github.com/biopelayo/K-CHOPORE)*
*Analysis performed: February 2026*
