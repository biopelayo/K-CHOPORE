# K-CHOPORE
### Keen Comprehensive High-throughput Omics Pipeline Organizer

Just like the iconic Asturian **cachopo**, K-CHOPORE is a layered and satisfying bioinformatics pipeline designed for analyzing **Nanopore direct RNA-seq** data with a focus on **epitranscriptomics** and **non-coding RNA biology**! Dive into the data feast, where each component works together to deliver high-quality, reproducible results.

---
<img width="1024" height="572" alt="image" src="https://github.com/user-attachments/assets/e66807c5-a076-4912-b40b-9fe9c7a1f327" />




## Overview
**K-CHOPORE** is an open-source pipeline for comprehensive analysis of Oxford Nanopore Technologies (ONT) **direct RNA sequencing** data. It handles every step from basecalling through epitranscriptomic modification detection, differential expression analysis, and multi-layered non-coding RNA characterization.

The pipeline integrates **Snakemake**, **Docker**, **Python**, and **R** with established bioinformatics tools to provide a **FAIR-compliant** (Findable, Accessible, Interoperable, and Reusable) workflow.

### Key Features
- **FAIR-Compliant**: Adheres to FAIR principles for reproducible workflows
- **Automated Workflow**: Snakemake-driven with conditional module execution
- **Containerized Environment**: Docker image with all dependencies pre-installed
- **Direct RNA-seq Optimized**: Splice-aware alignment, RNA modification detection, isoform analysis
- **Modular Design**: Enable/disable pipeline stages via config flags (`params.run_*` toggles)
- **Multi-tool Epitranscriptomics**: ELIGOS2, m6Anet, Nanocompore consensus with biotype stratification
- **lncRNA Discovery**: FEELnc + CPC2 + CPAT consensus classification with CANTATAdb cross-reference
- **Small RNA Analysis**: ShortStack + miRDeep-P2 for miRNA discovery from Illumina sRNA-seq
- **miRNA Target Prediction**: TargetFinder + psRNATarget + CleaveLand4 degradome validation
- **Multi-omics Integration**: WGCNA co-expression networks, cis-regulatory lncRNA detection, population variation context
- **POD5 Support**: Compatible with both legacy FAST5 and modern POD5 formats

---

### Current Experiment
- **Organism**: *Arabidopsis thaliana* (TAIR10 genome + AtRTDv2 annotation)
- **Sequencing**: ONT MinION, R9.4.1 flowcell, Direct RNA
- **Design**: 2x2 factorial — Genotype (WT vs *anac017-1* mutant) x Treatment (Control vs Antimycin A)
- **Samples**: 10 (WT×C: 3, WT×AA: 3, anac017-1×C: 3, anac017-1×AA: 1)
- **DESeq2 model**: Additive (`~ genotype + treatment`) due to unbalanced design
- **Results**:
  - 20,958 isoforms quantified across samples (FLAIR)
  - 435 differentially expressed genes by genotype (padj < 0.05, |LFC| > 1)
  - 266 differentially expressed genes by treatment (padj < 0.05, |LFC| > 1)

---

## Table of Contents
- [Pipeline Architecture](#pipeline-architecture)
- [v3.0 Module Overview](#v30-module-overview)
- [Integrated Tools](#integrated-tools)
- [Getting Started](#getting-started)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
- [Pipeline Stages](#pipeline-stages)
- [v3.0 Non-Coding RNA Modules](#v30-non-coding-rna-modules)
- [Customization](#customization)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)
- [Contributing](#contributing)
- [Contact](#contact)

---

## Pipeline Architecture

```
FAST5/POD5 (raw signal)                     Illumina sRNA-seq FASTQ
    |                                              |
    v                                              v
[1. BASECALLING] ---- Dorado / Guppy         [M2. SMALL RNA ANALYSIS]
    |                                              |--- Trim Galore (adapter trimming)
    v                                              |--- ShortStack (sRNA loci + miRNAs)
FASTQ (reads)                                      |--- miRDeep-P2 (novel miRNA prediction)
    |                                              |--- Annotation (miRBase + PmiREN)
    v                                              |
[2. READ FILTERING] -- NanoFilt                    v
    |                                         [M3. miRNA TARGET PREDICTION]
    v                                              |--- TargetFinder + psRNATarget
[3. READ QC] --------- NanoPlot + NanoComp         |--- CleaveLand4 (degradome validation)
    |                                              |--- DE anti-correlation analysis
    v                                              |
[4. ALIGNMENT] ------- Minimap2 (-ax splice)       |
    |                                              |
    v                                              |
Sorted BAM                                         |
    |                                              |
    +---> [5. ALIGNMENT QC] --- samtools + pycoQC  |
    |                                              |
    +---> [6. ISOFORM ANALYSIS]                    |
    |         |--- FLAIR (align/correct/collapse/quantify)
    |         |                                    |
    |         +---> [M1. lncRNA DISCOVERY]         |
    |                   |--- gffcompare (classification)
    |                   |--- TransDecoder (ORF filtering)
    |                   |--- FEELnc + CPC2 + CPAT (consensus)
    |                   |--- lncRNA DESeq2         |
    |                                              |
    +---> [7. EPITRANSCRIPTOMICS]                  |
    |         |--- ELIGOS2 + m6Anet + Nanocompore  |
    |         |--- [M4. ENHANCED EPITX]            |
    |                   |--- Multi-tool consensus  |
    |                   |--- Biotype stratification|
    |                                              |
    +---> [8. DIFFERENTIAL EXPRESSION] --- DESeq2  |
    |                                              |
    v                                              v
[9. AGGREGATE QC] ---- MultiQC        [M5. DATA INTEGRATION]
                                           |--- WGCNA co-expression networks
                                           |--- cis-regulatory lncRNA pairs
                                           |--- Population variation context
                                           |--- Integration report
```

---

## v3.0 Module Overview

K-CHOPORE v3.0 adds five new analysis modules for non-coding RNA biology, each independently toggleable via `config/config.yml`:

| Module | Toggle | Description | Data Required |
|--------|--------|-------------|---------------|
| **M1: lncRNA Discovery** | `params.run_lncrna` | Identifies long non-coding RNAs from FLAIR isoforms using a consensus of FEELnc, CPC2, and CPAT. Cross-references CANTATAdb for known lncRNAs. | DRS data (existing) |
| **M2: Small RNA Analysis** | `params.run_smallrna` | Discovers and annotates miRNAs from Illumina sRNA-seq using ShortStack and miRDeep-P2. | Illumina sRNA-seq FASTQ |
| **M3: miRNA Targets** | `params.run_mirna_targets` | Predicts miRNA targets with TargetFinder and psRNATarget, validates with CleaveLand4 degradome analysis. | Module 2 output + degradome FASTQ |
| **M4: Enhanced Epitx** | `params.run_epitx_enhanced` | Builds a multi-tool consensus of RNA modifications from ELIGOS2, m6Anet, and Nanocompore. Stratifies by biotype. | DRS data (existing) |
| **M5: Integration** | `params.run_integration` | WGCNA co-expression networks, cis-regulatory lncRNA detection, population variation context. | All module outputs |

### Recommended Execution Order

```
Phase A: Module 1 (lncRNA) + Module 4 (Epitx)    ← Uses existing DRS data
Phase B: Module 2 (sRNA) + Module 3 (Targets)     ← Requires Illumina sRNA-seq
Phase C: Module 5 (Integration)                     ← Requires all module outputs
```

---

## Integrated Tools

| Category | Tool | Purpose |
|----------|------|---------|
| **Basecalling** | Dorado | Modern ONT basecaller (recommended, supports direct RNA models) |
| | Guppy | Legacy ONT basecaller (CPU version included) |
| | Bonito | Research-grade ONT basecaller |
| **Read Filtering** | NanoFilt | Quality and length filtering for nanopore reads |
| **Read QC** | NanoPlot | Per-sample read quality visualization (length, quality, N50) |
| | NanoComp | Cross-sample comparative QC (violin plots, statistics) |
| **Alignment** | Minimap2 | Splice-aware long-read alignment (`-ax splice -uf -k14`) |
| **BAM Processing** | samtools | Sort, index, flagstat, stats for alignment QC |
| | Picard | BAM validation and genome normalization |
| **Isoform Analysis** | FLAIR | Full isoform pipeline: align, correct, collapse, quantify, diffExp, diffSplice |
| | StringTie2 | Long-read isoform assembly (with `-L` flag) |
| **Modification Detection** | ELIGOS2 | Epitranscriptomic signal from base-level error patterns |
| | m6Anet | m6A detection from Nanopolish eventalign signal data |
| | xPore | Differential RNA modification between conditions |
| | Nanocompore | Signal-level differential modification (k-mer level) |
| | Nanopolish | Signal-level eventalign (required by m6Anet, xPore, Nanocompore) |
| **lncRNA Discovery** | gffcompare | Isoform classification against reference annotation |
| | gffread | Sequence extraction from GTF/GFF coordinates |
| | TransDecoder | ORF prediction for coding potential filtering |
| | FEELnc | lncRNA classification (filter + codpot + classifier) |
| | CPC2 | Coding Potential Calculator 2 |
| | CPAT | Coding Potential Assessment Tool (plant model) |
| **Small RNA** | Trim Galore | Adapter trimming + quality filtering for Illumina sRNA-seq |
| | ShortStack | sRNA locus discovery + miRNA annotation (v4) |
| | miRDeep-P2 | Plant-specific novel miRNA prediction |
| **Target Prediction** | TargetFinder | Plant miRNA target prediction (complementarity-based) |
| | psRNATarget | miRNA target prediction (energy + accessibility model) |
| | CleaveLand4 | Degradome-seq validation of miRNA cleavage sites |
| **Integration** | WGCNA | Weighted Gene Co-expression Network Analysis (R) |
| | bcftools | Variant data processing for population context |
| | bedtools | Genomic interval operations (cis-regulation, TE intersection) |
| **Differential Expression** | DESeq2 | Differential gene/isoform expression with MA, volcano, PCA plots |
| **Aggregate QC** | MultiQC | Combined HTML report from all QC tools |
| **File Formats** | pod5 | Support for modern ONT POD5 format |
| | ont-fast5-api | Support for legacy FAST5 format |

---

## Getting Started

### Prerequisites

1. **Docker** (recommended for full reproducibility)
   ```bash
   docker --version
   ```

2. **Snakemake** (workflow management)
   ```bash
   pip install snakemake
   snakemake --version
   ```

3. **Git**
   ```bash
   git --version
   ```

4. **Python 3.8+**
   ```bash
   python --version
   ```

### Installation

**Step 1: Clone the repository**
```bash
git clone https://github.com/biopelayo/K-CHOPORE.git
cd K-CHOPORE
```

**Step 2: Build the Docker image**
```bash
sudo docker build -t k-chopore .
```

**Step 3: Configure the pipeline**

Edit `config/config.yml` to set your input paths, sample names, and enable/disable modules.

---

## Configuration

The main configuration file is `config/config.yml`. Key sections:

### Samples (2x2 Factorial Design)
Each sample has genotype, treatment, data type, and NAS directory metadata:
```yaml
samples:
  WT_C_R1:
    genotype: "WT"
    treatment: "C"
    data_type: "fastq"
    nas_dirs: ["WT_C_R1"]
    run_subdirs: ["no_sample/20220224_1639_MC-112869_FAR90122_f656439d"]
  WT_AA_R3:
    genotype: "WT"
    treatment: "AA"
    data_type: "fastq"
    nas_dirs: ["WT_AA_R3", "WT_AA_R3_2"]  # Multiple runs merged
    run_subdirs:
      - "no_sample/20220315_1614_MC-112869_FAR90098_c891954b"
      - "no_sample/20220316_1636_MC-112869_FAR90098_6e408c50"
```

**Note:** Sample IDs use underscores (`anac017_1`) while NAS folder names use hyphens (`anac017-1`). The pipeline converts underscores to hyphens automatically for FLAIR compatibility.

### Tool Settings (Direct RNA-seq optimized)
```yaml
tools:
  basecaller: "guppy"
  guppy_config_file: "/opt/ont-guppy-cpu/data/rna_r9.4.1_70bps_hac.cfg"
  minimap2_preset: "splice"                       # Splice-aware alignment
  minimap2_kmer_size: 14                          # k=14 for RNA
  minimap2_extra_flags: "--secondary=no --MD"
  nanofilt_min_quality: 7
  nanofilt_min_length: 200
  flair_support: 3
  eligos2_pval: 0.05
  eligos2_oddR: 2.5
  eligos2_esb: 0.1
  deseq2_padj_threshold: 0.05
  deseq2_lfc_threshold: 1.0
```

### Key Parameters
```yaml
params:
  threads: 40                 # Adjust for your server
  latency_wait: 60            # Helps with NFS/NTFS mount delays
  run_m6anet: false           # Requires FAST5 on local disk (~644 GB)

  # v3.0 Module toggles
  run_lncrna: true            # Module 1: lncRNA discovery
  run_smallrna: false          # Module 2: Small RNA analysis (Illumina)
  run_mirna_targets: false     # Module 3: miRNA target prediction
  run_epitx_enhanced: true     # Module 4: Enhanced epitranscriptomics
  run_integration: false       # Module 5: Data integration
```

### Annotation Bundle (v3.0)
Module-specific reference annotations are configured in the `annotations:` section:
```yaml
annotations:
  araport11_gff: "data/reference/annotations/Araport11_GFF3_genes_transposons.gff"
  cantatadb_bed: "data/reference/annotations/CANTATAdb_v3_arabidopsis.bed"
  mirbase_gff: "data/reference/annotations/miRBase_v22_ath.gff3"
  mirbase_mature_fa: "data/reference/annotations/miRBase_v22_ath_mature.fa"
  pmiren_fa: "data/reference/annotations/PmiREN_ath_miRNA.fa"
  te_annotation: "data/reference/annotations/TAIR10_TE_annotation.bed"
```

Download all annotations automatically:
```bash
bash scripts/download_annotations.sh
```

---

## Running the Pipeline

### Full pipeline via Docker
```bash
docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:latest \
    snakemake --cores 40 --latency-wait 60 --printshellcmds --keep-going
```

### Run a specific rule
```bash
# Only alignment for one sample
docker run --rm -v "$(pwd)":/workspace -w /workspace k-chopore:latest \
    snakemake --cores 12 results/sorted_bam/WT_C_R1_sorted.bam

# Only DESeq2 analysis
docker run --rm -v "$(pwd)":/workspace -w /workspace k-chopore:latest \
    Rscript scripts/run_deseq2.R results/flair/counts_matrix.tsv \
    results/deseq2/sample_sheet.tsv results/deseq2 0.05 1.0

# Only m6Anet for one sample (requires FAST5 on disk)
docker run --rm -v "$(pwd)":/workspace -w /workspace k-chopore:latest \
    snakemake --cores 12 results/m6anet/WT_C_R1/data.site_proba.csv
```

### Dry run (check what will be executed)
```bash
docker run --rm -v "$(pwd)":/workspace -w /workspace k-chopore:latest \
    snakemake -n --printshellcmds
```

---

## Pipeline Stages

### 1. Basecalling (Optional)
Converts raw FAST5/POD5 signal to FASTQ using **Dorado** (recommended) or **Guppy**.
- Dorado uses modern RNA models (`rna004_130bps_hac`)
- Supports both FAST5 and POD5 input formats

### 2. Read Filtering (NanoFilt)
Filters reads by quality score and length to remove low-quality data before alignment.

### 3. Read QC (NanoPlot + NanoComp)
- **NanoPlot**: Per-sample statistics and plots (read length distribution, quality scores, N50)
- **NanoComp**: Cross-sample comparison with violin plots for quality, length, and yield

### 4. Alignment (Minimap2, splice-aware)
Uses splice-aware alignment mode critical for direct RNA-seq:
- `-ax splice`: Splice-aware alignment preset
- `-uf`: Forward-strand alignment (ONT direct RNA reads are forward-stranded)
- `-k 14`: Kmer size optimized for RNA
- `--junc-bed`: Known junction annotation for guided alignment
- `--secondary=no`: Report only primary alignments
- `--MD`: Include MD tag for downstream modification detection

### 5. Alignment QC
- **samtools flagstat**: Mapping rates, duplicate counts, paired-end stats
- **samtools stats**: Comprehensive alignment statistics (insert size, coverage, etc.)
- **pycoQC**: Interactive HTML report from ONT sequencing summary

### 6. Isoform Analysis
**FLAIR** (primary):
- `flair align`: Initial alignment to reference
- `flair correct`: Correct splice sites using reference annotation
- `flair collapse`: Collapse reads into isoform models
- `flair quantify`: Quantify isoform expression across samples
- `flair diffExp`: Differential isoform expression between conditions
- `flair diffSplice`: Differential splicing events

**StringTie2** (alternative):
- Long-read mode (`-L` flag) for ONT isoform assembly
- Gene abundance estimation

### 7. Epitranscriptomic Modification Detection

**ELIGOS2**: Detects RNA modifications from base-level error patterns in direct RNA-seq data. Identifies m6A and other modifications by comparing error signals against an expected reference.

**m6Anet** (via Nanopolish):
1. `nanopolish index`: Links FAST5 signal data to FASTQ reads
2. `nanopolish eventalign`: Extracts signal-level alignment data
3. `m6anet dataprep`: Prepares eventalign data for inference
4. `m6anet inference`: Predicts m6A modification probability per DRACH site

**xPore**: Differential RNA modification detection between two conditions using nanopolish eventalign data. Identifies positions with significantly different modification rates.

### 8. Differential Expression (DESeq2)
Statistical differential expression analysis from FLAIR counts matrix using a **2x2 factorial design** (genotype x treatment):
- Full factorial model: `~ genotype + treatment + genotype:treatment` (when all groups have >=2 replicates)
- Additive model: `~ genotype + treatment` (automatically used when any group has <2 replicates)
- Extracts three contrasts: genotype effect, treatment effect, and interaction
- Generates MA plots, volcano plots, and PCA plot (colored by genotype, shaped by treatment)
- Reference levels: WT (genotype), C (treatment)

### 9. MultiQC Aggregation
Combines outputs from NanoPlot, samtools, pycoQC, and other tools into a single interactive HTML report.

---

## v3.0 Non-Coding RNA Modules

### Module 1: lncRNA Discovery (`run_lncrna: true`)

Identifies long non-coding RNAs from FLAIR-collapsed isoforms using a multi-tool consensus approach:

1. **`lncrna_filter_candidates`**: Merges all FLAIR isoform GTFs, classifies with gffcompare against Araport11, filters by length (>200 nt) and exon count, excludes known protein-coding transcripts.
2. **`lncrna_transdecoder`**: Predicts ORFs with TransDecoder. Transcripts with ORFs >100 aa are flagged as likely coding.
3. **`lncrna_feelnc`**: Three-step FEELnc analysis (filter → codpot → classifier) for lncRNA classification (lincRNA, antisense, intronic).
4. **`lncrna_cpc2`**: Coding Potential Calculator 2 scoring.
5. **`lncrna_cpat`**: Coding Potential Assessment Tool with Arabidopsis-trained model (threshold 0.39).
6. **`lncrna_consensus`**: Consensus classification requiring ≥2 of 3 tools (FEELnc + CPC2 + CPAT) to agree "non-coding" AND no large ORF. Cross-references CANTATAdb to mark known vs novel lncRNAs. Classifies TE-derived lncRNAs via bedtools intersection with TAIR10 TE annotation.
7. **`lncrna_deseq2`**: Differential expression analysis of lncRNAs using the same DESeq2 R script as the main pipeline.

**Key outputs:**
- `results/lncrna/lncrna_final.gtf` — Consensus lncRNA annotation
- `results/lncrna/lncrna_summary_report.tsv` — Classification and statistics
- `results/lncrna/deseq2/` — Differential lncRNA expression results

### Module 2: Small RNA Analysis (`run_smallrna: true`)

Processes Illumina small RNA-seq data to discover and annotate miRNAs:

1. **`srna_trim`**: Trim Galore adapter trimming with size selection (18-30 nt).
2. **`srna_shortstack`**: ShortStack v4 for sRNA locus discovery and MIRNA annotation from all samples jointly.
3. **`srna_mirdeep_p2`**: miRDeep-P2 for novel miRNA prediction with plant-specific criteria.
4. **`srna_annotate_known`**: Cross-references discovered miRNAs with miRBase v22 and PmiREN databases to classify known vs novel.
5. **`srna_counts_matrix`**: Generates miRNA count matrix and sRNA size distribution plots.

**Key outputs:**
- `results/smallrna/shortstack/Results.txt` — ShortStack locus table
- `results/smallrna/mirna_annotated.tsv` — Annotated miRNA catalog
- `results/smallrna/mirna_counts_matrix.tsv` — miRNA count matrix for downstream DE

### Module 3: miRNA Target Prediction (`run_mirna_targets: true`)

Predicts and validates miRNA-target interactions:

1. **`targets_prediction`**: TargetFinder for complementarity-based target prediction.
2. **`targets_psrnatarget`**: psRNATarget for energy + accessibility-based prediction.
3. **`degradome_validation`**: CleaveLand4 degradome-seq analysis to experimentally validate cleavage sites.
4. **`targets_integrate`**: Merges predictions from both tools with degradome evidence, cross-references with DESeq2 results to identify anti-correlated miRNA-target pairs, and flags lncRNA targets specifically.

**Key outputs:**
- `results/targets/target_evidence_table.tsv` — Comprehensive evidence table with prediction + validation scores
- `results/targets/mirna_target_network.tsv` — Network representation of miRNA-target interactions

### Module 4: Enhanced Epitranscriptomics (`run_epitx_enhanced: true`)

Builds a multi-tool consensus of RNA modifications:

1. **`nanocompore_run`**: Nanocompore signal-level differential modification analysis for each configured comparison (pairwise condition contrasts).
2. **`epitx_consensus`**: Parses and standardizes outputs from ELIGOS2, m6Anet, and Nanocompore into a common format. Identifies consensus sites (called by ≥2 tools). Stratifies modifications by biotype (mRNA vs lncRNA).
3. **`epitx_report`**: Generates an HTML summary report with modification counts per tool, biotype breakdown, and tool agreement statistics.

**Key outputs:**
- `results/epitx/modification_consensus.tsv` — High-confidence modification sites
- `results/epitx/modification_by_biotype.tsv` — Modifications stratified by transcript type
- `results/epitx/epitx_report.html` — Visual summary report

### Module 5: Data Integration (`run_integration: true`)

Integrates outputs from all modules into a systems-level view:

1. **`wgcna_network`**: WGCNA co-expression network analysis combining mRNA, lncRNA, and miRNA expression data. Auto-detects soft-thresholding power and correlates modules with genotype/treatment traits.
2. **`cis_regulation`**: Identifies lncRNAs within a configurable window (default 10 kb) of differentially expressed protein-coding genes. Classifies by orientation (sense/antisense) and position (upstream/downstream/overlapping).
3. **`population_context`**: Intersects ncRNA loci with 1001 Genomes Project SNP data to identify loci with high natural variation.
4. **`integration_report`**: Comprehensive HTML report combining WGCNA modules, cis-regulatory pairs, population variation, and miRNA-target networks.

**Key outputs:**
- `results/integration/wgcna_modules.tsv` — Co-expression network modules
- `results/integration/cis_lncrna_pairs.tsv` — Candidate cis-regulatory lncRNAs
- `results/integration/integration_report.html` — Final integrative summary

---

## Project Structure

```
K-CHOPORE/
├── config/
│   └── config.yml              # Pipeline configuration (samples + module toggles)
├── data/
│   ├── raw/                    # Input data (per-sample subdirectories)
│   │   ├── srna/               # Illumina sRNA-seq FASTQ (Module 2)
│   │   └── degradome/          # PARE/degradome FASTQ (Module 3)
│   └── reference/
│       ├── genome/             # TAIR10 reference genome FASTA + minimap2 index
│       ├── annotations/        # AtRTDv2 GTF + Araport11 + CANTATAdb + miRBase + PmiREN + TE
│       └── transcriptome/      # Transcriptome FASTA
├── results/                    # All pipeline outputs
│   ├── basecalls/              # Basecalling output (FASTQ + sequencing summaries)
│   ├── fastq_filtered/         # NanoFilt filtered reads
│   ├── nanoplot/               # NanoPlot QC reports (per sample)
│   ├── nanocomp/               # NanoComp cross-sample comparison
│   ├── sorted_bam/             # Sorted + indexed BAM files
│   ├── samtools_stats/         # Alignment statistics (flagstat + stats)
│   ├── quality_analysis/       # pycoQC HTML reports
│   ├── flair/                  # FLAIR isoform results + counts_matrix.tsv
│   ├── eligos/                 # ELIGOS2 modification results
│   ├── nanopolish/             # Nanopolish eventalign output
│   ├── m6anet/                 # m6Anet m6A predictions
│   ├── deseq2/                 # DESeq2 DE results + plots
│   ├── multiqc/                # MultiQC aggregate report
│   ├── lncrna/                 # [M1] lncRNA discovery outputs + lncRNA DESeq2
│   ├── smallrna/               # [M2] ShortStack + miRDeep-P2 + miRNA counts
│   ├── targets/                # [M3] Target prediction + degradome validation
│   ├── epitx/                  # [M4] Multi-tool modification consensus + report
│   └── integration/            # [M5] WGCNA + cis-regulation + integration report
├── scripts/
│   ├── run_deseq2.R            # DESeq2 factorial analysis (reused by lncRNA module)
│   ├── fix_eligos2.py          # ELIGOS2 compatibility fix (applied during Docker build)
│   ├── download_last_fast5.py  # FAST5 download utility (NAS -> local)
│   ├── transfer_v4_twostep.py  # Two-step data transfer (NAS -> local -> server)
│   ├── download_annotations.sh # [v3] Annotation bundle download + versioning
│   ├── lncrna_consensus.py     # [M1] Multi-tool lncRNA consensus classification
│   ├── annotate_mirna.py       # [M2] miRNA annotation (miRBase + PmiREN cross-ref)
│   ├── srna_counts.py          # [M2] sRNA count matrix + size distribution plots
│   ├── predict_targets.py      # [M3] TargetFinder wrapper for plant miRNA targets
│   ├── integrate_targets.py    # [M3] Target evidence integration + DE cross-ref
│   ├── epitx_consensus.py      # [M4] Multi-tool modification consensus builder
│   ├── epitx_report.py         # [M4] Epitranscriptomic summary HTML report
│   ├── run_wgcna.R             # [M5] WGCNA co-expression network analysis
│   ├── cis_regulation.py       # [M5] cis-regulatory lncRNA pair detection
│   └── integration_report.py   # [M5] Final integrative HTML report
├── tests/
│   ├── conftest.py             # Shared pytest fixtures
│   ├── test_config.yml         # Test configuration (toy data)
│   ├── test_config_validation.py  # Config structure validation tests
│   ├── test_snakefile_rules.py    # Rule presence + structure tests
│   ├── test_snakefile_dag.py      # DAG validation + rule count tests
│   ├── test_lncrna_module.py      # [M1] lncRNA consensus script tests
│   ├── test_smallrna_module.py    # [M2-M5] Module rule + script tests
│   └── toy_data/               # Minimal test data (FASTQ, FASTA, GTF, TSV)
├── Snakefile                   # Main Snakemake workflow (46 rules)
├── Dockerfile                  # Docker build file (all tools pre-installed)
├── run_pipeline.sh             # Pipeline launcher script
├── deploy_server.sh            # Server deployment script
├── requirements.txt            # Python dependencies
├── CHANGELOG.md                # Version history
└── README.md
```

---

## Troubleshooting

- **Minimap2 alignment issues**: Ensure you use `-ax splice -uf` for direct RNA data (not `map-ont`). The `-uf` flag is critical as direct RNA reads are in the forward orientation.
- **FLAIR underscore restriction**: FLAIR quantify does NOT allow underscores in sample ID, condition, or batch fields of the reads manifest. The pipeline automatically converts underscores to hyphens in the manifest while keeping file paths unchanged.
- **FLAIR correct fails**: Ensure GTF chromosome names match the reference genome FASTA headers. The pipeline strips `Chr` prefixes from FLAIR BED files to match TAIR10 numeric chromosome names.
- **ELIGOS2 CMH test failure**: Known rpy2/R compatibility issue in the Docker image causes ELIGOS2 to fail with `error testCMH` on all samples. The pipeline handles this gracefully by creating placeholder output files. This requires investigation of the R/rpy2 bridge in the Docker environment.
- **DESeq2 unbalanced design**: When any group in the 2x2 factorial has fewer than 2 replicates, the R script automatically falls back to an additive model (`~ genotype + treatment`), omitting the interaction term.
- **m6Anet disk requirements**: m6Anet requires FAST5 files on local disk (~644 GB for this experiment). Set `run_m6anet: false` unless FAST5 data fits on local storage.
- **Docker on NTFS**: Use `--latency-wait 60` for Snakemake when running Docker volumes on NTFS. Directory creation can be slow.
- **Guppy CPU basecalling**: Very slow without GPU (days per sample). Consider using `--keep-going` so other pipeline steps can proceed while basecalling runs.
- **Module 1 (lncRNA) no candidates**: If gffcompare finds no novel transcripts, check that FLAIR collapse ran successfully and produced isoform GTFs. Ensure Araport11 GFF is downloaded (`bash scripts/download_annotations.sh`).
- **Module 2 (sRNA) ShortStack slow**: ShortStack processes all samples jointly. For initial testing, use a subset via the `srna_samples` config section. Ensure trimmed reads are 18-30 nt.
- **Module 4 (Epitx) empty consensus**: The consensus requires ≥2 tools to agree. Check that at least ELIGOS2 and one other tool (m6Anet or Nanocompore) have run successfully.
- **Annotation downloads fail**: Some databases may change URLs. Check `scripts/download_annotations.sh` for current URLs and update if needed.
- **Check logs**: All rules write detailed logs to the `logs/` directory.

---

## Citations

If you use K-CHOPORE, please cite the integrated tools:

### Core Pipeline
- **Dorado/Guppy**: Oxford Nanopore Technologies basecaller software
- **NanoPlot/NanoComp/NanoFilt**: De Coster, W., et al. (2018). NanoPack. Bioinformatics, 34(15), 2666-2669
- **Minimap2**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100
- **pycoQC**: Leger, A., & Leonardi, T. (2019). pycoQC. Bioinformatics, 35(23), 5243-5245
- **FLAIR**: Tang, A. D., et al. (2020). Nature Communications, 11(1), 1438
- **StringTie2**: Kovaka, S., et al. (2019). Genome Biology, 20(1), 278
- **DESeq2**: Love, M. I., et al. (2014). Genome Biology, 15(12), 550
- **MultiQC**: Ewels, P., et al. (2016). Bioinformatics, 32(19), 3047-3048
- **samtools**: Danecek, P., et al. (2021). GigaScience, 10(2), giab008

### Epitranscriptomics
- **ELIGOS2**: Jenjaroenpun, P., et al. (2021). Nature Communications, 12(1), 7091
- **m6Anet**: Hendra, C., et al. (2022). Nature Methods, 19(12), 1590-1598
- **xPore**: Pratanwanich, P. N., et al. (2021). Nature Biotechnology, 39(11), 1394-1402
- **Nanocompore**: Leger, A., et al. (2021). Nature Communications, 12(1), 7025
- **Nanopolish**: Simpson, J. T., et al. (2017). Nature Methods, 14, 407-410

### lncRNA Discovery (Module 1)
- **FEELnc**: Wucher, V., et al. (2017). Nucleic Acids Research, 45(5), e21
- **CPC2**: Kang, Y. J., et al. (2017). Nucleic Acids Research, 45(W1), W12-W16
- **CPAT**: Wang, L., et al. (2013). Nucleic Acids Research, 41(6), e74
- **TransDecoder**: Haas, B., et al. TransDecoder. GitHub repository
- **gffcompare**: Pertea, G., & Pertea, M. (2020). F1000Research, 9, 304

### Small RNA Analysis (Module 2)
- **ShortStack**: Axtell, M. J. (2013). RNA, 19(6), 740-751
- **miRDeep-P2**: Kuang, Z., et al. (2019). Bioinformatics, 35(14), 2521-2523
- **Trim Galore**: Felix Krueger. Trim Galore. GitHub repository

### miRNA Target Prediction (Module 3)
- **TargetFinder**: Fahlgren, N., & Carrington, J. C. (2010). Methods in Molecular Biology, 592, 51-57
- **psRNATarget**: Dai, X., et al. (2018). Nucleic Acids Research, 46(W1), W49-W54
- **CleaveLand4**: Addo-Quaye, C., et al. (2009). Bioinformatics, 25(1), 130-131

### Data Integration (Module 5)
- **WGCNA**: Langfelder, P., & Horvath, S. (2008). BMC Bioinformatics, 9, 559

---

## Limitations

### Basecalling Bottlenecks
Basecalling with Guppy CPU is extremely slow (days per sample on a 40-thread server without GPU). The Docker image includes the CPU version; GPU support requires NVIDIA Docker runtime and significantly reduces processing time.

### Signal-Level Analysis Requirements
m6Anet, xPore, and Nanocompore require raw FAST5 files and Nanopolish eventalign, which generates large intermediate files. For this experiment, FAST5 data totals ~644 GB. Plan for sufficient local disk space or use SSHFS to mount remote storage.

### ELIGOS2 Known Issue
The ELIGOS2 CMH (Cochran-Mantel-Haenszel) statistical test fails in the current Docker environment due to an rpy2/R compatibility issue. The pipeline creates placeholder output files when this occurs. This is a known issue that requires rebuilding the Docker image with updated R/rpy2 packages.

### Unbalanced Factorial Design
The anac017-1 x AA group has only 1 replicate (R2 and R3 required GPU basecalling and were excluded). The DESeq2 R script detects this automatically and uses an additive model (`~ genotype + treatment`) instead of the full factorial, omitting the interaction term.

### Module Dependencies
Module 3 (miRNA Targets) requires Module 2 (Small RNA) outputs. Module 5 (Integration) requires all other modules. Enable modules in the recommended phase order (A → B → C) to ensure all dependencies are satisfied.

### Annotation Database Availability
The annotation download script (`scripts/download_annotations.sh`) fetches databases from external servers (TAIR, miRBase, PmiREN, CANTATAdb). Some URLs may change over time. Check the script if downloads fail and update URLs as needed.

---

## Contributing

We welcome contributions! To contribute:

1. Fork the repository
2. Create a feature branch
3. Make your changes and test them
4. Submit a pull request with a detailed description

---

## Contact and Support

- **Author**: Pelayo Gonzalez de Lena Rodriguez, MSc
- **Institution**: Cancer Epigenetics and Nanomedicine Lab | FINBA / Systems Biology Lab | University of Oviedo
- **Email**: pelayo.gonzalez@ispasturias.es
- **LinkedIn**: https://www.linkedin.com/in/biopelayo/
- **GitLab**: https://gitlab.com/bio.pelayo/

---

## Acknowledgments

We thank the Cancer Epigenetics and Nanomedicine Lab and the Systems Biology Lab at the University of Oviedo for their support.
