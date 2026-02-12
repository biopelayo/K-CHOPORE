# K-CHOPORE User Manual
## Keen Comprehensive High-throughput Omics Pipeline Organizer for Epitranscriptomics

**Version:** 2.0
**Author:** Pelayo Gonzalez de Lena Rodriguez, MSc
**Affiliation:** Cancer Epigenetics and Nanomedicine Lab | FINBA / Systems Biology Lab | University of Oviedo
**Purpose:** End-to-end analysis of Oxford Nanopore Technologies (ONT) Direct RNA Sequencing data with focus on epitranscriptomics, isoform discovery, and differential expression.

---

## Table of Contents

1. [What is K-CHOPORE?](#1-what-is-k-chopore)
2. [Prerequisites](#2-prerequisites)
3. [Installation](#3-installation)
4. [Understanding the Pipeline Architecture](#4-understanding-the-pipeline-architecture)
5. [Preparing Your Data](#5-preparing-your-data)
6. [Configuring the Pipeline](#6-configuring-the-pipeline)
7. [Running the Pipeline](#7-running-the-pipeline)
8. [Understanding Each Analysis Stage](#8-understanding-each-analysis-stage)
9. [Interpreting Results](#9-interpreting-results)
10. [Advanced Usage](#10-advanced-usage)
11. [Troubleshooting](#11-troubleshooting)
12. [Tool Reference](#12-tool-reference)

---

## 1. What is K-CHOPORE?

K-CHOPORE is a Dockerized Snakemake pipeline designed specifically for **ONT Direct RNA Sequencing** analysis. Unlike DNA or cDNA sequencing, direct RNA-seq reads the native RNA molecule, preserving chemical modifications (m6A, pseudouridine, etc.) in the raw electrical signal.

The pipeline covers the full analysis lifecycle:

```
Raw Signal (FAST5/POD5)
    |
    v
[Basecalling] --> FASTQ reads
    |
    v
[Quality Control] --> Filtered reads + QC reports
    |
    v
[Splice-aware Alignment] --> Sorted BAM files + alignment stats
    |
    +---> [Isoform Analysis] --> Novel isoforms, transcript quantification
    |
    +---> [Modification Detection] --> m6A sites, base modifications
    |
    +---> [Differential Expression] --> DE genes, volcano plots
    |
    v
[Aggregate QC Report] --> MultiQC HTML report
```

### What makes K-CHOPORE unique for direct RNA?

- Uses **splice-aware alignment** (`minimap2 -ax splice -uf -k14`) -- the `-uf` flag is critical because direct RNA reads are forward-stranded
- Integrates **signal-level modification detection** (m6Anet, ELIGOS2) that exploit the electrical signature of modified bases
- Provides **complete FLAIR isoform pipeline** (align, correct, collapse, quantify, diffExp, diffSplice) tailored for long reads
- All tools containerized in Docker for perfect reproducibility

---

## 2. Prerequisites

### Hardware Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU | 4 cores | 12+ cores |
| RAM | 16 GB | 32+ GB |
| Storage | 50 GB free | 200+ GB free |
| GPU | Not required | NVIDIA GPU for Dorado GPU basecalling |

### Software Requirements

| Software | Version | How to Install |
|----------|---------|----------------|
| Docker Desktop | 20.10+ | https://docs.docker.com/get-docker/ |
| Git | 2.30+ | https://git-scm.com/downloads |

That's it. Everything else is inside the Docker container.

### Verify Docker is Running

```bash
docker --version
docker ps
```

If Docker Desktop is installed but the daemon is not running, start it from your system tray or:

```bash
# Windows (PowerShell)
Start-Process "C:\Program Files\Docker\Docker\Docker Desktop.exe"

# macOS
open -a Docker

# Linux
sudo systemctl start docker
```

---

## 3. Installation

### Step 3.1 -- Clone the Repository

```bash
git clone https://gitlab.com/bio.pelayo/K-CHOPORE.git
cd K-CHOPORE
```

Or if you downloaded the ZIP:

```bash
unzip K-CHOPORE-master.zip
cd K-CHOPORE-master
```

### Step 3.2 -- Build the Docker Image

This downloads and compiles all bioinformatics tools. The first build takes **30-60 minutes** depending on your internet speed (the image is ~22 GB).

```bash
docker build -t k-chopore .
```

You will see verification output at the end:

```
=== K-CHOPORE Build Verification ===
Minimap2: 2.22-r1101
Samtools: samtools 1.19
Dorado: 0.8.0+acec121
StringTie: 2.2.1
NanoPlot: NanoPlot 1.46.2
NanoFilt: NanoFilt 2.8.0
m6anet: Version: 2.0.1
MultiQC: multiqc, version 1.33
Snakemake: 7.32.4
R: R version 4.1.2
=== K-CHOPORE Docker image built successfully! ===
```

### Step 3.3 -- Verify the Build

```bash
# Check the image exists
docker images k-chopore

# List all available pipeline rules
docker run --rm k-chopore snakemake --list
```

You should see 27 rules listed (from `all` to `xpore_diffmod`).

---

## 4. Understanding the Pipeline Architecture

### Module Toggle System

K-CHOPORE uses a **modular toggle** system. Each analysis stage can be independently enabled or disabled in `config/config.yml`:

```yaml
params:
  run_basecalling: false    # Only if starting from raw signal (FAST5/POD5)
  run_nanofilt: true        # Read quality filtering
  run_nanoplot: true        # Per-sample QC visualizations
  run_nanocomp: true        # Cross-sample comparison
  run_pycoqc: true          # Sequencing run QC (needs sequencing summary)
  run_flair: true           # Full isoform analysis pipeline
  run_stringtie: false      # Alternative isoform assembler
  run_eligos2: true         # Error-based modification detection
  run_m6anet: true          # m6A detection from nanopore signal
  run_xpore: false          # Differential modification (needs 2+ conditions)
  run_deseq2: true          # Differential gene expression
  run_multiqc: true         # Aggregate QC report
```

Set any module to `false` to skip it. The pipeline automatically adjusts its execution graph.

### Data Flow Between Stages

```
data/raw/fastq/{sample}.fastq
        |
        v
  [nanofilt] --> results/fastq_filtered/{sample}_filtered.fastq
        |
   +----|----+
   |    |    |
   v    v    v
[nanoplot] [nanocomp]  [minimap2 alignment]
                              |
                        results/mapped/{sample}.sam
                              |
                        [samtools sort/index]
                              |
                        results/sorted_bam/{sample}_sorted.bam
                              |
        +-----+-----+--------+--------+--------+
        |     |     |        |        |        |
        v     v     v        v        v        v
    [FLAIR] [StringTie] [ELIGOS2] [m6Anet] [pycoQC] [samtools stats]
        |                                              |
        v                                              v
  [FLAIR quantify]                               [MultiQC]
        |
   +----|----+
   |         |
   v         v
[DESeq2] [FLAIR diffExp/diffSplice]
```

---

## 5. Preparing Your Data

### Step 5.1 -- Organize Your Input Files

Create the following directory structure for your data:

```
my_project_data/
    raw/
        fastq/
            Sample1.fastq       # One FASTQ per sample
            Sample2.fastq
            Sample3.fastq
        fast5/
            Sample1/            # One folder per sample (for m6Anet)
                *.fast5
            Sample2/
                *.fast5
        pod5/                   # Modern ONT format (if using Dorado)
            Sample1/
                *.pod5
        summaries/
            Sample1_sequencing_summary.txt    # From MinKNOW
            Sample2_sequencing_summary.txt
        reads_manifest.tsv      # For FLAIR quantification
    reference/
        genome/
            genome.fasta        # Reference genome FASTA
        annotations/
            annotations.gtf     # Gene annotation GTF
        transcriptome/
            transcriptome.fa    # Transcriptome FASTA (optional, for xPore)
```

### Step 5.2 -- Input File Formats

**FASTQ files**: Standard FASTQ format from ONT basecalling. One file per sample. Can be gzipped.

```
@read_id_001 runid=abc123 ...
ACGUACGUACGUACGU...
+
!#$%&'()*+,-./01...
```

**Sequencing summary**: Tab-separated file produced by MinKNOW or Guppy. Required columns include `read_id`, `sequence_length_template`, `mean_qscore_template`. This file is needed for pycoQC.

**FLAIR reads manifest** (`reads_manifest.tsv`): Tab-separated file mapping samples to conditions for isoform quantification:

```
sample1_name    condition1    batch1    path/to/Sample1_filtered.fastq
sample2_name    condition1    batch1    path/to/Sample2_filtered.fastq
sample3_name    condition2    batch1    path/to/Sample3_filtered.fastq
```

**Reference genome**: FASTA format. For Arabidopsis, use TAIR10. For human, use GRCh38. Ensure chromosome names match the GTF annotation.

**GTF annotation**: Standard GTF/GFF format with gene and transcript features. This is used for splice junction guidance during alignment and isoform analysis.

### Step 5.3 -- Naming Convention

Sample names should follow a consistent pattern. The DESeq2 script automatically extracts conditions from sample names using the pattern:

```
{Condition}_{Replicate}
```

Examples:
- `WT_C_R1`, `WT_C_R2` --> condition = `WT_C`
- `MUT_T_R1`, `MUT_T_R2` --> condition = `MUT_T`
- `Control_R1`, `Control_R2` --> condition = `Control`

The last underscore-delimited segment is treated as the replicate identifier.

---

## 6. Configuring the Pipeline

### Step 6.1 -- Edit the Configuration File

Open `config/config.yml` in any text editor. This is the single file that controls the entire pipeline.

### Step 6.2 -- Define Your Samples

Replace the example samples with yours:

```yaml
samples:
  - "WT_R1"
  - "WT_R2"
  - "KO_R1"
  - "KO_R2"

conditions:
  WT_R1: "wildtype"
  WT_R2: "wildtype"
  KO_R1: "knockout"
  KO_R2: "knockout"
```

### Step 6.3 -- Set Input Paths

Point to your reference files. These paths are **relative to the workspace inside Docker** (`/workspace/`):

```yaml
input_files:
  fastq_dir: "data/raw/fastq"
  fast5_dir: "data/raw/fast5"
  pod5_dir: "data/raw/pod5"
  sequencing_summaries_dir: "data/raw/summaries"

  reference_genome: "data/reference/genome/GRCh38.fasta"
  reference_genome_mmi: "data/reference/genome/GRCh38.fasta.mmi"
  gtf_file: "data/reference/annotations/gencode.v44.annotation.gtf"

  reads_manifest: "data/raw/reads_manifest.tsv"
```

### Step 6.4 -- Configure Tool Parameters

Adjust parameters based on your experiment:

```yaml
tools:
  # Alignment (DO NOT change for direct RNA-seq)
  minimap2_preset: "splice"          # Always "splice" for RNA
  minimap2_kmer_size: 14             # k=14 optimal for RNA
  minimap2_extra_flags: "--secondary=no --MD"

  # Quality filtering -- adjust based on your data quality
  nanofilt_min_quality: 7            # Min Q-score (7-10 typical)
  nanofilt_min_length: 200           # Min read length in bp
  nanofilt_max_length: 0             # 0 = no maximum

  # Isoform analysis
  flair_support: 3                   # Min reads supporting an isoform

  # Modification detection
  eligos2_pval: 0.05                 # Significance threshold
  eligos2_oddR: 5                    # Odds ratio cutoff
  m6anet_num_iterations: 1000        # More = more accurate, slower

  # Differential expression
  deseq2_padj_threshold: 0.05       # Adjusted p-value cutoff
  deseq2_lfc_threshold: 1.0         # Log2 fold change cutoff
```

### Step 6.5 -- Enable/Disable Modules

Choose which analyses to run:

```yaml
params:
  threads: 12                        # Adjust to your CPU cores

  run_basecalling: false             # true only if starting from FAST5/POD5
  run_nanofilt: true
  run_nanoplot: true
  run_nanocomp: true
  run_pycoqc: true                   # Needs sequencing_summary files
  run_flair: true
  run_stringtie: false               # Alternative to FLAIR
  run_eligos2: true
  run_m6anet: true                   # Needs FAST5 files for signal data
  run_xpore: false                   # Needs 2+ conditions
  run_deseq2: true                   # Needs 2+ conditions
  run_multiqc: true
```

### Step 6.6 -- Common Configuration Scenarios

**Scenario A: Quick QC + alignment only (no modification detection)**
```yaml
run_nanofilt: true
run_nanoplot: true
run_flair: false
run_eligos2: false
run_m6anet: false
run_deseq2: false
```

**Scenario B: Full epitranscriptomic analysis**
```yaml
run_nanofilt: true
run_nanoplot: true
run_flair: true
run_eligos2: true
run_m6anet: true       # Requires FAST5 + Nanopolish
run_deseq2: true
run_multiqc: true
```

**Scenario C: Starting from raw signal (no pre-basecalled FASTQ)**
```yaml
run_basecalling: true   # Will use Dorado to produce FASTQ
# All other modules as desired
```

---

## 7. Running the Pipeline

### Step 7.1 -- Launch the Full Pipeline

Mount your data directory into the container and run Snakemake:

```bash
docker run --rm \
    -v /absolute/path/to/my_project_data:/workspace/data \
    -v /absolute/path/to/results:/workspace/results \
    k-chopore \
    snakemake --cores 12 --latency-wait 30
```

**Explanation of flags:**
- `--rm`: Remove the container after completion (image is preserved)
- `-v source:dest`: Mount host directories into the container
- `--cores 12`: Use 12 CPU threads (match your `params.threads`)
- `--latency-wait 30`: Wait 30 seconds for NFS/slow filesystems

### Step 7.2 -- Dry Run (Recommended First Step)

Before running the real analysis, always do a dry run to verify the pipeline will execute correctly:

```bash
docker run --rm \
    -v /absolute/path/to/my_project_data:/workspace/data \
    k-chopore \
    snakemake --cores 12 -n --printshellcmds
```

The `-n` flag means "dry run" -- it shows what *would* be executed without actually running anything. You should see a list of jobs and their shell commands.

### Step 7.3 -- Run Specific Modules Only

To run only specific targets:

```bash
# Only run QC (NanoFilt + NanoPlot)
docker run --rm -v ...:/workspace/data k-chopore \
    snakemake --cores 12 results/nanoplot/Sample1/NanoStats.txt

# Only run alignment for one sample
docker run --rm -v ...:/workspace/data k-chopore \
    snakemake --cores 12 results/sorted_bam/Sample1_sorted.bam

# Only run FLAIR isoform analysis
docker run --rm -v ...:/workspace/data k-chopore \
    snakemake --cores 12 results/flair/counts_matrix.tsv

# Only run DESeq2
docker run --rm -v ...:/workspace/data k-chopore \
    snakemake --cores 12 results/deseq2/deseq2_results.csv
```

### Step 7.4 -- Interactive Mode

For debugging or manual exploration:

```bash
docker run -it --rm \
    -v /absolute/path/to/my_project_data:/workspace/data \
    k-chopore \
    bash
```

This gives you a shell inside the container where you can run any tool manually:

```bash
# Inside the container:
minimap2 --version
samtools --version
NanoPlot --version
snakemake --list
```

### Step 7.5 -- Resume After Failure

If the pipeline fails partway through, Snakemake automatically resumes from where it stopped:

```bash
# Just run the same command again -- completed steps are skipped
docker run --rm -v ...:/workspace/data k-chopore \
    snakemake --cores 12 --latency-wait 30
```

To force re-execution of all steps:

```bash
docker run --rm -v ...:/workspace/data k-chopore \
    snakemake --cores 12 --forceall
```

---

## 8. Understanding Each Analysis Stage

### Stage 1: Basecalling (Optional)

**When to use:** Only if you have raw signal files (FAST5/POD5) and have not yet basecalled them.

**What it does:** Converts electrical signal to nucleotide sequences (FASTQ).

**Tools available:**
| Basecaller | When to Use |
|------------|-------------|
| Dorado | Modern data (POD5 format, R10.4.1 chemistry) |
| Guppy | Legacy data (FAST5 format, R9.4.1 chemistry) |
| Bonito | Research/alternative basecaller |

**Enable:** Set `run_basecalling: true` and configure:
```yaml
tools:
  basecaller: "dorado"
  dorado_model: "rna004_130bps_hac@v5.0.0"
```

**Output:** `results/basecalls/{sample}.fastq`

---

### Stage 2: Read Filtering (NanoFilt)

**What it does:** Removes low-quality reads and reads that are too short or too long.

**Why it matters:** ONT direct RNA reads can have variable quality. Filtering improves downstream alignment rates and reduces false positives in modification detection.

**Default thresholds:**
- Minimum quality score: Q7
- Minimum length: 200 bp
- Maximum length: unlimited

**Tuning guidance:**
- For high-depth experiments: increase to Q10 for cleaner data
- For low-depth experiments: reduce to Q5 to preserve more reads
- Direct RNA reads are typically 200-2000 bp; set max_length to 5000 to remove artifacts

**Input:** `data/raw/fastq/{sample}.fastq`
**Output:** `results/fastq_filtered/{sample}_filtered.fastq`

---

### Stage 3: Read Quality Control (NanoPlot + NanoComp)

**NanoPlot** generates per-sample QC visualizations:
- Read length distribution
- Quality score distribution
- Read length vs quality scatter plot
- Yield over time

**NanoComp** compares all samples side-by-side:
- Violin plots of read length
- Quality comparison
- Helps identify batch effects

**Output:**
- `results/nanoplot/{sample}/NanoStats.txt` -- Summary statistics
- `results/nanoplot/{sample}/` -- PNG/SVG/PDF plots
- `results/nanocomp/NanoComp-report.html` -- Interactive comparison report

---

### Stage 4: Splice-Aware Alignment (Minimap2)

**What it does:** Aligns filtered reads to the reference genome using splice-aware mode.

**Critical direct RNA flags:**
```
minimap2 -ax splice -uf -k14 --junc-bed annotation.gtf --secondary=no --MD
```

| Flag | Purpose |
|------|---------|
| `-ax splice` | Splice-aware alignment for RNA |
| `-uf` | Forward strand only (direct RNA reads are sense-strand) |
| `-k14` | k-mer size 14, optimized for RNA |
| `--junc-bed` | Use known splice junctions from GTF to improve alignment |
| `--secondary=no` | Report only primary alignments |
| `--MD` | Include MD tag for downstream modification tools (ELIGOS2) |

**Output:**
- `results/sorted_bam/{sample}_sorted.bam` -- Sorted, indexed BAM file
- `results/sorted_bam/{sample}_sorted.bam.bai` -- BAM index
- `results/samtools_stats/{sample}_flagstat.txt` -- Alignment summary
- `results/samtools_stats/{sample}_stats.txt` -- Detailed statistics

**Key metrics to check in flagstat:**
- Total reads mapped (should be >80% for good libraries)
- Supplementary alignments (chimeric reads)
- Properly paired (N/A for single-end direct RNA)

---

### Stage 5: Sequencing Run QC (pycoQC)

**What it does:** Generates an interactive HTML report from the MinKNOW sequencing summary, overlaid with alignment data.

**Requires:** Sequencing summary file from MinKNOW.

**Output:** `results/quality_analysis/pycoQC_output_{sample}.html`

**Report includes:**
- Pass/fail read counts
- Read length distribution over time
- Quality over time
- Channel activity heatmap
- Alignment rate per quality bin

---

### Stage 6: Isoform Analysis (FLAIR)

FLAIR runs a complete 6-step isoform analysis pipeline optimized for long reads:

**Step 6a -- FLAIR Align:** Aligns reads to genome and converts to BED12 format for splice junction identification.

**Step 6b -- FLAIR Correct:** Corrects splice sites using the reference GTF annotation. Misaligned splice junctions are snapped to the nearest annotated junction.

**Step 6c -- FLAIR Collapse:** Groups reads into isoform consensus sequences. Requires minimum support reads (default: 3). Produces:
- Isoform BED file
- Isoform FASTA sequences
- Isoform GTF annotation

**Step 6d -- FLAIR Quantify:** Counts reads per isoform per sample using the reads manifest. Produces a counts matrix (TSV).

**Step 6e -- FLAIR DiffExp:** Runs differential isoform expression between conditions (requires 2+ conditions).

**Step 6f -- FLAIR DiffSplice:** Identifies differential alternative splicing events (alternative 3' splice sites, alternative 5' splice sites, retained introns, etc.).

**Output:**
- `results/flair/{sample}_flair.collapse.isoforms.bed` -- Discovered isoforms
- `results/flair/{sample}_flair.collapse.isoforms.gtf` -- GTF for visualization
- `results/flair/counts_matrix.tsv` -- Isoform counts matrix
- `results/flair/diffExp/` -- Differential expression results
- `results/flair/diffSplice/` -- Differential splicing events

---

### Stage 7: RNA Modification Detection

K-CHOPORE offers three complementary modification detection approaches:

#### 7a -- ELIGOS2 (Error-based)

**Principle:** RNA modifications cause systematic basecalling errors. ELIGOS2 detects positions with significantly higher error rates than expected from sequence context alone.

**Strengths:** Does not require signal-level data; works from BAM files alone.

**Parameters:**
- `pval`: Significance threshold (default 0.05)
- `oddR`: Odds ratio cutoff (default 5)
- `esb`: Error of specific bases (default 0.2)

**Output:** `results/eligos/{sample}_eligos_output.txt`

#### 7b -- m6Anet (Signal-based m6A detection)

**Principle:** Uses the raw electrical signal from Nanopolish eventalign to detect m6A (N6-methyladenosine) at DRACH motifs using a trained neural network.

**Requirements:**
- FAST5 files (raw signal data)
- Nanopolish installed (for eventalign)

**Pipeline flow:**
1. `nanopolish index` -- Links FAST5 signals to reads
2. `nanopolish eventalign` -- Extracts per-position signal data
3. `m6anet dataprep` -- Formats signal data for inference
4. `m6anet inference` -- Predicts m6A probability per site

**Output:** `results/m6anet/{sample}/data.result.csv.gz` -- Per-site m6A probabilities

#### 7c -- xPore (Differential modification between conditions)

**Principle:** Compares nanopore signal distributions between conditions to identify differentially modified positions.

**Requirements:**
- Nanopolish eventalign data for all samples
- Two or more experimental conditions

**Output:** `results/xpore/diffmod.table`

---

### Stage 8: Differential Expression (DESeq2)

**What it does:** Tests for statistically significant differences in gene/isoform expression between conditions.

**Requirements:** At least 2 conditions with 2+ replicates each. Uses the FLAIR counts matrix as input.

**Output:**
- `results/deseq2/deseq2_results.csv` -- Full results table with log2FC, p-values
- `results/deseq2/MA_plot.pdf` -- Mean expression vs fold change
- `results/deseq2/volcano_plot.pdf` -- Significance vs fold change
- `results/deseq2/PCA_plot.pdf` -- Sample clustering
- `results/deseq2/normalized_counts.csv` -- Normalized expression values

**Single-condition mode:** If only one condition is present, DESeq2 outputs descriptive statistics (mean counts, standard deviation) instead.

---

### Stage 9: Aggregate QC Report (MultiQC)

**What it does:** Collects all QC outputs (NanoPlot stats, samtools stats, pycoQC) into a single interactive HTML report.

**Output:** `results/multiqc/multiqc_report.html`

Open this file in any web browser for a comprehensive overview of all samples.

---

## 9. Interpreting Results

### 9.1 -- Key Output Files Summary

| File | What It Contains |
|------|-----------------|
| `results/nanoplot/*/NanoStats.txt` | Read count, N50, mean quality per sample |
| `results/nanocomp/NanoComp-report.html` | Cross-sample comparison plots |
| `results/samtools_stats/*_flagstat.txt` | Alignment rate, mapped reads |
| `results/flair/*_flair.collapse.isoforms.bed` | Discovered transcript isoforms |
| `results/flair/counts_matrix.tsv` | Isoform expression matrix |
| `results/eligos/*_eligos_output.txt` | Predicted RNA modification sites |
| `results/m6anet/*/data.result.csv.gz` | m6A probabilities per DRACH motif |
| `results/deseq2/deseq2_results.csv` | Differentially expressed genes |
| `results/deseq2/volcano_plot.pdf` | Visual summary of DE analysis |
| `results/multiqc/multiqc_report.html` | Aggregate QC dashboard |

### 9.2 -- Quality Checkpoints

After the pipeline completes, verify these quality metrics:

**1. Read filtering (NanoFilt):**
- Check `results/nanoplot/*/NanoStats.txt`
- Expect >70% reads passing filter
- Mean quality should be Q8+

**2. Alignment rate:**
- Check `results/samtools_stats/*_flagstat.txt`
- Expect >80% mapping rate for good libraries
- If <60%, check: wrong reference genome? Low-quality reads? DNA contamination?

**3. Isoform discovery:**
- Check `results/flair/*_flair.collapse.isoforms.bed`
- Number of isoforms should be reasonable for your organism
- Compare with known annotation counts

**4. Differential expression:**
- Check `results/deseq2/PCA_plot.pdf`
- Replicates should cluster together
- Conditions should separate on PC1 or PC2

### 9.3 -- Interpreting m6Anet Results

The m6Anet output (`data.result.csv.gz`) contains columns:
- `transcript_id`: Transcript with the modification
- `transcript_position`: Position within the transcript
- `n_reads`: Number of reads covering this position
- `probability_modified`: Probability of m6A (0-1)
- `kmer`: The 5-mer sequence context

**Filtering recommendations:**
- `probability_modified > 0.9` for high-confidence m6A sites
- `n_reads >= 20` for reliable estimates
- Expect DRACH motif enrichment (D=A/G/U, R=A/G, H=A/C/U)

### 9.4 -- Interpreting DESeq2 Results

The DESeq2 output (`deseq2_results.csv`) contains columns:
- `gene`: Gene/isoform identifier
- `baseMean`: Average normalized count across samples
- `log2FoldChange`: Effect size (positive = upregulated in condition 2)
- `lfcSE`: Standard error of log2FC
- `stat`: Wald test statistic
- `pvalue`: Raw p-value
- `padj`: Benjamini-Hochberg adjusted p-value
- `significant`: "significant" or "not_significant" based on your thresholds

**Filtering for significant genes:**
- `padj < 0.05` AND `|log2FoldChange| > 1.0` (default thresholds)
- Adjust in config.yml if needed

---

## 10. Advanced Usage

### 10.1 -- Using a Custom Configuration File

```bash
docker run --rm \
    -v /path/to/data:/workspace/data \
    -v /path/to/my_custom_config.yml:/workspace/config/config.yml \
    k-chopore \
    snakemake --cores 12
```

### 10.2 -- Running on a Cluster (SLURM)

Create a Snakemake cluster profile:

```bash
docker run --rm \
    -v /path/to/data:/workspace/data \
    k-chopore \
    snakemake --cores 100 \
    --cluster "sbatch -p gpu -c {threads} --mem=32G -t 24:00:00" \
    --jobs 10
```

### 10.3 -- Using GPU for Basecalling

For Dorado GPU basecalling, pass through the GPU:

```bash
docker run --rm --gpus all \
    -v /path/to/data:/workspace/data \
    k-chopore \
    snakemake --cores 12
```

Make sure `nvidia-container-toolkit` is installed on the host.

### 10.4 -- Generating a Pipeline DAG Visualization

```bash
docker run --rm \
    -v /path/to/data:/workspace/data \
    -v /path/to/output:/workspace/dag_output \
    k-chopore \
    bash -c "snakemake --dag | dot -Tpdf > dag_output/pipeline_dag.pdf"
```

This requires graphviz (already in the container). Open the PDF to see the complete dependency graph.

### 10.5 -- Adjusting for Different Organisms

K-CHOPORE works with any organism. Adjust:

1. **Reference genome**: Download from Ensembl, NCBI, or TAIR
2. **GTF annotation**: Must match the reference genome chromosome names
3. **Dorado model**: Use the appropriate RNA model for your flowcell chemistry
4. **FLAIR support**: Reduce for low-coverage organisms, increase for high-coverage

### 10.6 -- Working with Multiple Flowcell Runs

If you sequenced the same sample across multiple flowcells, concatenate FASTQ files before running the pipeline:

```bash
cat run1/Sample1.fastq run2/Sample1.fastq > data/raw/fastq/Sample1.fastq
```

### 10.7 -- Downloading Bonito Models (Optional)

Bonito models are not included in the Docker image to save space. To download them:

```bash
docker run -v bonito_models:/root/.local/share/bonito k-chopore bonito download --models
```

The models are stored in a Docker volume and persist across container runs.

---

## 11. Troubleshooting

### Problem: "MissingInputException: Missing input files for rule index_genome"

**Cause:** The reference genome file is not found at the configured path.

**Solution:** Ensure your data directory is correctly mounted and the reference genome path in config.yml matches:
```bash
# Check files are visible inside the container
docker run --rm -v /path/to/data:/workspace/data k-chopore ls -la data/reference/genome/
```

### Problem: "No space left on device" during build

**Cause:** Docker image is ~22 GB; you need additional space for intermediate layers.

**Solution:**
```bash
# Clean Docker cache
docker system prune -a
# Ensure at least 50 GB free disk space
```

### Problem: NanoFilt removing too many reads

**Cause:** Quality threshold too strict for your data.

**Solution:** Lower the thresholds in config.yml:
```yaml
nanofilt_min_quality: 5    # Reduce from 7
nanofilt_min_length: 100   # Reduce from 200
```

### Problem: m6Anet fails with "Nanopolish not found"

**Cause:** Nanopolish compilation can fail during Docker build.

**Solution:** m6Anet requires signal-level data from Nanopolish eventalign. You can:
1. Install Nanopolish separately via conda on your host
2. Run Nanopolish eventalign outside Docker and provide the output
3. Use ELIGOS2 instead (works from BAM files, no signal data needed)

### Problem: DESeq2 outputs "single condition mode" warning

**Cause:** All your samples belong to the same condition. DESeq2 requires at least 2 conditions.

**Solution:** Add treatment/control samples, or update the `conditions` mapping in config.yml to correctly assign different conditions.

### Problem: FLAIR quantify fails

**Cause:** Missing or incorrectly formatted reads manifest.

**Solution:** Ensure `data/raw/reads_manifest.tsv` exists with the correct format:
```
SampleName\tCondition\tBatch\tPath/to/filtered.fastq
```
Paths should be relative to the workspace (`results/fastq_filtered/{sample}_filtered.fastq`).

### Problem: Pipeline runs too slowly

**Solutions:**
1. Increase threads: `params.threads: 24`
2. Run only necessary modules: disable unused modules
3. Use SSD storage for data/results
4. For basecalling: use Dorado with GPU (`--gpus all`)

---

## 12. Tool Reference

### Complete Tool Inventory

| Tool | Version | Category | Purpose |
|------|---------|----------|---------|
| Dorado | 0.8.0 | Basecalling | Modern ONT basecaller (GPU/CPU) |
| Guppy | 6.1.5 | Basecalling | Legacy ONT basecaller (CPU) |
| Bonito | 1.0.1 | Basecalling | Research basecaller |
| NanoFilt | 2.8.0 | QC | Read quality/length filtering |
| NanoPlot | 1.46.2 | QC | Read QC visualization |
| NanoComp | latest | QC | Cross-sample comparison |
| pycoQC | 2.5.2 | QC | Sequencing run QC |
| Minimap2 | 2.22 | Alignment | Splice-aware long-read aligner |
| Samtools | 1.19 | Alignment | BAM sort, index, stats |
| FLAIR | 2.0.0 | Isoforms | Full isoform analysis (6 steps) |
| StringTie2 | 2.2.1 | Isoforms | Alternative isoform assembler |
| ELIGOS2 | latest | Modifications | Error-based modification detection |
| m6Anet | 2.0.1 | Modifications | Signal-based m6A detection |
| xPore | latest | Modifications | Differential modification |
| Nanopolish | latest | Signal | Signal-level data extraction |
| DESeq2 | (R/Bioc) | Expression | Differential expression |
| MultiQC | 1.33 | Reporting | Aggregate QC reports |
| Snakemake | 7.32.4 | Workflow | Pipeline engine |

### Key Citations

If you use K-CHOPORE, please cite the relevant tools:

- **Minimap2:** Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100.
- **FLAIR:** Tang et al. (2020). Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. *Nature Communications*, 11, 1438.
- **ELIGOS2:** Jenjaroenpun et al. (2021). Decoding the epitranscriptional landscape from native RNA sequences. *Nucleic Acids Research*, 49(2), e7.
- **m6Anet:** Hendra et al. (2022). Detection of m6A from direct RNA sequencing using a multiple instance learning framework. *Nature Methods*, 19, 1590-1598.
- **xPore:** Pratanwanich et al. (2021). Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. *Nature Biotechnology*, 39, 1394-1402.
- **DESeq2:** Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
- **NanoPlot:** De Coster, W. et al. (2018). NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics*, 34(15), 2666-2669.
- **StringTie2:** Kovaka et al. (2019). Transcriptome assembly from long-read RNA-seq alignments with StringTie2. *Genome Biology*, 20, 278.
- **Snakemake:** Molder et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33.

---

## Quick Start Checklist

- [ ] Docker Desktop installed and running
- [ ] K-CHOPORE Docker image built (`docker build -t k-chopore .`)
- [ ] Reference genome FASTA prepared
- [ ] GTF annotation file prepared (matching chromosome names)
- [ ] FASTQ files named consistently (one per sample)
- [ ] `config/config.yml` edited with your sample names
- [ ] `config/config.yml` edited with correct file paths
- [ ] Module toggles set appropriately for your analysis
- [ ] `reads_manifest.tsv` created (if running FLAIR quantify)
- [ ] Sequencing summary files present (if running pycoQC)
- [ ] FAST5 directories organized per sample (if running m6Anet)
- [ ] Dry run successful (`snakemake -n`)
- [ ] Full pipeline launched (`snakemake --cores 12`)
- [ ] Results verified via MultiQC report

---

*K-CHOPORE -- Making ONT Direct RNA Epitranscriptomics Accessible*
