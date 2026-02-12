# K-CHOPORE
### Keen Comprehensive High-throughput Omics Pipeline Organizer

Just like the iconic Asturian **cachopo**, K-CHOPORE is a layered and satisfying bioinformatics pipeline designed for analyzing **Nanopore direct RNA-seq** data with a focus on **epitranscriptomics**! Dive into the data feast, where each component works together to deliver high-quality, reproducible results.

---

![image](https://github.com/user-attachments/assets/64bb5e93-3e63-4352-99a2-bd1ba350a670)


## Overview
**K-CHOPORE** is an open-source pipeline for comprehensive analysis of Oxford Nanopore Technologies (ONT) **direct RNA sequencing** data. It handles every step from basecalling through epitranscriptomic modification detection and differential expression analysis.

The pipeline integrates **Snakemake**, **Docker**, **Python**, and **R** with established bioinformatics tools to provide a **FAIR-compliant** (Findable, Accessible, Interoperable, and Reusable) workflow.

### Key Features
- **FAIR-Compliant**: Adheres to FAIR principles for reproducible workflows
- **Automated Workflow**: Snakemake-driven with conditional module execution
- **Containerized Environment**: Docker image with all dependencies pre-installed
- **Direct RNA-seq Optimized**: Splice-aware alignment, RNA modification detection, isoform analysis
- **Modular Design**: Enable/disable pipeline stages via config flags
- **Multi-tool Epitranscriptomics**: ELIGOS2, m6Anet, and xPore for RNA modification detection
- **POD5 Support**: Compatible with both legacy FAST5 and modern POD5 formats

---

## Table of Contents
- [Pipeline Architecture](#pipeline-architecture)
- [Integrated Tools](#integrated-tools)
- [Getting Started](#getting-started)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
- [Pipeline Stages](#pipeline-stages)
- [Customization](#customization)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)
- [Contributing](#contributing)
- [Contact](#contact)

---

## Pipeline Architecture

```
FAST5/POD5 (raw signal)
    |
    v
[1. BASECALLING] ---- Dorado / Guppy / Bonito
    |
    v
FASTQ (reads)
    |
    v
[2. READ FILTERING] -- NanoFilt (quality + length filtering)
    |
    v
[3. READ QC] --------- NanoPlot (per-sample) + NanoComp (cross-sample)
    |
    v
[4. ALIGNMENT] ------- Minimap2 (splice-aware: -ax splice -uf -k14)
    |
    v
Sorted BAM
    |
    +---> [5. ALIGNMENT QC] --- samtools flagstat/stats + pycoQC
    |
    +---> [6. ISOFORM ANALYSIS]
    |         |--- FLAIR (align -> correct -> collapse -> quantify -> diffExp/diffSplice)
    |         |--- StringTie2 (long-read isoform assembly, -L flag)
    |
    +---> [7. EPITRANSCRIPTOMIC MODIFICATION DETECTION]
    |         |--- ELIGOS2 (base-level error signal analysis)
    |         |--- m6Anet (Nanopolish eventalign -> dataprep -> inference)
    |         |--- xPore (differential modification between conditions)
    |
    +---> [8. DIFFERENTIAL EXPRESSION] --- DESeq2 (from FLAIR counts)
    |
    v
[9. AGGREGATE QC] ---- MultiQC (combined report from all tools)
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
| | Nanopolish | Signal-level eventalign (required by m6Anet and xPore) |
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
git clone https://github.com/pelayovic/K-CHOPORE.git
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

### Samples
```yaml
samples:
  - "WT_C_R1"
  - "WT_C_R2"

conditions:
  WT_C_R1: "control"
  WT_C_R2: "control"
```

### Tool Settings (Direct RNA-seq optimized)
```yaml
tools:
  basecaller: "dorado"
  dorado_model: "rna004_130bps_hac@v5.0.0"     # Direct RNA model
  minimap2_preset: "splice"                       # Splice-aware alignment
  minimap2_kmer_size: 14                          # k=14 for RNA
  minimap2_extra_flags: "--secondary=no --MD"
  nanofilt_min_quality: 7
  nanofilt_min_length: 200
```

### Module Toggles
```yaml
params:
  run_basecalling: false    # Set true if starting from FAST5/POD5
  run_nanofilt: true        # Read quality filtering
  run_nanoplot: true        # Per-sample QC plots
  run_nanocomp: true        # Cross-sample comparison
  run_pycoqc: true          # pycoQC report
  run_flair: true           # Isoform analysis
  run_stringtie: false      # StringTie2 (alternative to FLAIR)
  run_eligos2: true         # ELIGOS2 modification detection
  run_m6anet: true          # m6Anet m6A detection
  run_xpore: false          # xPore differential modification
  run_deseq2: true          # DESeq2 differential expression
  run_multiqc: true         # Aggregate QC report
```

---

## Running the Pipeline

### Full pipeline via Docker
```bash
sudo docker run -it --rm \
    -v /path/to/your/local/data:/workspace \
    k-chopore \
    snakemake --snakefile /workspace/Snakefile \
    --configfile /workspace/config/config.yml \
    --cores 12 --latency-wait 30 --printshellcmds
```

### Run a specific rule
```bash
# Only alignment
snakemake --cores 12 results/sorted_bam/WT_C_R1_sorted.bam

# Only NanoPlot QC
snakemake --cores 4 results/nanoplot/WT_C_R1/NanoStats.txt

# Only m6Anet
snakemake --cores 12 results/m6anet/WT_C_R1/data.result.csv.gz
```

### Dry run (check what will be executed)
```bash
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
Statistical differential expression analysis from FLAIR counts matrix:
- Generates results CSV with log2FC, p-values, adjusted p-values
- MA plot, volcano plot, and PCA plot
- Handles single-condition gracefully (descriptive stats only)

### 9. MultiQC Aggregation
Combines outputs from NanoPlot, samtools, pycoQC, and other tools into a single interactive HTML report.

---

## Project Structure

```
K-CHOPORE/
├── config/
│   └── config.yml              # Main pipeline configuration
├── data/
│   ├── raw/
│   │   ├── fastq/              # Input FASTQ files
│   │   ├── fast5/              # Raw FAST5 files (legacy)
│   │   ├── pod5/               # Raw POD5 files (modern)
│   │   └── summaries/          # ONT sequencing summaries
│   └── reference/
│       ├── genome/             # Reference genome FASTA
│       ├── annotations/        # GTF/BED annotation files
│       └── transcriptome/      # Transcriptome FASTA
├── results/                    # All pipeline outputs
│   ├── basecalls/              # Basecalling output
│   ├── fastq_filtered/         # NanoFilt filtered reads
│   ├── nanoplot/               # NanoPlot QC reports
│   ├── nanocomp/               # NanoComp comparison
│   ├── mapped/                 # SAM alignment files
│   ├── sorted_bam/             # Sorted + indexed BAM files
│   ├── samtools_stats/         # Alignment statistics
│   ├── quality_analysis/       # pycoQC HTML reports
│   ├── flair/                  # FLAIR isoform results
│   ├── stringtie/              # StringTie2 isoform results
│   ├── eligos/                 # ELIGOS2 modification results
│   ├── nanopolish/             # Nanopolish eventalign output
│   ├── m6anet/                 # m6Anet m6A predictions
│   ├── xpore/                  # xPore differential modification
│   ├── deseq2/                 # DESeq2 DE results + plots
│   └── multiqc/                # MultiQC aggregate report
├── scripts/
│   ├── run_deseq2.R            # DESeq2 analysis script
│   ├── alignment/              # Alignment helper scripts
│   └── pipelines/              # Legacy pipeline scripts
├── Snakefile                   # Main Snakemake workflow
├── Dockerfile                  # Docker build file
├── requirements.txt            # Python dependencies
└── README.md
```

---

## Troubleshooting

- **Minimap2 alignment issues**: Ensure you use `-ax splice -uf` for direct RNA data (not `map-ont`). The `-uf` flag is critical as direct RNA reads are in the forward orientation.
- **m6Anet fails**: Requires Nanopolish eventalign with `--signal-index --scale-events` flags. FAST5 files must be accessible.
- **FLAIR correct fails**: Ensure GTF chromosome names match the reference genome FASTA headers.
- **DESeq2 single condition**: The script handles gracefully when only one condition is present (outputs descriptive stats).
- **Memory issues**: For large datasets, increase Docker memory allocation and thread count in `config.yml`.
- **Check logs**: All rules write detailed logs to the `logs/` directory.

---

## Citations

If you use K-CHOPORE, please cite the integrated tools:

- **Dorado/Guppy**: Oxford Nanopore Technologies basecaller software
- **NanoPlot/NanoComp/NanoFilt**: De Coster, W., et al. (2018). NanoPack. Bioinformatics, 34(15), 2666-2669
- **Minimap2**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100
- **pycoQC**: Leger, A., & Leonardi, T. (2019). pycoQC. Bioinformatics, 35(23), 5243-5245
- **FLAIR**: Tang, A. D., et al. (2020). Nature Communications, 11(1), 1438
- **StringTie2**: Kovaka, S., et al. (2019). Genome Biology, 20(1), 278
- **ELIGOS2**: Jenjaroenpun, P., et al. (2021). Nature Communications, 12(1), 7091
- **m6Anet**: Hendra, C., et al. (2022). Nature Methods, 19(12), 1590-1598
- **xPore**: Pratanwanich, P. N., et al. (2021). Nature Biotechnology, 39(11), 1394-1402
- **Nanopolish**: Simpson, J. T., et al. (2017). Nature Methods, 14, 407-410
- **DESeq2**: Love, M. I., et al. (2014). Genome Biology, 15(12), 550
- **MultiQC**: Ewels, P., et al. (2016). Bioinformatics, 32(19), 3047-3048
- **Picard**: Broad Institute. Picard toolkit. GitHub repository
- **samtools**: Danecek, P., et al. (2021). GigaScience, 10(2), giab008

---

## Limitations

### Basecalling Bottlenecks
Basecalling with Dorado/Guppy can be compute-intensive, especially without GPU acceleration. The Docker image includes CPU versions; GPU support requires NVIDIA Docker runtime.

### Signal-Level Analysis Requirements
m6Anet and xPore require raw FAST5 files and Nanopolish eventalign, which generates large intermediate files. Plan for sufficient disk space.

### Single-Condition Support
DESeq2 differential expression requires at least two conditions. The pipeline gracefully degrades to descriptive statistics with a single condition.

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
