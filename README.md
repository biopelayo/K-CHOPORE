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

## Project Structure

```
K-CHOPORE/
├── config/
│   └── config.yml              # Pipeline configuration (10 samples)
├── data/
│   ├── raw/                    # Input data (per-sample subdirectories)
│   └── reference/
│       ├── genome/             # TAIR10 reference genome FASTA + minimap2 index
│       ├── annotations/        # AtRTDv2 GTF annotation
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
│   └── multiqc/                # MultiQC aggregate report
├── scripts/
│   ├── run_deseq2.R            # DESeq2 factorial analysis (handles unbalanced designs)
│   ├── fix_eligos2.py          # ELIGOS2 compatibility fix (applied during Docker build)
│   ├── download_last_fast5.py  # FAST5 download utility (NAS -> local)
│   └── transfer_v4_twostep.py  # Two-step data transfer (NAS -> local -> server)
├── Snakefile                   # Main Snakemake workflow
├── Dockerfile                  # Docker build file (22.7 GB image with all tools)
├── run_pipeline.sh             # Pipeline launcher script
├── deploy_server.sh            # Server deployment script
├── requirements.txt            # Python dependencies
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
Basecalling with Guppy CPU is extremely slow (days per sample on a 40-thread server without GPU). The Docker image includes the CPU version; GPU support requires NVIDIA Docker runtime and significantly reduces processing time.

### Signal-Level Analysis Requirements
m6Anet and xPore require raw FAST5 files and Nanopolish eventalign, which generates large intermediate files. For this experiment, FAST5 data totals ~644 GB. Plan for sufficient local disk space or use SSHFS to mount remote storage.

### ELIGOS2 Known Issue
The ELIGOS2 CMH (Cochran-Mantel-Haenszel) statistical test fails in the current Docker environment due to an rpy2/R compatibility issue. The pipeline creates placeholder output files when this occurs. This is a known issue that requires rebuilding the Docker image with updated R/rpy2 packages.

### Unbalanced Factorial Design
The anac017-1 x AA group has only 1 replicate (R2 and R3 required GPU basecalling and were excluded). The DESeq2 R script detects this automatically and uses an additive model (`~ genotype + treatment`) instead of the full factorial, omitting the interaction term.

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
