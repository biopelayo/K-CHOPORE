# K-CHOPORE v3.0 — Hands-On Tutorial

## Quick Start Guide for Workshop Participants

**Duration**: ~90 minutes
**Prerequisites**: Docker installed, Git, basic command line familiarity

---

## Table of Contents

1. [Setup and Installation](#1-setup-and-installation)
2. [Understanding the Pipeline Architecture](#2-understanding-the-pipeline-architecture)
3. [Configuration Walkthrough](#3-configuration-walkthrough)
4. [Exercise 1: Running the Core Pipeline (Dry Run)](#exercise-1-running-the-core-pipeline-dry-run)
5. [Exercise 2: Exploring the Test Suite](#exercise-2-exploring-the-test-suite)
6. [Exercise 3: lncRNA Discovery Module](#exercise-3-lncrna-discovery-module)
7. [Exercise 4: Epitranscriptomic Consensus](#exercise-4-epitranscriptomic-consensus)
8. [Exercise 5: Interpreting Results](#exercise-5-interpreting-results)
9. [Advanced: Adding Your Own Data](#advanced-adding-your-own-data)
10. [Troubleshooting FAQ](#troubleshooting-faq)

---

## 1. Setup and Installation

### 1.1 Clone the Repository

```bash
git clone https://github.com/biopelayo/K-CHOPORE.git
cd K-CHOPORE
```

### 1.2 Verify Prerequisites

```bash
# Check Docker
docker --version
# Expected: Docker version 24.x or newer

# Check Python
python3 --version
# Expected: Python 3.8+

# Check pip packages
pip install snakemake pyyaml pytest
```

### 1.3 Build the Docker Image (Optional for Full Runs)

```bash
# Full build takes ~30 min
docker build -t k-chopore:v3 .
```

> **Workshop note**: For the hands-on exercises, you do NOT need to build the full Docker image. We will use the test suite and script-level exercises that run natively with Python.

### 1.4 Install Python Dependencies for Testing

```bash
pip install pytest pyyaml
```

---

## 2. Understanding the Pipeline Architecture

K-CHOPORE v3.0 has **46 Snakemake rules** organized into a core pipeline and 5 optional modules:

```
Core Pipeline (22 rules):
  Basecalling → Filtering → QC → Alignment → Isoforms → Modifications → DESeq2 → MultiQC

Optional Modules (24 rules):
  M1: lncRNA Discovery      (7 rules)  ← Uses existing DRS data
  M2: Small RNA Analysis     (5 rules)  ← Requires Illumina sRNA-seq
  M3: miRNA Target Prediction (4 rules) ← Requires Module 2 + degradome
  M4: Enhanced Epitranscriptomics (3 rules) ← Uses existing DRS data
  M5: Data Integration       (4 rules)  ← Requires all modules
```

### Key Design Principles

- **Modular**: Each module can be enabled/disabled independently
- **Containerized**: All tools in a single Docker image
- **Tested**: 136 automated tests validate every component
- **FAIR-compliant**: Reproducible, documented, version-controlled

---

## 3. Configuration Walkthrough

Open `config/config.yml` in your favorite editor:

```bash
# Use any text editor
nano config/config.yml
# or
code config/config.yml
```

### Key Sections to Explore

**1. Samples** (line ~5): How samples are defined with metadata:
```yaml
samples:
  WT_C_R1:
    genotype: "WT"
    treatment: "C"
    data_type: "fastq"
```

**2. Module toggles** (search for `run_lncrna`):
```yaml
params:
  run_lncrna: true           # Module 1: ON
  run_smallrna: false         # Module 2: OFF (no Illumina data yet)
  run_mirna_targets: false    # Module 3: OFF
  run_epitx_enhanced: true    # Module 4: ON
  run_integration: false      # Module 5: OFF
```

**3. Tool parameters**: Every tool has configurable thresholds:
```yaml
tools:
  minimap2_preset: "splice"
  minimap2_kmer_size: 14
  flair_support: 3
  eligos2_pval: 0.05
```

> **Exercise**: Find the CPAT threshold for lncRNA classification. What value is it set to? (Answer: 0.39)

---

## Exercise 1: Running the Core Pipeline (Dry Run)

A dry run shows what Snakemake *would* execute without actually running anything:

```bash
# If you have Docker + Snakemake installed:
docker run --rm -v "$(pwd)":/workspace -w /workspace k-chopore:v3 \
    snakemake -n --printshellcmds 2>&1 | head -50

# Without Docker (just to see the DAG):
snakemake -n --printshellcmds --configfile config/config.yml 2>&1 | head -50
```

### What to Look For

1. How many jobs does Snakemake plan to execute?
2. Which rules appear first? (Hint: they should be basecalling or index creation)
3. Can you spot the conditional module targets?

> **Key concept**: Snakemake works backwards from targets — it figures out what input files are needed, then what rules produce those inputs, recursively.

---

## Exercise 2: Exploring the Test Suite

### 2.1 Run All Tests

```bash
python -m pytest tests/ -v
```

**Expected output**: 136 passed in ~3 seconds

### 2.2 Understanding Test Categories

```bash
# Config validation tests
python -m pytest tests/test_config_validation.py -v

# Snakefile structure tests
python -m pytest tests/test_snakefile_rules.py -v

# v3 integration tests with toy data
python -m pytest tests/test_v3_integration.py -v
```

### 2.3 Inspect Toy Data

```bash
# Look at the toy genome (2 small chromosomes)
cat tests/toy_data/reference/toy_genome.fasta

# Look at the annotation (3 genes)
cat tests/toy_data/annotations/toy_annotation.gtf

# Look at mock FLAIR counts
cat tests/toy_data/mock_results/flair/counts_matrix.tsv
```

> **Discussion question**: Why is toy data important for testing bioinformatics pipelines? What are the advantages over using real data for tests?

---

## Exercise 3: lncRNA Discovery Module

This is the heart of Module 1. Let's walk through the consensus logic step by step.

### 3.1 Understanding the Consensus Approach

The lncRNA consensus requires agreement from multiple tools:

```
Candidate transcript
        |
        v
[TransDecoder] ──> ORF > 100aa? ──> YES ──> EXCLUDE (likely coding)
        |
        NO (small/no ORF)
        |
        v
[FEELnc] ──> non-coding? ──> Vote +1
[CPC2]   ──> non-coding? ──> Vote +1
[CPAT]   ──> prob < 0.39? ──> Vote +1
        |
        v
Votes >= 2/3? ──> YES ──> CONSENSUS lncRNA
                  NO  ──> EXCLUDE (low confidence)
```

### 3.2 Run the Consensus Script with Mock Data

```bash
# Run lncRNA consensus with toy data
python scripts/lncrna_consensus.py \
    --transdecoder tests/toy_data/mock_results/lncrna/transdecoder/longest_orfs.pep \
    --feelnc-gtf tests/toy_data/mock_results/lncrna/feelnc/candidate_lncRNA.gtf \
    --feelnc-class tests/toy_data/mock_results/lncrna/feelnc/candidate_lncRNA.gtf \
    --cpc2 tests/toy_data/mock_results/lncrna/cpc2/cpc2_results.txt \
    --cpat tests/toy_data/mock_results/lncrna/cpat/cpat_results.tsv \
    --candidates-gtf tests/toy_data/mock_results/lncrna/candidates_merged.gtf \
    --counts-matrix tests/toy_data/mock_results/flair/counts_matrix.tsv \
    --cantatadb tests/toy_data/annotations/cantata_test.bed \
    --te-bed tests/toy_data/annotations/te_test.bed \
    --min-consensus 2 \
    --max-orf-aa 100 \
    --cpat-threshold 0.39 \
    --out-gtf /tmp/lncrna_final.gtf \
    --out-bed /tmp/lncrna_final.bed \
    --out-counts /tmp/lncrna_counts.tsv \
    --out-report /tmp/lncrna_report.tsv
```

### 3.3 Examine the Output

```bash
echo "=== lncRNA Report ==="
cat /tmp/lncrna_report.tsv

echo "=== lncRNA GTF ==="
cat /tmp/lncrna_final.gtf

echo "=== lncRNA Counts ==="
cat /tmp/lncrna_counts.tsv
```

### 3.4 Discussion Questions

1. How many lncRNAs passed the consensus filter?
2. Which transcripts were excluded and why?
3. Look at the `consensus_tools` attribute in the GTF — which tools agreed?
4. What would happen if you changed `--min-consensus` to 3? To 1?
5. Why do we need multiple tools instead of just one?

> **Hands-on challenge**: Try running with `--min-consensus 3` and observe which lncRNAs are lost. Then try `--cpat-threshold 0.5` — does MSTRG.1.1 become non-coding?

---

## Exercise 4: Epitranscriptomic Consensus

### 4.1 Understanding Multi-Tool Modification Detection

Each tool uses a different approach:
- **ELIGOS2**: Base-level error patterns (purely computational)
- **m6Anet**: Deep learning on signal data (requires Nanopolish)
- **Nanocompore**: Signal-level k-mer comparisons (differential)

Sites detected by >= 2 tools are high-confidence modifications.

### 4.2 Run the Epitx Consensus

```bash
python scripts/epitx_consensus.py \
    --eligos-dir tests/toy_data/mock_results/eligos \
    --m6anet-dir tests/toy_data/mock_results/m6anet \
    --min-tools 1 \
    --pval 0.05 \
    --out-consensus /tmp/epitx_consensus.tsv \
    --out-biotype /tmp/epitx_biotype.tsv
```

### 4.3 Examine the Results

```bash
echo "=== Modification Consensus ==="
cat /tmp/epitx_consensus.tsv

echo "=== Biotype Stratification ==="
cat /tmp/epitx_biotype.tsv
```

### 4.4 Discussion Questions

1. How many sites were detected by each tool?
2. Which genomic positions appear in both ELIGOS2 and m6Anet?
3. What does the biotype column tell us about modification patterns?
4. Why is the m6Anet threshold different from the ELIGOS2 threshold?
5. How would Nanocompore add value when available?

---

## Exercise 5: Interpreting Results

### 5.1 Connecting the Dots

In a real experiment, the modules connect:

```
lncRNAs (M1) ──┬──> Modified lncRNAs (M4)
                ├──> cis-regulatory pairs (M5)
                └──> miRNA targets (M3)

miRNAs (M2) ───┬──> Target genes (M3)
                └──> Co-expression networks (M5)

Modifications ──> Modified vs unmodified genes ──> Functional enrichment
```

### 5.2 Mock Integration Analysis

Let's manually connect some results:

```bash
# Which lncRNAs are differentially expressed?
echo "=== DE lncRNAs ==="
cat tests/toy_data/mock_results/deseq2/deseq2_genotype_results.csv | grep "MSTRG"

# Cross-reference with lncRNA consensus
echo "=== lncRNA consensus ==="
cat /tmp/lncrna_report.tsv
```

**Questions**:
1. Which MSTRG transcripts are both lncRNAs AND differentially expressed?
2. What biological story could connect an anac017-1 mutant to lncRNA dysregulation?
3. How would modification data on these lncRNAs change the interpretation?

### 5.3 Visualizing the Workflow

Create a Snakemake DAG visualization:
```bash
# If Snakemake and graphviz are installed:
snakemake --rulegraph | dot -Tpng > dag.png
```

---

## Advanced: Adding Your Own Data

### Step 1: Prepare Config

Edit `config/config.yml`:

```yaml
samples:
  my_sample_1:
    genotype: "WT"
    treatment: "control"
    data_type: "fastq"
    nas_dirs: ["my_sample_1"]
    run_subdirs: ["run1"]
```

### Step 2: Place Your Data

```
data/
  raw/
    my_sample_1/
      run1/
        *.fastq.gz
  reference/
    genome/
      my_genome.fasta
    annotations/
      my_annotation.gtf
```

### Step 3: Download Annotations (Arabidopsis)

```bash
bash scripts/download_annotations.sh
```

### Step 4: Enable Modules

```yaml
params:
  run_lncrna: true
  run_epitx_enhanced: true
  # Enable others as data becomes available
```

### Step 5: Run

```bash
docker run --rm -v "$(pwd)":/workspace -w /workspace k-chopore:v3 \
    snakemake --cores 40 --latency-wait 60 --printshellcmds --keep-going
```

---

## Troubleshooting FAQ

### Q: Tests fail with `ModuleNotFoundError: yaml`
```bash
pip install pyyaml
```

### Q: Docker build fails on NTFS
Add `--latency-wait 60` to Snakemake commands. NTFS mount propagation can be slow.

### Q: ELIGOS2 gives `error testCMH`
Known rpy2/R compatibility issue. The pipeline creates placeholder files. This doesn't affect m6Anet or Nanocompore.

### Q: FLAIR crashes with "underscore in ID"
FLAIR quantify doesn't allow underscores in sample IDs. The pipeline auto-converts to hyphens in the manifest.

### Q: How do I run just one module?
Target specific outputs:
```bash
# Only lncRNA discovery
snakemake --cores 12 results/lncrna/lncrna_final.gtf

# Only epitx consensus
snakemake --cores 12 results/epitx/modification_consensus.tsv
```

### Q: What does each test file cover?

| Test File | Tests | Focus |
|-----------|-------|-------|
| `test_config_validation.py` | 31 | Config structure, required keys, value ranges |
| `test_lncrna_module.py` | 12 | lncRNA consensus script unit tests |
| `test_smallrna_module.py` | 17 | All v3 module rules + configs |
| `test_snakefile_rules.py` | 11 | Rule consistency, log directives |
| `test_snakefile_dag.py` | 14 | DAG structure, rule counts |
| `test_minimap2_fix.py` | 8 | Alignment parameter validation |
| `test_v3_integration.py` | 38 | End-to-end with toy data |
| **Total** | **136** | |

---

## Workshop Summary

### What You Learned
1. How K-CHOPORE organizes a complex multi-tool bioinformatics pipeline
2. The Snakemake + Docker + config-driven approach to reproducible workflows
3. Multi-tool consensus strategies for lncRNA discovery and modification detection
4. How automated tests validate pipeline correctness
5. How to configure and extend the pipeline for your own data

### Key Takeaways
- **Consensus > single tool**: Using >=2 tools reduces false positives
- **Modularity is key**: Enable only what you need, add data sources incrementally
- **Test-driven development**: 136 tests catch regressions before they reach production
- **FAIR principles**: Docker + config + Git = reproducible science

### Next Steps
1. Run Phase A (lncRNA + Epitx) with your DRS data
2. Generate Illumina sRNA-seq for Phase B
3. Enable Module 5 integration after all phases complete
4. Contribute to the project: https://github.com/biopelayo/K-CHOPORE

---

**Happy analyzing!** 🧬🥩 (*K-CHOPORE — as layered as a cachopo!*)
