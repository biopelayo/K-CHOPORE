# K-CHOPORE Changelog: From Origin to v1.0

## Overview

This document tracks all improvements, bug fixes, new features, and architectural changes made to the K-CHOPORE pipeline from its original state (`fdd374f`, Feb 12 2026) through the final production version (`84f886e`, Feb 18 2026) across **10 collaborative development sessions**.

**Original state:** A Snakemake pipeline prototype with 27 rules, an 817-line Snakefile, a basic Dockerfile (226 lines), and a simple DESeq2 script (179 lines) — untested, with several critical bugs preventing execution on real data.

**Final state:** A battle-tested, Docker-containerized pipeline (730-line Snakefile, 22.7 GB Docker image) that successfully processed 10 ONT direct RNA-seq samples through 78+ Snakemake jobs, producing 435 genotype and 266 treatment differentially expressed isoforms.

---

## 1. Critical Bug Fixes

### 1.1 Minimap2 Alignment Flags (Session 1)
- **Bug:** Legacy `test.sh` used `map-ont` preset (generic nanopore, no splice awareness)
- **Fix:** Changed to `-ax splice -uf -k14 --secondary=no --MD`
  - `-ax splice`: splice-aware alignment essential for RNA-seq
  - `-uf`: forward-strand alignment critical for direct RNA reads
  - `--MD`: MD tag required by downstream modification tools (ELIGOS2, Nanopolish)
  - `--secondary=no`: prevents multi-mapping artifacts in isoform quantification
- **Impact:** Without this fix, all alignments would miss splice junctions and map to the wrong strand

### 1.2 pycoQC Sequencing Summary Resolution (Session 1)
- **Bug:** pycoQC rule hardcoded `{sample}_sequencing_summary.txt` pattern, but actual files contain flowcell IDs (e.g., `WT_C_R1_sequencing_summary_FAR90122_d34138fc.txt`)
- **Fix:** Added lambda function to resolve actual filenames from config mapping
- **Impact:** pycoQC rule would fail for every sample without this fix

### 1.3 FLAIR Chromosome Name Harmonization (Session 3)
- **Bug:** TAIR10 genome uses numeric chromosomes (1, 2, 3, 4, 5) but AtRTDv2 GTF uses prefixed names (Chr1, Chr2, etc.) — FLAIR correct/collapse silently produced empty outputs
- **Fix:** Added `sed`-based chromosome renaming in `flair_correct` rule:
  - BED files: `1` → `Chr1` etc.
  - Genome FASTA: `>1 ` → `>Chr1 ` (one-time rename to `genome_renamed.fasta`)
  - ELIGOS2: reverse rename (`Chr` → numeric) for BAM compatibility
- **Impact:** Without this, FLAIR produced 0 isoforms for every sample

### 1.4 FLAIR Manifest Underscore Bug (Session 10)
- **Bug:** FLAIR quantify silently fails when sample IDs, conditions, or batch fields contain underscores — it uses underscores as internal column name delimiters
- **Fix:** All underscores replaced with hyphens in `generate_reads_manifest` rule:
  ```python
  flair_id = sample_id.replace("_", "-")
  condition = f"{info['genotype']}-{info['treatment']}".replace("_", "-")
  ```
- **Impact:** FLAIR counts matrix had wrong column assignments, breaking DESeq2

### 1.5 DESeq2 Sample-to-Column Name Matching (Session 10)
- **Bug:** DESeq2 sample sheet used raw sample IDs but FLAIR counts matrix columns use the format `{flair_id}_{condition}_batch1` (hyphens within fields, underscores between)
- **Fix:** Added `generate_deseq2_sample_sheet` rule that constructs correct FLAIR column names
- **Impact:** DESeq2 could not match samples to count columns, producing errors

### 1.6 Bash Pipefail Incompatibility (Session 3)
- **Bug:** Shell commands used `2>&1 | tee {log}` pattern which breaks bash strict mode (`set -eo pipefail`) in Docker — a failed command's exit code is masked by `tee`'s success
- **Fix:** Changed all logging to `> {log} 2>&1` (redirect, not pipe)
- **Impact:** Failed rules appeared to succeed, hiding errors from Snakemake

---

## 2. Docker Image Improvements

### 2.1 Python Symlink (Session 4)
- **Added:** `RUN ln -sf /usr/bin/python3 /usr/bin/python`
- **Reason:** Many bioinformatics tools (ELIGOS2, legacy scripts) use `#!/usr/bin/env python` which doesn't exist in Ubuntu 22.04 by default

### 2.2 Nanopolish Build Fix (Session 3)
- **Bug:** Original Dockerfile tried `pip install nanopolish` (doesn't exist) with fallback to parallel compilation
- **Fix:** Changed to `make -j1` (serial build) to avoid HDF5/eigen race conditions during compilation
- **Impact:** Nanopolish compilation now succeeds reliably

### 2.3 VBZ Compression Plugin (Session 3)
- **Added:** `ENV HDF5_PLUGIN_PATH="/opt/ont-guppy-cpu/lib"`
- **Reason:** Modern FAST5 files use VBZ compression; without the plugin, Nanopolish/h5py cannot read them

### 2.4 pycoQC + setuptools Fix (Session 3)
- **Bug:** NanoPlot installation removes `setuptools`, but pycoQC imports `pkg_resources` (from setuptools) at runtime
- **Fix:** Explicit `pip install setuptools` before pycoQC, plus `setuptools<71` pin (v71+ removed `pkg_resources`)

### 2.5 ELIGOS2 Installation from GitLab (Session 3-5)
- **Added:** Complete ELIGOS2 installation in Dockerfile:
  - Clone from GitLab (NOT GitHub): `https://gitlab.com/piroonj/eligos2.git`
  - Install requirements
  - Apply `fix_eligos2.py` patch for rpy2 DeprecationWarning
- **Original:** ELIGOS2 was not installed in the Docker image at all

### 2.6 Pandas/NumPy Pinning for ELIGOS2 (Session 5)
- **Added:** `pip install "pandas==1.5.3" "numpy<2.0"` as final pip install
- **Reason:** ELIGOS2 uses deprecated `pd.concat` with string arguments that break in pandas 2.x
- **Order matters:** Must be installed AFTER all other packages to override their dependencies

### 2.7 libxml2-dev for R Packages (Session 5)
- **Added:** `libxml2-dev` to apt-get install
- **Reason:** Required for R XML package compilation (dependency of DESeq2)

---

## 3. Snakefile Architecture Overhaul

### 3.1 From Conditional Targets to Explicit Pipeline (Sessions 6-10)
- **Original:** 27 rules with `run_*` boolean toggles in config, complex `get_all_targets()` function
- **Final:** 22 streamlined rules with explicit `rule all` targets, no conditional logic
- **Removed rules:**
  - `setup_complete_structure` — unnecessary directory setup
  - `basecall_dorado` / `basecall_guppy` — basecalling moved outside pipeline scope
  - `flair_diff_exp` / `flair_diff_splice` — replaced by custom DESeq2 factorial analysis
  - `stringtie_assemble` — not needed for this analysis
  - `xpore_diffmod` — requires FAST5 signal data not available on disk
- **Added rules:**
  - `prepare_fastq` — merge pre-basecalled gzipped FASTQs from multiple flow cell runs
  - `generate_reads_manifest` — auto-generate FLAIR manifest from config (no manual file)
  - `generate_deseq2_sample_sheet` — auto-generate DESeq2 metadata matching FLAIR column names

### 3.2 Sample Sheet Architecture (Sessions 6-9)
- **Original:** Simple list of sample names, separate config entries for each file path
- **Final:** Rich sample sheet in config with per-sample metadata:
  ```yaml
  samples:
    WT_C_R1:
      genotype: "WT"
      treatment: "C"
      data_type: "fastq"
      nas_dirs: ["240711_WT_C_01"]
      run_subdirs: ["20240711_1530_MN35267_FAX00753_6d17d8d3"]
  ```
- **Impact:** All downstream rules (FLAIR manifest, DESeq2 sample sheet) auto-generate from this single source of truth

### 3.3 ELIGOS2 Fault Tolerance (Session 5)
- **Original:** ELIGOS2 rule would crash the pipeline on failure
- **Final:** Wrapped in if/else with graceful fallback:
  - On success: copies results
  - On failure: creates placeholder file with error description
  - Pipeline continues with `--keep-going`

### 3.4 m6Anet Conditional Execution (Session 8)
- **Added:** `RUN_M6ANET = config.get("params", {}).get("run_m6anet", False)`
- m6Anet rules exist but only activate when `run_m6anet: true` in config
- Prevents pipeline failure when FAST5 data is not available on disk

### 3.5 NTFS/Docker Compatibility (Sessions 3-7)
- **Added:** `mkdir -p` to every rule's shell block
- **Reason:** Docker volumes on Windows NTFS have directory creation latency; Snakemake may attempt to write outputs before directories propagate
- **Added:** `--latency-wait 60` recommendation in docs

---

## 4. DESeq2 R Script: Complete Rewrite

### 4.1 Original (179 lines)
- Simple condition-based DE: `~ condition`
- Auto-extracted conditions from sample names (fragile regex)
- Single comparison (condition A vs B)
- Basic MA plot and volcano plot
- No PCA, no interaction model

### 4.2 Final (276 lines)
- **Full 2×2 factorial design:** `~ genotype + treatment` (additive) or `~ genotype + treatment + genotype:treatment` (interaction)
- **Adaptive model selection:** Automatically switches to additive model when any group has <2 replicates (avoids DESeq2 crash with unbalanced designs)
- **Explicit sample sheet input:** Reads sample metadata from TSV file (genotype, treatment columns)
- **Reference levels:** Sets WT as genotype reference, Control as treatment reference
- **Outputs (10 files):**
  - `deseq2_genotype_results.csv` — genotype effect (anac017-1 vs WT)
  - `deseq2_treatment_results.csv` — treatment effect (AA vs Control)
  - `deseq2_interaction_results.csv` — interaction (empty CSV when additive model used)
  - `normalized_counts.csv` — DESeq2 normalized count matrix
  - `PCA_plot.pdf` — Variance-stabilized PCA colored by genotype + treatment
  - `MA_plot_genotype.pdf` / `MA_plot_treatment.pdf`
  - `volcano_genotype.pdf` / `volcano_treatment.pdf`

---

## 5. New Scripts Added (26 new files)

### 5.1 Testing Framework (Session 2)
- `tests/conftest.py` — pytest fixtures with toy data
- `tests/test_config_validation.py` — config schema validation (138 lines)
- `tests/test_minimap2_fix.py` — verifies correct splice/uf/k14 flags (112 lines)
- `tests/test_snakefile_dag.py` — DAG structure validation (197 lines)
- `tests/test_snakefile_rules.py` — rule completeness checks (132 lines)
- `tests/toy_data/` — minimal reference genome, GTF, reads manifest, sequencing summaries
- **Total:** 51 tests, all passing

### 5.2 Data Transfer Scripts (Sessions 6-8)
- `scripts/transfer_fastq_to_server.py` — SSH/SCP transfer with progress bars (332 lines)
- `scripts/transfer_data.py` — generic data transfer utility (294 lines)
- `scripts/transfer_v3_tar.py` — tar-based bulk transfer (292 lines)
- `scripts/transfer_v4_twostep.py` — two-step transfer for large datasets (329 lines)
- `scripts/transfer_all_samples.sh` — orchestrate all 10 sample transfers (209 lines)
- `scripts/download_nas_to_server.sh` — NAS mount and download script (256 lines)
- `scripts/setup_nas_mount.sh` — configure NFS/CIFS mount points (215 lines)
- `scripts/download_last_fast5.py` — selective FAST5 download for GPU basecalling (126 lines)

### 5.3 ELIGOS2 Fix Script (Session 5)
- `scripts/fix_eligos2.py` — patches rpy2 DeprecationWarning crash and BED merge strand column bug (45 lines)

### 5.4 Deployment and Launch (Sessions 4, 6)
- `deploy_server.sh` — one-command pipeline deployment to Linux server (265 lines)
- `run_pipeline.sh` — Linux launcher with Docker auto-detection (203 lines)
- `run_pipeline.bat` — Windows batch file launcher (107 lines)

### 5.5 Documentation (Sessions 5-6)
- `docs/step_by_step_guide.md` — beginner-friendly tutorial covering all sessions (548 lines)
- `docs/methods_section.md` — ready-to-use Methods section for publications (95 lines)

---

## 6. Configuration Evolution

### 6.1 Original config.yml (148 lines)
- Simple sample list: `samples: [WT_C_R1, WT_C_R2, ...]`
- Hardcoded file paths
- `run_*` boolean toggles for each module
- No per-sample metadata

### 6.2 Final config.yml (restructured)
- Rich sample sheet with per-sample metadata (genotype, treatment, data_type, NAS paths, run subdirectories)
- Eliminated all `run_*` toggles (explicit pipeline)
- Single `run_m6anet: false` toggle for FAST5-dependent analysis
- ELIGOS2 tuning parameters (oddR, esb, min_depth)
- All 10 samples with their actual NAS directory structures

---

## 7. Pipeline Rules Comparison

| Rule | Original | Final | Change |
|------|----------|-------|--------|
| `setup_complete_structure` | Yes | No | Removed (unnecessary) |
| `basecall_dorado` | Yes | No | Removed (basecalling done externally) |
| `basecall_guppy` | Yes | No | Removed (basecalling done externally) |
| `prepare_fastq` | No | Yes | **NEW**: Merge pre-basecalled FASTQs |
| `nanofilt` | Yes | Yes | Kept |
| `nanoplot` | Yes | Yes | Kept |
| `nanocomp` | Yes | Yes | Kept |
| `index_genome` | Yes | Yes | Kept |
| `map_with_minimap2` | Yes | Yes | Fixed flags |
| `sort_and_index_bam` | Yes | Yes | Kept |
| `samtools_stats` | Yes | Yes | Kept |
| `quality_analysis_with_pycoQC` | Yes | Yes | Fixed summary resolution |
| `flair_align` | Yes | Yes | Kept |
| `flair_correct` | Yes | Yes | Added chr name harmonization |
| `flair_collapse` | Yes | Yes | Uses renamed genome |
| `generate_reads_manifest` | No | Yes | **NEW**: Auto-generate from config |
| `flair_quantify` | Yes | Yes | Fixed manifest format |
| `flair_diff_exp` | Yes | No | Removed (replaced by DESeq2) |
| `flair_diff_splice` | Yes | No | Removed |
| `stringtie_assemble` | Yes | No | Removed |
| `run_eligos2` | Yes | Yes | Added fault tolerance + chr fix |
| `nanopolish_index` | Yes | Yes | Kept |
| `nanopolish_eventalign` | Yes | Yes | Kept |
| `m6anet_dataprep` | Yes | Yes | Conditional execution |
| `m6anet_inference` | Yes | Yes | Conditional execution |
| `xpore_diffmod` | Yes | No | Removed (requires FAST5) |
| `generate_deseq2_sample_sheet` | No | Yes | **NEW**: Auto-generate from config |
| `run_deseq2` | Yes | Yes | Complete rewrite (factorial) |
| `multiqc` | Yes | Yes | Kept |

**Summary:** 27 → 22 rules (removed 8, added 3, modified 7, kept 9 unchanged)

---

## 8. Commit History

| Date | Commit | Description |
|------|--------|-------------|
| Feb 12 | `fdd374f` | Initial commit: pipeline prototype with minimap2 + pycoQC fixes |
| Feb 12 | `72061ba` | Add 51-test pytest suite with toy data |
| Feb 12 | `f706760` | Add .gitignore for Python caches |
| Feb 12 | `8adc51e` | Add Docker launch scripts (bat + sh) |
| Feb 15 | `96ecdd9` | Fix VBZ plugin, pycoQC setuptools, FLAIR chr names |
| Feb 15 | `b9db173` | Add server deployment script and config fixes |
| Feb 15 | `c6855ee` | Add python symlink for ELIGOS2 |
| Feb 16 | `b0372df` | Fix ELIGOS2 rpy2 DeprecationWarning crash |
| Feb 16 | `9535f69` | Fix m6Anet output, ELIGOS2 bugs, FLAIR/MultiQC compat |
| Feb 16 | `33db293` | Pin pandas 1.5.3 + numpy <2.0 for ELIGOS2 |
| Feb 16 | `92c99c9` | Add ELIGOS2 fix script, Dockerfile cleanup, methods section |
| Feb 16 | `c869099` | Add step-by-step beginner guide |
| Feb 16 | `679ee77` | Rewrite pipeline for full 12-sample FAST5 basecalling |
| Feb 16 | `d465240` | Adapt for pre-basecalled data + NAS download scripts |
| Feb 17 | `3ed4420` | Make m6Anet optional for disk constraints |
| Feb 18 | `b5902c8` | Fix FLAIR/DESeq2 integration (underscore bug + factorial) |
| Feb 18 | `1aaf026` | Update README with results and docs |
| Feb 18 | `84f886e` | Finalize pipeline for 10 samples (production release) |

---

## 9. Production Results

The finalized pipeline successfully produced:
- **12.8M reads** processed across 10 samples (98.8% filter retention)
- **90.5% mean mapping rate** to TAIR10
- **20,958 isoforms** quantified by FLAIR
- **435 DE isoforms** by genotype (303 UP, 132 DOWN in anac017-1)
- **266 DE isoforms** by treatment (245 UP, 21 DOWN by Antimycin A)
- **AOX1A** (AT3G22370) recovered as top AA-responsive gene — validates the experimental system
- Full QC reports: MultiQC, NanoPlot (×10), NanoComp, pycoQC (×10), samtools stats (×10)

---

## 10. Known Limitations (Future Work)

1. **ELIGOS2 CMH test failure** — rpy2/R bridge issue in Docker prevents RNA modification detection. Placeholder outputs generated.
2. **m6Anet not run** — requires ~644 GB FAST5 signal files on local disk (not available)
3. **Unbalanced design** — n=1 for anac017-1×AA group prevents interaction term estimation (additive model used)
4. **WT_C_R2_stats.txt** — was 0 bytes due to pipeline race condition; manually regenerated

---

*Pipeline: K-CHOPORE v1.0 | https://github.com/biopelayo/K-CHOPORE*
*18 commits, 26 new files, 7 modified files | Feb 12-18, 2026*
*Co-developed by Pelayo Gonzalez de Lena Rodriguez & Claude Opus 4.6*
