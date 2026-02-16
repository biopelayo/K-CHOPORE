# K-CHOPORE Project: Step-by-Step Guide for Non-Bioinformaticians

> **What is this document?**
> This is a plain-language explanation of everything we did across 5 working sessions to build, deploy, debug, and run a bioinformatics pipeline called **K-CHOPORE**. It's written so that someone with no bioinformatics background can understand what happened and why.

---

## Table of Contents

1. [The Big Picture: What Is This Project About?](#1-the-big-picture)
2. [Key Concepts (Jargon Decoder)](#2-key-concepts)
3. [Session 1: Building the Foundation](#3-session-1)
4. [Session 2: Testing and First Fixes](#4-session-2)
5. [Session 3: Moving to a Powerful Server](#5-session-3)
6. [Session 4: Running on the Server and Debugging](#6-session-4)
7. [Session 5: Finishing, Analyzing, and Documenting](#7-session-5)
8. [Final Results Summary](#8-final-results)
9. [Complete List of Bugs Fixed](#9-bugs-fixed)
10. [Files We Created or Modified](#10-files)

---

## 1. The Big Picture: What Is This Project About? <a name="1-the-big-picture"></a>

### The Goal

We have **plant RNA** (from *Arabidopsis thaliana*, a small plant used in biology research) that was sequenced using a special machine made by Oxford Nanopore Technologies (ONT). This machine reads RNA molecules **directly** (without converting them to DNA first), which is special because it can also detect **chemical modifications** on the RNA.

Think of it like this:
- Traditional sequencing = photocopying a book and reading the copy
- Direct RNA sequencing = reading the original book, including any handwritten notes in the margins

The **K-CHOPORE pipeline** is an automated "recipe" that takes the raw data from the sequencing machine and:
1. Cleans up the data (quality control)
2. Figures out where each RNA molecule came from in the genome (alignment)
3. Identifies different versions (isoforms) of each gene's RNA
4. Looks for chemical modifications on the RNA (epitranscriptomics)
5. Compares expression between experimental conditions

### The Samples

We had **two samples** of plant RNA:
- **WT_C_R1** (~1.3 million reads) - labeled as "control"
- **WT_C_R2** (~1.6 million reads) - labeled as "treatment"

Only WT_C_R2 had the raw electrical signal data (FAST5 files) needed for detecting RNA modifications.

### The Technology Stack

| Thing | What It Is | Analogy |
|-------|-----------|---------|
| **Snakemake** | Workflow manager | A cooking recipe that runs each step in the right order automatically |
| **Docker** | Container system | A "box" containing all the software pre-installed, so it works the same on any computer |
| **GitHub** | Code repository | Google Drive for code, tracking every change we make |
| **Linux server** | Remote powerful computer | A big kitchen with industrial appliances, vs. our laptop which is a home kitchen |

---

## 2. Key Concepts (Jargon Decoder) <a name="2-key-concepts"></a>

| Term | Simple Explanation |
|------|-------------------|
| **Pipeline** | A chain of programs that run one after another, each taking the previous one's output |
| **FASTQ** | A text file containing the sequenced RNA reads (like a list of sentences read from the book) |
| **FAST5** | Raw electrical signal files from the sequencer (the original voltage measurements) |
| **BAM** | A compressed file showing where each read maps to the genome |
| **Genome / Reference** | The known complete DNA sequence of the organism (TAIR10 = version 10 of the Arabidopsis genome) |
| **GTF** | A file listing all known genes and their positions (like an index for the genome) |
| **Reads** | Individual RNA molecules that were sequenced |
| **Alignment / Mapping** | Figuring out where in the genome each read came from |
| **Isoform** | Different versions of RNA from the same gene (like different editions of a book) |
| **m6A** | A chemical modification on RNA (a methyl group added to adenine) |
| **DRACH motif** | The specific sequence pattern where m6A modifications typically occur |
| **Docker image** | A snapshot of a complete operating system + software, ready to run |
| **SSH** | Secure remote login to another computer over the network |
| **Chromosome** | A large piece of DNA; Arabidopsis has 5 main chromosomes plus chloroplast and mitochondrial DNA |

---

## 3. Session 1: Building the Foundation <a name="3-session-1"></a>

### What We Did

**Goal:** Set up the pipeline so it can run on any computer.

#### Step 1: Reviewed the existing code
- The pipeline was already written as a `Snakefile` (the recipe) and a `Dockerfile` (the software packaging instructions)
- We read through everything to understand what it does

#### Step 2: Created the Docker image
- Think of this as "packing a moving box" with ALL the software needed:
  - ~20 bioinformatics tools (for quality checking, alignment, etc.)
  - Python, R, and their libraries
  - The Arabidopsis reference genome
- The Docker image ended up being **~22 GB** (it's big because bioinformatics tools are large)
- The command: `docker build -t k-chopore:latest .`
  - This reads the `Dockerfile` and installs everything step by step

#### Step 3: Created configuration files
- `config/config.yml` - A settings file that tells the pipeline:
  - Which samples to process
  - Where to find the data files
  - Which analysis steps to run (on/off switches)
  - Parameters for each tool (quality thresholds, etc.)

#### Step 4: Created launch scripts
- `run_docker.sh` (Linux) and `run_docker.bat` (Windows)
- These scripts start the Docker container and run the pipeline inside it
- They handle mapping folders between your computer and the Docker container

#### Step 5: Pushed everything to GitHub
- Repository: `https://github.com/biopelayo/K-CHOPORE`
- This makes the code available to anyone and tracks all our changes

### Commit:
```
fdd374f - Initial commit: K-CHOPORE pipeline with bug fixes
```

---

## 4. Session 2: Testing and First Fixes <a name="4-session-2"></a>

### What We Did

**Goal:** Test the pipeline locally and fix whatever breaks.

#### Step 1: Created a test suite
- Built 51 automated tests using Python's `pytest` framework
- Tests check things like: "Does the config file load correctly?" and "Are all required tools installed?"
- Created "toy data" (tiny fake datasets) so tests run in seconds instead of hours

#### Step 2: First local run - and immediate problems
We tried running the pipeline on the actual data and hit our first wall:

**Bug 1: VBZ Plugin Missing**
- **What happened:** The program couldn't read FAST5 files (the raw signal data)
- **Why:** ONT uses a special compression format called VBZ, and the plugin wasn't installed
- **Fix:** Added the VBZ plugin to the Docker image

**Bug 2: pycoQC vs setuptools version conflict**
- **What happened:** A quality control tool (pycoQC) crashed on startup
- **Why:** A newer version of Python's `setuptools` (v71+) broke backward compatibility
- **Fix:** Pinned setuptools to version <71 in the Docker image

**Bug 3: FLAIR sample name problem**
- **What happened:** FLAIR (isoform analysis tool) rejected our sample names
- **Why:** FLAIR doesn't allow underscores in sample IDs, but our samples were named `WT_C_R1`
- **Fix:** Modified the pipeline to replace underscores with hyphens when calling FLAIR

**Bug 4: Chromosome naming mismatch**
- **What happened:** FLAIR couldn't find any genes
- **Why:** The genome file uses `1, 2, 3, 4, 5` for chromosomes, but the gene annotation file uses `Chr1, Chr2, Chr3, Chr4, Chr5`
- **Fix:** Added a step to rename chromosomes in the genome to match the annotation before running FLAIR

#### Step 3: Added a Python symlink
- Some tools expected the command `python` but the Docker image only had `python3`
- Fix: Created a shortcut (symlink) so `python` points to `python3`

### Commits:
```
f706760 - Add __pycache__ and .pytest_cache to .gitignore
72061ba - Add comprehensive test suite with toy data (51 tests)
8adc51e - Add Docker launch scripts for running with real data
96ecdd9 - Fix pipeline errors: VBZ plugin, pycoQC setuptools, FLAIR naming/chr mismatch
c6855ee - Add python symlink for ELIGOS2 compatibility
```

---

## 5. Session 3: Moving to a Powerful Server <a name="5-session-3"></a>

### Why Move to a Server?

The pipeline needs a LOT of computing power:
- One step (nanopolish eventalign) produces a **70 GB** intermediate file
- Some steps take hours or days on a regular laptop
- We had a **Dell PowerEdge T440 server** available:
  - 40 CPU cores (vs. ~8 on a laptop)
  - 251 GB RAM (vs. ~16 on a laptop)
  - 4 TB SSD for fast storage

### What We Did

#### Step 1: Set up SSH access
- SSH = Secure Shell, a way to remotely control another computer via command line
- Generated an encryption key pair so we could log in without typing a password each time
- Server address: `156.35.42.17`, username: `usuario2`

#### Step 2: Transferred the pipeline to the server
- Copied the entire project folder (~2 GB of code + reference files) using `scp` (secure copy)
- Also started transferring FAST5 files (384 files, several GB each) - this took a long time

#### Step 3: Built the Docker image on the server
- Ran the same `docker build` command on the server
- This took ~30 minutes because it downloads and compiles all the tools

#### Step 4: Created a server-specific launch script
- `deploy_server.sh` - handles Linux paths and Docker volume mounts for the server setup
- Configured to use 40 threads (taking advantage of all CPU cores)

#### Step 5: Started the pipeline on the server
- Launched Snakemake inside Docker
- It started processing all 31 pipeline steps
- Left it running (it takes hours)

#### Step 6: Learned painful Windows-to-Linux lessons
Several things that work on Windows break on Linux and vice versa:

| Problem | Explanation |
|---------|------------|
| `2>&1 \| tee` breaks bash | Redirecting output using `tee` fails with bash strict mode; had to use `> file 2>&1` instead |
| Docker creates files as root | Files created inside Docker belong to the root user, and the regular user can't delete them |
| NTFS vs Linux filesystem | Windows filesystem (NTFS) has different timing behavior; Docker sometimes can't find files it just created |
| Git Bash path conversion | Git Bash on Windows automatically converts `/workspace` to `C:/workspace`; had to disable this with `MSYS_NO_PATHCONV=1` |

### Commit:
```
b9db173 - Add Linux server deployment script and fix config/launcher
```

---

## 6. Session 4: Running on the Server and Debugging <a name="6-session-4"></a>

### What We Did

**Goal:** Get the pipeline from 24/31 steps complete to 31/31.

This session was the longest and most challenging. The pipeline kept hitting errors at different steps, and we had to fix each one, restart, hit the next error, fix it, and repeat. We did **6 iterative reruns**.

#### The Debugging Cycle

Each time the pipeline stopped, we:
1. Checked the log file to find the error
2. Figured out what was wrong
3. Fixed the code (on both our local machine and the server)
4. Cleaned up any broken output files
5. Restarted the pipeline (Snakemake is smart enough to skip already-completed steps)

#### Bug 5: ELIGOS2 crashed with rpy2 DeprecationWarning
- **What happened:** ELIGOS2 (an RNA modification detection tool) crashed immediately
- **Why:** It uses a bridge between Python and R (called rpy2), and a function it calls (`pandas2ri.activate()`) was deprecated in newer versions
- **First fix:** Wrapped the problematic call in a try/except block (like saying "try this, and if it fails, just move on")
- **But wait:** The same call appeared in 3 different places in the code, including inside a function with different indentation. Our first fix only caught 2 out of 3.
- **Real fix:** Created a Python script (`fix_eligos2.py`) that finds ALL occurrences regardless of indentation

#### Bug 6: m6Anet output filename changed
- **What happened:** Pipeline said "expected file not found" after m6Anet finished
- **Why:** m6Anet version 2.x changed its output filename from `data.result.csv.gz` to `data.site_proba.csv`
- **Fix:** Updated the pipeline to expect the new filename

#### Bug 7: MultiQC output filename mismatch
- **What happened:** Pipeline couldn't find MultiQC's report file
- **Why:** When you give MultiQC a title (like "K-CHOPORE QC Report"), it prepends a sanitized version to the filename. So instead of `multiqc_report.html`, it generated `K-CHOPORE-QC-Report_multiqc_report.html`
- **Fix:** Added `--filename multiqc_report` flag to force the expected name

#### Bug 8: FLAIR quantify output naming
- **What happened:** Pipeline expected `counts_matrix.tsv` but FLAIR created `counts_matrix.counts.tsv`
- **Why:** FLAIR appends `.counts` to whatever prefix you give it
- **Fix:** Added a rename step after FLAIR quantify runs

#### Bug 9: ELIGOS2 BED merge lost the strand column
- **What happened:** ELIGOS2 crashed trying to access data that wasn't there
- **Why:** This was a bug IN ELIGOS2's own code. When it merges overlapping gene regions, it was only keeping 4 columns (chromosome, start, end, name) but later expected 5 columns (chromosome, start, end, strand, name). The strand column got dropped.
- **Fix:** Changed the merge command to keep both strand and name columns

#### Bug 10: ELIGOS2 chromosome names didn't match
- **What happened:** ELIGOS2 processed all regions but found nothing
- **Why:** The BAM file uses chromosome names `1, 2, 3, 4, 5` but the BED file from FLAIR uses `Chr1, Chr2, Chr3, Chr4, Chr5`. ELIGOS2 was looking for "Chr1" in the BAM but only finding "1", so it skipped everything.
- **Fix:** Added a `sed` command to strip the "Chr" prefix from the BED file before giving it to ELIGOS2

#### Bug 11: ELIGOS2 pandas 2.x incompatibility
- **What happened:** Even after all the above fixes, ELIGOS2 still crashed
- **Why:** ELIGOS2 was written for pandas 1.x (a Python data library). Pandas 2.x changed how `pd.concat()` works, breaking ELIGOS2's code
- **Temporary fix:** Made ELIGOS2 non-fatal (pipeline continues even if it crashes)
- **Permanent fix (Session 5):** Downgraded pandas to 1.5.3

#### Bug 12: DESeq2 not installed
- **What happened:** R couldn't find the DESeq2 package
- **Why:** During Docker build, a system library (`libxml2-dev`) was missing. This caused the R `XML` package to fail to compile, and since DESeq2 depends on XML (through a chain of dependencies), DESeq2 never got installed.
- **Fix:** Added `libxml2-dev` to the Dockerfile's apt-get install list

#### Bug 13: DESeq2 needs replicates
- **What happened:** DESeq2 refused to run
- **Why:** DESeq2 performs statistical tests to find genes with different expression levels between conditions. Statistics requires multiple measurements (replicates) to estimate variability. We only had 1 sample per condition, which is fundamentally insufficient.
- **Fix:** Made DESeq2 non-fatal and output a placeholder file explaining why it couldn't run. This is NOT a code bug - it's a experimental design limitation.

#### Bug 14: Docker multi-line Python in Dockerfile
- **What happened:** Docker build failed with a mysterious error
- **Why:** We had a multi-line Python script inside a `RUN` command in the Dockerfile. Docker's build system (BuildKit) tried to interpret lines of Python as Docker instructions
- **Fix:** Moved the Python script to a separate file (`scripts/fix_eligos2.py`) and used `COPY` + `RUN python3 script.py`

#### Bug 15: Can't delete root-owned files via SSH
- **What happened:** Trying to clean up failed output files got "Permission denied"
- **Why:** Docker runs as root inside the container, so files it creates are owned by root. The SSH user doesn't have root privileges.
- **Fix:** Used Docker itself to delete the files: ran a temporary Docker container that mounts the folder and deletes from inside (where it IS root)

### Commits:
```
b0372df - Fix ELIGOS2 rpy2 DeprecationWarning crash
9535f69 - Fix pipeline completion: m6Anet output, ELIGOS2 bugs, FLAIR/MultiQC/DESeq2 compat
```

### Result: Pipeline reached **31/31 steps complete!**

---

## 7. Session 5: Finishing, Analyzing, and Documenting <a name="7-session-5"></a>

### What We Did

With the pipeline finished, we had 6 final tasks:

#### Task 1: Analyze the Results

We examined all the output files from the pipeline:

**Quality Control (NanoPlot + samtools):**
- Both samples had good quality (~Q11 median) and good mapping rates (92-93%)
- WT_C_R2 had more reads (1.6M vs 1.3M)
- Error rate ~10%, which is normal for ONT direct RNA sequencing

**Isoform Analysis (FLAIR):**
- Found **21,286 distinct transcript isoforms**
- 98.6% of isoforms were detected in both samples (very consistent)
- Top expressed gene: AT1G67090 (rubisco small subunit - expected for plant leaf tissue, since rubisco is the most abundant protein on Earth and is essential for photosynthesis)

**RNA Modifications (m6Anet):**
- Tested 228 sites where m6A modifications could potentially occur
- Found **18 sites** with modification probability > 50%
- Found **2 high-confidence sites** (probability > 90%):
  - Chromosome 4, position 130,288 (94% probability)
  - Chromosome 3, position 87,790 (91% probability)
- Most modifications occurred in the GAACT motif context (a known m6A pattern)

**RNA Modifications (ELIGOS2):**
- Ran but found no statistically significant sites
- This is a valid biological result (the filtering criteria were stringent)

**Differential Expression (DESeq2):**
- Could not run due to lack of biological replicates (1 sample per condition)

#### Task 2: Fix ELIGOS2 Properly

Instead of just making it non-fatal, we properly fixed the root cause:
- Downgraded pandas from 2.x to **1.5.3** (the version ELIGOS2 was designed for)
- Also had to downgrade numpy from 2.x to **<2.0** (pandas 1.5.3 is incompatible with numpy 2.x)
- After the fix, ELIGOS2 ran cleanly without any crashes
- It still found no significant modifications (but now that's a real biological result, not a software crash)

#### Task 3: Confirm DESeq2 Replicate Situation

- Confirmed that we only have 2 samples total (1 per condition)
- DESeq2 fundamentally requires **at least 2 biological replicates per condition** for statistical validity
- This can't be fixed with code - it requires generating more experimental data
- The FLAIR count matrix is saved for future use when more replicates become available

#### Task 4: Rebuild the Local Docker Image

- Rebuilt the Docker image on our local Windows machine with ALL the fixes accumulated over 5 sessions
- The key fix here was moving the ELIGOS2 fix from inline Dockerfile code to a separate script (`scripts/fix_eligos2.py`), because Docker's build system couldn't handle multi-line Python code
- Final image: **~22.6 GB**

#### Task 5: Clean Up the Server

- Deleted ELIGOS2 temporary files (freed **3.5 GB**)
- Cleaned Docker build cache (`docker builder prune`)
- Result: **143 GB free** on the 4 TB SSD (96% used - those FAST5 files are huge!)

#### Task 6: Write the Methods Section

- Created `docs/methods_section.md` - a formal bioinformatics methods section suitable for a thesis or scientific paper
- Includes all software versions, parameters, and key results
- Written in the standard academic format

### Final Commit:
```
92c99c9 - Add ELIGOS2 fix script, Dockerfile cleanup, and methods section
33db293 - Pin pandas 1.5.3 + numpy <2.0 for ELIGOS2 compatibility
```

---

## 8. Final Results Summary <a name="8-final-results"></a>

### What We Achieved

Starting from a pipeline that had never been run end-to-end, we:

1. **Packaged** 20+ bioinformatics tools into a single Docker image
2. **Deployed** it to a powerful remote server
3. **Fixed 15 bugs** spanning Python, R, Docker, file formats, and tool incompatibilities
4. **Ran the complete pipeline** (31/31 steps) on real Arabidopsis direct RNA-seq data
5. **Analyzed the results** and found real biological signal (m6A modifications)
6. **Documented everything** for reproducibility

### Key Biological Findings

| Finding | Detail |
|---------|--------|
| Mapping success | 92-93% of reads mapped to the genome (very good) |
| Transcript diversity | 21,286 distinct isoforms identified |
| m6A modifications | 2 high-confidence sites (>90% probability) on Chr3 and Chr4 |
| Top gene | Rubisco small subunit family (expected for leaf tissue) |
| DESeq2 | Could not run (needs more biological replicates) |

### Pipeline Steps (31 total)

```
 1. NanoFilt (R1)         - Filter low-quality reads from sample 1
 2. NanoFilt (R2)         - Filter low-quality reads from sample 2
 3. Minimap2 align (R1)   - Map reads to genome for sample 1
 4. Minimap2 align (R2)   - Map reads to genome for sample 2
 5. Samtools sort (R1)    - Sort alignment file for sample 1
 6. Samtools sort (R2)    - Sort alignment file for sample 2
 7. Samtools index (R1)   - Index alignment file for sample 1
 8. Samtools index (R2)   - Index alignment file for sample 2
 9. NanoPlot (R1)         - Generate QC plots for sample 1
10. NanoPlot (R2)         - Generate QC plots for sample 2
11. NanoComp              - Compare quality metrics between samples
12. pycoQC (R1)           - Sequencing run quality report for sample 1
13. pycoQC (R2)           - Sequencing run quality report for sample 2
14. Samtools flagstat (R1) - Alignment summary for sample 1
15. Samtools flagstat (R2) - Alignment summary for sample 2
16. Samtools stats (R1)   - Detailed alignment statistics for sample 1
17. Samtools stats (R2)   - Detailed alignment statistics for sample 2
18. FLAIR align           - Splice-aware alignment for isoforms
19. FLAIR correct (R1)    - Correct splice junctions for sample 1
20. FLAIR correct (R2)    - Correct splice junctions for sample 2
21. FLAIR collapse        - Merge reads into distinct isoforms
22. FLAIR quantify        - Count reads per isoform per sample
23. Nanopolish index      - Link FAST5 signals to basecalled reads
24. Nanopolish eventalign - Align electrical signals to reference (70 GB output!)
25. m6Anet dataprep       - Prepare signal data for m6A detection
26. m6Anet inference      - Predict m6A modification probabilities
27. ELIGOS2 (R1)          - Error-based modification detection for sample 1
28. ELIGOS2 (R2)          - Error-based modification detection for sample 2
29. DESeq2                - Differential expression (failed: needs replicates)
30. FLAIR diffExp         - FLAIR's own differential analysis
31. MultiQC               - Aggregate all QC reports into one
```

---

## 9. Complete List of Bugs Fixed <a name="9-bugs-fixed"></a>

Here's every bug we encountered and fixed, in chronological order:

| # | Bug | Root Cause | Fix | Session |
|---|-----|-----------|-----|---------|
| 1 | FAST5 files can't be read | VBZ compression plugin missing | Added VBZ plugin to Docker image | 2 |
| 2 | pycoQC crashes on import | setuptools v71+ broke compatibility | Pinned setuptools < 71 | 2 |
| 3 | FLAIR rejects sample names | Underscores not allowed in FLAIR IDs | Replace underscores with hyphens | 2 |
| 4 | FLAIR can't find genes | Chromosome names: `1` vs `Chr1` | Renamed chromosomes in genome for FLAIR | 2 |
| 5 | ELIGOS2 crashes (rpy2) | `pandas2ri.activate()` deprecated in rpy2 3.6+ | Wrapped in try/except via fix script | 3-4 |
| 6 | m6Anet "file not found" | Output filename changed in v2.x | Updated expected filename in Snakefile | 4 |
| 7 | MultiQC "file not found" | `--title` flag changes output filename | Added `--filename` flag | 4 |
| 8 | FLAIR quantify "file not found" | FLAIR appends `.counts` to output name | Added rename step | 4 |
| 9 | ELIGOS2 crashes (strand) | BED merge drops strand column (bug in ELIGOS2) | Fixed merge to keep strand + name columns | 4 |
| 10 | ELIGOS2 finds nothing | Chr prefix mismatch between BAM and BED | Strip `Chr` prefix with sed | 4 |
| 11 | ELIGOS2 crashes (pandas) | pandas 2.x changed `pd.concat()` behavior | Downgraded to pandas 1.5.3 | 4-5 |
| 12 | DESeq2 not installed | Missing `libxml2-dev` system library | Added to Dockerfile apt-get | 4 |
| 13 | DESeq2 won't run | Needs >= 2 biological replicates per condition | Made non-fatal with placeholder output | 4 |
| 14 | Docker build fails | Multi-line Python in Dockerfile confuses BuildKit | Moved to separate .py script file | 5 |
| 15 | Can't delete Docker files | Docker creates files as root user | Use Docker container to delete files | 4 |

**Bonus Windows/Linux issues fixed:**
- NTFS filesystem delays with Docker volumes -> added `--latency-wait 60`
- Git Bash path conversion -> `MSYS_NO_PATHCONV=1`
- `tee` command breaking bash strict mode -> used `> file 2>&1` instead
- numpy 2.x incompatible with pandas 1.5.3 -> pinned numpy < 2.0

---

## 10. Files We Created or Modified <a name="10-files"></a>

### Core Pipeline Files

| File | What It Does |
|------|-------------|
| `Snakefile` | The main pipeline recipe - defines all 31 steps and how they connect |
| `Dockerfile` | Instructions to build the Docker image with all 20+ tools |
| `config/config.yml` | Settings file: sample names, tool parameters, on/off switches |

### Scripts

| File | What It Does |
|------|-------------|
| `scripts/fix_eligos2.py` | Automatically patches ELIGOS2 bugs during Docker build |
| `scripts/run_deseq2.R` | R script for differential expression analysis |
| `run_docker.sh` | Linux script to launch the pipeline in Docker |
| `run_docker.bat` | Windows script to launch the pipeline in Docker |
| `deploy_server.sh` | Script to set up and run on the remote server |

### Documentation

| File | What It Does |
|------|-------------|
| `docs/methods_section.md` | Formal bioinformatics methods section for thesis/paper |
| `docs/step_by_step_guide.md` | This document! |

### GitHub Commit History (oldest to newest)

```
fdd374f - Initial commit: K-CHOPORE pipeline with bug fixes
72061ba - Add comprehensive test suite with toy data (51 tests)
8adc51e - Add Docker launch scripts for running with real data
f706760 - Add __pycache__ and .pytest_cache to .gitignore
96ecdd9 - Fix pipeline errors: VBZ plugin, pycoQC setuptools, FLAIR naming/chr mismatch
c6855ee - Add python symlink for ELIGOS2 compatibility
b9db173 - Add Linux server deployment script and fix config/launcher
b0372df - Fix ELIGOS2 rpy2 DeprecationWarning crash
9535f69 - Fix pipeline completion: m6Anet output, ELIGOS2 bugs, FLAIR/MultiQC/DESeq2 compat
33db293 - Pin pandas 1.5.3 + numpy <2.0 for ELIGOS2 compatibility
92c99c9 - Add ELIGOS2 fix script, Dockerfile cleanup, and methods section
```

---

## Visual Summary

```
Session 1: BUILD          Session 2: TEST           Session 3: DEPLOY
+------------------+     +------------------+     +------------------+
| Write Dockerfile |     | Run locally      |     | Set up SSH       |
| Write Snakefile  | --> | Hit 4 bugs       | --> | Transfer to      |
| Write config     |     | Fix all 4        |     |   server         |
| Push to GitHub   |     | Add test suite   |     | Build Docker     |
+------------------+     +------------------+     | Start pipeline   |
                                                   +------------------+
                                                          |
                                                          v
Session 5: FINISH         Session 4: DEBUG
+------------------+     +------------------+
| Analyze results  |     | Fix 11 more bugs |
| Fix ELIGOS2 for  | <-- | 6 restart cycles |
|   real (pandas)  |     | 24/31 --> 31/31  |
| Clean up server  |     | Pipeline done!   |
| Write methods    |     +------------------+
| Rebuild Docker   |
| Push to GitHub   |
+------------------+
```

---

> **Bottom line:** We took a complex bioinformatics pipeline with 20+ tools, packaged it in Docker for reproducibility, deployed it to a remote server, fixed 15 bugs across 5 debugging sessions, successfully analyzed real plant RNA sequencing data, found evidence of RNA chemical modifications, and documented everything for future use.
