# =============================================================
# K-CHOPORE Pipeline - Snakemake Workflow
# =============================================================
# Comprehensive pipeline for ONT Direct RNA Sequencing analysis
# including basecalling, QC, alignment, isoform analysis,
# epitranscriptomic modification detection, and differential expression.
#
# Created by Pelayo Gonzalez de Lena Rodriguez, MSc
# FPI Severo Ochoa Fellow
# Cancer Epigenetics and Nanomedicine Lab | FINBA
# Systems Biology Lab | University of Oviedo
# https://www.linkedin.com/in/biopelayo/
# https://gitlab.com/bio.pelayo/
# =============================================================

import os

# Load configuration
configfile: "config/config.yml"

# -------------------------------------------------------------
# Global variables from config
# -------------------------------------------------------------
SAMPLES = config["samples"]
THREADS = config["params"]["threads"]

# Reference files
REFERENCE_GENOME = config["input_files"]["reference_genome"]
REFERENCE_INDEX = config["input_files"]["reference_genome_mmi"]
GTF_FILE = config["input_files"]["gtf_file"]

# Output directories
OUT = config["output"]["path"]
LOGS = config["output"]["logs"]

# Tool settings
MINIMAP2_PRESET = config["tools"]["minimap2_preset"]
MINIMAP2_KMER = config["tools"]["minimap2_kmer_size"]
MINIMAP2_EXTRA = config["tools"]["minimap2_extra_flags"]

# Build sample-to-sequencing-summary mapping
# This resolves the mismatch between simple {sample}_sequencing_summary.txt patterns
# and actual filenames that may contain flowcell IDs (e.g., WT_C_R1_sequencing_summary_FAR90122_d34138fc.txt)
_seq_summaries = config["input_files"].get("sequencing_summaries", [])
SEQUENCING_SUMMARY = {}
for i, sample in enumerate(SAMPLES):
    if i < len(_seq_summaries):
        SEQUENCING_SUMMARY[sample] = _seq_summaries[i]
    else:
        # Fallback: use default naming pattern
        SEQUENCING_SUMMARY[sample] = f"data/raw/summaries/{sample}_sequencing_summary.txt"

# Print configuration for debugging
print(f"[K-CHOPORE] Samples: {SAMPLES}")
print(f"[K-CHOPORE] Reference genome: {REFERENCE_GENOME}")
print(f"[K-CHOPORE] Threads: {THREADS}")
print(f"[K-CHOPORE] Minimap2 preset: {MINIMAP2_PRESET} (k={MINIMAP2_KMER})")
print(f"[K-CHOPORE] Sequencing summaries: {SEQUENCING_SUMMARY}")

# -------------------------------------------------------------
# Helper: collect conditional targets for the 'all' rule
# -------------------------------------------------------------
def get_all_targets():
    targets = []

    # Always: directory structure
    targets.append("complete_structure_created.txt")

    # Always: genome index
    targets.append(REFERENCE_INDEX)

    # Basecalling (optional, if starting from FAST5/POD5)
    if config["params"].get("run_basecalling", False):
        targets.extend(
            expand("results/basecalls/{sample}.fastq", sample=SAMPLES)
        )

    # NanoFilt read filtering
    if config["params"].get("run_nanofilt", True):
        targets.extend(
            expand("results/fastq_filtered/{sample}_filtered.fastq", sample=SAMPLES)
        )

    # NanoPlot QC per sample
    if config["params"].get("run_nanoplot", True):
        targets.extend(
            expand("results/nanoplot/{sample}/NanoStats.txt", sample=SAMPLES)
        )

    # NanoComp comparison
    if config["params"].get("run_nanocomp", True):
        targets.append("results/nanocomp/NanoComp-report.html")

    # Alignment (always)
    targets.extend(
        expand("results/sorted_bam/{sample}_sorted.bam", sample=SAMPLES)
    )

    # Samtools stats
    targets.extend(
        expand("results/samtools_stats/{sample}_flagstat.txt", sample=SAMPLES)
    )
    targets.extend(
        expand("results/samtools_stats/{sample}_stats.txt", sample=SAMPLES)
    )

    # pycoQC
    if config["params"].get("run_pycoqc", True):
        targets.extend(
            expand("results/quality_analysis/pycoQC_output_{sample}.html", sample=SAMPLES)
        )

    # FLAIR isoform analysis
    if config["params"].get("run_flair", True):
        targets.extend(
            expand("results/flair/{sample}_flair.collapse.isoforms.bed", sample=SAMPLES)
        )
        targets.append("results/flair/counts_matrix.tsv")

    # StringTie2 isoform assembly
    if config["params"].get("run_stringtie", False):
        targets.extend(
            expand("results/stringtie/{sample}_stringtie.gtf", sample=SAMPLES)
        )

    # ELIGOS2 epitranscriptomic modification
    if config["params"].get("run_eligos2", True):
        targets.extend(
            expand("results/eligos/{sample}_eligos_output.txt", sample=SAMPLES)
        )

    # m6Anet modification detection
    if config["params"].get("run_m6anet", True):
        targets.extend(
            expand("results/m6anet/{sample}/data.result.csv.gz", sample=SAMPLES)
        )

    # xPore differential modification
    if config["params"].get("run_xpore", False):
        targets.append("results/xpore/diffmod.table")

    # DESeq2 differential expression
    if config["params"].get("run_deseq2", True):
        targets.append("results/deseq2/deseq2_results.csv")

    # MultiQC aggregate report
    if config["params"].get("run_multiqc", True):
        targets.append("results/multiqc/multiqc_report.html")

    return targets

# =============================================================
# RULE ALL - Master target
# =============================================================
rule all:
    input:
        get_all_targets()

# =============================================================
# RULE: Setup project directory structure
# =============================================================
rule setup_complete_structure:
    output:
        "complete_structure_created.txt"
    run:
        dirs = [
            "config", "data", "data/raw/fastq", "data/raw/fast5",
            "data/raw/pod5", "data/raw/summaries",
            "data/reference/genome", "data/reference/annotations",
            "data/reference/transcriptome",
            "docs", "envs", "logs", "notebooks", "publication",
            "results", "reviews", "scripts",
            "results/basecalls", "results/fastq_filtered",
            "results/mapped", "results/sorted_bam",
            "results/quality_analysis", "results/nanoplot",
            "results/nanocomp", "results/samtools_stats",
            "results/flair", "results/stringtie", "results/bambu",
            "results/eligos", "results/m6anet", "results/xpore",
            "results/nanopolish", "results/deseq2", "results/multiqc"
        ]
        for d in dirs:
            os.makedirs(d, exist_ok=True)
            print(f"[K-CHOPORE] Directory ready: {d}")
        with open(output[0], 'w') as f:
            f.write("K-CHOPORE project structure created successfully.")

# =============================================================
# BASECALLING (Optional - from FAST5/POD5 to FASTQ)
# =============================================================

# Rule: Basecall with Dorado (recommended for modern ONT data)
rule basecall_dorado:
    input:
        structure_created="complete_structure_created.txt",
        pod5_dir=config["input_files"]["pod5_dir"]
    output:
        fastq="results/basecalls/{sample}.fastq"
    params:
        model=config["tools"]["dorado_model"],
        dorado=config["tools"]["dorado_path"]
    threads: THREADS
    log:
        "logs/basecall_dorado_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Basecalling with Dorado for {wildcards.sample}..."
        {params.dorado} basecaller {params.model} \
            {input.pod5_dir}/{wildcards.sample}/ \
            --emit-fastq > {output.fastq} 2> {log}
        echo "[K-CHOPORE] Basecalling completed for {wildcards.sample}."
        """

# Rule: Basecall with Guppy (legacy support)
rule basecall_guppy:
    input:
        structure_created="complete_structure_created.txt",
        fast5_dir=config["input_files"]["fast5_dir"]
    output:
        fastq="results/basecalls/{sample}_guppy.fastq"
    params:
        guppy_cfg=config["tools"]["guppy_config_file"]
    threads: THREADS
    log:
        "logs/basecall_guppy_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Basecalling with Guppy for {wildcards.sample}..."
        guppy_basecaller \
            -i {input.fast5_dir}/{wildcards.sample}/ \
            -s results/basecalls/{wildcards.sample}_guppy/ \
            -c {params.guppy_cfg} \
            --num_callers {threads} \
            --compress_fastq 2>&1 | tee {log}
        cat results/basecalls/{wildcards.sample}_guppy/pass/*.fastq > {output.fastq}
        echo "[K-CHOPORE] Guppy basecalling completed for {wildcards.sample}."
        """

# =============================================================
# READ QC AND FILTERING
# =============================================================

# Rule: Filter reads with NanoFilt
rule nanofilt:
    input:
        fastq="data/raw/fastq/{sample}.fastq"
    output:
        filtered="results/fastq_filtered/{sample}_filtered.fastq"
    params:
        min_qual=config["tools"]["nanofilt_min_quality"],
        min_len=config["tools"]["nanofilt_min_length"],
        max_len=config["tools"]["nanofilt_max_length"]
    log:
        "logs/nanofilt_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Filtering reads with NanoFilt for {wildcards.sample}..."
        max_len_flag=""
        if [ {params.max_len} -gt 0 ]; then
            max_len_flag="--maxlength {params.max_len}"
        fi
        NanoFilt -q {params.min_qual} -l {params.min_len} $max_len_flag \
            < {input.fastq} > {output.filtered} 2> {log}
        echo "[K-CHOPORE] NanoFilt completed for {wildcards.sample}."
        """

# Rule: NanoPlot per-sample QC
rule nanoplot:
    input:
        fastq="results/fastq_filtered/{sample}_filtered.fastq"
    output:
        stats="results/nanoplot/{sample}/NanoStats.txt"
    params:
        outdir="results/nanoplot/{sample}",
        fmt=config["tools"]["nanoplot_format"]
    threads: 4
    log:
        "logs/nanoplot_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running NanoPlot QC for {wildcards.sample}..."
        NanoPlot --fastq {input.fastq} \
            --outdir {params.outdir} \
            --format {params.fmt} \
            --threads {threads} \
            --loglength \
            --title "{wildcards.sample} Read QC" \
            --plots dot kde 2>&1 | tee {log}
        echo "[K-CHOPORE] NanoPlot completed for {wildcards.sample}."
        """

# Rule: NanoComp comparative QC across all samples
rule nanocomp:
    input:
        fastqs=expand("results/fastq_filtered/{sample}_filtered.fastq", sample=SAMPLES)
    output:
        report="results/nanocomp/NanoComp-report.html"
    params:
        outdir="results/nanocomp",
        names=" ".join(SAMPLES)
    threads: 4
    log:
        "logs/nanocomp.log"
    shell:
        """
        echo "[K-CHOPORE] Running NanoComp across all samples..."
        NanoComp --fastq {input.fastqs} \
            --names {params.names} \
            --outdir {params.outdir} \
            --threads {threads} \
            --plot violin 2>&1 | tee {log}
        echo "[K-CHOPORE] NanoComp completed."
        """

# =============================================================
# ALIGNMENT
# =============================================================

# Rule: Index reference genome with Minimap2 (splice-aware for RNA)
rule index_genome:
    input:
        structure_created="complete_structure_created.txt",
        reference_genome=REFERENCE_GENOME
    output:
        reference_index=REFERENCE_INDEX
    log:
        "logs/index_genome.log"
    shell:
        """
        echo "[K-CHOPORE] Indexing reference genome for splice-aware RNA alignment..."
        minimap2 -d {output.reference_index} -k {k} {input.reference_genome} 2>&1 | tee {log}
        echo "[K-CHOPORE] Genome indexing completed."
        """.replace("{k}", str(MINIMAP2_KMER))

# Rule: Align reads with Minimap2 (splice-aware for direct RNA-seq)
rule map_with_minimap2:
    input:
        structure_created="complete_structure_created.txt",
        reference_index=REFERENCE_INDEX,
        fastq="results/fastq_filtered/{sample}_filtered.fastq",
        bed=GTF_FILE
    output:
        sam="results/mapped/{sample}.sam"
    params:
        preset=MINIMAP2_PRESET,
        kmer=MINIMAP2_KMER,
        extra=MINIMAP2_EXTRA
    threads: THREADS
    log:
        "logs/minimap2_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Aligning {wildcards.sample} with Minimap2 (splice-aware, direct RNA)..."
        minimap2 -ax {params.preset} \
            -k {params.kmer} \
            -uf \
            --junc-bed {input.bed} \
            {params.extra} \
            -t {threads} \
            {input.reference_index} \
            {input.fastq} > {output.sam} 2> {log}
        echo "[K-CHOPORE] Alignment completed for {wildcards.sample}."
        """

# Rule: Sort and index BAM files
rule sort_and_index_bam:
    input:
        sam="results/mapped/{sample}.sam"
    output:
        bam="results/sorted_bam/{sample}_sorted.bam",
        bai="results/sorted_bam/{sample}_sorted.bam.bai"
    threads: 4
    log:
        "logs/sort_index_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Sorting and indexing BAM for {wildcards.sample}..."
        samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        echo "[K-CHOPORE] BAM sorted and indexed for {wildcards.sample}."
        """

# Rule: Samtools flagstat and stats
rule samtools_stats:
    input:
        bam="results/sorted_bam/{sample}_sorted.bam",
        bai="results/sorted_bam/{sample}_sorted.bam.bai"
    output:
        flagstat="results/samtools_stats/{sample}_flagstat.txt",
        stats="results/samtools_stats/{sample}_stats.txt"
    log:
        "logs/samtools_stats_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Computing alignment statistics for {wildcards.sample}..."
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        samtools stats {input.bam} > {output.stats} 2>> {log}
        echo "[K-CHOPORE] Alignment stats completed for {wildcards.sample}."
        """

# =============================================================
# QUALITY CONTROL
# =============================================================

# Rule: pycoQC quality analysis
rule quality_analysis_with_pycoQC:
    input:
        summary=lambda wildcards: SEQUENCING_SUMMARY[wildcards.sample],
        bam="results/sorted_bam/{sample}_sorted.bam",
        bai="results/sorted_bam/{sample}_sorted.bam.bai"
    output:
        html="results/quality_analysis/pycoQC_output_{sample}.html"
    params:
        min_qual=config["tools"]["pycoqc_min_pass_qual"]
    log:
        "logs/pycoqc_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running pycoQC for {wildcards.sample}..."
        pycoQC -f {input.summary} \
            -a {input.bam} \
            -o {output.html} \
            --min_pass_qual {params.min_qual} 2>&1 | tee {log}
        echo "[K-CHOPORE] pycoQC completed for {wildcards.sample}."
        """

# =============================================================
# ISOFORM ANALYSIS - FLAIR
# =============================================================

# Rule: FLAIR align
rule flair_align:
    input:
        genome=REFERENCE_GENOME,
        fastq="results/fastq_filtered/{sample}_filtered.fastq"
    output:
        bed="results/flair/{sample}_flair.aligned.bed"
    params:
        outprefix="results/flair/{sample}_flair"
    threads: THREADS
    log:
        "logs/flair_align_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running FLAIR align for {wildcards.sample}..."
        flair align \
            -g {input.genome} \
            -r {input.fastq} \
            -o {params.outprefix} \
            --threads {threads} \
            --nvrna 2>&1 | tee {log}
        echo "[K-CHOPORE] FLAIR align completed for {wildcards.sample}."
        """

# Rule: FLAIR correct
rule flair_correct:
    input:
        bed="results/flair/{sample}_flair.aligned.bed",
        genome=REFERENCE_GENOME,
        gtf=GTF_FILE
    output:
        corrected_bed="results/flair/{sample}_flair_all_corrected.bed"
    params:
        outprefix="results/flair/{sample}_flair"
    threads: THREADS
    log:
        "logs/flair_correct_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running FLAIR correct for {wildcards.sample}..."
        flair correct \
            -q {input.bed} \
            -g {input.genome} \
            -f {input.gtf} \
            -o {params.outprefix} \
            --threads {threads} 2>&1 | tee {log}
        echo "[K-CHOPORE] FLAIR correct completed for {wildcards.sample}."
        """

# Rule: FLAIR collapse
rule flair_collapse:
    input:
        corrected_bed="results/flair/{sample}_flair_all_corrected.bed",
        genome=REFERENCE_GENOME,
        gtf=GTF_FILE,
        fastq="results/fastq_filtered/{sample}_filtered.fastq"
    output:
        isoforms_bed="results/flair/{sample}_flair.collapse.isoforms.bed",
        isoforms_fa="results/flair/{sample}_flair.collapse.isoforms.fa",
        isoforms_gtf="results/flair/{sample}_flair.collapse.isoforms.gtf"
    params:
        outprefix="results/flair/{sample}_flair.collapse",
        support=config["tools"]["flair_support"]
    threads: THREADS
    log:
        "logs/flair_collapse_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running FLAIR collapse for {wildcards.sample}..."
        flair collapse \
            -g {input.genome} \
            -r {input.fastq} \
            -q {input.corrected_bed} \
            -f {input.gtf} \
            -o {params.outprefix} \
            -s {params.support} \
            --threads {threads} 2>&1 | tee {log}
        echo "[K-CHOPORE] FLAIR collapse completed for {wildcards.sample}."
        """

# Rule: FLAIR quantify
rule flair_quantify:
    input:
        isoforms_fa=expand("results/flair/{sample}_flair.collapse.isoforms.fa", sample=SAMPLES),
        reads_manifest=config["input_files"]["reads_manifest"]
    output:
        counts="results/flair/counts_matrix.tsv"
    params:
        threads=THREADS
    log:
        "logs/flair_quantify.log"
    shell:
        """
        echo "[K-CHOPORE] Quantifying isoforms with FLAIR..."
        flair quantify \
            -r {input.reads_manifest} \
            -i {input.isoforms_fa[0]} \
            --tpm \
            --threads {params.threads} \
            -o results/flair/counts_matrix 2>&1 | tee {log}
        echo "[K-CHOPORE] FLAIR quantification completed."
        """

# Rule: FLAIR differential expression
rule flair_diff_exp:
    input:
        counts_matrix="results/flair/counts_matrix.tsv"
    output:
        diff_exp="results/flair/diffExp/genes_deseq2_sig.tsv"
    params:
        out_dir="results/flair/diffExp"
    log:
        "logs/flair_diffexp.log"
    shell:
        """
        echo "[K-CHOPORE] Running FLAIR differential expression..."
        flair diffExp \
            -q {input.counts_matrix} \
            -o {params.out_dir} 2>&1 | tee {log}
        echo "[K-CHOPORE] FLAIR diffExp completed."
        """

# Rule: FLAIR differential splicing
rule flair_diff_splice:
    input:
        isoforms_bed=expand("results/flair/{sample}_flair.collapse.isoforms.bed", sample=SAMPLES[0]),
        counts_matrix="results/flair/counts_matrix.tsv"
    output:
        diff_splice="results/flair/diffSplice/diffsplice.alt3.events.quant.tsv"
    params:
        out_dir="results/flair/diffSplice"
    log:
        "logs/flair_diffsplice.log"
    shell:
        """
        echo "[K-CHOPORE] Running FLAIR differential splicing..."
        flair diffSplice \
            -i {input.isoforms_bed} \
            -q {input.counts_matrix} \
            --test 2>&1 | tee {log}
        echo "[K-CHOPORE] FLAIR diffSplice completed."
        """

# =============================================================
# ISOFORM ANALYSIS - StringTie2 (Alternative)
# =============================================================

rule stringtie_assemble:
    input:
        bam="results/sorted_bam/{sample}_sorted.bam",
        bai="results/sorted_bam/{sample}_sorted.bam.bai",
        gtf=GTF_FILE
    output:
        gtf="results/stringtie/{sample}_stringtie.gtf"
    params:
        min_cov=config["tools"]["stringtie_min_cov"],
        min_tpm=config["tools"]["stringtie_min_tpm"]
    threads: THREADS
    log:
        "logs/stringtie_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running StringTie2 for {wildcards.sample}..."
        stringtie {input.bam} \
            -G {input.gtf} \
            -o {output.gtf} \
            -p {threads} \
            -L \
            -c {params.min_cov} \
            -A results/stringtie/{wildcards.sample}_gene_abund.tab 2>&1 | tee {log}
        echo "[K-CHOPORE] StringTie2 completed for {wildcards.sample}."
        """

# =============================================================
# EPITRANSCRIPTOMIC MODIFICATION DETECTION - ELIGOS2
# =============================================================

rule run_eligos2:
    input:
        bam="results/sorted_bam/{sample}_sorted.bam",
        bai="results/sorted_bam/{sample}_sorted.bam.bai",
        reference_genome=REFERENCE_GENOME,
        region_bed="results/flair/{sample}_flair.collapse.isoforms.bed"
    output:
        eligos_output="results/eligos/{sample}_eligos_output.txt"
    params:
        pval=config["tools"]["eligos2_pval"],
        oddR=config["tools"]["eligos2_oddR"],
        esb=config["tools"]["eligos2_esb"],
        outdir="results/eligos/{sample}"
    threads: THREADS
    log:
        "logs/eligos2_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running ELIGOS2 for {wildcards.sample}..."
        eligos2 rna_mod \
            -i {input.bam} \
            -reg {input.region_bed} \
            -ref {input.reference_genome} \
            -o {params.outdir} \
            --pval {params.pval} \
            --oddR {params.oddR} \
            --esb {params.esb} \
            --threads {threads} 2>&1 | tee {log}
        # Collect main output
        cp {params.outdir}/*_baseExt0.txt {output.eligos_output} 2>/dev/null || \
            touch {output.eligos_output}
        echo "[K-CHOPORE] ELIGOS2 completed for {wildcards.sample}."
        """

# =============================================================
# EPITRANSCRIPTOMIC MODIFICATION DETECTION - m6Anet
# =============================================================

# Step 1: Nanopolish index FAST5 files for signal-level access
rule nanopolish_index:
    input:
        fastq="results/fastq_filtered/{sample}_filtered.fastq",
        fast5_dir=config["input_files"]["fast5_dir"]
    output:
        index=touch("results/nanopolish/{sample}_index.done")
    log:
        "logs/nanopolish_index_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Indexing FAST5 for Nanopolish ({wildcards.sample})..."
        nanopolish index \
            -d {input.fast5_dir}/{wildcards.sample}/ \
            {input.fastq} 2>&1 | tee {log}
        echo "[K-CHOPORE] Nanopolish indexing completed for {wildcards.sample}."
        """

# Step 2: Nanopolish eventalign for signal-level data
rule nanopolish_eventalign:
    input:
        fastq="results/fastq_filtered/{sample}_filtered.fastq",
        bam="results/sorted_bam/{sample}_sorted.bam",
        bai="results/sorted_bam/{sample}_sorted.bam.bai",
        genome=REFERENCE_GENOME,
        index_done="results/nanopolish/{sample}_index.done"
    output:
        eventalign="results/nanopolish/{sample}_eventalign.txt"
    threads: THREADS
    log:
        "logs/nanopolish_eventalign_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running Nanopolish eventalign for {wildcards.sample}..."
        nanopolish eventalign \
            --reads {input.fastq} \
            --bam {input.bam} \
            --genome {input.genome} \
            --signal-index \
            --scale-events \
            --summary results/nanopolish/{wildcards.sample}_summary.txt \
            --threads {threads} > {output.eventalign} 2> {log}
        echo "[K-CHOPORE] Nanopolish eventalign completed for {wildcards.sample}."
        """

# Step 3: m6Anet dataprep - prepare data for m6A inference
rule m6anet_dataprep:
    input:
        eventalign="results/nanopolish/{sample}_eventalign.txt"
    output:
        dataprep_done=touch("results/m6anet/{sample}/dataprep.done")
    params:
        outdir="results/m6anet/{sample}"
    threads: config["tools"]["m6anet_num_processors"]
    log:
        "logs/m6anet_dataprep_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running m6Anet dataprep for {wildcards.sample}..."
        m6anet dataprep \
            --eventalign {input.eventalign} \
            --out_dir {params.outdir} \
            --n_processes {threads} 2>&1 | tee {log}
        echo "[K-CHOPORE] m6Anet dataprep completed for {wildcards.sample}."
        """

# Step 4: m6Anet inference - detect m6A modifications
rule m6anet_inference:
    input:
        dataprep_done="results/m6anet/{sample}/dataprep.done"
    output:
        result="results/m6anet/{sample}/data.result.csv.gz"
    params:
        indir="results/m6anet/{sample}",
        outdir="results/m6anet/{sample}",
        n_iters=config["tools"]["m6anet_num_iterations"],
        n_proc=config["tools"]["m6anet_num_processors"]
    log:
        "logs/m6anet_inference_{sample}.log"
    shell:
        """
        echo "[K-CHOPORE] Running m6Anet inference for {wildcards.sample}..."
        m6anet inference \
            --input_dir {params.indir} \
            --out_dir {params.outdir} \
            --n_processes {params.n_proc} \
            --num_iterations {params.n_iters} 2>&1 | tee {log}
        echo "[K-CHOPORE] m6Anet inference completed for {wildcards.sample}."
        """

# =============================================================
# EPITRANSCRIPTOMIC MODIFICATION - xPore (Differential)
# =============================================================

# xPore requires eventalign data from nanopolish for multiple conditions
rule xpore_diffmod:
    input:
        eventalign=expand("results/nanopolish/{sample}_eventalign.txt", sample=SAMPLES)
    output:
        table="results/xpore/diffmod.table"
    params:
        outdir="results/xpore",
        min_reads=config["tools"]["xpore_min_reads"]
    threads: THREADS
    log:
        "logs/xpore_diffmod.log"
    shell:
        """
        echo "[K-CHOPORE] Running xPore differential modification analysis..."
        # xPore requires a YAML config - generate it
        cat > results/xpore/xpore_config.yml << 'XPORE_EOF'
notes: K-CHOPORE xPore analysis
out: {params.outdir}
XPORE_EOF
        xpore diffmod \
            --config results/xpore/xpore_config.yml \
            --n_processes {threads} 2>&1 | tee {log}
        echo "[K-CHOPORE] xPore analysis completed."
        """

# =============================================================
# DIFFERENTIAL EXPRESSION - DESeq2
# =============================================================

rule run_deseq2:
    input:
        counts="results/flair/counts_matrix.tsv"
    output:
        results="results/deseq2/deseq2_results.csv",
        ma_plot="results/deseq2/MA_plot.pdf",
        volcano="results/deseq2/volcano_plot.pdf"
    params:
        padj=config["tools"]["deseq2_padj_threshold"],
        lfc=config["tools"]["deseq2_lfc_threshold"],
        outdir="results/deseq2"
    log:
        "logs/deseq2.log"
    shell:
        """
        echo "[K-CHOPORE] Running DESeq2 differential expression analysis..."
        Rscript scripts/run_deseq2.R \
            {input.counts} \
            {params.outdir} \
            {params.padj} \
            {params.lfc} 2>&1 | tee {log}
        echo "[K-CHOPORE] DESeq2 analysis completed."
        """

# =============================================================
# MULTIQC - Aggregate all QC reports
# =============================================================

rule multiqc:
    input:
        nanoplot=expand("results/nanoplot/{sample}/NanoStats.txt", sample=SAMPLES) if config["params"].get("run_nanoplot", True) else [],
        flagstat=expand("results/samtools_stats/{sample}_flagstat.txt", sample=SAMPLES),
        stats=expand("results/samtools_stats/{sample}_stats.txt", sample=SAMPLES),
        pycoqc=expand("results/quality_analysis/pycoQC_output_{sample}.html", sample=SAMPLES) if config["params"].get("run_pycoqc", True) else []
    output:
        report="results/multiqc/multiqc_report.html"
    params:
        outdir="results/multiqc"
    log:
        "logs/multiqc.log"
    shell:
        """
        echo "[K-CHOPORE] Aggregating QC reports with MultiQC..."
        multiqc results/ \
            -o {params.outdir} \
            --force \
            --title "K-CHOPORE QC Report" 2>&1 | tee {log}
        echo "[K-CHOPORE] MultiQC report generated."
        """
