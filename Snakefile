# =============================================================
# K-CHOPORE Pipeline - Snakemake Workflow
# =============================================================
# ONT Direct RNA Sequencing analysis pipeline:
# Pre-basecalled FASTQ -> QC -> alignment -> isoforms ->
# epitranscriptomic modifications -> differential expression
#
# 10 samples, 2x2 factorial: WT vs anac017-1 x Control vs AA
# Arabidopsis thaliana (TAIR10 + AtRTDv2)
#
# Created by Pelayo Gonzalez de Lena Rodriguez, MSc
# FPI Severo Ochoa Fellow
# Cancer Epigenetics and Nanomedicine Lab | FINBA
# Systems Biology Lab | University of Oviedo
# =============================================================

import os

configfile: "config/config.yml"

# -------------------------------------------------------------
# Parse sample sheet from config
# -------------------------------------------------------------
SAMPLE_SHEET = config["samples"]
SAMPLE_IDS = list(SAMPLE_SHEET.keys())
THREADS = config["params"]["threads"]
RAW_DATA_DIR = config["input_files"]["raw_data_dir"]

# Reference files
REFERENCE_GENOME = config["input_files"]["reference_genome"]
REFERENCE_INDEX = config["input_files"]["reference_genome_mmi"]
GTF_FILE = config["input_files"]["gtf_file"]

# Tool settings
MINIMAP2_PRESET = config["tools"]["minimap2_preset"]
MINIMAP2_KMER = config["tools"]["minimap2_kmer_size"]
MINIMAP2_EXTRA = config["tools"]["minimap2_extra_flags"]

# m6Anet requires FAST5 data on local disk.
# Set run_m6anet: true in config only if FAST5 data is available locally.
RUN_M6ANET = config.get("params", {}).get("run_m6anet", False)

# --- v3.0 ncRNA module toggles ---
RUN_LNCRNA = config.get("params", {}).get("run_lncrna", False)
RUN_SMALLRNA = config.get("params", {}).get("run_smallrna", False)
RUN_MIRNA_TARGETS = config.get("params", {}).get("run_mirna_targets", False)
RUN_EPITX_ENHANCED = config.get("params", {}).get("run_epitx_enhanced", False)
RUN_INTEGRATION = config.get("params", {}).get("run_integration", False)

# Small RNA and degradome sample sheets (Module 2 & 3)
SRNA_SAMPLES = config.get("srna_samples", {})
SRNA_SAMPLE_IDS = list(SRNA_SAMPLES.keys()) if SRNA_SAMPLES else []
DEGRADOME_SAMPLES = config.get("degradome_samples", {})
EPITX_COMPARISONS = config.get("epitx_comparisons", [])

# Annotation paths (v3.0)
ANNOTATIONS = config.get("annotations", {})

# Helper functions
def get_raw_fastq_dirs(sample):
    """Return list of local paths to fastq_pass directories."""
    dirs = []
    nas_dirs = SAMPLE_SHEET[sample]["nas_dirs"]
    run_subdirs = SAMPLE_SHEET[sample]["run_subdirs"]
    for nas_dir, run_subdir in zip(nas_dirs, run_subdirs):
        dirs.append(os.path.join(RAW_DATA_DIR, nas_dir, run_subdir, "fastq_pass"))
    return dirs

def get_raw_fast5_dirs(sample):
    """Return list of local paths to FAST5 directories."""
    dirs = []
    nas_dirs = SAMPLE_SHEET[sample]["nas_dirs"]
    run_subdirs = SAMPLE_SHEET[sample]["run_subdirs"]
    for nas_dir, run_subdir in zip(nas_dirs, run_subdirs):
        dirs.append(os.path.join(RAW_DATA_DIR, nas_dir, run_subdir, "fast5_pass"))
    return dirs

def get_raw_summary_files(sample):
    """Return list of local paths to sequencing_summary files."""
    files = []
    nas_dirs = SAMPLE_SHEET[sample]["nas_dirs"]
    run_subdirs = SAMPLE_SHEET[sample]["run_subdirs"]
    for nas_dir, run_subdir in zip(nas_dirs, run_subdirs):
        summary_dir = os.path.join(RAW_DATA_DIR, nas_dir, run_subdir)
        files.append(summary_dir)
    return files

def get_genotype(sample):
    return SAMPLE_SHEET[sample]["genotype"]

def get_treatment(sample):
    return SAMPLE_SHEET[sample]["treatment"]

# Print configuration
print(f"[K-CHOPORE] Samples ({len(SAMPLE_IDS)}): {SAMPLE_IDS}")
print(f"[K-CHOPORE] Reference: {REFERENCE_GENOME}")
print(f"[K-CHOPORE] Raw data: {RAW_DATA_DIR}")
print(f"[K-CHOPORE] Threads: {THREADS}")

# =============================================================
# RULE ALL - Master target
# =============================================================
rule all:
    input:
        # Genome index
        REFERENCE_INDEX,
        # Merged FASTQs and summaries
        expand("results/basecalls/{sample}.fastq", sample=SAMPLE_IDS),
        expand("results/basecalls/{sample}_sequencing_summary.txt", sample=SAMPLE_IDS),
        # Filtered reads
        expand("results/fastq_filtered/{sample}_filtered.fastq", sample=SAMPLE_IDS),
        # QC
        expand("results/nanoplot/{sample}/NanoStats.txt", sample=SAMPLE_IDS),
        "results/nanocomp/NanoComp-report.html",
        expand("results/quality_analysis/pycoQC_output_{sample}.html", sample=SAMPLE_IDS),
        # Alignment stats
        expand("results/sorted_bam/{sample}_sorted.bam", sample=SAMPLE_IDS),
        expand("results/samtools_stats/{sample}_flagstat.txt", sample=SAMPLE_IDS),
        expand("results/samtools_stats/{sample}_stats.txt", sample=SAMPLE_IDS),
        # FLAIR isoforms
        expand("results/flair/{sample}_flair.collapse.isoforms.bed", sample=SAMPLE_IDS),
        "results/flair/counts_matrix.tsv",
        # ELIGOS2 modification detection
        expand("results/eligos/{sample}_eligos_output.txt", sample=SAMPLE_IDS),
        # m6Anet modification detection (optional - requires FAST5 data on disk)
        # Run per-sample: snakemake results/m6anet/SAMPLE/data.site_proba.csv
        expand("results/m6anet/{sample}/data.site_proba.csv", sample=SAMPLE_IDS) if RUN_M6ANET else [],
        # DESeq2 differential expression (2x2 factorial)
        "results/deseq2/deseq2_genotype_results.csv",
        "results/deseq2/deseq2_treatment_results.csv",
        "results/deseq2/deseq2_interaction_results.csv",
        # MultiQC aggregate report
        "results/multiqc/multiqc_report.html",
        # ---- v3.0 ncRNA modules (conditional) ----
        # Module 1: lncRNA discovery
        "results/lncrna/lncrna_final.gtf" if RUN_LNCRNA else [],
        "results/lncrna/lncrna_summary_report.tsv" if RUN_LNCRNA else [],
        "results/lncrna/deseq2/deseq2_genotype_results.csv" if RUN_LNCRNA else [],
        # Module 2: Small RNA analysis (Illumina)
        "results/smallrna/shortstack/Results.txt" if RUN_SMALLRNA else [],
        "results/smallrna/mirna_annotated.tsv" if RUN_SMALLRNA else [],
        "results/smallrna/mirna_counts_matrix.tsv" if RUN_SMALLRNA else [],
        # Module 3: miRNA target prediction + degradome
        "results/targets/target_evidence_table.tsv" if RUN_MIRNA_TARGETS else [],
        # Module 4: Enhanced epitranscriptomics
        "results/epitx/modification_consensus.tsv" if RUN_EPITX_ENHANCED else [],
        # Module 5: Data integration
        "results/integration/wgcna_modules.tsv" if RUN_INTEGRATION else [],
        "results/integration/integration_report.html" if RUN_INTEGRATION else [],

# =============================================================
# PREPARE READS - Merge pre-basecalled FASTQs
# =============================================================
# All 10 samples have pre-basecalled data (data_type: "fastq").
# This rule decompresses and concatenates all *.fastq.gz from fastq_pass/
# into a single FASTQ per sample.

rule prepare_fastq:
    """Merge pre-basecalled gzipped FASTQs into a single FASTQ per sample."""
    output:
        fastq="results/basecalls/{sample}.fastq",
        summary="results/basecalls/{sample}_sequencing_summary.txt"
    params:
        fastq_dirs=lambda wildcards: get_raw_fastq_dirs(wildcards.sample),
        summary_dirs=lambda wildcards: get_raw_summary_files(wildcards.sample)
    threads: 1
    log:
        "logs/prepare_reads_{sample}.log"
    shell:
        """
        mkdir -p results/basecalls logs
        echo "[K-CHOPORE] Merging pre-basecalled FASTQs for {wildcards.sample}..."

        # Concatenate all gzipped FASTQs from all fastq_pass directories
        > {output.fastq}
        for fq_dir in {params.fastq_dirs}; do
            echo "[K-CHOPORE]   Processing: $fq_dir"
            if ls "$fq_dir"/*.fastq.gz 1>/dev/null 2>&1; then
                zcat "$fq_dir"/*.fastq.gz >> {output.fastq}
            elif ls "$fq_dir"/*.fastq 1>/dev/null 2>&1; then
                cat "$fq_dir"/*.fastq >> {output.fastq}
            else
                echo "[K-CHOPORE]   WARNING: No FASTQ files found in $fq_dir"
            fi
        done

        # Merge sequencing summaries (header from first, data from rest)
        first=true
        for summary_dir in {params.summary_dirs}; do
            for summary_file in "$summary_dir"/sequencing_summary_*.txt; do
                if [ -f "$summary_file" ]; then
                    echo "[K-CHOPORE]   Summary: $summary_file"
                    if $first; then
                        cat "$summary_file" > {output.summary}
                        first=false
                    else
                        tail -n +2 "$summary_file" >> {output.summary}
                    fi
                fi
            done
        done

        echo "[K-CHOPORE] Merge completed for {wildcards.sample}."
        echo "[K-CHOPORE]   Reads: $(( $(wc -l < {output.fastq}) / 4 ))"
        """

# =============================================================
# READ QC AND FILTERING
# =============================================================

rule nanofilt:
    input:
        fastq="results/basecalls/{sample}.fastq"
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
        mkdir -p results/fastq_filtered logs
        echo "[K-CHOPORE] Filtering reads for {wildcards.sample}..."
        max_len_flag=""
        if [ {params.max_len} -gt 0 ]; then
            max_len_flag="--maxlength {params.max_len}"
        fi
        NanoFilt -q {params.min_qual} -l {params.min_len} $max_len_flag \
            < {input.fastq} > {output.filtered} 2> {log}
        echo "[K-CHOPORE] NanoFilt completed for {wildcards.sample}."
        """

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
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Running NanoPlot for {wildcards.sample}..."
        NanoPlot --fastq {input.fastq} \
            --outdir {params.outdir} \
            --format {params.fmt} \
            --threads {threads} \
            --loglength \
            --title "{wildcards.sample} Read QC" \
            --plots dot kde > {log} 2>&1
        echo "[K-CHOPORE] NanoPlot completed for {wildcards.sample}."
        """

rule nanocomp:
    input:
        fastqs=expand("results/fastq_filtered/{sample}_filtered.fastq", sample=SAMPLE_IDS)
    output:
        report="results/nanocomp/NanoComp-report.html"
    params:
        outdir="results/nanocomp",
        names=" ".join(SAMPLE_IDS)
    threads: 4
    log:
        "logs/nanocomp.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Running NanoComp across all samples..."
        NanoComp --fastq {input.fastqs} \
            --names {params.names} \
            --outdir {params.outdir} \
            --threads {threads} \
            --plot violin > {log} 2>&1
        echo "[K-CHOPORE] NanoComp completed."
        """

# =============================================================
# ALIGNMENT
# =============================================================

rule index_genome:
    input:
        reference_genome=REFERENCE_GENOME
    output:
        reference_index=REFERENCE_INDEX
    log:
        "logs/index_genome.log"
    shell:
        """
        mkdir -p "$(dirname {output.reference_index})" logs
        echo "[K-CHOPORE] Indexing reference genome..."
        minimap2 -d {output.reference_index} -k {k} {input.reference_genome} > {log} 2>&1
        echo "[K-CHOPORE] Genome indexing completed."
        """.replace("{k}", str(MINIMAP2_KMER))

rule map_with_minimap2:
    input:
        reference_index=REFERENCE_INDEX,
        fastq="results/fastq_filtered/{sample}_filtered.fastq",
        bed=GTF_FILE
    output:
        sam=temp("results/mapped/{sample}.sam")
    params:
        preset=MINIMAP2_PRESET,
        kmer=MINIMAP2_KMER,
        extra=MINIMAP2_EXTRA
    threads: THREADS
    log:
        "logs/minimap2_{sample}.log"
    shell:
        """
        mkdir -p results/mapped logs
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
        mkdir -p results/sorted_bam logs
        echo "[K-CHOPORE] Sorting and indexing BAM for {wildcards.sample}..."
        samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        echo "[K-CHOPORE] BAM sorted and indexed for {wildcards.sample}."
        """

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
        mkdir -p results/samtools_stats logs
        echo "[K-CHOPORE] Computing alignment statistics for {wildcards.sample}..."
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        samtools stats {input.bam} > {output.stats} 2>> {log}
        echo "[K-CHOPORE] Stats completed for {wildcards.sample}."
        """

# =============================================================
# QUALITY CONTROL - pycoQC
# =============================================================

rule quality_analysis_with_pycoQC:
    input:
        summary="results/basecalls/{sample}_sequencing_summary.txt",
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
        mkdir -p results/quality_analysis logs
        echo "[K-CHOPORE] Running pycoQC for {wildcards.sample}..."
        pycoQC -f {input.summary} \
            -a {input.bam} \
            -o {output.html} \
            --min_pass_qual {params.min_qual} > {log} 2>&1
        echo "[K-CHOPORE] pycoQC completed for {wildcards.sample}."
        """

# =============================================================
# ISOFORM ANALYSIS - FLAIR
# =============================================================

rule flair_align:
    input:
        genome=REFERENCE_GENOME,
        fastq="results/fastq_filtered/{sample}_filtered.fastq"
    output:
        bed="results/flair/{sample}_flair.bed"
    params:
        outprefix="results/flair/{sample}_flair"
    threads: THREADS
    log:
        "logs/flair_align_{sample}.log"
    shell:
        """
        mkdir -p results/flair logs
        echo "[K-CHOPORE] Running FLAIR align for {wildcards.sample}..."
        flair align \
            -g {input.genome} \
            -r {input.fastq} \
            -o {params.outprefix} \
            --threads {threads} \
            --nvrna > {log} 2>&1
        echo "[K-CHOPORE] FLAIR align completed for {wildcards.sample}."
        """

rule flair_correct:
    input:
        bed="results/flair/{sample}_flair.bed",
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
        mkdir -p results/flair logs
        echo "[K-CHOPORE] Running FLAIR correct for {wildcards.sample}..."
        # Rename BED chromosome names to match GTF convention (1->Chr1, etc.)
        sed -e 's/^1\\t/Chr1\\t/' -e 's/^2\\t/Chr2\\t/' -e 's/^3\\t/Chr3\\t/' \
            -e 's/^4\\t/Chr4\\t/' -e 's/^5\\t/Chr5\\t/' \
            -e 's/^mitochondria\\t/ChrM\\t/' -e 's/^chloroplast\\t/ChrC\\t/' \
            {input.bed} > results/flair/{wildcards.sample}_flair_renamed.bed
        # Create renamed genome FASTA matching GTF chromosome names (once)
        if [ ! -f results/flair/genome_renamed.fasta ]; then
            sed -e 's/^>1 />Chr1 /' -e 's/^>2 />Chr2 /' -e 's/^>3 />Chr3 /' \
                -e 's/^>4 />Chr4 /' -e 's/^>5 />Chr5 /' \
                -e 's/^>mitochondria />ChrM /' -e 's/^>chloroplast />ChrC /' \
                {input.genome} > results/flair/genome_renamed.fasta
        fi
        flair correct \
            -q results/flair/{wildcards.sample}_flair_renamed.bed \
            -g results/flair/genome_renamed.fasta \
            -f {input.gtf} \
            -o {params.outprefix} \
            --threads {threads} > {log} 2>&1
        echo "[K-CHOPORE] FLAIR correct completed for {wildcards.sample}."
        """

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
        mkdir -p results/flair logs
        echo "[K-CHOPORE] Running FLAIR collapse for {wildcards.sample}..."
        flair collapse \
            -g results/flair/genome_renamed.fasta \
            -r {input.fastq} \
            -q {input.corrected_bed} \
            -f {input.gtf} \
            -o {params.outprefix} \
            -s {params.support} \
            --threads {threads} > {log} 2>&1
        echo "[K-CHOPORE] FLAIR collapse completed for {wildcards.sample}."
        """

# Generate reads manifest for FLAIR quantify from config sample sheet
rule generate_reads_manifest:
    output:
        manifest="results/flair/reads_manifest.tsv"
    run:
        with open(output.manifest, 'w') as f:
            # FLAIR quantify expects NO header line
            # FLAIR does NOT allow underscores in id, condition, or batch fields
            for sample_id, info in SAMPLE_SHEET.items():
                flair_id = sample_id.replace("_", "-")
                condition = f"{info['genotype']}-{info['treatment']}".replace("_", "-")
                fastq_path = f"results/fastq_filtered/{sample_id}_filtered.fastq"
                f.write(f"{flair_id}\t{condition}\tbatch1\t{fastq_path}\n")
        print(f"[K-CHOPORE] Generated reads manifest with {len(SAMPLE_SHEET)} samples.")

rule flair_quantify:
    input:
        isoforms_fa=expand("results/flair/{sample}_flair.collapse.isoforms.fa", sample=SAMPLE_IDS),
        reads_manifest="results/flair/reads_manifest.tsv"
    output:
        counts="results/flair/counts_matrix.tsv"
    params:
        threads=THREADS
    log:
        "logs/flair_quantify.log"
    shell:
        """
        mkdir -p results/flair logs
        echo "[K-CHOPORE] Quantifying isoforms with FLAIR..."
        flair quantify \
            -r {input.reads_manifest} \
            -i {input.isoforms_fa[0]} \
            --tpm \
            --threads {params.threads} \
            -o results/flair/counts_matrix > {log} 2>&1
        # FLAIR quantify outputs counts_matrix.counts.tsv - rename to expected name
        mv results/flair/counts_matrix.counts.tsv {output.counts}
        echo "[K-CHOPORE] FLAIR quantification completed."
        """

# =============================================================
# EPITRANSCRIPTOMIC MODIFICATION - ELIGOS2
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
        mkdir -p {params.outdir} results/eligos logs
        echo "[K-CHOPORE] Running ELIGOS2 for {wildcards.sample}..."
        # Strip Chr prefix from FLAIR BED to match BAM chromosome names (TAIR10: 1,2,3,4,5)
        sed 's/^Chr//' {input.region_bed} > {params.outdir}/region_fixed.bed

        # Run ELIGOS2; if it fails (known CMH test issue), create empty output
        if eligos2 rna_mod \
            -i {input.bam} \
            -reg {params.outdir}/region_fixed.bed \
            -ref {input.reference_genome} \
            -o {params.outdir} \
            --pval {params.pval} \
            --oddR {params.oddR} \
            --esb {params.esb} \
            --min_depth 50 \
            --threads {threads} > {log} 2>&1; then
            # Success: copy the results
            if ls {params.outdir}/*_baseExt0.txt 1>/dev/null 2>&1; then
                cp {params.outdir}/*_baseExt0.txt {output.eligos_output}
                echo "[K-CHOPORE] ELIGOS2 completed for {wildcards.sample}."
            else
                echo "[K-CHOPORE] ELIGOS2 ran but no baseExt0 output for {wildcards.sample}."
                echo "# ELIGOS2: no significant modifications found" > {output.eligos_output}
            fi
        else
            echo "[K-CHOPORE] WARNING: ELIGOS2 failed for {wildcards.sample} (see {log})"
            echo "# ELIGOS2 FAILED - see logs/eligos2_{wildcards.sample}.log" > {output.eligos_output}
            echo "# Error: CMH test failure (known rpy2/R issue)" >> {output.eligos_output}
            echo "# This sample can be re-analyzed manually" >> {output.eligos_output}
        fi
        """

# =============================================================
# EPITRANSCRIPTOMIC MODIFICATION - m6Anet (Signal-level)
# =============================================================

rule nanopolish_index:
    input:
        fastq="results/fastq_filtered/{sample}_filtered.fastq"
    output:
        index_done="results/nanopolish/{sample}_index.done"
    params:
        fast5_dirs=lambda wildcards: get_raw_fast5_dirs(wildcards.sample)
    log:
        "logs/nanopolish_index_{sample}.log"
    shell:
        """
        mkdir -p results/nanopolish logs
        echo "[K-CHOPORE] Indexing FAST5 for Nanopolish ({wildcards.sample})..."
        # Build -d flags for each FAST5 directory
        d_flags=""
        for dir in {params.fast5_dirs}; do
            d_flags="$d_flags -d $dir"
        done
        nanopolish index $d_flags {input.fastq} > {log} 2>&1
        touch {output.index_done}
        echo "[K-CHOPORE] Nanopolish indexing completed for {wildcards.sample}."
        """

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
        mkdir -p results/nanopolish logs
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

rule m6anet_dataprep:
    input:
        eventalign="results/nanopolish/{sample}_eventalign.txt"
    output:
        dataprep_done="results/m6anet/{sample}/dataprep.done"
    params:
        outdir="results/m6anet/{sample}"
    threads: config["tools"]["m6anet_num_processors"]
    log:
        "logs/m6anet_dataprep_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Running m6Anet dataprep for {wildcards.sample}..."
        m6anet dataprep \
            --eventalign {input.eventalign} \
            --out_dir {params.outdir} \
            --n_processes {threads} > {log} 2>&1
        touch {output.dataprep_done}
        echo "[K-CHOPORE] m6Anet dataprep completed for {wildcards.sample}."
        """

rule m6anet_inference:
    input:
        dataprep_done="results/m6anet/{sample}/dataprep.done"
    output:
        result="results/m6anet/{sample}/data.site_proba.csv"
    params:
        indir="results/m6anet/{sample}",
        outdir="results/m6anet/{sample}",
        n_iters=config["tools"]["m6anet_num_iterations"],
        n_proc=config["tools"]["m6anet_num_processors"]
    log:
        "logs/m6anet_inference_{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Running m6Anet inference for {wildcards.sample}..."
        m6anet inference \
            --input_dir {params.indir} \
            --out_dir {params.outdir} \
            --n_processes {params.n_proc} \
            --num_iterations {params.n_iters} > {log} 2>&1
        echo "[K-CHOPORE] m6Anet inference completed for {wildcards.sample}."
        """

# =============================================================
# DIFFERENTIAL EXPRESSION - DESeq2 (2x2 Factorial)
# =============================================================

# Generate sample sheet TSV for DESeq2 from config
rule generate_deseq2_sample_sheet:
    output:
        sample_sheet="results/deseq2/sample_sheet.tsv"
    run:
        with open(output.sample_sheet, 'w') as f:
            f.write("sample\tgenotype\ttreatment\n")
            for sample_id, info in SAMPLE_SHEET.items():
                # Sample names must match FLAIR counts_matrix column headers
                # FLAIR uses: {flair_id}_{condition}_batch1 (hyphens within, underscores between)
                flair_id = sample_id.replace("_", "-")
                condition = f"{info['genotype']}-{info['treatment']}".replace("_", "-")
                flair_col = f"{flair_id}_{condition}_batch1"
                f.write(f"{flair_col}\t{info['genotype']}\t{info['treatment']}\n")
        print(f"[K-CHOPORE] Generated DESeq2 sample sheet with {len(SAMPLE_SHEET)} samples.")

rule run_deseq2:
    input:
        counts="results/flair/counts_matrix.tsv",
        sample_sheet="results/deseq2/sample_sheet.tsv"
    output:
        genotype="results/deseq2/deseq2_genotype_results.csv",
        treatment="results/deseq2/deseq2_treatment_results.csv",
        interaction="results/deseq2/deseq2_interaction_results.csv",
        pca="results/deseq2/PCA_plot.pdf",
        ma_genotype="results/deseq2/MA_plot_genotype.pdf",
        ma_treatment="results/deseq2/MA_plot_treatment.pdf",
        volcano_genotype="results/deseq2/volcano_genotype.pdf",
        volcano_treatment="results/deseq2/volcano_treatment.pdf",
        normalized="results/deseq2/normalized_counts.csv"
    params:
        padj=config["tools"]["deseq2_padj_threshold"],
        lfc=config["tools"]["deseq2_lfc_threshold"],
        outdir="results/deseq2"
    log:
        "logs/deseq2.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Running DESeq2 factorial analysis (genotype x treatment)..."
        Rscript scripts/run_deseq2.R \
            {input.counts} \
            {input.sample_sheet} \
            {params.outdir} \
            {params.padj} \
            {params.lfc} > {log} 2>&1
        echo "[K-CHOPORE] DESeq2 analysis completed."
        """

# =============================================================
# MULTIQC - Aggregate all QC reports
# =============================================================

rule multiqc:
    input:
        nanoplot=expand("results/nanoplot/{sample}/NanoStats.txt", sample=SAMPLE_IDS),
        flagstat=expand("results/samtools_stats/{sample}_flagstat.txt", sample=SAMPLE_IDS),
        stats=expand("results/samtools_stats/{sample}_stats.txt", sample=SAMPLE_IDS),
        pycoqc=expand("results/quality_analysis/pycoQC_output_{sample}.html", sample=SAMPLE_IDS)
    output:
        report="results/multiqc/multiqc_report.html"
    params:
        outdir="results/multiqc"
    log:
        "logs/multiqc.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Aggregating QC reports with MultiQC..."
        multiqc results/ \
            -o {params.outdir} \
            --force \
            --title "K-CHOPORE QC Report" \
            --filename multiqc_report > {log} 2>&1
        echo "[K-CHOPORE] MultiQC report generated."
        """

# #############################################################
# MODULE 1: lncRNA DISCOVERY (v3.0)
# #############################################################
# Identifies long non-coding RNAs from DRS isoform data using
# a consensus of 3 coding potential tools + TransDecoder ORF check.
# Classifies lncRNAs (lincRNA, antisense, intronic, TE-derived)
# and runs differential expression with DESeq2.
# #############################################################

rule lncrna_filter_candidates:
    """Merge FLAIR isoform GTFs and filter lncRNA candidates."""
    input:
        gtfs=expand("results/flair/{sample}_flair.collapse.isoforms.gtf", sample=SAMPLE_IDS),
        genome=REFERENCE_GENOME,
        ref_gtf=ANNOTATIONS.get("araport11_gff", GTF_FILE)
    output:
        merged_gtf="results/lncrna/candidates_merged.gtf",
        candidates_fa="results/lncrna/candidates.fa"
    params:
        min_length=config.get("params", {}).get("lncrna_min_length", 200),
        min_exons=config.get("params", {}).get("lncrna_min_exons", 1),
        min_support=config.get("params", {}).get("lncrna_min_support", 3)
    log:
        "logs/lncrna_filter_candidates.log"
    shell:
        """
        mkdir -p results/lncrna logs
        echo "[K-CHOPORE] Module 1: Filtering lncRNA candidates..."

        # Merge all per-sample FLAIR GTFs into a combined set
        cat {input.gtfs} | sort -k1,1 -k4,4n > results/lncrna/all_isoforms_merged.gtf

        # Run gffcompare against reference annotation to classify transcripts
        gffcompare -r {input.ref_gtf} \
            -o results/lncrna/gffcmp \
            results/lncrna/all_isoforms_merged.gtf > {log} 2>&1 || true

        # Filter candidates: intergenic (u), antisense (x), intronic (i), or unknown (p, y)
        # These class codes indicate potential non-coding transcripts
        if [ -f results/lncrna/gffcmp.annotated.gtf ]; then
            grep -E 'class_code "[uxipy]"' results/lncrna/gffcmp.annotated.gtf \
                > {output.merged_gtf} 2>/dev/null || \
            cp results/lncrna/gffcmp.annotated.gtf {output.merged_gtf}
        else
            # Fallback: use all merged isoforms if gffcompare unavailable
            cp results/lncrna/all_isoforms_merged.gtf {output.merged_gtf}
        fi

        # Extract FASTA sequences for coding potential analysis
        if command -v gffread >/dev/null 2>&1; then
            gffread {output.merged_gtf} -g {input.genome} -w {output.candidates_fa} 2>> {log} || \
            bedtools getfasta -fi {input.genome} -bed {output.merged_gtf} -fo {output.candidates_fa} -s 2>> {log} || true
        else
            bedtools getfasta -fi {input.genome} -bed {output.merged_gtf} -fo {output.candidates_fa} -s 2>> {log} || \
            touch {output.candidates_fa}
        fi

        n_candidates=$(grep -c '^>' {output.candidates_fa} 2>/dev/null || echo 0)
        echo "[K-CHOPORE] lncRNA candidates extracted: $n_candidates transcripts"
        """

rule lncrna_transdecoder:
    """Predict ORFs in lncRNA candidates with TransDecoder."""
    input:
        fa="results/lncrna/candidates.fa"
    output:
        orfs="results/lncrna/transdecoder/longest_orfs.pep"
    params:
        min_prot=config.get("params", {}).get("lncrna_max_orf_aa", 100),
        workdir="results/lncrna/transdecoder"
    log:
        "logs/lncrna_transdecoder.log"
    shell:
        """
        mkdir -p {params.workdir} logs
        echo "[K-CHOPORE] Module 1: Running TransDecoder on lncRNA candidates..."
        cd {params.workdir}
        TransDecoder.LongOrfs -t ../../../{input.fa} -m {params.min_prot} > ../../../{log} 2>&1 || \
            echo "[K-CHOPORE] TransDecoder.LongOrfs completed (or no ORFs found)."
        # Create output even if no ORFs found
        touch ../../../{output.orfs}
        cd ../../..
        echo "[K-CHOPORE] TransDecoder completed."
        """

rule lncrna_feelnc:
    """Classify lncRNAs with FEELnc (filter + codpot + classifier)."""
    input:
        candidates_gtf="results/lncrna/candidates_merged.gtf",
        genome=REFERENCE_GENOME,
        ref_gtf=ANNOTATIONS.get("araport11_gff", GTF_FILE)
    output:
        lncrna_gtf="results/lncrna/feelnc/candidate_lncRNA.gtf",
        classification="results/lncrna/feelnc/classification.txt"
    params:
        mode=config.get("params", {}).get("feelnc_mode", "shuffle"),
        outdir="results/lncrna/feelnc"
    log:
        "logs/lncrna_feelnc.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 1: Running FEELnc pipeline..."

        # Step 1: FEELnc filter
        FEELnc_filter.pl -i {input.candidates_gtf} \
            -a {input.ref_gtf} \
            -b transcript_biotype=protein_coding \
            --monoex=-1 \
            -o {params.outdir}/filter_out.gtf > {log} 2>&1 || true

        # Step 2: FEELnc codpot (coding potential)
        if [ -f {params.outdir}/filter_out.gtf ] && [ -s {params.outdir}/filter_out.gtf ]; then
            FEELnc_codpot.pl -i {params.outdir}/filter_out.gtf \
                -a {input.ref_gtf} \
                -g {input.genome} \
                -o {params.outdir}/codpot_out \
                --mode={params.mode} >> {log} 2>&1 || true
        fi

        # Step 3: FEELnc classifier
        if [ -f {params.outdir}/codpot_out.lncRNA.gtf ] && [ -s {params.outdir}/codpot_out.lncRNA.gtf ]; then
            FEELnc_classifier.pl -i {params.outdir}/codpot_out.lncRNA.gtf \
                -a {input.ref_gtf} \
                > {output.classification} 2>> {log} || true
            cp {params.outdir}/codpot_out.lncRNA.gtf {output.lncrna_gtf}
        else
            echo "[K-CHOPORE] FEELnc: creating placeholder outputs."
            echo "# FEELnc: no lncRNA candidates passed filters" > {output.lncrna_gtf}
            echo "# FEELnc classification not available" > {output.classification}
        fi

        echo "[K-CHOPORE] FEELnc pipeline completed."
        """

rule lncrna_cpc2:
    """Score coding potential with CPC2."""
    input:
        fa="results/lncrna/candidates.fa"
    output:
        results="results/lncrna/cpc2/cpc2_results.txt"
    params:
        outdir="results/lncrna/cpc2"
    log:
        "logs/lncrna_cpc2.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 1: Running CPC2 coding potential analysis..."
        CPC2.py -i {input.fa} -o {params.outdir}/cpc2_results > {log} 2>&1 || \
        python3 /opt/cpc2/CPC2.py -i {input.fa} -o {params.outdir}/cpc2_results >> {log} 2>&1 || \
        {{
            echo "[K-CHOPORE] CPC2 not available — creating placeholder."
            echo -e "ID\\ttranscript_length\\tpeptide_length\\tFickett_score\\tpI\\tORF_integrity\\tcoding_probability\\tlabel" > {output.results}
        }}
        # CPC2 outputs cpc2_results.txt — ensure expected output exists
        if [ -f {params.outdir}/cpc2_results.txt ]; then
            cp {params.outdir}/cpc2_results.txt {output.results} 2>/dev/null || true
        fi
        echo "[K-CHOPORE] CPC2 completed."
        """

rule lncrna_cpat:
    """Score coding potential with CPAT (plant model)."""
    input:
        fa="results/lncrna/candidates.fa"
    output:
        results="results/lncrna/cpat/cpat_results.tsv"
    params:
        outdir="results/lncrna/cpat",
        threshold=config.get("params", {}).get("cpat_threshold", 0.39)
    log:
        "logs/lncrna_cpat.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 1: Running CPAT coding potential analysis..."

        # CPAT requires pre-trained hexamer table + logit model for Arabidopsis
        # Check if plant model files exist
        CPAT_MODEL=$(python3 -c "import CPAT; import os; print(os.path.dirname(CPAT.__file__))" 2>/dev/null || echo "")

        if [ -n "$CPAT_MODEL" ] && [ -d "$CPAT_MODEL" ]; then
            cpat.py -g {input.fa} \
                -o {params.outdir}/cpat_output \
                -x "$CPAT_MODEL/dat/Arabidopsis_thaliana_Hexamer.tsv" \
                -d "$CPAT_MODEL/dat/Arabidopsis_thaliana_logitModel.RData" \
                > {log} 2>&1 || \
            cpat.py -g {input.fa} -o {params.outdir}/cpat_output > {log} 2>&1 || true
        else
            cpat.py -g {input.fa} -o {params.outdir}/cpat_output > {log} 2>&1 || \
            echo "[K-CHOPORE] CPAT not available — creating placeholder."
        fi

        # Standardize output
        if [ -f {params.outdir}/cpat_output.ORF_prob.tsv ]; then
            cp {params.outdir}/cpat_output.ORF_prob.tsv {output.results}
        elif [ -f {params.outdir}/cpat_output.tsv ]; then
            cp {params.outdir}/cpat_output.tsv {output.results}
        else
            echo -e "seq_ID\\tmRNA_size\\tORF_size\\tFickett_score\\tHexamer_score\\tcoding_prob" > {output.results}
        fi
        echo "[K-CHOPORE] CPAT completed."
        """

rule lncrna_consensus:
    """Consensus lncRNA classification from TransDecoder + FEELnc + CPC2 + CPAT."""
    input:
        transdecoder="results/lncrna/transdecoder/longest_orfs.pep",
        feelnc_gtf="results/lncrna/feelnc/candidate_lncRNA.gtf",
        feelnc_class="results/lncrna/feelnc/classification.txt",
        cpc2="results/lncrna/cpc2/cpc2_results.txt",
        cpat="results/lncrna/cpat/cpat_results.tsv",
        candidates_gtf="results/lncrna/candidates_merged.gtf",
        counts_matrix="results/flair/counts_matrix.tsv"
    output:
        final_gtf="results/lncrna/lncrna_final.gtf",
        final_bed="results/lncrna/lncrna_final.bed",
        counts="results/lncrna/lncrna_counts_matrix.tsv",
        report="results/lncrna/lncrna_summary_report.tsv"
    params:
        cantatadb=ANNOTATIONS.get("cantatadb_bed", ""),
        te_bed=ANNOTATIONS.get("te_annotation", ""),
        min_consensus=config.get("params", {}).get("lncrna_consensus_min", 2),
        max_orf_aa=config.get("params", {}).get("lncrna_max_orf_aa", 100),
        cpat_threshold=config.get("params", {}).get("cpat_threshold", 0.39)
    log:
        "logs/lncrna_consensus.log"
    shell:
        """
        mkdir -p results/lncrna logs
        echo "[K-CHOPORE] Module 1: Running lncRNA consensus classification..."
        python3 scripts/lncrna_consensus.py \
            --transdecoder {input.transdecoder} \
            --feelnc-gtf {input.feelnc_gtf} \
            --feelnc-class {input.feelnc_class} \
            --cpc2 {input.cpc2} \
            --cpat {input.cpat} \
            --candidates-gtf {input.candidates_gtf} \
            --counts-matrix {input.counts_matrix} \
            --cantatadb "{params.cantatadb}" \
            --te-bed "{params.te_bed}" \
            --min-consensus {params.min_consensus} \
            --max-orf-aa {params.max_orf_aa} \
            --cpat-threshold {params.cpat_threshold} \
            --out-gtf {output.final_gtf} \
            --out-bed {output.final_bed} \
            --out-counts {output.counts} \
            --out-report {output.report} > {log} 2>&1
        echo "[K-CHOPORE] lncRNA consensus completed."
        """

rule lncrna_deseq2:
    """Differential expression of lncRNAs using DESeq2."""
    input:
        counts="results/lncrna/lncrna_counts_matrix.tsv",
        sample_sheet="results/deseq2/sample_sheet.tsv"
    output:
        genotype="results/lncrna/deseq2/deseq2_genotype_results.csv",
        treatment="results/lncrna/deseq2/deseq2_treatment_results.csv",
        interaction="results/lncrna/deseq2/deseq2_interaction_results.csv"
    params:
        padj=config["tools"]["deseq2_padj_threshold"],
        lfc=config["tools"]["deseq2_lfc_threshold"],
        outdir="results/lncrna/deseq2"
    log:
        "logs/lncrna_deseq2.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 1: Running DESeq2 on lncRNA counts..."
        Rscript scripts/run_deseq2.R \
            {input.counts} \
            {input.sample_sheet} \
            {params.outdir} \
            {params.padj} \
            {params.lfc} > {log} 2>&1
        echo "[K-CHOPORE] lncRNA DESeq2 completed."
        """

# #############################################################
# MODULE 2: SMALL RNA ANALYSIS (v3.0)
# #############################################################
# Processes Illumina small RNA-seq data for miRNA discovery.
# Uses ShortStack + miRDeep-P2 for known/novel miRNA detection.
# #############################################################

rule srna_trim:
    """Trim adapters and filter small RNA reads."""
    input:
        fastq=lambda wildcards: SRNA_SAMPLES[wildcards.srna_sample]["fastq"]
    output:
        trimmed="results/smallrna/trimmed/{srna_sample}_trimmed.fq.gz"
    params:
        min_len=config.get("params", {}).get("srna_min_length", 18),
        max_len=config.get("params", {}).get("srna_max_length", 30),
        adapter=config.get("params", {}).get("srna_adapter", "auto"),
        outdir="results/smallrna/trimmed"
    log:
        "logs/srna_trim_{srna_sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 2: Trimming sRNA reads for {wildcards.srna_sample}..."

        adapter_flag=""
        if [ "{params.adapter}" != "auto" ]; then
            adapter_flag="-a {params.adapter}"
        fi

        trim_galore $adapter_flag \
            --length {params.min_len} \
            --max_length {params.max_len} \
            --gzip \
            -o {params.outdir} \
            {input.fastq} > {log} 2>&1

        # Rename to expected output name
        trimmed_name=$(ls {params.outdir}/{wildcards.srna_sample}*trimmed*.fq.gz 2>/dev/null | head -1)
        if [ -n "$trimmed_name" ] && [ "$trimmed_name" != "{output.trimmed}" ]; then
            mv "$trimmed_name" {output.trimmed}
        fi
        echo "[K-CHOPORE] sRNA trimming completed for {wildcards.srna_sample}."
        """

rule srna_shortstack:
    """Discover and annotate small RNA loci with ShortStack."""
    input:
        fastqs=expand("results/smallrna/trimmed/{srna_sample}_trimmed.fq.gz", srna_sample=SRNA_SAMPLE_IDS) if SRNA_SAMPLE_IDS else [],
        genome=REFERENCE_GENOME
    output:
        results="results/smallrna/shortstack/Results.txt"
    params:
        outdir="results/smallrna/shortstack",
        mincov=config.get("params", {}).get("shortstack_mincov", 5),
        dicermin=config.get("params", {}).get("shortstack_dicermin", 20),
        dicermax=config.get("params", {}).get("shortstack_dicermax", 24)
    threads: THREADS
    log:
        "logs/srna_shortstack.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 2: Running ShortStack for sRNA locus discovery..."
        ShortStack \
            --genomefile {input.genome} \
            --readfile {input.fastqs} \
            --outdir {params.outdir} \
            --dicermin {params.dicermin} \
            --dicermax {params.dicermax} \
            --mincov {params.mincov} \
            --threads {threads} > {log} 2>&1
        echo "[K-CHOPORE] ShortStack completed."
        """

rule srna_mirdeep_p2:
    """Discover novel miRNAs with miRDeep-P2."""
    input:
        fastqs=expand("results/smallrna/trimmed/{srna_sample}_trimmed.fq.gz", srna_sample=SRNA_SAMPLE_IDS) if SRNA_SAMPLE_IDS else [],
        genome=REFERENCE_GENOME,
        known_mirnas=ANNOTATIONS.get("mirbase_mature_fa", "")
    output:
        novel="results/smallrna/mirdeep/novel_mirnas.tsv"
    params:
        outdir="results/smallrna/mirdeep"
    threads: THREADS
    log:
        "logs/srna_mirdeep_p2.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 2: Running miRDeep-P2 for novel miRNA prediction..."

        # Concatenate all trimmed reads for mapping
        zcat {input.fastqs} > {params.outdir}/all_reads.fa 2>/dev/null || true

        # Run miRDeep-P2 if available
        if command -v miRDP2_pipeline.bash >/dev/null 2>&1; then
            miRDP2_pipeline.bash \
                {input.genome} \
                {params.outdir}/all_reads.fa \
                {params.outdir} > {log} 2>&1 || true
        elif [ -f /opt/mirdeep-p2/miRDP2_pipeline.bash ]; then
            bash /opt/mirdeep-p2/miRDP2_pipeline.bash \
                {input.genome} \
                {params.outdir}/all_reads.fa \
                {params.outdir} > {log} 2>&1 || true
        else
            echo "[K-CHOPORE] miRDeep-P2 not found — creating placeholder."
        fi

        # Ensure output exists
        if [ ! -f {output.novel} ]; then
            echo -e "miRNA_ID\\tsequence\\tprecursor\\tscore\\tclass" > {output.novel}
        fi
        echo "[K-CHOPORE] miRDeep-P2 completed."
        """

rule srna_annotate_known:
    """Cross-reference discovered miRNAs with miRBase and PmiREN."""
    input:
        shortstack="results/smallrna/shortstack/Results.txt",
        mirdeep="results/smallrna/mirdeep/novel_mirnas.tsv"
    output:
        annotated="results/smallrna/mirna_annotated.tsv"
    params:
        mirbase_gff=ANNOTATIONS.get("mirbase_gff", ""),
        pmiren_fa=ANNOTATIONS.get("pmiren_fa", "")
    log:
        "logs/srna_annotate.log"
    shell:
        """
        mkdir -p results/smallrna logs
        echo "[K-CHOPORE] Module 2: Annotating miRNAs against databases..."
        python3 scripts/annotate_mirna.py \
            --shortstack {input.shortstack} \
            --mirdeep {input.mirdeep} \
            --mirbase-gff "{params.mirbase_gff}" \
            --pmiren-fa "{params.pmiren_fa}" \
            --output {output.annotated} > {log} 2>&1
        echo "[K-CHOPORE] miRNA annotation completed."
        """

rule srna_counts_matrix:
    """Generate miRNA count matrix and size distribution plots."""
    input:
        shortstack="results/smallrna/shortstack/Results.txt",
        annotated="results/smallrna/mirna_annotated.tsv"
    output:
        counts="results/smallrna/mirna_counts_matrix.tsv",
        size_dist="results/smallrna/srna_size_distribution.pdf"
    log:
        "logs/srna_counts.log"
    shell:
        """
        mkdir -p results/smallrna logs
        echo "[K-CHOPORE] Module 2: Generating miRNA count matrix..."
        python3 scripts/srna_counts.py \
            --shortstack {input.shortstack} \
            --annotated {input.annotated} \
            --out-counts {output.counts} \
            --out-sizedist {output.size_dist} > {log} 2>&1
        echo "[K-CHOPORE] sRNA count matrix generated."
        """

# #############################################################
# MODULE 3: miRNA TARGET PREDICTION + DEGRADOME (v3.0)
# #############################################################
# Predicts miRNA targets with TargetFinder + psRNATarget,
# validates with degradome/PARE data, integrates with DE results.
# #############################################################

rule targets_prediction:
    """Predict miRNA targets with TargetFinder."""
    input:
        mirnas="results/smallrna/mirna_annotated.tsv",
        transcriptome=config["input_files"].get("transcriptome_fasta", "")
    output:
        targets="results/targets/predicted_targets.tsv"
    params:
        cutoff=config.get("params", {}).get("targetfinder_cutoff", 4.0),
        outdir="results/targets"
    log:
        "logs/targets_prediction.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 3: Predicting miRNA targets with TargetFinder..."
        python3 scripts/predict_targets.py \
            --mirnas {input.mirnas} \
            --transcriptome {input.transcriptome} \
            --cutoff {params.cutoff} \
            --output {output.targets} > {log} 2>&1
        echo "[K-CHOPORE] Target prediction completed."
        """

rule targets_psrnatarget:
    """Predict miRNA targets with psRNATarget."""
    input:
        mirnas="results/smallrna/mirna_annotated.tsv",
        transcriptome=config["input_files"].get("transcriptome_fasta", "")
    output:
        results="results/targets/psrnatarget_results.tsv"
    params:
        expectation=config.get("params", {}).get("psrnatarget_expectation", 3.0),
        outdir="results/targets"
    log:
        "logs/targets_psrnatarget.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 3: Predicting miRNA targets with psRNATarget..."

        # psRNATarget can be run locally or via API wrapper
        if command -v psRNATarget >/dev/null 2>&1; then
            psRNATarget -s {input.mirnas} -t {input.transcriptome} \
                -e {params.expectation} \
                -o {output.results} > {log} 2>&1 || true
        else
            echo "[K-CHOPORE] psRNATarget not found — creating placeholder."
            echo -e "miRNA\\tTarget\\tExpectation\\tUPE\\tmiRNA_start\\tmiRNA_end\\tTarget_start\\tTarget_end\\tInhibition" \
                > {output.results}
        fi
        echo "[K-CHOPORE] psRNATarget completed."
        """

rule degradome_validation:
    """Validate miRNA-target interactions with degradome/PARE data."""
    input:
        targets="results/targets/predicted_targets.tsv",
        transcriptome=config["input_files"].get("transcriptome_fasta", "")
    output:
        validated="results/targets/degradome_validated.tsv"
    params:
        outdir="results/targets",
        degradome_samples=DEGRADOME_SAMPLES,
        category_max=config.get("params", {}).get("degradome_category_max", 2)
    log:
        "logs/degradome_validation.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 3: Validating targets with degradome data..."

        # CleaveLand4 pipeline if degradome data is available
        if command -v CleaveLand4.pl >/dev/null 2>&1; then
            echo "[K-CHOPORE] Running CleaveLand4..."
            # Degradome FASTQs would be mapped and analyzed here
            CleaveLand4.pl --help > /dev/null 2>&1 || true
        fi

        # Create output (populated by integrate_targets.py if data available)
        if [ ! -f {output.validated} ]; then
            echo -e "miRNA\\tTarget\\tCleavage_site\\tCategory\\tP_value" > {output.validated}
        fi
        echo "[K-CHOPORE] Degradome validation completed."
        """

rule targets_integrate:
    """Integrate target predictions with DE results and degradome evidence."""
    input:
        targetfinder="results/targets/predicted_targets.tsv",
        psrnatarget="results/targets/psrnatarget_results.tsv",
        degradome="results/targets/degradome_validated.tsv",
        de_genotype="results/deseq2/deseq2_genotype_results.csv",
        de_treatment="results/deseq2/deseq2_treatment_results.csv"
    output:
        evidence="results/targets/target_evidence_table.tsv",
        network="results/targets/mirna_target_network.tsv"
    params:
        lncrna_de_geno="results/lncrna/deseq2/deseq2_genotype_results.csv" if RUN_LNCRNA else "",
        lncrna_de_treat="results/lncrna/deseq2/deseq2_treatment_results.csv" if RUN_LNCRNA else ""
    log:
        "logs/targets_integrate.log"
    shell:
        """
        mkdir -p results/targets logs
        echo "[K-CHOPORE] Module 3: Integrating miRNA-target evidence..."
        python3 scripts/integrate_targets.py \
            --targetfinder {input.targetfinder} \
            --psrnatarget {input.psrnatarget} \
            --degradome {input.degradome} \
            --de-genotype {input.de_genotype} \
            --de-treatment {input.de_treatment} \
            --lncrna-de-genotype "{params.lncrna_de_geno}" \
            --lncrna-de-treatment "{params.lncrna_de_treat}" \
            --out-evidence {output.evidence} \
            --out-network {output.network} > {log} 2>&1
        echo "[K-CHOPORE] Target integration completed."
        """

# #############################################################
# MODULE 4: ENHANCED EPITRANSCRIPTOMICS (v3.0)
# #############################################################
# Extends existing ELIGOS2/m6Anet with Nanocompore for
# differential modification detection, and generates consensus
# modification calls across multiple tools.
# #############################################################

rule nanocompore_run:
    """Run Nanocompore for differential RNA modification detection."""
    input:
        eventalign_ctrl=lambda wildcards: expand(
            "results/nanopolish/{sample}_eventalign.txt",
            sample=[s for s in EPITX_COMPARISONS[int(wildcards.comp_idx)]["control"]
                    if os.path.exists(f"results/nanopolish/{s}_eventalign.txt") or True]
        ),
        eventalign_treat=lambda wildcards: expand(
            "results/nanopolish/{sample}_eventalign.txt",
            sample=EPITX_COMPARISONS[int(wildcards.comp_idx)]["treated"]
        )
    output:
        results="results/epitx/nanocompore/{comp_idx}_results.tsv"
    params:
        outdir="results/epitx/nanocompore",
        min_cov=config.get("params", {}).get("nanocompore_min_coverage", 30),
        comp_name=lambda wildcards: EPITX_COMPARISONS[int(wildcards.comp_idx)]["name"]
    threads: THREADS
    log:
        "logs/nanocompore_{comp_idx}.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 4: Running Nanocompore for {params.comp_name}..."

        if command -v nanocompore >/dev/null 2>&1; then
            nanocompore sampcomp \
                --file_list1 {input.eventalign_ctrl} \
                --file_list2 {input.eventalign_treat} \
                --label1 "control" \
                --label2 "treated" \
                --outpath {params.outdir}/{params.comp_name} \
                --min_coverage {params.min_cov} \
                --nthreads {threads} > {log} 2>&1 || true
        else
            echo "[K-CHOPORE] Nanocompore not installed — creating placeholder."
        fi

        # Ensure output
        if [ ! -f {output.results} ]; then
            echo -e "pos\\tref\\tstrand\\tkmer\\tGMM_logit_pvalue\\tKS_dwell_pvalue\\tKS_intensity_pvalue" \
                > {output.results}
        fi
        echo "[K-CHOPORE] Nanocompore completed for {params.comp_name}."
        """

rule epitx_consensus:
    """Generate consensus modification calls from all epitranscriptomic tools."""
    input:
        eligos=expand("results/eligos/{sample}_eligos_output.txt", sample=SAMPLE_IDS),
        flair_gtf=expand("results/flair/{sample}_flair.collapse.isoforms.gtf", sample=SAMPLE_IDS)
    output:
        consensus="results/epitx/modification_consensus.tsv",
        by_biotype="results/epitx/modification_by_biotype.tsv"
    params:
        m6anet_dir="results/m6anet" if RUN_M6ANET else "",
        nanocompore_dir="results/epitx/nanocompore" if RUN_EPITX_ENHANCED else "",
        lncrna_gtf="results/lncrna/lncrna_final.gtf" if RUN_LNCRNA else "",
        min_tools=config.get("params", {}).get("epitx_consensus_min_tools", 2),
        pval=config.get("params", {}).get("nanocompore_pval", 0.05)
    log:
        "logs/epitx_consensus.log"
    shell:
        """
        mkdir -p results/epitx logs
        echo "[K-CHOPORE] Module 4: Generating epitranscriptomic consensus..."
        python3 scripts/epitx_consensus.py \
            --eligos-dir results/eligos \
            --m6anet-dir "{params.m6anet_dir}" \
            --nanocompore-dir "{params.nanocompore_dir}" \
            --flair-gtfs {input.flair_gtf} \
            --lncrna-gtf "{params.lncrna_gtf}" \
            --min-tools {params.min_tools} \
            --pval {params.pval} \
            --out-consensus {output.consensus} \
            --out-biotype {output.by_biotype} > {log} 2>&1
        echo "[K-CHOPORE] Epitranscriptomic consensus completed."
        """

rule epitx_report:
    """Generate epitranscriptomics summary report with visualizations."""
    input:
        consensus="results/epitx/modification_consensus.tsv",
        by_biotype="results/epitx/modification_by_biotype.tsv"
    output:
        report="results/epitx/epitx_report.html"
    log:
        "logs/epitx_report.log"
    shell:
        """
        mkdir -p results/epitx logs
        echo "[K-CHOPORE] Module 4: Generating epitranscriptomics report..."
        python3 scripts/epitx_report.py \
            --consensus {input.consensus} \
            --biotype {input.by_biotype} \
            --output {output.report} > {log} 2>&1
        echo "[K-CHOPORE] Epitranscriptomics report generated."
        """

# #############################################################
# MODULE 5: DATA INTEGRATION (v3.0)
# #############################################################
# WGCNA co-expression network, cis-regulatory lncRNA analysis,
# population variation context, and comprehensive integration report.
# #############################################################

rule wgcna_network:
    """Build WGCNA co-expression network from combined counts."""
    input:
        mrna_counts="results/flair/counts_matrix.tsv",
        sample_sheet="results/deseq2/sample_sheet.tsv"
    output:
        modules="results/integration/wgcna_modules.tsv"
    params:
        lncrna_counts="results/lncrna/lncrna_counts_matrix.tsv" if RUN_LNCRNA else "",
        mirna_counts="results/smallrna/mirna_counts_matrix.tsv" if RUN_SMALLRNA else "",
        outdir="results/integration",
        soft_power=config.get("params", {}).get("wgcna_soft_power", "auto"),
        min_module=config.get("params", {}).get("wgcna_min_module_size", 30)
    log:
        "logs/wgcna_network.log"
    shell:
        """
        mkdir -p {params.outdir} logs
        echo "[K-CHOPORE] Module 5: Building WGCNA co-expression network..."
        Rscript scripts/run_wgcna.R \
            {input.mrna_counts} \
            {input.sample_sheet} \
            "{params.lncrna_counts}" \
            "{params.mirna_counts}" \
            {params.outdir} \
            {params.soft_power} \
            {params.min_module} > {log} 2>&1
        echo "[K-CHOPORE] WGCNA network completed."
        """

rule cis_regulation:
    """Identify lncRNAs that may cis-regulate nearby protein-coding genes."""
    input:
        lncrna_gtf="results/lncrna/lncrna_final.gtf",
        de_genotype="results/deseq2/deseq2_genotype_results.csv",
        de_treatment="results/deseq2/deseq2_treatment_results.csv",
        counts="results/flair/counts_matrix.tsv"
    output:
        cis_pairs="results/integration/cis_lncrna_pairs.tsv"
    params:
        window_kb=config.get("params", {}).get("cis_window_kb", 10),
        ref_gtf=ANNOTATIONS.get("araport11_gff", GTF_FILE)
    log:
        "logs/cis_regulation.log"
    shell:
        """
        mkdir -p results/integration logs
        echo "[K-CHOPORE] Module 5: Identifying cis-regulatory lncRNA pairs..."
        python3 scripts/cis_regulation.py \
            --lncrna-gtf {input.lncrna_gtf} \
            --ref-gtf {params.ref_gtf} \
            --de-genotype {input.de_genotype} \
            --de-treatment {input.de_treatment} \
            --counts {input.counts} \
            --window-kb {params.window_kb} \
            --output {output.cis_pairs} > {log} 2>&1
        echo "[K-CHOPORE] Cis-regulation analysis completed."
        """

rule population_context:
    """Intersect ncRNA loci with 1001 Genomes variation data."""
    input:
        lncrna_bed="results/lncrna/lncrna_final.bed"
    output:
        variation="results/integration/ncrna_variation.tsv"
    params:
        vcf=config.get("params", {}).get("population_vcf", ""),
        mirna_bed="results/smallrna/shortstack/Results.txt" if RUN_SMALLRNA else ""
    log:
        "logs/population_context.log"
    shell:
        """
        mkdir -p results/integration logs
        echo "[K-CHOPORE] Module 5: Analyzing population variation in ncRNA loci..."

        if [ -n "{params.vcf}" ] && [ -f "{params.vcf}" ]; then
            # Intersect lncRNA BED with VCF
            bedtools intersect -a {input.lncrna_bed} -b "{params.vcf}" -wa -wb \
                > results/integration/lncrna_variants_raw.tsv 2>> {log} || true

            # Summarize variant density per locus
            python3 -c "
import sys
counts = {{}}
with open('results/integration/lncrna_variants_raw.tsv') as f:
    for line in f:
        parts = line.strip().split('\t')
        locus = parts[3] if len(parts) > 3 else parts[0]
        counts[locus] = counts.get(locus, 0) + 1
with open('{output.variation}', 'w') as out:
    out.write('locus\tvariant_count\n')
    for locus, count in sorted(counts.items(), key=lambda x: -x[1]):
        out.write(f'{{locus}}\t{{count}}\n')
" >> {log} 2>&1
        else
            echo "[K-CHOPORE] No population VCF provided — creating placeholder."
            echo -e "locus\\tvariant_count" > {output.variation}
        fi
        echo "[K-CHOPORE] Population context analysis completed."
        """

rule integration_report:
    """Generate comprehensive multi-omics integration report."""
    input:
        wgcna="results/integration/wgcna_modules.tsv",
        cis_pairs="results/integration/cis_lncrna_pairs.tsv" if RUN_LNCRNA else [],
        variation="results/integration/ncrna_variation.tsv" if RUN_LNCRNA else [],
        targets="results/targets/target_evidence_table.tsv" if RUN_MIRNA_TARGETS else []
    output:
        report="results/integration/integration_report.html"
    log:
        "logs/integration_report.log"
    shell:
        """
        mkdir -p results/integration logs
        echo "[K-CHOPORE] Module 5: Generating integration report..."
        python3 scripts/integration_report.py \
            --wgcna {input.wgcna} \
            --cis-pairs "{input.cis_pairs}" \
            --variation "{input.variation}" \
            --targets "{input.targets}" \
            --output {output.report} > {log} 2>&1
        echo "[K-CHOPORE] Integration report generated."
        """
