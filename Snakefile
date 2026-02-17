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
