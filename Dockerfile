# ==============================================================
# K-CHOPORE Dockerfile
# Comprehensive ONT Direct RNA Sequencing Analysis Pipeline
# ==============================================================
# This Dockerfile builds the complete K-CHOPORE environment with
# all tools for basecalling, QC, alignment, isoform analysis,
# epitranscriptomic modification detection, and differential expression.
#
# Base image: Ubuntu 22.04
# Maintained by: Pelayo Gonzalez de Lena Rodriguez
# ==============================================================

FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Ensure 'python' command exists (many tools use #!/usr/bin/env python)
RUN ln -sf /usr/bin/python3 /usr/bin/python

LABEL maintainer="pelayovic"
LABEL description="K-CHOPORE: ONT Direct RNA-seq Analysis Pipeline"

# ==============================================================
# 1. System dependencies and essential tools
# ==============================================================
RUN apt-get update && apt-get install -y \
    curl wget git \
    jq \
    build-essential \
    gcc g++ \
    libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev \
    liblzma-dev libssl-dev libffi-dev libcurl4-openssl-dev \
    python3 python3-pip python3-venv python3-dev \
    default-jdk default-jre \
    r-base \
    bedtools \
    autoconf automake \
    libhdf5-dev \
    libxml2-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Install yq for YAML parsing
RUN curl -L https://github.com/mikefarah/yq/releases/download/v4.9.8/yq_linux_amd64 \
    -o /usr/bin/yq && chmod +x /usr/bin/yq

# ==============================================================
# 2. Samtools (alignment stats, sort, index, flagstat)
# ==============================================================
RUN wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 && \
    tar -xjf samtools-1.19.tar.bz2 && \
    cd samtools-1.19 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf samtools-1.19*

# ==============================================================
# 3. Minimap2 (splice-aware long-read aligner)
# ==============================================================
RUN apt-get update && apt-get install -y minimap2 && rm -rf /var/lib/apt/lists/*

# ==============================================================
# 4. Picard (BAM validation and manipulation)
# ==============================================================
RUN wget https://github.com/broadinstitute/picard/releases/download/2.25.7/picard.jar \
    -P /usr/local/bin/

# ==============================================================
# 5. Dorado (modern ONT basecaller for direct RNA)
# ==============================================================
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.0-linux-x64.tar.gz -P /opt/ && \
    tar -xzvf /opt/dorado-0.8.0-linux-x64.tar.gz -C /opt/ && \
    rm /opt/dorado-0.8.0-linux-x64.tar.gz && \
    ln -s /opt/dorado-0.8.0-linux-x64/bin/dorado /usr/local/bin/dorado
ENV PATH="/opt/dorado-0.8.0-linux-x64/bin:$PATH"

# ==============================================================
# 6. Guppy (legacy ONT basecaller, CPU version)
# ==============================================================
RUN rm -rf /opt/ont-guppy && \
    wget --no-check-certificate \
    https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.1.5_linux64.tar.gz \
    -O /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz && \
    tar -xvzf /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz -C /opt/ && \
    rm /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz
ENV PATH="/opt/ont-guppy-cpu/bin:$PATH"

# ==============================================================
# 7. StringTie2 (long-read isoform assembly)
# ==============================================================
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz && \
    tar -xzvf stringtie-2.2.1.Linux_x86_64.tar.gz && \
    cp stringtie-2.2.1.Linux_x86_64/stringtie /usr/local/bin/ && \
    rm -rf stringtie-2.2.1*

# ==============================================================
# 8. Nanopolish (signal-level analysis for m6Anet/xPore)
# ==============================================================
# Nanopolish must be compiled from source (no pip package available)
# Use serial build (-j1) to avoid race conditions between HDF5/eigen
# downloads and source compilation
RUN git clone --recursive https://github.com/jts/nanopolish.git /opt/nanopolish && \
    cd /opt/nanopolish && \
    make -j1 && \
    ln -s /opt/nanopolish/nanopolish /usr/local/bin/nanopolish && \
    ln -s /opt/nanopolish/scripts/nanopolish_makerange.py /usr/local/bin/nanopolish_makerange.py

# VBZ compression plugin for HDF5 (required to read VBZ-compressed FAST5 files)
# The plugin is already bundled with Guppy at /opt/ont-guppy-cpu/lib/
# Point HDF5_PLUGIN_PATH there so nanopolish (and its bundled HDF5) can find it
ENV HDF5_PLUGIN_PATH="/opt/ont-guppy-cpu/lib"

# ==============================================================
# 9. Python package manager and Snakemake
# ==============================================================
RUN pip install --upgrade pip && \
    pip install pulp==2.7.0 && \
    pip install --upgrade snakemake && \
    sed -i 's/list_solvers/listSolvers/' /usr/local/lib/python3.10/dist-packages/snakemake/__init__.py || true

# ==============================================================
# 10. Bonito (alternative ONT basecaller)
# ==============================================================
RUN pip install ont-bonito || true
# NOTE: Bonito models are large (~5GB). Download at runtime instead:
#   docker run -v bonito_models:/root/.local/share/bonito k-chopore bonito download --models

ENV PATH="/usr/local/bin/bonito:$PATH"

# ==============================================================
# 11. Python requirements (core pipeline)
# ==============================================================
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir setuptools wheel && \
    pip install --no-cache-dir -r /tmp/requirements.txt

# ==============================================================
# 12. NanoPlot, NanoComp, NanoFilt (ONT-specific QC tools)
# ==============================================================
RUN pip install --no-cache-dir \
    NanoPlot \
    NanoComp \
    NanoFilt \
    nanoget \
    nanomath

# pycoQC installed with --no-deps to avoid plotly version conflict
# (NanoPlot requires plotly>=6.1.1, pycoQC pins plotly==4.1.0)
# setuptools must be reinstalled after NanoPlot (which may remove it)
# because pycoQC imports pkg_resources at runtime
RUN pip install --no-cache-dir setuptools && \
    pip install --no-cache-dir --no-deps pycoQC

# ==============================================================
# 13. m6Anet (m6A modification detection from nanopore signal)
# ==============================================================
# m6anet 2.0.1 pins numpy==1.18.0 which is incompatible with Python 3.10.
# Install with --no-deps and provide compatible dependencies manually.
RUN pip install --no-cache-dir --no-deps m6anet && \
    pip install --no-cache-dir torch scikit-learn

# ==============================================================
# 14. xPore (differential RNA modification from nanopore)
# ==============================================================
RUN pip install --no-cache-dir xpore || \
    pip install --no-cache-dir --no-deps xpore

# ==============================================================
# 15. ELIGOS2 (RNA modification detection from error signatures)
# ==============================================================
COPY scripts/fix_eligos2.py /tmp/fix_eligos2.py
RUN git clone https://gitlab.com/piroonj/eligos2.git /home/eligos2 && \
    cd /home/eligos2 && \
    pip install --no-cache-dir -r requirements.txt 2>/dev/null || true && \
    chmod +x /home/eligos2/Scripts/*.py 2>/dev/null || true && \
    # Fix rpy2 DeprecationWarning + BED merge strand column bug \
    python3 /tmp/fix_eligos2.py

# ==============================================================
# 16. MultiQC (aggregate QC report)
# ==============================================================
RUN pip install --no-cache-dir multiqc

# ==============================================================
# 16. POD5 tools (modern ONT file format support)
# ==============================================================
RUN pip install --no-cache-dir pod5

# ==============================================================
# 17. Additional Python packages
# ==============================================================
RUN pip install --no-cache-dir \
    pysam \
    rpy2 \
    tabulate==0.9.0 \
    ont-fast5-api

# ==============================================================
# 18. R packages for differential expression
# ==============================================================
RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")' && \
    Rscript -e 'BiocManager::install(c("DESeq2", "tximport", "apeglm"), ask=FALSE)' && \
    Rscript -e 'install.packages(c("samplesizeCMH", "ggplot2", "pheatmap", "EnhancedVolcano"), repos="https://cloud.r-project.org")' || true

# ==============================================================
# 19. FLAIR (isoform analysis for long reads)
# ==============================================================
RUN pip install --no-cache-dir flair-brookslab

# ==============================================================
# v3.0 MODULE 1: lncRNA Discovery Tools
# ==============================================================

# gffcompare + gffread (transcript classification + sequence extraction)
RUN apt-get update && apt-get install -y gffread 2>/dev/null || \
    (git clone https://github.com/gpertea/gffread.git /opt/gffread && \
     cd /opt/gffread && make release && cp gffread /usr/local/bin/) && \
    (git clone https://github.com/gpertea/gffcompare.git /opt/gffcompare && \
     cd /opt/gffcompare && make release && cp gffcompare /usr/local/bin/) || true

# TransDecoder (ORF prediction to filter coding transcripts)
RUN git clone https://github.com/TransDecoder/TransDecoder.git /opt/TransDecoder && \
    chmod +x /opt/TransDecoder/TransDecoder.LongOrfs /opt/TransDecoder/TransDecoder.Predict || true
ENV PATH="/opt/TransDecoder:$PATH"

# FEELnc (lncRNA classification: filter + codpot + classifier)
RUN apt-get update && apt-get install -y cpanminus libdb-dev 2>/dev/null && \
    cpanm Bio::Perl Bio::DB::Fasta Parallel::ForkManager 2>/dev/null || true && \
    git clone https://github.com/tderber/FEELnc.git /opt/FEELnc || true
ENV FEELNCPATH="/opt/FEELnc"
ENV PATH="/opt/FEELnc/scripts:$PATH"
ENV PERL5LIB="/opt/FEELnc/lib:$PERL5LIB"

# CPC2 (coding potential calculator)
RUN git clone https://github.com/gao-lab/CPC2_standalone.git /opt/cpc2 && \
    cd /opt/cpc2 && chmod +x bin/CPC2.py 2>/dev/null || true
ENV PATH="/opt/cpc2/bin:$PATH"

# CPAT (coding potential assessment tool)
RUN pip install --no-cache-dir CPAT || true

# ==============================================================
# v3.0 MODULE 2: Small RNA Analysis Tools
# ==============================================================

# cutadapt + Trim Galore (adapter trimming for Illumina sRNA-seq)
RUN pip install --no-cache-dir cutadapt && \
    wget -q -O /usr/local/bin/trim_galore \
        "https://raw.githubusercontent.com/FelixKrueger/TrimGalore/master/trim_galore" && \
    chmod +x /usr/local/bin/trim_galore || true

# ShortStack v4 (sRNA locus annotation + MIRNA identification)
RUN pip install --no-cache-dir ShortStack || \
    (git clone https://github.com/MikeAxtworthy/ShortStack.git /opt/ShortStack && \
     cd /opt/ShortStack && pip install --no-cache-dir .) || true

# miRDeep-P2 (plant miRNA prediction)
RUN git clone https://github.com/mhbrennan/miRDeep-P2.git /opt/mirdeep-p2 || true
ENV PATH="/opt/mirdeep-p2:$PATH"

# ==============================================================
# v3.0 MODULE 3: miRNA Target Prediction Tools
# ==============================================================

# TargetFinder (plant miRNA target prediction)
RUN git clone https://github.com/carringtonlab/TargetFinder.git /opt/targetfinder && \
    chmod +x /opt/targetfinder/targetfinder.pl 2>/dev/null || true
ENV PATH="/opt/targetfinder:$PATH"

# CleaveLand4 (degradome / PARE analysis)
RUN pip install --no-cache-dir cleaveland4 || \
    (git clone https://github.com/MikeAxtworthy/CleaveLand4.git /opt/CleaveLand4 && \
     cd /opt/CleaveLand4 && pip install --no-cache-dir .) || true

# ==============================================================
# v3.0 MODULE 4: Enhanced Epitranscriptomics Tools
# ==============================================================

# Nanocompore (differential RNA modification from nanopore)
RUN pip install --no-cache-dir nanocompore || true

# matplotlib + seaborn + upsetplot (visualization for reports)
RUN pip install --no-cache-dir matplotlib seaborn upsetplot jinja2 || true

# ==============================================================
# v3.0 MODULE 5: Data Integration Tools
# ==============================================================

# WGCNA (weighted gene co-expression network analysis, R package)
RUN Rscript -e 'if (!requireNamespace("WGCNA", quietly=TRUE)) BiocManager::install("WGCNA", ask=FALSE)' || true

# bcftools (variant analysis for population context)
RUN apt-get update && apt-get install -y bcftools 2>/dev/null || \
    (wget -q https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
     tar -xjf bcftools-1.21.tar.bz2 && cd bcftools-1.21 && \
     ./configure --prefix=/usr/local && make -j$(nproc) && make install && \
     cd .. && rm -rf bcftools-1.21*) || true

# ==============================================================
# 19b. Pin pandas/numpy for ELIGOS2 compatibility
# ==============================================================
# ELIGOS2 uses pd.concat with strings (fixed in pandas 1.x, breaks in 2.x)
# Must be installed AFTER all other Python packages to override their deps
RUN pip install --no-cache-dir "pandas==1.5.3" "numpy<2.0"

# Ensure setuptools < 71 for pycoQC compatibility (must be last pip install)
RUN pip install --no-cache-dir "setuptools<71"

# ==============================================================
# 20. Set up workspace
# ==============================================================
WORKDIR /workspace

# Copy project files (data/ is mounted at runtime via docker -v)
COPY scripts /workspace/scripts
COPY config /workspace/config
COPY requirements.txt /workspace/requirements.txt
COPY Snakefile /workspace/Snakefile

# Create data directory structure (will be overridden by volume mount)
RUN mkdir -p /workspace/data/raw/fastq /workspace/data/raw/fast5 \
    /workspace/data/raw/pod5 /workspace/data/raw/summaries \
    /workspace/data/reference/genome /workspace/data/reference/annotations \
    /workspace/data/reference/transcriptome

# ==============================================================
# 21. Fix setuptools for pycoQC (pkg_resources removed in v71+)
# ==============================================================
# Must be the LAST pip install to prevent later packages upgrading it
RUN pip install --no-cache-dir "setuptools<71"

# ==============================================================
# 22. Environment configuration
# ==============================================================
ENV PYTHONPATH="/usr/local/lib/python3.10/dist-packages"
ENV PATH="/home/eligos2:/home/eligos2/Scripts:$PATH"
ENV CUDA_HOME=/usr/local/cuda

# ==============================================================
# Verify installation
# ==============================================================
RUN echo "=== K-CHOPORE Build Verification ===" && \
    echo "Minimap2: $(minimap2 --version 2>&1 || echo 'not found')" && \
    echo "Samtools: $(samtools --version 2>&1 | head -1 || echo 'not found')" && \
    echo "Dorado: $(dorado --version 2>&1 || echo 'not found')" && \
    echo "StringTie: $(stringtie --version 2>&1 || echo 'not found')" && \
    echo "Nanopolish: $(nanopolish --version 2>&1 || echo 'not found')" && \
    echo "NanoPlot: $(NanoPlot --version 2>&1 || echo 'not found')" && \
    echo "NanoFilt: $(NanoFilt --version 2>&1 || echo 'not found')" && \
    echo "m6anet: $(pip show m6anet 2>&1 | grep Version || echo 'not found')" && \
    echo "MultiQC: $(multiqc --version 2>&1 || echo 'not found')" && \
    echo "Snakemake: $(snakemake --version 2>&1 || echo 'not found')" && \
    echo "R: $(R --version 2>&1 | head -1 || echo 'not found')" && \
    echo "=== K-CHOPORE Docker image built successfully! ==="
