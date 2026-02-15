#!/bin/bash
# ==============================================================
# K-CHOPORE Server Deployment Script
# ==============================================================
# Run this on your Linux server to set up and launch the pipeline.
#
# PREREQUISITES:
#   - Docker installed and running (docker --version)
#   - Git installed (git --version)
#   - ~300 GB free disk space (data + intermediate files)
#
# USAGE:
#   1. Transfer your data to the server first (see DATA TRANSFER below)
#   2. Run: bash deploy_server.sh /path/to/workspace [threads]
#
# DATA TRANSFER (from Windows to server):
#   Use rsync or scp to transfer data. Example with rsync:
#
#   # From Windows (Git Bash), transfer data to server:
#   rsync -avP --progress E:/kchopore_workspace/data/ user@server:/data/kchopore/data/
#
#   # Or use scp (slower, no resume):
#   scp -r E:/kchopore_workspace/data/ user@server:/data/kchopore/data/
#
# The workspace directory should contain:
#   data/
#   +-- raw/
#   |   +-- fastq/
#   |   |   +-- WT_C_R1.fastq        (2.9 GB)
#   |   |   +-- WT_C_R2.fastq        (3.4 GB)
#   |   +-- fast5/
#   |   |   +-- WT_C_R2/             (73 GB, 384 files)
#   |   +-- summaries/
#   |       +-- WT_C_R1_sequencing_summary_FAR90122_d34138fc.txt  (455 MB)
#   |       +-- WT_C_R2_sequencing_summary_FAR91957_a56dafa5.txt  (572 MB)
#   +-- reference/
#       +-- genome/
#       |   +-- TAIR10_chr_all.fas.fasta  (116 MB)
#       +-- annotations/
#           +-- AtRTDv2_QUASI_19April2016.gtf  (55 MB)
#
# Total data to transfer: ~81 GB
# Total disk needed (data + results): ~300 GB
# ==============================================================

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Parse arguments
WORKSPACE="${1:-}"
THREADS="${2:-$(nproc)}"
IMAGE_NAME="k-chopore:latest"

if [ -z "$WORKSPACE" ]; then
    echo -e "${RED}ERROR: Please specify the workspace path.${NC}"
    echo ""
    echo "Usage: bash $0 /path/to/workspace [threads]"
    echo ""
    echo "The workspace should contain your data/ directory."
    echo "If starting fresh, data will be expected at /path/to/workspace/data/"
    exit 1
fi

echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN} K-CHOPORE Server Deployment${NC}"
echo -e "${CYAN}========================================${NC}"
echo ""
echo -e "Workspace:  ${YELLOW}$WORKSPACE${NC}"
echo -e "Threads:    ${YELLOW}$THREADS${NC}"
echo -e "Docker img: ${YELLOW}$IMAGE_NAME${NC}"
echo ""

# --- Step 1: System checks ---
echo -e "${GREEN}[1/6] Checking system requirements...${NC}"

# Docker
if ! command -v docker &> /dev/null; then
    echo -e "${RED}ERROR: Docker not found. Install Docker first:${NC}"
    echo "  curl -fsSL https://get.docker.com | sh"
    echo "  sudo usermod -aG docker \$USER"
    exit 1
fi
echo -e "  Docker: ${GREEN}$(docker --version)${NC}"

# Check Docker daemon
if ! docker info &> /dev/null; then
    echo -e "${RED}ERROR: Docker daemon not running. Start with: sudo systemctl start docker${NC}"
    exit 1
fi

# Git
if ! command -v git &> /dev/null; then
    echo -e "${RED}ERROR: Git not found. Install with: sudo apt install git${NC}"
    exit 1
fi
echo -e "  Git: ${GREEN}$(git --version)${NC}"

# Disk space
AVAIL_GB=$(df -BG "$WORKSPACE" 2>/dev/null | tail -1 | awk '{print $4}' | tr -d 'G' || echo "unknown")
echo -e "  Disk free: ${YELLOW}${AVAIL_GB} GB${NC}"
if [ "$AVAIL_GB" != "unknown" ] && [ "$AVAIL_GB" -lt 200 ]; then
    echo -e "  ${YELLOW}WARNING: Less than 200 GB free. Pipeline may need up to 300 GB.${NC}"
fi

# CPU
echo -e "  CPUs: ${YELLOW}$(nproc)${NC}"
echo -e "  RAM: ${YELLOW}$(free -h | awk '/Mem:/{print $2}')${NC}"
echo ""

# --- Step 2: Clone repo ---
echo -e "${GREEN}[2/6] Setting up workspace...${NC}"
mkdir -p "$WORKSPACE"

if [ ! -f "$WORKSPACE/Snakefile" ]; then
    echo "  Cloning K-CHOPORE repository..."
    git clone https://github.com/biopelayo/K-CHOPORE.git "$WORKSPACE/repo"
    # Copy pipeline files into workspace (keep data separate)
    cp "$WORKSPACE/repo/Snakefile" "$WORKSPACE/"
    cp "$WORKSPACE/repo/Dockerfile" "$WORKSPACE/"
    cp -r "$WORKSPACE/repo/config" "$WORKSPACE/"
    cp -r "$WORKSPACE/repo/scripts" "$WORKSPACE/"
    echo -e "  ${GREEN}Repository cloned and files copied.${NC}"
else
    echo -e "  ${GREEN}Snakefile already present, skipping clone.${NC}"
fi
echo ""

# --- Step 3: Check data ---
echo -e "${GREEN}[3/6] Checking data files...${NC}"
MISSING=0

check_file() {
    if [ ! -f "$WORKSPACE/$1" ]; then
        echo -e "  ${RED}MISSING: $1${NC}"
        MISSING=1
    else
        SIZE=$(du -sh "$WORKSPACE/$1" 2>/dev/null | cut -f1)
        echo -e "  ${GREEN}OK: $1 ($SIZE)${NC}"
    fi
}

check_dir() {
    if [ ! -d "$WORKSPACE/$1" ]; then
        echo -e "  ${RED}MISSING: $1/${NC}"
        MISSING=1
    else
        COUNT=$(find "$WORKSPACE/$1" -type f | wc -l)
        SIZE=$(du -sh "$WORKSPACE/$1" 2>/dev/null | cut -f1)
        echo -e "  ${GREEN}OK: $1/ ($COUNT files, $SIZE)${NC}"
    fi
}

check_file "data/reference/genome/TAIR10_chr_all.fas.fasta"
check_file "data/reference/annotations/AtRTDv2_QUASI_19April2016.gtf"
check_file "data/raw/fastq/WT_C_R1.fastq"
check_file "data/raw/fastq/WT_C_R2.fastq"
check_file "data/raw/summaries/WT_C_R1_sequencing_summary_FAR90122_d34138fc.txt"
check_file "data/raw/summaries/WT_C_R2_sequencing_summary_FAR91957_a56dafa5.txt"
check_dir  "data/raw/fast5/WT_C_R2"

echo ""
if [ "$MISSING" -eq 1 ]; then
    echo -e "${RED}Data files are missing! Transfer them first:${NC}"
    echo ""
    echo "  From your Windows machine (Git Bash):"
    echo "  rsync -avP E:/kchopore_workspace/data/ user@this-server:$WORKSPACE/data/"
    echo ""
    echo "  Or transfer individual components:"
    echo "  rsync -avP E:/kchopore_workspace/data/raw/fastq/ user@this-server:$WORKSPACE/data/raw/fastq/"
    echo "  rsync -avP E:/kchopore_workspace/data/raw/fast5/ user@this-server:$WORKSPACE/data/raw/fast5/"
    echo "  rsync -avP E:/kchopore_workspace/data/raw/summaries/ user@this-server:$WORKSPACE/data/raw/summaries/"
    echo "  rsync -avP E:/kchopore_workspace/data/reference/ user@this-server:$WORKSPACE/data/reference/"
    echo ""
    read -p "Continue anyway (build Docker image while waiting for data)? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# --- Step 4: Build Docker image ---
echo -e "${GREEN}[4/6] Building Docker image (this takes 15-30 minutes)...${NC}"

if docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo -e "  ${GREEN}Docker image '$IMAGE_NAME' already exists.${NC}"
    read -p "  Rebuild? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "  Building..."
        docker build -t "$IMAGE_NAME" "$WORKSPACE/" 2>&1 | tail -20
    fi
else
    echo "  Building..."
    docker build -t "$IMAGE_NAME" "$WORKSPACE/" 2>&1 | tail -20
fi

# Verify image
if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo -e "${RED}ERROR: Docker image build failed!${NC}"
    exit 1
fi
echo -e "  ${GREEN}Docker image ready: $(docker image inspect "$IMAGE_NAME" --format '{{.Size}}' | numfmt --to=iec 2>/dev/null || docker images "$IMAGE_NAME" --format '{{.Size}}')${NC}"
echo ""

# --- Step 5: Verify tools inside container ---
echo -e "${GREEN}[5/6] Verifying tools inside container...${NC}"
docker run --rm "$IMAGE_NAME" bash -c '
    echo "  minimap2: $(minimap2 --version 2>/dev/null || echo MISSING)"
    echo "  samtools: $(samtools --version 2>/dev/null | head -1 || echo MISSING)"
    echo "  nanopolish: $(nanopolish --version 2>/dev/null | head -1 || echo MISSING)"
    echo "  flair: $(flair --version 2>/dev/null || echo present)"
    echo "  NanoPlot: $(NanoPlot --version 2>/dev/null || echo MISSING)"
    echo "  snakemake: $(snakemake --version 2>/dev/null || echo MISSING)"
    echo "  HDF5_PLUGIN_PATH: $HDF5_PLUGIN_PATH"
    python3 -c "from pkg_resources import resource_filename; print(\"  pkg_resources: OK\")" 2>/dev/null || echo "  pkg_resources: MISSING"
'
echo ""

# --- Step 6: Pre-create directories and launch ---
echo -e "${GREEN}[6/6] Creating output directories...${NC}"
mkdir -p "$WORKSPACE/results" "$WORKSPACE/logs" "$WORKSPACE/.snakemake"
echo -e "  ${GREEN}Done.${NC}"
echo ""

# --- Ready to run ---
echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN} Deployment complete!${NC}"
echo -e "${CYAN}========================================${NC}"
echo ""
echo "To run the full pipeline:"
echo ""
echo -e "  ${YELLOW}# Interactive (see output live):${NC}"
echo "  docker run --rm \\"
echo "    -v $WORKSPACE:/workspace \\"
echo "    -w /workspace \\"
echo "    $IMAGE_NAME \\"
echo "    snakemake --snakefile Snakefile --configfile config/config.yml \\"
echo "    --cores $THREADS --latency-wait 60 --keep-going --rerun-incomplete"
echo ""
echo -e "  ${YELLOW}# Background with nohup (recommended for long runs):${NC}"
echo "  nohup docker run --rm \\"
echo "    -v $WORKSPACE:/workspace \\"
echo "    -w /workspace \\"
echo "    $IMAGE_NAME \\"
echo "    snakemake --snakefile Snakefile --configfile config/config.yml \\"
echo "    --cores $THREADS --latency-wait 60 --keep-going --rerun-incomplete \\"
echo "    > $WORKSPACE/logs/pipeline_\$(date +%Y%m%d_%H%M%S).log 2>&1 &"
echo ""
echo -e "  ${YELLOW}# Using the launcher script:${NC}"
echo "  bash $WORKSPACE/run_pipeline.sh $WORKSPACE $THREADS"
echo ""
echo -e "  ${YELLOW}# Dry run (check what will be executed):${NC}"
echo "  docker run --rm -v $WORKSPACE:/workspace -w /workspace $IMAGE_NAME \\"
echo "    snakemake -n --snakefile Snakefile --configfile config/config.yml --cores $THREADS"
echo ""
echo -e "  ${YELLOW}# Monitor running pipeline:${NC}"
echo "  tail -f $WORKSPACE/logs/pipeline_*.log"
echo "  docker ps    # check running containers"
echo ""
