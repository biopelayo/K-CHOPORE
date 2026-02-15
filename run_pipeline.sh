#!/bin/bash
# ==============================================================
# K-CHOPORE Pipeline Launcher
# ==============================================================
# This script runs the full K-CHOPORE pipeline inside Docker.
#
# USAGE:
#   ./run_pipeline.sh /path/to/your/data [threads]
#
# EXAMPLE:
#   ./run_pipeline.sh /mnt/e/kchopore_data 12
#   ./run_pipeline.sh D:/ONT_data 8
#
# Your data directory MUST contain:
#   data/
#   ├── raw/
#   │   ├── fastq/
#   │   │   ├── WT_C_R1.fastq
#   │   │   └── WT_C_R2.fastq
#   │   ├── fast5/           (needed for m6Anet/Nanopolish)
#   │   │   ├── WT_C_R1/
#   │   │   └── WT_C_R2/
#   │   └── summaries/
#   │       ├── WT_C_R1_sequencing_summary_FAR90122_d34138fc.txt
#   │       └── WT_C_R2_sequencing_summary_FAR91957_a56dafa5.txt
#   └── reference/
#       ├── genome/
#       │   └── TAIR10_chr_all.fas.fasta
#       └── annotations/
#           └── AtRTDv2_QUASI_19April2016.gtf
#
# The pipeline will create results/ and logs/ inside your data directory.
# ==============================================================

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
DATA_DIR="${1:-}"
THREADS="${2:-12}"
IMAGE_NAME="k-chopore:latest"

if [ -z "$DATA_DIR" ]; then
    echo -e "${RED}ERROR: Please specify the path to your data directory.${NC}"
    echo ""
    echo "Usage: $0 /path/to/your/data [threads]"
    echo "Example: $0 /mnt/e/kchopore_data 12"
    exit 1
fi

# Convert Windows paths to Unix if needed (for Git Bash / WSL)
DATA_DIR=$(echo "$DATA_DIR" | sed 's|\\|/|g')

# Verify data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo -e "${RED}ERROR: Data directory not found: $DATA_DIR${NC}"
    exit 1
fi

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN} K-CHOPORE Pipeline Launcher${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "Data directory: ${YELLOW}$DATA_DIR${NC}"
echo -e "Threads:        ${YELLOW}$THREADS${NC}"
echo -e "Docker image:   ${YELLOW}$IMAGE_NAME${NC}"
echo ""

# Verify Docker image exists
if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo -e "${RED}ERROR: Docker image '$IMAGE_NAME' not found.${NC}"
    echo "Build it first with: docker build -t k-chopore:latest ."
    exit 1
fi

# Check for required files
echo "Checking required files..."
MISSING=0

check_file() {
    if [ ! -f "$DATA_DIR/$1" ]; then
        echo -e "  ${RED}MISSING: $1${NC}"
        MISSING=1
    else
        SIZE=$(du -sh "$DATA_DIR/$1" 2>/dev/null | cut -f1)
        echo -e "  ${GREEN}OK: $1 ($SIZE)${NC}"
    fi
}

check_dir() {
    if [ ! -d "$DATA_DIR/$1" ]; then
        echo -e "  ${YELLOW}WARNING: Directory $1 not found (optional for some modules)${NC}"
    else
        COUNT=$(ls "$DATA_DIR/$1" 2>/dev/null | wc -l)
        echo -e "  ${GREEN}OK: $1/ ($COUNT items)${NC}"
    fi
}

# Reference files
check_file "data/reference/genome/TAIR10_chr_all.fas.fasta"
check_file "data/reference/annotations/AtRTDv2_QUASI_19April2016.gtf"

# Raw data
check_file "data/raw/fastq/WT_C_R1.fastq"
check_file "data/raw/fastq/WT_C_R2.fastq"

# Optional but recommended
check_file "data/raw/summaries/WT_C_R1_sequencing_summary_FAR90122_d34138fc.txt"
check_file "data/raw/summaries/WT_C_R2_sequencing_summary_FAR91957_a56dafa5.txt"
check_dir "data/raw/fast5"

echo ""
if [ "$MISSING" -eq 1 ]; then
    echo -e "${RED}Some required files are missing! Please check your data directory.${NC}"
    echo "See the header of this script for the expected directory structure."
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Create output directories on the host
mkdir -p "$DATA_DIR/results" "$DATA_DIR/logs"

# Run the pipeline
echo -e "${GREEN}Starting K-CHOPORE pipeline...${NC}"
echo ""

LOGFILE="$DATA_DIR/logs/pipeline_run_$(date +%Y%m%d_%H%M%S).log"

docker run --rm \
    -v "$DATA_DIR":/workspace \
    -w /workspace \
    "$IMAGE_NAME" \
    snakemake \
        --snakefile /workspace/Snakefile \
        --configfile /workspace/config/config.yml \
        --cores "$THREADS" \
        --latency-wait 60 \
        --printshellcmds \
        --reason \
        --keep-going \
        --rerun-incomplete \
        > >(tee "$LOGFILE") 2>&1

EXIT_CODE=$?

echo ""
if [ "$EXIT_CODE" -eq 0 ]; then
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN} Pipeline completed successfully!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo -e "Results: ${YELLOW}$DATA_DIR/results/${NC}"
    echo -e "Logs:    ${YELLOW}$DATA_DIR/logs/${NC}"
else
    echo -e "${RED}========================================${NC}"
    echo -e "${RED} Pipeline finished with errors (exit: $EXIT_CODE)${NC}"
    echo -e "${RED}========================================${NC}"
    echo -e "Check logs: ${YELLOW}$DATA_DIR/logs/${NC}"
fi

exit $EXIT_CODE
