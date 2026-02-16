#!/bin/bash
# ==============================================================
# K-CHOPORE Pipeline Launcher
# ==============================================================
# Runs the full K-CHOPORE pipeline inside Docker.
# Data is downloaded from NAS to local data/raw/ before running.
#
# USAGE:
#   ./run_pipeline.sh [threads]
#
# ARGUMENTS:
#   threads:  Number of CPU cores (default: 40)
#
# EXAMPLE:
#   ./run_pipeline.sh 40
#
# EXPECTED DIRECTORY STRUCTURE:
#   /home/usuario2/pelamovic/kchopore/
#   ├── Snakefile
#   ├── config/config.yml
#   ├── scripts/run_deseq2.R
#   ├── data/
#   │   ├── raw/            ← Downloaded from NAS
#   │   │   ├── WT_C_R1/    (fastq_pass/ + sequencing_summary)
#   │   │   ├── WT_C_R2/    ...
#   │   │   └── ...
#   │   └── reference/
#   │       ├── genome/TAIR10_chr_all.fas.fasta
#   │       └── annotations/AtRTDv2_QUASI_19April2016.gtf
#   └── results/             ← Pipeline outputs
# ==============================================================

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Parse arguments
THREADS="${1:-40}"
IMAGE_NAME="k-chopore:latest"
DATA_DIR="/home/usuario2/pelamovic/kchopore"

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN} K-CHOPORE Pipeline Launcher${NC}"
echo -e "${GREEN} 12 Samples | 2x2 Factorial Design${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo -e "Data directory: ${CYAN}$DATA_DIR${NC}"
echo -e "Threads:        ${CYAN}$THREADS${NC}"
echo -e "Docker image:   ${CYAN}$IMAGE_NAME${NC}"
echo ""

# -------------------------------------------------------------
# Preflight checks
# -------------------------------------------------------------

# Check Docker image
if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo -e "${RED}ERROR: Docker image '$IMAGE_NAME' not found.${NC}"
    echo "Build it first with: docker build -t k-chopore:latest ."
    exit 1
fi

# Check required files
echo "Checking required files..."
MISSING=0

check_file() {
    if [ ! -f "$DATA_DIR/$1" ]; then
        echo -e "  ${RED}MISSING: $1${NC}"
        MISSING=1
    else
        SIZE=$(du -sh "$DATA_DIR/$1" 2>/dev/null | cut -f1)
        echo -e "  ${GREEN}OK${NC}: $1 ($SIZE)"
    fi
}

# Reference files
check_file "data/reference/genome/TAIR10_chr_all.fas.fasta"
check_file "data/reference/annotations/AtRTDv2_QUASI_19April2016.gtf"

# Pipeline files
check_file "Snakefile"
check_file "config/config.yml"
check_file "scripts/run_deseq2.R"

# Check raw data directories
echo ""
echo "Checking raw data (downloaded from NAS)..."
RAW_DIRS=(
    "WT_C_R1" "WT_C_R2" "WT_C_R3"
    "WT_AA_R1" "WT_AA_R2" "WT_AA_R3"
    "anac017-1_C_R1" "anac017-1_C_R2" "anac017-1_C_R3"
    "anac017-1_AA_R1"
)
FAST5_DIRS=("anac017-1_AA_R2" "anac017-1_AA_R3" "anac017_AA_R2-2")

for dir in "${RAW_DIRS[@]}"; do
    if [ -d "$DATA_DIR/data/raw/$dir" ]; then
        SIZE=$(du -sh "$DATA_DIR/data/raw/$dir" 2>/dev/null | cut -f1)
        echo -e "  ${GREEN}OK${NC}: data/raw/$dir/ ($SIZE)"
    else
        echo -e "  ${RED}MISSING: data/raw/$dir/${NC}"
        MISSING=1
    fi
done

for dir in "${FAST5_DIRS[@]}"; do
    if [ -d "$DATA_DIR/data/raw/$dir" ]; then
        SIZE=$(du -sh "$DATA_DIR/data/raw/$dir" 2>/dev/null | cut -f1)
        echo -e "  ${GREEN}OK${NC}: data/raw/$dir/ ($SIZE) [FAST5-only]"
    else
        echo -e "  ${YELLOW}MISSING: data/raw/$dir/ (FAST5-only, needs basecalling)${NC}"
    fi
done

echo ""
if [ "$MISSING" -eq 1 ]; then
    echo -e "${YELLOW}Some data directories are missing.${NC}"
    echo "Transfer data with: python scripts/transfer_fastq_to_server.py"
    echo ""
fi

# Report disk space
echo -e "${YELLOW}Disk space:${NC}"
echo -n "  Server: "
df -h "$DATA_DIR" 2>/dev/null | tail -1 | awk '{print $4 " available of " $2}'
echo ""

# Create output directories
mkdir -p "$DATA_DIR/results" "$DATA_DIR/logs"

# Pre-create .snakemake directory (avoids Docker/NTFS race condition)
mkdir -p "$DATA_DIR/.snakemake/locks" \
         "$DATA_DIR/.snakemake/metadata" \
         "$DATA_DIR/.snakemake/incomplete" \
         "$DATA_DIR/.snakemake/shadow"

# -------------------------------------------------------------
# Launch pipeline in Docker
# -------------------------------------------------------------
echo -e "${GREEN}Starting K-CHOPORE pipeline (12 samples, $THREADS cores)...${NC}"
echo ""

LOGFILE="$DATA_DIR/logs/pipeline_run_$(date +%Y%m%d_%H%M%S).log"
echo -e "Log file: ${CYAN}$LOGFILE${NC}"
echo ""

docker run --rm \
    --name k-chopore-run \
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
        > "$LOGFILE" 2>&1

EXIT_CODE=$?

echo ""
if [ "$EXIT_CODE" -eq 0 ]; then
    echo -e "${GREEN}============================================${NC}"
    echo -e "${GREEN} Pipeline completed successfully!${NC}"
    echo -e "${GREEN}============================================${NC}"
    echo ""
    echo -e "Results: ${CYAN}$DATA_DIR/results/${NC}"
    echo -e "Log:     ${CYAN}$LOGFILE${NC}"
    echo ""
    echo "Key outputs:"
    echo "  results/deseq2/  - Differential expression (3 contrasts)"
    echo "  results/flair/   - Isoform analysis"
    echo "  results/eligos/  - RNA modifications (ELIGOS2)"
    echo "  results/multiqc/ - Aggregated QC report"
    echo ""
    echo "For m6Anet (requires FAST5 data on disk, ~55 GB/sample):"
    echo "  Transfer FAST5 for one sample, then run:"
    echo "  snakemake results/m6anet/SAMPLE/data.site_proba.csv --cores $THREADS"
else
    echo -e "${RED}============================================${NC}"
    echo -e "${RED} Pipeline finished with errors (exit: $EXIT_CODE)${NC}"
    echo -e "${RED}============================================${NC}"
    echo ""
    echo -e "Check log: ${CYAN}$LOGFILE${NC}"
    echo ""
    echo "To see the last 50 lines:"
    echo "  tail -50 $LOGFILE"
    echo ""
    echo "To resume from where it stopped:"
    echo "  $0 $THREADS"
fi

exit $EXIT_CODE
