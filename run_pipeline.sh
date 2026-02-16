#!/bin/bash
# ==============================================================
# K-CHOPORE Pipeline Launcher
# ==============================================================
# Runs the full K-CHOPORE pipeline inside Docker with NAS mount
# for reading FAST5 data directly from the Synology NAS.
#
# USAGE:
#   ./run_pipeline.sh <data_dir> <nas_mount> [threads]
#
# ARGUMENTS:
#   data_dir:   Path to K-CHOPORE project data directory
#   nas_mount:  Path where NAS is mounted via SSHFS (e.g., /mnt/nas)
#   threads:    Number of CPU cores (default: 40)
#
# EXAMPLE:
#   # First mount the NAS (one-time setup):
#   ./scripts/setup_nas_mount.sh usuario valmeilab.synology.me \
#       /HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS
#
#   # Then launch the pipeline:
#   ./run_pipeline.sh /home/usuario2/K-CHOPORE /mnt/nas 40
#
# DATA DIRECTORY MUST CONTAIN:
#   data/
#   └── reference/
#       ├── genome/
#       │   └── TAIR10_chr_all.fas.fasta
#       └── annotations/
#           └── AtRTDv2_QUASI_19April2016.gtf
#
# NAS MOUNT MUST CONTAIN (read-only via SSHFS):
#   WT_C_R1/        WT_C_R2/        WT_C_R3/
#   WT_AA_R1/       WT_AA_R2/       WT_AA_R3/       WT_AA_R3_2/
#   anac017-1_C_R1/ anac017-1_C_R2/ anac017-1_C_R3/
#   anac017-1_AA_R1/ anac017-1_AA_R2/ anac017-1_AA_R3/ anac017_AA_R2-2/
#
# OUTPUT (created in data_dir):
#   results/basecalls/    - Guppy basecalled FASTQs
#   results/filtered/     - NanoFilt quality-filtered reads
#   results/qc/           - NanoPlot, pycoQC reports
#   results/mapped/       - Minimap2 BAM alignments
#   results/flair/        - Isoform analysis
#   results/eligos2/      - RNA modification calls
#   results/m6anet/       - m6A modification predictions
#   results/deseq2/       - Factorial differential expression
#   results/multiqc/      - Aggregated QC report
# ==============================================================

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Parse arguments
DATA_DIR="${1:-}"
NAS_MOUNT="${2:-}"
THREADS="${3:-40}"
IMAGE_NAME="k-chopore:latest"

if [ -z "$DATA_DIR" ] || [ -z "$NAS_MOUNT" ]; then
    echo -e "${RED}ERROR: Missing required arguments.${NC}"
    echo ""
    echo "Usage: $0 <data_dir> <nas_mount> [threads]"
    echo ""
    echo "Arguments:"
    echo "  data_dir:   Path to K-CHOPORE project data"
    echo "  nas_mount:  Path where NAS is SSHFS-mounted (e.g., /mnt/nas)"
    echo "  threads:    CPU cores (default: 40)"
    echo ""
    echo "Example:"
    echo "  $0 /home/usuario2/K-CHOPORE /mnt/nas 40"
    exit 1
fi

# Convert Windows paths if needed (Git Bash / WSL)
DATA_DIR=$(echo "$DATA_DIR" | sed 's|\\|/|g')
NAS_MOUNT=$(echo "$NAS_MOUNT" | sed 's|\\|/|g')

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN} K-CHOPORE Pipeline Launcher${NC}"
echo -e "${GREEN} 12 Samples | 2x2 Factorial Design${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo -e "Data directory: ${CYAN}$DATA_DIR${NC}"
echo -e "NAS mount:      ${CYAN}$NAS_MOUNT${NC}"
echo -e "Threads:        ${CYAN}$THREADS${NC}"
echo -e "Docker image:   ${CYAN}$IMAGE_NAME${NC}"
echo ""

# -------------------------------------------------------------
# Preflight checks
# -------------------------------------------------------------

# Check data directory
if [ ! -d "$DATA_DIR" ]; then
    echo -e "${RED}ERROR: Data directory not found: $DATA_DIR${NC}"
    exit 1
fi

# Check NAS mount
if [ ! -d "$NAS_MOUNT" ]; then
    echo -e "${RED}ERROR: NAS mount point not found: $NAS_MOUNT${NC}"
    echo "Run scripts/setup_nas_mount.sh first to mount the NAS."
    exit 1
fi

if ! mountpoint -q "$NAS_MOUNT" 2>/dev/null; then
    echo -e "${YELLOW}WARNING: $NAS_MOUNT does not appear to be a mount point.${NC}"
    echo "If the NAS is mounted at a parent directory, this may still work."
    echo ""
fi

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

check_nas_dir() {
    if [ ! -d "$NAS_MOUNT/$1" ]; then
        echo -e "  ${RED}MISSING${NC}: NAS:$1/"
        MISSING=1
    else
        echo -e "  ${GREEN}OK${NC}: NAS:$1/"
    fi
}

# Reference files (on local disk)
check_file "data/reference/genome/TAIR10_chr_all.fas.fasta"
check_file "data/reference/annotations/AtRTDv2_QUASI_19April2016.gtf"

# Pipeline files
check_file "Snakefile"
check_file "config/config.yml"
check_file "scripts/run_deseq2.R"

# FAST5 folders on NAS
echo ""
echo "Checking FAST5 folders on NAS..."
NAS_DIRS=(
    "WT_C_R1" "WT_C_R2" "WT_C_R3"
    "WT_AA_R1" "WT_AA_R2" "WT_AA_R3" "WT_AA_R3_2"
    "anac017-1_C_R1" "anac017-1_C_R2" "anac017-1_C_R3"
    "anac017-1_AA_R1" "anac017-1_AA_R2" "anac017-1_AA_R3"
    "anac017_AA_R2-2"
)

for dir in "${NAS_DIRS[@]}"; do
    check_nas_dir "$dir"
done

echo ""
if [ "$MISSING" -eq 1 ]; then
    echo -e "${RED}Some required files or directories are missing!${NC}"
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
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
    -v "$NAS_MOUNT":/nas:ro \
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
    echo "  results/eligos2/ - RNA modifications"
    echo "  results/m6anet/  - m6A predictions"
    echo "  results/multiqc/ - Aggregated QC report"
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
    echo "  $0 $DATA_DIR $NAS_MOUNT $THREADS"
fi

exit $EXIT_CODE
