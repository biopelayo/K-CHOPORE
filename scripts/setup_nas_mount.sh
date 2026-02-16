#!/bin/bash
# =============================================================
# K-CHOPORE NAS Mount Setup Script
# =============================================================
# Mounts the Synology NAS via SSHFS so Docker can read FAST5
# data directly from the NAS without copying 900 GB locally.
#
# Run this ONCE on the server before launching the pipeline.
# Requires sudo for installing sshfs and creating mount point.
#
# Usage:
#   ./scripts/setup_nas_mount.sh <nas_user> <nas_host> <nas_path> [mount_point]
#
# Example:
#   ./scripts/setup_nas_mount.sh usuario valmeilab.synology.me \
#       /HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS \
#       /mnt/nas
#
# After mounting, verify with:
#   ls /mnt/nas/WT_C_R1/*.fast5 | head
# =============================================================

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# Parse arguments
NAS_USER="${1:-}"
NAS_HOST="${2:-}"
NAS_PATH="${3:-}"
MOUNT_POINT="${4:-/mnt/nas}"

if [ -z "$NAS_USER" ] || [ -z "$NAS_HOST" ] || [ -z "$NAS_PATH" ]; then
    echo -e "${RED}ERROR: Missing arguments.${NC}"
    echo ""
    echo "Usage: $0 <nas_user> <nas_host> <nas_path> [mount_point]"
    echo ""
    echo "Arguments:"
    echo "  nas_user:     SSH username for the NAS (e.g., 'usuario')"
    echo "  nas_host:     NAS hostname or IP (e.g., 'valmeilab.synology.me')"
    echo "  nas_path:     Path to FAST5 data on the NAS"
    echo "  mount_point:  Local mount point (default: /mnt/nas)"
    echo ""
    echo "Example:"
    echo "  $0 usuario valmeilab.synology.me \\"
    echo "      /HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS"
    exit 1
fi

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN} K-CHOPORE NAS Mount Setup${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo -e "NAS:         ${CYAN}${NAS_USER}@${NAS_HOST}:${NAS_PATH}${NC}"
echo -e "Mount point: ${CYAN}${MOUNT_POINT}${NC}"
echo ""

# -------------------------------------------------------------
# Step 1: Install sshfs if not present
# -------------------------------------------------------------
echo -e "${YELLOW}[Step 1/5] Checking sshfs...${NC}"

if command -v sshfs &> /dev/null; then
    echo -e "  ${GREEN}sshfs already installed: $(which sshfs)${NC}"
else
    echo -e "  ${YELLOW}Installing sshfs...${NC}"
    if command -v apt-get &> /dev/null; then
        sudo apt-get update -qq && sudo apt-get install -y -qq sshfs
    elif command -v yum &> /dev/null; then
        sudo yum install -y sshfs
    elif command -v dnf &> /dev/null; then
        sudo dnf install -y sshfs
    else
        echo -e "  ${RED}ERROR: Could not find package manager to install sshfs.${NC}"
        echo "  Please install sshfs manually and re-run this script."
        exit 1
    fi
    echo -e "  ${GREEN}sshfs installed successfully.${NC}"
fi

# -------------------------------------------------------------
# Step 2: Check FUSE is available
# -------------------------------------------------------------
echo -e "${YELLOW}[Step 2/5] Checking FUSE support...${NC}"

if [ -e /dev/fuse ]; then
    echo -e "  ${GREEN}/dev/fuse exists.${NC}"
else
    echo -e "  ${YELLOW}Loading fuse module...${NC}"
    sudo modprobe fuse 2>/dev/null || true
    if [ -e /dev/fuse ]; then
        echo -e "  ${GREEN}/dev/fuse now available.${NC}"
    else
        echo -e "  ${RED}WARNING: /dev/fuse not found. SSHFS may not work.${NC}"
    fi
fi

# Enable allow_other in fuse config (needed for Docker to read the mount)
FUSE_CONF="/etc/fuse.conf"
if [ -f "$FUSE_CONF" ]; then
    if grep -q "^#user_allow_other" "$FUSE_CONF"; then
        echo -e "  ${YELLOW}Enabling user_allow_other in $FUSE_CONF...${NC}"
        sudo sed -i 's/^#user_allow_other/user_allow_other/' "$FUSE_CONF"
        echo -e "  ${GREEN}Enabled.${NC}"
    elif grep -q "^user_allow_other" "$FUSE_CONF"; then
        echo -e "  ${GREEN}user_allow_other already enabled.${NC}"
    else
        echo -e "  ${YELLOW}Adding user_allow_other to $FUSE_CONF...${NC}"
        echo "user_allow_other" | sudo tee -a "$FUSE_CONF" > /dev/null
        echo -e "  ${GREEN}Added.${NC}"
    fi
fi

# -------------------------------------------------------------
# Step 3: Create mount point
# -------------------------------------------------------------
echo -e "${YELLOW}[Step 3/5] Creating mount point...${NC}"

if mountpoint -q "$MOUNT_POINT" 2>/dev/null; then
    echo -e "  ${GREEN}$MOUNT_POINT is already mounted.${NC}"
    echo -e "  ${YELLOW}Unmounting first...${NC}"
    sudo umount "$MOUNT_POINT" || fusermount -u "$MOUNT_POINT"
fi

sudo mkdir -p "$MOUNT_POINT"
echo -e "  ${GREEN}Mount point ready: $MOUNT_POINT${NC}"

# -------------------------------------------------------------
# Step 4: Mount NAS via SSHFS
# -------------------------------------------------------------
echo -e "${YELLOW}[Step 4/5] Mounting NAS via SSHFS...${NC}"
echo -e "  ${CYAN}You will be prompted for the NAS SSH password.${NC}"
echo ""

sshfs "${NAS_USER}@${NAS_HOST}:${NAS_PATH}" "$MOUNT_POINT" \
    -o allow_other \
    -o ro \
    -o reconnect \
    -o ServerAliveInterval=15 \
    -o ServerAliveCountMax=3 \
    -o cache=yes \
    -o kernel_cache \
    -o large_read \
    -o max_read=65536

echo ""
echo -e "  ${GREEN}NAS mounted successfully at $MOUNT_POINT${NC}"

# -------------------------------------------------------------
# Step 5: Verify FAST5 data is accessible
# -------------------------------------------------------------
echo -e "${YELLOW}[Step 5/5] Verifying FAST5 data access...${NC}"
echo ""

# Expected sample folders (from config)
EXPECTED_DIRS=(
    "WT_C_R1" "WT_C_R2" "WT_C_R3"
    "WT_AA_R1" "WT_AA_R2" "WT_AA_R3" "WT_AA_R3_2"
    "anac017-1_C_R1" "anac017-1_C_R2" "anac017-1_C_R3"
    "anac017-1_AA_R1" "anac017-1_AA_R2" "anac017-1_AA_R3"
    "anac017_AA_R2-2"
)

FOUND=0
MISSING=0

for dir in "${EXPECTED_DIRS[@]}"; do
    if [ -d "$MOUNT_POINT/$dir" ]; then
        # Count FAST5 files
        N_FAST5=$(find "$MOUNT_POINT/$dir" -name "*.fast5" -maxdepth 2 2>/dev/null | head -100 | wc -l)
        SIZE=$(du -sh "$MOUNT_POINT/$dir" 2>/dev/null | cut -f1 || echo "?")
        echo -e "  ${GREEN}OK${NC}: $dir/ (${N_FAST5}+ fast5 files, ~${SIZE})"
        FOUND=$((FOUND + 1))
    else
        echo -e "  ${RED}MISSING${NC}: $dir/"
        MISSING=$((MISSING + 1))
    fi
done

echo ""
echo -e "  Found: ${GREEN}$FOUND${NC} / ${#EXPECTED_DIRS[@]} expected directories"

if [ "$MISSING" -gt 0 ]; then
    echo -e "  ${RED}WARNING: $MISSING directories missing.${NC}"
    echo "  Check that NAS_PATH is correct and contains the FAST5 folders."
fi

# Disk space report
echo ""
echo -e "${YELLOW}Disk space on server:${NC}"
df -h / | head -2
echo ""

# Summary
echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN} NAS Mount Setup Complete${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo "Next steps:"
echo "  1. Launch the pipeline with:"
echo ""
echo -e "     ${CYAN}./run_pipeline.sh /path/to/kchopore_data $MOUNT_POINT 40${NC}"
echo ""
echo "  2. To unmount the NAS later:"
echo ""
echo -e "     ${CYAN}sudo umount $MOUNT_POINT${NC}"
echo "     or: fusermount -u $MOUNT_POINT"
echo ""
echo "  3. To remount after a reboot, re-run this script."
echo ""
