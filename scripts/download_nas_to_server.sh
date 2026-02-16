#!/bin/bash
# ==============================================================
# K-CHOPORE NAS Data Transfer Script
# ==============================================================
# Downloads FASTQ and sequencing summary data from the Synology
# NAS via the File Station API and uploads to the server via SSH.
#
# ARCHITECTURE: Windows -> NAS API (download) -> SSH -> Server
# This bypasses the NAS IP auto-block on the server by routing
# downloads through the Windows machine.
#
# USAGE (run from Windows/Git Bash):
#   ./scripts/download_nas_to_server.sh
#
# PREREQUISITES:
#   - curl, ssh available (Git Bash on Windows)
#   - SSH access to server (usuario2@156.35.42.17)
#   - Synology NAS API access (pelayo@valmeilab.synology.me:5001)
#
# DATA TRANSFERRED:
#   For 11 pre-basecalled samples:
#     - fastq_pass/*.fastq.gz files (~1.5 GB/sample compressed)
#     - sequencing_summary_*.txt (~500 MB/sample)
#   For 2 FAST5-only samples:
#     - fast5/*.fast5 files
#     - sequencing_summary_*.txt
#   Total: ~25 GB compressed FASTQs + ~6 GB summaries
# ==============================================================

set -euo pipefail

# Configuration
NAS_HOST="valmeilab.synology.me"
NAS_PORT="5001"
NAS_USER="pelayo"
NAS_PASS="Pelamovic39@"
NAS_BASE="/HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS"

SERVER_HOST="156.35.42.17"
SERVER_USER="usuario2"
SERVER_BASE="/home/usuario2/pelamovic/kchopore/data/raw"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# ---- Helper functions ----

nas_login() {
    local sid
    sid=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/auth.cgi" \
        --data-urlencode "api=SYNO.API.Auth" \
        --data-urlencode "version=6" \
        --data-urlencode "method=login" \
        --data-urlencode "account=${NAS_USER}" \
        --data-urlencode "passwd=${NAS_PASS}" \
        --data-urlencode "format=sid" 2>/dev/null \
        | grep -o '"sid":"[^"]*"' | cut -d'"' -f4)
    if [ -z "$sid" ]; then
        echo -e "${RED}ERROR: Failed to login to NAS API.${NC}" >&2
        exit 1
    fi
    echo "$sid"
}

nas_list_files() {
    local folder_path="$1"
    local sid="$2"
    local limit="${3:-2000}"
    curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi" \
        --data-urlencode "api=SYNO.FileStation.List" \
        --data-urlencode "version=2" \
        --data-urlencode "method=list" \
        --data-urlencode "folder_path=${folder_path}" \
        --data-urlencode "additional=[\"size\"]" \
        --data-urlencode "limit=${limit}" \
        --data-urlencode "_sid=${sid}" 2>/dev/null
}

nas_download_file() {
    local file_path="$1"
    local sid="$2"
    curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.Download&version=2&method=download&path=$(python3 -c "import urllib.parse; print(urllib.parse.quote('${file_path}', safe=''))")&mode=download&_sid=${sid}" 2>/dev/null
}

download_folder_to_server() {
    local nas_folder="$1"
    local server_dest="$2"
    local sid="$3"
    local file_pattern="${4:-}"  # optional: "fastq.gz" or "fast5" or "txt"

    echo -e "  ${CYAN}Listing files in: ${nas_folder}${NC}"

    # Get file list
    local file_list
    file_list=$(nas_list_files "$nas_folder" "$sid")

    local total
    total=$(echo "$file_list" | python3 -c "import sys,json; print(json.load(sys.stdin)['data']['total'])" 2>/dev/null || echo "0")

    if [ "$total" = "0" ]; then
        echo -e "  ${YELLOW}WARNING: No files found in $nas_folder${NC}"
        return
    fi

    echo -e "  ${GREEN}Found $total files${NC}"

    # Create destination directory on server
    ssh ${SERVER_USER}@${SERVER_HOST} "mkdir -p '$server_dest'" 2>/dev/null

    # Download each file and stream to server
    local count=0
    local offset=0
    local batch_size=100

    while [ "$offset" -lt "$total" ]; do
        local batch
        batch=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi" \
            --data-urlencode "api=SYNO.FileStation.List" \
            --data-urlencode "version=2" \
            --data-urlencode "method=list" \
            --data-urlencode "folder_path=${nas_folder}" \
            --data-urlencode "additional=[\"size\"]" \
            --data-urlencode "offset=${offset}" \
            --data-urlencode "limit=${batch_size}" \
            --data-urlencode "_sid=${sid}" 2>/dev/null)

        # Extract file paths and sizes
        local files
        files=$(echo "$batch" | python3 -c "
import sys, json
data = json.load(sys.stdin)
for f in data['data']['files']:
    if not f['isdir']:
        print(f['path'] + '|' + f['name'] + '|' + str(f['additional']['size']))
" 2>/dev/null)

        while IFS='|' read -r fpath fname fsize; do
            [ -z "$fpath" ] && continue

            # Apply file pattern filter
            if [ -n "$file_pattern" ]; then
                case "$fname" in
                    *${file_pattern}*) ;;
                    *) continue ;;
                esac
            fi

            count=$((count + 1))
            local fsize_mb=$((fsize / 1048576))
            echo -ne "  [${count}] ${fname} (${fsize_mb} MB)... "

            # Download from NAS and stream to server via SSH
            curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.Download&version=2&method=download&path=$(python3 -c "import urllib.parse; print(urllib.parse.quote('${fpath}', safe=''))")&mode=download&_sid=${sid}" 2>/dev/null \
                | ssh ${SERVER_USER}@${SERVER_HOST} "cat > '${server_dest}/${fname}'" 2>/dev/null

            echo -e "${GREEN}OK${NC}"
        done <<< "$files"

        offset=$((offset + batch_size))
    done

    echo -e "  ${GREEN}Transferred $count files to ${server_dest}${NC}"
}

# ==============================================================
# MAIN
# ==============================================================

echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN} K-CHOPORE NAS Data Transfer${NC}"
echo -e "${GREEN} NAS -> Windows -> Server${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""

# Login to NAS
echo -e "${YELLOW}[Step 1] Logging in to Synology NAS...${NC}"
SID=$(nas_login)
echo -e "  ${GREEN}Session established.${NC}"
echo ""

# Verify server is reachable
echo -e "${YELLOW}[Step 2] Testing server connectivity...${NC}"
ssh ${SERVER_USER}@${SERVER_HOST} "echo SERVER_OK && mkdir -p '${SERVER_BASE}'" 2>/dev/null
echo -e "  ${GREEN}Server reachable.${NC}"
echo ""

# Define samples and their data
declare -A SAMPLES
# format: SAMPLES["nas_dir"]="run_subdir|data_type"
SAMPLES["WT_C_R1"]="no_sample/20220224_1639_MC-112869_FAR90122_f656439d|fastq"
SAMPLES["WT_C_R2"]="no_sample/20220323_1648_MC-112869_FAR91957_0263d5a0|fastq"
SAMPLES["WT_C_R3"]="no_sample/20220301_1656_MC-112869_FAR90189_a35665b9|fastq"
SAMPLES["WT_AA_R1"]="20220221_1734_MC-112869_FAR92050_00e32dc0|fastq"
SAMPLES["WT_AA_R2"]="no_sample/20220316_1902_MC-112869_FAR91120_624c5f7d|fastq"
SAMPLES["WT_AA_R3"]="no_sample/20220315_1614_MC-112869_FAR90098_c891954b|fastq"
SAMPLES["WT_AA_R3_2"]="no_sample/20220316_1636_MC-112869_FAR90098_6e408c50|fastq"
SAMPLES["anac017-1_C_R1"]="no_sample/20220330_1936_MC-112869_FAR92015_6ac8671e|fastq"
SAMPLES["anac017-1_C_R2"]="no_sample/20220404_1440_MC-112869_FAR91811_a135af49|fastq"
SAMPLES["anac017-1_C_R3"]="no_sample/20220405_1610_MC-112869_FAR91811_d5fba51e|fastq"
SAMPLES["anac017-1_AA_R1"]="no_sample/20220328_1730_MC-112869_FAR90074_970aa1f5|fastq"
SAMPLES["anac017-1_AA_R2"]="no_sample/20220930_1148_MC-112869_FAR81498_9fdc7765|fast5"
SAMPLES["anac017-1_AA_R3"]="no_sample/20220929_1527_MC-112869_FAU09642_5f85f08b|fast5"
SAMPLES["anac017_AA_R2-2"]="no_sample/20221215_1232_MC-112869_FAT21048_d754e27f|fast5"

TOTAL_SAMPLES=${#SAMPLES[@]}
CURRENT=0

echo -e "${YELLOW}[Step 3] Transferring data for $TOTAL_SAMPLES sample directories...${NC}"
echo ""

for nas_dir in "${!SAMPLES[@]}"; do
    IFS='|' read -r run_subdir data_type <<< "${SAMPLES[$nas_dir]}"
    CURRENT=$((CURRENT + 1))
    NAS_RUN_PATH="${NAS_BASE}/${nas_dir}/${run_subdir}"
    SERVER_RUN_PATH="${SERVER_BASE}/${nas_dir}/${run_subdir}"

    echo -e "${GREEN}[${CURRENT}/${TOTAL_SAMPLES}] ${nas_dir} (${data_type})${NC}"

    if [ "$data_type" = "fastq" ]; then
        # Download fastq_pass/ files
        echo -e "  ${YELLOW}Downloading fastq_pass...${NC}"
        download_folder_to_server "${NAS_RUN_PATH}/fastq_pass" "${SERVER_RUN_PATH}/fastq_pass" "$SID"

        # Download sequencing_summary
        echo -e "  ${YELLOW}Downloading sequencing_summary...${NC}"
        download_folder_to_server "${NAS_RUN_PATH}" "${SERVER_RUN_PATH}" "$SID" "sequencing_summary"
    else
        # FAST5-only: download fast5/ directory
        echo -e "  ${YELLOW}Downloading fast5...${NC}"
        download_folder_to_server "${NAS_RUN_PATH}/fast5" "${SERVER_RUN_PATH}/fast5" "$SID"

        # Download sequencing_summary
        echo -e "  ${YELLOW}Downloading sequencing_summary...${NC}"
        download_folder_to_server "${NAS_RUN_PATH}" "${SERVER_RUN_PATH}" "$SID" "sequencing_summary"
    fi

    echo ""
done

# Final verification
echo -e "${YELLOW}[Step 4] Verifying data on server...${NC}"
ssh ${SERVER_USER}@${SERVER_HOST} "du -sh ${SERVER_BASE}/*/ 2>/dev/null"

echo ""
echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN} Data transfer complete!${NC}"
echo -e "${GREEN}============================================${NC}"
echo ""
echo "Next steps:"
echo "  1. SSH to server:  ssh ${SERVER_USER}@${SERVER_HOST}"
echo "  2. Run pipeline:   cd /home/usuario2/pelamovic/kchopore && ./run_pipeline.sh"
echo ""
