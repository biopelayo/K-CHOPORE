#!/bin/bash
# Transfer all sample data from NAS to server
# Run from Windows/Git Bash
set -euo pipefail

NAS_HOST="valmeilab.synology.me"
NAS_PORT="5001"
NAS_BASE="/HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS"
SERVER="usuario2@156.35.42.17"
SERVER_BASE="/home/usuario2/pelamovic/kchopore/data/raw"

# Login to NAS
echo "[$(date +%H:%M:%S)] Logging in to NAS..."
SID=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/auth.cgi" \
    --data-urlencode "api=SYNO.API.Auth" \
    --data-urlencode "version=6" \
    --data-urlencode "method=login" \
    --data-urlencode "account=pelayo" \
    --data-urlencode "passwd=Pelamovic39@" \
    --data-urlencode "format=sid" 2>/dev/null \
    | grep -o '"sid":"[^"]*"' | cut -d'"' -f4)
echo "[$(date +%H:%M:%S)] SID: ${SID:0:10}..."

url_encode_path() {
    local path="$1"
    echo "$path" | sed 's|/|%2F|g' | sed 's| |%20|g' | sed 's|-|%2D|g'
}

download_file_to_server() {
    local nas_path="$1"
    local server_path="$2"
    local encoded
    encoded=$(url_encode_path "$nas_path")
    curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.Download&version=2&method=download&path=${encoded}&mode=download&_sid=${SID}" 2>/dev/null \
        | ssh $SERVER "cat > '${server_path}'" 2>/dev/null
}

transfer_sample_fastq() {
    local nas_dir="$1"
    local run_subdir="$2"
    local nas_folder="${NAS_BASE}/${nas_dir}/${run_subdir}/fastq_pass"
    local server_folder="${SERVER_BASE}/${nas_dir}/${run_subdir}/fastq_pass"

    # Create dir on server
    ssh $SERVER "mkdir -p '${server_folder}'" 2>/dev/null

    # List files
    local encoded_folder
    encoded_folder=$(url_encode_path "$nas_folder")
    local total offset=0 batch_size=200

    # Get total count first
    total=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.List&version=2&method=list&folder_path=${encoded_folder}&limit=1&_sid=${SID}" 2>/dev/null \
        | grep -o '"total":[0-9]*' | cut -d: -f2)

    echo "    fastq_pass: $total files"

    local count=0
    while [ "$offset" -lt "$total" ]; do
        # Get batch of filenames
        local batch
        batch=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.List&version=2&method=list&folder_path=${encoded_folder}&offset=${offset}&limit=${batch_size}&_sid=${SID}" 2>/dev/null)

        # Extract file paths - use simple grep/sed since python3 may not work in Git Bash
        local files
        files=$(echo "$batch" | grep -o '"path":"[^"]*"' | cut -d'"' -f4 | grep -v '/$')

        for fpath in $files; do
            local fname
            fname=$(basename "$fpath")
            count=$((count + 1))

            if [ $((count % 50)) -eq 0 ] || [ "$count" -eq 1 ]; then
                echo "    [${count}/${total}] ${fname}"
            fi

            download_file_to_server "$fpath" "${server_folder}/${fname}"
        done

        offset=$((offset + batch_size))
    done
    echo "    Done: ${count} files transferred"
}

transfer_sample_summary() {
    local nas_dir="$1"
    local run_subdir="$2"
    local nas_folder="${NAS_BASE}/${nas_dir}/${run_subdir}"
    local server_folder="${SERVER_BASE}/${nas_dir}/${run_subdir}"

    ssh $SERVER "mkdir -p '${server_folder}'" 2>/dev/null

    # List files in run directory, find sequencing_summary
    local encoded_folder
    encoded_folder=$(url_encode_path "$nas_folder")
    local file_list
    file_list=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.List&version=2&method=list&folder_path=${encoded_folder}&limit=20&_sid=${SID}" 2>/dev/null)

    local summary_paths
    summary_paths=$(echo "$file_list" | grep -o '"path":"[^"]*sequencing_summary[^"]*"' | cut -d'"' -f4)

    for fpath in $summary_paths; do
        local fname
        fname=$(basename "$fpath")
        echo "    summary: ${fname}"
        download_file_to_server "$fpath" "${server_folder}/${fname}"
    done
}

transfer_sample_fast5() {
    local nas_dir="$1"
    local run_subdir="$2"
    local nas_folder="${NAS_BASE}/${nas_dir}/${run_subdir}/fast5"
    local server_folder="${SERVER_BASE}/${nas_dir}/${run_subdir}/fast5"

    ssh $SERVER "mkdir -p '${server_folder}'" 2>/dev/null

    local encoded_folder
    encoded_folder=$(url_encode_path "$nas_folder")
    local total offset=0 batch_size=200

    total=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.List&version=2&method=list&folder_path=${encoded_folder}&limit=1&_sid=${SID}" 2>/dev/null \
        | grep -o '"total":[0-9]*' | cut -d: -f2)

    echo "    fast5: $total files"

    local count=0
    while [ "$offset" -lt "$total" ]; do
        local batch
        batch=$(curl -sk "https://${NAS_HOST}:${NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.List&version=2&method=list&folder_path=${encoded_folder}&offset=${offset}&limit=${batch_size}&_sid=${SID}" 2>/dev/null)

        local files
        files=$(echo "$batch" | grep -o '"path":"[^"]*"' | cut -d'"' -f4 | grep -v '/$')

        for fpath in $files; do
            local fname
            fname=$(basename "$fpath")
            count=$((count + 1))

            if [ $((count % 20)) -eq 0 ] || [ "$count" -eq 1 ]; then
                echo "    [${count}/${total}] ${fname}"
            fi

            download_file_to_server "$fpath" "${server_folder}/${fname}"
        done

        offset=$((offset + batch_size))
    done
    echo "    Done: ${count} files transferred"
}

# ---- MAIN ----
echo ""
echo "=== K-CHOPORE Data Transfer: NAS -> Windows -> Server ==="
echo ""

# Pre-basecalled samples (fastq_pass + summary)
echo "=== PRE-BASECALLED SAMPLES (FASTQ) ==="
declare -a FASTQ_SAMPLES=(
    "WT_C_R1|no_sample/20220224_1639_MC-112869_FAR90122_f656439d"
    "WT_C_R2|no_sample/20220323_1648_MC-112869_FAR91957_0263d5a0"
    "WT_C_R3|no_sample/20220301_1656_MC-112869_FAR90189_a35665b9"
    "WT_AA_R1|20220221_1734_MC-112869_FAR92050_00e32dc0"
    "WT_AA_R2|no_sample/20220316_1902_MC-112869_FAR91120_624c5f7d"
    "WT_AA_R3|no_sample/20220315_1614_MC-112869_FAR90098_c891954b"
    "WT_AA_R3_2|no_sample/20220316_1636_MC-112869_FAR90098_6e408c50"
    "anac017-1_C_R1|no_sample/20220330_1936_MC-112869_FAR92015_6ac8671e"
    "anac017-1_C_R2|no_sample/20220404_1440_MC-112869_FAR91811_a135af49"
    "anac017-1_C_R3|no_sample/20220405_1610_MC-112869_FAR91811_d5fba51e"
    "anac017-1_AA_R1|no_sample/20220328_1730_MC-112869_FAR90074_970aa1f5"
)

IDX=0
TOTAL_FASTQ=${#FASTQ_SAMPLES[@]}
for entry in "${FASTQ_SAMPLES[@]}"; do
    IFS='|' read -r nas_dir run_subdir <<< "$entry"
    IDX=$((IDX + 1))
    echo ""
    echo "[$(date +%H:%M:%S)] [${IDX}/${TOTAL_FASTQ}] ${nas_dir}"
    transfer_sample_fastq "$nas_dir" "$run_subdir"
    transfer_sample_summary "$nas_dir" "$run_subdir"
done

# FAST5-only samples (fast5 + summary)
echo ""
echo "=== FAST5-ONLY SAMPLES (need basecalling) ==="
declare -a FAST5_SAMPLES=(
    "anac017-1_AA_R2|no_sample/20220930_1148_MC-112869_FAR81498_9fdc7765"
    "anac017-1_AA_R3|no_sample/20220929_1527_MC-112869_FAU09642_5f85f08b"
    "anac017_AA_R2-2|no_sample/20221215_1232_MC-112869_FAT21048_d754e27f"
)

IDX=0
TOTAL_FAST5=${#FAST5_SAMPLES[@]}
for entry in "${FAST5_SAMPLES[@]}"; do
    IFS='|' read -r nas_dir run_subdir <<< "$entry"
    IDX=$((IDX + 1))
    echo ""
    echo "[$(date +%H:%M:%S)] [${IDX}/${TOTAL_FAST5}] ${nas_dir} (FAST5)"
    transfer_sample_fast5 "$nas_dir" "$run_subdir"
    transfer_sample_summary "$nas_dir" "$run_subdir"
done

echo ""
echo "[$(date +%H:%M:%S)] === TRANSFER COMPLETE ==="
echo ""
ssh $SERVER "du -sh ${SERVER_BASE}/*/" 2>/dev/null
echo ""
ssh $SERVER "df -h /" 2>/dev/null
