#!/usr/bin/env python3
"""
K-CHOPORE NAS Data Transfer Script (v2 - batch mode)
Downloads each sample's data from NAS to local temp dir,
then uploads the whole sample to server via a single scp/rsync.

Usage: python scripts/transfer_data.py [sample_index]
  sample_index: optional, 0-based index to start from (for resuming)
"""
import os
import sys
import json
import subprocess
import urllib.request
import urllib.parse
import ssl
import time
import shutil
import tempfile

# Configuration
NAS_HOST = "valmeilab.synology.me"
NAS_PORT = "5001"
NAS_USER = "pelayo"
NAS_PASS = "Pelamovic39@"
NAS_BASE = "/HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS"

SERVER = "usuario2@156.35.42.17"
SERVER_BASE = "/home/usuario2/pelamovic/kchopore/data/raw"

# Local temp directory for downloads
LOCAL_TEMP = os.path.join(tempfile.gettempdir(), "kchopore_transfer")

# SSL context (ignore self-signed cert)
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# All samples
ALL_SAMPLES = [
    # (nas_dir, run_subdir, data_type)
    ("WT_C_R1", "no_sample/20220224_1639_MC-112869_FAR90122_f656439d", "fastq"),
    ("WT_C_R2", "no_sample/20220323_1648_MC-112869_FAR91957_0263d5a0", "fastq"),
    ("WT_C_R3", "no_sample/20220301_1656_MC-112869_FAR90189_a35665b9", "fastq"),
    ("WT_AA_R1", "20220221_1734_MC-112869_FAR92050_00e32dc0", "fastq"),
    ("WT_AA_R2", "no_sample/20220316_1902_MC-112869_FAR91120_624c5f7d", "fastq"),
    ("WT_AA_R3", "no_sample/20220315_1614_MC-112869_FAR90098_c891954b", "fastq"),
    ("WT_AA_R3_2", "no_sample/20220316_1636_MC-112869_FAR90098_6e408c50", "fastq"),
    ("anac017-1_C_R1", "no_sample/20220330_1936_MC-112869_FAR92015_6ac8671e", "fastq"),
    ("anac017-1_C_R2", "no_sample/20220404_1440_MC-112869_FAR91811_a135af49", "fastq"),
    ("anac017-1_C_R3", "no_sample/20220405_1610_MC-112869_FAR91811_d5fba51e", "fastq"),
    ("anac017-1_AA_R1", "no_sample/20220328_1730_MC-112869_FAR90074_970aa1f5", "fastq"),
    ("anac017-1_AA_R2", "no_sample/20220930_1148_MC-112869_FAR81498_9fdc7765", "fast5"),
    ("anac017-1_AA_R3", "no_sample/20220929_1527_MC-112869_FAU09642_5f85f08b", "fast5"),
    ("anac017_AA_R2-2", "no_sample/20221215_1232_MC-112869_FAT21048_d754e27f", "fast5"),
]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def nas_api_get(params, sid=None):
    """Make a GET request to Synology File Station API."""
    if sid:
        params["_sid"] = sid
    query = urllib.parse.urlencode(params, safe="[]")
    url = f"https://{NAS_HOST}:{NAS_PORT}/webapi/entry.cgi?{query}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, context=ctx, timeout=60) as resp:
        return json.loads(resp.read())


def nas_login():
    """Login to Synology and return session ID."""
    params = {
        "api": "SYNO.API.Auth",
        "version": "6",
        "method": "login",
        "account": NAS_USER,
        "passwd": NAS_PASS,
        "format": "sid",
    }
    query = urllib.parse.urlencode(params)
    url = f"https://{NAS_HOST}:{NAS_PORT}/webapi/auth.cgi?{query}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, context=ctx, timeout=30) as resp:
        result = json.loads(resp.read())
    if not result.get("success"):
        raise Exception(f"NAS login failed: {result}")
    return result["data"]["sid"]


def nas_list_all_files(folder_path, sid):
    """List ALL files in a NAS folder (paginated)."""
    all_files = []
    offset = 0
    total = None
    while total is None or offset < total:
        data = nas_api_get({
            "api": "SYNO.FileStation.List",
            "version": "2",
            "method": "list",
            "folder_path": folder_path,
            "additional": '["size"]',
            "offset": str(offset),
            "limit": "500",
        }, sid=sid)
        if not data.get("success"):
            log(f"  WARNING: List failed for {folder_path}: {data}")
            break
        total = data["data"]["total"]
        for f in data["data"]["files"]:
            if not f["isdir"]:
                all_files.append({
                    "path": f["path"],
                    "name": f["name"],
                    "size": f["additional"]["size"],
                })
        offset += 500
    return all_files


def nas_download_file(nas_path, local_path, sid):
    """Download a single file from NAS to local disk."""
    encoded = urllib.parse.quote(nas_path, safe="")
    url = f"https://{NAS_HOST}:{NAS_PORT}/webapi/entry.cgi?api=SYNO.FileStation.Download&version=2&method=download&path={encoded}&mode=download&_sid={sid}"
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, context=ctx, timeout=120) as resp:
        with open(local_path, 'wb') as f:
            while True:
                chunk = resp.read(65536)
                if not chunk:
                    break
                f.write(chunk)


def scp_to_server(local_dir, server_dir):
    """Upload a local directory to the server using scp -r."""
    # Create parent dir on server
    subprocess.run(
        ["ssh", "-o", "ConnectTimeout=30", SERVER, f"mkdir -p '{server_dir}'"],
        capture_output=True, timeout=60
    )
    # Use scp -r to upload entire directory
    result = subprocess.run(
        ["scp", "-r", "-o", "ConnectTimeout=30", "-q",
         local_dir, f"{SERVER}:{server_dir}/"],
        capture_output=True, text=True, timeout=3600  # 1h per sample max
    )
    return result.returncode == 0


def transfer_sample(nas_dir, run_subdir, data_type, sid):
    """Download one sample from NAS to local, then upload to server."""
    nas_run = f"{NAS_BASE}/{nas_dir}/{run_subdir}"
    local_sample = os.path.join(LOCAL_TEMP, nas_dir, run_subdir)

    # Determine what to download
    if data_type == "fastq":
        folders_to_download = [("fastq_pass", f"{nas_run}/fastq_pass")]
    else:
        folders_to_download = [("fast5", f"{nas_run}/fast5")]

    total_files = 0
    total_bytes = 0

    for folder_name, nas_folder in folders_to_download:
        local_folder = os.path.join(local_sample, folder_name)
        os.makedirs(local_folder, exist_ok=True)

        log(f"  Listing {folder_name}...")
        files = nas_list_all_files(nas_folder, sid)
        log(f"  Found {len(files)} files in {folder_name}")

        for i, f in enumerate(files, 1):
            local_path = os.path.join(local_folder, f["name"])
            if os.path.exists(local_path) and os.path.getsize(local_path) == f["size"]:
                # Skip already downloaded
                total_files += 1
                total_bytes += f["size"]
                continue

            if i == 1 or i % 50 == 0 or i == len(files):
                mb = f["size"] / (1024*1024)
                log(f"  Downloading [{i}/{len(files)}] {f['name']} ({mb:.1f} MB)")

            try:
                nas_download_file(f["path"], local_path, sid)
                total_files += 1
                total_bytes += f["size"]
            except Exception as e:
                log(f"  ERROR downloading {f['name']}: {e}")

    # Also download sequencing_summary
    log(f"  Downloading sequencing_summary...")
    os.makedirs(local_sample, exist_ok=True)
    run_files = nas_list_all_files(nas_run, sid)
    for f in run_files:
        if "sequencing_summary" in f["name"]:
            local_path = os.path.join(local_sample, f["name"])
            mb = f["size"] / (1024*1024)
            log(f"  Summary: {f['name']} ({mb:.1f} MB)")
            try:
                nas_download_file(f["path"], local_path, sid)
                total_files += 1
                total_bytes += f["size"]
            except Exception as e:
                log(f"  ERROR downloading summary: {e}")

    gb = total_bytes / (1024**3)
    log(f"  Downloaded {total_files} files ({gb:.2f} GB) locally")

    # Upload to server in one go
    server_dest = f"{SERVER_BASE}/{nas_dir}"
    log(f"  Uploading to server: {server_dest}")
    # Upload the nas_dir folder (which contains run_subdir/fastq_pass etc)
    local_nas_dir = os.path.join(LOCAL_TEMP, nas_dir)
    success = scp_to_server(local_nas_dir, SERVER_BASE)
    if success:
        log(f"  Upload complete!")
    else:
        log(f"  WARNING: Upload may have had issues. Check server.")

    # Clean local temp for this sample to save disk space
    shutil.rmtree(local_nas_dir, ignore_errors=True)
    log(f"  Local temp cleaned.")

    return total_files, total_bytes


def main():
    start_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0

    print("=" * 60, flush=True)
    print(" K-CHOPORE NAS Data Transfer (v2 - batch mode)", flush=True)
    print(" Download locally, then scp to server", flush=True)
    print(f" Local temp: {LOCAL_TEMP}", flush=True)
    print("=" * 60, flush=True)
    print(flush=True)

    os.makedirs(LOCAL_TEMP, exist_ok=True)

    log("Logging in to NAS...")
    sid = nas_login()
    log("Session established.")
    print(flush=True)

    grand_files = 0
    grand_bytes = 0

    for i, (nas_dir, run_subdir, data_type) in enumerate(ALL_SAMPLES):
        if i < start_idx:
            log(f"[{i+1}/{len(ALL_SAMPLES)}] SKIPPED {nas_dir}")
            continue

        log(f"[{i+1}/{len(ALL_SAMPLES)}] {nas_dir} ({data_type})")

        # Re-login every 3 samples
        if i > 0 and i % 3 == 0:
            log("  Refreshing NAS session...")
            try:
                sid = nas_login()
            except:
                time.sleep(5)
                sid = nas_login()

        try:
            nf, nb = transfer_sample(nas_dir, run_subdir, data_type, sid)
            grand_files += nf
            grand_bytes += nb
        except Exception as e:
            log(f"  FAILED: {e}")
            log(f"  You can resume with: python transfer_data.py {i}")

        print(flush=True)

    grand_gb = grand_bytes / (1024**3)
    print("=" * 60, flush=True)
    log(f"TRANSFER COMPLETE: {grand_files} files, {grand_gb:.2f} GB")
    print("=" * 60, flush=True)

    # Verify on server
    log("Checking server:")
    result = subprocess.run(
        ["ssh", "-o", "ConnectTimeout=30", SERVER,
         f"du -sh {SERVER_BASE}/*/ 2>/dev/null; echo '---'; df -h /"],
        capture_output=True, text=True, timeout=30
    )
    print(result.stdout, flush=True)


if __name__ == "__main__":
    main()
