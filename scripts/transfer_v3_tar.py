#!/usr/bin/env python3
"""
K-CHOPORE: Transfer data from NAS to server via tar stream.
NAS API -> download to local -> tar | ssh cat > server

Much faster than per-file SCP because it uses a single SSH connection.

Usage: python scripts/transfer_v3_tar.py [start_index]
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

LOCAL_TEMP = os.path.join(tempfile.gettempdir(), "kchopore_transfer")

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# All 14 samples
ALL_SAMPLES = [
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


def nas_login():
    params = urllib.parse.urlencode({
        "api": "SYNO.API.Auth", "version": "6", "method": "login",
        "account": NAS_USER, "passwd": NAS_PASS, "format": "sid",
    })
    url = f"https://{NAS_HOST}:{NAS_PORT}/webapi/auth.cgi?{params}"
    with urllib.request.urlopen(urllib.request.Request(url), context=ctx, timeout=30) as resp:
        result = json.loads(resp.read())
    if not result.get("success"):
        raise Exception(f"NAS login failed: {result}")
    return result["data"]["sid"]


def nas_list_all_files(folder_path, sid):
    all_files = []
    offset = 0
    total = None
    while total is None or offset < total:
        params = urllib.parse.urlencode({
            "api": "SYNO.FileStation.List", "version": "2", "method": "list",
            "folder_path": folder_path, "additional": '["size"]',
            "offset": str(offset), "limit": "500", "_sid": sid,
        }, safe='[]')
        url = f"https://{NAS_HOST}:{NAS_PORT}/webapi/entry.cgi?{params}"
        with urllib.request.urlopen(urllib.request.Request(url), context=ctx, timeout=60) as resp:
            data = json.loads(resp.read())
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
    encoded = urllib.parse.quote(nas_path, safe="")
    url = (f"https://{NAS_HOST}:{NAS_PORT}/webapi/entry.cgi?"
           f"api=SYNO.FileStation.Download&version=2&method=download"
           f"&path={encoded}&mode=download&_sid={sid}")
    with urllib.request.urlopen(urllib.request.Request(url), context=ctx, timeout=300) as resp:
        with open(local_path, 'wb') as f:
            while True:
                chunk = resp.read(65536)
                if not chunk:
                    break
                f.write(chunk)


def tar_upload_to_server(local_dir, server_dest):
    """Upload directory to server using tar pipe (much faster than scp -r)."""
    # Create parent dir on server
    subprocess.run(
        ["ssh", "-o", "ConnectTimeout=30", SERVER, f"mkdir -p '{server_dest}'"],
        capture_output=True, timeout=60
    )
    # Tar the local directory and pipe through SSH
    # tar cf - -C parent_dir dirname | ssh server "tar xf - -C dest"
    parent = os.path.dirname(local_dir)
    dirname = os.path.basename(local_dir)

    cmd = f'tar cf - -C "{parent}" "{dirname}" | ssh -o ConnectTimeout=30 {SERVER} "tar xf - -C \'{server_dest}\'"'
    log(f"  Tar-piping {dirname} to server...")
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=7200
    )
    if result.returncode != 0:
        log(f"  Tar pipe error: {result.stderr[:300]}")
    return result.returncode == 0


def download_and_upload_sample(nas_dir, run_subdir, data_type, sid):
    """Download one sample from NAS, then tar-pipe to server."""
    nas_run = f"{NAS_BASE}/{nas_dir}/{run_subdir}"
    local_sample = os.path.join(LOCAL_TEMP, nas_dir, run_subdir)
    total_files = 0
    total_bytes = 0

    # Determine what to download
    if data_type == "fastq":
        # Download fastq_pass
        folder = f"{nas_run}/fastq_pass"
        local_folder = os.path.join(local_sample, "fastq_pass")
        os.makedirs(local_folder, exist_ok=True)
        log(f"  Listing fastq_pass...")
        files = nas_list_all_files(folder, sid)
        log(f"  Found {len(files)} FASTQ files")

        for i, f in enumerate(files, 1):
            local_path = os.path.join(local_folder, f["name"])
            if os.path.exists(local_path) and os.path.getsize(local_path) == f["size"]:
                total_files += 1
                total_bytes += f["size"]
                continue
            if i == 1 or i % 100 == 0 or i == len(files):
                mb = f["size"] / (1024*1024)
                log(f"  Downloading [{i}/{len(files)}] {f['name']} ({mb:.1f} MB)")
            try:
                nas_download_file(f["path"], local_path, sid)
                total_files += 1
                total_bytes += f["size"]
            except Exception as e:
                log(f"  ERROR: {f['name']}: {e}")

    else:
        # Download fast5
        folder = f"{nas_run}/fast5"
        local_folder = os.path.join(local_sample, "fast5")
        os.makedirs(local_folder, exist_ok=True)
        log(f"  Listing fast5...")
        files = nas_list_all_files(folder, sid)
        log(f"  Found {len(files)} FAST5 files")

        for i, f in enumerate(files, 1):
            local_path = os.path.join(local_folder, f["name"])
            if os.path.exists(local_path) and os.path.getsize(local_path) == f["size"]:
                total_files += 1
                total_bytes += f["size"]
                continue
            if i == 1 or i % 20 == 0 or i == len(files):
                mb = f["size"] / (1024*1024)
                log(f"  Downloading [{i}/{len(files)}] {f['name']} ({mb:.1f} MB)")
            try:
                nas_download_file(f["path"], local_path, sid)
                total_files += 1
                total_bytes += f["size"]
            except Exception as e:
                log(f"  ERROR: {f['name']}: {e}")

    # Download sequencing_summary
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
                log(f"  ERROR summary: {e}")

    gb = total_bytes / (1024**3)
    log(f"  Downloaded {total_files} files ({gb:.2f} GB)")

    # Upload via tar pipe
    local_nas_dir = os.path.join(LOCAL_TEMP, nas_dir)
    success = tar_upload_to_server(local_nas_dir, SERVER_BASE)
    if success:
        log(f"  Upload complete!")
    else:
        log(f"  WARNING: Upload may have had issues.")

    # Clean local temp
    shutil.rmtree(local_nas_dir, ignore_errors=True)
    log(f"  Local temp cleaned.")
    return total_files, total_bytes


def main():
    start_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0

    print("=" * 60, flush=True)
    print(" K-CHOPORE Data Transfer v3 (tar pipe)", flush=True)
    print(f" NAS -> Windows -> tar | ssh -> Server", flush=True)
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
            nf, nb = download_and_upload_sample(nas_dir, run_subdir, data_type, sid)
            grand_files += nf
            grand_bytes += nb
        except Exception as e:
            log(f"  FAILED: {e}")
            import traceback
            traceback.print_exc()
            log(f"  Resume with: python scripts/transfer_v3_tar.py {i}")

        print(flush=True)

    grand_gb = grand_bytes / (1024**3)
    print("=" * 60, flush=True)
    log(f"TRANSFER COMPLETE: {grand_files} files, {grand_gb:.2f} GB")
    print("=" * 60, flush=True)

    # Verify on server
    log("Checking server:")
    result = subprocess.run(
        ["ssh", "-o", "ConnectTimeout=30", SERVER,
         f"du -sh {SERVER_BASE}/*/ 2>/dev/null; echo '---'; df -h / /media/usuario2/ssd4TB"],
        capture_output=True, text=True, timeout=30
    )
    print(result.stdout, flush=True)


if __name__ == "__main__":
    main()
