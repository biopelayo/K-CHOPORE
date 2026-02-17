#!/usr/bin/env python3
"""
K-CHOPORE Transfer v4: Two-step approach (sleep-resistant)

Step 1: Download ALL samples from NAS to local Windows disk
Step 2: Upload ALL at once to server via tar pipe

If computer sleeps during step 1, just re-run (has resume).
If computer sleeps during step 2, just re-run step 2 (data is local).

Usage:
  python scripts/transfer_v4_twostep.py download [start_index]
  python scripts/transfer_v4_twostep.py upload
  python scripts/transfer_v4_twostep.py status
"""
import os
import sys
import json
import subprocess
import urllib.request
import urllib.parse
import ssl
import time
import glob

# Configuration
NAS_HOST = "valmeilab.synology.me"
NAS_PORT = "5001"
NAS_USER = "pelayo"
NAS_PASS = "Pelamovic39@"
NAS_BASE = "/HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS"

SERVER = "usuario2@156.35.42.17"
SERVER_BASE = "/home/usuario2/pelamovic/kchopore/data/raw"

# Local storage for downloaded data (moved to E: for space)
LOCAL_DATA = os.path.join("E:", os.sep, "kchopore_nas_data")

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

ALL_SAMPLES = [
    ("WT_C_R1",         "no_sample/20220224_1639_MC-112869_FAR90122_f656439d", "fastq"),
    ("WT_C_R2",         "no_sample/20220323_1648_MC-112869_FAR91957_0263d5a0", "fastq"),
    ("WT_C_R3",         "no_sample/20220301_1656_MC-112869_FAR90189_a35665b9", "fastq"),
    ("WT_AA_R1",        "20220221_1734_MC-112869_FAR92050_00e32dc0",           "fastq"),
    ("WT_AA_R2",        "no_sample/20220316_1902_MC-112869_FAR91120_624c5f7d", "fastq"),
    ("WT_AA_R3",        "no_sample/20220315_1614_MC-112869_FAR90098_c891954b", "fastq"),
    ("WT_AA_R3_2",      "no_sample/20220316_1636_MC-112869_FAR90098_6e408c50", "fastq"),
    ("anac017-1_C_R1",  "no_sample/20220330_1936_MC-112869_FAR92015_6ac8671e", "fastq"),
    ("anac017-1_C_R2",  "no_sample/20220404_1440_MC-112869_FAR91811_a135af49", "fastq"),
    ("anac017-1_C_R3",  "no_sample/20220405_1610_MC-112869_FAR91811_d5fba51e", "fastq"),
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


def do_download(start_idx=0):
    """Step 1: Download all samples from NAS to local disk."""
    print("=" * 60, flush=True)
    print(f" STEP 1: Download from NAS to local disk", flush=True)
    print(f" Local: {LOCAL_DATA}", flush=True)
    print("=" * 60, flush=True)

    os.makedirs(LOCAL_DATA, exist_ok=True)

    log("Logging in to NAS...")
    sid = nas_login()
    log("OK")

    for i, (nas_dir, run_subdir, data_type) in enumerate(ALL_SAMPLES):
        if i < start_idx:
            log(f"[{i+1}/14] SKIP {nas_dir}")
            continue

        log(f"[{i+1}/14] {nas_dir} ({data_type})")
        nas_run = f"{NAS_BASE}/{nas_dir}/{run_subdir}"
        local_sample = os.path.join(LOCAL_DATA, nas_dir, run_subdir)

        # Re-login every 4 samples
        if i > 0 and i % 4 == 0:
            try:
                sid = nas_login()
            except:
                time.sleep(5)
                sid = nas_login()

        # Download data folder
        if data_type == "fastq":
            folder = f"{nas_run}/fastq_pass"
            local_folder = os.path.join(local_sample, "fastq_pass")
        else:
            folder = f"{nas_run}/fast5"
            local_folder = os.path.join(local_sample, "fast5")

        os.makedirs(local_folder, exist_ok=True)
        files = nas_list_all_files(folder, sid)
        log(f"  {len(files)} files to download")

        downloaded = 0
        skipped = 0
        for j, f in enumerate(files, 1):
            local_path = os.path.join(local_folder, f["name"])
            # Skip if already downloaded and correct size
            if os.path.exists(local_path) and os.path.getsize(local_path) == f["size"]:
                skipped += 1
                continue
            if j == 1 or j % 100 == 0 or j == len(files):
                log(f"  [{j}/{len(files)}] {f['name']}")
            try:
                nas_download_file(f["path"], local_path, sid)
                downloaded += 1
            except Exception as e:
                log(f"  ERROR {f['name']}: {e}")

        # Download sequencing_summary
        os.makedirs(local_sample, exist_ok=True)
        run_files = nas_list_all_files(nas_run, sid)
        for f in run_files:
            if "sequencing_summary" in f["name"]:
                local_path = os.path.join(local_sample, f["name"])
                if not (os.path.exists(local_path) and os.path.getsize(local_path) == f["size"]):
                    log(f"  Summary: {f['name']}")
                    nas_download_file(f["path"], local_path, sid)

        log(f"  Done: {downloaded} new, {skipped} cached")
        print(flush=True)

    log("ALL DOWNLOADS COMPLETE")
    do_status()


def do_status():
    """Show what's downloaded locally."""
    print(flush=True)
    print("=" * 60, flush=True)
    print(f" LOCAL DATA STATUS: {LOCAL_DATA}", flush=True)
    print("=" * 60, flush=True)

    total_size = 0
    total_files = 0
    for i, (nas_dir, run_subdir, data_type) in enumerate(ALL_SAMPLES):
        local_sample = os.path.join(LOCAL_DATA, nas_dir, run_subdir)
        if data_type == "fastq":
            data_dir = os.path.join(local_sample, "fastq_pass")
            ext = "*.fastq.gz"
        else:
            data_dir = os.path.join(local_sample, "fast5")
            ext = "*.fast5"

        if os.path.isdir(data_dir):
            files = glob.glob(os.path.join(data_dir, ext))
            size = sum(os.path.getsize(f) for f in files)
            total_size += size
            total_files += len(files)
            # Check summary
            summaries = glob.glob(os.path.join(local_sample, "sequencing_summary*"))
            status = "OK" if len(files) > 0 and len(summaries) > 0 else "INCOMPLETE"
            print(f"  [{status:>10}] {nas_dir}: {len(files)} files ({size/(1024**3):.2f} GB) + {len(summaries)} summary", flush=True)
        else:
            print(f"  [   MISSING] {nas_dir}", flush=True)

    print(f"\n  TOTAL: {total_files} files, {total_size/(1024**3):.2f} GB", flush=True)


def do_upload():
    """Step 2: Upload all local data to server via tar pipe."""
    print("=" * 60, flush=True)
    print(f" STEP 2: Upload to server via tar pipe", flush=True)
    print("=" * 60, flush=True)

    if not os.path.isdir(LOCAL_DATA):
        log("ERROR: No local data found. Run 'download' first.")
        return

    # Create destination on server
    subprocess.run(
        ["ssh", "-o", "ConnectTimeout=30", SERVER, f"mkdir -p '{SERVER_BASE}'"],
        capture_output=True, timeout=60
    )

    # Upload each sample individually (allows resume if interrupted)
    for i, (nas_dir, run_subdir, data_type) in enumerate(ALL_SAMPLES):
        local_nas_dir = os.path.join(LOCAL_DATA, nas_dir)
        if not os.path.isdir(local_nas_dir):
            log(f"[{i+1}/14] SKIP {nas_dir} (not downloaded)")
            continue

        # Check if already on server
        result = subprocess.run(
            ["ssh", "-o", "ConnectTimeout=30", SERVER,
             f"test -d '{SERVER_BASE}/{nas_dir}' && du -sh '{SERVER_BASE}/{nas_dir}' 2>/dev/null | cut -f1"],
            capture_output=True, text=True, timeout=30
        )
        server_size = result.stdout.strip()
        if server_size and not server_size.startswith("0"):
            # Rough check: compare number of files
            if data_type == "fastq":
                local_count = len(glob.glob(os.path.join(local_nas_dir, "**", "*.fastq.gz"), recursive=True))
            else:
                local_count = len(glob.glob(os.path.join(local_nas_dir, "**", "*.fast5"), recursive=True))

            result2 = subprocess.run(
                ["ssh", "-o", "ConnectTimeout=30", SERVER,
                 f"find '{SERVER_BASE}/{nas_dir}' -name '*.fastq.gz' -o -name '*.fast5' 2>/dev/null | wc -l"],
                capture_output=True, text=True, timeout=30
            )
            server_count = int(result2.stdout.strip()) if result2.stdout.strip().isdigit() else 0

            if server_count >= local_count and local_count > 0:
                log(f"[{i+1}/14] SKIP {nas_dir} (already on server: {server_size}, {server_count} files)")
                continue

        log(f"[{i+1}/14] Uploading {nas_dir}...")

        # Remove partial upload if exists
        subprocess.run(
            ["ssh", "-o", "ConnectTimeout=30", SERVER, f"rm -rf '{SERVER_BASE}/{nas_dir}'"],
            capture_output=True, timeout=60
        )

        # Tar pipe upload
        cmd = f'tar cf - -C "{LOCAL_DATA}" "{nas_dir}" | ssh -o ConnectTimeout=30 {SERVER} "tar xf - -C \'{SERVER_BASE}\'"'
        t0 = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=7200)
        elapsed = time.time() - t0

        if result.returncode == 0:
            log(f"  Done ({elapsed:.0f}s)")
        else:
            log(f"  ERROR: {result.stderr[:200]}")
            log(f"  Re-run: python scripts/transfer_v4_twostep.py upload")
            return

    log("ALL UPLOADS COMPLETE")

    # Verify
    log("Verifying on server...")
    result = subprocess.run(
        ["ssh", "-o", "ConnectTimeout=30", SERVER,
         f"du -sh {SERVER_BASE}/*/ 2>/dev/null; echo '---'; df -h /"],
        capture_output=True, text=True, timeout=30
    )
    print(result.stdout, flush=True)


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python scripts/transfer_v4_twostep.py download [start_index]")
        print("  python scripts/transfer_v4_twostep.py upload")
        print("  python scripts/transfer_v4_twostep.py status")
        return

    cmd = sys.argv[1]
    if cmd == "download":
        start = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        do_download(start)
    elif cmd == "upload":
        do_upload()
    elif cmd == "status":
        do_status()
    else:
        print(f"Unknown command: {cmd}")


if __name__ == "__main__":
    main()
