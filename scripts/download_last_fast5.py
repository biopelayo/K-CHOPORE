#!/usr/bin/env python3
"""Download the remaining FAST5 files for anac017_AA_R2-2 to D: drive."""
import os
import sys
import json
import urllib.request
import urllib.parse
import ssl
import time

NAS_HOST = "valmeilab.synology.me"
NAS_PORT = "5001"
NAS_USER = "pelayo"
NAS_PASS = "Pelamovic39@"
NAS_BASE = "/HTData_and_DBs/NGS/Nanopore/2022_Arabidopsis_AA_anac017_DRS"

LOCAL_DATA = os.path.join("D:", os.sep, "kchopore_nas_data")

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

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

def main():
    nas_dir = "anac017_AA_R2-2"
    run_subdir = "no_sample/20221215_1232_MC-112869_FAT21048_d754e27f"

    log(f"Downloading remaining FAST5 for {nas_dir} to D:")
    log("Logging in to NAS...")
    sid = nas_login()
    log("OK")

    nas_run = f"{NAS_BASE}/{nas_dir}/{run_subdir}"
    local_sample = os.path.join(LOCAL_DATA, nas_dir, run_subdir)

    # Download fast5 folder
    folder = f"{nas_run}/fast5"
    local_folder = os.path.join(local_sample, "fast5")
    os.makedirs(local_folder, exist_ok=True)

    files = nas_list_all_files(folder, sid)
    log(f"  {len(files)} files total on NAS")

    downloaded = 0
    skipped = 0
    for j, f in enumerate(files, 1):
        local_path = os.path.join(local_folder, f["name"])
        if os.path.exists(local_path) and os.path.getsize(local_path) == f["size"]:
            skipped += 1
            continue
        if j == 1 or j % 10 == 0 or j == len(files):
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
    log("ALL DONE")

if __name__ == "__main__":
    main()
