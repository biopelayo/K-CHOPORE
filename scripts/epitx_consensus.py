#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — Epitranscriptomic Consensus
=============================================
Standardizes outputs from ELIGOS2, m6Anet, and Nanocompore into a
common format, generates consensus modification calls (sites called
by ≥N tools), and stratifies by biotype (mRNA vs lncRNA).
"""

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Epitranscriptomic consensus analysis")
    parser.add_argument("--eligos-dir", required=True, help="ELIGOS2 output directory")
    parser.add_argument("--m6anet-dir", default="", help="m6Anet output directory")
    parser.add_argument("--nanocompore-dir", default="", help="Nanocompore output directory")
    parser.add_argument("--flair-gtfs", nargs="+", help="FLAIR isoform GTF files")
    parser.add_argument("--lncrna-gtf", default="", help="lncRNA final GTF (optional)")
    parser.add_argument("--min-tools", type=int, default=2, help="Min tools for consensus")
    parser.add_argument("--pval", type=float, default=0.05, help="P-value threshold")
    parser.add_argument("--out-consensus", required=True, help="Output consensus TSV")
    parser.add_argument("--out-biotype", required=True, help="Output biotype stratification TSV")
    return parser.parse_args()


def parse_eligos(eligos_dir, pval_thresh):
    """Parse ELIGOS2 outputs. Returns list of modification site dicts."""
    sites = []
    if not os.path.isdir(eligos_dir):
        return sites

    for filepath in glob.glob(os.path.join(eligos_dir, "*_eligos_output.txt")):
        sample = os.path.basename(filepath).replace("_eligos_output.txt", "")
        with open(filepath) as f:
            for line in f:
                if line.startswith("#") or line.startswith("chrom"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                try:
                    pval = float(parts[5]) if parts[5] != "NA" else 1.0
                except (ValueError, IndexError):
                    pval = 1.0

                if pval < pval_thresh:
                    sites.append({
                        "chr": parts[0],
                        "pos": int(parts[1]),
                        "strand": parts[3] if len(parts) > 3 else ".",
                        "mod_type": "unknown",
                        "score": pval,
                        "tool": "ELIGOS2",
                        "sample": sample
                    })
    return sites


def parse_m6anet(m6anet_dir, pval_thresh):
    """Parse m6Anet site probability outputs."""
    sites = []
    if not m6anet_dir or not os.path.isdir(m6anet_dir):
        return sites

    for filepath in glob.glob(os.path.join(m6anet_dir, "*/data.site_proba.csv")):
        sample = os.path.basename(os.path.dirname(filepath))
        with open(filepath) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    prob = float(row.get("probability_modified", 0))
                except (ValueError, TypeError):
                    prob = 0

                if prob > (1 - pval_thresh):  # m6Anet uses probability, not p-value
                    sites.append({
                        "chr": row.get("transcript_id", "").split("|")[0] if "|" in row.get("transcript_id", "") else "",
                        "pos": int(row.get("transcript_position", 0)),
                        "strand": ".",
                        "mod_type": "m6A",
                        "score": prob,
                        "tool": "m6Anet",
                        "sample": sample
                    })
    return sites


def parse_nanocompore(nanocompore_dir, pval_thresh):
    """Parse Nanocompore results."""
    sites = []
    if not nanocompore_dir or not os.path.isdir(nanocompore_dir):
        return sites

    for filepath in glob.glob(os.path.join(nanocompore_dir, "*_results.tsv")):
        comparison = os.path.basename(filepath).replace("_results.tsv", "")
        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                try:
                    pval = float(row.get("GMM_logit_pvalue", row.get("KS_dwell_pvalue", 1)))
                except (ValueError, TypeError):
                    pval = 1.0

                if pval < pval_thresh:
                    sites.append({
                        "chr": row.get("ref", row.get("chr", "")),
                        "pos": int(row.get("pos", 0)),
                        "strand": row.get("strand", "."),
                        "mod_type": "differential",
                        "score": pval,
                        "tool": "Nanocompore",
                        "sample": comparison
                    })
    return sites


def parse_lncrna_regions(lncrna_gtf):
    """Parse lncRNA GTF to get genomic regions for biotype classification."""
    regions = []
    if not lncrna_gtf or not os.path.isfile(lncrna_gtf):
        return regions
    import re
    with open(lncrna_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            try:
                regions.append({
                    "chr": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6] if len(parts) > 6 else "."
                })
            except (ValueError, IndexError):
                pass
    return regions


def site_in_lncrna(site, lncrna_regions):
    """Check if a modification site falls within a lncRNA region."""
    for region in lncrna_regions:
        if (site["chr"] == region["chr"] and
            region["start"] <= site["pos"] <= region["end"]):
            return True
    return False


def main():
    args = parse_args()

    print("[K-CHOPORE] Module 4: Epitranscriptomic Consensus Analysis")
    print(f"[K-CHOPORE] Minimum tools for consensus: {args.min_tools}")
    print(f"[K-CHOPORE] P-value threshold: {args.pval}")

    # Parse all tool outputs
    eligos_sites = parse_eligos(args.eligos_dir, args.pval)
    print(f"[K-CHOPORE] ELIGOS2 significant sites: {len(eligos_sites)}")

    m6anet_sites = parse_m6anet(args.m6anet_dir, args.pval)
    print(f"[K-CHOPORE] m6Anet significant sites: {len(m6anet_sites)}")

    nanocompore_sites = parse_nanocompore(args.nanocompore_dir, args.pval)
    print(f"[K-CHOPORE] Nanocompore significant sites: {len(nanocompore_sites)}")

    all_sites = eligos_sites + m6anet_sites + nanocompore_sites
    print(f"[K-CHOPORE] Total significant sites across tools: {len(all_sites)}")

    # Build site index: (chr, pos, strand) -> list of tools
    site_tools = defaultdict(lambda: {"tools": set(), "mod_types": set(), "best_score": 1.0, "samples": set()})

    for site in all_sites:
        key = (site["chr"], site["pos"], site["strand"])
        site_tools[key]["tools"].add(site["tool"])
        site_tools[key]["mod_types"].add(site["mod_type"])
        site_tools[key]["samples"].add(site["sample"])
        if site["score"] < site_tools[key]["best_score"]:
            site_tools[key]["best_score"] = site["score"]

    # Apply consensus filter
    consensus_sites = {}
    for key, info in site_tools.items():
        if len(info["tools"]) >= args.min_tools:
            consensus_sites[key] = info

    print(f"\n[K-CHOPORE] === Consensus Results ===")
    print(f"[K-CHOPORE] Sites called by >={args.min_tools} tools: {len(consensus_sites)}")

    # Tool overlap statistics
    tool_counts = defaultdict(int)
    for info in site_tools.values():
        for tool in info["tools"]:
            tool_counts[tool] += 1
    for tool, count in sorted(tool_counts.items()):
        print(f"[K-CHOPORE]   {tool}: {count} sites")

    # Parse lncRNA regions for biotype classification
    lncrna_regions = parse_lncrna_regions(args.lncrna_gtf)

    # Write consensus output
    with open(args.out_consensus, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "chr", "pos", "strand", "mod_types", "n_tools",
            "tools", "best_pvalue", "samples"
        ])
        for (chrom, pos, strand), info in sorted(consensus_sites.items()):
            writer.writerow([
                chrom, pos, strand,
                ",".join(info["mod_types"]),
                len(info["tools"]),
                ",".join(sorted(info["tools"])),
                f"{info['best_score']:.2e}",
                ",".join(sorted(info["samples"]))
            ])

    # Write biotype stratification (all sites, not just consensus)
    biotype_stats = {"mRNA": 0, "lncRNA": 0, "unknown": 0}

    with open(args.out_biotype, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "chr", "pos", "strand", "biotype", "n_tools", "tools", "best_pvalue"
        ])
        for (chrom, pos, strand), info in sorted(site_tools.items()):
            site_dict = {"chr": chrom, "pos": pos, "strand": strand}
            if lncrna_regions and site_in_lncrna(site_dict, lncrna_regions):
                biotype = "lncRNA"
            elif chrom:
                biotype = "mRNA"
            else:
                biotype = "unknown"

            biotype_stats[biotype] += 1
            writer.writerow([
                chrom, pos, strand, biotype,
                len(info["tools"]),
                ",".join(sorted(info["tools"])),
                f"{info['best_score']:.2e}"
            ])

    print(f"\n[K-CHOPORE] Biotype distribution:")
    for bt, count in biotype_stats.items():
        print(f"[K-CHOPORE]   {bt}: {count} sites")

    print(f"\n[K-CHOPORE] Output:")
    print(f"[K-CHOPORE]   Consensus: {args.out_consensus}")
    print(f"[K-CHOPORE]   Biotype: {args.out_biotype}")
    print("[K-CHOPORE] Epitranscriptomic consensus completed.")


if __name__ == "__main__":
    main()
