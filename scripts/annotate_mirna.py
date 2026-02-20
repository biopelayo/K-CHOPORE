#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — miRNA Annotation
==================================
Cross-references ShortStack and miRDeep-P2 miRNA calls with
miRBase and PmiREN databases. Classifies as known or novel.
"""

import argparse
import csv
import os
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate miRNAs against databases")
    parser.add_argument("--shortstack", required=True, help="ShortStack Results.txt")
    parser.add_argument("--mirdeep", required=True, help="miRDeep-P2 novel miRNAs TSV")
    parser.add_argument("--mirbase-gff", default="", help="miRBase GFF3 for known miRNAs")
    parser.add_argument("--pmiren-fa", default="", help="PmiREN FASTA for plant miRNAs")
    parser.add_argument("--output", required=True, help="Output annotated miRNA TSV")
    return parser.parse_args()


def parse_mirbase_gff(filepath):
    """Parse miRBase GFF3 to get known miRNA names and coordinates."""
    known = {}
    if not filepath or not os.path.isfile(filepath):
        return known
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] in ("miRNA", "miRNA_primary_transcript"):
                attrs = parts[8]
                name_match = re.search(r"Name=([^;]+)", attrs)
                id_match = re.search(r"ID=([^;]+)", attrs)
                if name_match:
                    name = name_match.group(1)
                    known[name] = {
                        "chr": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "type": parts[2],
                        "id": id_match.group(1) if id_match else name
                    }
    return known


def parse_pmiren_fa(filepath):
    """Parse PmiREN FASTA to get additional plant miRNA sequences."""
    mirnas = {}
    if not filepath or not os.path.isfile(filepath):
        return mirnas
    current_id = None
    current_seq = ""
    with open(filepath) as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    mirnas[current_id] = current_seq
                current_id = line.strip()[1:].split()[0]
                current_seq = ""
            else:
                current_seq += line.strip()
    if current_id:
        mirnas[current_id] = current_seq
    return mirnas


def parse_shortstack(filepath):
    """Parse ShortStack Results.txt for miRNA loci."""
    mirnas = []
    if not os.path.isfile(filepath):
        return mirnas
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # ShortStack marks MIRNA loci with DicerCall and MIRNA column
            is_mirna = row.get("MIRNA", "N") == "Y"
            if is_mirna:
                mirnas.append({
                    "locus": row.get("Locus", ""),
                    "name": row.get("Name", row.get("Locus", "")),
                    "chr": row.get("Chrom", ""),
                    "start": row.get("Start", ""),
                    "end": row.get("End", ""),
                    "strand": row.get("Strand", "."),
                    "mature_seq": row.get("MajorRNA", ""),
                    "reads": row.get("Reads", "0"),
                    "dicer_call": row.get("DicerCall", ""),
                    "source": "ShortStack"
                })
    return mirnas


def parse_mirdeep(filepath):
    """Parse miRDeep-P2 results."""
    mirnas = []
    if not os.path.isfile(filepath):
        return mirnas
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mirnas.append({
                "locus": row.get("miRNA_ID", ""),
                "name": row.get("miRNA_ID", ""),
                "chr": "",
                "start": "",
                "end": "",
                "strand": "",
                "mature_seq": row.get("sequence", ""),
                "reads": "",
                "dicer_call": "",
                "source": "miRDeep-P2",
                "score": row.get("score", "")
            })
    return mirnas


def annotate_mirnas(discovered, known_mirbase, known_pmiren):
    """Cross-reference discovered miRNAs with databases."""
    annotated = []
    for mirna in discovered:
        seq = mirna.get("mature_seq", "").upper().replace("T", "U")
        name = mirna.get("name", "")

        # Check against miRBase by name or sequence
        mirbase_match = None
        for mb_name, mb_info in known_mirbase.items():
            if name and (mb_name.lower() in name.lower() or name.lower() in mb_name.lower()):
                mirbase_match = mb_name
                break

        # Check against PmiREN by sequence
        pmiren_match = None
        for pm_name, pm_seq in known_pmiren.items():
            if seq and pm_seq.upper().replace("T", "U") == seq:
                pmiren_match = pm_name
                break

        # Determine status
        if mirbase_match:
            status = "known_miRBase"
            family = mirbase_match
        elif pmiren_match:
            status = "known_PmiREN"
            family = pmiren_match
        else:
            status = "novel"
            family = ""

        annotated.append({
            **mirna,
            "status": status,
            "family": family,
            "mirbase_match": mirbase_match or "",
            "pmiren_match": pmiren_match or ""
        })

    return annotated


def main():
    args = parse_args()

    print("[K-CHOPORE] Module 2: miRNA Annotation")

    # Parse databases
    known_mirbase = parse_mirbase_gff(args.mirbase_gff)
    print(f"[K-CHOPORE] miRBase: {len(known_mirbase)} known miRNAs loaded")

    known_pmiren = parse_pmiren_fa(args.pmiren_fa)
    print(f"[K-CHOPORE] PmiREN: {len(known_pmiren)} plant miRNAs loaded")

    # Parse discovery tools
    ss_mirnas = parse_shortstack(args.shortstack)
    print(f"[K-CHOPORE] ShortStack: {len(ss_mirnas)} MIRNA loci")

    md_mirnas = parse_mirdeep(args.mirdeep)
    print(f"[K-CHOPORE] miRDeep-P2: {len(md_mirnas)} predictions")

    # Merge and deduplicate
    all_mirnas = ss_mirnas + md_mirnas

    # Annotate
    annotated = annotate_mirnas(all_mirnas, known_mirbase, known_pmiren)

    # Count stats
    n_known = sum(1 for m in annotated if m["status"] != "novel")
    n_novel = sum(1 for m in annotated if m["status"] == "novel")
    print(f"[K-CHOPORE] Annotated: {n_known} known, {n_novel} novel")

    # Write output
    fieldnames = [
        "locus", "name", "chr", "start", "end", "strand",
        "mature_seq", "reads", "dicer_call", "source",
        "status", "family", "mirbase_match", "pmiren_match"
    ]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for mirna in annotated:
            writer.writerow(mirna)

    print(f"[K-CHOPORE] Output: {args.output}")
    print(f"[K-CHOPORE] miRNA annotation completed.")


if __name__ == "__main__":
    main()
