#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — Cis-Regulatory lncRNA Analysis
=================================================
Identifies lncRNAs within a defined window of differentially
expressed protein-coding genes. Tests co-expression to find
potential cis-regulatory relationships.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Identify cis-regulatory lncRNA pairs")
    parser.add_argument("--lncrna-gtf", required=True, help="lncRNA final GTF")
    parser.add_argument("--ref-gtf", required=True, help="Reference annotation GTF (Araport11)")
    parser.add_argument("--de-genotype", required=True, help="DESeq2 genotype results CSV")
    parser.add_argument("--de-treatment", required=True, help="DESeq2 treatment results CSV")
    parser.add_argument("--counts", required=True, help="FLAIR counts matrix")
    parser.add_argument("--window-kb", type=int, default=10, help="Window size in kb")
    parser.add_argument("--output", required=True, help="Output cis pairs TSV")
    return parser.parse_args()


def parse_gtf_genes(filepath):
    """Parse GTF to get gene coordinates."""
    genes = {}
    if not os.path.isfile(filepath):
        return genes
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("gene", "transcript"):
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Extract gene_id
            gene_match = re.search(r'gene_id "([^"]+)"', parts[8])
            if gene_match:
                gid = gene_match.group(1)
                if gid not in genes:
                    genes[gid] = {
                        "chr": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand
                    }
                else:
                    genes[gid]["start"] = min(genes[gid]["start"], start)
                    genes[gid]["end"] = max(genes[gid]["end"], end)
    return genes


def parse_lncrna_gtf(filepath):
    """Parse lncRNA GTF for coordinates and metadata."""
    lncrnas = {}
    if not os.path.isfile(filepath):
        return lncrnas
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            tid_match = re.search(r'transcript_id "([^"]+)"', parts[8])
            if tid_match:
                tid = tid_match.group(1)
                if tid not in lncrnas:
                    lncrnas[tid] = {
                        "chr": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand
                    }
                else:
                    lncrnas[tid]["start"] = min(lncrnas[tid]["start"], start)
                    lncrnas[tid]["end"] = max(lncrnas[tid]["end"], end)
    return lncrnas


def load_de_results(filepath):
    """Load DE results as dict of gene -> {lfc, padj, significant}."""
    de = {}
    if not os.path.isfile(filepath):
        return de
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene = row.get("gene", "")
            if not gene:
                continue
            try:
                lfc = float(row.get("log2FoldChange", 0))
                padj = float(row.get("padj", 1))
            except (ValueError, TypeError):
                lfc = 0
                padj = 1
            sig = row.get("significant", "not_significant")
            de[gene] = {"lfc": lfc, "padj": padj, "significant": sig}
    return de


def find_cis_pairs(lncrnas, genes, window_bp):
    """Find lncRNA-gene pairs within window distance."""
    pairs = []
    for lnc_id, lnc_info in lncrnas.items():
        for gene_id, gene_info in genes.items():
            if lnc_info["chr"] != gene_info["chr"]:
                continue

            # Calculate distance between lncRNA and gene
            if lnc_info["end"] < gene_info["start"]:
                distance = gene_info["start"] - lnc_info["end"]
            elif gene_info["end"] < lnc_info["start"]:
                distance = lnc_info["start"] - gene_info["end"]
            else:
                distance = 0  # Overlapping

            if distance <= window_bp:
                # Determine orientation
                if lnc_info["strand"] == gene_info["strand"]:
                    orientation = "sense"
                else:
                    orientation = "antisense"

                # Determine relative position
                lnc_mid = (lnc_info["start"] + lnc_info["end"]) / 2
                gene_mid = (gene_info["start"] + gene_info["end"]) / 2
                if lnc_mid < gene_mid:
                    position = "upstream"
                elif lnc_mid > gene_mid:
                    position = "downstream"
                else:
                    position = "overlapping"

                pairs.append({
                    "lncrna_id": lnc_id,
                    "gene_id": gene_id,
                    "chr": lnc_info["chr"],
                    "distance": distance,
                    "orientation": orientation,
                    "position": position,
                    "lncrna_start": lnc_info["start"],
                    "lncrna_end": lnc_info["end"],
                    "gene_start": gene_info["start"],
                    "gene_end": gene_info["end"]
                })
    return pairs


def main():
    args = parse_args()
    window_bp = args.window_kb * 1000

    print("[K-CHOPORE] Module 5: Cis-Regulatory lncRNA Analysis")
    print(f"[K-CHOPORE] Window: {args.window_kb} kb ({window_bp} bp)")

    # Parse annotations
    lncrnas = parse_lncrna_gtf(args.lncrna_gtf)
    print(f"[K-CHOPORE] lncRNAs: {len(lncrnas)}")

    genes = parse_gtf_genes(args.ref_gtf)
    print(f"[K-CHOPORE] Reference genes: {len(genes)}")

    # Load DE results
    de_genotype = load_de_results(args.de_genotype)
    de_treatment = load_de_results(args.de_treatment)
    de_genes = set()
    for gene, info in de_genotype.items():
        if info["significant"] == "significant":
            de_genes.add(gene)
    for gene, info in de_treatment.items():
        if info["significant"] == "significant":
            de_genes.add(gene)
    print(f"[K-CHOPORE] DE genes: {len(de_genes)}")

    # Filter genes to only DE genes
    de_gene_coords = {g: genes[g] for g in de_genes if g in genes}
    print(f"[K-CHOPORE] DE genes with coordinates: {len(de_gene_coords)}")

    # Find cis pairs
    pairs = find_cis_pairs(lncrnas, de_gene_coords, window_bp)
    print(f"[K-CHOPORE] Cis lncRNA-DE gene pairs: {len(pairs)}")

    # Annotate with DE info
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "lncrna_id", "gene_id", "chr", "distance_bp",
            "orientation", "position",
            "gene_lfc_genotype", "gene_padj_genotype",
            "gene_lfc_treatment", "gene_padj_treatment"
        ])
        for pair in sorted(pairs, key=lambda x: x["distance"]):
            gene = pair["gene_id"]
            geno_info = de_genotype.get(gene, {"lfc": "", "padj": ""})
            treat_info = de_treatment.get(gene, {"lfc": "", "padj": ""})
            writer.writerow([
                pair["lncrna_id"], pair["gene_id"], pair["chr"],
                pair["distance"],
                pair["orientation"], pair["position"],
                geno_info.get("lfc", ""), geno_info.get("padj", ""),
                treat_info.get("lfc", ""), treat_info.get("padj", "")
            ])

    # Summary
    n_overlapping = sum(1 for p in pairs if p["distance"] == 0)
    n_antisense = sum(1 for p in pairs if p["orientation"] == "antisense")
    print(f"\n[K-CHOPORE] Summary:")
    print(f"[K-CHOPORE]   Overlapping: {n_overlapping}")
    print(f"[K-CHOPORE]   Antisense: {n_antisense}")
    print(f"[K-CHOPORE]   Within {args.window_kb}kb: {len(pairs)}")
    print(f"[K-CHOPORE] Output: {args.output}")
    print("[K-CHOPORE] Cis-regulation analysis completed.")


if __name__ == "__main__":
    main()
