#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — miRNA-Target Integration
==========================================
Merges TargetFinder + psRNATarget predictions, validates with
degradome data, and cross-references with differential expression
to find anti-correlated miRNA-target pairs.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="Integrate miRNA-target evidence")
    parser.add_argument("--targetfinder", required=True, help="TargetFinder results TSV")
    parser.add_argument("--psrnatarget", required=True, help="psRNATarget results TSV")
    parser.add_argument("--degradome", required=True, help="Degradome validated TSV")
    parser.add_argument("--de-genotype", required=True, help="DESeq2 genotype results CSV")
    parser.add_argument("--de-treatment", required=True, help="DESeq2 treatment results CSV")
    parser.add_argument("--lncrna-de-genotype", default="", help="lncRNA DESeq2 genotype CSV")
    parser.add_argument("--lncrna-de-treatment", default="", help="lncRNA DESeq2 treatment CSV")
    parser.add_argument("--out-evidence", required=True, help="Output evidence table TSV")
    parser.add_argument("--out-network", required=True, help="Output network TSV")
    return parser.parse_args()


def load_tsv(filepath):
    """Load a TSV file as list of dicts."""
    rows = []
    if not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        return rows
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def load_csv(filepath):
    """Load a CSV file as list of dicts."""
    rows = []
    if not filepath or not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        return rows
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def parse_de_results(filepath):
    """Parse DESeq2 results into dict of gene -> {lfc, padj, significant}."""
    de = {}
    rows = load_csv(filepath)
    for row in rows:
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


def main():
    args = parse_args()

    print("[K-CHOPORE] Module 3: miRNA-Target Integration")

    # Load predictions
    tf_results = load_tsv(args.targetfinder)
    ps_results = load_tsv(args.psrnatarget)
    deg_results = load_tsv(args.degradome)

    print(f"[K-CHOPORE] TargetFinder predictions: {len(tf_results)}")
    print(f"[K-CHOPORE] psRNATarget predictions: {len(ps_results)}")
    print(f"[K-CHOPORE] Degradome validated: {len(deg_results)}")

    # Load DE results
    de_genotype = parse_de_results(args.de_genotype)
    de_treatment = parse_de_results(args.de_treatment)
    de_lnc_genotype = parse_de_results(args.lncrna_de_genotype) if args.lncrna_de_genotype else {}
    de_lnc_treatment = parse_de_results(args.lncrna_de_treatment) if args.lncrna_de_treatment else {}

    # Merge all DE
    all_de = {**de_genotype, **de_treatment, **de_lnc_genotype, **de_lnc_treatment}

    # Build unified target evidence table
    # Key: (miRNA, target) -> evidence dict
    evidence = defaultdict(lambda: {
        "targetfinder": False, "tf_score": "",
        "psrnatarget": False, "ps_expectation": "",
        "degradome": False, "deg_category": "",
        "target_lfc_genotype": "", "target_padj_genotype": "",
        "target_lfc_treatment": "", "target_padj_treatment": "",
        "is_lncrna_target": False,
        "evidence_score": 0
    })

    # TargetFinder
    for row in tf_results:
        mirna = row.get("miRNA", row.get("query", ""))
        target = row.get("Target", row.get("target", ""))
        score = row.get("Score", row.get("score", ""))
        if mirna and target:
            key = (mirna, target)
            evidence[key]["targetfinder"] = True
            evidence[key]["tf_score"] = score
            evidence[key]["evidence_score"] += 1

    # psRNATarget
    for row in ps_results:
        mirna = row.get("miRNA", row.get("query", ""))
        target = row.get("Target", row.get("target", ""))
        exp = row.get("Expectation", row.get("expectation", ""))
        if mirna and target:
            key = (mirna, target)
            evidence[key]["psrnatarget"] = True
            evidence[key]["ps_expectation"] = exp
            evidence[key]["evidence_score"] += 1

    # Degradome
    for row in deg_results:
        mirna = row.get("miRNA", row.get("query", ""))
        target = row.get("Target", row.get("target", ""))
        cat = row.get("Category", row.get("category", ""))
        if mirna and target:
            key = (mirna, target)
            evidence[key]["degradome"] = True
            evidence[key]["deg_category"] = cat
            evidence[key]["evidence_score"] += 2  # Degradome counts double

    # Cross with DE
    for (mirna, target), ev in evidence.items():
        if target in de_genotype:
            ev["target_lfc_genotype"] = de_genotype[target]["lfc"]
            ev["target_padj_genotype"] = de_genotype[target]["padj"]
        if target in de_treatment:
            ev["target_lfc_treatment"] = de_treatment[target]["lfc"]
            ev["target_padj_treatment"] = de_treatment[target]["padj"]
        if target in de_lnc_genotype or target in de_lnc_treatment:
            ev["is_lncrna_target"] = True

    print(f"[K-CHOPORE] Unique miRNA-target pairs: {len(evidence)}")

    # Count evidence levels
    n_predicted_only = sum(1 for e in evidence.values() if not e["degradome"])
    n_degradome_val = sum(1 for e in evidence.values() if e["degradome"])
    n_both_tools = sum(1 for e in evidence.values() if e["targetfinder"] and e["psrnatarget"])
    n_lncrna = sum(1 for e in evidence.values() if e["is_lncrna_target"])

    print(f"[K-CHOPORE] Predicted only: {n_predicted_only}")
    print(f"[K-CHOPORE] Degradome validated: {n_degradome_val}")
    print(f"[K-CHOPORE] Both tools agree: {n_both_tools}")
    print(f"[K-CHOPORE] lncRNA targets: {n_lncrna}")

    # Write evidence table
    with open(args.out_evidence, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "miRNA", "target", "evidence_score",
            "targetfinder", "tf_score", "psrnatarget", "ps_expectation",
            "degradome", "deg_category",
            "target_lfc_genotype", "target_padj_genotype",
            "target_lfc_treatment", "target_padj_treatment",
            "is_lncrna_target"
        ])
        for (mirna, target), ev in sorted(evidence.items(), key=lambda x: -x[1]["evidence_score"]):
            writer.writerow([
                mirna, target, ev["evidence_score"],
                ev["targetfinder"], ev["tf_score"],
                ev["psrnatarget"], ev["ps_expectation"],
                ev["degradome"], ev["deg_category"],
                ev["target_lfc_genotype"], ev["target_padj_genotype"],
                ev["target_lfc_treatment"], ev["target_padj_treatment"],
                ev["is_lncrna_target"]
            ])

    # Write network (for visualization tools like Cytoscape)
    with open(args.out_network, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["source", "target", "interaction", "weight"])
        for (mirna, target), ev in evidence.items():
            interaction = "validated" if ev["degradome"] else "predicted"
            writer.writerow([mirna, target, interaction, ev["evidence_score"]])

    print(f"[K-CHOPORE] Evidence table: {args.out_evidence}")
    print(f"[K-CHOPORE] Network table: {args.out_network}")
    print("[K-CHOPORE] Target integration completed.")


if __name__ == "__main__":
    main()
