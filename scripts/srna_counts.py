#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — Small RNA Count Matrix Generator
===================================================
Generates miRNA count matrix from ShortStack output and
produces sRNA size distribution plots.
"""

import argparse
import csv
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Generate sRNA count matrix")
    parser.add_argument("--shortstack", required=True, help="ShortStack Results.txt")
    parser.add_argument("--annotated", required=True, help="Annotated miRNA TSV")
    parser.add_argument("--out-counts", required=True, help="Output count matrix TSV")
    parser.add_argument("--out-sizedist", required=True, help="Output size distribution PDF")
    return parser.parse_args()


def parse_shortstack_counts(filepath):
    """Extract per-sample counts from ShortStack Results.txt."""
    count_data = {}
    sample_names = []

    if not os.path.isfile(filepath):
        return count_data, sample_names

    with open(filepath) as f:
        header = f.readline().strip().split("\t")

        # ShortStack includes per-sample read counts as columns
        # Identify sample columns (they contain count data)
        name_col = None
        mirna_col = None
        count_cols = []

        for i, col in enumerate(header):
            if col in ("Name", "Locus"):
                name_col = i
            elif col == "MIRNA":
                mirna_col = i
            elif col.endswith(".bam") or col.endswith("_trimmed"):
                count_cols.append((i, col))
                sample_names.append(col.replace(".bam", "").replace("_trimmed", ""))

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(i for i, _ in count_cols) if count_cols else 0:
                continue

            locus_name = parts[name_col] if name_col is not None else parts[0]

            # Only include MIRNA loci
            if mirna_col is not None and parts[mirna_col] != "Y":
                continue

            counts = {}
            for col_idx, col_name in count_cols:
                try:
                    counts[col_name] = int(float(parts[col_idx]))
                except (ValueError, IndexError):
                    counts[col_name] = 0

            count_data[locus_name] = counts

    return count_data, sample_names


def generate_size_distribution(filepath, output_pdf):
    """Generate sRNA size distribution plot."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        sizes = {}
        if os.path.isfile(filepath):
            with open(filepath) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    dicer = row.get("DicerCall", row.get("dicer_call", ""))
                    try:
                        size = int(dicer)
                        sizes[size] = sizes.get(size, 0) + int(float(row.get("Reads", row.get("reads", 1))))
                    except (ValueError, TypeError):
                        pass

        if sizes:
            fig, ax = plt.subplots(figsize=(8, 5))
            x = sorted(sizes.keys())
            y = [sizes[s] for s in x]
            ax.bar(x, y, color="#2196F3", edgecolor="white")
            ax.set_xlabel("sRNA size (nt)")
            ax.set_ylabel("Read count")
            ax.set_title("K-CHOPORE: Small RNA Size Distribution")
            ax.set_xticks(range(min(x), max(x) + 1))
            plt.tight_layout()
            plt.savefig(output_pdf, dpi=150)
            plt.close()
            print(f"[K-CHOPORE] Size distribution plot saved: {output_pdf}")
        else:
            # Create minimal PDF placeholder
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No sRNA size data available",
                    ha="center", va="center", transform=ax.transAxes)
            plt.savefig(output_pdf)
            plt.close()

    except ImportError:
        print("[K-CHOPORE] matplotlib not available — skipping size distribution plot")
        # Create empty file
        with open(output_pdf, "w") as f:
            f.write("# Size distribution plot requires matplotlib\n")


def main():
    args = parse_args()

    print("[K-CHOPORE] Module 2: Generating sRNA count matrix")

    # Parse ShortStack counts
    count_data, sample_names = parse_shortstack_counts(args.shortstack)
    print(f"[K-CHOPORE] MIRNA loci with counts: {len(count_data)}")
    print(f"[K-CHOPORE] Samples detected: {len(sample_names)}")

    # Write count matrix
    with open(args.out_counts, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["miRNA_locus"] + sample_names)
        for locus, counts in sorted(count_data.items()):
            row = [locus] + [counts.get(s, 0) for s in sample_names]
            writer.writerow(row)

    print(f"[K-CHOPORE] Count matrix: {args.out_counts}")

    # Generate size distribution
    generate_size_distribution(args.shortstack, args.out_sizedist)

    print("[K-CHOPORE] sRNA count matrix generation completed.")


if __name__ == "__main__":
    main()
