#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — miRNA Target Prediction Wrapper
=================================================
Wraps TargetFinder for plant miRNA target prediction.
Reads mature miRNA sequences from the annotated miRNA table
and predicts targets against the transcriptome.
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description="Predict miRNA targets")
    parser.add_argument("--mirnas", required=True, help="Annotated miRNA TSV")
    parser.add_argument("--transcriptome", required=True, help="Transcriptome FASTA")
    parser.add_argument("--cutoff", type=float, default=4.0, help="TargetFinder score cutoff")
    parser.add_argument("--output", required=True, help="Output targets TSV")
    return parser.parse_args()


def load_mirna_sequences(filepath):
    """Load mature miRNA sequences from annotated TSV."""
    mirnas = {}
    if not os.path.isfile(filepath):
        return mirnas
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row.get("name", row.get("locus", ""))
            seq = row.get("mature_seq", "")
            if name and seq:
                # Convert DNA to RNA if needed
                seq = seq.upper().replace("T", "U")
                mirnas[name] = seq
    return mirnas


def run_targetfinder(mirna_name, mirna_seq, transcriptome, cutoff):
    """Run TargetFinder for a single miRNA. Returns list of targets."""
    targets = []

    # Check if TargetFinder is available
    targetfinder_cmd = None
    for cmd in ["targetfinder.pl", "TargetFinder", "/opt/targetfinder/targetfinder.pl"]:
        try:
            subprocess.run([cmd, "--help"], capture_output=True, timeout=5)
            targetfinder_cmd = cmd
            break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue

    if not targetfinder_cmd:
        return targets

    try:
        result = subprocess.run(
            [targetfinder_cmd, "-s", mirna_seq, "-d", transcriptome, "-c", str(cutoff)],
            capture_output=True, text=True, timeout=300
        )
        # Parse TargetFinder output
        for line in result.stdout.strip().split("\n"):
            if line and not line.startswith("#"):
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    targets.append({
                        "miRNA": mirna_name,
                        "target": parts[0],
                        "score": parts[1] if len(parts) > 1 else "",
                        "alignment": parts[2] if len(parts) > 2 else ""
                    })
    except (subprocess.TimeoutExpired, Exception) as e:
        print(f"[K-CHOPORE] TargetFinder failed for {mirna_name}: {e}", file=sys.stderr)

    return targets


def main():
    args = parse_args()

    print("[K-CHOPORE] Module 3: miRNA Target Prediction")
    print(f"[K-CHOPORE] Score cutoff: {args.cutoff}")

    # Load miRNA sequences
    mirnas = load_mirna_sequences(args.mirnas)
    print(f"[K-CHOPORE] miRNAs loaded: {len(mirnas)}")

    if not mirnas:
        print("[K-CHOPORE] No miRNA sequences found — creating empty output.")
        with open(args.output, "w") as f:
            f.write("miRNA\ttarget\tscore\talignment\n")
        return

    # Check if transcriptome exists
    if not os.path.isfile(args.transcriptome):
        print(f"[K-CHOPORE] WARNING: Transcriptome not found: {args.transcriptome}")
        with open(args.output, "w") as f:
            f.write("miRNA\ttarget\tscore\talignment\n")
        return

    # Run predictions for each miRNA
    all_targets = []
    for i, (name, seq) in enumerate(mirnas.items()):
        if (i + 1) % 10 == 0 or i == 0:
            print(f"[K-CHOPORE] Processing miRNA {i+1}/{len(mirnas)}: {name}")
        targets = run_targetfinder(name, seq, args.transcriptome, args.cutoff)
        all_targets.extend(targets)

    print(f"[K-CHOPORE] Total predicted targets: {len(all_targets)}")

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["miRNA", "target", "score", "alignment"],
                                delimiter="\t")
        writer.writeheader()
        for target in all_targets:
            writer.writerow(target)

    print(f"[K-CHOPORE] Output: {args.output}")
    print("[K-CHOPORE] Target prediction completed.")


if __name__ == "__main__":
    main()
