#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — lncRNA Consensus Classification
=================================================
Integrates results from TransDecoder, FEELnc, CPC2, and CPAT to
produce a consensus set of high-confidence lncRNAs.

Classification categories:
  - lincRNA: intergenic long non-coding RNA
  - antisense: antisense to protein-coding gene
  - intronic: within intron of protein-coding gene
  - TE-derived: overlapping transposable element
  - known: matches CANTATAdb entry

Consensus rule: transcript must be called "non-coding" by ≥N of 3
tools (FEELnc, CPC2, CPAT) AND must NOT have ORF >100aa by TransDecoder.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description="lncRNA consensus classification")
    parser.add_argument("--transdecoder", required=True, help="TransDecoder longest_orfs.pep")
    parser.add_argument("--feelnc-gtf", required=True, help="FEELnc lncRNA GTF")
    parser.add_argument("--feelnc-class", required=True, help="FEELnc classification file")
    parser.add_argument("--cpc2", required=True, help="CPC2 results file")
    parser.add_argument("--cpat", required=True, help="CPAT results TSV")
    parser.add_argument("--candidates-gtf", required=True, help="Candidates merged GTF")
    parser.add_argument("--counts-matrix", required=True, help="FLAIR counts_matrix.tsv")
    parser.add_argument("--cantatadb", default="", help="CANTATAdb BED file (optional)")
    parser.add_argument("--te-bed", default="", help="TE annotation BED (optional)")
    parser.add_argument("--min-consensus", type=int, default=2, help="Min tools agreeing non-coding")
    parser.add_argument("--max-orf-aa", type=int, default=100, help="TransDecoder ORF cutoff (aa)")
    parser.add_argument("--cpat-threshold", type=float, default=0.39, help="CPAT coding probability cutoff")
    parser.add_argument("--out-gtf", required=True, help="Output final lncRNA GTF")
    parser.add_argument("--out-bed", required=True, help="Output final lncRNA BED")
    parser.add_argument("--out-counts", required=True, help="Output lncRNA counts matrix")
    parser.add_argument("--out-report", required=True, help="Output summary report TSV")
    return parser.parse_args()


def parse_transdecoder(filepath, max_orf_aa):
    """Return set of transcript IDs with ORF > max_orf_aa."""
    coding_ids = set()
    if not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        return coding_ids
    current_id = None
    current_len = 0
    with open(filepath) as f:
        for line in f:
            if line.startswith(">"):
                # Save previous
                if current_id and current_len > max_orf_aa:
                    coding_ids.add(current_id)
                # Parse new header: >transcript_id.p1 ...
                header = line.strip().split()[0][1:]
                # Remove .p suffix
                current_id = re.sub(r"\.p\d+$", "", header)
                current_len = 0
            else:
                current_len += len(line.strip())
    # Last entry
    if current_id and current_len > max_orf_aa:
        coding_ids.add(current_id)
    return coding_ids


def parse_feelnc_ids(gtf_path):
    """Return set of transcript IDs called lncRNA by FEELnc."""
    ids = set()
    if not os.path.isfile(gtf_path) or os.path.getsize(gtf_path) == 0:
        return ids
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            # Extract transcript_id from GTF attributes
            match = re.search(r'transcript_id "([^"]+)"', line)
            if match:
                ids.add(match.group(1))
    return ids


def parse_feelnc_classification(filepath):
    """Return dict of transcript_id -> classification type."""
    classifications = {}
    if not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        return classifications
    with open(filepath) as f:
        for line in f:
            if line.startswith("#") or line.startswith("isBest"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                tid = parts[1] if len(parts) > 1 else parts[0]
                lnc_type = parts[2] if len(parts) > 2 else "unknown"
                classifications[tid] = lnc_type
    return classifications


def parse_cpc2(filepath):
    """Return set of transcript IDs called noncoding by CPC2."""
    noncoding_ids = set()
    if not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        return noncoding_ids
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            label = row.get("label", "").strip().lower()
            tid = row.get("ID", row.get("#ID", "")).strip()
            if label == "noncoding" and tid:
                noncoding_ids.add(tid)
    return noncoding_ids


def parse_cpat(filepath, threshold):
    """Return set of transcript IDs called noncoding by CPAT (prob < threshold)."""
    noncoding_ids = set()
    if not os.path.isfile(filepath) or os.path.getsize(filepath) == 0:
        return noncoding_ids
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # CPAT output may have different column names
            prob = None
            for col in ["coding_prob", "Coding_prob", "coding_probability"]:
                if col in row:
                    try:
                        prob = float(row[col])
                    except (ValueError, TypeError):
                        pass
                    break
            tid = row.get("seq_ID", row.get("ID", "")).strip()
            if prob is not None and prob < threshold and tid:
                noncoding_ids.add(tid)
    return noncoding_ids


def parse_cantatadb(filepath):
    """Return set of lncRNA coordinates from CANTATAdb for known annotation."""
    known = set()
    if not filepath or not os.path.isfile(filepath):
        return known
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                known.add(parts[3])  # lncRNA ID
    return known


def parse_candidate_gtf(filepath):
    """Parse candidates GTF, return dict of transcript_id -> (chr, start, end, strand, attrs)."""
    transcripts = {}
    if not os.path.isfile(filepath):
        return transcripts
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            match = re.search(r'transcript_id "([^"]+)"', parts[8])
            if match:
                tid = match.group(1)
                if tid not in transcripts:
                    transcripts[tid] = {
                        "chr": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "attrs": parts[8],
                        "lines": []
                    }
                transcripts[tid]["lines"].append(line.strip())
                # Update coordinates to span full transcript
                transcripts[tid]["start"] = min(transcripts[tid]["start"], int(parts[3]))
                transcripts[tid]["end"] = max(transcripts[tid]["end"], int(parts[4]))
    return transcripts


def extract_lncrna_counts(counts_matrix_path, lncrna_ids, output_path):
    """Extract rows matching lncRNA IDs from FLAIR counts matrix."""
    written = 0
    with open(counts_matrix_path) as fin, open(output_path, "w") as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            row_id = line.split("\t")[0]
            # Match if the row ID contains any lncRNA transcript ID
            for lid in lncrna_ids:
                if lid in row_id or row_id in lid:
                    fout.write(line)
                    written += 1
                    break
    return written


def main():
    args = parse_args()

    print("[K-CHOPORE] lncRNA Consensus Classification")
    print(f"[K-CHOPORE] Minimum consensus: {args.min_consensus}/3 tools")
    print(f"[K-CHOPORE] Max ORF size: {args.max_orf_aa} aa")
    print(f"[K-CHOPORE] CPAT threshold: {args.cpat_threshold}")

    # Parse all tool results
    coding_orfs = parse_transdecoder(args.transdecoder, args.max_orf_aa)
    print(f"[K-CHOPORE] TransDecoder: {len(coding_orfs)} transcripts with ORF >{args.max_orf_aa}aa")

    feelnc_ids = parse_feelnc_ids(args.feelnc_gtf)
    feelnc_classes = parse_feelnc_classification(args.feelnc_class)
    print(f"[K-CHOPORE] FEELnc: {len(feelnc_ids)} lncRNA candidates")

    cpc2_noncoding = parse_cpc2(args.cpc2)
    print(f"[K-CHOPORE] CPC2: {len(cpc2_noncoding)} noncoding calls")

    cpat_noncoding = parse_cpat(args.cpat, args.cpat_threshold)
    print(f"[K-CHOPORE] CPAT: {len(cpat_noncoding)} noncoding calls")

    cantata_known = parse_cantatadb(args.cantatadb)
    print(f"[K-CHOPORE] CANTATAdb: {len(cantata_known)} known lncRNAs")

    # Parse candidate transcripts
    candidates = parse_candidate_gtf(args.candidates_gtf)
    print(f"[K-CHOPORE] Candidates: {len(candidates)} transcripts")

    # Apply consensus filter
    consensus_lncrnas = {}
    stats = defaultdict(int)

    for tid, info in candidates.items():
        # Check TransDecoder: exclude if coding ORF found
        if tid in coding_orfs:
            stats["excluded_coding_orf"] += 1
            continue

        # Count tool votes for "non-coding"
        votes = 0
        tool_calls = []

        if tid in feelnc_ids:
            votes += 1
            tool_calls.append("FEELnc")
        if tid in cpc2_noncoding:
            votes += 1
            tool_calls.append("CPC2")
        if tid in cpat_noncoding:
            votes += 1
            tool_calls.append("CPAT")

        if votes >= args.min_consensus:
            # Classify type
            lnc_type = feelnc_classes.get(tid, "unclassified")
            if lnc_type in ("", "unclassified", "unknown"):
                lnc_type = "lincRNA"  # Default for intergenic

            # Check if known in CANTATAdb
            is_known = "known" if tid in cantata_known else "novel"

            consensus_lncrnas[tid] = {
                **info,
                "type": lnc_type,
                "known_status": is_known,
                "votes": votes,
                "tools": ",".join(tool_calls)
            }
            stats[f"type_{lnc_type}"] += 1
            stats[f"status_{is_known}"] += 1
        else:
            stats["excluded_low_consensus"] += 1

    print(f"\n[K-CHOPORE] === Consensus Results ===")
    print(f"[K-CHOPORE] Total lncRNAs passing consensus: {len(consensus_lncrnas)}")
    for key, val in sorted(stats.items()):
        print(f"[K-CHOPORE]   {key}: {val}")

    # Write output GTF
    with open(args.out_gtf, "w") as f:
        f.write(f"# K-CHOPORE v3.0 lncRNA consensus classification\n")
        f.write(f"# Consensus: {args.min_consensus}/3 tools (FEELnc, CPC2, CPAT)\n")
        f.write(f"# Total lncRNAs: {len(consensus_lncrnas)}\n")
        for tid, info in consensus_lncrnas.items():
            for line in info["lines"]:
                # Append lncRNA metadata to GTF attributes
                extra = f' lncrna_type "{info["type"]}"; known_status "{info["known_status"]}"; consensus_tools "{info["tools"]}";'
                f.write(line + extra + "\n")

    # Write output BED
    with open(args.out_bed, "w") as f:
        for tid, info in consensus_lncrnas.items():
            f.write(f"{info['chr']}\t{info['start']-1}\t{info['end']}\t{tid}\t{info['votes']}\t{info['strand']}\n")

    # Extract lncRNA counts from FLAIR matrix
    n_counts = extract_lncrna_counts(
        args.counts_matrix, set(consensus_lncrnas.keys()), args.out_counts
    )
    print(f"[K-CHOPORE] lncRNA counts matrix: {n_counts} features extracted")

    # Write summary report
    with open(args.out_report, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "transcript_id", "chr", "start", "end", "strand",
            "lncrna_type", "known_status", "consensus_votes", "tools"
        ])
        for tid, info in sorted(consensus_lncrnas.items()):
            writer.writerow([
                tid, info["chr"], info["start"], info["end"], info["strand"],
                info["type"], info["known_status"], info["votes"], info["tools"]
            ])

    print(f"\n[K-CHOPORE] lncRNA consensus classification completed.")
    print(f"[K-CHOPORE] Output files:")
    print(f"[K-CHOPORE]   GTF: {args.out_gtf}")
    print(f"[K-CHOPORE]   BED: {args.out_bed}")
    print(f"[K-CHOPORE]   Counts: {args.out_counts}")
    print(f"[K-CHOPORE]   Report: {args.out_report}")


if __name__ == "__main__":
    main()
