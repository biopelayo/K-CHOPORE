#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — Multi-Omics Integration Report
=================================================
Generates a comprehensive HTML report summarizing results from
all K-CHOPORE modules: WGCNA co-expression, cis-regulatory pairs,
population variation, and miRNA-target networks.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser(description="Generate integration report")
    parser.add_argument("--wgcna", required=True, help="WGCNA modules TSV")
    parser.add_argument("--cis-pairs", default="", help="Cis lncRNA pairs TSV")
    parser.add_argument("--variation", default="", help="ncRNA variation TSV")
    parser.add_argument("--targets", default="", help="Target evidence table TSV")
    parser.add_argument("--output", required=True, help="Output HTML report")
    return parser.parse_args()


def load_tsv(filepath):
    """Load TSV as list of dicts."""
    rows = []
    if not filepath or not os.path.isfile(filepath):
        return rows
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def main():
    args = parse_args()

    print("[K-CHOPORE] Module 5: Generating Integration Report")

    # Load all data
    wgcna_data = load_tsv(args.wgcna)
    cis_data = load_tsv(args.cis_pairs)
    var_data = load_tsv(args.variation)
    target_data = load_tsv(args.targets)

    print(f"[K-CHOPORE] WGCNA genes: {len(wgcna_data)}")
    print(f"[K-CHOPORE] Cis pairs: {len(cis_data)}")
    print(f"[K-CHOPORE] Variant loci: {len(var_data)}")
    print(f"[K-CHOPORE] Target pairs: {len(target_data)}")

    # WGCNA statistics
    module_counts = defaultdict(int)
    for row in wgcna_data:
        mod = row.get("module_color", row.get("module_number", "unassigned"))
        module_counts[mod] += 1

    # Cis statistics
    n_overlapping = sum(1 for r in cis_data if r.get("distance_bp", "0") == "0")
    n_antisense = sum(1 for r in cis_data if r.get("orientation", "") == "antisense")

    # Target statistics
    n_validated = sum(1 for r in target_data if r.get("degradome", "False") == "True")
    n_lncrna_targets = sum(1 for r in target_data if r.get("is_lncrna_target", "False") == "True")

    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>K-CHOPORE v3.0 — Multi-Omics Integration Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 40px; background: #fafafa; color: #333; }}
        h1 {{ color: #2E7D32; border-bottom: 3px solid #2E7D32; padding-bottom: 10px; }}
        h2 {{ color: #388E3C; margin-top: 30px; }}
        h3 {{ color: #43A047; }}
        table {{ border-collapse: collapse; margin: 15px 0; width: 100%; max-width: 700px; }}
        th {{ background: #2E7D32; color: white; padding: 10px 15px; text-align: left; }}
        td {{ padding: 8px 15px; border-bottom: 1px solid #ddd; }}
        tr:nth-child(even) {{ background: #f5f5f5; }}
        .stat-box {{ background: white; border: 1px solid #ddd; border-radius: 8px; padding: 20px; margin: 10px; display: inline-block; min-width: 180px; text-align: center; }}
        .stat-number {{ font-size: 36px; font-weight: bold; color: #2E7D32; }}
        .stat-label {{ font-size: 14px; color: #666; }}
        .section {{ background: white; border-radius: 8px; padding: 20px; margin: 20px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
        .note {{ background: #E8F5E9; border-left: 4px solid #4CAF50; padding: 10px 15px; margin: 15px 0; }}
        .warning {{ background: #FFF3E0; border-left: 4px solid #FF9800; padding: 10px 15px; margin: 15px 0; }}
        footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #999; font-size: 12px; }}
    </style>
</head>
<body>
    <h1>K-CHOPORE v3.0 — Multi-Omics Integration Report</h1>
    <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
    <p>Arabidopsis thaliana | 2&times;2 factorial: WT vs anac017-1 &times; Control vs Antimycin A</p>

    <div>
        <div class="stat-box">
            <div class="stat-number">{len(module_counts)}</div>
            <div class="stat-label">WGCNA Modules</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{len(wgcna_data)}</div>
            <div class="stat-label">Genes in Network</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{len(cis_data)}</div>
            <div class="stat-label">Cis lncRNA-Gene Pairs</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{len(target_data)}</div>
            <div class="stat-label">miRNA-Target Pairs</div>
        </div>
    </div>

    <div class="section">
        <h2>1. WGCNA Co-expression Network</h2>
        <p>Weighted gene co-expression network analysis identified co-regulated
        gene modules across all RNA types (mRNA, lncRNA, miRNA).</p>

        <h3>Module Sizes</h3>
        <table>
            <tr><th>Module</th><th>Genes</th></tr>
"""

    for mod, count in sorted(module_counts.items(), key=lambda x: -x[1]):
        html += f"            <tr><td>{mod}</td><td>{count}</td></tr>\n"

    html += f"""        </table>
    </div>

    <div class="section">
        <h2>2. Cis-Regulatory lncRNA Pairs</h2>
        <p>lncRNAs located within 10 kb of differentially expressed protein-coding genes,
        representing potential cis-regulatory relationships.</p>

        <table>
            <tr><th>Metric</th><th>Count</th></tr>
            <tr><td>Total cis pairs</td><td>{len(cis_data)}</td></tr>
            <tr><td>Overlapping lncRNA-gene</td><td>{n_overlapping}</td></tr>
            <tr><td>Antisense orientation</td><td>{n_antisense}</td></tr>
        </table>
"""

    if cis_data:
        html += """
        <h3>Top Cis Pairs (by proximity)</h3>
        <table>
            <tr><th>lncRNA</th><th>Gene</th><th>Distance (bp)</th><th>Orientation</th></tr>
"""
        for row in cis_data[:20]:
            html += f"""            <tr>
                <td>{row.get('lncrna_id', '')}</td>
                <td>{row.get('gene_id', '')}</td>
                <td>{row.get('distance_bp', '')}</td>
                <td>{row.get('orientation', '')}</td>
            </tr>\n"""
        html += "        </table>\n"

    html += f"""    </div>

    <div class="section">
        <h2>3. miRNA-Target Network</h2>
        <table>
            <tr><th>Metric</th><th>Count</th></tr>
            <tr><td>Total miRNA-target pairs</td><td>{len(target_data)}</td></tr>
            <tr><td>Degradome-validated</td><td>{n_validated}</td></tr>
            <tr><td>lncRNA targets</td><td>{n_lncrna_targets}</td></tr>
        </table>
    </div>

    <div class="section">
        <h2>4. Population Variation Context</h2>
        <p>ncRNA loci intersected with 1001 Genomes variation data.</p>
        <table>
            <tr><th>Metric</th><th>Count</th></tr>
            <tr><td>Loci with variants</td><td>{len(var_data)}</td></tr>
        </table>
"""

    if var_data:
        html += """
        <h3>Most Variable Loci</h3>
        <table>
            <tr><th>Locus</th><th>Variants</th></tr>
"""
        for row in var_data[:15]:
            html += f"""            <tr><td>{row.get('locus', '')}</td><td>{row.get('variant_count', '')}</td></tr>\n"""
        html += "        </table>\n"

    html += f"""    </div>

    <div class="note">
        <strong>Methods note:</strong> This integration report combines results from
        K-CHOPORE v3.0 modules: Module 1 (lncRNA discovery), Module 2 (small RNA analysis),
        Module 3 (miRNA target prediction), Module 4 (epitranscriptomics), and
        Module 5 (data integration). All analyses use the same 2&times;2 factorial
        experimental design.
    </div>

    <footer>
        K-CHOPORE v3.0 — ONT Direct RNA-seq Analysis Pipeline<br>
        Integration Module (Module 5)<br>
        Universidad de Oviedo | FINBA | {datetime.now().year}
    </footer>
</body>
</html>"""

    with open(args.output, "w") as f:
        f.write(html)

    print(f"[K-CHOPORE] Integration report: {args.output}")
    print("[K-CHOPORE] Integration report completed.")


if __name__ == "__main__":
    main()
