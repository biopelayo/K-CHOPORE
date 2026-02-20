#!/usr/bin/env python3
"""
K-CHOPORE v3.0 — Epitranscriptomics Report Generator
=====================================================
Generates an HTML summary report with modification counts per tool,
per biotype, per condition, and tool overlap visualization.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser(description="Generate epitranscriptomics report")
    parser.add_argument("--consensus", required=True, help="Consensus TSV")
    parser.add_argument("--biotype", required=True, help="Biotype stratification TSV")
    parser.add_argument("--output", required=True, help="Output HTML report")
    return parser.parse_args()


def load_tsv(filepath):
    """Load TSV as list of dicts."""
    rows = []
    if not os.path.isfile(filepath):
        return rows
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def generate_html_report(consensus_data, biotype_data, output_path):
    """Generate HTML report."""

    # Compute statistics
    total_consensus = len(consensus_data)

    # Tool counts from consensus
    tool_counts = defaultdict(int)
    for row in consensus_data:
        tools = row.get("tools", "").split(",")
        for tool in tools:
            if tool.strip():
                tool_counts[tool.strip()] += 1

    # Biotype counts
    biotype_counts = defaultdict(int)
    biotype_tool_counts = defaultdict(lambda: defaultdict(int))
    for row in biotype_data:
        bt = row.get("biotype", "unknown")
        biotype_counts[bt] += 1
        tools = row.get("tools", "").split(",")
        for tool in tools:
            if tool.strip():
                biotype_tool_counts[bt][tool.strip()] += 1

    # N-tool distribution from all sites
    ntool_dist = defaultdict(int)
    for row in biotype_data:
        try:
            n = int(row.get("n_tools", 1))
        except ValueError:
            n = 1
        ntool_dist[n] += 1

    total_all_sites = len(biotype_data)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>K-CHOPORE v3.0 — Epitranscriptomics Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 40px; background: #fafafa; color: #333; }}
        h1 {{ color: #1565C0; border-bottom: 3px solid #1565C0; padding-bottom: 10px; }}
        h2 {{ color: #1976D2; margin-top: 30px; }}
        table {{ border-collapse: collapse; margin: 15px 0; width: 100%; max-width: 600px; }}
        th {{ background: #1565C0; color: white; padding: 10px 15px; text-align: left; }}
        td {{ padding: 8px 15px; border-bottom: 1px solid #ddd; }}
        tr:nth-child(even) {{ background: #f5f5f5; }}
        .stat-box {{ background: white; border: 1px solid #ddd; border-radius: 8px; padding: 20px; margin: 10px 0; display: inline-block; min-width: 200px; }}
        .stat-number {{ font-size: 36px; font-weight: bold; color: #1565C0; }}
        .stat-label {{ font-size: 14px; color: #666; }}
        .note {{ background: #FFF3E0; border-left: 4px solid #FF9800; padding: 10px 15px; margin: 15px 0; }}
        footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #999; font-size: 12px; }}
    </style>
</head>
<body>
    <h1>K-CHOPORE v3.0 — Epitranscriptomics Report</h1>
    <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>

    <div>
        <div class="stat-box">
            <div class="stat-number">{total_all_sites}</div>
            <div class="stat-label">Total modification sites</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{total_consensus}</div>
            <div class="stat-label">Consensus sites (≥2 tools)</div>
        </div>
        <div class="stat-box">
            <div class="stat-number">{len(tool_counts)}</div>
            <div class="stat-label">Tools with detections</div>
        </div>
    </div>

    <h2>Sites per Tool</h2>
    <table>
        <tr><th>Tool</th><th>Sites</th></tr>
"""
    for tool, count in sorted(tool_counts.items(), key=lambda x: -x[1]):
        html += f"        <tr><td>{tool}</td><td>{count}</td></tr>\n"

    html += """    </table>

    <h2>Tool Agreement Distribution</h2>
    <table>
        <tr><th>Number of Tools</th><th>Sites</th></tr>
"""
    for n in sorted(ntool_dist.keys()):
        html += f"        <tr><td>{n} tool(s)</td><td>{ntool_dist[n]}</td></tr>\n"

    html += """    </table>

    <h2>Biotype Distribution</h2>
    <table>
        <tr><th>Biotype</th><th>Sites</th></tr>
"""
    for bt, count in sorted(biotype_counts.items(), key=lambda x: -x[1]):
        html += f"        <tr><td>{bt}</td><td>{count}</td></tr>\n"

    html += """    </table>

    <h2>Sites per Tool × Biotype</h2>
    <table>
        <tr><th>Biotype</th><th>Tool</th><th>Sites</th></tr>
"""
    for bt in sorted(biotype_tool_counts.keys()):
        for tool, count in sorted(biotype_tool_counts[bt].items(), key=lambda x: -x[1]):
            html += f"        <tr><td>{bt}</td><td>{tool}</td><td>{count}</td></tr>\n"

    html += f"""    </table>

    <div class="note">
        <strong>Note:</strong> Consensus sites are those detected by ≥2 independent tools.
        This multi-tool approach reduces false positives inherent in any single method.
    </div>

    <footer>
        K-CHOPORE v3.0 — ONT Direct RNA-seq Analysis Pipeline<br>
        Epitranscriptomics Module (Module 4)<br>
        Universidad de Oviedo | FINBA
    </footer>
</body>
</html>"""

    with open(output_path, "w") as f:
        f.write(html)


def main():
    args = parse_args()

    print("[K-CHOPORE] Module 4: Generating epitranscriptomics report...")

    consensus_data = load_tsv(args.consensus)
    biotype_data = load_tsv(args.biotype)

    print(f"[K-CHOPORE] Consensus sites: {len(consensus_data)}")
    print(f"[K-CHOPORE] Total sites (all tools): {len(biotype_data)}")

    generate_html_report(consensus_data, biotype_data, args.output)

    print(f"[K-CHOPORE] Report saved: {args.output}")
    print("[K-CHOPORE] Epitranscriptomics report completed.")


if __name__ == "__main__":
    main()
