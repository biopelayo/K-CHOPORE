#!/usr/bin/env bash
# =============================================================
# K-CHOPORE v3.0 — Download Reference Annotations
# =============================================================
# Downloads and versions all reference annotations required by
# the ncRNA expansion modules (lncRNA, sRNA, epitranscriptomics,
# targets, and integration).
#
# Usage:
#   bash scripts/download_annotations.sh [target_dir]
#
# Default target: data/reference/annotations/
# =============================================================

set -euo pipefail

TARGET_DIR="${1:-data/reference/annotations}"
mkdir -p "$TARGET_DIR"

echo "[K-CHOPORE] ================================================="
echo "[K-CHOPORE] Downloading reference annotations for K-CHOPORE v3"
echo "[K-CHOPORE] Target directory: $TARGET_DIR"
echo "[K-CHOPORE] ================================================="

# ------------------------------------------------------------------
# 1. Araport11 GFF3 (comprehensive gene annotation with lncRNAs + TEs)
# ------------------------------------------------------------------
ARAPORT11_FILE="$TARGET_DIR/Araport11_GFF3_genes_transposons.gff"
if [ ! -f "$ARAPORT11_FILE" ]; then
    echo "[K-CHOPORE] Downloading Araport11 GFF3..."
    wget -q -O "${ARAPORT11_FILE}.gz" \
        "https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.current.gff.gz" \
        2>/dev/null || \
    wget -q -O "${ARAPORT11_FILE}.gz" \
        "https://arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.Mar2021.gff.gz" \
        2>/dev/null || \
    echo "[K-CHOPORE] WARNING: Araport11 auto-download failed. Please download manually from TAIR."
    if [ -f "${ARAPORT11_FILE}.gz" ]; then
        gunzip -f "${ARAPORT11_FILE}.gz"
        echo "[K-CHOPORE] Araport11 GFF3 downloaded: $(wc -l < "$ARAPORT11_FILE") lines"
    fi
else
    echo "[K-CHOPORE] Araport11 GFF3 already present."
fi

# ------------------------------------------------------------------
# 2. CANTATAdb v3 — Known lncRNA coordinates for Arabidopsis
# ------------------------------------------------------------------
CANTATA_FILE="$TARGET_DIR/CANTATAdb_v3_arabidopsis.bed"
if [ ! -f "$CANTATA_FILE" ]; then
    echo "[K-CHOPORE] Downloading CANTATAdb v3 lncRNA annotations..."
    # CANTATAdb may require manual download from http://cantata.amu.edu.pl/
    # Attempt direct download; create placeholder if unavailable
    wget -q -O "$CANTATA_FILE" \
        "http://cantata.amu.edu.pl/download/Arabidopsis_thaliana_lncRNAs.bed" \
        2>/dev/null || \
    {
        echo "[K-CHOPORE] WARNING: CANTATAdb auto-download failed."
        echo "[K-CHOPORE] Please download manually from http://cantata.amu.edu.pl/"
        echo "# CANTATAdb v3 Arabidopsis lncRNAs — download manually" > "$CANTATA_FILE"
        echo "# Source: http://cantata.amu.edu.pl/" >> "$CANTATA_FILE"
    }
else
    echo "[K-CHOPORE] CANTATAdb BED already present."
fi

# ------------------------------------------------------------------
# 3. miRBase v22 — Known miRNAs for Arabidopsis thaliana
# ------------------------------------------------------------------
MIRBASE_GFF="$TARGET_DIR/miRBase_v22_ath.gff3"
MIRBASE_FA="$TARGET_DIR/miRBase_v22_ath_mature.fa"
if [ ! -f "$MIRBASE_GFF" ]; then
    echo "[K-CHOPORE] Downloading miRBase v22 GFF3..."
    wget -q -O "$TARGET_DIR/miRBase_v22_all.gff3" \
        "https://www.mirbase.org/download/CURRENT/genomes/ath.gff3" \
        2>/dev/null || \
    wget -q -O "$TARGET_DIR/miRBase_v22_all.gff3" \
        "https://mirbase.org/download/ath.gff3" \
        2>/dev/null || \
    echo "[K-CHOPORE] WARNING: miRBase GFF3 auto-download failed."
    if [ -f "$TARGET_DIR/miRBase_v22_all.gff3" ]; then
        mv "$TARGET_DIR/miRBase_v22_all.gff3" "$MIRBASE_GFF"
        echo "[K-CHOPORE] miRBase GFF3 downloaded: $(wc -l < "$MIRBASE_GFF") lines"
    fi
else
    echo "[K-CHOPORE] miRBase GFF3 already present."
fi

if [ ! -f "$MIRBASE_FA" ]; then
    echo "[K-CHOPORE] Downloading miRBase mature miRNA sequences..."
    wget -q -O "$TARGET_DIR/mature_all.fa.gz" \
        "https://www.mirbase.org/download/CURRENT/mature.fa.gz" \
        2>/dev/null || \
    wget -q -O "$TARGET_DIR/mature_all.fa.gz" \
        "https://mirbase.org/download/mature.fa.gz" \
        2>/dev/null || true
    if [ -f "$TARGET_DIR/mature_all.fa.gz" ]; then
        # Extract only Arabidopsis thaliana entries (ath-miR*)
        gunzip -f "$TARGET_DIR/mature_all.fa.gz"
        awk '/^>ath-/{p=1; print; next} /^>/{p=0} p' "$TARGET_DIR/mature_all.fa" > "$MIRBASE_FA"
        rm -f "$TARGET_DIR/mature_all.fa"
        echo "[K-CHOPORE] miRBase mature extracted: $(grep -c '^>' "$MIRBASE_FA") Arabidopsis miRNAs"
    fi
else
    echo "[K-CHOPORE] miRBase mature FA already present."
fi

# ------------------------------------------------------------------
# 4. PmiREN — Plant miRNA Encyclopedia
# ------------------------------------------------------------------
PMIREN_FILE="$TARGET_DIR/PmiREN_ath_miRNA.fa"
if [ ! -f "$PMIREN_FILE" ]; then
    echo "[K-CHOPORE] Downloading PmiREN Arabidopsis miRNAs..."
    wget -q -O "$PMIREN_FILE" \
        "https://www.pmiren.com/ftp/Arabidopsis_thaliana/AT_miR.fa" \
        2>/dev/null || \
    {
        echo "[K-CHOPORE] WARNING: PmiREN auto-download failed."
        echo "[K-CHOPORE] Please download manually from https://www.pmiren.com/"
        echo "> placeholder_pmiren" > "$PMIREN_FILE"
    }
else
    echo "[K-CHOPORE] PmiREN FASTA already present."
fi

# ------------------------------------------------------------------
# 5. TAIR10 TE annotation (Transposable Elements)
# ------------------------------------------------------------------
TE_FILE="$TARGET_DIR/TAIR10_TE_annotation.bed"
if [ ! -f "$TE_FILE" ]; then
    echo "[K-CHOPORE] Downloading TAIR10 TE annotation..."
    wget -q -O "$TARGET_DIR/TAIR10_Transposable_Elements.txt" \
        "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.txt" \
        2>/dev/null || \
    wget -q -O "$TARGET_DIR/TAIR10_TE.gff3" \
        "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_TE.gff3" \
        2>/dev/null || \
    echo "[K-CHOPORE] WARNING: TAIR10 TE auto-download failed."
    # Convert to BED format if GFF3 was downloaded
    if [ -f "$TARGET_DIR/TAIR10_TE.gff3" ]; then
        awk -F'\t' '$3 == "transposable_element" {print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' \
            "$TARGET_DIR/TAIR10_TE.gff3" > "$TE_FILE"
        echo "[K-CHOPORE] TE BED created: $(wc -l < "$TE_FILE") elements"
    elif [ -f "$TARGET_DIR/TAIR10_Transposable_Elements.txt" ]; then
        # Tab-delimited: chr, start, end, id, family, superfamily, strand
        tail -n +2 "$TARGET_DIR/TAIR10_Transposable_Elements.txt" | \
            awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7}' > "$TE_FILE"
        echo "[K-CHOPORE] TE BED created: $(wc -l < "$TE_FILE") elements"
    fi
else
    echo "[K-CHOPORE] TAIR10 TE BED already present."
fi

# ------------------------------------------------------------------
# Generate checksums
# ------------------------------------------------------------------
echo "[K-CHOPORE] Generating checksums..."
cd "$TARGET_DIR"
md5sum *.gff *.gff3 *.bed *.fa *.fasta 2>/dev/null | sort > CHECKSUMS.md5 || true
cd - > /dev/null

echo ""
echo "[K-CHOPORE] ================================================="
echo "[K-CHOPORE] Annotation download complete."
echo "[K-CHOPORE] Files in $TARGET_DIR:"
ls -lh "$TARGET_DIR/" 2>/dev/null || true
echo "[K-CHOPORE] ================================================="
echo ""
echo "[K-CHOPORE] NOTE: Some databases may require manual download."
echo "[K-CHOPORE] Check any WARNING messages above."
echo "[K-CHOPORE] CANTATAdb: http://cantata.amu.edu.pl/"
echo "[K-CHOPORE] PmiREN: https://www.pmiren.com/"
echo "[K-CHOPORE] TAIR: https://www.arabidopsis.org/"
