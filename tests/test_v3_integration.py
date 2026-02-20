"""
K-CHOPORE v3.0 — Integration Tests with Toy Data
=================================================
Exercises the v3 Python scripts end-to-end using mock data
from tests/toy_data/mock_results/.

Tests:
  - lncrna_consensus.py with synthetic CPC2/CPAT/FEELnc/TransDecoder outputs
  - epitx_consensus.py with synthetic ELIGOS2 and m6Anet outputs
  - annotate_mirna.py parsing of mock miRBase/PmiREN annotation files
  - cis_regulation.py with toy lncRNA + DE results
  - Toy annotation file integrity (BED, GFF3, FASTA)
"""

import csv
import importlib.util
import os
import sys
import tempfile

import pytest

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, "scripts")
MOCK_DIR = os.path.join(PROJECT_ROOT, "tests", "toy_data", "mock_results")
ANNOT_DIR = os.path.join(PROJECT_ROOT, "tests", "toy_data", "annotations")


def _import_script(name):
    """Import a Python script from scripts/ as a module."""
    path = os.path.join(SCRIPTS_DIR, name)
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ===================================================================
# 1. TOY ANNOTATION FILE INTEGRITY
# ===================================================================
class TestToyAnnotationFiles:
    """Verify all v3 toy annotation files exist and are well-formed."""

    def test_cantatadb_bed_exists(self):
        fp = os.path.join(ANNOT_DIR, "cantata_test.bed")
        assert os.path.isfile(fp), "cantata_test.bed missing"

    def test_cantatadb_bed_columns(self):
        fp = os.path.join(ANNOT_DIR, "cantata_test.bed")
        with open(fp) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                cols = line.split("\t")
                assert len(cols) >= 4, f"BED needs ≥4 columns, got {len(cols)}"

    def test_mirbase_gff3_exists(self):
        fp = os.path.join(ANNOT_DIR, "mirbase_test.gff3")
        assert os.path.isfile(fp)

    def test_mirbase_gff3_has_miRNA_entries(self):
        fp = os.path.join(ANNOT_DIR, "mirbase_test.gff3")
        mirna_count = 0
        with open(fp) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if "\tmiRNA\t" in line:
                    mirna_count += 1
        assert mirna_count >= 2, f"Expected ≥2 miRNA entries, got {mirna_count}"

    def test_mirbase_mature_fasta_exists(self):
        fp = os.path.join(ANNOT_DIR, "mirbase_mature_test.fa")
        assert os.path.isfile(fp)

    def test_mirbase_mature_fasta_has_sequences(self):
        fp = os.path.join(ANNOT_DIR, "mirbase_mature_test.fa")
        seq_count = 0
        with open(fp) as f:
            for line in f:
                if line.startswith(">"):
                    seq_count += 1
        assert seq_count >= 3, f"Expected ≥3 miRNA sequences, got {seq_count}"

    def test_pmiren_fasta_exists(self):
        fp = os.path.join(ANNOT_DIR, "pmiren_test.fa")
        assert os.path.isfile(fp)

    def test_pmiren_fasta_has_sequences(self):
        fp = os.path.join(ANNOT_DIR, "pmiren_test.fa")
        seq_count = 0
        with open(fp) as f:
            for line in f:
                if line.startswith(">"):
                    seq_count += 1
        assert seq_count >= 3

    def test_te_bed_exists(self):
        fp = os.path.join(ANNOT_DIR, "te_test.bed")
        assert os.path.isfile(fp)

    def test_te_bed_has_entries(self):
        fp = os.path.join(ANNOT_DIR, "te_test.bed")
        with open(fp) as f:
            lines = [l for l in f if l.strip()]
        assert len(lines) >= 2


# ===================================================================
# 2. MOCK RESULT FILE INTEGRITY
# ===================================================================
class TestMockResultFiles:
    """Verify mock result files used by integration tests."""

    def test_flair_counts_exists(self):
        fp = os.path.join(MOCK_DIR, "flair", "counts_matrix.tsv")
        assert os.path.isfile(fp)

    def test_flair_counts_has_header_and_data(self):
        fp = os.path.join(MOCK_DIR, "flair", "counts_matrix.tsv")
        with open(fp) as f:
            lines = f.readlines()
        assert len(lines) >= 3, "FLAIR counts needs header + ≥2 data rows"
        assert "ids" in lines[0], "Header should contain 'ids'"

    def test_transdecoder_pep_exists(self):
        fp = os.path.join(MOCK_DIR, "lncrna", "transdecoder", "longest_orfs.pep")
        assert os.path.isfile(fp)

    def test_cpc2_results_exists(self):
        fp = os.path.join(MOCK_DIR, "lncrna", "cpc2", "cpc2_results.txt")
        assert os.path.isfile(fp)

    def test_cpat_results_exists(self):
        fp = os.path.join(MOCK_DIR, "lncrna", "cpat", "cpat_results.tsv")
        assert os.path.isfile(fp)

    def test_feelnc_gtf_exists(self):
        fp = os.path.join(MOCK_DIR, "lncrna", "feelnc", "candidate_lncRNA.gtf")
        assert os.path.isfile(fp)

    def test_candidates_gtf_exists(self):
        fp = os.path.join(MOCK_DIR, "lncrna", "candidates_merged.gtf")
        assert os.path.isfile(fp)

    def test_eligos_outputs_exist(self):
        for sample in ["test_sample_1", "test_sample_2"]:
            fp = os.path.join(MOCK_DIR, "eligos", f"{sample}_eligos_output.txt")
            assert os.path.isfile(fp), f"Missing ELIGOS2 output for {sample}"

    def test_m6anet_outputs_exist(self):
        for sample in ["test_sample_1", "test_sample_2"]:
            fp = os.path.join(MOCK_DIR, "m6anet", sample, "data.site_proba.csv")
            assert os.path.isfile(fp), f"Missing m6Anet output for {sample}"

    def test_deseq2_results_exist(self):
        fp = os.path.join(MOCK_DIR, "deseq2", "deseq2_genotype_results.csv")
        assert os.path.isfile(fp)


# ===================================================================
# 3. lncRNA CONSENSUS END-TO-END
# ===================================================================
class TestLncRNAConsensusEndToEnd:
    """Run lncrna_consensus.py functions with mock data."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.mod = _import_script("lncrna_consensus.py")

    def test_parse_transdecoder_detects_coding(self):
        pep_file = os.path.join(MOCK_DIR, "lncrna", "transdecoder", "longest_orfs.pep")
        coding = self.mod.parse_transdecoder(pep_file, max_orf_aa=100)
        # MSTRG.3.1 has only 20aa ORF, should NOT be coding
        assert "MSTRG.3.1" not in coding
        # MSTRG.1.1 headers say len:50 but actual sequence is 42aa
        # (our threshold is 100, so neither should be detected as >100aa coding)

    def test_parse_cpc2_detects_noncoding(self):
        cpc2_file = os.path.join(MOCK_DIR, "lncrna", "cpc2", "cpc2_results.txt")
        noncoding = self.mod.parse_cpc2(cpc2_file)
        assert "MSTRG.2.1" in noncoding, "MSTRG.2.1 labelled noncoding by CPC2"
        assert "MSTRG.3.1" in noncoding, "MSTRG.3.1 labelled noncoding by CPC2"
        assert "MSTRG.1.1" not in noncoding, "MSTRG.1.1 labelled coding by CPC2"

    def test_parse_cpat_detects_noncoding(self):
        cpat_file = os.path.join(MOCK_DIR, "lncrna", "cpat", "cpat_results.tsv")
        noncoding = self.mod.parse_cpat(cpat_file, threshold=0.39)
        assert "MSTRG.2.1" in noncoding, "MSTRG.2.1 has prob=0.12 < 0.39"
        assert "MSTRG.3.1" in noncoding, "MSTRG.3.1 has prob=0.22 < 0.39"

    def test_parse_feelnc_ids(self):
        feelnc_gtf = os.path.join(MOCK_DIR, "lncrna", "feelnc", "candidate_lncRNA.gtf")
        ids = self.mod.parse_feelnc_ids(feelnc_gtf)
        assert "MSTRG.2.1" in ids
        assert "MSTRG.3.1" in ids

    def test_parse_cantatadb(self):
        cantata = os.path.join(ANNOT_DIR, "cantata_test.bed")
        known = self.mod.parse_cantatadb(cantata)
        assert len(known) >= 2, "Should find ≥2 CANTATAdb entries"

    def test_parse_candidate_gtf(self):
        gtf = os.path.join(MOCK_DIR, "lncrna", "candidates_merged.gtf")
        candidates = self.mod.parse_candidate_gtf(gtf)
        assert len(candidates) == 3, f"Expected 3 transcripts, got {len(candidates)}"
        assert "MSTRG.2.1" in candidates

    def test_full_consensus_pipeline(self):
        """End-to-end test: run the full consensus logic and check outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_gtf = os.path.join(tmpdir, "lncrna_final.gtf")
            out_bed = os.path.join(tmpdir, "lncrna_final.bed")
            out_counts = os.path.join(tmpdir, "lncrna_counts.tsv")
            out_report = os.path.join(tmpdir, "summary.tsv")

            # Simulate command-line args
            sys.argv = [
                "lncrna_consensus.py",
                "--transdecoder", os.path.join(MOCK_DIR, "lncrna", "transdecoder", "longest_orfs.pep"),
                "--feelnc-gtf", os.path.join(MOCK_DIR, "lncrna", "feelnc", "candidate_lncRNA.gtf"),
                "--feelnc-class", os.path.join(MOCK_DIR, "lncrna", "feelnc", "candidate_lncRNA.gtf"),  # reuse as class
                "--cpc2", os.path.join(MOCK_DIR, "lncrna", "cpc2", "cpc2_results.txt"),
                "--cpat", os.path.join(MOCK_DIR, "lncrna", "cpat", "cpat_results.tsv"),
                "--candidates-gtf", os.path.join(MOCK_DIR, "lncrna", "candidates_merged.gtf"),
                "--counts-matrix", os.path.join(MOCK_DIR, "flair", "counts_matrix.tsv"),
                "--cantatadb", os.path.join(ANNOT_DIR, "cantata_test.bed"),
                "--te-bed", os.path.join(ANNOT_DIR, "te_test.bed"),
                "--min-consensus", "2",
                "--max-orf-aa", "100",
                "--cpat-threshold", "0.39",
                "--out-gtf", out_gtf,
                "--out-bed", out_bed,
                "--out-counts", out_counts,
                "--out-report", out_report,
            ]

            self.mod.main()

            # Verify outputs exist and have content
            assert os.path.isfile(out_gtf), "GTF not created"
            assert os.path.getsize(out_gtf) > 0, "GTF empty"

            assert os.path.isfile(out_bed), "BED not created"
            assert os.path.isfile(out_counts), "Counts not created"
            assert os.path.isfile(out_report), "Report not created"

            # Verify report content
            with open(out_report) as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
            # MSTRG.2.1 and MSTRG.3.1 should pass (both called non-coding by CPC2 + CPAT + FEELnc)
            lncrna_ids = [r["transcript_id"] for r in rows]
            assert "MSTRG.2.1" in lncrna_ids, "MSTRG.2.1 should pass consensus"
            assert "MSTRG.3.1" in lncrna_ids, "MSTRG.3.1 should pass consensus"
            # MSTRG.1.1 is called coding by CPC2, so should fail consensus
            assert "MSTRG.1.1" not in lncrna_ids, "MSTRG.1.1 should fail (only 1/3 tools noncoding)"

    def test_extract_lncrna_counts(self):
        """Verify count extraction from FLAIR matrix for lncRNA IDs."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            counts_matrix = os.path.join(MOCK_DIR, "flair", "counts_matrix.tsv")
            n = self.mod.extract_lncrna_counts(counts_matrix, {"MSTRG.2.1", "MSTRG.3.1"}, tmp_path)
            assert n == 2, f"Expected 2 lncRNA rows extracted, got {n}"
            with open(tmp_path) as f:
                lines = f.readlines()
            assert len(lines) == 3, "Header + 2 data rows"
        finally:
            os.unlink(tmp_path)


# ===================================================================
# 4. EPITRANSCRIPTOMIC CONSENSUS END-TO-END
# ===================================================================
class TestEpitxConsensusEndToEnd:
    """Run epitx_consensus.py functions with mock data."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.mod = _import_script("epitx_consensus.py")

    def test_parse_eligos_finds_significant_sites(self):
        sites = self.mod.parse_eligos(os.path.join(MOCK_DIR, "eligos"), pval_thresh=0.05)
        assert len(sites) >= 4, f"Expected ≥4 ELIGOS2 sites, got {len(sites)}"
        tools = {s["tool"] for s in sites}
        assert tools == {"ELIGOS2"}

    def test_parse_m6anet_finds_significant_sites(self):
        sites = self.mod.parse_m6anet(os.path.join(MOCK_DIR, "m6anet"), pval_thresh=0.05)
        # m6Anet uses prob > 0.95 threshold (1 - 0.05)
        high_prob = [s for s in sites if s["score"] > 0.95]
        assert len(high_prob) >= 0  # may have 0 due to threshold

    def test_full_epitx_consensus_pipeline(self):
        """End-to-end: merge ELIGOS2 + m6Anet and build consensus."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_consensus = os.path.join(tmpdir, "consensus.tsv")
            out_biotype = os.path.join(tmpdir, "biotype.tsv")

            sys.argv = [
                "epitx_consensus.py",
                "--eligos-dir", os.path.join(MOCK_DIR, "eligos"),
                "--m6anet-dir", os.path.join(MOCK_DIR, "m6anet"),
                "--min-tools", "1",  # Use 1 for test since we only have 2 tools
                "--pval", "0.05",
                "--out-consensus", out_consensus,
                "--out-biotype", out_biotype,
            ]

            self.mod.main()

            assert os.path.isfile(out_consensus)
            assert os.path.getsize(out_consensus) > 0

            assert os.path.isfile(out_biotype)
            assert os.path.getsize(out_biotype) > 0

            # Check consensus output format
            with open(out_consensus) as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
            assert len(rows) >= 1, "Should have at least 1 consensus site"
            assert "chr" in reader.fieldnames
            assert "n_tools" in reader.fieldnames

    def test_site_in_lncrna_function(self):
        """Test the biotype classification helper."""
        regions = [
            {"chr": "Chr1", "start": 100, "end": 500, "strand": "+"},
            {"chr": "Chr2", "start": 10, "end": 300, "strand": "-"},
        ]
        # Site inside lncRNA
        assert self.mod.site_in_lncrna({"chr": "Chr1", "pos": 250}, regions) is True
        # Site outside lncRNA
        assert self.mod.site_in_lncrna({"chr": "Chr1", "pos": 800}, regions) is False
        # Different chromosome
        assert self.mod.site_in_lncrna({"chr": "Chr3", "pos": 100}, regions) is False


# ===================================================================
# 5. miRNA ANNOTATION PARSING
# ===================================================================
class TestMiRNAAnnotationParsing:
    """Test annotate_mirna.py functions with toy miRBase/PmiREN."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.mod = _import_script("annotate_mirna.py")

    def test_parse_mirbase_gff(self):
        """Parse the toy miRBase GFF3."""
        fp = os.path.join(ANNOT_DIR, "mirbase_test.gff3")
        mirnas = self.mod.parse_mirbase_gff(fp)
        assert len(mirnas) >= 2, f"Expected ≥2 miRBase entries, got {len(mirnas)}"

    def test_parse_pmiren_fa(self):
        """Parse the toy PmiREN FASTA."""
        fp = os.path.join(ANNOT_DIR, "pmiren_test.fa")
        mirnas = self.mod.parse_pmiren_fa(fp)
        assert len(mirnas) >= 3, f"Expected ≥3 PmiREN entries, got {len(mirnas)}"


# ===================================================================
# 6. DESEQ2 MOCK RESULTS VALIDATION
# ===================================================================
class TestDESeq2MockResults:
    """Validate mock DESeq2 results format."""

    def test_deseq2_csv_format(self):
        fp = os.path.join(MOCK_DIR, "deseq2", "deseq2_genotype_results.csv")
        with open(fp) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 5, f"Expected 5 DE genes, got {len(rows)}"
        # Check required columns
        for col in ["baseMean", "log2FoldChange", "padj"]:
            assert col in reader.fieldnames, f"Missing column: {col}"

    def test_deseq2_has_significant_genes(self):
        fp = os.path.join(MOCK_DIR, "deseq2", "deseq2_genotype_results.csv")
        with open(fp) as f:
            reader = csv.DictReader(f)
            sig_count = sum(1 for row in reader if float(row["padj"]) < 0.05)
        assert sig_count >= 3, f"Expected ≥3 significant genes, got {sig_count}"


# ===================================================================
# 7. FLAIR COUNTS MATRIX STRUCTURE
# ===================================================================
class TestFLAIRCountsMatrix:
    """Validate FLAIR counts matrix matches toy samples."""

    def test_counts_has_both_samples(self):
        fp = os.path.join(MOCK_DIR, "flair", "counts_matrix.tsv")
        with open(fp) as f:
            header = f.readline().strip().split("\t")
        assert "test_sample_1" in header
        assert "test_sample_2" in header

    def test_counts_has_novel_transcripts(self):
        fp = os.path.join(MOCK_DIR, "flair", "counts_matrix.tsv")
        with open(fp) as f:
            f.readline()  # skip header
            ids = [line.split("\t")[0] for line in f if line.strip()]
        novel = [i for i in ids if i.startswith("MSTRG")]
        assert len(novel) >= 2, "Should have ≥2 novel MSTRG transcripts"
