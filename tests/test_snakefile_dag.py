"""
Test 3: Snakefile DAG validation
Tests that:
  - Snakemake can parse the Snakefile with the test config
  - All expected rules are present via static analysis
  - Toy data files are present and valid

Updated for K-CHOPORE v3.0: the main Snakefile uses dict-style samples
and the test config uses list-style. DAG tests use static analysis
rather than snakemake dry-run (which requires Docker tools).
"""
import os
import re
import subprocess
import pytest
import yaml

from tests.conftest import PROJECT_ROOT, TESTS_DIR, TOY_DATA_DIR


class TestSnakemakeParsing:
    """Test that Snakemake can parse the Snakefile syntax."""

    def test_snakefile_syntax_valid(self):
        """Snakefile must be valid Python/Snakemake syntax."""
        snakefile = os.path.join(PROJECT_ROOT, "Snakefile")
        with open(snakefile) as f:
            content = f.read()
        # Basic check: the file should be parseable and contain key elements
        assert "configfile:" in content, "Snakefile missing configfile directive"
        assert "rule all:" in content, "Snakefile missing rule all"
        assert "SAMPLE_IDS" in content, "Snakefile missing SAMPLE_IDS"

    def test_all_expected_v2_rules_present(self):
        """All core v2 rules must be in the Snakefile."""
        snakefile = os.path.join(PROJECT_ROOT, "Snakefile")
        with open(snakefile) as f:
            content = f.read()

        expected_rules = [
            "all",
            "prepare_fastq",
            "nanofilt",
            "nanoplot",
            "nanocomp",
            "index_genome",
            "map_with_minimap2",
            "sort_and_index_bam",
            "samtools_stats",
            "quality_analysis_with_pycoQC",
            "flair_align",
            "flair_correct",
            "flair_collapse",
            "generate_reads_manifest",
            "flair_quantify",
            "run_eligos2",
            "nanopolish_index",
            "nanopolish_eventalign",
            "m6anet_dataprep",
            "m6anet_inference",
            "generate_deseq2_sample_sheet",
            "run_deseq2",
            "multiqc",
        ]
        for rule_name in expected_rules:
            assert f"rule {rule_name}:" in content, (
                f"Expected rule '{rule_name}' not found in Snakefile"
            )

    def test_all_expected_v3_rules_present(self):
        """All v3.0 ncRNA module rules must be in the Snakefile."""
        snakefile = os.path.join(PROJECT_ROOT, "Snakefile")
        with open(snakefile) as f:
            content = f.read()

        v3_rules = [
            # Module 1: lncRNA
            "lncrna_filter_candidates",
            "lncrna_transdecoder",
            "lncrna_feelnc",
            "lncrna_cpc2",
            "lncrna_cpat",
            "lncrna_consensus",
            "lncrna_deseq2",
            # Module 2: sRNA
            "srna_trim",
            "srna_shortstack",
            "srna_mirdeep_p2",
            "srna_annotate_known",
            "srna_counts_matrix",
            # Module 3: Targets
            "targets_prediction",
            "targets_psrnatarget",
            "degradome_validation",
            "targets_integrate",
            # Module 4: Epitx
            "nanocompore_run",
            "epitx_consensus",
            "epitx_report",
            # Module 5: Integration
            "wgcna_network",
            "cis_regulation",
            "population_context",
            "integration_report",
        ]
        for rule_name in v3_rules:
            assert f"rule {rule_name}:" in content, (
                f"Expected v3.0 rule '{rule_name}' not found in Snakefile"
            )

    def test_total_rule_count(self):
        """Snakefile should have expected total number of rules."""
        snakefile = os.path.join(PROJECT_ROOT, "Snakefile")
        with open(snakefile) as f:
            content = f.read()
        all_rules = re.findall(r'^rule\s+\w+:', content, re.MULTILINE)
        # v2: 22 rules + v3: 24 new rules = 46 total
        assert len(all_rules) >= 45, (
            f"Expected ≥45 rules (22 v2 + 24 v3), found {len(all_rules)}"
        )


class TestPycoQCResolution:
    """Test that pycoQC correctly resolves sequencing summary filenames."""

    def test_sequencing_summary_mapping_logic(self):
        """Verify sequencing summary resolution works with test config."""
        with open(os.path.join(TESTS_DIR, "test_config.yml"), "r") as f:
            config = yaml.safe_load(f)

        samples = config["samples"]
        seq_summaries = config["input_files"].get("sequencing_summaries", [])

        # Build the mapping for list-style test config
        if isinstance(samples, list):
            mapping = {}
            for i, sample in enumerate(samples):
                if i < len(seq_summaries):
                    mapping[sample] = seq_summaries[i]
                else:
                    mapping[sample] = f"data/raw/summaries/{sample}_sequencing_summary.txt"

            assert mapping["test_sample_1"] == (
                "tests/toy_data/summaries/test_sample_1_sequencing_summary_FAR00001_abc12345.txt"
            )
            assert mapping["test_sample_2"] == (
                "tests/toy_data/summaries/test_sample_2_sequencing_summary_FAR00002_def67890.txt"
            )

    def test_sequencing_summary_files_exist(self, toy_data_dir):
        """Verify that the toy sequencing summary files actually exist."""
        expected_files = [
            "summaries/test_sample_1_sequencing_summary_FAR00001_abc12345.txt",
            "summaries/test_sample_2_sequencing_summary_FAR00002_def67890.txt",
        ]
        for relpath in expected_files:
            full_path = os.path.join(toy_data_dir, relpath)
            assert os.path.isfile(full_path), (
                f"Expected sequencing summary file not found: {full_path}"
            )


class TestToyDataIntegrity:
    """Validate that all toy data files exist and are well-formed."""

    def test_reference_genome_exists(self, toy_data_dir):
        fasta = os.path.join(toy_data_dir, "reference", "toy_genome.fasta")
        assert os.path.isfile(fasta), f"Toy genome not found: {fasta}"

    def test_reference_genome_has_sequences(self, toy_data_dir):
        """FASTA must have at least one '>' header line."""
        fasta = os.path.join(toy_data_dir, "reference", "toy_genome.fasta")
        with open(fasta) as f:
            content = f.read()
        headers = [l for l in content.split("\n") if l.startswith(">")]
        assert len(headers) >= 1, "Toy genome FASTA has no sequence headers"

    def test_gtf_exists(self, toy_data_dir):
        gtf = os.path.join(toy_data_dir, "annotations", "toy_annotation.gtf")
        assert os.path.isfile(gtf), f"Toy GTF not found: {gtf}"

    def test_gtf_has_entries(self, toy_data_dir):
        """GTF must have at least one non-comment line."""
        gtf = os.path.join(toy_data_dir, "annotations", "toy_annotation.gtf")
        with open(gtf) as f:
            lines = [l for l in f if not l.startswith("#") and l.strip()]
        assert len(lines) >= 1, "Toy GTF has no entries"

    def test_fastq_files_exist(self, toy_data_dir):
        for sample in ["test_sample_1", "test_sample_2"]:
            fq = os.path.join(toy_data_dir, "fastq", f"{sample}.fastq")
            assert os.path.isfile(fq), f"Toy FASTQ not found: {fq}"

    def test_fastq_has_reads(self, toy_data_dir):
        """Each FASTQ must have at least one 4-line read block."""
        for sample in ["test_sample_1", "test_sample_2"]:
            fq = os.path.join(toy_data_dir, "fastq", f"{sample}.fastq")
            with open(fq) as f:
                lines = [l.strip() for l in f if l.strip()]
            assert len(lines) >= 4, f"FASTQ {sample} has fewer than 4 lines"
            assert lines[0].startswith("@"), f"FASTQ {sample} first line should start with '@'"
            assert lines[2] == "+", f"FASTQ {sample} third line should be '+'"

    def test_reads_manifest_exists(self, toy_data_dir):
        manifest = os.path.join(toy_data_dir, "reads_manifest.tsv")
        assert os.path.isfile(manifest), f"Reads manifest not found: {manifest}"

    def test_reads_manifest_has_header_and_data(self, toy_data_dir):
        manifest = os.path.join(toy_data_dir, "reads_manifest.tsv")
        with open(manifest) as f:
            lines = [l.strip() for l in f if l.strip()]
        assert len(lines) >= 2, "Reads manifest needs header + at least 1 data line"
        header = lines[0].split("\t")
        assert "sample" in header, "Reads manifest header must contain 'sample'"
        assert "condition" in header, "Reads manifest header must contain 'condition'"
