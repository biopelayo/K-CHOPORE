"""
Test 3: Snakefile DAG validation
Tests that:
  - Snakemake can parse the Snakefile with the test config
  - The DAG resolves correctly (dry run)
  - All expected rules are present
  - pycoQC input resolution works with the actual filenames
"""
import os
import shutil
import subprocess
import pytest
import yaml

from tests.conftest import PROJECT_ROOT, TESTS_DIR, TOY_DATA_DIR


def _setup_toy_inputs():
    """
    Copy toy data into the locations expected by the Snakefile rules.
    The nanofilt rule expects: data/raw/fastq/{sample}.fastq
    The pycoQC rule uses the config sequencing_summaries paths (already in tests/toy_data).
    """
    # Create the data/raw/fastq directory and copy toy FASTQs
    raw_fastq = os.path.join(PROJECT_ROOT, "data", "raw", "fastq")
    os.makedirs(raw_fastq, exist_ok=True)
    for sample in ["test_sample_1", "test_sample_2"]:
        src = os.path.join(TOY_DATA_DIR, "fastq", f"{sample}.fastq")
        dst = os.path.join(raw_fastq, f"{sample}.fastq")
        if not os.path.exists(dst):
            shutil.copy2(src, dst)


class TestSnakemakeParsing:
    """Test that Snakemake can parse the Snakefile and build the DAG."""

    def test_snakemake_dry_run(self):
        """Snakemake --dry-run should succeed with test config."""
        _setup_toy_inputs()
        result = subprocess.run(
            [
                "python", "-m", "snakemake",
                "--snakefile", os.path.join(PROJECT_ROOT, "Snakefile"),
                "--configfile", os.path.join(TESTS_DIR, "test_config.yml"),
                "--dry-run",
                "--quiet",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0, (
            f"Snakemake dry run failed!\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_snakemake_list_rules(self):
        """Snakemake --list should show all expected rules."""
        result = subprocess.run(
            [
                "python", "-m", "snakemake",
                "--snakefile", os.path.join(PROJECT_ROOT, "Snakefile"),
                "--configfile", os.path.join(TESTS_DIR, "test_config.yml"),
                "--list",
            ],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        assert result.returncode == 0, (
            f"Snakemake --list failed!\n"
            f"STDERR:\n{result.stderr}"
        )
        output = result.stdout
        expected_rules = [
            "all",
            "setup_complete_structure",
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
            "flair_quantify",
            "run_eligos2",
            "nanopolish_index",
            "nanopolish_eventalign",
            "m6anet_dataprep",
            "m6anet_inference",
            "run_deseq2",
            "multiqc",
        ]
        for rule_name in expected_rules:
            assert rule_name in output, (
                f"Expected rule '{rule_name}' not found in Snakemake --list output"
            )


class TestPycoQCResolution:
    """Test that pycoQC correctly resolves sequencing summary filenames."""

    def test_sequencing_summary_mapping_logic(self):
        """Verify the SEQUENCING_SUMMARY dict is built correctly."""
        with open(os.path.join(TESTS_DIR, "test_config.yml"), "r") as f:
            config = yaml.safe_load(f)

        samples = config["samples"]
        seq_summaries = config["input_files"].get("sequencing_summaries", [])

        # Build the mapping the same way the Snakefile does
        mapping = {}
        for i, sample in enumerate(samples):
            if i < len(seq_summaries):
                mapping[sample] = seq_summaries[i]
            else:
                mapping[sample] = f"data/raw/summaries/{sample}_sequencing_summary.txt"

        # Verify each sample maps to its specific file
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
