"""
Test 1: Config validation
Ensures both the main and test configs have all required keys,
correct types, and consistent values.
"""
import os
import pytest
import yaml


# ----------------------------------------------------------------
# Required config structure (keys that MUST exist)
# ----------------------------------------------------------------
REQUIRED_TOP_KEYS = ["samples", "conditions", "input_files", "output", "tools", "params"]

REQUIRED_INPUT_FILES_KEYS = [
    "fastq_dir", "fast5_dir", "pod5_dir",
    "reference_genome", "reference_genome_mmi", "gtf_file",
    "reads_manifest",
]

REQUIRED_TOOLS_KEYS = [
    "minimap2_preset", "minimap2_kmer_size", "minimap2_extra_flags",
    "nanofilt_min_quality", "nanofilt_min_length",
    "flair_support",
    "eligos2_pval", "eligos2_oddR", "eligos2_esb",
    "m6anet_num_iterations", "m6anet_num_processors",
    "deseq2_padj_threshold", "deseq2_lfc_threshold",
]

REQUIRED_PARAMS_KEYS = [
    "threads",
    "run_basecalling", "run_nanofilt", "run_nanoplot", "run_nanocomp",
    "run_pycoqc", "run_flair", "run_eligos2", "run_m6anet", "run_deseq2",
    "run_multiqc",
]


class TestConfigStructure:
    """Validate that configs have all required keys."""

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_top_level_keys(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for key in REQUIRED_TOP_KEYS:
            assert key in config, f"Missing top-level key '{key}' in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_input_files_keys(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for key in REQUIRED_INPUT_FILES_KEYS:
            assert key in config["input_files"], (
                f"Missing input_files.{key} in {config_fixture}"
            )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_tools_keys(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for key in REQUIRED_TOOLS_KEYS:
            assert key in config["tools"], f"Missing tools.{key} in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_params_keys(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for key in REQUIRED_PARAMS_KEYS:
            assert key in config["params"], f"Missing params.{key} in {config_fixture}"


class TestConfigValues:
    """Validate that config values are correct types and within expected ranges."""

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_samples_is_nonempty_list(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        assert isinstance(config["samples"], list), "samples must be a list"
        assert len(config["samples"]) > 0, "samples list must not be empty"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_conditions_match_samples(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for sample in config["samples"]:
            assert sample in config["conditions"], (
                f"Sample '{sample}' not found in conditions mapping"
            )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_minimap2_preset_is_splice(self, config_fixture, request):
        """Critical: direct RNA-seq MUST use 'splice' preset, NOT 'map-ont'."""
        config = request.getfixturevalue(config_fixture)
        preset = config["tools"]["minimap2_preset"]
        assert preset == "splice", (
            f"minimap2_preset must be 'splice' for direct RNA-seq, got '{preset}'"
        )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_minimap2_kmer_is_14(self, config_fixture, request):
        """k=14 is recommended for direct RNA."""
        config = request.getfixturevalue(config_fixture)
        kmer = config["tools"]["minimap2_kmer_size"]
        assert kmer == 14, f"minimap2_kmer_size should be 14 for RNA, got {kmer}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_minimap2_extra_flags_contain_MD(self, config_fixture, request):
        """--MD is needed for downstream modification tools."""
        config = request.getfixturevalue(config_fixture)
        extra = config["tools"]["minimap2_extra_flags"]
        assert "--MD" in extra, (
            f"minimap2_extra_flags must contain '--MD' for modification detection, got '{extra}'"
        )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_threads_positive(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        threads = config["params"]["threads"]
        assert isinstance(threads, int) and threads > 0, (
            f"threads must be a positive integer, got {threads}"
        )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_sequencing_summaries_match_samples(self, config_fixture, request):
        """Each sample should have a corresponding sequencing summary."""
        config = request.getfixturevalue(config_fixture)
        summaries = config["input_files"].get("sequencing_summaries", [])
        samples = config["samples"]
        if summaries:
            assert len(summaries) == len(samples), (
                f"Number of sequencing_summaries ({len(summaries)}) "
                f"doesn't match number of samples ({len(samples)})"
            )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_run_flags_are_booleans(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for key, val in config["params"].items():
            if key.startswith("run_"):
                assert isinstance(val, bool), (
                    f"params.{key} should be boolean, got {type(val).__name__}: {val}"
                )
