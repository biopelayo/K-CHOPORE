"""
Test 1: Config validation
Ensures both the main and test configs have all required keys,
correct types, and consistent values.

Updated for K-CHOPORE v3.0: includes v2 + v3 module config checks.
"""
import os
import pytest
import yaml


# ----------------------------------------------------------------
# Required config structure (keys that MUST exist)
# ----------------------------------------------------------------
# Note: main config uses dict-style samples (v2), test config uses list (v1 legacy).
# We test both patterns.

REQUIRED_TOOLS_KEYS = [
    "minimap2_preset", "minimap2_kmer_size", "minimap2_extra_flags",
    "nanofilt_min_quality", "nanofilt_min_length",
    "flair_support",
    "eligos2_pval", "eligos2_oddR", "eligos2_esb",
    "m6anet_num_iterations", "m6anet_num_processors",
    "deseq2_padj_threshold", "deseq2_lfc_threshold",
]

# v3.0 module toggle flags
V3_MODULE_TOGGLES = [
    "run_lncrna", "run_smallrna", "run_mirna_targets",
    "run_epitx_enhanced", "run_integration",
]

# v3.0 lncRNA parameters
V3_LNCRNA_PARAMS = [
    "lncrna_min_length", "lncrna_consensus_min", "cpat_threshold",
]


class TestConfigStructure:
    """Validate that configs have all required keys."""

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_top_level_has_samples(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        assert "samples" in config, f"Missing 'samples' in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_top_level_has_input_files(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        assert "input_files" in config, f"Missing 'input_files' in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_top_level_has_tools(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        assert "tools" in config, f"Missing 'tools' in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_top_level_has_params(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        assert "params" in config, f"Missing 'params' in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_tools_keys(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for key in REQUIRED_TOOLS_KEYS:
            assert key in config["tools"], f"Missing tools.{key} in {config_fixture}"


class TestConfigValues:
    """Validate that config values are correct types and within expected ranges."""

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_samples_is_nonempty(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        samples = config["samples"]
        if isinstance(samples, dict):
            assert len(samples) > 0, "samples dict must not be empty"
        elif isinstance(samples, list):
            assert len(samples) > 0, "samples list must not be empty"
        else:
            pytest.fail(f"samples must be dict or list, got {type(samples)}")

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
            f"minimap2_extra_flags must contain '--MD', got '{extra}'"
        )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_threads_positive(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        threads = config["params"]["threads"]
        assert isinstance(threads, int) and threads > 0, (
            f"threads must be a positive integer, got {threads}"
        )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_run_flags_are_booleans(self, config_fixture, request):
        config = request.getfixturevalue(config_fixture)
        for key, val in config["params"].items():
            if key.startswith("run_"):
                assert isinstance(val, bool), (
                    f"params.{key} should be boolean, got {type(val).__name__}: {val}"
                )


class TestV3ModuleConfig:
    """Validate K-CHOPORE v3.0 module configuration."""

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_v3_module_toggles_present(self, config_fixture, request):
        """All v3.0 module toggle flags must exist as booleans."""
        config = request.getfixturevalue(config_fixture)
        params = config.get("params", {})
        for toggle in V3_MODULE_TOGGLES:
            assert toggle in params, f"Missing params.{toggle} in {config_fixture}"
            assert isinstance(params[toggle], bool), (
                f"params.{toggle} must be bool in {config_fixture}"
            )

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_v3_annotations_section(self, config_fixture, request):
        """Annotations section must exist with required keys."""
        config = request.getfixturevalue(config_fixture)
        assert "annotations" in config, f"Missing 'annotations' in {config_fixture}"
        ann = config["annotations"]
        for key in ["araport11_gff", "cantatadb_bed", "mirbase_gff", "pmiren_fa", "te_annotation"]:
            assert key in ann, f"Missing annotations.{key} in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_v3_lncrna_params(self, config_fixture, request):
        """lncRNA module parameters must be present."""
        config = request.getfixturevalue(config_fixture)
        params = config.get("params", {})
        for key in V3_LNCRNA_PARAMS:
            assert key in params, f"Missing params.{key} in {config_fixture}"

    @pytest.mark.parametrize("config_fixture", ["main_config", "test_config"])
    def test_v3_epitx_comparisons(self, config_fixture, request):
        """epitx_comparisons must be a list."""
        config = request.getfixturevalue(config_fixture)
        assert "epitx_comparisons" in config, f"Missing 'epitx_comparisons' in {config_fixture}"
        comps = config["epitx_comparisons"]
        assert isinstance(comps, list), "epitx_comparisons must be a list"

    def test_main_config_lncrna_defaults(self, main_config):
        """Main config lncRNA defaults should be sensible."""
        params = main_config["params"]
        assert params["lncrna_min_length"] >= 200
        assert 0 < params["cpat_threshold"] < 1
        assert params["lncrna_consensus_min"] in (1, 2, 3)
