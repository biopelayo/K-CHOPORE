"""
Test Suite: Module 2 — Small RNA Analysis
Tests for small RNA processing, miRNA annotation, and count matrix generation.
"""
import os
import pytest
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))


class TestSmallRNAScripts:
    """Test sRNA analysis helper scripts."""

    def test_annotate_mirna_exists(self, project_root):
        """annotate_mirna.py must exist."""
        script = os.path.join(project_root, "scripts", "annotate_mirna.py")
        assert os.path.isfile(script)

    def test_srna_counts_exists(self, project_root):
        """srna_counts.py must exist."""
        script = os.path.join(project_root, "scripts", "srna_counts.py")
        assert os.path.isfile(script)

    def test_annotate_mirna_importable(self):
        """annotate_mirna.py must have no syntax errors."""
        import importlib.util
        script_path = os.path.join(
            os.path.dirname(__file__), "..", "scripts", "annotate_mirna.py"
        )
        spec = importlib.util.spec_from_file_location("annotate_mirna", script_path)
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
        except SystemExit:
            pass

    def test_parse_mirbase_gff_empty(self, tmp_path):
        """miRBase parser handles empty file."""
        from annotate_mirna import parse_mirbase_gff
        result = parse_mirbase_gff(str(tmp_path / "nonexistent.gff3"))
        assert result == {}

    def test_parse_pmiren_fa(self, tmp_path):
        """PmiREN parser correctly reads FASTA."""
        from annotate_mirna import parse_pmiren_fa
        fa = tmp_path / "test.fa"
        fa.write_text(">ath-miR156a\nUGACAGAAGAGAGUGAGCAC\n>ath-miR172a\nAGAAUCUUGAUGAUGCUGCAU\n")
        result = parse_pmiren_fa(str(fa))
        assert len(result) == 2
        assert "ath-miR156a" in result
        assert result["ath-miR156a"] == "UGACAGAAGAGAGUGAGCAC"


class TestSmallRNASnakefileRules:
    """Test that sRNA rules are properly defined."""

    def test_srna_rules_present(self, snakefile_content):
        """All 5 sRNA rules must be defined."""
        expected = [
            "srna_trim",
            "srna_shortstack",
            "srna_mirdeep_p2",
            "srna_annotate_known",
            "srna_counts_matrix"
        ]
        for rule in expected:
            assert f"rule {rule}:" in snakefile_content, f"Rule '{rule}' missing"

    def test_srna_toggle(self, snakefile_content):
        """RUN_SMALLRNA toggle must exist."""
        assert "RUN_SMALLRNA" in snakefile_content

    def test_srna_sample_parsing(self, snakefile_content):
        """SRNA_SAMPLES must be parsed from config."""
        assert "SRNA_SAMPLES" in snakefile_content
        assert "SRNA_SAMPLE_IDS" in snakefile_content


class TestModuleTargetsRules:
    """Test Module 3 target prediction rules."""

    def test_target_rules_present(self, snakefile_content):
        """Module 3 rules must be defined."""
        expected = [
            "targets_prediction",
            "targets_psrnatarget",
            "degradome_validation",
            "targets_integrate"
        ]
        for rule in expected:
            assert f"rule {rule}:" in snakefile_content, f"Rule '{rule}' missing"

    def test_target_scripts_exist(self, project_root):
        """Target prediction scripts must exist."""
        for script in ["predict_targets.py", "integrate_targets.py"]:
            path = os.path.join(project_root, "scripts", script)
            assert os.path.isfile(path), f"scripts/{script} not found"


class TestModuleEpitxRules:
    """Test Module 4 epitranscriptomics rules."""

    def test_epitx_rules_present(self, snakefile_content):
        """Module 4 rules must be defined."""
        expected = [
            "epitx_consensus",
            "epitx_report"
        ]
        for rule in expected:
            assert f"rule {rule}:" in snakefile_content, f"Rule '{rule}' missing"

    def test_epitx_scripts_exist(self, project_root):
        """Epitx scripts must exist."""
        for script in ["epitx_consensus.py", "epitx_report.py"]:
            path = os.path.join(project_root, "scripts", script)
            assert os.path.isfile(path), f"scripts/{script} not found"


class TestModuleIntegrationRules:
    """Test Module 5 integration rules."""

    def test_integration_rules_present(self, snakefile_content):
        """Module 5 rules must be defined."""
        expected = [
            "wgcna_network",
            "cis_regulation",
            "population_context",
            "integration_report"
        ]
        for rule in expected:
            assert f"rule {rule}:" in snakefile_content, f"Rule '{rule}' missing"

    def test_integration_scripts_exist(self, project_root):
        """Integration scripts must exist."""
        for script in ["run_wgcna.R", "cis_regulation.py", "integration_report.py"]:
            path = os.path.join(project_root, "scripts", script)
            assert os.path.isfile(path), f"scripts/{script} not found"


class TestAllModuleConfigs:
    """Test that all v3.0 module configs are present."""

    def test_module_toggles_exist(self, main_config):
        """All module toggle flags must be present."""
        params = main_config.get("params", {})
        expected_toggles = [
            "run_lncrna", "run_smallrna", "run_mirna_targets",
            "run_epitx_enhanced", "run_integration"
        ]
        for toggle in expected_toggles:
            assert toggle in params, f"Missing params.{toggle}"
            assert isinstance(params[toggle], bool), (
                f"params.{toggle} must be boolean, got {type(params[toggle])}"
            )

    def test_annotations_section_exists(self, main_config):
        """Annotations section must be present in config."""
        assert "annotations" in main_config, "Missing 'annotations' section in config"
        annotations = main_config["annotations"]
        expected_keys = [
            "araport11_gff", "cantatadb_bed", "mirbase_gff",
            "mirbase_mature_fa", "pmiren_fa", "te_annotation"
        ]
        for key in expected_keys:
            assert key in annotations, f"Missing annotations.{key}"

    def test_srna_samples_section_exists(self, main_config):
        """srna_samples section must be present (even if empty)."""
        assert "srna_samples" in main_config, "Missing 'srna_samples' section"

    def test_epitx_comparisons_exist(self, main_config):
        """epitx_comparisons must be defined."""
        assert "epitx_comparisons" in main_config, "Missing 'epitx_comparisons' section"
        comps = main_config["epitx_comparisons"]
        assert isinstance(comps, list), "epitx_comparisons must be a list"
        if comps:
            assert "name" in comps[0], "Each comparison needs a 'name'"
            assert "control" in comps[0], "Each comparison needs 'control' samples"
            assert "treated" in comps[0], "Each comparison needs 'treated' samples"
