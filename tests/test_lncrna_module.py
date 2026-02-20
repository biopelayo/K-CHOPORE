"""
Test Suite: Module 1 — lncRNA Discovery
Tests for the lncRNA consensus classification pipeline.
"""
import os
import pytest
import sys

# Add scripts to path for import testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))


class TestLncRNAConsensusScript:
    """Test the lncrna_consensus.py helper functions."""

    def test_script_exists(self, project_root):
        """lncrna_consensus.py must exist."""
        script = os.path.join(project_root, "scripts", "lncrna_consensus.py")
        assert os.path.isfile(script), "scripts/lncrna_consensus.py not found"

    def test_script_imports(self):
        """Script must be importable (no syntax errors)."""
        import importlib.util
        script_path = os.path.join(
            os.path.dirname(__file__), "..", "scripts", "lncrna_consensus.py"
        )
        spec = importlib.util.spec_from_file_location("lncrna_consensus", script_path)
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
        except SystemExit:
            pass  # argparse may call sys.exit when no args provided

    def test_parse_transdecoder_empty(self, tmp_path):
        """TransDecoder parser handles empty file gracefully."""
        from lncrna_consensus import parse_transdecoder
        empty_file = tmp_path / "empty.pep"
        empty_file.write_text("")
        result = parse_transdecoder(str(empty_file), 100)
        assert result == set()

    def test_parse_transdecoder_with_orfs(self, tmp_path):
        """TransDecoder parser correctly identifies long ORFs."""
        from lncrna_consensus import parse_transdecoder
        pep_file = tmp_path / "test.pep"
        # Create a protein with >100 aa
        long_protein = "M" + "A" * 150
        short_protein = "M" + "A" * 50
        pep_file.write_text(
            f">transcript1.p1 type:complete\n{long_protein}\n"
            f">transcript2.p1 type:complete\n{short_protein}\n"
        )
        result = parse_transdecoder(str(pep_file), 100)
        assert "transcript1" in result
        assert "transcript2" not in result

    def test_parse_cpc2_noncoding(self, tmp_path):
        """CPC2 parser correctly identifies noncoding transcripts."""
        from lncrna_consensus import parse_cpc2
        cpc2_file = tmp_path / "cpc2.txt"
        cpc2_file.write_text(
            "ID\ttranscript_length\tpeptide_length\tFickett_score\tpI\tORF_integrity\tcoding_probability\tlabel\n"
            "tx1\t500\t0\t0.5\t7.0\tno\t0.1\tnoncoding\n"
            "tx2\t1000\t200\t0.9\t6.0\tyes\t0.95\tcoding\n"
            "tx3\t800\t50\t0.3\t8.0\tno\t0.05\tnoncoding\n"
        )
        result = parse_cpc2(str(cpc2_file))
        assert "tx1" in result
        assert "tx2" not in result
        assert "tx3" in result

    def test_parse_cpat(self, tmp_path):
        """CPAT parser correctly identifies noncoding (prob < threshold)."""
        from lncrna_consensus import parse_cpat
        cpat_file = tmp_path / "cpat.tsv"
        cpat_file.write_text(
            "seq_ID\tmRNA_size\tORF_size\tFickett_score\tHexamer_score\tcoding_prob\n"
            "tx1\t500\t0\t0.5\t-0.3\t0.1\n"
            "tx2\t1000\t600\t0.9\t0.8\t0.95\n"
            "tx3\t800\t100\t0.3\t-0.5\t0.2\n"
        )
        result = parse_cpat(str(cpat_file), 0.39)
        assert "tx1" in result
        assert "tx2" not in result
        assert "tx3" in result


class TestLncRNASnakefileRules:
    """Test that lncRNA rules are properly defined in the Snakefile."""

    def test_lncrna_rules_present(self, snakefile_content):
        """All 7 lncRNA rules must be defined."""
        expected_rules = [
            "lncrna_filter_candidates",
            "lncrna_transdecoder",
            "lncrna_feelnc",
            "lncrna_cpc2",
            "lncrna_cpat",
            "lncrna_consensus",
            "lncrna_deseq2"
        ]
        for rule in expected_rules:
            assert f"rule {rule}:" in snakefile_content, (
                f"Rule '{rule}' not found in Snakefile"
            )

    def test_lncrna_toggle_exists(self, snakefile_content):
        """RUN_LNCRNA toggle variable must exist."""
        assert "RUN_LNCRNA" in snakefile_content

    def test_lncrna_conditional_targets(self, snakefile_content):
        """rule all must include conditional lncRNA targets."""
        assert 'if RUN_LNCRNA' in snakefile_content
        assert 'lncrna_final.gtf' in snakefile_content


class TestLncRNAConfig:
    """Test lncRNA config parameters."""

    def test_lncrna_params_in_main_config(self, main_config):
        """Main config must have lncRNA parameters."""
        params = main_config.get("params", {})
        assert "run_lncrna" in params, "Missing params.run_lncrna"
        assert "lncrna_min_length" in params, "Missing params.lncrna_min_length"
        assert "lncrna_consensus_min" in params, "Missing params.lncrna_consensus_min"
        assert "cpat_threshold" in params, "Missing params.cpat_threshold"

    def test_lncrna_min_length_valid(self, main_config):
        """lncRNA minimum length must be >= 200."""
        min_len = main_config["params"]["lncrna_min_length"]
        assert min_len >= 200, f"lncrna_min_length should be >= 200, got {min_len}"

    def test_cpat_threshold_range(self, main_config):
        """CPAT threshold must be between 0 and 1."""
        threshold = main_config["params"]["cpat_threshold"]
        assert 0 < threshold < 1, f"cpat_threshold must be (0,1), got {threshold}"
