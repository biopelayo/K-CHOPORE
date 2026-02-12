"""
Test 4: Snakefile rule content validation
Static analysis of the Snakefile to check for common issues
without actually running the pipeline.
"""
import re
import pytest


class TestRuleConsistency:
    """Check that Snakefile rules follow expected patterns."""

    def test_all_rules_have_log_directive(self, snakefile_content):
        """Every rule (except 'all' and 'setup') should have a log directive."""
        # Extract rule names
        rules = re.findall(r'^rule\s+(\w+):', snakefile_content, re.MULTILINE)
        skip_rules = {"all", "setup_complete_structure"}

        for rule_name in rules:
            if rule_name in skip_rules:
                continue
            # Find the rule block and check for 'log:'
            pattern = rf'rule\s+{rule_name}:.*?(?=\nrule\s|\Z)'
            match = re.search(pattern, snakefile_content, re.DOTALL)
            assert match, f"Could not find rule block for '{rule_name}'"
            block = match.group(0)
            assert "log:" in block, (
                f"Rule '{rule_name}' is missing a 'log:' directive"
            )

    def test_all_shell_rules_have_echo(self, snakefile_content):
        """Every shell rule should have K-CHOPORE echo messages for logging."""
        rules = re.findall(r'^rule\s+(\w+):', snakefile_content, re.MULTILINE)
        skip_rules = {"all", "setup_complete_structure"}

        for rule_name in rules:
            if rule_name in skip_rules:
                continue
            pattern = rf'rule\s+{rule_name}:.*?(?=\nrule\s|\Z)'
            match = re.search(pattern, snakefile_content, re.DOTALL)
            if match and 'shell:' in match.group(0):
                assert "[K-CHOPORE]" in match.group(0), (
                    f"Rule '{rule_name}' shell block should contain [K-CHOPORE] logging"
                )

    def test_no_hardcoded_thread_counts(self, snakefile_content):
        """Shell blocks should not hardcode thread counts; use {threads}."""
        # Find all shell blocks
        shell_blocks = re.findall(
            r'shell:\s*"""(.*?)"""', snakefile_content, re.DOTALL
        )
        for i, block in enumerate(shell_blocks):
            # Check for -t followed by a hardcoded number (not {threads})
            hardcoded = re.findall(r'-t\s+\d+(?!\s*\{)', block)
            assert len(hardcoded) == 0, (
                f"Shell block {i+1} has hardcoded thread count: {hardcoded}. "
                "Use {{threads}} instead."
            )

    def test_rule_all_calls_get_all_targets(self, snakefile_content):
        """Rule 'all' should use the get_all_targets() function."""
        assert "get_all_targets()" in snakefile_content, (
            "Rule 'all' should use get_all_targets() for dynamic target generation"
        )

    def test_minimap2_rule_uses_params(self, snakefile_content):
        """map_with_minimap2 should use params, not hardcoded values."""
        pattern = r'rule\s+map_with_minimap2:.*?(?=\nrule\s|\Z)'
        match = re.search(pattern, snakefile_content, re.DOTALL)
        assert match, "Could not find rule 'map_with_minimap2'"
        block = match.group(0)
        assert "{params.preset}" in block, (
            "map_with_minimap2 should use {params.preset} from config"
        )
        assert "{params.kmer}" in block, (
            "map_with_minimap2 should use {params.kmer} from config"
        )
        assert "{params.extra}" in block, (
            "map_with_minimap2 should use {params.extra} from config"
        )

    def test_pycoqc_uses_lambda_for_summary(self, snakefile_content):
        """pycoQC rule should use lambda to resolve summary paths (not hardcoded pattern)."""
        pattern = r'rule\s+quality_analysis_with_pycoQC:.*?(?=\nrule\s|\Z)'
        match = re.search(pattern, snakefile_content, re.DOTALL)
        assert match, "Could not find rule 'quality_analysis_with_pycoQC'"
        block = match.group(0)
        assert "lambda wildcards" in block, (
            "pycoQC rule should use 'lambda wildcards: SEQUENCING_SUMMARY[wildcards.sample]' "
            "instead of a hardcoded filename pattern"
        )
        # Make sure the old hardcoded pattern is NOT there
        assert "{sample}_sequencing_summary.txt" not in block, (
            "pycoQC rule still has hardcoded '{sample}_sequencing_summary.txt' pattern. "
            "Should use lambda for flexible filename resolution."
        )


class TestConditionalExecution:
    """Test that the conditional module toggles work correctly."""

    def test_get_all_targets_respects_run_flags(self):
        """Verify the get_all_targets function uses config run_* flags."""
        # We test this by reading the Snakefile and checking the function
        from tests.conftest import PROJECT_ROOT
        import os
        with open(os.path.join(PROJECT_ROOT, "Snakefile")) as f:
            content = f.read()

        # Extract the get_all_targets function
        pattern = r'def get_all_targets\(\):.*?return targets'
        match = re.search(pattern, content, re.DOTALL)
        assert match, "Could not find get_all_targets() function"
        func_body = match.group(0)

        # Check that it references the key run flags
        expected_checks = [
            "run_basecalling",
            "run_nanofilt",
            "run_nanoplot",
            "run_nanocomp",
            "run_pycoqc",
            "run_flair",
            "run_eligos2",
            "run_m6anet",
            "run_deseq2",
            "run_multiqc",
        ]
        for flag in expected_checks:
            assert flag in func_body, (
                f"get_all_targets() should check '{flag}' config flag"
            )
