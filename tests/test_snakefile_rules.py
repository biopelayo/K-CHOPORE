"""
Test 4: Snakefile rule content validation
Static analysis of the Snakefile to check for common issues
without actually running the pipeline.

Updated for K-CHOPORE v3.0: adapted to v2/v3 Snakefile architecture.
"""
import re
import pytest


class TestRuleConsistency:
    """Check that Snakefile rules follow expected patterns."""

    def test_all_rules_have_log_directive(self, snakefile_content):
        """Every rule with shell: (except 'all' and 'run:') should have a log directive."""
        rules = re.findall(r'^rule\s+(\w+):', snakefile_content, re.MULTILINE)
        skip_rules = {"all", "setup_complete_structure", "generate_reads_manifest",
                       "generate_deseq2_sample_sheet"}

        for rule_name in rules:
            if rule_name in skip_rules:
                continue
            # Find rule block
            pattern = rf'rule\s+{rule_name}:.*?(?=\nrule\s|\n#\s*#{{5,}}|\Z)'
            match = re.search(pattern, snakefile_content, re.DOTALL)
            assert match, f"Could not find rule block for '{rule_name}'"
            block = match.group(0)
            # Only check shell rules (not run: rules)
            if "shell:" in block:
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
            pattern = rf'rule\s+{rule_name}:.*?(?=\nrule\s|\n#\s*#{{5,}}|\Z)'
            match = re.search(pattern, snakefile_content, re.DOTALL)
            if match and 'shell:' in match.group(0):
                assert "[K-CHOPORE]" in match.group(0), (
                    f"Rule '{rule_name}' shell block should contain [K-CHOPORE] logging"
                )

    def test_no_hardcoded_thread_counts(self, snakefile_content):
        """Shell blocks should not hardcode thread counts; use {threads}."""
        shell_blocks = re.findall(
            r'shell:\s*"""(.*?)"""', snakefile_content, re.DOTALL
        )
        for i, block in enumerate(shell_blocks):
            hardcoded = re.findall(r'-t\s+\d+(?!\s*\{)', block)
            assert len(hardcoded) == 0, (
                f"Shell block {i+1} has hardcoded thread count: {hardcoded}. "
                "Use {{threads}} instead."
            )

    def test_rule_all_has_conditional_targets(self, snakefile_content):
        """Rule 'all' should use conditional targets for v3.0 modules."""
        # In v2/v3, rule all directly lists targets with conditionals
        assert "rule all:" in snakefile_content, "Rule 'all' not found"
        assert "if RUN_LNCRNA" in snakefile_content, (
            "Rule 'all' should conditionally include lncRNA targets"
        )
        assert "if RUN_M6ANET" in snakefile_content or "RUN_M6ANET" in snakefile_content, (
            "Snakefile should reference RUN_M6ANET toggle"
        )

    def test_minimap2_rule_uses_params(self, snakefile_content):
        """map_with_minimap2 should use params, not hardcoded values."""
        pattern = r'rule\s+map_with_minimap2:.*?(?=\nrule\s|\n#\s*#{{5,}}|\Z)'
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

    def test_minimap2_uses_uf_flag(self, snakefile_content):
        """Critical: minimap2 for direct RNA must use -uf (forward strand)."""
        pattern = r'rule\s+map_with_minimap2:.*?(?=\nrule\s|\n#\s*#{{5,}}|\Z)'
        match = re.search(pattern, snakefile_content, re.DOTALL)
        assert match, "Could not find rule 'map_with_minimap2'"
        block = match.group(0)
        assert "-uf" in block, (
            "map_with_minimap2 MUST use -uf flag for direct RNA-seq (forward strand)"
        )


class TestConditionalExecution:
    """Test that v3.0 conditional module toggles are defined."""

    def test_module_toggle_variables_defined(self, snakefile_content):
        """All v3.0 module toggle variables must be defined."""
        expected_vars = [
            "RUN_LNCRNA", "RUN_SMALLRNA", "RUN_MIRNA_TARGETS",
            "RUN_EPITX_ENHANCED", "RUN_INTEGRATION"
        ]
        for var in expected_vars:
            assert var in snakefile_content, (
                f"Module toggle variable '{var}' not found in Snakefile"
            )

    def test_conditional_targets_in_rule_all(self, snakefile_content):
        """rule all should include conditional targets for each module."""
        # Check that conditional patterns exist
        assert 'if RUN_LNCRNA else []' in snakefile_content
        assert 'if RUN_SMALLRNA else []' in snakefile_content
        assert 'if RUN_MIRNA_TARGETS else []' in snakefile_content
        assert 'if RUN_EPITX_ENHANCED else []' in snakefile_content
        assert 'if RUN_INTEGRATION else []' in snakefile_content

    def test_srna_sample_parsing(self, snakefile_content):
        """Small RNA sample list must be parsed from config."""
        assert "SRNA_SAMPLES" in snakefile_content
        assert "SRNA_SAMPLE_IDS" in snakefile_content

    def test_annotation_paths_parsed(self, snakefile_content):
        """Annotation paths must be parsed from config."""
        assert "ANNOTATIONS" in snakefile_content


class TestV3RuleCount:
    """Test that all expected v3.0 rules exist."""

    def test_module1_rule_count(self, snakefile_content):
        """Module 1 should have 7 lncRNA rules."""
        lncrna_rules = re.findall(r'^rule\s+lncrna_\w+:', snakefile_content, re.MULTILINE)
        assert len(lncrna_rules) == 7, (
            f"Expected 7 lncRNA rules, found {len(lncrna_rules)}: {lncrna_rules}"
        )

    def test_module2_rule_count(self, snakefile_content):
        """Module 2 should have 5 sRNA rules."""
        srna_rules = re.findall(r'^rule\s+srna_\w+:', snakefile_content, re.MULTILINE)
        assert len(srna_rules) == 5, (
            f"Expected 5 sRNA rules, found {len(srna_rules)}: {srna_rules}"
        )

    def test_module3_rule_count(self, snakefile_content):
        """Module 3 should have 4 target rules."""
        target_rules = re.findall(r'^rule\s+(targets_|degradome_)\w+:', snakefile_content, re.MULTILINE)
        assert len(target_rules) == 4, (
            f"Expected 4 target rules, found {len(target_rules)}: {target_rules}"
        )

    def test_module4_has_epitx_rules(self, snakefile_content):
        """Module 4 should have epitx consensus and report rules."""
        assert "rule epitx_consensus:" in snakefile_content
        assert "rule epitx_report:" in snakefile_content

    def test_module5_rule_count(self, snakefile_content):
        """Module 5 should have 4 integration rules."""
        integration_rules = [
            "wgcna_network", "cis_regulation",
            "population_context", "integration_report"
        ]
        for rule in integration_rules:
            assert f"rule {rule}:" in snakefile_content, (
                f"Missing Module 5 rule: {rule}"
            )
