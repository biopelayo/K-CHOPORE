"""
Test 2: Minimap2 flags validation
Ensures that BOTH the Snakefile AND the legacy test.sh use the correct
minimap2 flags for ONT direct RNA-seq alignment.

Critical flags:
  -ax splice   (splice-aware alignment for RNA)
  -uf          (forward strand for direct RNA reads)
  -k14         (kmer size recommended for RNA)
  --MD          (MD tag needed for modification detection)
  --secondary=no  (no secondary alignments)

WRONG:  -ax map-ont  (DNA preset, NOT suitable for direct RNA)
"""
import re
import pytest


class TestSnakefileMinimap2:
    """Validate minimap2 flags in the main Snakefile."""

    def test_snakefile_uses_splice_preset(self, snakefile_content):
        """Snakefile must use -ax splice, not -ax map-ont."""
        assert "-ax {params.preset}" in snakefile_content or "-ax splice" in snakefile_content, (
            "Snakefile minimap2 rule should use splice preset"
        )
        # Ensure map-ont is NOT used anywhere in the Snakefile
        # (except maybe comments)
        lines = snakefile_content.split("\n")
        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            assert "map-ont" not in stripped, (
                f"Snakefile line {i} contains 'map-ont' which is wrong for direct RNA: {stripped}"
            )

    def test_snakefile_uses_uf_flag(self, snakefile_content):
        """Snakefile must use -uf for forward-stranded direct RNA reads."""
        assert "-uf" in snakefile_content, (
            "Snakefile minimap2 rule must include '-uf' for direct RNA"
        )

    def test_snakefile_has_junc_bed(self, snakefile_content):
        """Snakefile should use --junc-bed for splice junction guidance."""
        assert "--junc-bed" in snakefile_content, (
            "Snakefile minimap2 rule should include '--junc-bed' for splice-aware alignment"
        )


class TestLegacyTestSh:
    """Validate minimap2 flags in the legacy test.sh script."""

    def test_no_map_ont_in_test_sh(self, legacy_test_sh_content):
        """Legacy test.sh must NOT use -ax map-ont."""
        lines = legacy_test_sh_content.split("\n")
        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            assert "map-ont" not in stripped, (
                f"test.sh line {i} still uses 'map-ont': {stripped}\n"
                "This is WRONG for direct RNA-seq. Must use '-ax splice -uf -k14'"
            )

    def _get_minimap2_command_lines(self, content):
        """Extract lines that actually invoke minimap2 binary (not echo/debug lines)."""
        lines = content.split("\n")
        cmd_lines = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            if stripped.startswith("echo "):
                continue
            # Only lines that call the minimap2 binary (contain path or direct call)
            if "minimap2 " in stripped and "-a" in stripped:
                cmd_lines.append(stripped)
        return cmd_lines

    def test_uses_splice_preset(self, legacy_test_sh_content):
        """Legacy test.sh must use -ax splice."""
        cmd_lines = self._get_minimap2_command_lines(legacy_test_sh_content)
        assert len(cmd_lines) > 0, "No minimap2 command invocations found in test.sh"
        for line in cmd_lines:
            assert "-ax splice" in line, (
                f"minimap2 command should use '-ax splice', found: {line}"
            )

    def test_uses_uf_flag(self, legacy_test_sh_content):
        """Legacy test.sh must use -uf for forward-stranded reads."""
        cmd_lines = self._get_minimap2_command_lines(legacy_test_sh_content)
        for line in cmd_lines:
            assert "-uf" in line, (
                f"minimap2 command should include '-uf', found: {line}"
            )

    def test_uses_k14(self, legacy_test_sh_content):
        """Legacy test.sh must use -k14 for RNA reads."""
        cmd_lines = self._get_minimap2_command_lines(legacy_test_sh_content)
        for line in cmd_lines:
            assert "-k14" in line or "-k 14" in line, (
                f"minimap2 command should include '-k14', found: {line}"
            )

    def test_uses_MD_tag(self, legacy_test_sh_content):
        """Legacy test.sh must use --MD for modification detection."""
        cmd_lines = self._get_minimap2_command_lines(legacy_test_sh_content)
        for line in cmd_lines:
            assert "--MD" in line, (
                f"minimap2 command should include '--MD', found: {line}"
            )
