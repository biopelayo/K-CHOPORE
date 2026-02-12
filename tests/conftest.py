"""
Pytest fixtures for K-CHOPORE test suite.
Provides shared paths, config loading, and toy data references.
"""
import os
import pytest
import yaml

# ----------------------------------------------------------------
# Paths
# ----------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
TOY_DATA_DIR = os.path.join(TESTS_DIR, "toy_data")


@pytest.fixture
def project_root():
    """Return absolute path to the project root."""
    return PROJECT_ROOT


@pytest.fixture
def tests_dir():
    """Return absolute path to the tests directory."""
    return TESTS_DIR


@pytest.fixture
def toy_data_dir():
    """Return absolute path to toy data directory."""
    return TOY_DATA_DIR


@pytest.fixture
def main_config_path():
    """Return path to the main production config."""
    return os.path.join(PROJECT_ROOT, "config", "config.yml")


@pytest.fixture
def test_config_path():
    """Return path to the test config."""
    return os.path.join(TESTS_DIR, "test_config.yml")


@pytest.fixture
def main_config(main_config_path):
    """Load and return the main production config as a dict."""
    with open(main_config_path, "r") as f:
        return yaml.safe_load(f)


@pytest.fixture
def test_config(test_config_path):
    """Load and return the test config as a dict."""
    with open(test_config_path, "r") as f:
        return yaml.safe_load(f)


@pytest.fixture
def snakefile_path():
    """Return path to the main Snakefile."""
    return os.path.join(PROJECT_ROOT, "Snakefile")


@pytest.fixture
def snakefile_content(snakefile_path):
    """Return the full Snakefile content as a string."""
    with open(snakefile_path, "r") as f:
        return f.read()


@pytest.fixture
def legacy_test_sh_path():
    """Return path to the legacy test.sh script."""
    return os.path.join(PROJECT_ROOT, "scripts", "pipelines", "test.sh")


@pytest.fixture
def legacy_test_sh_content(legacy_test_sh_path):
    """Return the full legacy test.sh content as a string."""
    with open(legacy_test_sh_path, "r") as f:
        return f.read()
