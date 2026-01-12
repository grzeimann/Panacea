# Ensure the package in src/ is importable when running tests without installation
from __future__ import annotations

import sys
from pathlib import Path
import pytest

# Prepend the repository's src directory to sys.path
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

# Provide a small synthetic dataset for tests
from tests.fixtures.sample_data import write_sample_dataset  # noqa: E402


@pytest.fixture(scope="session")
def sample_dataset(tmp_path_factory):
    """Create a tiny synthetic dataset once per test session and return a dict of paths.

    Returns a mapping with keys: bias, flat, arc, twilight, science.
    """
    base = tmp_path_factory.mktemp("sample_data")
    paths = write_sample_dataset(base)
    return paths
