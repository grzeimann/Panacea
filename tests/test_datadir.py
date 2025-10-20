from __future__ import annotations

import os
from pathlib import Path


def _guess_baseraw() -> Path:
    """Return a plausible base raw-data directory for tests.

    Priority:
    - PANACEA_BASERAW environment variable if set.
    - Repository-local LRS2/ directory (contains example layout in this repo).
    """
    env = os.getenv("PANACEA_BASERAW")
    if env:
        return Path(env).expanduser().resolve()
    repo_root = Path(__file__).resolve().parents[1]
    return repo_root / "LRS2"


def test_baseraw_has_night_data():
    """Ensure there is at least one night's worth of data available.

    This is a lightweight presence test. It checks that the base raw-data
    directory exists and has at least one subdirectory (date/program folder)
    or at least one tarball somewhere beneath it.

    You can override the default base directory by setting PANACEA_BASERAW.
    """
    baseraw = _guess_baseraw()
    assert baseraw.exists() and baseraw.is_dir(), f"Base raw directory not found: {baseraw}"

    # Consider either subdirectories (e.g., YYYYMMDD or program folders)
    # or tarballs as evidence of data presence.
    subdirs = [p for p in baseraw.iterdir() if p.is_dir()]
    tarballs = list(baseraw.rglob("*.tar"))

    assert (len(subdirs) > 0) or (len(tarballs) > 0), (
        f"No data found under {baseraw}. Expected at least one subdirectory or .tar file."
    )
