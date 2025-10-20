# Ensure the package in src/ is importable when running tests without installation
from __future__ import annotations

import sys
from pathlib import Path

# Prepend the repository's src directory to sys.path
_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
