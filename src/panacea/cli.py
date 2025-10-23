"""Package entry point for Panacea.

Decoupled from the legacy full_lrs2_reduction script so that help and argument parsing
work without the previous monolithic module.
"""
from __future__ import annotations

import sys


def main():
    """
    Entry point used by the console script.
    """
    from . import run_panacea as runner

    # This will print help and raise SystemExit on "-h" as expected by tests.
    # Delegate without passing argv; runner.main() internally reads sys.argv.
    runner.main()
