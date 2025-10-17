"""Package CLI entry point for Panacea.

Decoupled from the legacy full_lrs2_reduction script. This delegates to the
package's modular CLI (panacea.run_panacea) so that help and argument parsing
work without the monolithic module.
"""
from __future__ import annotations

import sys


def main():
    """Entry point used by the console script.

    Delegates to the modular CLI implemented in ``panacea.run_panacea`` so that
    ``-h`` and argument parsing do not depend on the legacy script slated for
    removal.
    """
    from . import run_panacea as runner

    # This will print help and raise SystemExit on "-h" as expected by tests.
    runner.main(sys.argv[1:])
