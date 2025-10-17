"""Package CLI delegator for Panacea.

Currently delegates to the legacy full_lrs2_reduction parser to preserve
behavior while refactoring proceeds.
"""
from __future__ import annotations

import sys


def main():
    """Entry point used by the console script.

    For now, simply invoke the legacy argparse parser so that ``-h`` and other
    options behave identically to the monolithic script.
    """
    import full_lrs2_reduction as legacy

    # This will print help and raise SystemExit on "-h" as expected by tests.
    legacy.parser.parse_args(sys.argv[1:])
