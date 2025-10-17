import pytest

import full_lrs2_reduction as panacea


def test_imports():
    # Basic import smoke test
    assert hasattr(panacea, "parser")


def test_cli_help_exit():
    # argparse -h should trigger SystemExit (exit code 0)
    with pytest.raises(SystemExit) as exc:
        panacea.parser.parse_args(["-h"])  # noqa: F841
    assert exc.value.code == 0
