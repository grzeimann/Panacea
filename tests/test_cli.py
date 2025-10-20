import sys
import pytest


def test_cli_main_exists():
    import panacea.cli as cli
    assert callable(cli.main)


def test_cli_help_exit_zero(monkeypatch):
    import panacea.cli as cli
    # Simulate `panacea-lrs2 -h`
    monkeypatch.setattr(sys, "argv", ["panacea-lrs2", "-h"])  # type: ignore[attr-defined]
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0
