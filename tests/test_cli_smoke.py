from __future__ import annotations

import sys
import pytest


def test_cli_help_exits_zero(capsys, monkeypatch):
    # Simulate calling the console script with -h and ensure it exits cleanly
    monkeypatch.setenv("PYTHONWARNINGS", "ignore")
    monkeypatch.setenv("NO_COLOR", "1")
    monkeypatch.setenv("TERM", "xterm")
    monkeypatch.setenv("COLUMNS", "100")

    monkeypatch.setattr(sys, "argv", ["panacea-lrs2", "-h"])  # show help

    from panacea.cli import main

    with pytest.raises(SystemExit) as exc:
        main()

    # Help should exit with code 0
    assert exc.value.code == 0

    out, err = capsys.readouterr()
    # Expect some help text mentioning options or usage
    assert ("-d" in out or "--date" in out or "usage" in out.lower())
