from __future__ import annotations

import sys
import pytest


def test_pipeline_smoke_resources_and_wiring(capsys, monkeypatch):
    """Run the CLI with --smoke-test to verify packaged resources and wiring.

    This does not access raw data or perform heavy I/O. It ensures that the
    console script is callable and that required packaged configuration files
    are discoverable for all default sides.
    """
    # Simulate calling the console script with the hidden smoke-test flag
    monkeypatch.setenv("PYTHONWARNINGS", "ignore")
    monkeypatch.setenv("NO_COLOR", "1")
    monkeypatch.setenv("TERM", "xterm")
    monkeypatch.setenv("COLUMNS", "100")

    # Default sides in the CLI are uv,orange,red,farred; rely on defaults here
    monkeypatch.setattr(
        sys, "argv", ["panacea-lrs2", "--smoke-test"]
    )  # early exit path

    from panacea.cli import main

    with pytest.raises(SystemExit) as exc:
        main()

    # Smoke test should exit with code 0
    assert exc.value.code == 0

    out, err = capsys.readouterr()
    assert "smoke test" in out.lower()
