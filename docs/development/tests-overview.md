# Tests overview (what we cover and where)

This page explains where the project’s tests live and what the most important ones validate. It’s intended to help users, reviewers (e.g., JOSS), and new contributors understand how to run and interpret the suite.

Where tests live
- All tests are under the repository’s tests/ directory.
- They are designed to be fast and data‑free except where noted; you can run them locally with:
  - make test  (or)  pytest -q
  - For coverage: coverage run -m pytest -q && coverage report -m

Key tests and what they check
- CLI and wiring
  - tests/test_cli.py — console entry panacea-lrs2 exists; -h exits 0
  - tests/test_cli_smoke.py — simulates panacea-lrs2 -h and checks help text
  - tests/test_pipeline_smoke.py — uses the hidden --smoke-test path to verify that packaged resources are discoverable for all default sides without touching raw data
- Parsers and resources
  - tests/test_parsers.py — validates simple parsing against packaged config files (DAR tables, skylines) and the read_arc_lines helper returns an Astropy Table with expected columns
  - tests/test_resources.py — ensures packaged resource files are found via importlib.resources and that read_arc_lines parses minimal valid content
- Calibration and utilities
  - tests/test_calibration_utils.py — unit tests for safe_division, build_weight_matrix, robust_polyfit, and response‑continuum estimation
- Sample data sanity (for reviewers)
  - tests/test_sample_data.py — uses a tiny synthetic dataset created on the fly (via tests/fixtures/) to assert that expected FITS files (bias, flat, arc, twilight, science) exist with small shapes and basic headers (INSTRUME=LRS2, CHANNEL present, EXPTIME set). It also checks that the arc has columns significantly brighter than median (a simple spectral‑line heuristic).

How this relates to the sample dataset (~400 MB)
- For a realistic end‑to‑end run, see Quickstart → “Run on sample data (~400 MB)”. Those instructions are meant for users or reviewers who want to try a real reduction quickly.
- The tests/test_sample_data.py file is intentionally small and synthetic to keep CI fast; it is not the same as the ~400 MB dataset used in the Quickstart run. The two complement each other:
  - tests/test_sample_data.py proves readers and I/O paths work on representative FITS structures
  - The Quickstart sample dataset demonstrates full pipeline behavior on real data

How to run only selected tests
- Single file: pytest tests/test_parsers.py -q
- Single test by name: pytest -k test_read_arc_lines_stringio_basic -q
- Stop on first failure (handy when debugging): pytest -x -q

See also
- Developer quickstart (env, hooks, CI parity): ../development/dev-quickstart.md
- Quickstart sample dataset and run: ../getting-started/quickstart.md
