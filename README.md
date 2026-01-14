# Panacea: Automatic LRS2 data-reduction pipeline for HET

[![CI](https://github.com/grzeimann/Panacea/actions/workflows/python-tests.yml/badge.svg)](https://github.com/grzeimann/Panacea/actions/workflows/python-tests.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

Panacea reduces LRS2 raw CCD frames into science‑ready spectra and data cubes for the Hobby–Eberly Telescope (HET). It powers the daily automated reductions on TACC and can be installed locally for development or small analyses. For a structured overview, see docs/index.md.

- Getting started: docs/getting-started/installation.md · docs/getting-started/quickstart.md
- TACC users: docs/tacc/overview.md · docs/tacc/running.md
- Data products: docs/data-products/overview.md

## Install (conda + pip)

```bash
conda env create -f environment.yml
conda activate panacea
pip install .[dev]
```

Compatibility note: The environment is pinned to numpy<2.0 with astropy<6 (see environment.yml / pyproject.toml). If you need NumPy 2.x, upgrade Astropy to ≥6 and test locally.

## Quickstart (CLI)

```bash
# Show help and verify installation
panacea-lrs2 -h

# Smoke test: verify packaged resources without raw data
panacea-lrs2 --smoke-test

# Example run (adjust to your data)
panacea-lrs2 -d 20181108 -s uv
```

More examples and options: docs/user-guide/cli.md and docs/user-guide/examples.md.

## Data and configuration
- Required configuration files (e.g., line lists, DAR tables, responses, fplane.txt) are bundled with the package and resolved via importlib.resources.
- Raw data layout: the CLI expects an LRS2/TACC‑like directory structure. Point to your data root using `--baseraw`.

## TACC pipeline users
If you use the facility pipeline on TACC, start with docs/tacc/overview.md.

## Troubleshooting
- Import error mentioning `numpy.in1d` or Astropy: ensure your environment follows the pins above.
- Need dev commands (pre‑commit, local CI, coverage)? See docs/development/dev-quickstart.md.
- More help: docs/faq.md.

## Citing and License
- Citation: CITATION.cff (see also docs/citation/citation.md)
- License: BSD‑3‑Clause (LICENSE)

## Contributing
We welcome issues and pull requests. Please read docs/community/contributing.md. For a concise developer setup and CI‑equivalent checks, see docs/development/dev-quickstart.md.
