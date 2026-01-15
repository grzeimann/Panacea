# Panacea: Automatic LRS2 data-reduction pipeline for HET

[![CI](https://github.com/grzeimann/Panacea/actions/workflows/python-tests.yml/badge.svg)](https://github.com/grzeimann/Panacea/actions/workflows/python-tests.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18250411.svg)](https://doi.org/10.5281/zenodo.18250411)

Panacea reduces LRS2 raw CCD frames into science‑ready spectra and data cubes for the Hobby–Eberly Telescope (HET). It powers the daily automated reductions on TACC and can be installed locally for development or small analyses. For a structured overview, see [Project documentation index](docs/index.md).

- Getting started: [Installation](docs/getting-started/installation.md) · [Quickstart](docs/getting-started/quickstart.md)
- TACC users: [Overview](docs/tacc/overview.md) · [Running on TACC](docs/tacc/running.md)
- Data products: [Overview](docs/data-products/overview.md)

## Install (conda + pip)

```bash
conda env create -f environment.yml
conda activate panacea
pip install .[dev]
```

Compatibility note: The environment is pinned to numpy<2.0 with astropy<6 (see [environment.yml](environment.yml) / [pyproject.toml](pyproject.toml)). If you need NumPy 2.x, upgrade Astropy to ≥6 and test locally.

## Quickstart (CLI)

```bash
# Show help and verify installation
panacea-lrs2 -h

# Smoke test: verify packaged resources without raw data
panacea-lrs2 --smoke-test

# Example run (adjust to your data)
panacea-lrs2 -d 20181108 -s uv
```

More examples and options: [CLI user guide](docs/user-guide/cli.md) and [Examples](docs/user-guide/examples.md).

## Data and configuration
- Required configuration files (e.g., line lists, DAR tables, responses, fplane.txt) are bundled with the package and resolved via importlib.resources.
- Raw data layout: the CLI expects an LRS2/TACC‑like directory structure. Point to your data root using `--baseraw`.

## TACC pipeline users
If you use the facility pipeline on TACC, start with [TACC overview](docs/tacc/overview.md).

## Troubleshooting
- Import error mentioning `numpy.in1d` or Astropy: ensure your environment follows the pins above.
- Need dev commands (pre‑commit, local CI, coverage)? See [Dev quickstart](docs/development/dev-quickstart.md). For exact git workflow commands (stage/commit/push), see [Stage, commit, push](docs/development/dev-quickstart.md#stage-commit-push-exact-commands).
- More help: [FAQ](docs/faq.md).

## Citing and License
- Citation: [CITATION.cff](CITATION.cff) (see also [Citation docs](docs/citation/citation.md))
- DOI (concept, stable): https://doi.org/10.5281/zenodo.18250411
- DOI (version v1.0.1): https://doi.org/10.5281/zenodo.18250412
- License: BSD‑3‑Clause ([LICENSE](LICENSE))

## Contributing
We welcome issues and pull requests. Please read [Contributing guide](docs/community/contributing.md). For a concise developer setup and CI‑equivalent checks, see [Dev quickstart](docs/development/dev-quickstart.md).
