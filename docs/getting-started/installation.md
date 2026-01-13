# Installation

This page covers local installation and dependencies. For TACC usage, see ../tacc/overview.md.

## Dependencies
- Python 3.9+
- NumPy, SciPy, Astropy, Matplotlib, PyYAML, tqdm, requests, scikit-learn (installed automatically from PyPI)

## Create a conda environment and install
```bash
conda env create -f environment.yml  # run from the repository root
conda activate panacea
# Install Panacea (use --no-deps because conda provided most packages)
pip install . --no-deps
```

Notes:
- Some steps require significant memory/storage; for full-scale reductions consider running on TACC (see ../tacc/overview.md).


See also: [Quickstart](../getting-started/quickstart.md), [TACC Overview](../tacc/overview.md)
