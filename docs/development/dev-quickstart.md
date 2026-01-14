# Developer Quickstart (local CI, hooks, and reproducible env)

This page collects the copy‑paste commands for creating a clean development environment, running the same checks as CI, and setting up pre‑commit hooks.

If you already have a `panacea` conda environment

```bash
conda deactivate
conda env remove -n panacea
```

Get the package

```bash
git clone https://github.com/grzeimann/Panacea.git
cd Panacea
```

Create the environment and activate it

```bash
conda env create -f environment.yml
conda activate panacea
```

Install in developer mode

```bash
pip install .[dev]
```

Initial run tests (sanity checks)

```bash
panacea-lrs2 -h
panacea-lrs2 --smoke-test
```

If you are going to make code changes

```bash
pre-commit clean
pre-commit install --hook-type pre-push
make ci
pre-commit run --all-files

# Stage and commit tracked changes
git add -u
git commit -m "message"

# If pre-commit hooks modify files on commit (EOF/trailing whitespace, ruff --fix), re‑stage and commit those auto‑fixes
git add -u && git commit -m "Apply pre-commit auto-fixes"

# Push to your remote
git push
```

Notes
- The pre‑push hook runs tests under coverage in an isolated environment pinned to compatible versions (numpy<2.0 with astropy<6). If it’s too slow for your workflow, you can skip once with `git push --no-verify` (not recommended routinely).
- Makefile targets mirror CI:
  - `make lint` → ruff check
  - `make test` → pytest -q
  - `make coverage` → coverage run/report/xml
  - `make ci` → lint + coverage (closest to CI)
