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

Verify docs build (Sphinx)

If you recently hit a Sphinx error like "Linkify enabled but not installed", we fixed this by adding linkify-it-py to the docs extras and enabling MyST linkify only when available. Verify locally with:

```bash
# Ensure docs extras (includes linkify-it-py and optional themes) are installed
pip install -U .[docs]

# Build the HTML docs
sphinx-build -b html docs docs/_build/html

# Open the result (macOS example)
open docs/_build/html/index.html  # use xdg-open on Linux
```

Switching Sphinx themes (quick ways)

Option A: One-off build with an environment variable

```bash
# Use Furo
SPHINX_THEME=furo sphinx-build -b html docs docs/_build/html
# Use PyData
SPHINX_THEME=pydata_sphinx_theme sphinx-build -b html docs docs/_build/html
# Use Book theme
SPHINX_THEME=sphinx_book_theme sphinx-build -b html docs docs/_build/html
# Use Read the Docs theme
SPHINX_THEME=sphinx_rtd_theme sphinx-build -b html docs docs/_build/html
```

Option B: Make it the default (edit docs/conf.py)
- Set `html_theme = "pydata_sphinx_theme"` (or another) near the bottom where theme selection occurs, or export `SPHINX_THEME` in your shell/profile.
- The config auto-detects installed themes and falls back gracefully.

Linkify backend status (optional check)
- During the build, you should see a log line like:
  - [Panacea docs] MyST linkify backend: ENABLED or DISABLED
- The generated HTML also includes a meta tag recording the status. View page source and look for:
  - <meta name="panacea-linkify" content="enabled"> or "disabled"
- If DISABLED, this is okay: docs still build and render; only auto-linking of bare URLs is skipped. To enable, install docs extras: `pip install .[docs]`.

Notes
- The pre‑push hook runs tests under coverage in an isolated environment pinned to compatible versions (numpy<2.0 with astropy<6). If it’s too slow for your workflow, you can skip once with `git push --no-verify` (not recommended routinely).
- Makefile targets mirror CI:
  - `make lint` → ruff check
  - `make test` → pytest -q
  - `make coverage` → coverage run/report/xml
  - `make ci` → lint + coverage (closest to CI)


## Stage, commit, push (exact commands)

These commands are safe to copy‑paste. They cover the common cases and a few handy fixes.

Staging and committing tracked changes only

```bash
# Stage modifications and deletions to already-tracked files
git add -u

# Optional: review what will be committed
git status -s
git diff --staged

# Commit
git commit -m "Describe your change"
```

Including new, untracked files (additions)

```bash
# Stage everything, including new files and deletions
git add -A

# Optional: review
git status -s
git diff --staged

# Commit
git commit -m "Describe your change"
```

Pushing your branch

```bash
# Push current branch to origin
git push

# If this is a new local branch tracking a remote for the first time
git push -u origin $(git branch --show-current)
```

If pre-commit modifies files during commit

```bash
# When hooks auto-fix (EOF newline, trailing whitespace, ruff --fix), re-stage and commit the fixes
git add -u && git commit -m "Apply pre-commit auto-fixes"
```

Helpful one-offs

```bash
# See unstaged vs staged changes
git status -s
git diff          # unstaged
git diff --staged # staged

# Unstage everything (keep working tree changes)
git restore --staged .

# Amend the last commit message or include newly-staged changes
git commit --amend --no-edit    # keep the same message

# Bypass hooks once (not recommended routinely)
git push --no-verify

# Clean ignored build artifacts (safe preview first)
git clean -ndX   # preview removal of ignored files
git clean -fdX   # actually remove ignored files
```

Tip: If your push is blocked by the pre-push test hook and you’ve already run the full local CI (`make ci`), you can skip once with `git push --no-verify`. Prefer fixing failures so that hooks remain helpful.

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

Verify docs build (Sphinx)

If you recently hit a Sphinx error like "Linkify enabled but not installed", we fixed this by adding linkify-it-py to the docs extras and enabling MyST linkify only when available. Verify locally with:

```bash
# Ensure docs extras (includes linkify-it-py and optional themes) are installed
pip install -U .[docs]

# Build the HTML docs
sphinx-build -b html docs docs/_build/html

# Open the result (macOS example)
open docs/_build/html/index.html  # use xdg-open on Linux
```

Switching Sphinx themes (quick ways)

Option A: One-off build with an environment variable

```bash
# Use Furo
SPHINX_THEME=furo sphinx-build -b html docs docs/_build/html
# Use PyData
SPHINX_THEME=pydata_sphinx_theme sphinx-build -b html docs docs/_build/html
# Use Book theme
SPHINX_THEME=sphinx_book_theme sphinx-build -b html docs docs/_build/html
# Use Read the Docs theme
SPHINX_THEME=sphinx_rtd_theme sphinx-build -b html docs docs/_build/html
```

Option B: Make it the default (edit docs/conf.py)
- Set `html_theme = "pydata_sphinx_theme"` (or another) near the bottom where theme selection occurs, or export `SPHINX_THEME` in your shell/profile.
- The config auto-detects installed themes and falls back gracefully.

Linkify backend status (optional check)
- During the build, you should see a log line like:
  - [Panacea docs] MyST linkify backend: ENABLED or DISABLED
- The generated HTML also includes a meta tag recording the status. View page source and look for:
  - <meta name="panacea-linkify" content="enabled"> or "disabled"
- If DISABLED, this is okay: docs still build and render; only auto-linking of bare URLs is skipped. To enable, install docs extras: `pip install .[docs]`.

Notes
- The pre‑push hook runs tests under coverage in an isolated environment pinned to compatible versions (numpy<2.0 with astropy<6). If it’s too slow for your workflow, you can skip once with `git push --no-verify` (not recommended routinely).
- Makefile targets mirror CI:
  - `make lint` → ruff check
  - `make test` → pytest -q
  - `make coverage` → coverage run/report/xml
  - `make ci` → lint + coverage (closest to CI)


## Stage, commit, push (exact commands)

These commands are safe to copy‑paste. They cover the common cases and a few handy fixes.

Staging and committing tracked changes only

```bash
# Stage modifications and deletions to already-tracked files
git add -u

# Optional: review what will be committed
git status -s
git diff --staged

# Commit
git commit -m "Describe your change"
```

Including new, untracked files (additions)

```bash
# Stage everything, including new files and deletions
git add -A

# Optional: review
git status -s
git diff --staged

# Commit
git commit -m "Describe your change"
```

Pushing your branch

```bash
# Push current branch to origin
git push

# If this is a new local branch tracking a remote for the first time
git push -u origin $(git branch --show-current)
```

If pre-commit modifies files during commit

```bash
# When hooks auto-fix (EOF newline, trailing whitespace, ruff --fix), re-stage and commit the fixes
git add -u && git commit -m "Apply pre-commit auto-fixes"
```

Helpful one-offs

```bash
# See unstaged vs staged changes
git status -s
git diff          # unstaged
git diff --staged # staged

# Unstage everything (keep working tree changes)
git restore --staged .

# Amend the last commit message or include newly-staged changes
git commit --amend --no-edit    # keep the same message

# Bypass hooks once (not recommended routinely)
git push --no-verify

# Clean ignored build artifacts (safe preview first)
git clean -ndX   # preview removal of ignored files
git clean -fdX   # actually remove ignored files
```

Tip: If your push is blocked by the pre-push test hook and you’ve already run the full local CI (`make ci`), you can skip once with `git push --no-verify`. Prefer fixing failures so that hooks remain helpful.


## Troubleshooting: "nothing added to commit but untracked files present"

This message means you ran `git add -u` (stages only already-tracked files), but your working tree only has untracked files (e.g., build artifacts like `build/` or `src/panacea.egg-info/`). Use one of the following depending on your intent:

1) You want to include brand‑new files in the commit
```bash
git add -A
git status -s
git diff --staged
git commit -m "Include new files"
```

2) You do NOT want to commit generated artifacts; just clean them
```bash
# Preview what would be removed (ignored files only)
git clean -ndX
# Remove ignored files (safe: respects .gitignore)
git clean -fdX
```

3) Files were previously tracked and keep showing up despite .gitignore
```bash
# Stop tracking, but keep them on disk
git rm -r --cached build/ src/panacea.egg-info/
# Commit the change so Git stops tracking those paths
git commit -m "Stop tracking build artifacts"
```

4) Quick status to see what's actually changed
```bash
git status -s
```

After cleaning or staging appropriately, try again with your normal flow:
```bash
git add -u && git commit -m "Prepare v1.0.0: bump version, update changelog, and citation"
```
