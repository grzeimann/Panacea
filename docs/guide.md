# Documentation Guide

This folder contains the Panacea documentation split from the long README into topic-focused pages. All links are relative so the docs render correctly on GitHub and remain compatible with a future docs site (e.g., MkDocs).

## Structure
- index.md — landing page and entry points
- getting-started/
  - installation.md — local install and dependencies
  - quickstart.md — first-run CLI examples
- tacc/
  - overview.md — accounts, access, and where to find products on TACC
  - running.md — interactive and batch usage on TACC
- user-guide/
  - cli.md — panacea-lrs2 CLI options and examples
  - configuration.md — config files, responses, fplane, line lists
- data-products/
  - overview.md — definitions and directory layout of outputs
- algorithms/
  - overview.md — high-level algorithms and references
- citation/
  - citation.md — how to cite Panacea (see CITATION.cff at the repository root: https://github.com/grzeimann/Panacea/blob/HEAD/CITATION.cff)
- community/
  - contributing.md — contribution guide
  - code_of_conduct.md — community standards
- release_info/
  - changelog.md — link to ../../CHANGELOG.md
- images/ — images referenced from docs and README

## Conventions
- One H1 (#) per page.
- Keep pages focused; prefer new pages when content gets long.
- Use relative links (./ or ../) within docs/.
- Add a short "See also" section at the bottom to connect related pages.
- Add language specifiers to fenced code blocks (```bash, ```text, etc.).

## Linking and compatibility
- Relative links are intentionally used so the same Markdown pages work on:
  - GitHub (links resolve relative to the current file)
  - MkDocs (links resolve within the docs/ tree)
  - Sphinx, when using MyST Markdown (links resolve as document references as long as the target files are part of the Sphinx source)
- For Sphinx adoption without converting formats, see migration notes: ./migration/sphinx.md

## Editing and contributions
- Keep filenames lowercase-with-hyphens to preserve URL stability.
- When adding a new page, link it from index.md and include a short "See also".
- Prefer incremental PRs that add or move a small number of pages at a time.
- If you change paths or filenames, verify all inbound links are updated.

See also: [Home](./index.md), [Installation](./getting-started/installation.md), [TACC Overview](./tacc/overview.md)

## Building the documentation

You can build the docs in two ways. Use Sphinx for the full API reference (autodoc), or MkDocs for a quick brochure site. Both can be run locally.

### Option A: Build with Sphinx (recommended for API pages)
This uses docs/conf.py and will render the API reference via autodoc. Output goes to _build/html.

1) Create/activate an environment and install doc deps
- Using pip only:
  - python -m pip install --upgrade pip
  - pip install .[docs]
- Using conda (optional) then pip for extras:
  - conda env update -f environment.yml  # creates/updates env
  - python -m pip install --upgrade pip
  - pip install .[docs]

2) Build HTML into _build/html
- Clean old build (optional):
  - rm -rf _build/html
- Build:
  - sphinx-build -b html docs _build/html
  - If you run into stale caches, use: sphinx-build -E -a -b html docs _build/html

3) View locally
- Open _build/html/index.html in your browser.

Notes
- docs/conf.py already adds ../src to sys.path, so autodoc can import panacea from a source checkout. If import errors persist, install the package in editable mode: pip install -e .
- If you prefer a theme like sphinx_rtd_theme, install it (auto‑detected in conf.py when available): pip install sphinx-rtd-theme.

### Option B: Build with MkDocs (quick brochure site)
This will build the navigation defined in mkdocs.yml. Sphinx-only blocks like ```{eval-rst} with autodoc will not execute under MkDocs; they render as code blocks.

1) Install MkDocs
- python -m pip install mkdocs
  - (Optional) Add a theme: pip install mkdocs-material

2) Build or serve
- Build to the default site/ directory:
  - mkdocs build
- Or build to _build/mkdocs:
  - mkdocs build -d _build/mkdocs
- Live-reload preview at http://127.0.0.1:8000/:
  - mkdocs serve

3) View locally
- Open site/index.html (or _build/mkdocs/index.html if you set -d).

Common issues and fixes
- ImportError during Sphinx autodoc: Ensure the environment is active and either rely on docs/conf.py’s sys.path tweak or run pip install -e . to make panacea importable.
- Extensions missing: If myst-parser or furo aren’t found, make sure you ran pip install .[docs].
- Stale API pages: Clean the build or pass -E -a to sphinx-build as shown above.
