# JOSS Submission Readiness

This page tracks Panacea’s readiness for submission to the Journal of Open Source Software (JOSS). It summarizes the JOSS checklist, what’s already done, and what remains.

Status last updated: 2026-01-13

References
- JOSS Author Guide: https://joss.theoj.org/about
- Submission Guidelines: https://joss.readthedocs.io/en/latest/submitting.html
- Review Checklist: https://joss.readthedocs.io/en/latest/review_checklist.html

## 1) Core Submission Artifacts

- Software repository is public and archived at a versioned release with DOI (e.g., Zenodo)
  - Status: TODO
  - Notes: Set up Zenodo GitHub integration; create a release (e.g., v1.0.0) and mint DOI. Add badge to README and CITATION.cff.

- Open-source license present and OSI-approved
  - Status: DONE
  - Where: LICENSE (MIT)

- Author list and affiliations
  - Status: DONE (single author)
  - Where: paper/paper.md
  - Notes: Confirm institutional naming per JOSS and ORCID correctness.

- Versioned release associated with the paper
  - Status: IN PROGRESS
  - Notes: docs/release_info/changelog.md contains [1.0.0] planned tag; publish a matching Git tag and GitHub Release before submission.

- Machine-readable citation metadata (CITATION.cff)
  - Status: DONE
  - Where: CITATION.cff

## 2) Paper Requirements (short, 250–1000 words)

- Title and summary describing the software’s purpose and domain
  - Status: DONE
  - Where: paper/paper.md (Summary)

- Statement of need: who the target audience is and what problem it solves
  - Status: DONE
  - Where: paper/paper.md (Statement of Need)

- State of the field / prior work with references
  - Status: DONE
  - Where: paper/paper.md (State of the Field)

- Software description: key functionality, algorithms, and implementation
  - Status: DONE
  - Where: paper/paper.md (Methods and Implementation)

- Validation/Performance: evidence it works for stated claims
  - Status: DONE
  - Where: paper/paper.md (Validation)

- Acknowledgements and references with DOIs when available
  - Status: IN PROGRESS
  - Where: paper/paper.bib
  - Notes: Ensure all references have DOIs where possible; run a DOI audit and add missing DOIs.

## 3) Documentation & Usability

- Installation instructions
  - Status: DONE
  - Where: docs/getting-started/installation.md, README.md

- Quickstart / usage examples
  - Status: DONE
  - Where: docs/getting-started/quickstart.md, docs/user-guide/cli.md

- API documentation (as appropriate) and explanation of core concepts
  - Status: DONE
  - Where: docs/api (utils, routine, astrometry, io), docs/algorithms/overview.md, docs/user-guide/examples.md
  - Notes: Expanded API coverage using Sphinx autodoc blocks embedded in MyST Markdown; added concrete examples in the User Guide.

- Community guidelines (how to contribute, report issues, and code of conduct)
  - Status: DONE
  - Where: docs/community/contributing.md, docs/community/code_of_conduct.md

## 4) Tests and Continuous Integration

- Automated tests that exercise major functionality
  - Status: IN PROGRESS
  - Where: tests/
  - Notes: Synthetic sample dataset available (tests/fixtures/sample_data.py). Tests cover core calibration utilities, parsers, packaged resources, and environment/data-presence checks. Suite runs under CI with coverage; continue expanding end-to-end pipeline tests.

- Continuous integration (CI) running tests
  - Status: DONE
  - Notes: GitHub Actions workflow (.github/workflows/python-tests.yml) runs Ruff lint (pinned 0.14.x) and the test suite under coverage.py on Python 3.11. Coverage XML is generated and uploaded as an artifact. A separate workflow checks external links with Lychee, and MkDocs deploys the brochure site to GitHub Pages on docs changes.

  - Canonical local commands: For exact, copy‑paste steps to create a clean environment and run lint/tests/coverage, see README.md sections “Reproducible dev environment and test commands” and “Run CI checks locally (pre‑commit + Makefile)”. For exact git stage/commit/push commands, see “Stage, commit, push (exact commands)” in docs/development/dev-quickstart.md. This readiness page summarizes and links only.

### What are “linting” and “coverage reporting”?
- Linting: Automated checks that look for common mistakes, inconsistencies, and style issues in code (e.g., unused imports, undefined names, dead code, formatting). Linters catch problems early and keep the codebase readable and consistent. In this project we use ruff for fast linting. Typical local command: `ruff check .`.
- Coverage reporting: A measure of how much of the code is executed while the tests run. It helps reveal untested paths and guide where to add tests. Typical local commands: `pytest --cov=panacea --cov-report=term-missing` (summary in terminal) or add `--cov-report=xml` to generate an XML file for CI services.

## 5) Functionality & Installability (reviewer checks)

- The software installs and runs as documented
  - Status: IN PROGRESS
  - Notes: Verify pip install from a clean environment using environment.yml/pyproject.toml; ensure panacea-lrs2 CLI works without TACC-specific paths when using bundled configs.

- Statement of need is addressed and functionality is sufficiently documented
  - Status: DONE

- Appropriate tests exist and pass
  - Status: IN PROGRESS

## 6) Archiving & Metadata

- Archive the exact release submitted (e.g., Zenodo) and include the DOI in the paper and README
  - Status: TODO
  - Notes: After making the submission tag, create Zenodo snapshot, update README badge and paper/paper.md references if required by JOSS.

- Ensure authorship and affiliations are correct and complete
  - Status: DONE (confirm at submission)

## Action Items Before Submission

- [x] Add/expand tests to cover CLI and core utilities; ensure tests pass locally and in CI.
- [x] Add GitHub Actions workflow for tests (see .github/workflows/python-tests.yml).
- [x] Add linting to CI (ruff) — see .github/workflows/python-tests.yml.
- [x] Add coverage reporting in CI and publish summary/artifact.
  - Implemented: GitHub Actions runs tests under `coverage run -m pytest -q`, then generates reports with `coverage xml -o coverage.xml` and `coverage report -m`. The `coverage.xml` is uploaded as an artifact per Python version.
- [x] Add link checking in CI (Lychee) and fix broken links (see .github/workflows/link-check.yml, lychee.toml).
- [x] Ensure docs build without warnings in MkDocs; deploy site from main (see .github/workflows/mkdocs-deploy.yml).
- [ ] Publish v1.0.0 (or appropriate) release; ensure CHANGELOG and pyproject.toml version match.
- [ ] Enable Zenodo archiving; mint DOI on release; add DOI badge to README and CITATION.cff.
- [ ] Audit paper.bib for DOIs; update entries.
- [ ] Sanity-check installation from a clean environment and basic CLI run with bundled configs.

## Exact walkthrough: Cut release, archive on Zenodo, submit to JOSS

This section gives a precise, copy‑pasteable sequence to: 1) prepare a tagged release, 2) archive it on Zenodo and capture a DOI, and 3) submit to JOSS with the correct metadata.

Prerequisites
- You have push rights on GitHub and a Zenodo account (https://zenodo.org) with GitHub integration enabled.
- Local tests pass and docs build cleanly. See Dev Quickstart for commands.

Step 1 — Prepare the release candidate branch and version
```bash
# Create a release candidate branch
git checkout -b rc/v1.0.0

# Update project version and changelog (edit files)
# - pyproject.toml: set [project].version = "1.0.0"
# - docs/release_info/changelog.md: ensure a 1.0.0 section with today’s date and summary
# - paper/paper.md: confirm date and content
# - CITATION.cff: ensure version and metadata are current (DOI added in Step 3)

# Commit your edits
git add pyproject.toml docs/release_info/changelog.md paper/paper.md CITATION.cff
git commit -m "Prepare v1.0.0: bump version, update changelog and paper metadata"
```

Step 2 — Verify from a clean environment (matches CI)
```bash
# Optional: fresh env
conda env remove -n panacea || true
conda env create -f environment.yml
conda activate panacea
pip install .[dev]

# Lint + tests + coverage
ruff check .
coverage run -m pytest -q && coverage report -m && coverage xml -o coverage.xml

# CLI smoke test
panacea-lrs2 -h
panacea-lrs2 --smoke-test

# Docs build (HTML)
pip install -U .[docs]
sphinx-build -b html docs docs/_build/html
```

Step 3 — Tag and publish the GitHub Release
```bash
# Tag the exact commit that passed verification
git tag -a v1.0.0 -m "Panacea v1.0.0"

git push origin rc/v1.0.0
git push origin v1.0.0
```
- Create a GitHub Release from tag v1.0.0 (on the web UI) and paste the 1.0.0 notes from docs/release_info/changelog.md. Publish the release.
- Wait for GitHub Actions to finish and show green.

Step 4 — Archive on Zenodo and mint DOI
1) Enable GitHub repository in Zenodo (one‑time)
- Log in to https://zenodo.org with your ORCID or GitHub.
- Go to Account → GitHub and flip the switch to enable archiving for grzeimann/Panacea.

2) Trigger the Zenodo snapshot
- After enabling, re‑publish the GitHub Release (edit → publish) if Zenodo didn’t capture it automatically.
- Zenodo will create:
  - A concept DOI (stable across versions)
  - A version‑specific DOI for v1.0.0

3) Complete Zenodo metadata for this release
- Title: Panacea v1.0.0
- Upload type: Software
- Creators: include full names and ORCID for the author(s)
- License: BSD‑3‑Clause (matches LICENSE)
- Related identifiers: add the GitHub release URL (isSupplementTo) and the software repository URL (isSupplementTo/compiles)
- Keywords: astronomy; spectroscopy; IFU; HET; LRS2; pipeline
- Save and publish. Copy both DOIs (concept and version).

Step 5 — Add DOI everywhere in the repository
```bash
# README.md: add a DOI badge using the concept DOI
# Example badge (replace with your concept DOI):
# [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

# CITATION.cff: set doi: 10.5281/zenodo.XXXXXXX (concept DOI)
# paper/paper.md: add the archive DOI in the header if desired by JOSS, or cite in text

git add README.md CITATION.cff paper/paper.md
git commit -m "Add Zenodo DOI badge and citation metadata"
```

Step 6 — Run a DOI audit on references and rebuild docs (optional but recommended)
```bash
# Dry-run
python scripts/doi_audit.py paper/paper.bib
# Apply updates (creates paper/paper.bib.bak)
python scripts/doi_audit.py paper/paper.bib --write
# Rebuild docs
pip install -U .[docs]
sphinx-build -b html docs docs/_build/html
```
Notes
- The audit uses Crossref; only high‑confidence matches are written automatically.

Step 7 — Submit to JOSS
A) Quick checks with Whedon (optional but useful)
```bash
# In a GitHub issue on your repo, post these comments for Whedon (the JOSS bot):
# @whedon generate pdf
# @whedon check references
```
- Whedon will build a draft PDF of paper/paper.md and report any reference issues.

B) Submit via JOSS web form
- Open https://joss.theoj.org/papers/new
- Fill in:
  - Repository URL: https://github.com/grzeimann/Panacea
  - Version: 1.0.0 (must match the Git tag)
  - Archive DOI: the Zenodo version‑specific DOI for v1.0.0
  - Title: matches paper/paper.md
  - Authors with ORCIDs and affiliations
  - Keywords and summary
- Submit and note the tracking issue URL.

Step 8 — After acceptance (later)
- Update README.md badges with the JOSS DOI.
- Add the JOSS article DOI to CITATION.cff and docs/citation page.

Troubleshooting
- Zenodo didn’t capture my release: ensure the repo is enabled in Zenodo, then edit and re‑publish the GitHub Release to retrigger.
- DOI badge shows 404: concept DOI may take a minute to propagate. Verify the badge URL matches Zenodo’s provided snippet.
- JOSS reference checks fail: fix DOIs in paper/paper.bib and re‑run Whedon.

## Nice-to-haves (not strictly required but helpful)

- [x] Add small sample dataset or fixtures for tests to avoid large downloads (see tests/fixtures/sample_data.py and pytest fixture in tests/conftest.py).
- [x] Increase API documentation breadth and add more examples in the User Guide.
  - Implemented: Added API pages for utils, routine, astrometry; expanded examples in a new User Guide page (user-guide/examples.md). Updated mkdocs.yml navigation accordingly.
- [ ] Add a Contributor Guide section mapping modules to responsibilities for reviewers.


## Addendum: Role of AI in This Project

Status last updated: 2026-01-13

This project was originally developed and authored 100% by a human (Greg Zeimann). In preparation for the JOSS submission during January 2026, limited AI assistance was used as follows:

- Junie (an autonomous programmer by JetBrains, powered by the GPT-5-2025-08-07 model) assisted with small, well-scoped tasks:
  - CI coverage integration (pytest-cov) and coverage artifact upload in GitHub Actions.
  - Addition of targeted unit tests for core calibration utilities.
  - This documentation addendum clarifying the role of AI.
  All AI-generated changes were reviewed and accepted by the human author before inclusion.

- PyCharm docstring autocomplete was occasionally used to suggest or complete short docstring phrases. These were treated as editor assistance (typing/completion) and were reviewed/edited by the human author.

Important notes:
- No scientific claims, algorithmic designs, or core implementation logic were introduced solely by AI without human validation.
- The core algorithms, design decisions, and initial implementation of Panacea remain human-authored.
- Paper text and results are human-written unless otherwise noted; any AI involvement is limited to the items listed above and is subject to human review.
