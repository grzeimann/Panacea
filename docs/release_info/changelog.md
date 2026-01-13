# Changelog

How we write future entries
- Heading format: [X.Y.Z] - YYYY-MM-DD
- Start with a 2–5 line summary of why the release matters for users.
- Use simple sections only when they help: Added, Changed, Fixed, Docs, Dev.
- One bullet per item; link issues/PRs when relevant (e.g., #123).
- Call out any breaking changes clearly.

Example
[unreleased]
- Short bullets about things merged since the last tag.

## [1.0.0] - 2025-11-11 — Initial public release for JOSS
This tag marks the first stable, citable Panacea release submitted to the Journal of Open Source Software (JOSS). It completes the refactor from a single monolithic script into a maintainable, modular Python package with basic tests, uniform documentation, and an open path for collaboration within the HET community.

### Highlights
- Refactor: split the legacy monolithic pipeline into cohesive modules (ccd, fiber, sky, wavelength, trace, routine, astrometry, utils, io) with clear responsibilities and imports.
- CLI: panacea-lrs2 wraps the quicklook reduction with concise, practical options for nightly use.
- Reproducibility: bundles required configuration resources (line lists, DAR tables, responses, fiber locations, fplane) inside the package for consistent results.

### Added
- Initial test coverage in tests/ to exercise parsers, resource loaders, and core calibration utilities.
- Packaged configuration data under panacea/lrs2_config loaded via importlib.resources.

### Docs
- Consolidated documentation site with user guides (installation, CLI, configuration, data products), API reference, and FAQ.
- More uniform, concise docstrings across modules; improved README with links to contributing and citation.

## [0.0.5] - 2019-01-01 — Start of automated processing (legacy CLI)
This internal milestone marks the first automated nightly reductions using the original command-line driver, full_reduction.py — a single ~2000-line monolithic script that orchestrated the pipeline end to end.
- Automatic processing began on 2019-01-01 using the legacy CLI.
- Established end-to-end quicklook reductions with minimal configuration.

---

See also: [Citation](../citation/citation.md), [Contributing](../community/contributing.md)
