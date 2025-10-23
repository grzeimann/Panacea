# Changelog

A lightweight, human-first summary of what changed between releases. Keep it short, friendly, and useful. The newest release goes on top.

How we write future entries
- Heading format: [X.Y.Z] - YYYY-MM-DD
- Start with a 2–5 line summary of why the release matters for users.
- Use simple sections only when they help: Added, Changed, Fixed, Docs, Dev.
- One bullet per item; link issues/PRs when relevant (e.g., #123).
- Call out any breaking changes clearly.

Example
[unreleased]
- Short bullets about things merged since the last tag.

## [0.9.2] - 2025-11-15
- One-line highlights here.

## [0.9.1] - 2025-10-23 — Initial JOSS submission
This is the first public-ready cut submitted to JOSS. The focus is on making Panacea easy to review and reproduce: clearer docs, a small but meaningful test slice, and simple command-line ergonomics for day-to-day use.

### Highlights
- New CLI option `--baseraw` in panacea-lrs2 (run_panacea.py) to point at a raw-data root without touching code.
- Tests: parser/unit coverage for `read_arc_lines`; integrity checks for packaged DAR and skyline tables.
- Docs: README touch-ups (link CONTRIBUTING and CODE OF CONDUCT, clarify optional pyhetdex features) and a clearer, concise description of the automatic IFU-detection workflow in the paper.

### Changed
- Polished inline comments and internal docs to align with the JOSS checklist and improve readability.

### Fixed
- No functional bugs addressed in this cut; this release focuses on documentation, tests, and packaging hygiene.

