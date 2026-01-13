# Using Sphinx with Markdown (MyST)

This repo’s documentation is written in plain Markdown with relative links (./ and ../). It renders correctly on GitHub and with MkDocs. If you prefer Sphinx, you can keep the same Markdown files and relative links by enabling the MyST parser.

## Why relative links are fine
- GitHub: relative links are resolved based on the current file’s path.
- MkDocs: relative links are resolved against the docs/ directory and work as expected.
- Sphinx + MyST: relative Markdown links between project files are supported and will resolve as document references, provided the linked files are included in Sphinx’s source tree (e.g., under docs/).

For maximum robustness in Sphinx, you may also use Sphinx roles like `{doc}` and section labels, but this is optional if your links are simple `./file.md` or `../path/file.md` references.

## Minimal Sphinx setup (keeping Markdown)

1) Install Sphinx and MyST
```bash
pip install sphinx myst-parser
# Optional: if you enable MyST's 'linkify' feature, install its backend:
pip install linkify-it-py
# Optional: to get a persistent sidebar/navigation like the docs site, install the RTD theme:
pip install sphinx_rtd_theme
```

2) Create docs/conf.py (example)
```python
# docs/conf.py
project = 'Panacea Documentation'
extensions = [
    'myst_parser',
]
# Treat .md files as source
source_suffix = {
    '.md': 'markdown',
    '.rst': 'restructuredtext',
}
# Root document (our docs/ already has index.md)
master_doc = 'index'
# Optional: add anchors to headings so #section links work
myst_heading_anchors = 3
# Optional: enable extra MyST features if desired
myst_enable_extensions = [
    'linkify',  # autolink URLs in text
]
# If using relative image paths, Sphinx resolves them like Markdown
html_theme = 'alabaster'  # or any installed Sphinx theme
```

3) Ensure your index
- You already have `docs/index.md`. Sphinx will use it as the root (`master_doc = 'index'`).
- Make sure all pages you want included are reachable from `index.md` via links, or list them in a toctree (see below).

4) (Optional) Add an index toctree if you want sidebar navigation
Add this block somewhere in `docs/index.md` to drive Sphinx nav (MkDocs will ignore it). NOTE: This is an example shown as a literal code block so it does not affect your build from this page.

```text
{toctree}
:hidden:
:maxdepth: 2

getting-started/installation.md
getting-started/quickstart.md
tacc/overview.md
tacc/running.md
user-guide/cli.md
user-guide/configuration.md
data-products/overview.md
algorithms/overview.md
faq.md
community/contributing.md
community/code_of_conduct.md
release_info/changelog.md
citation/citation.md
```

5) Build
```bash
sphinx-build -b html docs/ _build/html
```
Open `_build/html/index.html` in your browser.

## Notes and caveats
- Relative links like `./tacc/overview.md` or `../data-products/overview.md` will resolve under MyST so long as the target files are in the Sphinx source (`docs/`).
- For cross-page section anchors (e.g., `overview.md#where-to-find-data-products`), set `myst_heading_anchors` (as above) so Sphinx creates stable IDs for headings.
- For advanced cross-references, prefer MyST/Sphinx roles, e.g., `[{doc}](/path)` → `[TACC Overview]({doc}`tacc/overview`)` or labeled headings with `{#my-label}` and use `{ref}` roles. This increases resilience to file renames.
- If you later mix `.rst` and `.md`, both are supported with the `source_suffix` mapping shown above.
- Our existing MkDocs config (`mkdocs.yml`) continues to work independently; you can keep both builders side-by-side.

### Optional: API docs from source code
- Sphinx: Already configured. We enabled `autodoc`, `autosummary`, and `napoleon` in `docs/conf.py` and added `docs/api/index.md`. Build with `sphinx-build -b html docs/ _build/html` to generate API pages under the sidebar "API Reference".
- MkDocs: If you prefer API under MkDocs, install `mkdocstrings[python]` and add an `API` section to `mkdocs.yml`. Example snippet:

```yaml
plugins:
  - mkdocstrings:
      handlers:
        python:
          options:
            docstring_style: google
nav:
  - API:
      - panacea: api/panacea.md
```

Then create a page like `docs/api/panacea.md` with:

```markdown
# panacea

::: panacea
```

See also: [Docs Guide](../guide.md), [Home](../index.md)
