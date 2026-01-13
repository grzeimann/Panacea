# docs/conf.py
import os
import sys

# Make the package importable for autodoc (../src)
ROOT = os.path.abspath(os.path.join(__file__, '..', '..'))
SRC = os.path.join(ROOT, 'src')
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

project = 'Panacea Documentation'
extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
]
# Autodoc defaults: include public members only (no private/special) and show undocumented public members
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    # Do not include private or special members globally (omit keys so they are excluded)
    # 'private-members': False,
    # 'special-members': False,
    'inherited-members': True,
    'show-inheritance': True,
}
# Generate summary pages for modules/classes automatically
autosummary_generate = True
# Treat .md files as source
source_suffix = {
    '.md': 'markdown',
    '.rst': 'restructuredtext',
}
# Root document (our docs/ already has index.md)
master_doc = 'index'
# Exclude stale autosummary .rst stubs that can cause duplicate indexing
exclude_patterns = [
    'api/panacea.*.rst',
]
# Optional: add anchors to headings so #section links work
myst_heading_anchors = 3
# Optional: enable extra MyST features if desired (linkify is optional)
myst_enable_extensions = []
try:
    import linkify_it  # provided by linkify-it-py
except Exception:
    pass  # leave linkify disabled if backend isn't installed
else:
    myst_enable_extensions.append('linkify')  # autolink URLs in text

# Theme: prefer sphinx_rtd_theme (persistent sidebar), fallback to alabaster
try:
    import sphinx_rtd_theme  # noqa: F401
except Exception:
    html_theme = 'alabaster'
    html_theme_options = {}
else:
    html_theme = 'sphinx_rtd_theme'
    html_theme_options = {
        'navigation_depth': 4,
        'collapse_navigation': False,
        'sticky_navigation': True,
    }