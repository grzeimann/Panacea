# docs/conf.py
import os
import sys

# Make the package importable for autodoc (../src)
ROOT = os.path.abspath(os.path.join(__file__, "..", ".."))
SRC = os.path.join(ROOT, "src")
if os.path.isdir(SRC) and SRC not in sys.path:
    sys.path.insert(0, SRC)

# Project metadata and titles
project = "Panacea"
html_title = "Panacea"
author = "Greg Zeimann"
release = "1.0.1"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
]
# Autodoc defaults: include public members only (no private/special) and show undocumented public members
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    # Do not include private or special members globally (omit keys so they are excluded)
    # 'private-members': False,
    # 'special-members': False,
    "inherited-members": True,
    "show-inheritance": True,
}
# Generate summary pages for modules/classes automatically
autosummary_generate = True
# Treat .md files as source
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}
# Root document (our docs/ already has index.md)
master_doc = "index"
# Exclude stale autosummary .rst stubs that can cause duplicate indexing
exclude_patterns = [
    "api/panacea.*.rst",
]
# Optional: add anchors to headings so #section links work
myst_heading_anchors = 3
# Optional: enable extra MyST features if desired (linkify is optional)
myst_enable_extensions = []
# Enable MyST linkify only if the backend is available and record status
LINKIFY_ENABLED = False
try:
    pass  # provided by linkify-it-py
except Exception:
    LINKIFY_ENABLED = False
else:
    myst_enable_extensions.append("linkify")  # autolink URLs in text
    LINKIFY_ENABLED = True

# Expose status in HTML metadata for easy inspection and debugging
html_meta = {
    "panacea-linkify": "enabled" if LINKIFY_ENABLED else "disabled",
}

# Print a concise build-time note (helps when scanning Sphinx logs)
print(
    f"[Panacea docs] MyST linkify backend: {'ENABLED' if LINKIFY_ENABLED else 'DISABLED'}"
)

# Theme selection with optional override via SPHINX_THEME
# Supported values: furo, pydata_sphinx_theme, sphinx_book_theme, sphinx_rtd_theme, alabaster
_requested_theme = os.environ.get("SPHINX_THEME", "").strip()
_supported = {
    "furo": "furo",
    "pydata_sphinx_theme": "pydata_sphinx_theme",
    "sphinx_book_theme": "sphinx_book_theme",
    "sphinx_rtd_theme": "sphinx_rtd_theme",
    "alabaster": "alabaster",
}

html_theme_options = {}


def _try_set_theme(name: str) -> bool:
    try:
        __import__(name)
        globals()["html_theme"] = name
        return True
    except Exception:
        return False


# If user requested a specific theme via env var, try that first
if _requested_theme and _requested_theme in _supported:
    if not _try_set_theme(_supported[_requested_theme]):
        print(
            f"[Panacea docs] Requested theme '{_requested_theme}' not installed; falling back."
        )

# If not set by the override above, choose a reasonable default chain
if "html_theme" not in globals():
    if _try_set_theme("furo"):
        html_theme_options = {}
    elif _try_set_theme("pydata_sphinx_theme"):
        html_theme_options = {
            "navigation_depth": 4,
            "show_toc_level": 2,
        }
    elif _try_set_theme("sphinx_book_theme"):
        html_theme_options = {
            "show_toc_level": 2,
        }
    elif _try_set_theme("sphinx_rtd_theme"):
        html_theme_options = {
            "navigation_depth": 4,
            "collapse_navigation": False,
            "sticky_navigation": True,
        }
    else:
        html_theme = "alabaster"
        html_theme_options = {}

# Code block color styles (helps ensure good contrast in light/dark)
pygments_style = "sphinx"
pygments_dark_style = "native"
