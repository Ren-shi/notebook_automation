"""Sphinx configuration for the physim documentation.

Build: ``pip install -r docs/requirements.txt`` with the extension installed (``maturin develop --release``),
then ``sphinx-build -W -b html docs docs/_build/html``. The example notebooks are copied from ``notebooks/``
and executed during the build; the Rust API docs (``cargo doc --no-deps``) are copied to ``rust/`` by CI.
"""

import shutil
from pathlib import Path

import physim

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent

project = "physim"
author = "The physim authors"
copyright = "2026, The physim authors"
release = physim.__version__
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "myst_nb",
]
myst_enable_extensions = ["dollarmath", "amsmath", "deflist"]
myst_heading_anchors = 3
nb_execution_mode = "force"
nb_execution_timeout = 600
nb_execution_raise_on_error = True
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "requirements.txt"]

autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": False, "show-inheritance": False}
autodoc_typehints = "description"

html_theme = "furo"
html_title = f"physim {release}"


def _copy_notebooks(app):
    """Copy the example and validation notebooks into the source tree so Sphinx can include them."""
    for src, dest in ((ROOT / "notebooks", HERE / "examples"), (ROOT / "notebooks" / "validation", HERE / "validation")):
        dest.mkdir(exist_ok=True)
        for nb in sorted(src.glob("*.ipynb")):
            shutil.copy2(nb, dest / nb.name)


def setup(app):
    app.connect("builder-inited", _copy_notebooks)
