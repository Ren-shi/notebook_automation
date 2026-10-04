"""Write the validation notebooks in notebooks/validation/, one per physics-register page.

    python scripts/make_validation_notebooks.py

Each notebook runs the checks of physim.nuclear.validation for its capabilities, shows the results table and plots
every comparison. The docs build executes them, so the plots always show the current code against the reference
files in tests/reference/nuclear/.
"""

from __future__ import annotations

from pathlib import Path

import nbformat

OUT = Path(__file__).resolve().parents[1] / "notebooks" / "validation"

PAGES = {
    "nuclear-data": ("Validation: atomic masses and material data", ["Atomic masses and material data"],
                     "Nuclear masses against LISE++'s mass table, and Q-values against published values."),
    "kinematics": ("Validation: two-body kinematics", ["Two-body reaction kinematics"],
                   "Lab energies against LISE++'s kinematics calculator, and the exact classical limit."),
    "stopping": ("Validation: stopping, range and straggling",
                 ["Stopping power and range", "Energy and angular straggling"],
                 "Stopping powers, ranges and straggling against LISE++ (ATIMA 1.4), NIST's CSDA ranges and Bohr's "
                 "formula."),
    "rutherford": ("Validation: Rutherford scattering",
                   ["Rutherford cross section", "Distance of closest approach and validity checks",
                    "Coulomb trajectories"],
                   "Cross sections against LISE++ and Geiger and Marsden's 1913 data; orbit geometry against LISE++ "
                   "and the exact Coulomb deflection."),
    "detectors": ("Validation: detector solid angles", ["Detector solid angles and response"],
                  "Solid angles against exact closed forms."),
    "rates-and-events": ("Validation: count rates and Monte Carlo spectra",
                         ["Count rates and beam time", "Monte Carlo spectra"],
                         "Rates against flux × target × cross section × solid angle, and simulated peaks against "
                         "the analytic expectation."),
}


def notebook(title: str, capabilities: list, intro: str) -> nbformat.NotebookNode:
    nb = nbformat.v4.new_notebook()
    cells = [
        nbformat.v4.new_markdown_cell(
            f"# {title}\n\n{intro}\n\nEach check is defined in `physim.nuclear.validation`; the same checks run in "
            "the test suite (`tests/python/test_nuclear_validation.py`), and the physics register may show ✅ only "
            "for checks that pass. 🟡 means the reference has been requested but not produced yet (see "
            "`tests/reference/nuclear/pending/`)."),
        nbformat.v4.new_code_cell(
            "import matplotlib.pyplot as plt\n"
            "from IPython.display import Markdown\n"
            "from physim.nuclear import validation\n\n"
            f"capabilities = {capabilities!r}\n"
            "results = [r for c in capabilities for r in validation.run(capability=c)]\n"
            "Markdown(validation.report(results))"),
        nbformat.v4.new_markdown_cell("## Comparisons\n\nRed dots are physim, open squares the reference."),
        nbformat.v4.new_code_cell(
            "shown = [r for r in results if r.plot]\n"
            "fig, axes = plt.subplots(len(shown), 1, figsize=(6.5, 3.2 * len(shown)), squeeze=False)\n"
            "for ax, r in zip(axes[:, 0], shown):\n"
            "    validation.plot(r, ax=ax)\n"
            "fig.tight_layout()"),
    ]
    nb.cells = cells
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    return nb


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, (title, caps, intro) in PAGES.items():
        nbformat.write(notebook(title, caps, intro), OUT / f"{name}.ipynb")
        print("wrote", OUT / f"{name}.ipynb")


if __name__ == "__main__":
    main()
