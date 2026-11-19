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
    "coulex": ("Validation: Coulomb excitation and γ rays",
               ["Coulomb excitation, first-order semiclassical", "γ-ray Doppler shift and broadening"],
               "First-order Coulomb excitation against its exact results (GOSIA requested), and the Doppler "
               "formula against the Lorentz boost."),
    "rates-and-events": ("Validation: count rates and Monte Carlo spectra",
                         ["Count rates and beam time", "Monte Carlo spectra"],
                         "Rates against flux × target × cross section × solid angle, and simulated peaks against "
                         "the analytic expectation."),
    "fusion-evaporation": ("Validation: fusion-evaporation",
                           ["Fusion-evaporation residues from a channel table",
                            "Fusion cross section (Bass barrier, Wong)", "Evaporation Monte Carlo (Hauser–Feshbach)",
                            "Neutron detectors and channel selection"],
                           "The compound nucleus against the hand calculation; fusion cross sections against PACE4 "
                           "4.36.0, measured excitation functions and the model's exact results; evaporation against "
                           "PACE4, GEMINI++ 3.3.0, measured ¹⁶⁵Ho(α,xn) cross sections and exact results; neutron "
                           "detectors against closed forms, the binomial expectation and ENDF/B-VIII.0 carbon."),
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


def reproduction_notebook() -> nbformat.NotebookNode:
    """The reproduced fusion-evaporation experiment (backlog 75, D10): the four excitation functions and the
    comparisons, read from the saved result (rerunning it takes about 1.5 minutes: scripts/reproduce_usman2020.py)."""
    nb = nbformat.v4.new_notebook()
    nb.cells = [
        nbformat.v4.new_markdown_cell(
            "# Reproduced experiment: ⁴He + ¹⁶⁵Ho, Usman et al. 2020\n\n"
            "A. R. Usman et al., Nucl. Instrum. Methods B 469, 42 (2020) (EXFOR E2650): ¹⁶⁵Ho(α,xn)Tm by stacked-foil "
            "activation at the RIKEN AVF cyclotron, 26 Ho foils from 49.6 to 13.2 MeV, reproduced end to end through "
            "the planner. For each foil the `fusion_tm165` example is set to the foil's energy, its channel table "
            "filled from physim's own fusion and evaporation (\"Fill from physim\"), the residues simulated, and the "
            "production cross sections inferred from the residues made (`physim.nuclear.reproduction`). Kormazeva et "
            "al. 2022 and Gadkari et al. 1997 are the cross-checks; PACE4 4.36.0 ran at the same energies. The "
            "validation check `reproduction_usman2020` reruns the reproduction and judges it; this page reads the "
            "saved result, `tests/reference/nuclear/reproduction_usman2020.csv` (`scripts/reproduce_usman2020.py`)."),
        nbformat.v4.new_code_cell(
            "from physim.nuclear import reproduction as rp\n\n"
            "fig = rp.figure()"),
        nbformat.v4.new_markdown_cell(
            "physim's band is its statistical error. Judged above 1.1 V_B: the peaks of 2n, 3n and 4n (energy spread "
            "+ 1 MeV), the strongest channel at each foil where the measurement resolves it, and the strongest three "
            "channels within 30% of the band of the three measurements up to each channel's peak. Above the peaks "
            "the measured tails are fed by pre-equilibrium emission, which neither physim nor PACE4 models: "
            "reported, not judged."),
        nbformat.v4.new_code_cell(
            "from IPython.display import Markdown\n\n"
            "c = rp.compare(rp.saved(), rp.measured())\n"
            "lines = ['| Comparison | Tolerance | physim | measured | Result |', '|---|---|---|---|---|']\n"
            "for x in c['peaks']:\n"
            "    res = ('pass' if x['passed'] else 'fail') if x['judged'] else 'reported'\n"
            "    lines.append(f\"| {x['channel']} peak | ±{x['tolerance_mev']:.1f} MeV | {x['physim_mev']:.1f} MeV | \"\n"
            "                 f\"{x['measured_mev']:.1f} MeV | {res} |\")\n"
            "s = [x for x in c['strongest'] if x['judged']]\n"
            "lines.append(f\"| strongest channel | the same | {sum(x['passed'] for x in s)} of {len(s)} foils | | \"\n"
            "             f\"{'pass' if c['ok']['strongest'] else 'fail'} |\")\n"
            "m = [x for x in c['magnitudes'] if x['judged']]\n"
            "w = max(m, key=lambda x: abs(x['deviation']))\n"
            "lines.append(f\"| magnitudes ({len(m)} points) | 30% of the band | worst {w['deviation']:+.0%} \"\n"
            "             f\"({w['channel']}, {w['e_lab_mev']:g} MeV) | | {'pass' if c['ok']['magnitudes'] else 'fail'} |\")\n"
            "for x in c['ratio_3n_4n']:\n"
            "    lines.append(f\"| σ(3n)/σ(4n), {x['ref']} {x['e_lab_mev']:g} MeV | reported | {x['physim']:.3f} | \"\n"
            "                 f\"{x['measured']:.3f} | pre-equilibrium |\")\n"
            "Markdown('\\n'.join(lines) + '\\n\\n' + rp.reported(c) + '.')"),
    ]
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    return nb


def lifetimes_notebook() -> nbformat.NotebookNode:
    """Lifetime measurements (backlog 73, wave 2D): the generic page, plus the ⁹⁸Ru attenuation curves and the ⁸³Y
    line shape with AhKin's and the measured one."""
    nb = notebook(
        "Validation: lifetime measurements (RDDS and DSAM)",
        ["Lifetimes: recoil distance (RDDS)", "Lifetimes: Doppler-shift attenuation (DSAM)"],
        "The plunger's decay curves and Bateman's equations against their closed forms, the deorientation in vacuum "
        "against Radeck et al.'s Table I, the recoil transport against TRIM, the DSAM attenuation against its exact "
        "result for stopping proportional to speed, and DSAM line shapes against APCAD's and AhKin's examples "
        "(shape only, not a lifetime fit). The DDCM against NaPyTau, the closure tests and the three published "
        "experiments wait for simulated spectra and the lifetime analysis (wave 3C). The two line-shape comparisons "
        "fail their tolerances; why is on the register page `docs/physics-register/lifetimes.md`.")
    nb.cells += [
        nbformat.v4.new_markdown_cell(
            "## Deorientation in vacuum, ⁹⁸Ru 2⁺₁\n\nRadeck et al.'s Table I parameters: the attenuation of the "
            "unshifted γ rays, G_k at the stopper, and of the shifted ones, G̃_k averaged over the decays in flight "
            "(their Eq. 7), against the target-to-stopper distance."),
        nbformat.v4.new_code_cell(
            "r = next(x for x in results if x.name == 'deorientation_ru98')\n"
            "fig, ax = plt.subplots(figsize=(6.5, 3.2))\n"
            "for k, colour in ((2, '#d62728'), (4, '#1f77b4')):\n"
            "    c = r.plot['curves'][k]\n"
            "    ax.plot(r.plot['x'], c['unshifted'], color=colour, label=f'G{k} (unshifted)')\n"
            "    ax.plot(r.plot['x'], c['shifted'], '--', color=colour, label=f'G̃{k} (shifted)')\n"
            "ax.set_xlabel('distance (µm)')\n"
            "ax.set_ylabel('attenuation')\n"
            "ax.legend(fontsize=8)\n"
            "fig.tight_layout()"),
        nbformat.v4.new_markdown_cell(
            "## The ⁸³Y line shape\n\nAhKin's example at its τ = 0.1313 ps: the measured line, AhKin's best line "
            "and physim's with AhKin's inputs (its stopping table, cascade, detector and resolution)."),
        nbformat.v4.new_code_cell(
            "r = next(x for x in results if x.name == 'dsam_lineshape_vs_ahkin')\n"
            "fig, ax = plt.subplots(figsize=(6.5, 3.2))\n"
            "ax.step(r.plot['x'], r.plot['reference'], where='mid', color='#444444', label='measured')\n"
            "ax.plot(r.plot['x'], r.plot['ahkin'], color='#1f77b4', label='AhKin')\n"
            "ax.plot(r.plot['x'], r.plot['physim'], color='#d62728', label='physim')\n"
            "ax.set_xlabel('energy (keV)')\n"
            "ax.set_ylabel('counts per channel')\n"
            "ax.legend(fontsize=8)\n"
            "fig.tight_layout()"),
    ]
    return nb


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, (title, caps, intro) in PAGES.items():
        nbformat.write(notebook(title, caps, intro), OUT / f"{name}.ipynb")
        print("wrote", OUT / f"{name}.ipynb")
    nbformat.write(reproduction_notebook(), OUT / "reproduction-usman2020.ipynb")
    print("wrote", OUT / "reproduction-usman2020.ipynb")
    nbformat.write(lifetimes_notebook(), OUT / "lifetimes.ipynb")
    print("wrote", OUT / "lifetimes.ipynb")


if __name__ == "__main__":
    main()
