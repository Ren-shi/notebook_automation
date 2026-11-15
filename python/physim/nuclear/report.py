"""The beam-time report and the data exports: what a student attaches to a proposal or takes to their own analysis.

::

    from physim.nuclear import Experiment, report

    rep = report.build(Experiment.example("oxygen_on_lead_array"), seed=1)
    rep.write("my-report")        # report.html, CSV tables, figures (PNG and PDF), setup file, events.root
    rep.write("my-report", pdf=True)   # and report.pdf, printed by a browser found on this computer

or from a terminal::

    python -m physim.nuclear.report my_setup.toml -o my-report --seed 1

The HTML report is self-contained (figures embedded) and has print styles: the browser's "Print → Save as PDF"
gives the PDF version from the same source. The same setup file and seed always give the same numbers.
"""

from __future__ import annotations

import base64
import csv
import datetime as _dt
import html
import io
import math
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import numpy as np

from . import data
from .detectors import Array
from .experiment import Experiment
from .paper import JOURNALS, WIDTHS
from .planner import Planner
from .rutherford import Rutherford

#: Data shipped with physim that the numbers depend on (see physim/nuclear/data/SOURCES.md).
DATA_SOURCES = [
    "Atomic masses: AME2020 (W. J. Huang et al., Chin. Phys. C 45, 030002; M. Wang et al., Chin. Phys. C 45, "
    "030003 (2021)).",
    "Isotopic compositions: NIST Atomic Weights and Isotopic Compositions (Coursey et al., 2015; Berglund and "
    "Wieser 2009).",
    "Element densities and mean excitation energies: NIST SRD 126 (Hubbell and Seltzer, 2004).",
    "Proton and alpha stopping powers: NIST PSTAR/ASTAR (ICRU Reports 49 and 90).",
]

#: Papers behind the models, for the reference list of a proposal.
REFERENCES = [
    "E. Rutherford, Phil. Mag. 21, 669 (1911): the scattering cross section.",
    "H. Geiger and E. Marsden, Phil. Mag. 25, 604 (1913): its first test (used to validate physim).",
    "ICRU Report 49 (1993) and Report 90 (2016): stopping powers of protons and alpha particles.",
    "W. Brandt and M. Kitagawa, Phys. Rev. B 25, 5631 (1982): effective charge of heavy ions.",
    "J. F. Ziegler, J. P. Biersack and U. Littmark, The Stopping and Range of Ions in Solids (1985): nuclear "
    "stopping (ZBL).",
    "N. Bohr, Phil. Mag. 30, 581 (1915); J. Lindhard and M. Scharff, K. Dan. Vidensk. Selsk. Mat.-Fys. Medd. 27, "
    "15 (1953): energy straggling.",
    "V. L. Highland, Nucl. Instrum. Methods 129, 497 (1975): multiple scattering.",
]

#: Which register capabilities each part of the report uses.
MODELS = ("Atomic masses and material data", "Two-body reaction kinematics", "Stopping power and range",
          "Energy and angular straggling", "Rutherford cross section",
          "Distance of closest approach and validity checks", "Detector solid angles and response",
          "Count rates and beam time", "Monte Carlo spectra")


#: The register capabilities of a fusion-evaporation run. The report picks the validation suite's checks by
#: capability, never by name, so a check added to one of these rows appears in the report without a change here.
#: A capability added later whose name speaks of fusion, evaporation, residues or neutrons is taken too.
FUSION_MODELS = ("Fusion-evaporation residues from a channel table", "Fusion cross section (Bass barrier, Wong)",
                 "Evaporation Monte Carlo (Hauser–Feshbach)", "Neutron detectors and channel selection")
_FUSION_WORDS = ("fusion", "evaporation", "residue", "neutron")

#: Where physim's evaporation was seen to leave too much spin at the end of the chain: chains that start below
#: this spin give PACE4's channels (docs/theory/evaporation.md, "High spin").
HIGH_SPIN_HBAR = 40
#: Quasi-fission is also seen below the Z₁Z₂ ≈ 1600 guide, as low as ³²S + ²⁰⁸Pb (Z₁Z₂ = 1312;
#: docs/theory/fusion.md): from here the report notes it.
QUASI_FISSION_NOTE_Z1Z2 = 1300

#: The documented limits of the fusion-evaporation models (backlog item 75 "Left out", docs/theory/fusion.md,
#: evaporation.md, neutrons.md), each with the results it affects: "fusion" (σ_fus), "residues" (the channel
#: table and its rates), "selection" (the channel-selection matrix). :func:`fusion_limits` adds whether this run
#: falls where the limit bites.
LIMITS = (
    {"key": "coupled_channels", "affects": ("fusion",),
     "title": "No coupled-channels fusion below the barrier",
     "text": "One parabolic barrier (Bass, Wong) without coupling to the nuclei's excited states or to transfer. "
             "Below the barrier this underestimates fusion, by up to orders of magnitude, so σ_fus there is a "
             "lower bound; near it (up to V_B + ħω) coupled channels can change it noticeably.",
     "source": "docs/theory/fusion.md (Where the model stops being valid); backlog 75, Left out"},
    {"key": "incomplete_fusion", "affects": ("fusion", "residues"),
     "title": "No incomplete fusion and no pre-equilibrium emission",
     "text": "Complete fusion only, and every particle is evaporated from an equilibrated compound nucleus. Above "
             "about 10 MeV/u break-up of the projectile and pre-equilibrium emission take a growing share: "
             "complete fusion and the low-multiplicity channels are then overestimated, the high-energy "
             "particles missing.",
     "source": "docs/theory/fusion.md; docs/theory/evaporation.md (Limits); backlog 75, Left out"},
    {"key": "quasi_fission", "affects": ("fusion",),
     "title": "No quasi-fission",
     "text": "σ_fus is the capture cross section. Above Z₁Z₂ ≈ 1600 the system often re-separates before a "
             "compound nucleus forms, and capture overestimates fusion, by orders of magnitude for the heaviest "
             "systems; it is also seen in some lighter systems with heavy or deformed targets (³²S + ²⁰⁸Pb, "
             "Z₁Z₂ = 1312).",
     "source": "docs/theory/fusion.md; backlog 75, Left out"},
    {"key": "high_spin", "affects": ("residues",),
     "title": "The evaporation at high spin",
     "text": "physim's evaporation has no statistical E2 γ rays, so for heavy beams at high spin it leaves too much "
             "energy in rotation and emits too few neutrons (⁴⁸Ca + ¹²⁴Sn: 4n and 5n strongest instead of 6n "
             "and 5n; ³²S + ⁹²Mo off by up to 10× against PACE4). Chains that start below 40 ħ agree. Its "
             "neutrons are 10–25% softer than PACE4's at every spin.",
     "source": "docs/theory/evaporation.md (What agrees and what does not); backlog 75, Status"},
    {"key": "selection", "affects": ("selection",),
     "title": "Channel-selection simplifications still left out",
     "text": "Protons and α particles lose no energy in the target on their way out; the particle detectors tell "
             "protons from α particles perfectly (a ΔE-only silicon cannot, for punch-through particles); γ-ray "
             "cross-talk between neutron cells is left out (a γ ray that fires a cell is either called a γ ray or "
             "one false neutron).",
     "source": "docs/theory/neutrons.md (Limits); backlog 75, Done for the neutron detectors"},
)


#: The reaction as the report names it.
_REACTION_WORDS = {"elastic": "elastic scattering", "coulex": "Coulomb excitation (and elastic scattering)",
                   "fusion-evaporation": "fusion-evaporation (the particle detectors see the elastic scattering)"}


def _commit() -> str:
    """The git commit of the physim source, if it is a git checkout; otherwise "unknown"."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=Path(__file__).resolve().parent,
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def _fmt(x, digits: int = 4) -> str:
    if x is None:
        return "—"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, float) and not math.isfinite(x):
        return "∞" if x > 0 else "—"
    return f"{x:.{digits}g}"


def _hours(seconds) -> str:
    if seconds is None:
        return "—"
    if not math.isfinite(seconds):
        return "never"
    if seconds >= 360:
        return f"{seconds / 3600:.3g} h"
    return f"{seconds / 60:.3g} min" if seconds >= 60 else f"{seconds:.3g} s"


def _table(head, rows) -> str:
    e = html.escape
    out = ["<table><thead><tr>" + "".join(f"<th>{e(h)}</th>" for h in head) + "</tr></thead><tbody>"]
    for r in rows:
        out.append("<tr>" + "".join(f"<td>{e(str(c))}</td>" for c in r) + "</tr>")
    return "\n".join(out + ["</tbody></table>"])


@dataclass
class Report:
    """Everything in a beam-time report, as data. Build it with :func:`build`."""

    experiment: Experiment
    seed: int
    events: int
    #: physim version, commit, date, data sources.
    meta: dict
    #: One row per detector (see :meth:`detector_table`).
    detectors: list
    #: One row per expected peak.
    peaks: list
    #: Lab angle grid and energies of ejectiles and recoils per target nuclide.
    kinematics: dict
    energy_loss: dict
    #: (level, text) pairs, most important first.
    warnings: list
    #: One row per model used: capability, tool status, literature status.
    register: list
    planner: Planner = field(repr=False, default=None)
    #: Coulomb excitation (:meth:`physim.nuclear.planner.Planner.gamma`), or None for elastic setups.
    excitation: Optional[dict] = None
    #: The journal style of the figures (a key of :data:`physim.nuclear.paper.JOURNALS`) and their width.
    journal: str = "physical_review"
    width: str = "single"
    #: Fusion-evaporation, or None for other reactions: {"summary" (σ_fus, the barrier, ℓ; the planner's
    #: fusion summary), "e_cm_mev", "residues" (the channel table with rates), "selection" (the channel-selection
    #: matrix, with a neutron wall), "limits" (:func:`fusion_limits`: each documented limit and whether this run
    #: falls where it bites), "checks" (:func:`fusion_checks`, with ``validation``)}.
    fusion: Optional[dict] = None

    # -- numbers --------------------------------------------------------------------------------------------------

    def numbers(self) -> dict:
        """Every number in the report (no dates or versions): equal for the same setup and seed."""
        out = self._numbers()
        if self.fusion is not None:
            f = self.fusion
            out["fusion"] = {"summary": f["summary"], "residues": (f["residues"] or {}).get("rows"),
                             "selection": None if f["selection"] is None else f["selection"]["matrix"],
                             "flags": [lim["flag"] for lim in f["limits"]]}
        return out

    def _numbers(self) -> dict:
        return {"detectors": self.detectors, "peaks": self.peaks,
                "kinematics": {k: [None if not math.isfinite(x) else x for x in np.asarray(v, dtype=float).tolist()]
                               for k, v in self.kinematics.items()},
                "energy_loss": self.energy_loss["layers"], "warnings": self.warnings,
                "excitation": None if self.excitation is None else {
                    k: v for k, v in self.excitation.items() if k not in ("theta_cm", "probability")}}

    # -- CSV ------------------------------------------------------------------------------------------------------

    def write_csv(self, directory: Union[str, Path]) -> list:
        """detectors.csv, peaks.csv, kinematics.csv, energy_loss.csv and strips.csv (units in the headers)."""
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        paths = []

        def write(name, rows):
            p = d / name
            with open(p, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            paths.append(p)

        write("detectors.csv", self.detectors)
        write("peaks.csv", self.peaks)
        keys = list(self.kinematics)
        write("kinematics.csv", [{k: _round(self.kinematics[k][i]) for k in keys}
                                 for i in range(len(self.kinematics[keys[0]]))])
        write("energy_loss.csv", [{"layer": r["layer"], "material": r["material"],
                                   "thickness_mg_cm2": r["thickness_mg_cm2"], "energy_in_MeV": r["energy_in_mev"],
                                   "energy_out_MeV": r["energy_out_mev"], "loss_MeV": r["loss_mev"],
                                   "straggling_fwhm_MeV": r["straggling_fwhm_mev"]}
                                  for r in self.energy_loss["layers"]])
        strips = self.planner.rates()["strips"]
        write("strips.csv", [{"detector": name, "segment_i": seg[0], "segment_j": seg[1], "rate_per_s": rate}
                             for name, segs in strips.items() for seg, rate in segs.items()])
        if self.excitation is not None and self.excitation["doppler"]:
            write("gamma.csv", [{"particle_detector": r["particle_detector"], "gamma_detector": r["gamma_detector"],
                                 "emitter": r["emitter"], "E0_keV": r["e0_kev"], "mean_keV": r["mean_kev"],
                                 "shift_keV": r["shift_kev"], "doppler_fwhm_keV": r["doppler_fwhm_kev"],
                                 "fwhm_keV": r["fwhm_kev"], "beta": r["beta_mean"]}
                                for r in self.excitation["doppler"]])
        return paths

    # -- figures --------------------------------------------------------------------------------------------------

    def figures(self) -> dict:
        """The report's matplotlib figures (geometry, coverage, kinematics, spectra), in the report's journal
        style and at its column width (:mod:`physim.nuclear.paper`)."""
        from . import paper

        figs = {name: paper.figure(self.planner, name, self.journal, self.width)
                for name in ("geometry", "coverage", "kinematics")}
        figs["spectra"] = paper.figure(self.planner, "detector_spectra", self.journal, self.width,
                                       events=self.events, seed=self.seed)
        return figs

    def _save(self, fig, target, fmt: str, **kw) -> None:
        import matplotlib

        from . import paper

        with matplotlib.rc_context(paper.rc(self.journal)):  # the font embedding is decided when saving
            fig.savefig(target, format=fmt, **kw)

    def write_figures(self, directory: Union[str, Path]) -> list:
        """Each figure as PNG (600 dpi) and PDF, at the journal's column width."""
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        paths = []
        for name, fig in self.figures().items():
            for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):
                p = d / f"{name}.{ext}"
                self._save(fig, p, ext, **kw)
                paths.append(p)
        return paths

    # -- HTML -----------------------------------------------------------------------------------------------------

    def parts(self, dpi: int = 200) -> dict:
        """The report's sections as HTML fragments, for :meth:`to_html` and for the report of the experiment
        (:func:`physim.nuclear.logbook.experiment_html`), which holds them next to the runs and their analysis:
        {"warnings", "setup", "setup_figures" (its two figures alone), "detectors", "kinematics", "peaks", "excitation", "excitation_plain" (under a
        smaller heading, without the level schemes), "fusion" and "fusion_plain" (fusion-evaporation: σ_fus, the
        residues and the channel selection, each with its warnings, checks and limits beside it; "" for other
        reactions), "energy_loss", "spectra", "models" (with the fusion-evaporation checks and limits for such a
        run), "references", "meta"}.
        The figures are embedded as PNG at ``dpi``."""
        exp = self.experiment
        e = html.escape
        figs = {}
        for name, fig in self.figures().items():
            buf = io.BytesIO()
            self._save(fig, buf, "png", dpi=dpi)
            figs[name] = base64.b64encode(buf.getvalue()).decode()
        b = exp.beam
        t = exp.target

        def table(head, rows):
            out = ["<table><thead><tr>" + "".join(f"<th>{e(h)}</th>" for h in head) + "</tr></thead><tbody>"]
            for r in rows:
                out.append("<tr>" + "".join(f"<td>{e(str(c))}</td>" for c in r) + "</tr>")
            return "\n".join(out + ["</tbody></table>"])

        level_class = {"error": "err", "warning": "warn", "note": "note"}
        warn_html = "".join(f'<li class="{level_class[lv]}">{e(text)}</li>' for lv, text in self.warnings) \
            or "<li>None.</li>"
        setup_rows = [
            ("Beam", f"{b.nuclide}, {b.energy} ({_fmt(b.energy_mev)} MeV), {b.current}"
                     + (f", charge {b.charge_state}+" if b.charge_state else "")
                     + (f", spread {b.energy_spread}" if b.energy_spread else "")
                     + (f", spot {b.spot_size}" if b.spot_size else "")),
            ("Target", f"{t.material}, {t.thickness}" + (f", tilted {t.tilt}" if t.tilt else "")),
        ]
        if t.backing is not None:
            setup_rows.append(("Backing", f"{t.backing.material}, {t.backing.thickness}"))
        setup_rows.append(("Run", f"{exp.run.beam_time}" + (f", {exp.run.counts_wanted} counts wanted per detector"
                                                            if exp.run.counts_wanted else "")))
        det_rows = []
        for setup_det, row in zip(exp.detectors, self.detectors):
            det_rows.append((row["detector"], setup_det.shape, f"{row['theta_min_deg']:.1f}–{row['theta_max_deg']:.1f}",
                             f"{row['phi_min_deg']:.1f}–{row['phi_max_deg']:.1f}", _fmt(row["solid_angle_msr"]),
                             _fmt(row["dsigma_domega_lab_mb_sr"]), _fmt(row["rate_per_s"]),
                             f"{_fmt(row['mc_rate_per_s'])} ± {_fmt(row['mc_rate_error_per_s'], 2)}",
                             _fmt(row["counts_in_run"]), _hours(row["beam_time_s"])))
        peak_rows = [(p["detector"], p["channel"], p["particle"] + (" (2nd)" if p["branch"] == 2 else ""),
                      f"{p['mean_MeV']:.4f}", f"{p['fwhm_keV']:.1f}", _fmt(p["rate_per_s"]),
                      "yes" if p["above_threshold"] else "no") for p in self.peaks]
        loss_rows = [(r["layer"], r["material"], _fmt(r["thickness_mg_cm2"]), f"{r['energy_in_mev']:.4f}",
                      f"{r['energy_out_mev']:.4f}", f"{r['loss_mev'] * 1e3:.1f}", f"{r['straggling_fwhm_mev'] * 1e3:.1f}")
                     for r in self.energy_loss["layers"]]
        mark = {"pass": "✅", "pending": "🟡", "fail": "❌", "none": "—", "unchecked": "see register",
                "known": "❌ known failure"}
        reg_rows = [(r["capability"], mark[r["tool"]], mark[r["literature"]]) for r in self.register]
        m = self.meta
        what = {"all": "counts", "excitations": "Coulomb-excitation events",
                "coincidences": "particle–γ coincidences from Coulomb excitation"}[self.planner.rates()["measured"]]
        counts_note = (f"Counts in the run and beam time are for {what}; the beam time is for "
                       f"{exp.run.counts_wanted} of them (relative statistical error "
                       f"{100 / math.sqrt(exp.run.counts_wanted):.1f}%)." if exp.run.counts_wanted else "")
        r = self.planner.rates()
        if r["measured"] == "coincidences" and r["gamma_efficiency_typical"]:
            counts_note += (f" The γ-ray detectors are taken to catch {r['gamma_efficiency']:.2%} of the γ rays in "
                            "the full-energy peak. This comes from the typical response of such crystals, not from "
                            "a calibration of these detectors: the coincidence counts and the beam time scale with "
                            "the efficiency you measure.")
        setup_figures = f"""<figure><img alt="Setup in 3D" src="data:image/png;base64,{figs['geometry']}">
<figcaption>Beam, target and detectors (mm).</figcaption></figure>
<figure><img alt="Coverage" src="data:image/png;base64,{figs['coverage']}">
<figcaption>Where each detector sits in (θ, φ) as seen from the target.</figcaption></figure>"""
        return {
            "warnings": f'<ul class="warnings">{warn_html}</ul>',
            "setup": f"""{table(("", ""), setup_rows)}
<p>Reaction: {e(_REACTION_WORDS.get(exp.reaction, exp.reaction))}. Beam: {_fmt(self.planner.rates()['particles_per_second'])} particles
per second.</p>
{setup_figures}""",
            "setup_figures": setup_figures,
            "detectors": f"""{table(("Detector", "Shape", "θ (deg)", "φ (deg)", "Ω (msr)", "dσ/dΩ lab (mb/sr)", "Rate (1/s)",
        "Monte Carlo (1/s)", "Counts in run", "Beam time"), det_rows)}
<p>Rates count particles measured above threshold, from the target and its backing (scattered beam and recoils).
dσ/dΩ is the Rutherford cross section of the scattered beam on the main target nuclide at the detector's mean
angle, at the beam energy in the middle of the target. {e(counts_note)}</p>""",
            "kinematics": f"""<figure><img alt="Kinematics" src="data:image/png;base64,{figs['kinematics']}">
<figcaption>Lab energy against lab angle; the shaded bands are the detectors' angular coverage.</figcaption></figure>""",
            "peaks": table(("Detector", "From", "Particle", "Mean (MeV)", "FWHM (keV)", "Rate (1/s)",
                            "Above threshold"), peak_rows),
            "excitation": self._excitation_html(),
            "excitation_plain": self._excitation_html(level=3, schemes=False),
            "energy_loss": table(("Layer", "Material", "Thickness (mg/cm²)", "Beam in (MeV)", "Beam out (MeV)",
                                  "Loss (keV)", "Straggling FWHM (keV)"), loss_rows),
            "spectra": f"""<figure><img alt="Spectra" src="data:image/png;base64,{figs['spectra']}">
<figcaption>Measured-energy spectra, counts per bin in the planned beam time ({self.events:,} events, seed
{self.seed}).</figcaption></figure>""",
            "models": f"""<p>Status in physim's physics register (✅ checked by an automated test, 🟡 reference requested, — not
applicable; "see register": the tool reference files are only in physim's source tree, so the comparison was not run
here). Assumptions: elastic Rutherford scattering; straight tracks; Gaussian straggling and resolution;
detector faces without inter-strip effects. See the register for each model's range of validity.</p>
{table(("Model", "vs tools", "vs literature"), reg_rows)}""" + (
                "\n" + self._fusion_models_html() if self.fusion is not None else ""),
            "fusion": self._fusion_html(),
            "fusion_plain": self._fusion_html(level=3),
            "references": f"""<ul>{''.join(f'<li>{e(x)}</li>' for x in m['data'])}</ul>
<ul>{''.join(f'<li>{e(x)}</li>' for x in REFERENCES)}</ul>""",
            "meta": (f"physim {e(m['version'])} (commit {e(m['commit'])}) · {e(m['date'])} · Monte Carlo seed "
                     f"{self.seed}, {self.events:,} events"),
        }

    def to_html(self) -> str:
        """The report as one self-contained HTML page."""
        exp = self.experiment
        e = html.escape
        p = self.parts()
        m = self.meta
        return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{e(exp.title)}: beam-time report</title>
<style>
body {{ font-family: system-ui, -apple-system, "Segoe UI", sans-serif; max-width: 60rem; margin: 2rem auto;
  padding: 0 1rem; color: #1b1b1b; line-height: 1.45; }}
h1 {{ font-size: 1.6rem; margin-bottom: 0.2rem; }} h2 {{ font-size: 1.2rem; margin-top: 2rem;
  border-bottom: 1px solid #ccc; }}
.sub {{ color: #555; margin-top: 0; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.86rem; margin: 0.6rem 0; }}
th, td {{ border-bottom: 1px solid #ddd; padding: 0.25rem 0.4rem; text-align: left; vertical-align: top; }}
th {{ background: #f3f3f3; }}
ul.warnings li {{ margin: 0.2rem 0; }} .err {{ color: #a40000; font-weight: 600; }} .warn {{ color: #8a4b00; }}
.note {{ color: #444; }}
img {{ max-width: 100%; }} figure {{ margin: 1rem 0; }} figcaption {{ font-size: 0.85rem; color: #555; }}
.meta {{ font-size: 0.8rem; color: #555; }}
@media print {{ body {{ margin: 0; max-width: none; }} h2 {{ break-after: avoid; }}
  figure, table {{ break-inside: avoid; }} }}
</style></head><body>
<h1>{e(exp.title)}</h1>
<p class="sub">Beam-time report · {p['meta']}</p>
{f'<p>{e(exp.description)}</p>' if exp.description else ''}

<h2>Warnings</h2>
{p['warnings']}

<h2>Setup</h2>
{p['setup']}

<h2>Detectors: coverage, cross sections, rates and beam time</h2>
{p['detectors']}

<h2>Kinematics</h2>
{p['kinematics']}

<h2>Expected peaks</h2>
{p['peaks']}

{p['excitation']}
{p['fusion'] and p['fusion'] + '''
'''}<h2>Energy loss</h2>
{p['energy_loss']}

<h2>Simulated spectra</h2>
{p['spectra']}

<h2>Models and their validation</h2>
{p['models']}

<h2>Data and references</h2>
{p['references']}
<p class="meta">Generated by physim {e(m['version'])}, commit {e(m['commit'])}, on {e(m['date'])}. The same
setup file and seed reproduce every number.</p>
</body></html>
"""

    def _excitation_html(self, level: int = 2, schemes: bool = True) -> str:
        """The Coulomb-excitation section, under a heading of ``level``, with the level schemes (``schemes``)."""
        x = self.excitation
        if x is None:
            return ""
        e = html.escape
        s = x["state"]
        rows = "".join(f"<tr><td>{e(k)}</td><td>{_fmt(v)}</td></tr>" for k, v in x["rates"].items())
        dop = "".join(f"<tr><td>{e(r['particle_detector'])}</td><td>{e(r['gamma_detector'])}</td>"
                      f"<td>{r['mean_kev']:.2f}</td><td>{r['shift_kev']:+.2f}</td><td>{r['fwhm_kev']:.2f}</td></tr>"
                      for r in x["doppler"])
        safe = ("all angles" if x["max_safe_angle"] >= 180 else f"CM angles up to {x['max_safe_angle']:.0f}°")
        drawn = ""
        if schemes:
            try:
                from .logbook import level_schemes_html

                drawn = level_schemes_html(self.planner)
            except Exception:  # noqa: BLE001 - the report stands without the figure
                drawn = ""
        sub = level + 1
        return f"""<h{level}>Coulomb excitation</h{level}>
<p>{e(s['excite'].capitalize())} excited to {s['energy_kev']:g} keV (B({e(s['multipolarity'])}↑) =
{s['b_up_e2fm']:.4g} e² fm^{2 * int(s['multipolarity'][1])}), first-order semiclassical theory: adiabaticity
ξ = {x['xi']:.2f}, Sommerfeld parameter η = {x['eta']:.1f}, total cross section {x['total_mb']:.3g} mb. Safe
(Cline's criterion, closest approach ≥ {x['safe_distance_fm']:.1f} fm) at {safe}.</p>
{drawn}
<table><thead><tr><th>Detector</th><th>Excitation events (1/s)</th></tr></thead><tbody>{rows}</tbody></table>
{(f'<h{sub}>γ rays</h{sub}><table><thead><tr><th>Particle detector</th><th>γ detector</th><th>E_γ (keV)</th>'
  '<th>Doppler shift (keV)</th><th>FWHM (keV)</th></tr></thead><tbody>' + dop + '</tbody></table>') if dop else ''}
"""

    # -- fusion-evaporation ---------------------------------------------------------------------------------------

    _MARK = {"pass": "✅", "fail": "❌", "known failure": "❌", "pending": "🟡", "not run here": "—"}

    def _attached(self, where: str) -> str:
        """What goes beside a fusion-evaporation result: this run's warnings from the limits that cover it, the
        checks of its models (a known failure said plainly, with the suite's numbers for the system nearest this
        run), and the limits that apply to it but are not reached here."""
        e = html.escape
        f = self.fusion
        out = []
        for lim in f["limits"]:
            if where in lim["affects"] and lim["flag"]:
                word = "Warning" if lim["level"] == "warning" else "Note"
                cls = "warn" if lim["level"] == "warning" else "note"
                out.append(f'<p class="{cls} limit-flag">⚠ <b>{word} ({e(lim["title"].lower())}):</b> '
                           f"{e(lim['flag'])} <span class=\"note\">{e(lim['text'])}</span></p>")
        checks = [c for c in f["checks"] if c["affects"] == where]
        if checks:
            items = []
            for c in checks:
                who = c["checker"] + (f" {c['version']}" if c["version"] not in ("—", "not stated") else "")
                near = c["nearest"]
                near_txt = ""
                if near:
                    near_txt = (f" Nearest compared system: {near['system']}"
                                + (" (yours)" if near["same"] else " (your beam)" if near["same_beam"] else "")
                                + (f": {near['quote']}" if near["quote"] else "") + ".")
                if c["status"] == "known failure":
                    items.append(f'<li class="warn"><b>Known failure:</b> against {e(who)} this model is outside '
                                 f"its tolerance ({e(c['tolerance'])}) over the systems compared. The register "
                                 f"documents the disagreement; read this number with it.{e(near_txt)}</li>")
                elif c["status"] == "fail":
                    items.append(f'<li class="err"><b>Fails</b> against {e(who)} ({e(c["tolerance"])}).'
                                 f"{e(near_txt)}</li>")
                elif c["status"] == "pass":
                    items.append(f"<li>✅ Passes against {e(who)} (tolerance: {e(c['tolerance'])}).{e(near_txt)}</li>")
                else:
                    items.append(f"<li>{self._MARK[c['status']]} {e(who)}: {e(c['status'])}.</li>")
            out.append('<p class="checked-by">Checked against:</p><ul class="checks">' + "".join(items) + "</ul>")
        quiet = [lim for lim in f["limits"] if where in lim["affects"] and not lim["flag"]]
        if quiet:
            out.append('<p class="note limits">Limits (not reached by this run): '
                       + " ".join(f"<b>{e(lim['title'])}.</b> {e(lim['text'])}" for lim in quiet) + "</p>")
        return "\n".join(out)

    def _fusion_html(self, level: int = 2) -> str:
        """The fusion-evaporation section under a heading of ``level``: σ_fus, the residues and the channel
        selection, each with its warnings, checks and limits beside it."""
        f = self.fusion
        if f is None:
            return ""
        e = html.escape
        sub = level + 1
        out = [f"<h{level}>Fusion and evaporation</h{level}>"]
        s = f["summary"]
        if s:
            ratio = f"{f['e_cm_mev'] / s['barrier_cm_mev']:.3f}" if f["e_cm_mev"] else "—"
            rows = [("Compound nucleus", f"{s['compound']} ({s['beam']} + {s['target']})"),
                    ("Fusion barrier V_B (Bass 1980)",
                     f"{s['barrier_cm_mev']:.2f} MeV CM, {s['barrier_lab_mev']:.2f} MeV lab, ħω = "
                     f"{s['hbar_omega_mev']:.2f} MeV, R_B = {s['barrier_radius_fm']:.2f} fm"),
                    ("Beam in the target",
                     f"{s['e_in_mev']:.2f} → {s['e_out_mev']:.2f} MeV lab; E_cm / V_B = {ratio} at the front"),
                    ("σ_fus at the beam energy (Wong)", f"{_fmt(s['sigma_mb'])} mb"),
                    ("σ_fus averaged over the target", f"{_fmt(s['sigma_thick_mb'])} mb"),
                    ("Angular momentum", f"⟨ℓ⟩ = {s['mean_l']:.1f} ħ, up to {s['max_l']:.0f} ħ (critical "
                                         f"ℓ {s['critical_l']})"),
                    ("Excitation energy E*", f"{s['excitation_mev']:.2f} MeV ({s['excitation_range_mev'][0]:.2f}–"
                                             f"{s['excitation_range_mev'][1]:.2f} through the target)")]
            out += [f"<h{sub}>Fusion cross section</h{sub}>", _table(("Quantity", "Value"), rows),
                    self._attached("fusion")]
        r = f["residues"]
        if r:
            rows = [(x["residue"], x["channel"], _fmt(x["cross_section_mb"]), x["source"], _fmt(x["rate_per_s"]),
                     _fmt(x["residues_in_run"]), _fmt(x.get("entry_spin")), _hours(x.get("beam_time_s")))
                    for x in r["rows"]]
            own = [x for x in r["rows"] if not str(x.get("source", "")).startswith("physim")]
            note = ""
            if own and len(own) == len(r["rows"]):
                srcs = sorted({str(x.get("source") or "the setup") for x in own})
                note = (f"<p>The cross sections are the setup's own (source: {e(', '.join(srcs))}): physim's "
                        "evaporation did not make them, so its checks below bear on them only if the table is filled "
                        "from physim's evaporation.</p>")
            out += [f"<h{sub}>Residues</h{sub}>",
                    _table(("Residue", "Channel", "σ (mb)", "Source", "Residues (1/s)", "In the run",
                            "Entry spin (ħ)", "Beam time"), rows), note, self._attached("residues")]
        sel = f["selection"]
        if sel:
            rows = []
            for i, name in enumerate(sel["true"]):
                row = sel["matrix"][i]
                j = sel["assigned"].index(name) if name in sel["assigned"] else None
                right = row[j] if j is not None else None
                rows.append((name, "—" if right is None else f"{100 * right:.1f}%",
                             f"{100 * (1 - (right or 0)):.1f}%"))
            out += [f"<h{sub}>Channel selection</h{sub}>",
                    f"<p>The share of each channel's events that the neutron wall and the particle detectors sort "
                    f"into that channel (the planner's channel-selection Monte Carlo, {sel['n_events']:,} events, "
                    f"seed {sel['seed']}); \"—\" where no event can show it.</p>",
                    _table(("True channel", "Assigned to itself", "Elsewhere or unassigned"), rows),
                    self._attached("selection")]
        elif self.experiment.neutron_wall is None:
            out.append("<p>No neutron wall in the setup: no channel-selection matrix.</p>")
        return "\n".join(x for x in out if x)

    def _fusion_models_html(self) -> str:
        """Every fusion-evaporation check that applies, with its checker, version, tolerance, outcome and the
        suite's numbers for the system nearest this run; then every documented limit and whether this run
        reaches it."""
        f = self.fusion
        if f is None:
            return ""
        out = []
        if f["checks"]:
            rows = []
            for c in f["checks"]:
                near = c["nearest"]
                if near:
                    cmp_ = f"{near['system']}" + (" (your system)" if near["same"] else
                                                  " (your beam)" if near["same_beam"] else
                                                  (f" → {near['compound']}" if near["compound"] else ""))
                    cmp_ += f": {near['quote']}" if near["quote"] else f" (all systems: {c['detail']})"
                else:
                    cmp_ = c["detail"] or "—"
                rows.append((c["capability"], c["kind"], c["checker"], c["version"], c["tolerance"],
                             f"{self._MARK[c['status']]} {c['status']}", _fmt(c["worst"], 3), cmp_))
            known = any(c["status"] == "known failure" for c in f["checks"])
            out.append("<p><b>Fusion-evaporation checks.</b> Every check of the validation suite on the models this "
                       "run uses, picked by model. \"Nearest\" is the compared system closest to yours in beam and "
                       "compound nucleus, with the suite's numbers for it (physim against the reference, as the "
                       "suite writes them); where a check compares no systems, its overall result."
                       + (" A known failure is a disagreement the register documents, with the tolerance kept: "
                          "the numbers it bears on are marked where they appear." if known else "") + "</p>")
            out.append(_table(("Model", "Kind", "Checker", "Version", "Tolerance", "Outcome",
                               "Worst (tolerance's units)", "Nearest to your system, or overall"), rows))
        rows = [(lim["title"], ", ".join(lim["affects"]),
                 (("⚠ " + lim["flag"]) if lim["flag"] else "not reached"), lim["source"]) for lim in f["limits"]]
        out.append("<p><b>Limits of the fusion-evaporation models</b> (backlog item 75 and the theory notes), and "
                   "whether this run reaches them. Each is also given beside the numbers it affects.</p>")
        out.append(_table(("Limit", "Affects", "This run", "Documented in"), rows))
        return "\n".join(out)

    def write(self, directory: Union[str, Path], figures: bool = True, root: Optional[bool] = None,
              pdf: bool = False) -> list:
        """report.html (the report: this plan, and every number with its explanation, in one page,
        :func:`physim.nuclear.logbook.experiment_html`), ``report.pdf`` with ``pdf`` where a browser of this
        computer can print it, the CSV tables, the figures (PNG and PDF), the setup file and, with ``uproot``
        installed (or ``root=True``), ``events.root`` with the simulated events and spectra
        (:mod:`physim.nuclear.rootio`)."""
        from . import logbook

        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        paths = [d / "report.html"]
        page = logbook.experiment_html(self.planner, plan=self)
        paths[0].write_text(page, encoding="utf-8")
        if pdf:
            data = logbook.html_to_pdf(page)
            if data:
                (d / "report.pdf").write_bytes(data)
                paths.append(d / "report.pdf")
        paths += self.write_csv(d)
        if figures:
            paths += self.write_figures(d / "figures")
        setup = d / "setup.toml"
        setup.write_text(self.experiment.to_toml(), encoding="utf-8")
        paths.append(setup)
        if root is None:
            try:
                import uproot  # noqa: F401

                root = True
            except ImportError:
                root = False
        if root:
            from .rootio import write_root

            ev = self.planner.spectra(events=self.events, seed=self.seed)["events"]
            paths.append(write_root(ev, d / "events.root", self.experiment))
        return paths


def _round(x):
    return float(f"{x:.10g}") if isinstance(x, (float, np.floating)) and math.isfinite(x) else x


# -- fusion-evaporation: the checks that apply and the limits ---------------------------------------------------


def fusion_capabilities(capabilities) -> list:
    """The register capabilities that bear on a fusion-evaporation run, in the register's order:
    :data:`FUSION_MODELS` and any later one whose name speaks of fusion, evaporation, residues or neutrons."""
    return [c for c in capabilities if c in FUSION_MODELS or any(w in c.lower() for w in _FUSION_WORDS)]


def _affects(capability: str) -> str:
    """Which fusion-evaporation result a capability's checks bear on: "residues", "selection" or "fusion"."""
    c = capability.lower()
    if "residue" in c or "evaporation monte carlo" in c or ("evaporation" in c and "fusion" not in c):
        return "residues"
    if "neutron" in c or "selection" in c:
        return "selection"
    return "fusion"


def known_failures() -> frozenset:
    """The checks the validation suite's tests hold as known failures: disagreements the register documents, with
    their tolerances unchanged (``KNOWN_FAILURES`` in ``tests/python/test_nuclear_validation.py``). Read from the
    source tree, the one place they are listed, so the report follows the suite; empty where the tests are not
    installed (a wheel), and then a failing check is reported as a failure."""
    import ast

    from . import validation as v

    places = []
    ref = v.reference_dir()
    if ref is not None and len(Path(ref).parents) >= 2:
        places.append(Path(ref).parents[1] / "python" / "test_nuclear_validation.py")
    places += [p / "tests" / "python" / "test_nuclear_validation.py" for p in Path(__file__).resolve().parents]
    for path in places:
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (OSError, SyntaxError, ValueError):
            continue
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "KNOWN_FAILURES" for t in node.targets):
                try:
                    return frozenset(ast.literal_eval(node.value))
                except ValueError:
                    return frozenset()
    return frozenset()


_SUPERSCRIPTS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")
_NUCLIDE = r"(?<![\w.])(\d{1,3}[A-Z][a-z]?)(?![a-z\d])"
_PAIR = re.compile(_NUCLIDE + r"\s*\+\s*" + _NUCLIDE)
_LEAD = re.compile(r"(?:^|[;:]\s*)" + _NUCLIDE + r"(?=[:\s(])")


def _plain(text: str) -> str:
    """Text with superscript mass numbers written out (¹²C → 12C) and α as ⁴He, so systems can be compared."""
    return re.sub(r"(?<![\w])α(?=\s*\+\s*\d{1,3}[A-Z])", "4He", text.translate(_SUPERSCRIPTS))


def _nuclide(name: str):
    try:
        return data.nuclide(name)
    except Exception:  # noqa: BLE001 - not a nuclide: a word that looks like one
        return None


def _systems(text: str) -> list:
    """The (beam, target) systems a check's text names, in order, as "12C + 150Nd" pairs; where it names none,
    the beams that lead its parts ("12C: …") with the target None."""
    out = []
    for b, t in _PAIR.findall(text):
        if _nuclide(b) is not None and _nuclide(t) is not None and (b, t) not in out:
            out.append((b, t))
    if not out:
        for b in _LEAD.findall(text):
            if _nuclide(b) is not None and (b, None) not in out:
                out.append((b, None))
    return out


def _distance(system, beam, compound) -> float:
    """How far a compared system is from the run's: in beam mass and compound-nucleus mass, on a log scale."""
    b = _nuclide(system[0])
    d = abs(math.log(b.A / beam.A))
    if system[1] is not None:
        d += abs(math.log((b.A + _nuclide(system[1]).A) / compound.A))
    return d


def _quote(detail: str, system, systems) -> str:
    """The parts of a check's detail about one compared system, with the label of the list each part is in: the
    suite's own numbers for that system, nothing worked out here. Empty if the detail does not split by system."""
    def about(text, sys_) -> bool:
        """Whether a piece of text is about a system: it names the pair, or names the beam alone (not as a
        target) where no other compared system has that beam."""
        bb, tt = sys_
        if tt is not None and re.search(rf"(?<![\w.]){bb}\s*\+\s*{tt}(?![a-z\d])", text):
            return True
        return sum(1 for s in systems if s[0] == bb) == 1 and bool(
            re.search(rf"(?<!\+\s)(?<!\+)(?<![\w.]){bb}(?![a-z\d])(?!\s*\+)", text))

    def names_any(text) -> bool:
        return any(about(text, s) for s in systems)

    out, label = [], ""  # (label, text) pairs; pieces under one label are joined
    for piece in (p.strip() for p in re.split(r";\s*", detail)):
        head, sep, body = piece.partition(": ")
        if sep and not names_any(head):
            label = head
            items = body.split(", ")
            if sum(names_any(i) for i in items) >= 2:  # a list over the systems: keep this system's items
                mine = [i for i in items if about(i, system)]
                if mine:
                    out.append((head, ", ".join(mine)))
                continue
            if about(body, system):
                out.append((head, body))
            continue
        if about(piece, system):
            out.append((label, piece))
    joined = []
    for lab, text in out:
        if joined and joined[-1][0] == lab:
            joined[-1][1].append(text)
        else:
            joined.append((lab, [text]))
    return "; ".join((f"{lab}: " if lab else "") + ", ".join(texts) for lab, texts in joined)


def _checker(result) -> tuple:
    """(checker, version) of a check: the tool and its version for tool checks ("PACE4", "4.36.0"); for
    literature, data and exact results the reference itself and "—" (a publication, or exact)."""
    ref = result.reference
    if result.kind == "tool":
        m = re.match(r"\s*([A-Za-z][\w+\-]*(?:\+\+)?)\s+v?(\d+(?:\.\d+)+\w*)", ref)
        if m:
            return m.group(1), m.group(2)
        return ref.split(",")[0].split(":")[0].split("(")[0].strip(), "not stated"
    return ref.split(":")[0].strip() if ":" in ref[:120] else ref, "—"


def fusion_checks(results, experiment, have_refs: bool = True, known=None) -> list:
    """One row per check of the fusion-evaporation capabilities among ``results``
    (:func:`physim.nuclear.validation.run`): the model, the checker and its version, the tolerance, the outcome
    ("pass", "fail", "known failure", "pending", or "not run here" for a tool check whose reference files are not
    installed), the worst deviation, the suite's detail, and the compared system nearest the run's in beam and
    compound nucleus with the suite's numbers for it (``nearest``, None where the check compares no systems)."""
    from . import validation as v

    known = known_failures() if known is None else known
    beam = data.nuclide((experiment.beam.Z, experiment.beam.A))
    target = data.nuclide(_main_nuclide(experiment))
    compound = data.nuclide((beam.Z + target.Z, beam.A + target.A))
    caps = fusion_capabilities(v.CAPABILITIES)
    rows, seen = [], set()
    for r in results:
        if r.capability not in caps:
            continue
        seen.add(r.name)
        rows.append(_check_row(r, beam, target, compound, known))
    if not have_refs:  # the tool checks that could not run here are listed, not hidden
        for c in v.CHECKS:
            if c.capability in caps and c.name not in seen:
                r = v.Result(c.capability, c.kind, c.name, c.reference, c.tolerance, "unchecked")
                rows.append(_check_row(r, beam, target, compound, known))
    order = {c: i for i, c in enumerate(caps)}
    return sorted(rows, key=lambda x: order[x["capability"]])


def _check_row(r, beam, target, compound, known) -> dict:
    status = {"pass": "pass", "pending": "pending", "unchecked": "not run here"}.get(
        r.status, "known failure" if r.name in known else "fail")
    checker, version = _checker(r)
    detail = _plain(r.detail)
    systems = _systems(_plain(r.reference) + "; " + detail)
    nearest = None
    if systems:
        best = min(systems, key=lambda s: _distance(s, beam, compound))
        same = best[0] == beam.name and best[1] == target.name
        same_beam = best[0] == beam.name
        cn = (None if best[1] is None else
              data.nuclide((_nuclide(best[0]).Z + _nuclide(best[1]).Z, _nuclide(best[0]).A + _nuclide(best[1]).A)).name)
        nearest = {"system": best[0] + (f" + {best[1]}" if best[1] else ""), "same": same,
                   "same_beam": same_beam, "compound": cn,
                   "quote": _quote(detail, best, systems)}
    return {"capability": r.capability, "affects": _affects(r.capability), "kind": r.kind, "name": r.name,
            "checker": checker, "version": version, "reference": r.reference, "tolerance": r.tolerance,
            "status": status, "worst": None if not math.isfinite(r.worst) else _round(r.worst),
            "detail": r.detail, "nearest": nearest}


def _main_nuclide(experiment) -> tuple:
    """The (Z, A) of the target's main nuclide, the one the beam fuses with."""
    from .rates import stack

    return max(stack(experiment)[0].nuclides.items(), key=lambda kv: kv[1])[0]


def fusion_limits(experiment, fusion: Optional[dict], residues: Optional[dict] = None,
                  selection: Optional[dict] = None) -> list:
    """The documented limits (:data:`LIMITS`), each with ``flag``: None where this run is outside the region the
    limit covers, else the plain warning for this run (and ``level`` "warning" or "note"). ``fusion`` is the
    planner's fusion summary (:meth:`physim.nuclear.planner.Planner.fusion` without its arrays)."""
    from .exit_channels import fusion_of

    b = experiment.beam
    out = []
    f = fusion or {}
    fus = None
    try:
        fus = fusion_of(experiment)
    except Exception:  # noqa: BLE001 - no fusion model for this setup: nothing is flagged on its numbers
        pass
    rows = (residues or {}).get("rows", [])
    from_physim = any(str(r.get("source", "")).startswith("physim") for r in rows)
    for lim in LIMITS:
        flag, level = None, "warning"
        key = lim["key"]
        if key == "coupled_channels" and fus is not None and f:
            vb, hw = f["barrier_cm_mev"], f["hbar_omega_mev"]
            e_in, e_out = float(fus.e_cm(f["e_in_mev"])), float(fus.e_cm(max(f["e_out_mev"], 1e-9)))
            if e_in < vb:
                flag = (f"The beam is below the fusion barrier (E_cm = {e_in:.1f} MeV against V_B = {vb:.1f} MeV): "
                        f"σ_fus = {_fmt(f['sigma_mb'], 3)} mb is a lower bound, possibly by orders of magnitude.")
            elif e_out < vb:
                flag = (f"The beam drops below the fusion barrier inside the target (E_cm {e_in:.1f} → "
                        f"{e_out:.1f} MeV, V_B = {vb:.1f} MeV): the back of the target's share of σ_fus is a lower "
                        "bound.")
            elif e_out < vb + hw:
                flag = (f"The beam is near the fusion barrier (E_cm − V_B = {e_out - vb:.1f} MeV at the back of the "
                        f"target, within ħω = {hw:.1f} MeV): coupled channels can change σ_fus noticeably.")
        elif key == "incomplete_fusion":
            per_u = b.energy_mev / b.A
            from .fusion import INCOMPLETE_FUSION_MEV_PER_U as lim_u

            if per_u > lim_u:
                flag = (f"At {per_u:.1f} MeV/u, above about {lim_u:g} MeV/u: incomplete fusion and pre-equilibrium "
                        "emission take a sizeable share, so σ_fus and the low-multiplicity channels are "
                        "overestimated.")
        elif key == "quasi_fission" and f:
            from .fusion import QUASI_FISSION_Z1Z2 as lim_z

            t = data.nuclide(_main_nuclide(experiment))
            z1z2 = b.Z * t.Z
            if z1z2 > lim_z:
                flag = (f"Z₁Z₂ = {z1z2} is above about {lim_z}: quasi-fission competes, and σ_fus (capture) "
                        "overestimates fusion, often by orders of magnitude.")
            elif z1z2 >= QUASI_FISSION_NOTE_Z1Z2:
                flag, level = (f"Z₁Z₂ = {z1z2}: below the {lim_z} guide, but quasi-fission is seen from about "
                               f"{QUASI_FISSION_NOTE_Z1Z2} with heavy or deformed targets; σ_fus may be high."), "note"
        elif key == "high_spin" and f and f.get("max_l") is not None:
            lmax = f["max_l"]
            if lmax >= HIGH_SPIN_HBAR:
                where = ("the channel table comes from physim's evaporation" if from_physim else
                         "the channel table is the setup's own, so this matters if you fill it from physim's "
                         "evaporation")
                flag = (f"The fusion spins reach {lmax:.0f} ħ, above the {HIGH_SPIN_HBAR} ħ below which physim's "
                        f"evaporation agrees with PACE4: expect too few neutrons per event and the wrong strongest "
                        f"channel ({where}).")
                level = "warning" if from_physim else "note"
        elif key == "selection" and selection is not None:
            labels = [str(x) for x in selection.get("true", [])]
            has_p = any(re.search(r"p", s) for s in labels)
            has_a = any(re.search(r"a|α", s) for s in labels)
            if has_p and has_a:
                flag, level = ("Your channels have both protons and α particles: their separation is taken to be "
                               "perfect, so the p and α rows of the matrix are an upper bound on what the particle "
                               "detectors can tell apart."), "note"
        out.append(dict(lim, flag=flag, level=level if flag else None))
    return out


def _fusion_data(planner, experiment) -> Optional[dict]:
    """The fusion-evaporation numbers of the report: the fusion summary, the residue table and, with a neutron
    wall, the channel selection; None for other reactions."""
    if experiment.fusion is None:
        return None
    summary = planner._plain_fusion()
    try:
        residues = planner.residues()
    except Exception:  # noqa: BLE001 - the report stands without the table; the warnings say why
        residues = None
    selection = None
    if experiment.neutron_wall is not None:
        try:
            selection = planner.channel_selection()
        except Exception:  # noqa: BLE001 - as above
            selection = None
    e_cm = None
    if summary:
        try:
            from .exit_channels import fusion_of

            e_cm = float(fusion_of(experiment).e_cm(summary["e_in_mev"]))
        except Exception:  # noqa: BLE001
            e_cm = None
    return {"summary": summary, "e_cm_mev": e_cm, "residues": residues, "selection": selection,
            "limits": fusion_limits(experiment, summary, residues, selection), "checks": []}


def build(experiment: Experiment, seed: int = 1, events: int = 200_000, validation: bool = True,
          journal: str = "physical_review", width: str = "single") -> Report:
    """Compute everything in the report. ``validation`` runs the register checks to show each model's status (about
    1.5 s); without it the status is left out. ``journal`` and ``width`` set the figures' style and column width
    (:data:`physim.nuclear.paper.JOURNALS`)."""
    from .. import __version__
    from . import paper

    paper.width_mm(journal, width)  # a wrong choice fails here, not after the simulation

    p = Planner(experiment)
    r = p._rates()
    ev = p.spectra(events=events, seed=seed)["events"]
    beam_e = experiment.beam.energy_mev
    ion = data.nuclide((experiment.beam.Z, experiment.beam.A)).name
    main = max(r.layers[0].nuclides.items(), key=lambda kv: kv[1])[0]
    from .rates import beam_energy_at

    e_mid = float(beam_energy_at(experiment, 0, [r.layers[0].thickness / 2], r.layers)[0])
    ruth = Rutherford(ion, data.nuclide(main).name, e_mid)
    detectors = []
    for g, row in zip(Array.from_experiment(experiment), p.rates()["rows"]):
        lo, hi = row["theta_range"]
        plo, phi = g.phi_range()
        mc, mc_err = ev.rate(g.name)
        sigma = ruth.cross_section_lab(g.mean_theta())[0]
        detectors.append({
            "detector": g.name, "theta_min_deg": _round(lo), "theta_max_deg": _round(hi),
            "theta_mean_deg": _round(g.mean_theta()), "phi_min_deg": _round(plo), "phi_max_deg": _round(phi),
            "solid_angle_msr": _round(row["solid_angle_msr"]),
            "dsigma_domega_lab_mb_sr": _round(float(sigma)) if np.isfinite(sigma) else None,
            "rate_per_s": _round(row["rate_per_s"]), "rate_all_per_s": _round(row["rate_all_per_s"]),
            "mc_rate_per_s": _round(mc), "mc_rate_error_per_s": _round(mc_err),
            "counts_in_run": _round(row["counts_in_run"]),
            "beam_time_s": _round(row["beam_time_s"]) if row["beam_time_s"] is not None else None,
            "relative_error": _round(row["relative_error"]),
        })
    peaks = []
    for g in r.array:
        for pk in r.peaks(g.name):
            peaks.append({"detector": g.name, "channel": pk.channel, "particle": pk.particle, "branch": pk.branch,
                          "mean_MeV": _round(pk.mean), "fwhm_keV": _round(pk.fwhm * 1e3),
                          "rate_per_s": _round(pk.rate), "above_threshold": bool(pk.visible)})
    kin = p.kinematics()
    grid = np.arange(0.0, 180.5, 1.0)
    kinematics = {"theta_lab_deg": grid}
    for c in kin["curves"]:
        th, en = np.asarray(c["theta"]), np.asarray(c["energy"])
        # A recoil left at rest (θ* = 180°) has no defined direction: drop points past the maximum angle.
        ok = th <= c["max_angle"] + 1e-9
        th, en = th[ok], en[ok]
        # The curve is monotonic in lab angle on either side of its maximum angle; where both sides reach an angle
        # (double-valued kinematics) the table keeps the higher-energy solution.
        k = int(np.argmax(th))
        vals = np.full(grid.shape, np.nan)
        for part in (slice(0, k + 1), slice(k, None)):
            t, y = th[part], en[part]
            if len(t) < 2:
                continue
            order = np.argsort(t)
            v = np.interp(grid, t[order], y[order], left=np.nan, right=np.nan)
            vals = np.fmax(vals, v)
        kinematics[f"E_{c['label'].replace(' ', '_')}_MeV"] = vals
    status = []
    fusion = _fusion_data(p, experiment)
    if validation:
        from . import validation as v

        have_refs = v.reference_dir() is not None
        models = list(MODELS) + (fusion_capabilities(v.CAPABILITIES) if fusion is not None else [])
        # Only the checks of the models this report shows: the fusion-evaporation ones take most of the time.
        results = [r for cap in models for r in v.run(capability=cap, kind=None if have_refs else "literature")]
        if fusion is not None:
            fusion["checks"] = fusion_checks(results, experiment, have_refs)
        for cap in models:
            row = {"capability": cap}
            for kind in ("tool", "literature"):
                mine = [x.status for x in results if x.capability == cap and x.kind == kind]
                if kind == "tool" and not have_refs and any(c.capability == cap and c.kind == "tool"
                                                            for c in v.CHECKS):
                    row[kind] = "unchecked"  # the tool reference files are not installed with physim
                else:
                    row[kind] = ("fail" if "fail" in mine else "pending" if "pending" in mine
                                 else "pass" if mine else "none")
                if row[kind] == "fail" and cap not in MODELS and fusion is not None:
                    failed = {c["name"] for c in fusion["checks"] if c["capability"] == cap and c["kind"] == kind
                              and c["status"] in ("fail", "known failure")}
                    if failed and all(c["status"] == "known failure" for c in fusion["checks"]
                                      if c["name"] in failed):
                        row[kind] = "known"  # every failure here is one the register documents
            status.append(row)
    meta = {"version": __version__, "commit": _commit(), "date": _dt.date.today().isoformat(),
            "data": DATA_SOURCES}
    excitation = p.gamma() if experiment.excitation is not None else None
    return Report(experiment, seed, events, meta, detectors, peaks, kinematics, p.energy_loss(),
                  [(w.level, w.text) for w in p.warnings()], status, p, excitation, journal, width, fusion)


def main(argv: Optional[list] = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(prog="python -m physim.nuclear.report",
                                 description="Write a beam-time report (HTML, CSV, figures) for a setup file.")
    ap.add_argument("setup", help="setup file (.toml), or the name of an example setup")
    ap.add_argument("-o", "--output", default="report", help="output folder (default: report)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--events", type=int, default=200_000)
    ap.add_argument("--no-figures", action="store_true", help="skip the PNG/PDF figure files")
    ap.add_argument("--journal", default="physical_review", choices=list(JOURNALS),
                    help="journal style of the figures (default: physical_review)")
    ap.add_argument("--width", default="single", choices=list(WIDTHS),
                    help="column width of the figures (default: single; 'middle' only where the journal has one)")
    ap.add_argument("--no-root", action="store_true", help="skip events.root (written when uproot is installed)")
    ap.add_argument("--pdf", action="store_true",
                    help="also write report.pdf, printed by Edge, Chrome or Chromium if one is installed")
    args = ap.parse_args(argv)
    from .experiment import example_names

    exp = Experiment.example(args.setup) if args.setup in example_names() else Experiment.load(args.setup)
    import matplotlib

    matplotlib.use("Agg")
    try:
        rep = build(exp, seed=args.seed, events=args.events, journal=args.journal, width=args.width)
    except ValueError as err:
        ap.error(str(err))
    paths = rep.write(args.output, figures=not args.no_figures, root=False if args.no_root else None,
                      pdf=args.pdf)
    print(f"wrote {len(paths)} files to {args.output}/ (open report.html)")
    return 0


__all__ = ["DATA_SOURCES", "REFERENCES", "Report", "build", "main"]


if __name__ == "__main__":
    raise SystemExit(main())
