"""The beam-time report and the data exports: what a student attaches to a proposal or takes to their own analysis.

::

    from physim.nuclear import Experiment, report

    rep = report.build(Experiment.example("oxygen_on_lead_array"), seed=1)
    rep.write("my-report")        # report.html, CSV tables, figures (PNG and PDF), setup file, events.root

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

    # -- numbers --------------------------------------------------------------------------------------------------

    def numbers(self) -> dict:
        """Every number in the report (no dates or versions): equal for the same setup and seed."""
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

    def to_html(self) -> str:
        """The report as one self-contained HTML page."""
        exp = self.experiment
        e = html.escape
        figs = {}
        for name, fig in self.figures().items():
            buf = io.BytesIO()
            self._save(fig, buf, "png", dpi=200)
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
        mark = {"pass": "✅", "pending": "🟡", "fail": "❌", "none": "—", "unchecked": "see register"}
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
<p class="sub">Beam-time report · physim {e(m['version'])} (commit {e(m['commit'])}) · {e(m['date'])} ·
Monte Carlo seed {self.seed}, {self.events:,} events</p>
{f'<p>{e(exp.description)}</p>' if exp.description else ''}

<h2>Warnings</h2>
<ul class="warnings">{warn_html}</ul>

<h2>Setup</h2>
{table(("", ""), setup_rows)}
<p>Reaction: {e(exp.reaction)} scattering. Beam: {_fmt(self.planner.rates()['particles_per_second'])} particles
per second.</p>
<figure><img alt="Setup in 3D" src="data:image/png;base64,{figs['geometry']}">
<figcaption>Beam, target and detectors (mm).</figcaption></figure>
<figure><img alt="Coverage" src="data:image/png;base64,{figs['coverage']}">
<figcaption>Where each detector sits in (θ, φ) as seen from the target.</figcaption></figure>

<h2>Detectors: coverage, cross sections, rates and beam time</h2>
{table(("Detector", "Shape", "θ (deg)", "φ (deg)", "Ω (msr)", "dσ/dΩ lab (mb/sr)", "Rate (1/s)",
        "Monte Carlo (1/s)", "Counts in run", "Beam time"), det_rows)}
<p>Rates count particles measured above threshold, from the target and its backing (scattered beam and recoils).
dσ/dΩ is the Rutherford cross section of the scattered beam on the main target nuclide at the detector's mean
angle, at the beam energy in the middle of the target. {e(counts_note)}</p>

<h2>Kinematics</h2>
<figure><img alt="Kinematics" src="data:image/png;base64,{figs['kinematics']}">
<figcaption>Lab energy against lab angle; the shaded bands are the detectors' angular coverage.</figcaption></figure>

<h2>Expected peaks</h2>
{table(("Detector", "From", "Particle", "Mean (MeV)", "FWHM (keV)", "Rate (1/s)", "Above threshold"), peak_rows)}

{self._excitation_html()}
<h2>Energy loss</h2>
{table(("Layer", "Material", "Thickness (mg/cm²)", "Beam in (MeV)", "Beam out (MeV)", "Loss (keV)",
        "Straggling FWHM (keV)"), loss_rows)}

<h2>Simulated spectra</h2>
<figure><img alt="Spectra" src="data:image/png;base64,{figs['spectra']}">
<figcaption>Measured-energy spectra, counts per bin in the planned beam time ({self.events:,} events, seed
{self.seed}).</figcaption></figure>

<h2>Models and their validation</h2>
<p>Status in physim's physics register (✅ checked by an automated test, 🟡 reference requested, — not
applicable; "see register": the tool reference files are only in physim's source tree, so the comparison was not run
here). Assumptions: elastic Rutherford scattering; straight tracks; Gaussian straggling and resolution;
detector faces without inter-strip effects. See the register for each model's range of validity.</p>
{table(("Model", "vs tools", "vs literature"), reg_rows)}

<h2>Data and references</h2>
<ul>{''.join(f'<li>{e(s)}</li>' for s in m['data'])}</ul>
<ul>{''.join(f'<li>{e(s)}</li>' for s in REFERENCES)}</ul>
<p class="meta">Generated by physim {e(m['version'])}, commit {e(m['commit'])}, on {e(m['date'])}. The same
setup file and seed reproduce every number.</p>
</body></html>
"""

    def _excitation_html(self) -> str:
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
        return f"""<h2>Coulomb excitation</h2>
<p>{e(s['excite'].capitalize())} excited to {s['energy_kev']:g} keV (B({e(s['multipolarity'])}↑) =
{s['b_up_e2fm']:.4g} e² fm^{2 * int(s['multipolarity'][1])}), first-order semiclassical theory: adiabaticity
ξ = {x['xi']:.2f}, Sommerfeld parameter η = {x['eta']:.1f}, total cross section {x['total_mb']:.3g} mb. Safe
(Cline's criterion, closest approach ≥ {x['safe_distance_fm']:.1f} fm) at {safe}.</p>
<table><thead><tr><th>Detector</th><th>Excitation events (1/s)</th></tr></thead><tbody>{rows}</tbody></table>
{('<h3>γ rays</h3><table><thead><tr><th>Particle detector</th><th>γ detector</th><th>E_γ (keV)</th>'
  '<th>Doppler shift (keV)</th><th>FWHM (keV)</th></tr></thead><tbody>' + dop + '</tbody></table>') if dop else ''}
"""

    def write(self, directory: Union[str, Path], figures: bool = True, root: Optional[bool] = None) -> list:
        """report.html, the CSV tables, the figures (PNG and PDF), the setup file and, with ``uproot`` installed (or
        ``root=True``), ``events.root`` with the simulated events and spectra (:mod:`physim.nuclear.rootio`)."""
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        paths = [d / "report.html"]
        paths[0].write_text(self.to_html(), encoding="utf-8")
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
    if validation:
        from . import validation as v

        have_refs = v.reference_dir() is not None
        results = v.run() if have_refs else v.run(kind="literature")
        for cap in MODELS:
            row = {"capability": cap}
            for kind in ("tool", "literature"):
                mine = [x.status for x in results if x.capability == cap and x.kind == kind]
                if kind == "tool" and not have_refs and any(c.capability == cap and c.kind == "tool"
                                                            for c in v.CHECKS):
                    row[kind] = "unchecked"  # the tool reference files are not installed with physim
                else:
                    row[kind] = ("fail" if "fail" in mine else "pending" if "pending" in mine
                                 else "pass" if mine else "none")
            status.append(row)
    meta = {"version": __version__, "commit": _commit(), "date": _dt.date.today().isoformat(),
            "data": DATA_SOURCES}
    excitation = p.gamma() if experiment.excitation is not None else None
    return Report(experiment, seed, events, meta, detectors, peaks, kinematics, p.energy_loss(),
                  [(w.level, w.text) for w in p.warnings()], status, p, excitation, journal, width)


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
    args = ap.parse_args(argv)
    from .experiment import example_names

    exp = Experiment.example(args.setup) if args.setup in example_names() else Experiment.load(args.setup)
    import matplotlib

    matplotlib.use("Agg")
    try:
        rep = build(exp, seed=args.seed, events=args.events, journal=args.journal, width=args.width)
    except ValueError as err:
        ap.error(str(err))
    paths = rep.write(args.output, figures=not args.no_figures, root=False if args.no_root else None)
    print(f"wrote {len(paths)} files to {args.output}/ (open report.html)")
    return 0


__all__ = ["DATA_SOURCES", "REFERENCES", "Report", "build", "main"]


if __name__ == "__main__":
    raise SystemExit(main())
