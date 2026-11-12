"""The report of an experiment (backlog item 69): one page with the setup, the plan, every run and its spectra, the
gates, the Doppler correction, the analysis, every number explained, and the method.

There is one report. What used to be three pages (the beam-time report of the setup, the record of the experiment
and the run record) are its sections, so the plan, the spectra the runs gave and the numbers' explanations are
read together::

    from physim.nuclear import logbook

    page = logbook.experiment_html(planner, result)          # the report: one self-contained HTML page
    pdf = logbook.experiment_pdf(planner, result)            # the same as a PDF, printed by a browser found here
    data = logbook.experiment_zip(planner, result)           # report.html, report.pdf, the setup, the tables,
                                                             # and the experiment folder with every run
    for x in logbook.run_explanations(planner, planner.run):  # each with its formula, numbers and meaning
        assert x.recompute() == x.value
"""

from __future__ import annotations

import base64
import html
import io
import math
import zipfile
from pathlib import Path
from typing import Optional

import numpy as np

from .record import DOCS, Explanation, _explanation_html, _g, record_css


def _time(seconds: float) -> str:
    if seconds is None or not math.isfinite(seconds):
        return "—"
    if seconds >= 3600:
        return f"{seconds / 3600:.3g} h"
    return f"{seconds / 60:.3g} min" if seconds >= 60 else f"{seconds:.3g} s"


# -- the numbers the runs add ----------------------------------------------------------------------------------------

def run_explanations(planner, run=None, detector: Optional[str] = None) -> list:
    """The explanations of a run's numbers, each with its formula, the numbers substituted, what it means and a
    rule that recomputes it: a detector's measured rate, its counts over the whole run, its busy fraction, the live
    fraction, the coincidences in the peak, each gate's share of the particles, and predicted against measured."""
    run = run or planner.run
    if run is None or run.kind == "source":
        return []
    s = run.summary
    exp = run.experiment
    real, duration = run.real_s, run.duration_s
    dets = list(s.get("detectors", []))
    d = detector or (dets[0] if dets else None)
    out = []
    if d is not None:
        n = int(s["counts"].get(d, 0))
        out.append(Explanation(
            "run_rate", f"Measured rate of {d}", n / real, "1/s", "R = N / t_real",
            f"R = {n} / {real:.6g} s = {_g(n / real)} /s",
            f"How many particles {d} counted per second in the part of run {run.number} simulated event by event, "
            "where each event counts once.", "Unweighted events: each one stands for one reaction in the beam time; "
            "particles below the threshold are not recorded.", DOCS + "physics-register/rates-and-events.html",
            {"n": n, "t": real}, lambda n, t: n / t))
        out.append(Explanation(
            "run_counts", f"Counts in {d} over the whole run", n * duration / real, "", "N_run = N_real × t_run / t_real",
            f"N_run = {n} × {duration:.6g} s / {real:.6g} s = {_g(n * duration / real)}",
            f"The real part scaled to the run's {_time(duration)}; its histograms are redrawn with the full run's Poisson "
            "fluctuations.", "The beam and the setup stay the same over the whole run.",
            DOCS + "physics-register/rates-and-events.html", {"n": n, "t_run": duration, "t": real},
            lambda n, t_run, t: n * t_run / t))
        dead = _seconds(exp.run.dead_time)
        rate = n / real
        out.append(Explanation(
            "run_busy", f"Busy fraction of {d}", rate * dead / (1 + rate * dead), "", "b = R τ / (1 + R τ)",
            f"b = {_g(rate)} /s × {dead:.3g} s / (1 + {_g(rate)} /s × {dead:.3g} s) = {_g(rate * dead / (1 + rate * dead))}",
            "The share of the time this detector's electronics are busy and miss a particle.",
            "Non-paralysable dead time τ per count, the same for every count.", DOCS + "theory/coulex.html",
            {"r": rate, "tau": dead}, lambda r, tau: r * tau / (1 + r * tau)))
        peak = s.get("peak_coincidences", {})
        k = sum(v for key, v in peak.items() if key.split("|")[0] == d)
        if s.get("gammas"):
            out.append(Explanation(
                "run_coincidences", f"Coincidences in the full-energy peak with {d}", k / real, "1/s",
                "R_γ = N_peak / t_real", f"R_γ = {k} / {real:.6g} s = {_g(k / real)} /s",
                f"γ rays that left their full energy in a crystal while {d} counted the particle of the same event.",
                "Add-back and suppression as the setup has them; random coincidences are not in this count.",
                DOCS + "theory/coulex.html", {"k": k, "t": real}, lambda k, t: k / t))
    live = s.get("live_fraction")
    if live is not None:
        total = sum(s["counts"].values()) / real
        out.append(Explanation(
            "run_live", "Live fraction of the run", live, "", "L = 1 / (1 + τ Σ R)",
            f"L = 1 / (1 + {_seconds(exp.run.dead_time):.3g} s × (Σ R)) = {_g(live)}",
            "The share of the run the acquisition could take an event; every count is scaled by it.",
            "Σ R includes the particle singles and the γ-ray singles from the rates.", DOCS + "theory/coulex.html",
            {"live": live, "particles": total}, lambda live, particles: live))
    ev = None
    for g in planner.gates():
        if ev is None:
            ev = run.events()
        passed = int(g.mask(ev, exp).sum())
        total = len(ev)
        out.append(Explanation(
            f"gate_{g.name}", f"Gate “{g.name}”", passed / total if total else 0.0, "", "f = N_gate / N_all",
            f"f = {passed} / {total} = {_g(passed / total if total else 0.0)}",
            f"The share of the run's counted particles the gate keeps: {g.describe()}.",
            "A group is chosen by energy against the kinematic lines, as an experimentalist would, not by the "
            "simulation's knowledge of the reaction.", DOCS + "planner-app.html",
            {"n": passed, "total": total}, lambda n, total: n / total if total else 0.0))
    for r in planner.compare(run)["rows"]:
        if r["measured"] is None or not r["uncertainty"] or r["predicted"] is None:
            continue
        pull = (r["measured"] - r["predicted"]) / r["uncertainty"]
        out.append(Explanation(
            f"compare_{r['quantity']}_{r['what']}".replace(" ", "_"), f"{r['quantity']}, {r['what']}: measured "
            "against predicted", pull, "σ", "z = (measured − predicted) / σ_measured",
            f"z = ({_g(r['measured'])} − {_g(r['predicted'])}) / {_g(r['uncertainty'])} = {pull:+.2f}",
            "How far the run is from the Plan's prediction, in standard deviations of the run's statistics; within "
            "±3 they agree.", "The prediction is the Plan's at the moment the run was started.",
            DOCS + "physics-register/rates-and-events.html",
            {"m": r["measured"], "p": r["predicted"], "s": r["uncertainty"]}, lambda m, p, s: (m - p) / s))
    return out


def _seconds(q) -> float:
    if q is None:
        return 0.0
    from .quantity import Quantity

    return Quantity.parse(q).to("s") if isinstance(q, str) else float(q.to("s"))


# -- the report of the experiment ------------------------------------------------------------------------------------

def setup_history(planner) -> list:
    """How the setup changed from run to run, and since the last run: [{"run", "label", "changes"}], the first
    run's changes being those from... nothing (its setup is the start)."""
    out, previous = [], None
    for s in planner.runs():
        from .experiment import Experiment

        setup = Experiment.load(Path(s["_folder"]) / "setup.toml").to_dict()
        if previous is not None:
            from .runs import physics_changes

            out.append({"run": s["number"], "label": s["label"], "changes": physics_changes(previous, setup)})
        else:
            out.append({"run": s["number"], "label": s["label"], "changes": []})
        previous = setup
    if previous is not None:
        from .runs import physics_changes

        out.append({"run": None, "label": "the current setup",
                    "changes": physics_changes(previous, planner.experiment.to_dict())})
    return out


class _CachedGate:
    """A gate whose mask over a run's particles is worked out once: the report draws many spectra under it."""

    def __init__(self, gate):
        self.gate = gate
        self._masks: dict = {}

    @property
    def name(self) -> str:
        return self.gate.name

    @property
    def detectors(self):
        return self.gate.detectors

    def describe(self) -> str:
        return self.gate.describe()

    def mask(self, events, experiment) -> np.ndarray:
        key = id(events)
        if key not in self._masks:
            self._masks[key] = (events, self.gate.mask(events, experiment))
        return self._masks[key][1]


#: The width of a bin of the report's γ-ray spectra, keV.
GAMMA_BIN_KEV = 4.0


def gamma_lines(experiment) -> list:
    """The γ rays the excited nucleus can emit, for marking on a corrected spectrum: [{"energy_kev", "label"}],
    from its level scheme (E0 transitions emit none), or the one transition of the reaction without a scheme."""
    exc = experiment.excitation
    if exc is None:
        return []
    scheme = experiment.levels.get("target" if exc.excite == "target" else "beam")
    if scheme is None or not scheme.transitions:
        return [{"energy_kev": 1e3 * exc.energy_mev, "label": "the transition"}]
    out = []
    for t in scheme.transitions:
        if t.multipolarity.strip() == "E0":
            continue
        up, low = scheme.levels[t.initial], scheme.levels[t.final]
        out.append({"energy_kev": float(t.energy.value), "label": f"{up.label} → {low.label}",
                    "from_kev": float(up.energy.value)})
    return sorted(out, key=lambda x: x["energy_kev"])


def _correction_of(run, result) -> str:
    settings = dict(getattr(result, "settings", {}) or {})
    correction = settings.get("correction", "emitter")
    if correction not in ("projectile", "recoil"):
        exc = run.experiment.excitation
        correction = "recoil" if exc is not None and exc.excite == "target" else "projectile"
    return correction


def _gamma_set(run, correction: str, gate=None) -> dict:
    """The γ-ray spectra of a run under one gate: {"detectors": {name: {"raw", "corrected"}}, "total": {...}},
    each a spectrum of :func:`~physim.nuclear.dataviews.gamma_spectrum` to the top of the level scheme."""
    from . import dataviews as dv

    g = run.gammas()
    lo, hi = dv.gamma_range(g, run.experiment)
    bins = max(int(round(1e3 * (hi - lo) / GAMMA_BIN_KEV)), 50)
    out = {"detectors": {}, "total": None}
    for n in g.detector_names():
        out["detectors"][n] = {
            "raw": dv.gamma_spectrum(run, n, correction="off", gate=gate, randoms="none", bins=bins, range=(lo, hi)),
            "corrected": dv.gamma_spectrum(run, n, correction=correction, gate=gate, randoms="none", bins=bins,
                                           range=(lo, hi))}
    if out["detectors"]:
        first = next(iter(out["detectors"].values()))
        out["total"] = {k: {"edges": first[k]["edges"],
                            "counts": np.sum([s[k]["counts"] for s in out["detectors"].values()], axis=0)}
                        for k in ("raw", "corrected")}
    return out


def run_views(planner, run, result=None) -> dict:
    """What the report draws of a run: the spectra of every particle detector and, for every γ-ray detector, the
    raw spectrum and the Doppler-corrected one (:mod:`physim.nuclear.dataviews`), to the top of the level scheme,
    with the correction and the gate the analysis used (``result``) or, without an analysis, the correction for
    the nucleus that was excited and no gate. Returns {"particles": {name: spectrum}, "gammas": {name: {"raw":
    spectrum, "corrected": spectrum}}, "total": the γ-ray detectors summed, "particles_all": the particle
    detectors' spectra without the gate (None without a gate), "correction", "gate", "source": the crystals'
    spectra of a source run}."""
    from . import dataviews as dv

    if run.kind == "source":
        sr = run.source()
        return {"particles": {}, "particles_all": None, "gammas": {}, "total": None, "correction": None,
                "gate": None, "source": {n: {"counts": sr.spectra[n], "edges": sr.edges} for n in sr.names()}}
    settings = dict(getattr(result, "settings", {}) or {})
    correction = _correction_of(run, result)
    gate = _gate_of(planner, settings)
    cached = _CachedGate(gate) if gate is not None else None
    ev = run.events()
    out = {"particles": {d: dv.particle_spectrum(run, d, gate=cached) for d in ev.detectors}, "gammas": {},
           "particles_all": ({d: dv.particle_spectrum(run, d) for d in ev.detectors} if gate is not None else None),
           "total": None, "correction": correction, "gate": gate, "source": None}
    if run.gammas() is not None:
        both = _gamma_set(run, correction, cached)
        out["gammas"], out["total"] = both["detectors"], both["total"]
    return out


def _gate_of(planner, settings: dict):
    """The named gate of the Data tab an analysis was run with, if its settings are one gate's (the analysis keeps
    the gate's detectors, rings, group and energy window, not its name)."""
    def norm(v):
        return list(v) if isinstance(v, (list, tuple)) else v

    for g in planner.gates():
        gs = g.settings()
        if all(norm(settings.get(k)) == norm(getattr(gs, k))
               for k in ("detectors", "rings", "particle_gate", "particle_energy")):
            return g
    return None


def _png(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()


def _gate_words(gate) -> str:
    return (f"the particles of gate “{gate.name}” ({gate.describe()})" if gate is not None
            else "every counted particle")


def run_spectra_png(planner, run, result=None, dpi: int = 110, views: Optional[dict] = None) -> tuple:
    """The spectra of a run's detectors as one PNG (bytes), and its caption: a panel per particle detector with
    its measured-energy spectrum (a source run: a panel per crystal). The γ-ray spectra of a beam run have their
    own figure (:func:`run_gamma_png`). Needs matplotlib."""
    from matplotlib.figure import Figure

    v = views if views is not None else run_views(planner, run, result)
    panels = []
    if v["source"]:
        for n, s in v["source"].items():
            panels.append((n, [("", s["edges"], s["counts"])], "γ-ray energy (keV)"))
    else:
        everything = v.get("particles_all") or v["particles"]
        for n, s in everything.items():
            lines = [("every counted particle", s["edges"], s["counts"])]
            if v["gate"] is not None:
                lines.append((f"gate “{v['gate'].name}”", v["particles"][n]["edges"], v["particles"][n]["counts"]))
            panels.append((n, lines, "energy (MeV)"))
    if not panels:
        return b"", ""
    cols = min(3, len(panels))
    rows = -(-len(panels) // cols)
    fig = Figure(figsize=(3.4 * cols, 2.4 * rows), dpi=dpi)
    axes = fig.subplots(rows, cols, squeeze=False)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for k, (ax, (title, lines, xlabel)) in enumerate(zip(axes.flat, panels)):
        populated = []
        for i, (label, edges, counts) in enumerate(lines):
            y, edges = np.asarray(counts, dtype=float), np.asarray(edges, dtype=float)
            if not np.any(y > 0):
                continue
            ax.stairs(np.where(y > 0, y, np.nan), edges, label=label, color=("#0072B2", "#D55E00")[i % 2],
                      lw=0.9, fill=i > 0, alpha=1.0 if i == 0 else 0.55)
            populated.append(edges[:-1][y > 0])
        if populated:
            ax.set_yscale("log")
            if not v["source"]:  # the particles fill a small part of the range from zero
                lo, hi = min(float(p.min()) for p in populated), max(float(p.max()) for p in populated)
                ax.set_xlim(lo - 0.15 * (hi - lo) - 1e-9, hi + 0.15 * (hi - lo) + 1e-9)
        else:
            ax.text(0.5, 0.5, "no counts", transform=ax.transAxes, ha="center", va="center", fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("counts", fontsize=8)
        ax.tick_params(labelsize=7)
        if k == 0 and len(lines) > 1:
            ax.legend(fontsize=6, frameon=False, loc="lower center")
    fig.tight_layout()
    if v["source"]:
        caption = f"Run {run.number}: the counts each crystal recorded from the source, over {_time(run.duration_s)}."
    else:
        caption = (f"Run {run.number}: the measured energy in each particle detector, counts per bin in the real "
                   f"part of the run ({_time(run.real_s)}), for every counted particle"
                   + (f", and filled, for {_gate_words(v['gate'])}, which the analysis took" if v["gate"] else "")
                   + ".")
    return _png(fig), caption


def run_gamma_png(planner, run, result=None, dpi: int = 110, views: Optional[dict] = None) -> tuple:
    """The γ-ray spectra of a beam run as one PNG (bytes), and its caption: all γ-ray detectors summed, then each
    detector, as measured and Doppler-corrected over each other, to the top of the level scheme, with the γ
    rays of the scheme marked on the corrected spectrum (:func:`gamma_lines`). (b"", "") without γ rays."""
    from matplotlib.figure import Figure

    v = views if views is not None else run_views(planner, run, result)
    if not v["gammas"]:
        return b"", ""
    lines = gamma_lines(run.experiment)
    panels = [("All γ-ray detectors", v["total"])] + list(v["gammas"].items())
    cols = 1 if len(panels) <= 5 else 2  # five panels in a column fill a page
    rows = -(-len(panels) // cols)
    fig = Figure(figsize=(7.6, 1.62 * rows + 0.75), dpi=dpi)
    axes = fig.subplots(rows, cols, squeeze=False, sharex=True)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    x0, x1 = 1e3 * float(v["total"]["raw"]["edges"][0]), 1e3 * float(v["total"]["raw"]["edges"][-1])
    for k, (ax, (title, s)) in enumerate(zip(axes.flat, panels)):
        top = 1.0
        for i, (key, label) in enumerate((("raw", "as measured"), ("corrected", f"corrected for the {v['correction']}"))):
            y = np.asarray(s[key]["counts"], dtype=float)
            top = max(top, float(y.max()) if len(y) else 1.0)
            ax.stairs(np.where(y > 0, y, np.nan), 1e3 * np.asarray(s[key]["edges"], dtype=float), label=label,
                      color=("#0072B2", "#D55E00")[i], lw=0.8 if i == 0 else 1.1)
        ax.set_yscale("log")
        ax.set_ylim(0.5, top * 4)
        level, last = 0, -1e9
        for ln in lines:
            ax.axvline(ln["energy_kev"], color="#555", lw=0.5, ls=":")
            if k < cols:  # the energies of the lines, above the first row; close ones on two levels
                level = (level + 1) % 3 if ln["energy_kev"] - last < 0.03 * (x1 - x0) else 0
                last = ln["energy_kev"]
                ax.text(ln["energy_kev"], 1.02 + 0.13 * level, f"{ln['energy_kev']:.0f}", fontsize=6, ha="center",
                        va="bottom", color="#333", transform=ax.get_xaxis_transform(), clip_on=False)
        ax.set_xlim(x0, x1)
        ax.text(0.995, 0.93, title, transform=ax.transAxes, ha="right", va="top", fontsize=8)
        ax.set_ylabel(f"counts / {GAMMA_BIN_KEV:g} keV", fontsize=7)
        ax.tick_params(labelsize=7)
        if k == 0:
            ax.legend(fontsize=7, frameon=False, loc="upper left", ncols=2)
    for ax in axes[-1]:
        ax.set_xlabel("γ-ray energy (keV)", fontsize=8)
    fig.tight_layout()
    fig.subplots_adjust(top=1 - 0.5 / fig.get_figheight(), hspace=0.12)
    marked = (" The dotted lines, with their energies in keV above the top panel, mark the γ rays of the level "
              "scheme: "
              + ", ".join(f"{ln['energy_kev']:.1f} keV ({ln['label']})" for ln in lines) + "." if lines else "")
    caption = (f"Run {run.number}: the γ-ray spectra in the real part of the run ({_time(run.real_s)}), as measured "
               f"and Doppler-corrected for the {v['correction']}, in coincidence with {_gate_words(v['gate'])}; "
               "random coincidences not added; add-back and suppression as the setup has them. The correction "
               "gathers the line that the emitter's motion spreads over the measured spectrum." + marked)
    return _png(fig), caption


def run_rings_png(planner, run, dpi: int = 110) -> tuple:
    """The measured energy against ring (or strip) of every particle detector of a beam run as one PNG, with the
    kinematic line of each group and the window of every gate that has one drawn on it, and its caption."""
    from matplotlib.colors import LogNorm
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    from . import dataviews as dv

    ev = run.events()
    gates = list(planner.gates())
    names = list(ev.detectors)
    if not names:
        return b"", ""
    cols = min(3, len(names))
    rows = -(-len(names) // cols)
    fig = Figure(figsize=(3.6 * cols, 3.0 * rows), dpi=dpi)
    axes = fig.subplots(rows, cols, squeeze=False)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    colours = ("#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
    drawn, undrawn, outside = [], [], []
    for ax, name in zip(axes.flat, names):
        view = dv.energy_vs_ring(run, name)
        h = np.asarray(view["counts"], dtype=float).T
        n_rings = len(view["rings"])
        edges = np.asarray(view["edges"], dtype=float)
        lo, hi = float(edges[0]), float(edges[-1])
        if h.max() > 0:
            ax.pcolormesh(np.arange(n_rings + 1) + 0.5, edges, np.where(h > 0, h, np.nan),
                          norm=LogNorm(vmin=max(1.0, h.max() / 300), vmax=max(h.max(), 2)), cmap="Blues",
                          rasterized=True)
            filled = np.flatnonzero(h.sum(axis=1) > 0)
            lo, hi = float(edges[filled[0]]), float(edges[filled[-1] + 1])
            lo, hi = lo - 0.12 * (hi - lo), hi + 0.12 * (hi - lo)
        for ln in view["lines"]:
            ax.plot(np.asarray(ln["rings"]) + 1, ln["energy"], lw=0.9,
                    color="#000000" if ln["group"] == "elastic" else "#D55E00",
                    ls="-" if ln["group"] == "elastic" else "--")
        for i, gate in enumerate(gates):
            if gate.detectors is not None and name not in gate.detectors:
                continue
            if gate.rings is None and gate.energy is None:
                continue
            r0, r1 = (gate.rings[0], gate.rings[-1]) if gate.rings is not None else (0, n_rings - 1)
            e0, e1 = gate.energy if gate.energy is not None else (lo + 0.02 * (hi - lo), hi - 0.02 * (hi - lo))
            if e1 < lo or e0 > hi:  # its window lies outside the energies the run recorded: said, not drawn
                if gate.name not in outside:
                    outside.append(gate.name)
                continue
            ax.add_patch(Rectangle((r0 + 0.5, e0), r1 - r0 + 1, e1 - e0, fill=False, lw=1.2,
                                   ec=colours[i % len(colours)]))
            ax.text(r0 + 0.6, min(e1, hi - 0.01 * (hi - lo)), gate.name, fontsize=6, color=colours[i % len(colours)],
                    va="top" if e1 >= hi - 0.05 * (hi - lo) else "bottom")
            if gate.name not in drawn:
                drawn.append(gate.name)
        ax.set_title(name, fontsize=9)
        ax.set_xlabel("ring or strip", fontsize=8)
        ax.set_ylabel("measured energy (MeV)", fontsize=8)
        ax.set_xlim(0.5, n_rings + 0.5)
        ax.set_ylim(lo, hi)
        ax.tick_params(labelsize=7)
    undrawn = [g.name for g in gates if g.name not in drawn and g.name not in outside]
    fig.tight_layout()
    caption = (f"Run {run.number}: the measured energy against ring (or strip) in each particle detector, on a "
               "logarithmic scale of counts, over the energies the run recorded. The solid black lines are where "
               "elastic scattering puts each group, the dashed orange ones where the excited group falls.")
    if drawn:
        caption += " The boxes are the gates with a window: " + ", ".join(f"“{n}”" for n in drawn) + "."
    if outside:
        caption += (" The window of " + ", ".join(f"“{n}”" for n in outside) + " lies outside these energies: "
                    "no particle of the run passes it.")
    if undrawn:
        caption += (" Not drawn as a box (they take whole detectors or a group by its kinematic line): "
                    + ", ".join(f"“{n}”" for n in undrawn) + ".")
    return _png(fig), caption


def run_gates_png(planner, run, result=None, dpi: int = 110) -> tuple:
    """What each gate does to the γ-ray spectrum of a beam run: the Doppler-corrected spectrum of all γ-ray
    detectors with every counted particle and under each named gate, as one PNG, its caption, and the rows
    [{"gate", "what", "particles", "gammas", "in_peak"}] of the table beside it ("in_peak": γ rays within 1.5%
    of the transition energy after the correction). (b"", "", []) without γ rays."""
    from matplotlib.figure import Figure

    g = run.gammas()
    if g is None:
        return b"", "", []
    exp = run.experiment
    correction = _correction_of(run, result)
    ev = run.events()
    e0 = 1e3 * g.energy_mev
    sets = [("every counted particle", None)] + [(gt.name, _CachedGate(gt)) for gt in planner.gates()]
    fig = Figure(figsize=(7.6, 3.4), dpi=dpi)
    ax = fig.subplots()
    colours = ("#000000", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#0072B2")
    rows, top = [], 1.0
    for i, (label, gate) in enumerate(sets):
        s = _gamma_set(run, correction, gate)["total"]["corrected"]
        y, edges = np.asarray(s["counts"], dtype=float), 1e3 * np.asarray(s["edges"], dtype=float)
        top = max(top, float(y.max()))
        ax.stairs(np.where(y > 0, y, np.nan), edges, label=label, color=colours[i % len(colours)],
                  lw=1.0 if i == 0 else 0.9)
        centres = (edges[:-1] + edges[1:]) / 2
        passed = ev["counted"] if gate is None else (gate.mask(ev, exp) & ev["counted"])
        rows.append({"gate": label, "what": "" if gate is None else gate.describe(),
                     "particles": int(np.sum(passed)), "gammas": float(y.sum()),
                     "in_peak": float(y[np.abs(centres - e0) <= 0.015 * e0].sum())})
    ax.axvspan(0.985 * e0, 1.015 * e0, color="#999", alpha=0.15, lw=0)
    ax.set_yscale("log")
    ax.set_ylim(0.5, top * 3)
    ax.set_xlabel("γ-ray energy (keV)", fontsize=8)
    ax.set_ylabel(f"counts / {GAMMA_BIN_KEV:g} keV", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, frameon=False, ncols=2)
    fig.tight_layout()
    caption = (f"Run {run.number}: the γ-ray spectrum of all detectors, Doppler-corrected for the {correction}, in "
               "coincidence with every counted particle and with the particles of each gate. The grey band is "
               f"{e0:.0f} keV ± 1.5%, the window of the last column of the table.")
    return _png(fig), caption, rows


def _figure(png: bytes, caption: str, alt: str) -> str:
    return _img(png, alt, caption) if png else ""


def _run_pieces(planner, s: dict, result, gated: bool) -> dict:
    """What the report takes from one run, read once (a long run is gigabytes in memory and half a minute to
    read): {"figures": its detectors' spectra, the energy against ring with the gates, and the γ-ray spectra as
    measured and Doppler-corrected, as HTML, or a line saying why not; "gates": the γ-ray spectrum under each
    gate with its table (``gated``), as HTML; "seen": {"correction", "widths"} for the Doppler section}."""
    def e(text) -> str:
        return html.escape(str(text), quote=False)

    n = s["number"]
    out = {"figures": "", "gates": "", "seen": None}
    try:
        run = planner.run if planner.run is not None and planner.run.number == n else planner.folder.load_run(n)
    except Exception as err:  # noqa: BLE001 - the report stands without the run
        out["figures"] = f"<p>The run could not be read: {e(err)}.</p>"
        return out
    try:
        views = run_views(planner, run, result)
        png, caption = run_spectra_png(planner, run, result, views=views)
        parts = [_figure(png, caption, f"Spectra of run {n}")]
        if run.kind != "source":
            png, caption = run_rings_png(planner, run)
            parts.append(_figure(png, caption, f"Energy against ring of run {n}"))
            png, caption = run_gamma_png(planner, run, result, views=views)
            parts.append(_figure(png, caption, f"γ-ray spectra of run {n}"))
        out["figures"] = "".join(parts)
    except ImportError:
        out["figures"] = "<p>The spectra are not drawn: matplotlib is not installed.</p>"
    except Exception as err:  # noqa: BLE001 - the report stands without the figures
        out["figures"] = f"<p>The spectra could not be drawn: {e(err)}.</p>"
    if run.kind == "source":
        return out
    if gated:
        try:
            png, caption, rows = run_gates_png(planner, run, result)
            if png:
                out["gates"] = (f"<h3>Run {n} under each gate</h3>" + _figure(png, caption, f"Run {n} under each gate")
                                + "<table><thead><tr><th>Gate</th><th>Particles that pass</th>"
                                "<th>γ rays in coincidence</th><th>Of these, at the transition</th></tr></thead>"
                                "<tbody>"
                                + "".join(f"<tr><td>{e(r['gate'])}</td><td>{r['particles']:,}</td>"
                                          f"<td>{r['gammas']:,.0f}</td><td>{r['in_peak']:,.0f}</td></tr>"
                                          for r in rows) + "</tbody></table>")
        except ImportError:
            out["gates"] = "<p>The gated spectra are not drawn: matplotlib is not installed.</p>"
        except Exception as err:  # noqa: BLE001 - the report stands without the figure
            out["gates"] = f"<p>Run {n}: the gated spectra could not be drawn: {e(err)}.</p>"
    try:
        from .alignment import peak_widths

        g = run.gammas()
        out["seen"] = {"correction": _correction_of(run, result),
                       "widths": peak_widths(g, run.experiment) if g is not None else {}}
    except Exception:  # noqa: BLE001 - the Doppler section then reads the run itself
        out["seen"] = None
    return out


def _img(png: bytes, alt: str, caption: str) -> str:
    return (f'<figure><img alt="{html.escape(alt)}" style="max-width:100%" '
            f'src="data:image/png;base64,{base64.b64encode(png).decode("ascii")}"><figcaption>'
            f"{html.escape(caption, quote=False)}</figcaption></figure>")


def level_schemes(planner) -> list:
    """The level schemes the report draws: [{"role", "table", "source"}], the setup's own (from ENSDF or the
    user, :meth:`~physim.nuclear.planner.Planner.levels`) and, when the nucleus that is excited has none, the
    two levels the reaction is set up with (the ground state and the excited state)."""
    info = planner.levels()
    out = [{"role": role, "table": info[role], "source": "setup"} for role in ("target", "beam") if info.get(role)]
    exc = planner.experiment.excitation
    if exc is not None:
        role = "target" if exc.excite == "target" else "beam"
        if not info.get(role):
            from .analysis import ex_nuclide

            e = 1e3 * exc.energy_mev
            out.insert(0, {"role": role, "source": "reaction", "table": {
                "nuclide": ex_nuclide(planner.experiment), "reference": "", "notes": [],
                "levels": [{"n": 0, "energy_kev": 0.0, "jpi": "ground state"},
                           {"n": 1, "energy_kev": e, "jpi": "excited state"}],
                "transitions": [{"from": 1, "to": 0, "energy_kev": e, "multipolarity": exc.multipolarity,
                                 "branching": 1.0}], "matrix_elements": []}})
    return out


def level_scheme_png(table: dict, mark_kev: Optional[float] = None, dpi: int = 130, most: int = 24) -> bytes:
    """A level scheme as a PNG: a line per level with its spin and parity on the left and its energy on the
    right, an arrow per γ-ray transition with its energy (the ``most`` strongest branches if there are more), and
    the level at ``mark_kev`` (the state the reaction excites) drawn heavier. Needs matplotlib."""
    from matplotlib.figure import Figure

    levels, transitions = table["levels"], list(table["transitions"])
    if len(transitions) > most:
        transitions = sorted(transitions, key=lambda t: -(t["branching"] or 0.0))[:most]
        transitions.sort(key=lambda t: (t["from"], t["to"]))
    n = max(len(transitions), 1)
    top = max((lv["energy_kev"] for lv in levels), default=1.0) or 1.0
    fig = Figure(figsize=(min(max(3.6, 1.6 + 0.42 * n), 9.0), 3.8), dpi=dpi)
    ax = fig.add_subplot()
    pad = max(1.8, 0.25 * (n + 1))  # room beside the lines for the spins (left) and the energies (right)
    last = None
    for lv in levels:
        e = lv["energy_kev"]
        marked = mark_kev is not None and abs(e - mark_kev) < 0.5
        ax.plot([0, n + 1], [e, e], color="#0072B2", lw=2.2 if marked else 1.1)
        # Labels of levels closer than the text is high are moved apart.
        y = e if last is None or e - last > 0.045 * top else last + 0.045 * top
        last = y
        ax.text(-0.04 * pad, y, lv["jpi"] or "", ha="right", va="center", fontsize=7)
        ax.text(n + 1 + 0.04 * pad, y, f"{e:g} keV", ha="left", va="center", fontsize=7)
    for k, t in enumerate(transitions, start=1):
        e_top, e_bottom = levels[t["from"]]["energy_kev"], levels[t["to"]]["energy_kev"]
        ax.annotate("", xy=(k, e_bottom), xytext=(k, e_top),
                    arrowprops=dict(arrowstyle="-|>", color="#D55E00", lw=1.0, shrinkA=0, shrinkB=0))
        ax.text(k, e_top + 0.012 * top, f"{t['energy_kev']:.0f}", rotation=90, ha="center", va="bottom", fontsize=6,
                color="#D55E00")
    ax.set_xlim(-0.55 * pad, n + 1 + 0.45 * pad)
    ax.set_ylim(-0.04 * top, 1.16 * top)
    ax.set_axis_off()  # every level carries its energy, so no axis
    ax.set_title(table.get("nuclide") or "", fontsize=9)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.03)
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()


def level_schemes_html(planner) -> str:
    """The level schemes as figures with their captions, for the reports; empty without any."""
    exc = planner.experiment.excitation
    parts = []
    for sc in level_schemes(planner):
        t = sc["table"]
        mine = exc is not None and sc["role"] == ("target" if exc.excite == "target" else "beam")
        mark = 1e3 * exc.energy_mev if mine else None
        if sc["source"] == "reaction":
            lam = exc.multipolarity
            caption = (f"{t['nuclide']} (the {sc['role']}): the two levels the reaction is set up with, the ground "
                       f"state and the state at {mark:g} keV that is excited ({lam}, B({lam}↑) = "
                       f"{exc.b_up_e2fm:.4g} e²fm^{2 * int(lam[1])}). No level scheme has been looked up for this "
                       "nucleus (Physics, in the app), so no other level is drawn.")
        else:
            shown = min(len(t["transitions"]), 24)
            caption = (f"{t['nuclide']} (the {sc['role']}): {len(t['levels'])} levels and {shown} "
                       + (f"of its {len(t['transitions'])} " if shown < len(t["transitions"]) else "")
                       + "γ-ray transitions, labelled with their energies in keV"
                       + (f"; the state at {mark:g} keV that the reaction excites is drawn heavier" if mine else "")
                       + "." + (f" {t['reference']}" if t.get("reference") else ""))
        try:
            parts.append(_img(level_scheme_png(t, mark), f"Level scheme of {t['nuclide']}", caption))
        except ImportError:
            parts.append("<p>The level scheme is not drawn: matplotlib is not installed.</p>")
            break
        except Exception as err:  # noqa: BLE001 - the report stands without the figure
            parts.append(f"<p>The level scheme could not be drawn: {html.escape(str(err), quote=False)}.</p>")
    return "".join(parts)


def analysis_png(result, dpi: int = 130) -> bytes:
    """The spectrum an analysis fitted, with the fit: the gated, Doppler-corrected counts, the Gaussian on its
    line within the fit window, the ±3σ band the area is counted in, and the random coincidences expected.
    Empty for a result without its spectrum (one made before spectra were kept). Needs matplotlib."""
    sp = getattr(result, "spectrum", None)
    if not sp:
        return b""
    from matplotlib.figure import Figure

    from .response import FWHM_PER_SIGMA

    edges, counts = 1e3 * np.asarray(sp["edges"], dtype=float), np.asarray(sp["counts"], dtype=float)
    rnd = np.asarray(sp.get("randoms") or [], dtype=float)
    # The fit works on bins a quarter of the resolution wide; drawn four to a bin, a peak of few counts shows.
    k = 4 if len(counts) >= 40 else 1
    m = len(counts) // k * k
    counts, edges = counts[:m].reshape(-1, k).sum(axis=1), edges[:m + 1:k]
    rnd = rnd[:m].reshape(-1, k).sum(axis=1) if rnd.size >= m else rnd[:0]
    fit = result.fit
    fig = Figure(figsize=(6.4, 3.2), dpi=dpi)
    ax = fig.add_subplot()
    ax.stairs(counts, edges, color="#0072B2", lw=0.9, label=f"corrected for the {sp.get('correction', 'emitter')}")
    if rnd.size and rnd.any():
        ax.stairs(rnd, edges, color="#888888", lw=0.8, ls="--", label="random coincidences expected")
    if fit.sigma > 0 and math.isfinite(fit.centroid) and (fit.height > 0 or fit.level > 0):
        half = float(result.settings.get("fit_half_width", 4.0)) * FWHM_PER_SIGMA * fit.sigma
        x = np.linspace(fit.centroid - half, fit.centroid + half, 400)
        ax.plot(1e3 * x, k * fit.curve(x), color="#D55E00", lw=1.3, label="fit: Gaussian on a line")
        ax.axvspan(1e3 * (fit.centroid - 3 * fit.sigma), 1e3 * (fit.centroid + 3 * fit.sigma), color="#D55E00",
                   alpha=0.10, lw=0, label="±3σ: the area")
    if sp.get("energy"):
        ax.axvline(1e3 * sp["energy"], color="#444444", lw=0.7, ls=":", label=f"{1e3 * sp['energy']:.1f} keV")
    ax.set_xlabel("corrected γ-ray energy (keV)", fontsize=8)
    ax.set_ylabel(f"counts per {edges[1] - edges[0]:.2g} keV in the run", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()


def search_png(result, dpi: int = 130) -> bytes:
    """The peak search of an analysis: the spectrum corrected for the projectile and the one corrected for the
    recoil, each over its whole range, with the peaks found marked and those on a transition named. Empty for a
    result without a search (one made before the peaks were searched for). Needs matplotlib."""
    search = getattr(result, "search", None)
    if not search:
        return b""
    from matplotlib.figure import Figure

    fig = Figure(figsize=(6.6, 5.4), dpi=dpi)
    for k, which in enumerate(("projectile", "recoil")):
        sp = search[which]
        ax = fig.add_subplot(2, 1, k + 1)
        edges, counts = 1e3 * np.asarray(sp["edges"], dtype=float), np.asarray(sp["counts"], dtype=float)
        m = len(counts) // 2 * 2
        counts, edges = counts[:m].reshape(-1, 2).sum(axis=1), edges[:m + 1:2]
        ax.stairs(np.where(counts > 0, counts, np.nan), edges, color="#0072B2", lw=0.7)
        centres = (edges[:-1] + edges[1:]) / 2
        top = max(float(np.nanmax(counts)) if len(counts) else 1.0, 1.0)
        for pk in sp["peaks"]:
            near = np.abs(centres - pk["energy_kev"]) <= max(pk["fwhm_kev"], edges[1] - edges[0])
            y = max(float(np.nanmax(counts[near])) if near.any() else 1.0, 1.0)
            ax.plot([pk["energy_kev"]], [1.7 * y], marker="v", ms=4, color="#D55E00" if pk["label"] else "#888888",
                    ls="none")
            if pk["label"]:
                name = pk["label"].split(" ", 1)[1] if pk.get("nuclide") and " " in pk["label"] else pk["label"]
                ax.annotate(name, (pk["energy_kev"], 2.4 * y), rotation=55, fontsize=6.5, ha="left", va="bottom")
        ax.set_yscale("log")
        ax.set_ylim(0.6, 60 * top)
        ax.set_xlim(0, edges[-1])
        ax.set_title(f"corrected for the {which}" + (f": the lines of {sp['nuclide']}" if sp.get("nuclide") else ""),
                     fontsize=8, loc="left")
        ax.set_ylabel("counts per bin", fontsize=8)
        ax.tick_params(labelsize=7)
        if k == 1:
            ax.set_xlabel("corrected γ-ray energy (keV)", fontsize=8)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()


def line_png(line: dict, fit_half_width: float = 4.0, dpi: int = 130) -> bytes:
    """One integrated peak (an entry of an analysis result's ``lines``), as :func:`analysis_png` draws the first:
    the gated, corrected counts, the Gaussian on its line, the ±3σ band of the area and the random coincidences
    expected. Needs matplotlib."""
    from matplotlib.figure import Figure

    from .response import FWHM_PER_SIGMA

    sp, fit = line["spectrum"], line["fit"]
    edges, counts = 1e3 * np.asarray(sp["edges"], dtype=float), np.asarray(sp["counts"], dtype=float)
    rnd = np.asarray(sp.get("randoms") or [], dtype=float)
    k = 4 if len(counts) >= 40 else 1
    m = len(counts) // k * k
    counts, edges = counts[:m].reshape(-1, k).sum(axis=1), edges[:m + 1:k]
    rnd = rnd[:m].reshape(-1, k).sum(axis=1) if rnd.size >= m else rnd[:0]
    fig = Figure(figsize=(6.4, 3.2), dpi=dpi)
    ax = fig.add_subplot()
    ax.stairs(counts, edges, color="#0072B2", lw=0.9, label=f"corrected for the {line['correction']}")
    if rnd.size and rnd.any():
        ax.stairs(rnd, edges, color="#888888", lw=0.8, ls="--", label="random coincidences expected")
    sigma, centre = float(fit.get("sigma") or 0.0), float(fit.get("centroid") or 0.0)
    if sigma > 0 and math.isfinite(centre) and (fit.get("height", 0) > 0 or fit.get("level", 0) > 0):
        half = fit_half_width * FWHM_PER_SIGMA * sigma
        x = np.linspace(centre - half, centre + half, 400)
        curve = fit["height"] * np.exp(-0.5 * ((x - centre) / sigma) ** 2) + fit["level"] + fit["slope"] * (
            x - fit["origin"])
        ax.plot(1e3 * x, k * curve, color="#D55E00", lw=1.3, label="fit: Gaussian on a line")
        ax.axvspan(1e3 * (centre - 3 * sigma), 1e3 * (centre + 3 * sigma), color="#D55E00", alpha=0.10, lw=0,
                   label="±3σ: the area")
    ax.axvline(line["energy_kev"], color="#444444", lw=0.7, ls=":", label=f"{line['energy_kev']:.1f} keV")
    ax.set_xlabel("corrected γ-ray energy (keV)", fontsize=8)
    ax.set_ylabel(f"counts per {edges[1] - edges[0]:.2g} keV in the run", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    return buf.getvalue()


def _search_html(result, e) -> str:
    """The peaks the search found, as a figure and one table per corrected spectrum."""
    search = getattr(result, "search", None)
    if not search:
        return ""
    parts = ["<h4>The peaks found</h4><p>Each Doppler-corrected spectrum is searched over its whole range. A "
             "nucleus's lines are sharp in the spectrum corrected for it and smeared in the other, so each nucleus "
             "shows its peaks in its own spectrum. Every peak found is fitted; one that sits on a transition of "
             "the level scheme is named.</p>"]
    try:
        png = search_png(result)
    except Exception:  # noqa: BLE001 - the tables stand without the figure
        png = b""
    if png:
        parts.append(_img(png, "The peaks found", "The spectrum corrected for the projectile (top) and for the "
                          "recoil (bottom), random coincidences subtracted, with the peaks the search found: "
                          "orange on a transition of the level scheme, grey on none."))
    for which in ("projectile", "recoil"):
        sp = search[which]
        parts.append(f"<p><b>Corrected for the {which}</b>" + (f" (the lines of {e(sp['nuclide'])})"
                                                                 if sp.get("nuclide") else "")
                     + f", in coincidence with {e(sp['gate'])}: "
                     + (f"{len(sp['peaks'])} peak{'s' if len(sp['peaks']) != 1 else ''} found.</p>"
                        if sp["peaks"] else "no peak found.</p>"))
        if sp["peaks"]:
            parts.append("<table><thead><tr><th>Peak at (keV)</th><th>FWHM (keV)</th><th>Area (counts)</th>"
                         "<th>Significance</th><th>Transition</th></tr></thead><tbody>" + "".join(
                             f"<tr><td>{pk['energy_kev']:.1f}</td><td>{pk['fwhm_kev']:.1f}</td>"
                             f"<td>{pk['area']:,.0f} ± {pk['area_unc']:,.0f}</td><td>{pk['significance']:.0f} σ</td>"
                             "<td>" + (e(f"{pk['label']} ({pk['line_kev']:.1f} keV)") if pk["label"] else
                                       (e(f"none (smeared {pk['near']})") if pk.get("near") else "none"))
                             + "</td></tr>" for pk in sp["peaks"]) + "</tbody></table>")
    return "".join(parts)


def _lines_html(result, e) -> str:
    """The integrated peaks side by side, from the counts to the yield, and the figure of each peak after the
    first (which :func:`analysis_html` draws)."""
    lines = getattr(result, "lines", None)
    if not lines:
        return ""
    rows = [("Corrected for", lambda ln: f"the {ln['correction']}"),
            ("In coincidence with", lambda ln: ln["gate"]),
            ("Found by the search", lambda ln: (f"at {ln['found']['energy_kev']:.1f} keV, "
                                                f"{ln['found']['significance']:.0f} σ") if ln.get("found")
             else "no: fitted at the transition energy"),
            ("Fitted centre and FWHM", lambda ln: f"{1e3 * ln['fit']['centroid']:.1f} keV, "
                                                  f"{1e3 * 2.3548 * ln['fit']['sigma']:.1f} keV"),
            ("Counts over the line, ±3σ", lambda ln: f"{ln['fit']['area']:,.1f} ± {ln['fit']['area_unc']:,.1f}"),
            ("Random coincidences in the window", lambda ln: f"{ln['randoms']:,.1f}"),
            ("Area", lambda ln: f"{ln['area']:,.1f} ± {ln['area_unc']:,.1f}"),
            ("Full-energy efficiency", lambda ln: f"{100 * ln['efficiency']:.3f}%"),
            ("Angular correlation", lambda ln: f"{ln['correlation']:.3f}"),
            ("Share of the excitations the gate keeps", lambda ln: f"{100 * ln['gate_efficiency']:.1f}%"),
            ("γ rays emitted (the yield)", lambda ln: f"{ln['yield']:,.4g}"),
            ("Per particle in the same rings", lambda ln: f"{ln['per_particle']:.4g}"),
            ("Expected from the level scheme, all orders", lambda ln: "" if ln.get("expected_per_particle") is None
             else f"{ln['expected_per_particle']:.4g} (measured / expected {ln['ratio']:.3f})"),
            ("Also inside the window", lambda ln: ", ".join(ln.get("blends") or []))]
    parts = ["<h4>The peaks integrated</h4><p>The peak of each nucleus in its own corrected spectrum: a Gaussian "
             "on a straight line is fitted, and the area is the counts within ±3σ of its centre over that line, "
             "less the random coincidences expected there. The area over the efficiency and the angular "
             "correlation is the number of γ rays emitted.</p>"]
    for ln in lines[1:]:
        try:
            png = line_png(ln, float(result.settings.get("fit_half_width", 4.0)))
        except Exception:  # noqa: BLE001
            png = b""
        if png:
            parts.append(_img(png, f"The fitted peak of {ln['nuclide']}", (
                f"The {ln['label']} peak at {ln['energy_kev']:.1f} keV: in coincidence with {ln['gate']}, "
                f"Doppler-corrected for the {ln['correction']}, summed over the γ-ray detectors. Area "
                f"{ln['area']:.1f} ± {ln['area_unc']:.1f} counts after {ln['randoms']:.1f} random coincidences "
                "are subtracted.")))
    parts.append("<table><thead><tr><th>Quantity</th>" + "".join(f"<th>{e(ln['label'])}, {ln['energy_kev']:.1f} keV"
                                                                  "</th>" for ln in lines) + "</tr></thead><tbody>")
    for name, value in rows:
        cells = [value(ln) for ln in lines]
        if any(cells):
            parts.append(f"<tr><td>{e(name)}</td>" + "".join(f"<td>{e(c)}</td>" for c in cells) + "</tr>")
    parts.append("</tbody></table>")
    return "".join(parts)


def analysis_html(result, heading: str = "") -> str:
    """One analysis in full for the reports: the result, the fitted spectrum, the chain step by step, the
    uncertainty budget, the shape, and the notes."""
    def e(text) -> str:
        return html.escape(str(text), quote=False)

    parts = [heading]
    if result.b_e2fm4 == result.b_e2fm4:
        parts.append(f"<p>B(E2↑) = {result.b_e2fm4:.4g} e²fm⁴ ± {100 * result.total_unc:.1f}% "
                     f"(statistical {100 * result.statistical:.1f}%, systematic {100 * result.systematic:.1f}%), "
                     f"from {e(result.statistics_from)}; that is {result.b_e2b2:.4g} e²b² or {result.b_wu:.3g} "
                     f"Weisskopf units. The value put in was {result.truth_e2fm4:.4g} e²fm⁴ "
                     f"(pull {result.pull:+.2f}).</p>")
    else:
        parts.append("<p>B(E2↑) could not be determined: the analysis found no peak "
                     f"(the value put in was {result.truth_e2fm4:.4g} e²fm⁴).</p>")
    parts.append(_search_html(result, e))
    try:
        png = analysis_png(result)
    except ImportError:
        png = b""
    except Exception as err:  # noqa: BLE001 - the report stands without the figure
        png = b""
        parts.append(f"<p>The fitted spectrum could not be drawn: {e(err)}.</p>")
    if png:
        fit = result.fit
        parts.append(_img(png, "The fitted peak", (
            f"The γ-ray spectrum the analysis fitted: gated on the particles, Doppler-corrected for the "
            f"{result.spectrum.get('correction', 'emitter')}, summed over the γ-ray detectors. Peak at "
            f"{1e3 * fit.centroid:.1f} keV, FWHM {fit.fwhm_kev:.1f} keV, area {fit.area:.1f} ± {fit.area_unc:.1f} "
            f"counts over the line, of which {fit.randoms:.1f} random coincidences are subtracted.")))
    parts.append(_lines_html(result, e))
    parts.append("<h4>The chain, step by step</h4><ol>"
                 + "".join(f"<li><b>{e(st['step'])}</b>: {e(st['text'])}</li>" for st in result.steps) + "</ol>")
    if result.b_e2fm4 == result.b_e2fm4:
        parts.append("<h4>The uncertainty budget</h4><table><thead><tr><th>Source</th><th>Relative (%)</th></tr>"
                     f"</thead><tbody><tr><td>statistical</td><td>{100 * result.statistical:.2f}</td></tr>"
                     + "".join(f"<tr><td>{e(k)}</td><td>{100 * v:.2f}</td></tr>" for k, v in result.budget.items())
                     + f"<tr><td><b>total</b></td><td><b>{100 * result.total_unc:.2f}</b></td></tr></tbody></table>")
        if result.whole_run_statistical is not None:
            parts.append(f"<p>The whole run, scaled from its real part, would give a statistical uncertainty of "
                         f"{100 * result.whole_run_statistical:.1f}%.</p>")
        sh = result.shape or {}
        if sh:
            rows = [("β₂", f"{sh['beta2']:.3f}"), ("B(E2↑)", f"{result.b_wu:.3g} W.u."),
                    ("Q₀ of a rigid rotor", f"{sh['q0_efm2']:.3g} e fm²"),
                    ("Q_s(2⁺) of that rotor", f"{sh['qs_2plus_efm2']:.3g} e fm²")]
            if "e4_over_e2" in sh:
                rows.append(("E(4⁺)/E(2⁺)", f"{sh['e4_over_e2']:.3f}"))
            parts.append("<h4>The shape of the nucleus</h4><table><thead><tr><th>Quantity</th><th>Value</th></tr>"
                         "</thead><tbody>" + "".join(f"<tr><td>{e(q)}</td><td>{e(v)}</td></tr>" for q, v in rows)
                         + "</tbody></table>" + (f"<p>{e(sh['note'])}</p>" if sh.get("note") else ""))
        parts.append(f"<p>Beam time: {result.counts_per_shift:.3g} counts in the peak per 8 h shift; "
                     f"{result.hours_for_precision:.3g} h for a statistical precision of "
                     f"{100 * float(result.settings.get('wanted_precision', 0.05)):.3g}%.</p>")
    if result.notes:
        parts.append("<ul>" + "".join(f"<li>{e(n)}</li>" for n in result.notes) + "</ul>")
    return "".join(parts)


def alignment_html(a: dict) -> str:
    """The alignment check kept with a run, in a paragraph."""
    o = a["overlay"]
    text = (f"Alignment check: with the target assumed {a['offset_mm']:+g} mm along the beam from its place, the "
            f"corrected peak shifts by {o['shift_kev']:+.2f} keV and its FWHM goes from {o['fwhm_true_kev']:.1f} to "
            f"{o['fwhm_assumed_kev']:.1f} keV.")
    f = a.get("fit")
    if f:
        text += (f" The fit of the centroids against ring puts the target {f['offset_mm']:+.2f} ± "
                 f"{f['uncertainty_mm']:.2f} mm from where the analysis assumed it (truth {f['truth_mm']:+g} mm).")
    return f"<p>{html.escape(text, quote=False)}</p>"


def multistep_html(m: dict) -> str:
    """The multi-step solution (:meth:`~physim.nuclear.planner.Planner.multistep`) as tables: the γ rays in the
    run with all orders and with first order, and the prolate–oblate comparison."""
    def e(text) -> str:
        return html.escape(str(text), quote=False)

    t, dets = m["beam_time_s"], list(m["detectors"])
    parts = [f"<p>{e(m['nuclide'])}: the coupled equations on its level scheme, so excitation in several steps, "
             "the reorientation effect and interference between paths are included. Counts in the planned beam "
             f"time ({_time(t)}), with all orders / with first order.</p>"
             "<table><thead><tr><th>γ ray</th><th>keV</th>" + "".join(f"<th>{e(d)}</th>" for d in dets)
             + "</tr></thead><tbody>"
             + "".join(f"<tr><td>{e(x['label'])}</td><td>{x['energy_kev']:.1f}</td>"
                       + "".join(f"<td>{_g(x['all_orders'][d] * t)} / {_g(x['first_order'][d] * t)}</td>"
                                 for d in dets) + "</tr>" for x in m["transitions"]) + "</tbody></table>"]
    sh = m.get("shapes")
    if sh:
        q = sh["q_efm2"]
        parts.append(f"<p>Prolate or oblate: the 2⁺ state with Q(2⁺) = {q['prolate']:+.1f} (prolate rotor), 0 and "
                     f"{q['oblate']:+.1f} e fm² (oblate). The planned run "
                     + ("can" if sh["separable"] else "cannot") + " tell prolate from oblate: over all rings the "
                     f"two differ by {sh['total_difference_sigma']:.1f} standard deviations of the counts.</p>")
    if m.get("notes"):
        parts.append("<ul>" + "".join(f"<li>{e(n)}</li>" for n in m["notes"][:4]) + "</ul>")
    return "".join(parts)


def doppler_html(planner, analyses: Optional[dict] = None, seen: Optional[dict] = None) -> str:
    """The Doppler correction of the experiment for the report: which nucleus is corrected for and why, the
    formula, where each quantity in it comes from, the geometry it uses, the speeds and shifts it deals with for
    every pair of particle and γ-ray detector, and, for each beam run, the settings its analysis used and the
    width of the corrected peak in each crystal (``seen``: {run number: {"correction", "widths"}} if already
    worked out; the runs are read here otherwise). Empty without Coulomb excitation and γ-ray detectors."""
    from . import dataviews as dv

    exp = planner.experiment
    exc = exp.excitation
    if exc is None or not exp.gamma_detectors:
        return ""

    def e(text) -> str:
        return html.escape(str(text), quote=False)

    emitter = "recoil" if exc.excite == "target" else "projectile"
    who = "the target nucleus, which recoils" if emitter == "recoil" else "the beam nucleus, which is scattered"
    parts = [f"<p>{e(dv.correction_reason(exp, emitter))} The spectra and the analysis correct for {who}.</p>",
             '<p class="f">E₀ = E_measured × (1 − β cos α) / √(1 − β²)</p>',
             "<p>with β the emitter's speed in units of c and α the angle between its direction and the "
             "direction of the γ ray. Neither is measured directly; each event's are worked out from the "
             "segment of the particle detector and the crystal that fired:</p><ol>"
             "<li><b>The particle's direction</b> is the centre of the ring and sector (or the two strips) it "
             "hit, seen from the target at the position the analysis assumes.</li>"
             "<li><b>The orbit.</b> Two-body kinematics at the beam energy in the middle of the target, with the "
             f"excitation energy of {1e3 * exc.energy_mev:g} keV taken from the motion, gives the CM angle of "
             "the collision from that direction, whichever of the two nuclei was detected.</li>"
             "<li><b>The emitter.</b> The same kinematics gives the direction and the energy of the nucleus "
             "that emits: the detected particle itself, or its partner on the other side of the beam.</li>"
             "<li><b>Its speed β</b> is that energy less what the nucleus loses on its way out of the target "
             "from the middle of it.</li>"
             "<li><b>The γ ray's direction</b> is the centre of the crystal that fired. A crystal sees the "
             "emitter over its whole opening angle, so the corrected line keeps a width: the larger β and the "
             "crystal's opening angle, the wider.</li></ol>"]
    tgt = exp.target.position_mm
    parts.append(f"<h3>The geometry it uses</h3><p>The target stands at {tgt:+g} mm along the beam from the "
                 "chamber's centre." if abs(tgt) > 1e-9 else
                 "<h3>The geometry it uses</h3><p>The target stands at the chamber's centre.")
    parts.append(" The table gives each crystal as the correction sees it.</p>"
                 "<table><thead><tr><th>Crystal</th><th>θ (deg)</th><th>φ (deg)</th><th>Distance (mm)</th>"
                 "<th>Radius (mm)</th><th>Opening half-angle (deg)</th></tr></thead><tbody>")
    for i, gd in enumerate(exp.gamma_detectors):
        for label, centre, radius in gd.elements():
            c = np.asarray(centre, dtype=float)
            d = float(np.linalg.norm(c))
            theta = math.degrees(math.acos(max(-1.0, min(1.0, c[2] / d))))
            phi = math.degrees(math.atan2(c[1], c[0]))
            name = (gd.name or f"G{i + 1}") + (f" {label}" if label else "")
            parts.append(f"<tr><td>{e(name)}</td><td>{theta:.1f}</td><td>{phi:.1f}</td><td>{d:.1f}</td>"
                         f"<td>{radius:.1f}</td><td>{math.degrees(math.atan2(radius, d)):.1f}</td></tr>")
    parts.append("</tbody></table>")
    try:
        rows = planner.gamma().get("doppler", [])
    except Exception:  # noqa: BLE001 - the report stands without the table
        rows = []
    if rows:
        parts.append("<h3>The speeds and the shifts</h3><p>For the scattered particle reaching each particle "
                     "detector and the γ ray each γ-ray detector: the emitter's mean speed, the mean energy the "
                     "crystal measures, its shift from the transition energy, the width the motion gives the "
                     "measured line, and its width with the detector's resolution.</p><table><thead><tr>"
                     "<th>Particle detector</th><th>γ-ray detector</th><th>β</th><th>Measured (keV)</th>"
                     "<th>Shift (keV)</th><th>Doppler FWHM (keV)</th><th>FWHM (keV)</th></tr></thead><tbody>"
                     + "".join(f"<tr><td>{e(r['particle_detector'])}</td><td>{e(r['gamma_detector'])}</td>"
                               f"<td>{r['beta_mean']:.4f}</td><td>{r['mean_kev']:.1f}</td>"
                               f"<td>{r['shift_kev']:+.1f}</td><td>{r['doppler_fwhm_kev']:.1f}</td>"
                               f"<td>{r['fwhm_kev']:.1f}</td></tr>" for r in rows) + "</tbody></table>")
    beam_runs = [s for s in planner.runs() if s["kind"] != "source"]
    if beam_runs:
        from .alignment import peak_widths
        from .response import FWHM_PER_SIGMA

        parts.append("<h3>In each run</h3><table><thead><tr><th>Run</th><th>Corrected for</th>"
                     "<th>Target assumed</th><th>Particles</th><th>Corrected line, FWHM per crystal (keV)</th>"
                     "</tr></thead><tbody>")
        for s in beam_runs:
            r = (analyses or {}).get(s["number"])
            settings = dict(getattr(r, "settings", {}) or {})
            offset = float(settings.get("target_offset_mm", 0.0) or 0.0)
            try:
                if seen is not None and s["number"] in seen:
                    widths, correction = seen[s["number"]]["widths"], seen[s["number"]]["correction"]
                else:
                    run = planner.run if planner.run is not None and planner.run.number == s["number"] \
                        else planner.folder.load_run(s["number"])
                    g = run.gammas()
                    widths = peak_widths(g, run.experiment) if g is not None else {}
                    correction = _correction_of(run, r)
                lines = ", ".join(f"{k} {1e3 * FWHM_PER_SIGMA * v:.1f}" for k, v in widths.items()) or "no γ rays"
            except Exception as err:  # noqa: BLE001 - the report stands without the row's widths
                lines, correction = f"not read ({err})", emitter
            gate = _gate_of(planner, settings) if r is not None else None
            parts.append(f"<tr><td>{s['number']}</td><td>the {e(correction)}"
                         + ("" if r is not None else " (no analysis yet: the default)") + "</td><td>"
                         + ("where it is" if abs(offset) < 1e-9 else f"{offset:+g} mm from where it is")
                         + f"</td><td>{e(_gate_words(gate))}</td><td>{e(lines)}</td></tr>")
        parts.append("</tbody></table><p>The width of the corrected line is measured on the run's spectrum of "
                     "each crystal (its width at half maximum); it is never taken narrower than three times the "
                     "detector's resolution.</p>")
    return "".join(parts)


def _not_reached(experiment) -> list:
    """The levels of the excited nucleus's scheme that the simulated events never populate
    (:func:`~physim.nuclear.gamma.unreached_levels`): with first-order events those no matrix element joins to
    the ground state and no populated level decays to. [(label, energy in keV)]."""
    try:
        from .gamma import cascade_excitation, unreached_levels

        missing = unreached_levels(experiment)
        if not missing:
            return []
        levels = cascade_excitation(experiment).scheme.levels
        return [(levels[n].label, float(levels[n].energy.value)) for n in missing]
    except Exception:  # noqa: BLE001 - a note the report can do without
        return []


def pdf_browser() -> Optional[str]:
    """The program that prints a page to PDF: a Chromium-based browser (Edge, Chrome or Chromium) found on this
    computer, or the one ``$PHYSIM_PDF_BROWSER`` names. None if there is none."""
    import os
    import shutil

    named = os.environ.get("PHYSIM_PDF_BROWSER")
    if named:
        return named if Path(named).is_file() or shutil.which(named) else None
    for name in ("msedge", "chrome", "google-chrome", "chromium", "chromium-browser", "microsoft-edge"):
        found = shutil.which(name)
        if found:
            return found
    roots = [os.environ.get(k) for k in ("PROGRAMFILES(X86)", "PROGRAMFILES", "LOCALAPPDATA")]
    tails = (r"Microsoft\Edge\Application\msedge.exe", r"Google\Chrome\Application\chrome.exe",
             r"Chromium\Application\chrome.exe")
    for root in roots:
        for tail in tails:
            if root and (Path(root) / tail).is_file():
                return str(Path(root) / tail)
    for path in ("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
                 "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
                 "/Applications/Chromium.app/Contents/MacOS/Chromium"):
        if Path(path).is_file():
            return path
    return None


def html_to_pdf(page: str, timeout: float = 240.0) -> Optional[bytes]:
    """A self-contained HTML page as a PDF (bytes), printed by a browser of this computer without a window
    (:func:`pdf_browser`), with the page's own print styles. None when there is no such browser or it fails: the
    page itself prints to PDF from any browser."""
    import subprocess
    import tempfile

    browser = pdf_browser()
    if browser is None:
        return None
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
        src, out = Path(tmp) / "report.html", Path(tmp) / "report.pdf"
        src.write_text(page, encoding="utf-8")
        cmd = [browser, "--headless", "--disable-gpu", "--no-first-run", "--disable-extensions",
               f"--user-data-dir={Path(tmp) / 'profile'}", "--no-pdf-header-footer",
               f"--print-to-pdf={out}", src.as_uri()]
        import time

        started = time.monotonic()
        try:
            subprocess.run(cmd, check=False, timeout=timeout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.TimeoutExpired):
            return None
        # The browser's launcher may return before the page is printed: wait for the file, and until it has
        # stopped growing. A minute without a file means the browser gave none.
        size, steady = -1, 0
        while steady < 3 and time.monotonic() - started < timeout:
            now = out.stat().st_size if out.is_file() else -1
            steady = steady + 1 if now == size and now > 0 else 0
            if now < 0 and time.monotonic() - started > 60:
                break
            size = now
            time.sleep(0.3)
        try:
            return out.read_bytes() if size > 0 else None
        except OSError:
            return None


def run_analyses(planner, result=None) -> dict:
    """{run number: Result} for every beam run of the experiment that has been analysed: ``result`` for the
    current run if given, else what :meth:`~physim.nuclear.planner.Planner.run_analysis` finds (this session's,
    or the one kept in the run's folder)."""
    out = {}
    for s in planner.runs():
        if s["kind"] == "source":
            continue
        current = planner.run is not None and planner.run.number == s["number"]
        if current and result is not None:
            out[s["number"]] = result
            continue
        try:
            r = planner.run_analysis(planner.run if current else planner.folder.load_run(s["number"]))
        except Exception:  # noqa: BLE001 - a run that cannot be read has no analysis to show
            r = None
        if r is not None:
            out[s["number"]] = r
    return out


def _plan_report(planner, events: int, seed: int, journal: str, width: str, validation: bool):
    """The beam-time report of the current setup (:class:`physim.nuclear.report.Report`), or None with the
    reason it could not be built."""
    try:
        from .report import build

        return build(planner.experiment, seed=seed, events=events, validation=validation, journal=journal,
                     width=width), ""
    except Exception as err:  # noqa: BLE001 - the report stands without the plan's figures
        return None, str(err)


#: The sections of the report, in order: (anchor, title).
SECTIONS = (("summary", "Summary"), ("setup", "The setup"), ("plan", "The plan: what this setup should give"),
            ("runs", "The runs and their spectra"), ("gates", "The gates"),
            ("doppler", "The Doppler correction"), ("analysis", "The analysis"),
            ("numbers", "Every number, explained"), ("method", "The data, the method and what is left out"))


def experiment_html(planner, result=None, scene_png: Optional[str] = None, plan=None, plan_events: int = 50_000,
                    seed: int = 1, journal: str = "physical_review", width: str = "single",
                    validation: bool = False) -> str:
    """The report of the experiment, one self-contained page with everything in it: a summary; the setup with
    its warnings, level schemes and how it changed between runs; the plan (what the beam-time report holds:
    detectors, rates and beam time, kinematics, expected peaks, Coulomb excitation, energy loss, predicted
    spectra); every run with its summary and its spectra (the particle detectors, the energy against ring with
    the gates drawn on it, every γ-ray detector as measured and Doppler-corrected to the top of the level
    scheme); the gates and what each does to the γ-ray spectrum; the Doppler correction and its parameters; the
    analysis of every analysed run, the alignment and the multi-step solution; every number explained (what the
    run record holds); and the data, the method and what the simulation leaves out.

    ``result`` is the analysis of the current run (a :class:`~physim.nuclear.analysis.Result`; the kept one is
    used without it). ``plan`` is the setup's :class:`~physim.nuclear.report.Report` if already built; otherwise
    it is built here with ``plan_events`` Monte Carlo events, ``seed`` and the figures' ``journal`` style and
    ``width``, and ``validation`` runs the register's checks for the models' table."""
    from .planner import REGISTER
    from .record import LEFT_OUT, analysis_explanations, explanations, method_html, provenance

    def e(text) -> str:  # text, not attributes: quotes stay as they are
        return html.escape(str(text), quote=False)

    exp = planner.experiment
    why_no_plan = ""
    if plan is None:
        plan, why_no_plan = _plan_report(planner, plan_events, seed, journal, width, validation)
    pp = None
    if plan is not None:
        try:
            pp = plan.parts(dpi=150)
        except Exception as err:  # noqa: BLE001 - the report stands without the plan's figures
            why_no_plan = str(err)
    runs = planner.runs()
    analyses = run_analyses(planner, result)
    gates = planner.gates()
    pieces = {s["number"]: _run_pieces(planner, s, analyses.get(s["number"]), bool(gates)) for s in runs}
    b, t = exp.beam, exp.target
    sec = {key: [] for key, _ in SECTIONS}

    # -- 1. summary -------------------------------------------------------------------------------------------
    what = (f"Coulomb excitation of the {exp.excitation.excite} to {1e3 * exp.excitation.energy_mev:g} keV "
            f"({exp.excitation.multipolarity})" if exp.excitation else "elastic scattering")
    beam_runs = [s for s in runs if s["kind"] != "source"]
    items = [f"{e(b.nuclide)} at {e(str(b.energy))} on {e(t.material)} ({e(str(t.thickness))}): {e(what)}.",
             (f"{len(runs)} run{'s' * (len(runs) != 1)}: {len(beam_runs)} with beam, {len(runs) - len(beam_runs)} "
              "with a source." if runs else "No run yet: the report holds the setup and the plan.")]
    for n, r in analyses.items():
        if r.b_e2fm4 == r.b_e2fm4:
            items.append(f"Run {n}: B(E2↑) = {r.b_e2fm4:.4g} e²fm⁴ ± {100 * r.total_unc:.1f}% "
                         f"({r.b_wu:.3g} W.u.; the value put in was {r.truth_e2fm4:.4g} e²fm⁴, pull "
                         f"{r.pull:+.2f}).")
        else:
            items.append(f"Run {n}: the analysis found no peak.")
    if beam_runs and not analyses:
        items.append("No run has been analysed yet.")
    sec["summary"].append("<ul>" + "".join(f"<li>{x}</li>" for x in items) + "</ul>")

    # -- 2. the setup -----------------------------------------------------------------------------------------
    if exp.description:
        sec["setup"].append(f"<p>{e(exp.description)}</p>")
    if planner.folder is not None:
        sec["setup"].append(f"<p>Folder: <code>{e(str(planner.experiment_folder))}</code></p>")
    sec["setup"].append("<ul>"
                        f"<li>Beam: {e(b.nuclide)} at {e(str(b.energy))}, {e(str(b.current))}"
                        + (f", charge state {b.charge_state}+" if b.charge_state else "") + ".</li>"
                        f"<li>Target: {e(t.material)}, {e(str(t.thickness))}.</li>"
                        f"<li>Particle detectors: {e(', '.join(d.name or f'D{i + 1}' for i, d in enumerate(exp.detectors)))}.</li>"
                        f"<li>γ-ray detectors: {e(', '.join(g.name or f'G{i + 1}' for i, g in enumerate(exp.gamma_detectors))) or 'none'}.</li>"
                        f"<li>Reaction: {e(what)}.</li>"
                        f"<li>Beam time planned: {e(str(exp.run.beam_time))}"
                        + (f", for {exp.run.counts_wanted} counts per detector" if exp.run.counts_wanted else "")
                        + ".</li></ul>")
    if pp:
        sec["setup"].append(pp["setup_figures"])
        sec["setup"].append("<h3>Warnings</h3>" + pp["warnings"])
    elif scene_png:
        sec["setup"].append(f'<p><img src="data:image/png;base64,{scene_png}" alt="the scene" '
                            'style="max-width:100%"></p>')
    schemes = level_schemes_html(planner)
    if schemes:
        sec["setup"].append("<h3>The level schemes</h3>" + schemes)
    history = setup_history(planner)
    if history:
        sec["setup"].append("<h3>How the setup changed between runs</h3><ul>")
        for h in history:
            which = f"run {h['run']} ({e(h['label'])})" if h["run"] is not None else e(h["label"])
            if h["run"] is not None and h is history[0]:
                sec["setup"].append(f"<li>{which}: the setup as it started.</li>")
            else:
                sec["setup"].append(f"<li>{which}: " + (e("; ".join(h["changes"])) if h["changes"] else
                                                        "no change to the physics") + "</li>")
        sec["setup"].append("</ul>")

    # -- 3. the plan ------------------------------------------------------------------------------------------
    sec["plan"].append("<p>What the current setup should give, computed from the setup and not from the runs: "
                       "the runs are compared with it in the next section.</p>")
    try:
        pl = planner.plan()
    except Exception:  # noqa: BLE001 - the report stands without the plan
        pl = None
    if pl:
        sec["plan"].append("<h3>Beam time and load</h3>"
                           f"<p>Beam time needed: {_time(pl['beam_time_needed_s'])}"
                           + (f" ({e(str(pl['limiting_detector']))}, for {pl['counts_wanted']} {e(pl['measured'])})"
                              if pl["limiting_detector"] else "") + f"; planned: {_time(pl['beam_time_s'])}.</p>")
        sec["plan"].append("<table><thead><tr><th>Detector</th><th>Rate (1/s)</th><th>Busy</th>"
                           "<th>Coincidences (1/s)</th><th>Safe rings</th></tr></thead><tbody>"
                           + "".join(f"<tr><td>{e(d['detector'])}</td><td>{_g(d['rate_per_s'])}</td>"
                                     f"<td>{100 * d['dead_time_fraction']:.2g}%</td><td>{_g(d['coincidence_per_s'])}</td>"
                                     f"<td>{len(d['safe_rings'])} of {d['rings']}</td></tr>" for d in pl["detectors"])
                           + "</tbody></table>")
        if pl.get("other_peak"):
            op = pl["other_peak"]
            sec["plan"].append(
                f"<p>{e(op['nuclide'])} is excited too (its level scheme is in the setup). Its {e(op['label'])} peak "
                f"at {op['energy_kev']:.1f} keV, in coincidence: "
                + "; ".join(f"{e(d['detector'])} {_g(d['other_coincidence_per_s'])} /s" for d in pl["detectors"]
                            if d.get("other_coincidence_per_s") is not None)
                + ". The counts wanted and the beam time are for the nucleus of the reaction.</p>")
    if pp:
        sec["plan"] += ["<h3>Detectors: coverage, cross sections, rates and beam time</h3>", pp["detectors"],
                        "<h3>Kinematics</h3>", pp["kinematics"], "<h3>Expected peaks</h3>", pp["peaks"],
                        pp["excitation_plain"], "<h3>Energy loss</h3>", pp["energy_loss"],
                        "<h3>Simulated spectra</h3>", pp["spectra"]]
    elif why_no_plan:
        sec["plan"].append(f"<p>The plan's tables and figures could not be built: {e(why_no_plan)}.</p>")

    # -- 4. the runs ------------------------------------------------------------------------------------------
    if not runs:
        sec["runs"].append("<p>No run yet.</p>")
    elif exp.excitation is not None and exp.gamma_detectors:
        from .gamma import excitation_order

        missing = _not_reached(exp)
        first_order = excitation_order(exp) == "first"
        sec["runs"].append(
            "<p>Each beam run is drawn three times: the energy in each particle detector, the energy against "
            "ring with the gates on it, and the γ-ray spectra as measured and Doppler-corrected.</p>"
            + ("<p>The simulated events of this setup "
               + ("are first order (<code>order = \"first\"</code>): each level is excited in one step from the "
                  "ground state, and a level reached only in two steps is not populated"
                  if first_order else "do not populate every level of the scheme: no chain of matrix elements "
                                      "joins these to the ground state")
               + ". Their γ rays are not in these spectra: "
               + e(", ".join(f"{lab} at {en:.0f} keV" for lab, en in missing)) + ".</p>"
               if missing else ""))
    for s in runs:
        sec["runs"].append(f"<h3>Run {s['number']}: {e(s['describe'])}</h3><ul>"
                           f"<li>Folder <code>{e(s['label'])}</code>; seed {s.get('seed')}; finished "
                           f"{e(str(s.get('finished', ''))[:19])} (UTC).</li>")
        if s["kind"] == "source":
            src = s.get("source", {})
            sec["runs"].append(f"<li>{e(str(src.get('nuclide')))} source of {e(str(src.get('activity')))} at "
                               f"{e(str(src.get('position') or 'the target'))}; counts per crystal: "
                               + e(", ".join(f"{k} {v:,.0f}" for k, v in s.get("gamma_counts", {}).items()))
                               + ".</li>")
        else:
            sec["runs"].append(f"<li>{_time(s['real_s'])} simulated event by event of {_time(s['duration_s'])}"
                               + ("; the rest scaled" if s.get("scaled") else "")
                               + f"; {s.get('particles', 0):,} particles and {s.get('gammas', 0):,} γ rays "
                               "stored.</li><li>Counts (rate): "
                               + e(", ".join(f"{d} {c:,} ({s['rates_per_s'].get(d, 0):.4g}/s)"
                                             for d, c in s.get("counts", {}).items())) + ".</li>")
            if s.get("live_fraction") is not None:
                sec["runs"].append(f"<li>Live fraction {s['live_fraction']:.4f}.</li>")
            if s.get("assumed_offset_mm") is not None:
                sec["runs"].append(f"<li>Alignment run: the target assumed {s['assumed_offset_mm']:+g} mm off.</li>")
        sec["runs"].append("<li>" + ("Taken with another setup than the current one." if s["stale"]
                                     else "Taken with the current setup.") + "</li></ul>")
        sec["runs"].append(pieces[s["number"]]["figures"])
    if planner.run is not None and planner.run.kind != "source":
        c = planner.compare()
        if c["rows"]:
            sec["runs"].append(f"<h3>Run {c['run']}: predicted against measured</h3><table><thead><tr>"
                               "<th>Quantity</th><th></th><th>Predicted</th><th>Measured</th><th></th></tr></thead>"
                               "<tbody>"
                               + "".join(f"<tr><td>{e(r['quantity'])}</td><td>{e(r['what'])}</td>"
                                         f"<td>{_g(r['predicted'])}</td><td>{_g(r['measured'])} ± "
                                         f"{_g(r['uncertainty'])}</td><td>"
                                         f"{ {True: 'agree', False: 'differ', None: ''}[r['agree']] }</td></tr>"
                                         for r in c["rows"]) + "</tbody></table>")

    # -- 5. the gates -----------------------------------------------------------------------------------------
    if gates:
        sec["gates"].append("<p>A gate is a condition on the particles; a γ-ray spectrum under a gate holds the "
                            "γ rays in coincidence with the particles that pass it.</p><ul>"
                            + "".join(f"<li><b>{e(g.name)}</b>: {e(g.describe())}"
                                      + (f" ({e(g.note)})" if g.note else "") + ".</li>" for g in gates) + "</ul>")
        sec["gates"] += [pieces[s["number"]]["gates"] for s in beam_runs]
    else:
        sec["gates"].append("<p>No gates: the spectra take every counted particle. A gate is made on the Data "
                            "tab, from the energy-against-ring view or by naming detectors and rings.</p>")

    # -- 6. the Doppler correction ----------------------------------------------------------------------------
    try:
        dop = doppler_html(planner, analyses, {n: x["seen"] for n, x in pieces.items() if x["seen"] is not None})
    except Exception as err:  # noqa: BLE001 - the report stands without the section
        dop = f"<p>The Doppler correction could not be described: {e(err)}.</p>"
    sec["doppler"].append(dop or "<p>No excited state or no γ-ray detector: there is no γ ray to correct.</p>")

    # -- 7. the analysis --------------------------------------------------------------------------------------
    for s in runs:
        r = analyses.get(s["number"])
        if r is not None:
            sec["analysis"].append(analysis_html(r, f"<h3>Run {s['number']}: from the peak area to B(E2)</h3>"))
        if s["kind"] != "source":
            try:
                al = planner.kept_alignment(planner.run if planner.run is not None
                                            and planner.run.number == s["number"]
                                            else planner.folder.load_run(s["number"]))
            except Exception:  # noqa: BLE001 - the report stands without it
                al = None
            if al:
                sec["analysis"].append((f"<h3>Run {s['number']}: the alignment</h3>" if r is None else "")
                                       + alignment_html(al))
    if not analyses and result is not None:  # the analysis of a simulated sample, without a run
        sec["analysis"].append(analysis_html(result, "<h3>The simulated sample</h3>"))
    elif not analyses:
        sec["analysis"].append("<p>No run has been analysed yet: take a beam run and press Analyse on the "
                               "Analysis step. The analysis is then kept with the run and appears here.</p>")
    m = planner.kept_multistep()
    if m and m.get("available"):
        sec["analysis"].append("<h3>All orders: multi-step excitation and reorientation</h3>" + multistep_html(m))

    # -- 8. every number, explained ---------------------------------------------------------------------------
    sec["numbers"].append("<p>Each number with its formula, the formula with this experiment's numbers in it, "
                          "what it means and what it assumes.</p><h3>The plan's numbers</h3>")
    try:
        sec["numbers"] += [_explanation_html(x) for x in explanations(planner)]
    except Exception as err:  # noqa: BLE001 - the report stands without them
        sec["numbers"].append(f"<p>Not available: {e(err)}.</p>")
    xs = run_explanations(planner)
    if xs:
        sec["numbers"].append(f"<h3>Run {planner.run.number}: its counters, its gates, and the plan against what "
                              "it measured</h3>")
        sec["numbers"] += [_explanation_html(x) for x in xs]
    current = analyses.get(planner.run.number) if planner.run is not None else result
    if current is not None:
        try:
            ax = analysis_explanations(planner, current)
        except Exception:  # noqa: BLE001 - the report stands without them
            ax = []
        if ax:
            sec["numbers"].append("<h3>From the peak area to B(E2) and the shape</h3>")
            sec["numbers"] += [_explanation_html(x) for x in ax]

    # -- 9. the data, the method and what is left out ---------------------------------------------------------
    sec["method"].append("<h3>The nuclear data and where it comes from</h3>"
                         "<table><thead><tr><th>Quantity</th><th>Value</th><th>Source</th></tr></thead><tbody>"
                         + "".join(f"<tr><td>{e(q)}</td><td>{e(str(v))}</td><td>{e(src)}</td></tr>"
                                   for q, v, src in provenance(exp)) + "</tbody></table>"
                         "<p>Sources: <em>ensdf</em> is the evaluated data as read; <em>derived</em> is computed "
                         "from ENSDF's values by physim; <em>assumed</em> is a typical value, not a measurement; "
                         "<em>user</em> was given in the setup.</p>")
    if exp.excitation is not None:
        sec["method"].append("<h3>The method, step by step</h3>" + method_html())
    sec["method"].append("<h3>What the simulation leaves out</h3><ul>"
                         + "".join(f"<li>{e(x)}</li>" for x in LEFT_OUT) + "</ul>")
    if pp:
        sec["method"] += ["<h3>Models and their validation</h3>", pp["models"],
                          "<h3>Data and references</h3>", pp["references"]]
    sec["method"].append("<h3>The physics behind it</h3><ul>"
                         + "".join(f'<li><a href="{DOCS}{e(path)}.html">{e(tab.replace("_", " "))}</a></li>'
                                   for tab, path in REGISTER.items()) + "</ul>")

    body = [f"<h1>{e(planner.name)}: the report</h1>"]
    if pp:
        body.append(f'<p class="sub">{pp["meta"]}</p>')
    body.append('<ol class="toc">' + "".join(f'<li><a href="#{key}">{e(title)}</a></li>' for key, title in SECTIONS)
                + "</ol>")
    for i, (key, title) in enumerate(SECTIONS, 1):
        body.append(f'<h2 id="{key}">{i}. {e(title)}</h2>')
        body += sec[key]
    return ("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>"
            f"{e(planner.name)}: the report</title><style>{record_css()}</style></head><body>"
            + "\n".join(body) + "</body></html>")


def experiment_pdf(planner, result=None, page: Optional[str] = None, **options) -> Optional[bytes]:
    """The report of the experiment as a PDF: :func:`experiment_html` (or ``page``, if already made) printed by
    a browser of this computer (:func:`html_to_pdf`). None when no such browser is found."""
    return html_to_pdf(page if page is not None else experiment_html(planner, result, **options))


def experiment_zip(planner, result=None, folder: bool = True, pdf: bool = True, tables: bool = True,
                   figures: bool = False, plan_events: int = 50_000, seed: int = 1,
                   journal: str = "physical_review", width: str = "single", validation: bool = False) -> bytes:
    """The experiment as a zip: ``report.html`` (the report, :func:`experiment_html`), ``report.pdf`` (the same,
    when a browser of this computer can print it; ``pdf``), ``setup.toml``, the plan's tables as CSV under
    ``tables/`` (``tables``) and its figures as PNG and PDF under ``figures/`` (``figures``), and the experiment
    folder with every run under ``experiment/`` (``folder``)."""
    import tempfile

    plan, _ = _plan_report(planner, plan_events, seed, journal, width, validation)
    page = experiment_html(planner, result, plan=plan, plan_events=plan_events, seed=seed, journal=journal,
                           width=width, validation=validation)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("report.html", page)
        if pdf:
            data = html_to_pdf(page)
            if data:
                z.writestr("report.pdf", data)
        z.writestr("setup.toml", planner.to_toml())
        if plan is not None and (tables or figures):
            with tempfile.TemporaryDirectory() as tmp:
                if tables:
                    for p in plan.write_csv(Path(tmp) / "tables"):
                        z.write(p, "tables/" + p.name)
                if figures:
                    for p in plan.write_figures(Path(tmp) / "figures"):
                        z.write(p, "figures/" + p.name)
        if folder and planner.folder is not None:
            root = planner.experiment_folder
            for p in sorted(root.rglob("*")):
                if p.is_file():
                    z.write(p, "experiment/" + p.relative_to(root).as_posix())
    return buf.getvalue()


__all__ = ["GAMMA_BIN_KEV", "SECTIONS", "alignment_html", "analysis_html", "analysis_png", "line_png", "search_png", "doppler_html",
           "experiment_html", "experiment_pdf", "experiment_zip", "gamma_lines", "html_to_pdf",
           "level_scheme_png", "level_schemes", "level_schemes_html", "multistep_html", "pdf_browser",
           "run_analyses", "run_explanations", "run_gamma_png", "run_gates_png", "run_rings_png",
           "run_spectra_png", "run_views", "setup_history"]
