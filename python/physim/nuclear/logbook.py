"""The report of an experiment (backlog item 69): one page that leads with the results and keeps everything else
behind them.

The report opens with what the experiment found (:data:`RESULTS`): the result with its uncertainties, its method
and its runs; the result beside the other routes to the same quantity (the comparison of
:mod:`physim.nuclear.compare` over the rows of the results ledger, :mod:`physim.nuclear.results`); the
comparisons still to run; what limits the result; and what it means. A simulated run analysed back is *projected*:
it is closed against the value its simulation was given and set beside how well the quantity is already known,
never pulled against a measured or evaluated value as if it tested it.

Behind them come the appendices (:data:`APPENDICES`): the setup, the plan, every run and its spectra, the gates,
the Doppler correction, the analysis, every number explained, and the method. They are what used to be three
pages (the beam-time report of the setup, the record of the experiment and the run record), kept whole::

    from physim.nuclear import logbook

    page = logbook.experiment_html(planner, result)          # the report: one self-contained HTML page
    pdf = logbook.experiment_pdf(planner, result)            # the same as a PDF, printed by a browser found here
    data = logbook.experiment_zip(planner, result)           # report.html, report.pdf, the setup, the tables,
                                                             # and the experiment folder with every run
    for x in logbook.run_explanations(planner, planner.run):  # each with its formula, numbers and meaning
        assert x.recompute() == x.value
"""

from __future__ import annotations

import ast
import base64
import html
import io
import math
import re
import zipfile
from pathlib import Path
from typing import Optional

import numpy as np

from . import compare, nucleus, observables, results
from .record import DOCS, Explanation, record_css


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
            f"R = {n} / {real:.6g} s = {number(n / real)} /s",
            f"How many particles {d} counted per second in the part of run {run.number} simulated event by event, "
            "where each event counts once.", "Unweighted events: each one stands for one reaction in the beam time; "
            "particles below the threshold are not recorded.", DOCS + "physics-register/rates-and-events.html",
            {"n": n, "t": real}, lambda n, t: n / t))
        out.append(Explanation(
            "run_counts", f"Counts in {d} over the whole run", n * duration / real, "", "N_run = N_real × t_run / t_real",
            f"N_run = {n} × {duration:.6g} s / {real:.6g} s = {number(n * duration / real)}",
            f"The real part scaled to the run's {_time(duration)}; its histograms are redrawn with the full run's Poisson "
            "fluctuations.", "The beam and the setup stay the same over the whole run.",
            DOCS + "physics-register/rates-and-events.html", {"n": n, "t_run": duration, "t": real},
            lambda n, t_run, t: n * t_run / t))
        dead = _seconds(exp.run.dead_time)
        rate = n / real
        out.append(Explanation(
            "run_busy", f"Busy fraction of {d}", rate * dead / (1 + rate * dead), "", "b = R τ / (1 + R τ)",
            f"b = {number(rate)} /s × {dead:.3g} s / (1 + {number(rate)} /s × {dead:.3g} s) = "
            f"{number(rate * dead / (1 + rate * dead))}",
            "The share of the time this detector's electronics are busy and miss a particle.",
            "Non-paralysable dead time τ per count, the same for every count.", DOCS + "theory/coulex.html",
            {"r": rate, "tau": dead}, lambda r, tau: r * tau / (1 + r * tau)))
        peak = s.get("peak_coincidences", {})
        k = sum(v for key, v in peak.items() if key.split("|")[0] == d)
        if s.get("gammas"):
            out.append(Explanation(
                "run_coincidences", f"Coincidences in the full-energy peak with {d}", k / real, "1/s",
                "R_γ = N_peak / t_real", f"R_γ = {k} / {real:.6g} s = {number(k / real)} /s",
                f"γ rays that left their full energy in a crystal while {d} counted the particle of the same event.",
                "Add-back and suppression as the setup has them; random coincidences are not in this count.",
                DOCS + "theory/coulex.html", {"k": k, "t": real}, lambda k, t: k / t))
    live = s.get("live_fraction")
    if live is not None:
        total = sum(s["counts"].values()) / real
        out.append(Explanation(
            "run_live", "Live fraction of the run", live, "", "L = 1 / (1 + τ Σ R)",
            f"L = 1 / (1 + {_seconds(exp.run.dead_time):.3g} s × (Σ R)) = {number(live)}",
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
            f"f = {passed} / {total} = {number(passed / total if total else 0.0)}",
            f"The share of the run's counted particles the gate keeps: {g.describe()}.",
            "A group is chosen by energy against the kinematic lines, as an experimentalist would, not by the "
            "simulation's knowledge of the reaction.", DOCS + "planner-app.html",
            {"n": passed, "total": total}, lambda n, total: n / total if total else 0.0))
    for r in planner.compare(run)["rows"]:
        if r["measured"] is None or not r["uncertainty"] or r["predicted"] is None:
            continue
        name = f"compare_{r['quantity']}_{r['what']}".replace(" ", "_")
        if r["agree"] is None:  # shown side by side and not compared (a line under a trigger): no pull
            out.append(Explanation(
                name, f"{r['quantity']}, {r['what']}: measured beside predicted", r["measured"], r["unit"],
                "the run's number beside the Plan's: the two are not compared",
                f"measured {plus_minus(r['measured'], r['uncertainty'])}; predicted {number(r['predicted'])}",
                "What the run counted, beside what the Plan gives for the same line. The two are not the same "
                "quantity (the table of predicted against measured says what differs), so no difference in "
                "standard deviations is given.", "The prediction is the Plan's at the moment the run was started.",
                DOCS + "physics-register/rates-and-events.html", {"m": r["measured"]}, lambda m: m))
            continue
        pull = (r["measured"] - r["predicted"]) / r["uncertainty"]
        out.append(Explanation(
            name, f"{r['quantity']}, {r['what']}: measured against predicted", pull, "σ",
            "z = (measured − predicted) / σ_measured",
            f"z = ({number(r['measured'])} − {number(r['predicted'])}) / {number(r['uncertainty'])} = {pull:+.2f}",
            "How far the run is from the Plan's prediction, in standard deviations of the run's counting statistics; "
            "within ±3 they agree. The Plan's own accuracy is not in σ, so with many counts a difference of a percent "
            "already is many σ.", "The prediction is the Plan's at the moment the run was started.",
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
        from .runs import taken_setup

        setup = taken_setup(s["_folder"], s.get("kind")).to_dict()
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
    if run.kind == "fusion":  # no particle detectors, and γ rays kept per exit channel: the figures below do not fit
        out["figures"] = ("<p>The spectra of a fusion run are not drawn in the report yet: it has no particle "
                          "detectors, and its γ rays are kept per exit channel. The run against the plan is in "
                          "the table at the end of this appendix.</p>")
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
                       f"{number(exc.b_up_e2fm)} e²fm{str(2 * int(lam[1])).translate(_SUPERSCRIPT)}). No level "
                       "scheme has been looked up for this "
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
    line within the fit window (with its tails, where it has them), the ±3σ band, and the random coincidences expected.
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
        ax.plot(1e3 * x, k * fit.curve(x), color="#D55E00", lw=1.3, label=f"fit: {fit.shape}")
        ax.axvspan(1e3 * (fit.centroid - 3 * fit.sigma), 1e3 * (fit.centroid + 3 * fit.sigma), color="#D55E00",
                   alpha=0.10, lw=0, label="±3σ: the area" if fit.tails is None and fit.below is None else "±3σ")
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
    the gated, corrected counts, the fitted peak on its background, the ±3σ band and the random
    coincidences expected. Needs matplotlib."""
    from matplotlib.figure import Figure

    from .analysis import PeakFit
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
        peak = PeakFit.of(fit)
        ax.plot(1e3 * x, k * peak.curve(x), color="#D55E00", lw=1.3, label=f"fit: {peak.shape}")
        ax.axvspan(1e3 * (centre - 3 * sigma), 1e3 * (centre + 3 * sigma), color="#D55E00", alpha=0.10, lw=0,
                   label="±3σ: the area" if peak.tails is None and peak.below is None else "±3σ")
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
            ("Counts over the background", lambda ln: f"{ln['fit']['area']:,.1f} ± {ln['fit']['area_unc']:,.1f}"),
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
             "on a straight line is fitted, and the area is the counts within ±3σ of its centre over that line. "
             "Where the peak has tails it is the area of the peak fitted with them (the tails of a simulation of the "
             "setup, if the counts are too few to fix them). Where the fit window reaches the line's own Compton edge, "
             "the background is that line and the continuum below the peak that the detectors' response gives. The "
             "random coincidences expected within ±3σ are subtracted. The area over the efficiency and the angular "
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


def analysis_html(result, heading: str = "", key=None) -> str:
    """One analysis in full for the reports: the result, the fitted spectrum, the chain step by step, the
    uncertainty budget, the shape, and the notes. ``key`` is the transition's key (:func:`run_key`): Weisskopf
    units are given for the decay, B(E2)↓, and only where the key holds the spins of both levels
    (:func:`decay_of`); without them the excitation's value stands in e²fm⁴ alone."""
    def e(text) -> str:
        return html.escape(str(text), quote=False)

    parts = [heading]
    if result.b_e2fm4 == result.b_e2fm4:
        decay = decay_of(result.b_e2fm4, key)
        parts.append(f"<p>B(E2↑) = {number(result.b_e2fm4)} e²fm⁴ ± {100 * result.total_unc:.1f}% "
                     f"(statistical {100 * result.statistical:.1f}%, systematic {100 * result.systematic:.1f}%), "
                     f"from {e(result.statistics_from)}; that is {number(result.b_e2b2)} e²b² ↑. "
                     + (f"For the decay, B(E2↓) = B(E2↑) / {number(decay[2])} = {number(decay[0])} e²fm⁴ = "
                        f"{number(decay[1], 3)} Weisskopf units (↓). " if decay is not None else
                        "Weisskopf units are quoted for the decay, B(E2↓), which needs the spins of the two "
                        "levels: they are not known for this transition, so only the ↑ value is given. ")
                     + f"The value the simulation was given is {number(result.truth_e2fm4)} "
                     f"e²fm⁴ ↑: the analysis closes on it at {result.pull:+.2f}σ. This is a simulated run "
                     "analysed back, so the closure checks the analysis, not the nucleus.</p>")
    else:
        parts.append("<p>B(E2↑) could not be determined: the analysis found no peak "
                     f"(the value the simulation was given is {number(result.truth_e2fm4)} e²fm⁴).</p>")
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
                 + "".join(f"<li><b>{e(st['step'])}</b>: {e(plain_numbers(st['text']))}</li>" for st in result.steps)
                 + "</ol>")
    if result.b_e2fm4 == result.b_e2fm4:
        parts.append("<h4>The uncertainty budget</h4><table><thead><tr><th>Source</th><th>Relative (%)</th></tr>"
                     f"</thead><tbody><tr><td>statistical</td><td>{100 * result.statistical:.2f}</td></tr>"
                     + "".join(f"<tr><td>{e(k)}</td><td>{100 * v:.2f}</td></tr>" for k, v in result.budget.items())
                     + f"<tr><td><b>total</b></td><td><b>{100 * result.total_unc:.2f}</b></td></tr></tbody></table>")
        if getattr(result, "budget_notes", None):
            parts.append("<ul>" + "".join(f"<li>{e(n[0].upper() + n[1:])}.</li>" for n in result.budget_notes) + "</ul>")
        if result.whole_run_statistical is not None:
            parts.append(f"<p>The whole run, scaled from its real part, would give a statistical uncertainty of "
                         f"{100 * result.whole_run_statistical:.1f}%.</p>")
        sh = result.shape or {}
        if sh:
            rows = [("β₂", f"{sh['beta2']:.3f}"),
                    ("B(E2↓)", f"{number(decay[1], 3)} W.u. ↓") if decay is not None else
                    ("B(E2↑)", f"{number(result.b_e2fm4)} e²fm⁴"),
                    ("Q₀ of a rigid rotor", f"{sh['q0_efm2']:.3g} e fm²"),
                    ("Q_s(2⁺) of that rotor", f"{sh['qs_2plus_efm2']:.3g} e fm²")]
            if "e4_over_e2" in sh:
                rows.append(("E(4⁺)/E(2⁺)", f"{sh['e4_over_e2']:.3f}"))
            parts.append("<h4>The shape of the nucleus</h4><table><thead><tr><th>Quantity</th><th>Value</th></tr>"
                         "</thead><tbody>" + "".join(f"<tr><td>{e(q)}</td><td>{e(v)}</td></tr>" for q, v in rows)
                         + "</tbody></table>" + (f"<p>{e(sh['note'])}</p>" if sh.get("note") else ""))
        parts.append(f"<p>Beam time: {number(result.counts_per_shift, 3)} counts in the peak per 8 h shift; "
                     f"{number(result.hours_for_precision, 3)} h for a statistical precision of "
                     f"{100 * float(result.settings.get('wanted_precision', 0.05)):.3g}%.</p>")
    if result.notes:
        parts.append("<ul>" + "".join(f"<li>{e(note_words(n))}</li>" for n in result.notes) + "</ul>")
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


def _shapes_html(sh: dict, beam_time_s: float) -> str:
    """Prolate against oblate in a paragraph: what is compared (the γ rays of the 2⁺ state's decay, summed over
    all rings in the planned beam time, with the quadrupole moment of a prolate and of an oblate rotor), the
    difference in percent, and that the number of standard deviations counts only the statistics."""
    q = sh["q_efm2"]
    text = (f"Prolate or oblate: the 2⁺ state with Q(2⁺) = {q['prolate']:+.1f} (prolate rotor), 0 and "
            f"{q['oblate']:+.1f} e fm² (oblate). ")
    totals = sh.get("totals") or {}
    sigma = sh["total_difference_sigma"]
    if totals.get("prolate") and totals.get("oblate") and totals.get("spherical"):
        share = 100 * abs(totals["prolate"] - totals["oblate"]) / totals["spherical"]
        text += (f"Compared are the γ rays of its decay summed over all rings in the planned beam time "
                 f"({_time(beam_time_s)}): {number(totals['prolate'])} for the prolate moment against "
                 f"{number(totals['oblate'])} for the oblate one, a difference of {number(share, 2)} % that comes "
                 "from the reorientation effect. ")
    else:
        text += "Compared are the γ rays of its decay summed over all rings in the planned beam time. "
    text += ("The planned run " + ("can" if sh["separable"] else "cannot") + " tell prolate from oblate by its "
             f"counts: the two differ by {number(sigma, 3)} standard deviations of the counting statistics. That "
             "figure counts the statistics only" + (", and with so many γ rays it says little" if sigma >= 10 else "")
             + ": how well the yield is known in absolute terms is not in it. The efficiency, the target and B(E2) "
             "itself enter the yield (the uncertainty budget of the result), and a different B(E2) moves the yield "
             "as a different Q(2⁺) does.")
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
                       + "".join(f"<td>{number(x['all_orders'][d] * t)} / {number(x['first_order'][d] * t)}</td>"
                                 for d in dets) + "</tr>" for x in m["transitions"]) + "</tbody></table>"]
    sh = m.get("shapes")
    if sh:
        parts.append(_shapes_html(sh, t))
    if m.get("notes"):
        parts.append("<ul>" + "".join(f"<li>{e(n)}</li>" for n in m["notes"][:4]) + "</ul>")
    return "".join(parts)


def _doppler_depth(experiment) -> str:
    """At what depth of the target the Doppler correction takes the reaction (for the nucleus the ``[reaction]``
    excites, :func:`~physim.nuclear.gamma_events.reaction_depths`), in words."""
    from .gamma_events import reaction_depths
    from .rates import stack

    layers = stack(experiment)
    depths = reaction_depths(experiment, layers)
    if depths is None:
        return ("It is the middle of the target: the "
                "excitation hardly changes through this target, or this is a lifetime measurement.")
    lo, hi = min(depths.values()), max(depths.values())
    value = f"{lo:.3g}" if f"{lo:.3g}" == f"{hi:.3g}" else f"{lo:.3g} to {hi:.3g}"
    return (f"It is one for each particle detector and kind of particle: the mean depth of the excitations that the "
            f"setup's cross sections give there, here {value} of {layers[0].thickness:.3g} mg/cm². A thick target "
            "excites mostly near its front, where the beam is fastest. (For the other nucleus of the collision it is "
            "the middle of the target, where its excitation is taken.)")


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
    depth = _doppler_depth(exp)
    parts = [f"<p>{e(dv.correction_reason(exp, emitter))} The spectra and the analysis correct for {who}.</p>",
             '<p class="f">E₀ = E_measured × (1 − β cos α) / √(1 − β²)</p>',
             "<p>with β the emitter's speed in units of c and α the angle between its direction and the "
             "direction of the γ ray. Neither is measured directly; each event's are worked out from the "
             "segment of the particle detector and the crystal that fired:</p><ol>"
             "<li><b>The particle's direction</b> is the mean direction, in the ring and sector (or the two strips) "
             "it hit, of the particles of its kind (scattered beam or recoil), seen from the target at the position "
             "the analysis assumes.</li>"
             f"<li><b>The depth</b> of the reaction. {e(depth)}</li>"
             "<li><b>The orbit.</b> Two-body kinematics at the beam energy there, with the "
             f"excitation energy of {1e3 * exc.energy_mev:g} keV taken from the motion, gives the CM angle of "
             "the collision from that direction, whichever of the two nuclei was detected.</li>"
             "<li><b>The emitter.</b> The same kinematics gives the direction and the energy of the nucleus "
             "that emits: the detected particle itself, or its partner on the other side of the beam.</li>"
             "<li><b>Its speed β</b> is that energy less what the nucleus loses on its way out of the target "
             "from that depth.</li>"
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
                     width=width, planner=planner), ""
    except Exception as err:  # noqa: BLE001 - the report stands without the plan's figures
        return None, str(err)


# -- numbers and names as the report prints them --------------------------------------------------------------------

_SUPERSCRIPT = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

#: The units of the ledger's rows, of the B(Eλ) values and of the matrix elements, as the report prints them.
UNIT_WORDS = {"e2fm2": "e²fm²", "e2fm4": "e²fm⁴", "e2fm6": "e²fm⁶", "e2b": "e²b", "e2b2": "e²b²", "e2b3": "e²b³",
              "efm": "e fm", "efm2": "e fm²", "efm3": "e fm³", "eb": "e b", "mu_N": "μ_N", "mu_N2": "μ_N²"}

#: What stands for an uncertainty that cannot be shown as a number: NaN, infinite or negative.
NOT_AVAILABLE = "n/a"


def unit_words(unit: str) -> str:
    """A unit of the ledger as the report prints it ("e2fm4" → "e²fm⁴")."""
    return UNIT_WORDS.get(unit, unit)


def _rounded(x: float, digits: int) -> tuple:
    """(mantissa, exponent of ten) of ``x`` rounded to ``digits`` significant figures: 9.99996e12 at four figures
    is (1.0, 13), not (10.0, 12)."""
    mantissa, exponent = f"{x:.{max(digits, 1) - 1}e}".split("e")
    return float(mantissa), int(exponent)


def number(x, digits: int = 4) -> str:
    """A number as the report prints it: ``digits`` significant figures, thousands grouped, and no exponent form
    for ordinary sizes (from 0.0001 up to ten million), so 437100 is "437,100" and not "4.371e+05".
    Beyond those sizes it is a power of ten ("3.2 × 10¹²"). A whole number of counts is printed whole."""
    if x is None:
        return "—"
    if isinstance(x, (int, np.integer)) and not isinstance(x, bool):
        return f"{int(x):,}"
    x = float(x)
    if not math.isfinite(x):
        return "—"
    if x == 0:
        return "0"
    mantissa, exponent = _rounded(x, digits)
    if -4 <= exponent < 7:
        places = digits - 1 - exponent
        text = f"{round(x, places):,.{max(places, 0)}f}"
        return text.rstrip("0").rstrip(".") if "." in text else text
    text = f"{mantissa:.{max(digits, 1) - 1}f}"
    text = text.rstrip("0").rstrip(".") if "." in text else text  # never the zero of "10"
    return f"{text} × 10{str(exponent).translate(_SUPERSCRIPT)}"


def matched(value, *uncertainties, digits: int = 2) -> tuple:
    """A value and its uncertainties as texts with the same number of decimals, and the power of ten that
    multiplies all of them: ``(value, *uncertainties, factor)``, where ``factor`` is "" for ordinary sizes and
    "× 10ⁿ" for a value below 0.0001 or from ten million on. Whoever prints the texts prints the factor
    (:func:`plus_minus`, :func:`asymmetric`, :func:`row_value` do).

    The decimals show the largest uncertainty to ``digits`` significant figures and keep at least one figure of
    the smallest: an uncertainty that is not zero is never printed as 0, so (919, 0.3, 64) is ("919.0", "0.3",
    "64.0", ""). There is no exponent form: (745, 130) is ("745", "130", ""), not "1.3e+02". An uncertainty that
    is missing (None) comes back as "—", one that is exactly zero as "0", and one that is no number that can be
    shown (NaN, infinite, negative) as :data:`NOT_AVAILABLE`. Without any uncertainty to show, the value has four
    significant figures."""
    def kind(u) -> str:
        if u is None:
            return "missing"
        u = float(u)
        if not math.isfinite(u) or u < 0:
            return "bad"
        return "zero" if u == 0 else "shown"

    kinds = [kind(u) for u in uncertainties]
    blank = {"missing": "—", "bad": NOT_AVAILABLE, "zero": "0"}
    if value is None or not math.isfinite(float(value)):
        return ("—",) + tuple("—" for _ in uncertainties) + ("",)
    sizes = [float(u) for u, k in zip(uncertainties, kinds) if k == "shown"]
    if not sizes:
        return (number(value),) + tuple(blank[k] for k in kinds) + ("",)
    value = float(value)

    def exponent(x: float) -> int:
        return math.floor(math.log10(abs(x)))

    scale = exponent(value) if value != 0 and not 1e-4 <= abs(value) < 1e7 else 0
    # Two figures of the largest uncertainty, and at least one of the smallest.
    places = max(digits - 1 - exponent(max(sizes)), -exponent(min(sizes))) + scale
    if scale and abs(round(value / 10.0 ** scale, max(places, 0))) >= 10:  # 9.99996 rounds up to 10: one power on
        scale, places = scale + 1, places + 1

    def text(x, least: int) -> str:
        return f"{round(float(x) / 10.0 ** scale, max(places, least)):,.{max(places, 0)}f}"

    # An uncertainty of hundreds keeps its two figures ("130"); the value beside it is not rounded past the units.
    return ((text(value, 0),) + tuple(text(u, places) if k == "shown" else blank[k]
                                      for u, k in zip(uncertainties, kinds))
            + (f"× 10{str(scale).translate(_SUPERSCRIPT)}" if scale else "",))


def _scaled(body: str, factor: str) -> str:
    """A value with its uncertainties under the power of ten that multiplies them all: "(5.0 ± 1.0) × 10⁻⁵"."""
    return f"({body}) {factor}" if factor else body


def plus_minus(value, *uncertainties, digits: int = 2) -> str:
    """"919 ± 13 ± 63": a value and its uncertainties to matching decimals (:func:`matched`), under their power
    of ten where they have one ("(2.50 ± 0.10) × 10⁻¹²")."""
    *parts, factor = matched(value, *uncertainties, digits=digits)
    if not any(u is not None for u in uncertainties):
        return parts[0]
    return _scaled(" ± ".join(parts), factor)


def asymmetric(value, plus, minus, digits: int = 2) -> str:
    """"2.09 +0.12 −0.10": a value with an upper and a lower uncertainty to matching decimals, under their power
    of ten where they have one."""
    v, up, down, factor = matched(value, plus, minus, digits=digits)
    return _scaled(f"{v} +{up} −{down}", factor)


def percent(x: Optional[float]) -> str:
    """A relative uncertainty in percent as the report prints it: "±7.0 %"."""
    if x is None or not math.isfinite(x):
        return "—"
    return f"±{x:.1f} %" if x >= 0.1 else f"±{number(x, 2)} %"


def nuclide_words(nuclide: str) -> str:
    """ "62Ni" → "⁶²Ni"."""
    m = re.fullmatch(r"(\d+)([A-Za-z]+)", str(nuclide))
    return f"{m.group(1).translate(_SUPERSCRIPT)}{m.group(2)}" if m else str(nuclide)


def level_words(key) -> str:
    """A level's key (:class:`~physim.nuclear.results.LevelKey` or its text) as the report prints it: "2⁺₁" for
    a level anchored to ENSDF, "2⁺ (1173 keV)" for one known by its energy, "2⁺ (calculated, no. 1)" for a
    calculation's own ordinal."""
    k = key if isinstance(key, results.LevelKey) else results.LevelKey.parse(str(key))
    jpi = k.jpi.replace("+", "⁺").replace("-", "⁻")
    if k.anchored:
        return jpi + str(k.ordinal).translate(_SUBSCRIPT)
    if k.ordinal is not None:
        return f"{jpi} (calculated, no. {k.ordinal})"
    return jpi if not k.energy_kev else f"{jpi} ({k.energy_kev:.0f} keV)"


def quantity_words(observable: str, key, direction: Optional[str] = None) -> str:
    """What a row is a value of, as the report prints it, with the direction of a B stated: "B(E2; 0⁺₁ → 2⁺₁)↑"
    for the excitation, "B(E2; 2⁺₁ → 0⁺₁)↓" for the decay, "τ(2⁺₁)", "Q_s(2⁺₁)", "E(2⁺₁)"."""
    o = observables.get(observable)
    if o.belongs_to == "transition":
        t = key if isinstance(key, results.TransitionKey) else results.TransitionKey.parse(str(key))
        direction = direction or o.direction
        first, second = (t.final, t.initial) if direction == observables.UP else (t.initial, t.final)
        return (f"B({t.multipolarity}; {level_words(first)} → {level_words(second)})"
                f"{observables.DIRECTIONS[direction]}")
    return f"{o.label}({level_words(key)})"


def row_value(row, canonical: bool = False, unit: bool = True, direction: bool = False) -> str:
    """A row's value with its uncertainties and unit, to matching decimals: "919 ± 13 (stat) ± 63 (sys) e²fm⁴",
    "884 ± 13 e²fm⁴", "2.09 +0.12 −0.10 ps", or the bare value for a row without an uncertainty. A value in
    Weisskopf units always carries its ↓; with ``direction`` any other B carries its ↑ or ↓ too. With
    ``canonical`` the row is first brought to its observable's convention; without ``unit`` the unit is left
    out (a table that names it in its heading)."""
    if canonical:
        c = compare.canonical(row)
        value, name, stat, sys, plus, minus, arrow = c.value, c.unit, c.stat, c.sys, c.plus, c.minus, c.direction
    else:
        value, name, stat, sys, arrow = row.value, row.unit, row.stat, row.sys, row.direction
        plus, minus = (row.plus, row.minus) if row.asymmetric else (row.unc(), row.unc())
    tail = ""
    if unit:
        tail = " " + unit_words(name)
        if name == observables.WEISSKOPF:
            tail += " ↓"
        elif direction and arrow:
            tail += " " + observables.DIRECTIONS[arrow]
    if row.asymmetric:
        return asymmetric(value, plus, minus) + tail
    if stat is not None and sys is not None:
        v, s, y, factor = matched(value, stat, sys)
        return _scaled(f"{v} ± {s} (stat) ± {y} (sys)", factor) + tail
    if plus is None:
        return number(value) + tail
    return plus_minus(value, plus) + (" (stat)" if stat is not None else " (sys)" if sys is not None else "") + tail


def _run_words(label) -> str:
    """A run's folder name as the report calls the run: "03-beam-24h" → "run 3"."""
    m = re.match(r"(\d+)-", str(label or ""))
    return f"run {int(m.group(1))}" if m else str(label or "a run")


def source_words(row, here: Optional[str] = None, by_id: Optional[dict] = None) -> str:
    """Where a row comes from, in a few words: "this experiment, run 3" (``here`` is this experiment's folder
    name), "62Ni-on-194Pt-3, run 1", "Allmond et al. 2014", "JUN45, f5pg9", "own measurement, 2026-11-20". A row
    converted from one of ``by_id`` (the nucleus's rows) is said with that row: "Kenn et al. 2001, DSAM (line
    shape), converted from τ = 2.01 ± 0.12 ps"."""
    o = row.origin or {}
    kind = o.get("kind")
    if kind == "run":
        where = "this experiment" if here is not None and o.get("experiment") == here else str(o.get("experiment"))
        several = len(o.get("runs") or [])
        if several > 1:  # a result of a series of runs (a lifetime) names them all
            return f"{where}, {several} runs"
        return f"{where}, {_run_words(o['run'])}" if o.get("run") else where
    if kind == "reference":
        return " ".join(str(x) for x in (o.get("authors"), o.get("year")) if x) or "a reference"
    if kind == "shell model":
        return ", ".join(str(x) for x in (o.get("interaction"), o.get("space")) if x) or "a shell-model calculation"
    if kind == "own measurement":
        return "own measurement" + (f", {o['date']}" if o.get("date") else "")
    made_from = (by_id or {}).get(o.get("from")) if kind == "conversion" else None
    if made_from is not None:
        return (f"{source_words(made_from, here)}, {results.METHODS[made_from.method]}, converted from "
                f"{observables.get(made_from.observable).label} = {row_value(made_from, direction=True)}")
    if kind == "conversion":
        return "converted from another value" + (f" ({row.detail})" if row.detail else "")
    return "not stated"


def note_words(note: str) -> str:
    """A note of an analysis as the report prints it: settings that were written into it as a Python dict
    ("{'correction': 'projectile'}") become words ("correction = projectile"), and a B(E2) given in Weisskopf
    units says that it is the ↓ value."""
    def words(m) -> str:
        try:
            d = ast.literal_eval(m.group(0))
        except (ValueError, SyntaxError):
            return m.group(0)
        if not isinstance(d, dict):
            return m.group(0)

        def name(k) -> str:
            return " ".join(str(x) for x in k) if isinstance(k, tuple) else str(k).replace("_", " ")

        return ", ".join(f"{name(k)} = {v}" for k, v in d.items()) or "nothing"

    return plain_numbers(wu_marked(re.sub(r"\{[^{}]*\}", words, str(note))))


def plain_numbers(text: str) -> str:
    """Text in which numbers written in exponent form ("4.381e+05", "1.910e-01") are written as the report
    writes numbers (:func:`number`), with the figures they had."""
    def plain(m) -> str:
        return number(float(m.group(0)), max(len(m.group(1).replace(".", "")), 1))

    return re.sub(r"(?<![\w.])(\d(?:\.\d+)?)e([+-]\d+)(?!\w)", plain, str(text))


def wu_marked(text: str) -> str:
    """Text in which a B(E2) quoted in Weisskopf units carries its direction: "B(E2) = 12.13 W.u." becomes
    "B(E2)↓ = 12.13 W.u." (a value in Weisskopf units is always the decay's)."""
    return re.sub(r"B\((E\d)\)( = [^;,()]*? W\.u\.)", r"B(\1)↓\2", str(text))


def _explanation_html(x: Explanation) -> str:
    """One explained number as a box (as :func:`physim.nuclear.record._explanation_html`), its value printed
    without an exponent form for ordinary sizes."""
    e = html.escape
    unit = f" {e(x.unit)}" if x.unit else ""
    return (f'<div class="x" id="{e(x.key)}"><h3>{e(x.title)}: <span class="v">{number(x.value)}{unit}</span></h3>'
            f'<div class="f">{e(x.formula)}</div><div class="f">{e(plain_numbers(x.substituted))}</div>'
            f'<div class="m">{e(x.meaning)}</div>'
            f'<div class="a">Assumptions: {e(x.assumptions)} <a href="{e(x.reference)}">Theory</a></div></div>')


# -- the results: the rows of the ledger, the comparison, what is still to run --------------------------------------

#: What a value in Weisskopf units says, for "What it means".
WEISSKOPF_MEANING = ("The Weisskopf unit is the strength one proton moving alone would give. Tens of Weisskopf "
                     "units mean many nucleons move together: a collective state. A few or less is a "
                     "single-particle transition.")

#: The rule of the verdicts, in one line (the numbers are :data:`physim.nuclear.compare.AGREE_BELOW` and
#: :data:`~physim.nuclear.compare.DIFFER_ABOVE`).
RULE = (f"The rule, fixed beforehand: the pull is the difference over the two uncertainties in quadrature; below "
        f"{compare.AGREE_BELOW:g} the two agree, from {compare.AGREE_BELOW:g} to {compare.DIFFER_ABOVE:g} there is "
        f"tension, above {compare.DIFFER_ABOVE:g} they differ. A value without an uncertainty gets its ratio and no "
        "verdict. A projected value is closed against the value its simulation was given, by the same rule, and "
        "is not pulled against the reference: a simulated run cannot test a measured value.")

_BOX = ("border:1px solid #D5D9D6;border-left:4px solid #0072B2;border-radius:0 6px 6px 0;padding:10px 16px;"
        "margin:10px 0;background:rgba(127,127,127,0.08);-webkit-print-color-adjust:exact;"
        "print-color-adjust:exact")
_TAG = ("display:inline-block;font:600 0.62em monospace;letter-spacing:0.04em;color:#8a4b00;"
        "border:1px solid #8a4b00;border-radius:4px;padding:0 6px;vertical-align:middle;margin-left:0.5em")
_PLAIN_LIST = "list-style:none;padding-left:0"


def _esc(text) -> str:
    return html.escape(str(text), quote=False)


def run_key(planner, summary: Optional[dict]):
    """The key of the transition the ``[reaction]`` of a run's own setup names (``summary`` is the run's
    summary; None: the current setup), or None for a setup that is not a Coulomb excitation or cannot be read."""
    from .runs import taken_setup

    try:
        exp = taken_setup(summary["_folder"], summary.get("kind")) if summary else planner.experiment
        return results.key_of_reaction(exp)
    except (OSError, ValueError, TypeError, KeyError):
        return None


def analysis_row(planner, summary: Optional[dict], result) -> tuple:
    """(row, why): the row a kept Coulomb-excitation analysis stands for, made in memory and not written
    (:func:`physim.nuclear.result_rows.coulex_row`, as the planner writes it after an analysis), or None with the
    reason no row of the results ledger can be made for it. ``summary`` is the run's summary (None for the
    analysis of a simulated sample without a run)."""
    from . import result_rows
    from .runs import taken_setup

    if result.b_e2fm4 != result.b_e2fm4:
        return None, "the analysis found no peak, so there is no B(E2) to enter"
    try:
        exp = taken_setup(summary["_folder"], summary.get("kind")) if summary else planner.experiment
    except (OSError, ValueError, TypeError, KeyError) as err:
        return None, f"the setup the run was taken with could not be read ({err})"
    if exp.excitation is None:
        return None, "the setup the run was taken with is not a Coulomb excitation"
    if exp.excitation.multipolarity != "E2":
        return None, (f"the results ledger has no observable for a B({exp.excitation.multipolarity}) yet: it "
                      "holds B(E2)")
    if result.truth_e2fm4 != result.truth_e2fm4:
        return None, ("the analysis kept no value that the simulation was given, and a projected row is closed "
                      "against that value")
    here = planner.experiment_folder.name if planner.folder is not None else planner.name
    try:
        return result_rows.coulex_row(result, exp, here, summary["label"] if summary else "the simulated sample"), ""
    except (ValueError, TypeError, KeyError) as err:
        return None, f"the row could not be made ({err})"


def experiment_rows(planner, runs: Optional[list] = None, analyses: Optional[dict] = None, result=None,
                    problems: Optional[dict] = None) -> list:
    """The rows this experiment has produced that hold now: those of its ledger (``results.jsonl`` in its
    folder), and a row made in memory for every kept analysis the ledger does not hold as it stands.

    Nothing is written: the planner writes a row at each analysis and backfills an old experiment when it is
    opened; the report only reads. ``runs`` and ``analyses`` are :meth:`Planner.runs` and :func:`run_analyses`
    if already at hand; ``result`` is the analysis of a simulated sample where there is no run. An analysis no
    row can be made for is said in ``problems`` if a dict is given: {run number (None: the simulated sample):
    the reason}."""
    runs = planner.runs() if runs is None else runs
    analyses = run_analyses(planner, result) if analyses is None else analyses
    rows = []
    folder = planner.experiment_folder
    if folder is not None:
        try:
            rows = results.of_experiment(folder)
        except (OSError, ValueError):
            rows = []
    rows = results.current(rows)
    by_number = {s["number"]: s for s in runs}
    pending = [(n, by_number.get(n), r) for n, r in analyses.items()]
    if not analyses and result is not None:
        pending.append((None, None, result))
    for n, summary, r in pending:
        label = summary["label"] if summary else None
        kept = [x for x in rows if x.observable == "B(E2)" and x.method == "coulex-yield"
                and x.origin.get("run") == label]
        b = r.b_e2fm4
        if any(abs(x.value - b) <= 1e-6 * abs(b) for x in kept):
            continue
        made, why = analysis_row(planner, summary, r)
        if made is not None:  # the ledger's row of this run, if any, is of an earlier analysis: this one holds
            rows = [x for x in rows if x not in kept] + [made]
        elif problems is not None:
            problems[n] = why
    return rows


def with_lifetime_rows(planner, analysis: Optional[dict], own: list) -> list:
    """``own`` (:func:`experiment_rows`) with the rows of a lifetime analysis where the ledger does not hold them
    as they stand: τ of the level and the B(E2)↓ that depends on it, made in memory and not written
    (:func:`physim.nuclear.result_rows.lifetime_rows`, as the planner writes them after an analysis). ``analysis``
    is the dict of :meth:`~physim.nuclear.planner.Planner.analyse_lifetimes`. ``own`` comes back as it is without
    an analysis, or where the setup holds no τ for the level to be the value the simulation was given."""
    from . import result_rows

    if not analysis or not analysis.get("available"):
        return own
    here = planner.experiment_folder.name if planner.folder is not None else planner.name
    try:
        made = result_rows.lifetime_rows(analysis, here)
    except (ValueError, TypeError, KeyError):
        made = []
    if not made:
        return own
    tau = made[0]
    standing = [r for r in own if r.observable == "tau" and r.key == tau.key and r.method == tau.method]
    if any(abs(r.value - tau.value) <= 1e-6 * abs(tau.value) for r in standing):
        return own  # the ledger holds this analysis
    try:
        made += result_rows.lifetime_rows(analysis, here, depends_on=tau)
    except (ValueError, TypeError, KeyError):
        pass
    old = {r.id for r in standing}
    return [r for r in own if r.id not in old and not (old & set(r.depends_on))] + made


def nucleus_rows(planner, nuclide: str, own: list = ()) -> list:
    """The rows of a nucleus that hold now, which the report compares: those of the experiments beside this one
    and of the nucleus's own points in the nuclei home (:func:`physim.nuclear.results.of_nuclide`,
    :func:`~physim.nuclear.results.nuclei_home`), with this experiment's rows as ``own`` has them (which may
    hold rows that are not written, and leaves out rows a newer analysis replaced)."""
    folder = planner.experiment_folder
    name = results.nuclide_name(nuclide)
    try:
        rows = results.of_nuclide(name, folder.parent if folder is not None else None)
    except (OSError, ValueError):
        rows = []
    if folder is not None:  # this experiment's rows are taken from ``own``
        rows = [r for r in rows if not (r.origin.get("kind") == "run" and r.origin.get("experiment") == folder.name)]
    seen = {r.id for r in rows}
    # What holds now is decided over all the nucleus's rows, before a quantity is picked: a corrected data point
    # may have moved to another level or transition.
    return results.current(rows + [r for r in own if r.nuclide == name and r.id not in seen])


def observable_rows(planner, nuclide: str, observable: str, key: str, own: list = ()) -> tuple:
    """(view, rows, by_id, why not): what the nucleus page shows of one observable of one level or transition
    (:func:`physim.nuclear.nucleus.observable_view`), from the rows the report has (:func:`nucleus_rows`, with this
    experiment's ``own``), and its rows: those that hold and the B(E2) or τ values that others give through the
    level's γ ray. ``by_id`` are the nucleus's rows, which say what a converted row was made from. Nothing is
    written. Where the rows cannot be compared (a B↓ between levels of unknown spin) the view holds the routes
    alone and the reason comes back in words."""
    ledger = nucleus_rows(planner, nuclide, own)
    rows, via = nucleus.with_converted(nuclide, ledger)
    group = results.select(rows, observable, key)
    by_id = {r.id: r for r in ledger}
    try:
        return nucleus.observable_view(nuclide, observable, key, ledger=ledger), group, by_id, ""
    except ValueError as err:
        routes = nucleus.routes(nuclide, observable, key, group, via, rows)
        return {"comparison": None, "routes": routes}, group, by_id, str(err)


def evidence_missing(rows: list, outside: bool = False, routes: list = ()) -> list:
    """The independent evidence for one observable of one level or transition, kind by kind: [{"evidence",
    "label", "have" (the rows that hold now), "missing" (what to say when there are none)}]. With ``outside`` the
    result of an outside program (GOSIA on the same yields) is asked for too. ``routes`` are those of the nucleus
    view (:func:`physim.nuclear.nucleus.routes`): the ones known from the literature only are named."""
    now = [r for r in results.current(list(rows)) if r.independent]
    literature = ", ".join(x if x[:2].isupper() or not x[0].isalpha() else x[0].lower() + x[1:]
                           for x in (r["label"] for r in routes if r["means"] == "literature"))
    kinds = [("evaluated", "Evaluated value", "No evaluated value (ENSDF adopted) is entered: without it there is "
              "no reference, and nothing here says how well the value is already known."),
             ("measured", "Measured points", "No measured point is entered: neither a published value (by any "
              "route above" + (f", or by one known from the literature only: {literature}" if literature else "")
              + ") nor one of your own."),
             ("theory", "Shell model", "No shell-model value is entered for this nucleus.")]
    if outside:
        kinds.append(("outside program", "Outside program", "No result of an outside program (GOSIA on the same "
                      "yields) has been read back."))
    return [{"evidence": k, "label": label, "have": [r for r in now if r.evidence == k], "missing": missing}
            for k, label, missing in kinds]


def closure_words(p) -> str:
    """What a projected row's closure says, in a sentence."""
    given = (f"the value the simulation was given ({number(p.truth)} "
             f"{unit_words(observables.get(p.row.observable).unit)})")
    if p.closure is None:
        return f"The result has no uncertainty, so it cannot be closed against {given}."
    how = f"{abs(p.closure):.1f}σ ({p.closure:+.2f}σ)"
    if p.closure_verdict == "agree":
        return f"The analysis gives back {given} within {how}: the chain from the spectra to the value closes."
    if p.closure_verdict == "tension":
        return (f"The analysis is {how} from {given}: a tension. A second run, or the same run analysed with "
                "another gate, shows whether it is a fluctuation.")
    return (f"The analysis is {how} from {given}: it does not close. Something in the chain from the spectra to "
            "the value is biased for this setup; the result is not to be read as the setup's reach until it does.")


def known_words(p, c, by_id: Optional[dict] = None) -> str:
    """ "This setup reaches ±x %; the value is known to ±y %; the best single measurement to ±z %." """
    reach = (f"This setup reaches {percent(p.reach_percent)}" if p.reach_percent is not None
             else "This setup's reach is not stated (the result has no uncertainty)")
    if p.known_percent is not None:
        known = (f"the value is known to {percent(p.known_percent)} "
                 f"({source_words(c.reference.row, by_id=by_id)}, the reference)")
    elif c.reference is not None:
        known = "the reference has no uncertainty, so how well the value is known is not stated"
    else:
        known = "no evaluated value is entered, so how well the value is already known is not stated here"
    if p.best_measured_percent is not None:
        best = (f"the best single measurement reaches {percent(p.best_measured_percent)} "
                f"({source_words(c.best_measured.row, by_id=by_id)})")
    else:
        best = "no measured point with an uncertainty is entered"
    return f"{reach}; {known}; {best}."


def comparison_html(c, here: Optional[str] = None, mine: Optional[str] = None, rule: bool = True,
                    by_id: Optional[dict] = None) -> str:
    """A comparison (:func:`physim.nuclear.compare.compare`) as the report prints it: the independent evidence
    against the reference with difference, pull and verdict; then, apart, the projected rows with their closure
    against the value their simulation was given and what they reach, set beside how well the value is known;
    and the rule in one line (``rule``; a second comparison on the same page leaves it out). ``here`` is this
    experiment's folder name; ``mine`` the id of the row the report leads with (its sentence comes first). A row
    that was brought to the comparison's convention says what it was given as, and one that was converted from
    another observable (a B(E2) from a τ) what it was made from, where ``by_id`` (the nucleus's rows) has it."""
    e = _esc
    unit = unit_words(c.unit)

    def source(row) -> str:
        return source_words(row, here, by_id)

    parts = []
    if c.entries:
        ref = c.reference
        parts.append("<h3>Independent evidence</h3><p>"
                     + (f"Measured, evaluated and calculated values, against the reference "
                        f"({e(source(ref.row))}, {e(results.EVIDENCE[ref.evidence].lower())})."
                        if ref is not None else
                        "Measured, evaluated and calculated values. No evaluated value is entered, so there is no "
                        "reference: they are listed without differences.") + "</p>")
        rows = []
        for x in c.entries:
            value = e(plus_minus(x.value, x.plus) if x.plus == x.minus else asymmetric(x.value, x.plus, x.minus))
            given = f"<br><small>given as {e(row_value(x.row, direction=True))}</small>" if x.converted else ""
            if x.is_reference:
                against, says = "the reference", ""
            elif ref is None:
                against, says = "—", ""
            elif x.pull is None:
                why = "another evaluation" if compare.averages(x.row) else "no uncertainty"
                against = "—" if x.ratio is None else f"× {x.ratio:.2f} ({why}: no pull)"
                says = ""
            else:
                against = f"{x.difference_percent:+.1f} %, pull {x.pull:+.2f}"
                says = compare.VERDICTS[x.verdict].lower()
            if x.independent_of_reference is False:
                says += ("; " if says else "") + "not independent of the reference"
            rows.append(f"<tr><td>{e(results.EVIDENCE[x.evidence])}</td><td>{e(source(x.row))}</td>"
                        f"<td>{e(x.label)}" + (f"<br><small>{e(x.detail)}</small>" if x.detail else "")
                        + f"</td><td>{value}{given}</td><td>{e(percent(x.uncertainty_percent))}</td>"
                        f"<td>{e(against)}</td><td>{e(says)}</td></tr>")
        parts.append("<table><thead><tr><th>Evidence</th><th>From</th><th>Method</th>"
                     f"<th>{e(c.label)} ({e(unit)})</th><th>Known to</th><th>Against the reference</th>"
                     "<th>Verdict</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")
    else:
        parts.append("<h3>Independent evidence</h3><p>None is entered for this quantity: no evaluated value, no "
                     "measured point, no calculation. What is missing is listed in the next section.</p>")
    if c.projected:
        parts.append("<h3>Projected: simulated runs analysed back</h3><p>Each is closed against the value its "
                     "simulation was given, and its uncertainty is what the setup would reach. None is pulled "
                     "against the reference.</p>")
        rows = []
        for p in c.projected:
            closure = "—" if p.closure is None else (
                f"{p.closure:+.2f}σ: " + {"agree": "closes", "tension": "tension", "differ": "does not close"}[
                    p.closure_verdict])
            rows.append(f"<tr><td>{e(source(p.row))}</td><td>{e(p.label)}"
                        + (f"<br><small>{e(p.detail)}</small>" if p.detail else "")
                        + f"</td><td>{e(row_value(p.row, canonical=True, unit=False))}"
                        + (f"<br><small>given as {e(row_value(p.row, direction=True))}</small>" if p.converted
                           else "") + "</td>"
                        f"<td>{e(number(p.truth))}</td><td>{e(closure)}</td>"
                        f"<td>{e(percent(p.reach_percent))}</td></tr>")
        parts.append("<table><thead><tr><th>Simulated run</th><th>Method</th>"
                     f"<th>{e(c.label)} ({e(unit)})</th><th>Given to the simulation</th><th>Closure</th>"
                     "<th>Reaches</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")
        lead = next((p for p in c.projected if p.id == mine), c.projected[-1])
        parts.append(f"<p>{e(known_words(lead, c, by_id))}</p>")
    if rule:
        parts.append(f"<p><small>{e(RULE)}</small></p>")
    return "".join(parts)


def to_run_html(routes: list, rows: list, here: Optional[str] = None, unanalysed: Optional[list] = None,
                outside: bool = False, analyse_words: str = "", by_id: Optional[dict] = None) -> str:
    """The comparisons still to run for one observable of one level or transition, as three lists: the experiments
    among the ``routes`` of the nucleus view (:func:`physim.nuclear.nucleus.routes`, the states the nucleus page
    shows), the runs of this experiment that are not analysed (``unanalysed``: their numbers), and the independent
    evidence that is missing. ``rows`` are the rows of the observable, converted ones among them (``by_id`` says
    what each was made from). A route is ticked (☑) only when an independent value by it exists; one with nothing
    but simulated runs is marked projected (◐), one without a row is open (☐), and one that does not apply to this
    nucleus is crossed out (✗); the reason is said with it."""
    e = _esc
    held = {r.id: r for r in rows}

    def item(mark: str, title: str, text: str) -> str:
        return f"<li>{mark} <b>{e(title)}.</b> {e(text)}</li>"

    def have(found: list) -> str:
        return "; ".join(f"{source_words(r, here, by_id)} ({row_value(r, direction=True)}"
                         f"{', projected' if r.evidence == results.PROJECTED else ''})" for r in found)

    parts = [f'<h3>Other routes to it</h3><ul style="{_PLAIN_LIST}">']
    for r in routes:
        if r["kind"] is None:  # a calculation, the evaluation, a literature method: in the evidence below
            continue
        mine = [held[i] for i in r["ids"] if i in held]
        done = [x for x in mine if x.independent]
        projected = [x for x in mine if x.evidence == results.PROJECTED]
        how = (r["note"][0].upper() + r["note"][1:] + ". ") if r["note"] else ""
        kind = "" if r["kind_label"] == r["label"] else f" ({r['kind_label']})"
        says = f" For this nucleus: {r['reason'].rstrip('.')}." if r["reason"] else ""
        doubt = says if r["applies"] is not True else ""  # beside a value, only what speaks against the route
        if r["status"] == nucleus.DONE:
            parts.append(item("☑", r["label"], f"{how}Done: {have(done)}."
                              + (f" Also projected: {have(projected)}." if projected else "") + doubt))
        elif r["status"] == nucleus.ONLY_PROJECTED:
            parts.append(item("◐", r["label"], f"{how}Projected: a simulated run, no measurement yet: "
                                               f"{have(projected)}.{doubt}"))
        elif r["applies"] is False:
            parts.append(item("✗", r["label"], f"{how}Does not apply: {r['reason'].rstrip('.')}."))
        elif r["live"]:
            parts.append(item("☐", r["label"], f"{how}No value by this route yet. physim has this kind of "
                                               f"experiment{kind}.{says}"))
        else:
            parts.append(item("☐", r["label"], f"{how}No value by this route yet, and physim does not have "
                                               f"this kind of experiment yet{kind}."))
    parts.append("</ul>")
    if unanalysed is not None:
        parts.append(f'<h3>Runs of this experiment</h3><ul style="{_PLAIN_LIST}">')
        if unanalysed:
            many = len(unanalysed) != 1
            names = ", ".join(str(n) for n in unanalysed[:-1]) + (" and " if many else "") + str(unanalysed[-1])
            parts.append(item("☐", f"Run{'s' * many} {names}", f"Taken and not analysed{analyse_words}: "
                                                              "nothing yet shows that the result repeats."))
        else:
            parts.append(item("☑", "Every beam run is analysed", "No run waits for its analysis."))
        parts.append("</ul>")
    parts.append(f'<h3>Independent evidence</h3><ul style="{_PLAIN_LIST}">')
    for x in evidence_missing(rows, outside, routes):
        parts.append(item("☑", x["label"], f"Entered: {have(x['have'])}.") if x["have"]
                     else item("☐", x["label"], x["missing"]))
    parts.append("</ul>")
    return "".join(parts)


def budget_html(row, extra: Optional[list] = None) -> str:
    """A row's uncertainty budget as bars, the largest term first: the statistical uncertainty and each term of
    the systematic one in percent of the value, and the total. ``extra`` are lines to print under it."""
    e = _esc
    size = abs(row.value)
    terms = [(name, 100.0 * v / size) for name, v in (row.budget or {}).items() if v is not None and size]
    if row.stat is not None and size:
        terms.append(("statistics", 100.0 * row.stat / size))
    if not terms:
        return ""
    terms.sort(key=lambda t: -t[1])
    top = max(t[1] for t in terms) or 1.0
    bar = ("height:0.75em;background:#0072B2;-webkit-print-color-adjust:exact;print-color-adjust:exact;"
           "min-width:1px;width:{:.0f}%")
    out = ['<table style="width:100%;max-width:640px"><thead><tr><th>Source</th><th style="width:55%"></th>'
           "<th>Of the value</th></tr></thead><tbody>"]
    for name, pc in terms:
        out.append(f"<tr><td>{e(name[0].upper() + name[1:])}</td><td><div style=\"{bar.format(100 * pc / top)}\">"
                   f"</div></td><td>{pc:.1f} %</td></tr>")
    total = row.relative
    if total is not None:
        out.append(f"<tr><td><b>Total</b></td><td></td><td><b>{100 * total:.1f} %</b></td></tr>")
    out.append("</tbody></table>")
    if row.stat is not None and row.sys is not None and size:
        stat, sys = 100 * row.stat / size, 100 * row.sys / size
        out.append("<p>" + (f"The systematic uncertainties ({sys:.1f} %) set the precision, not the counts "
                            f"({stat:.1f} %): more beam time alone would not improve it." if sys > stat else
                            f"The counts ({stat:.1f} %) set the precision, not the systematic uncertainties "
                            f"({sys:.1f} %): more beam time would improve it.") + "</p>")
    for line in extra or []:
        out.append(f"<p>{e(line)}</p>")
    return "".join(out)


def _headline_html(row, runs_words: str, method_words: str, of: str = "") -> str:
    """The box the report opens with: the quantity, its value with its uncertainties, the method and the runs
    in a sentence, and, for a projected row, plainly that it comes from a simulated run."""
    e = _esc
    several = len(row.origin.get("runs") or []) > 1
    out = [f'<div style="{_BOX}"><div style="font-size:1.3em">'
           f"{e(quantity_words(row.observable, row.key, row.direction))} = {e(row_value(row))}"]
    if row.evidence == results.PROJECTED:
        out.append(f'<span style="{_TAG}">PROJECTED</span>')
    out.append(f"</div><div style=\"margin-top:4px\"><b>Method:</b> {e(method_words)}. <b>From:</b> "
               f"{e(runs_words)}.")
    if row.evidence == results.PROJECTED:
        reach = f" ({percent(100 * row.relative)})" if row.relative is not None else ""
        out.append((" <b>Projected from simulated runs:</b> the runs were simulated and analysed back" if several
                    else " <b>Projected from a simulated run:</b> the run was simulated and analysed back")
                   + f", so this is the precision the setup would reach{e(reach)}, not a measurement"
                   + (f" of {e(of)}" if of else "") + ".")
    out.append("</div></div>")
    return "".join(out)


def decay_of(b_up: float, key) -> Optional[tuple]:
    """(B(E2)↓ in e²fm⁴, the same in Weisskopf units, B(E2)↑ over B(E2)↓) of a transition from its B(E2)↑ in
    e²fm⁴: the ratio from the spins of the two levels (:func:`physim.nuclear.observables.b_up_over_down`), the
    Weisskopf unit from the mass number (:func:`~physim.nuclear.observables.weisskopf`). None where ``key`` (a
    :class:`~physim.nuclear.results.TransitionKey` or its text) is missing or does not hold both spins: no ratio
    is assumed."""
    if key is None:
        return None
    try:
        t = key if isinstance(key, results.TransitionKey) else results.TransitionKey.parse(str(key))
        if t.initial.spin is None or t.final.spin is None:
            return None
        ratio = observables.b_up_over_down(t.final.spin, t.initial.spin, t.multipolarity)
        down = b_up / ratio
        return down, observables.weisskopf(down, t.nuclide, t.multipolarity), ratio
    except (ValueError, TypeError, ZeroDivisionError):
        return None


def decay_words(row) -> str:
    """A B(E2) row's decay value in a sentence: B(E2)↓ in e²fm⁴ and in Weisskopf units, which are always ↓
    (:func:`decay_of`). Empty where the spins of the two levels are not known."""
    if row.observable != "B(E2)":
        return ""
    try:
        up = compare.canonical(row).value
    except ValueError:
        return ""
    decay = decay_of(up, row.key)
    if decay is None:
        return ""
    return (f"For the decay, {quantity_words('B(E2)', row.key, observables.DOWN)} = {number(decay[0])} e²fm⁴ = "
            f"{number(decay[1], 3)} W.u. ↓"
            + (f"; in other units the excitation value is {number(up * 1e-4)} e²b² ↑" if row.direction
               == observables.UP else "") + ".")


def _pick(rows: list, current_run: Optional[str]):
    """The row the report leads with among an experiment's rows of one quantity: the current run's, else the
    newest."""
    if not rows:
        return None
    mine = [r for r in rows if current_run is not None and r.origin.get("run") == current_run]
    return max(mine or rows, key=lambda r: r.created)


def _run_describe(summary: Optional[dict]) -> str:
    if not summary:
        return "a simulated sample of events, without a run"
    text = str(summary.get("describe") or "")
    return f"run {summary['number']}" + (f" ({text[0].lower() + text[1:]})" if text else "")


def weisskopf_shown(x: Explanation, key) -> Optional[Explanation]:
    """The explanation of B(E2) in Weisskopf units as the report shows it: for the decay (the analysis holds
    the excitation's value over the Weisskopf unit; a value in Weisskopf units is always the decay's). None
    where the transition's ``key`` does not hold the spins of both levels: no ratio is assumed."""
    b_up, b_w = x.inputs["b"], x.inputs["b_w"]
    decay = decay_of(b_up, key)
    if decay is None:
        return None
    down, wu, ratio = decay
    return Explanation(
        x.key, "B(E2)↓ in Weisskopf units", wu, "W.u. ↓",
        "B(E2)↓ = B(E2)↑ × (2J_lower + 1) / (2J_upper + 1);  B_W(E2) = (1/4π) (3/5)² (1.2 A^⅓ fm)⁴ e²;  B↓ / B_W",
        f"B(E2)↓ = {number(b_up)} / {number(ratio)} = {number(down)} e²fm⁴;  B_W = {number(b_w)} e²fm⁴;  "
        f"{number(down)} / {number(b_w)} = {number(wu, 3)} W.u. ↓",
        x.meaning + " Weisskopf units are quoted for the decay (↓), as ENSDF quotes them.", x.assumptions,
        x.reference, {"b": b_up, "ratio": ratio, "b_w": b_w}, lambda b, ratio, b_w: b / ratio / b_w)


#: The sections that lead the report, in order: (anchor, title). A section with nothing to say is left out
#: (an elastic or fusion-evaporation experiment has no result in the ledger: it opens with what it gives).
RESULTS = (("result", "The result"), ("routes", "Beside the other routes"),
           ("to-run", "Comparisons still to run"), ("limits", "What limits it"), ("meaning", "What it means"))
#: The appendices, in order: (anchor, title). The report as it was before it led with the results, kept whole.
APPENDICES = (("setup", "The setup"), ("plan", "The plan: what this setup should give"),
              ("runs", "The runs and their spectra"), ("gates", "The gates"),
              ("doppler", "The Doppler correction"), ("analysis", "The analysis"),
              ("numbers", "Every number, explained"), ("method", "The data, the method and what is left out"))
#: Every section of the report, in order: the results, numbered, then the appendices, lettered.
SECTIONS = RESULTS + APPENDICES


def reaction_words(planner) -> str:
    """The reaction of the setup in a line, for the head of the report: Coulomb excitation with the state,
    fusion-evaporation with the studied channel, or elastic scattering; and the lifetime method where the setup
    measures a lifetime."""
    exp = planner.experiment
    if exp.excitation is not None:
        what = (f"Coulomb excitation of the {exp.excitation.excite} to {1e3 * exp.excitation.energy_mev:g} keV "
                f"({exp.excitation.multipolarity})")
    elif exp.fusion is not None:
        what = "fusion-evaporation"
        try:
            ch = exp.fusion.channels[planner.studied_channel()]
            what += f", studied channel {ch.particles} to {ch.residue_name}"
        except (TypeError, IndexError, AttributeError):
            pass
    else:
        what = "elastic scattering"
    if exp.lifetime is not None:
        what += ("; lifetime measurement with a plunger (RDDS)" if exp.lifetime.method == "rdds"
                 else "; lifetime measurement from line shapes (DSAM)")
    return what


def reaction_b_source(experiment) -> str:
    """Where the B(Eλ↑) of the setup's ``[reaction]`` comes from: the source of the level scheme's matrix
    element for that transition when the scheme has one that gives the same B (ENSDF's value, read or derived,
    as the Physics step puts it into the ``[reaction]``); "user" when the setup holds a value of its own."""
    exc = experiment.excitation
    try:
        scheme = experiment.levels.get("target" if exc.excite == "target" else "beam")
        if scheme is None:
            return "user"
        for i, lev in enumerate(scheme.levels):
            if not i or abs(lev.energy.value - 1e3 * exc.energy_mev) >= 1.0:
                continue
            me, b = scheme.matrix_element(0, i, exc.multipolarity), scheme.b(0, i, exc.multipolarity)
            if me is None or not b:
                continue
            if abs(b - exc.b_up_e2fm) <= 5e-3 * exc.b_up_e2fm:
                return me.source + (f" ({me.note})" if me.note else "") + ": the level scheme's matrix element"
            return f"user (the level scheme's matrix element gives {number(b)}; {me.source})"
    except (ValueError, TypeError, AttributeError):
        pass
    return "user"


def predicted_html(c: dict) -> str:
    """Predicted against measured (:meth:`~physim.nuclear.planner.Planner.compare`) as a table: the plan's
    number, what the run measured with its counting uncertainty to matching decimals, the difference in percent
    and in that uncertainty, and whether the two agree; then what is compared, and why a row differs."""
    e = _esc
    rows = []
    for r in c["rows"]:
        p, m, u = r["predicted"], r["measured"], r["uncertainty"]
        measured = "not measured" if m is None else plus_minus(m, u)
        difference = "—"
        if p and m is not None and r["agree"] is not None:  # no difference where the two are not comparable
            difference = f"{100 * (m - p) / p:+.1f} %" + (f", {(m - p) / u:+.1f}σ" if u else "")
        says = {True: "agree", False: "differ"}.get(r["agree"], "" if m is None else "shown, not compared")
        rows.append(f"<tr><td>{e(r['quantity'])}</td><td>{e(r['what'])}</td><td>{e(number(p))}</td>"
                    f"<td>{e(measured)}</td><td>{e(difference)}</td><td>{says}</td></tr>")
    out = ["<table><thead><tr><th>Quantity</th><th></th><th>Predicted by the plan</th><th>Measured in the run"
           "</th><th>Difference</th><th></th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"]
    if any(r["agree"] is not None for r in c["rows"]):
        out.append("<p>The plan's number was computed from the setup when the run was started; the run's is "
                   "counted in its event-by-event part. σ is the run's counting uncertainty alone, and a row "
                   "differs beyond 3σ."
                   + (" A row that differs: the plan's own accuracy is not in σ (it takes one typical energy, "
                      "angle and depth where the events have their spread), so with many counts a difference of "
                      "a percent already is many σ. The percent column says how large the difference is."
                      if any(r["agree"] is False for r in c["rows"]) else "") + "</p>")
    for note in c.get("notes") or []:
        out.append(f"<p>{e(note)}</p>")
    return "".join(out)


def _lifetime_level(lifetime: Optional[dict], name: Optional[str] = None) -> Optional[tuple]:
    """(nuclide, level key, the plan's row) of the level a lifetime experiment studies: the one called ``name``
    in the plan, else the plan's first. None where the plan names no nucleus or level."""
    plan = (lifetime or {}).get("plan") or {}
    levels = plan.get("levels") or []
    level = next((x for x in levels if x.get("level") == name), levels[0] if levels and name is None else None)
    if level is None or not plan.get("nucleus"):
        return None
    m = re.match(r"(\d+)(/2)?([+-])?", str(level["level"]))
    spin = None if m is None else int(m.group(1)) / (2.0 if m.group(2) else 1.0)
    parity = None if m is None else {"+": 1, "-": -1}.get(m.group(3))
    try:
        return plan["nucleus"], results.level_key(plan["nucleus"], level["energy_kev"], spin, parity), level
    except (ValueError, KeyError):
        return None


def _counted(pl: dict) -> str:
    """What the plan's counts wanted are counts of, in words."""
    return {"all": "counts", "excitations": "Coulomb-excitation events",
            "coincidences": "particle–γ coincidences"}.get(pl.get("measured"), str(pl.get("measured")))


def _here(planner) -> Optional[str]:
    return planner.experiment_folder.name if planner.folder is not None else None


def _other_results_html(own: list, lead: list) -> str:
    """This experiment's rows of other quantities than the ones it leads with, one line each."""
    ids = {r.id for r in lead}
    keys = {(r.observable, r.key) for r in lead}
    rest = [r for r in own if r.id not in ids and (r.observable, r.key) not in keys]
    if not rest:
        return ""
    return ("<p>This experiment's other results:</p><ul>" + "".join(
        f"<li>{_esc(quantity_words(r.observable, r.key, r.direction))} = {_esc(row_value(r))}, "
        f"{_esc(results.METHODS[r.method].lower())}"
        + (", projected" if r.evidence == results.PROJECTED else "") + "</li>" for r in rest) + "</ul>")


def _routes_section(sec: dict, planner, nuclide: str, observable: str, key: str, own: list, lead=None,
                    unanalysed: Optional[list] = None, outside: bool = False, to_run: bool = True) -> Optional[object]:
    """Fill "Beside the other routes" and (with ``to_run``) "Comparisons still to run" for one quantity, from the
    view the nucleus page shows (:func:`observable_rows`); returns the comparison."""
    e = _esc
    here = _here(planner)
    o = observables.get(observable)
    view, rows, by_id, why = observable_rows(planner, nuclide, observable, key, own)
    c = view["comparison"]
    what = quantity_words(observable, key)
    made = any(r.method == "converted" for r in rows)
    if c is not None:
        sec["routes"].append(f"<p>Every value of {e(what)} of {e(nuclide_words(nuclide))} that is entered: those of "
                             "the experiments beside this one and the nucleus's own data points, in one "
                             f"convention ({e(o.convention)})."
                             + (" A value that says it is converted was made from the level's "
                                f"{'B(E2)' if observable == 'tau' else 'τ'}, through its γ ray: it is the same "
                                "measurement said in another way, not one more." if made else "") + "</p>"
                             + comparison_html(c, here, getattr(lead, "id", None), rule=to_run, by_id=by_id))
    elif why:
        sec["routes"].append(f"<p>The values of {e(what)} could not be set beside each other: {e(why)}.</p>")
    if to_run:
        sec["to-run"].append(f"<p>What would test or sharpen {e(what)} of {e(nuclide_words(nuclide))}, and is not "
                             "done.</p>" + to_run_html(view["routes"], rows, here, unanalysed, outside,
                                                       by_id=by_id))
    return c


def _coulex_results(sec: dict, planner, runs: list, analyses: dict, result, pl) -> None:
    """The leading sections of a Coulomb-excitation experiment: B(E2) of the ``[reaction]``'s transition."""
    from .record import analysis_explanations

    e = _esc
    exp = planner.experiment
    problems: dict = {}
    own = experiment_rows(planner, runs, analyses, result, problems)
    try:
        key, nuclide = results.key_of_reaction(exp).text, results.reaction_nuclide(exp)
    except ValueError:
        key, nuclide = None, None
    mine = [r for r in own if r.observable == "B(E2)" and r.evidence == results.PROJECTED]
    mine = [r for r in mine if r.key == key] or mine
    by_number = {s["number"]: s for s in runs}
    current = by_number.get(planner.run.number, {}).get("label") if planner.run is not None else None
    head = _pick(mine, current)
    beam_runs = [s for s in runs if s["kind"] != "source"]
    unanalysed = [s["number"] for s in beam_runs if s["number"] not in analyses] if beam_runs else None
    # An analysed run (or sample) the ledger can hold no row for: say which, and why.
    analysed = dict(analyses) if analyses else ({None: result} if result is not None else {})
    without = []
    for n, why in problems.items():
        r, which = analysed[n], "The simulated sample" if n is None else f"Run {n}"
        without.append(f"{which}: the analysis found no peak, so B(E2) could not be determined."
                       if r.b_e2fm4 != r.b_e2fm4 else
                       f"{which} is analysed (B(E2↑) = {number(r.b_e2fm4)} e²fm⁴, in appendix F), but no row of the "
                       f"results ledger can be made for it: {why}.")
    if head is None:
        if analysed:
            sec["result"].append("<p><b>No result to lead with.</b> " + e(" ".join(without) or
                                 "The analysed runs have no B(E2) row in the results ledger.") + "</p>")
        elif beam_runs:
            sec["result"].append(f"<p><b>No result yet.</b> {len(beam_runs)} beam run{'s' * (len(beam_runs) != 1)} "
                                 "taken and none analysed: analyse one on the Analysis step, and its B(E2) leads "
                                 "this report.</p>")
        else:
            sec["result"].append("<p><b>No result yet.</b> No run has been taken with this setup: the report "
                                 "holds what the setup should give (appendix B).</p>")
        if pl:
            sec["result"].append(f"<p>The plan: {_time(pl['beam_time_needed_s'])} of beam for "
                                 f"{pl['counts_wanted']} {e(_counted(pl))}"
                                 + (f" in {e(str(pl['limiting_detector']))}" if pl["limiting_detector"] else "")
                                 + f"; {_time(pl['beam_time_s'])} are planned.</p>")
        if key is not None:
            _routes_section(sec, planner, nuclide, "B(E2)", key, own, None, unanalysed, outside=True)
        return
    nuclide, key = head.nuclide, head.key
    summary = next((s for s in runs if s["label"] == head.origin.get("run")), None)
    r = analyses.get(summary["number"]) if summary else (result if not analyses else None)
    method = results.METHODS[head.method] + (f", {head.detail}" if head.detail else "")
    where = _run_describe(summary)
    xs = {}
    if r is not None and head.method == "coulex-yield":
        all_orders = any(st["step"] == "all orders" for st in r.steps)
        method = ("Coulomb excitation, the γ-ray yield of the line "
                  + ("per scattered particle in the same rings" if r.normalisation == "rutherford"
                     else "against a reference transition")
                  + (", solved with the level scheme's coupled equations (all orders)" if all_orders
                     else ", analysed in first order"))
        where += f"; the statistics are those of {r.statistics_from}"
        try:
            xs = {x.key: x for x in analysis_explanations(planner, r)}
        except Exception:  # noqa: BLE001 - the report stands without the meanings
            xs = {}
    sec["result"].append(_headline_html(head, where, method, nuclide_words(nuclide)))
    c = _routes_section(sec, planner, nuclide, "B(E2)", key, own, head, unanalysed, outside=True)
    lead = next((p for p in c.projected if p.id == head.id), None) if c is not None else None
    if lead is not None:
        sec["result"].append(f"<p>{e(closure_words(lead))}</p>")
    down = decay_words(head)
    if down:
        sec["result"].append(f"<p>{e(down)}</p>")
    sec["result"].append(_other_results_html(own, [head]))
    if without:
        sec["result"].append(f"<p>{e(' '.join(without))}</p>")
    # What limits it: the budget, then what the value assumes.
    extra = []
    if r is not None:
        extra += [n[0].upper() + n[1:] + "." for n in (getattr(r, "budget_notes", None) or [])]
        if r.whole_run_statistical is not None:
            extra.append(f"The statistics are those of the run's event-by-event part; the whole run, scaled from "
                         f"it, would give {100 * r.whole_run_statistical:.1f} %.")
    sec["limits"].append(budget_html(head, extra))
    assumes = [xs["b_e2"].assumptions] if "b_e2" in xs else []
    if r is not None and getattr(r, "gate_warning", None):
        assumes.append(str(r.gate_warning))
    if assumes:
        sec["limits"].append("<p>Of the limits of the method, these bear on this value:</p><ul>"
                             + "".join(f"<li>{e(note_words(a))}</li>" for a in assumes) + "</ul>")
    if sec["limits"] and any(sec["limits"]):
        sec["limits"].append("<p>What the simulation leaves out altogether is in appendix H; the analysis step by "
                             "step, with its notes, in appendix F.</p>")
    # What it means: the meanings the explanations of the chain already hold.
    items = []
    if "b_e2" in xs:
        items.append(f"<li><b>{e(quantity_words('B(E2)', key, head.direction))}.</b> {e(xs['b_e2'].meaning)}</li>")
    if down:
        items.append(f"<li><b>In Weisskopf units.</b> {e(down)} {e(WEISSKOPF_MEANING)}</li>")
    for k in ("beta2", "q0", "lifetime"):
        if k in xs:
            x = xs[k]
            items.append(f"<li><b>{e(x.title)}: {e(number(x.value))}{' ' + e(x.unit) if x.unit else ''}.</b> "
                         f"{e(x.meaning)}</li>")
    sh = (r.shape or {}) if r is not None else {}
    if sh.get("e4_over_e2") is not None or sh.get("note"):
        items.append("<li><b>The rotor reading.</b> "
                     + (f"E(4⁺)/E(2⁺) = {sh['e4_over_e2']:.2f}. " if sh.get("e4_over_e2") is not None else "")
                     + e(note_words(sh.get("note") or "")) + "</li>")
    if items:
        sec["meaning"].append("<ul>" + "".join(items) + "</ul>")


def _lifetime_results(sec: dict, planner, runs: list, analyses: dict, plan) -> None:
    """The leading sections of a lifetime experiment: τ of the studied level, and the B(E2) that follows."""
    e = _esc
    exp = planner.experiment
    lifetime = getattr(plan, "lifetime", None) if plan is not None else None
    a = (lifetime or {}).get("analysis") or {}
    # The ledger's rows, and the analysis of this report as rows where the ledger does not hold it (not written).
    own = with_lifetime_rows(planner, a, experiment_rows(planner, runs, analyses))
    taus = [r for r in own if r.observable == "tau" and r.evidence == results.PROJECTED]
    head = _pick(taus, None)
    found = _lifetime_level(lifetime, a.get("level")) or _lifetime_level(lifetime)
    how = "a plunger (recoil distance, RDDS)" if exp.lifetime.method == "rdds" else "line shapes (DSAM)"
    after = "fusion-evaporation" if exp.fusion is not None else "Coulomb excitation"
    beam_runs = [s for s in runs if s["kind"] != "source"]
    if head is None:
        if a.get("available"):  # a τ without a row: the setup holds no τ for the level to close against
            sec["result"].append(
                f'<div style="{_BOX}"><div style="font-size:1.3em">τ({e(a["level"])}) = '
                f'{e(plus_minus(a["tau_ps"], a["stat_ps"]))} (stat) ps <span style="{_TAG}">PROJECTED</span></div>'
                f"<div style=\"margin-top:4px\"><b>Method:</b> {e(how)} after {after}, {e(a['method'])} analysis. "
                f"<b>From:</b> {len(beam_runs)} runs. <b>Projected from simulated runs:</b> the precision the "
                f"setup would reach, not a measurement. In all ± {e(number(a['total_ps'], 2))} ps.</div></div>")
        else:
            why = a.get("reason") or ("No run has been taken with this setup." if not beam_runs else
                                      "The plan with the lifetime analysis could not be built.")
            sec["result"].append(f"<p><b>No lifetime yet.</b> {e(why)}</p>")
        if found is not None:
            _routes_section(sec, planner, found[0], "tau", found[1].text, own)
        return
    nuclide, key = head.nuclide, head.key
    n_runs = len(head.origin.get("runs") or beam_runs)
    distances = len(a.get("points") or [])
    where = f"{n_runs} run{'s' * (n_runs != 1)}" + (f" at {distances} plunger distances" if distances else "")
    method = f"{how} after {after}" + (f", {head.detail}" if head.detail else "")
    sec["result"].append(_headline_html(head, where, method, nuclide_words(nuclide)))
    c = _routes_section(sec, planner, nuclide, "tau", key, own, head)
    lead = next((p for p in c.projected if p.id == head.id), None) if c is not None else None
    if lead is not None:
        sec["result"].append(f"<p>{e(closure_words(lead))}</p>")
    # The B that follows from τ: the experiment's own row where it has one, else what the analysis gives.
    bs = [r for r in own if r.observable == "B(E2)" and r.evidence == results.PROJECTED
          and (head.id in r.depends_on or r.origin.get("analysis") == "lifetime")]
    b_row = _pick(bs, None)
    wu_words = ""
    if b_row is not None:
        wu_words = decay_words(b_row)
        decay = decay_of(compare.canonical(b_row).value, b_row.key)
        sec["result"].append("<p>The B(E2) that follows from τ, with the γ-ray energy, the branching and the "
                             "conversion of the level scheme: "
                             f"{e(quantity_words('B(E2)', b_row.key, b_row.direction))} = {e(row_value(b_row))}"
                             + (f", that is {e(number(decay[1], 3))} W.u. ↓" if decay is not None
                                and b_row.direction == observables.DOWN else "") + "."
                             + (f" {e(wu_words)}" if b_row.direction == observables.UP and wu_words else "")
                             + "</p>")
        sec["routes"].append(f"<h3>{e(quantity_words('B(E2)', b_row.key))}, which follows from τ</h3>")
        _routes_section(sec, planner, nuclide, "B(E2)", b_row.key, own, b_row, to_run=False)
    elif a.get("b_down"):
        lam = a.get("b_multipole", "E2")
        level = found[2] if found is not None else {}
        b, unc = a["b_down"]
        text = (f"B({lam}; {level_words(key)} → {level.get('final') or 'the level below'})↓ = "
                f"{plus_minus(b, unc)} {unit_words(a.get('b_unit') or '')}")
        if lam == "E2":
            try:
                wu_words = f"{text} = {number(observables.weisskopf(b, nuclide), 3)} W.u. ↓."
            except (ValueError, TypeError):
                wu_words = ""
        sec["result"].append(f"<p>The B({e(lam)}) that follows from τ, with the γ-ray energy, the branching and "
                             f"the conversion of the level scheme: {e(text)}.</p>")
    sec["result"].append(_other_results_html(own, [head] + ([b_row] if b_row is not None else [])))
    if a.get("compare") and (c is None or not c.entries):
        sec["routes"].append(
            "<h3>What the setup and ENSDF hold for this level</h3><p>The ledger has no independent value of τ "
            "for this level yet. These are the values the report read before it: the level scheme's is the one "
            "the simulation was given, so it is the closure above and no test.</p><table><thead><tr>"
            f"<th>From</th><th>τ (ps)</th><th>B({e(a.get('b_multipole', 'E2'))}↓)</th></tr></thead><tbody>" + "".join(
                f"<tr><td>{e(x['source'])}</td><td>{e(number(x['tau_ps']))}</td>"
                f"<td>{e(number(x['b_down']))}</td></tr>" for x in a["compare"]) + "</tbody></table>")
    sec["limits"].append(budget_html(head))
    # The limits that bear on τ; those of the plan's own numbers (one typical recoil, analytic spectra) stay
    # with the plan.
    reached = [lim for lim in (lifetime or {}).get("limits", [])
               if lim["flag"] and lim["key"] not in ("typical_recoil", "coulex_spectra")]
    failing = [x for x in (lifetime or {}).get("checks", []) if x["status"] in ("known failure", "fail")]
    if reached or failing:
        sec["limits"].append("<p>Of the limits of the lifetime models, these bear on this value (those of the "
                             "plan's own numbers are with the plan, appendix B):</p><ul>"
                             + "".join(f"<li><b>{e(lim['title'])}.</b> {e(lim['flag'])}</li>" for lim in reached)
                             + "".join(f"<li><b>{e(x['name'])}</b> ({e(x['status'])} against {e(x['checker'])}): "
                                       f"{e(x['detail'])}</li>" for x in failing) + "</ul>")
    if wu_words:
        sec["meaning"].append(f"<ul><li><b>In Weisskopf units.</b> {e(wu_words)} {e(WEISSKOPF_MEANING)}</li></ul>")


def _fusion_results(sec: dict, planner, runs: list, plan) -> None:
    """What a fusion-evaporation experiment gives: the fusion cross section and the residues with their rates.
    It has no result in the ledger yet."""
    e = _esc
    f = getattr(plan, "fusion", None) if plan is not None else None
    if not f:
        sec["result"].append("<p>The plan with the fusion cross section and the residues could not be built: see "
                             "appendix B.</p>")
        return
    s = f.get("summary")
    if s:
        ratio = f", E_cm / V_B = {f['e_cm_mev'] / s['barrier_cm_mev']:.2f}" if f.get("e_cm_mev") else ""
        sec["result"].append(
            f"<p>{e(nuclide_words(s['beam']))} + {e(nuclide_words(s['target']))} fuse to "
            f"{e(nuclide_words(s['compound']))} at {s['excitation_mev']:.1f} MeV of excitation{e(ratio)}: "
            f"σ<sub>fus</sub> = {e(number(s['sigma_thick_mb']))} mb averaged over the target, with "
            f"⟨ℓ⟩ = {s['mean_l']:.1f} ħ.</p>")
    rows = list((f.get("residues") or {}).get("rows") or [])
    if rows:
        studied = None
        try:
            studied = planner.experiment.fusion.channels[planner.studied_channel()].particles
        except (TypeError, IndexError, AttributeError):
            pass
        top = sorted(rows, key=lambda x: -(x.get("rate_per_s") or 0.0))[:6]
        sec["result"].append(
            "<table><thead><tr><th>Residue</th><th>Channel</th><th>σ (mb)</th><th>Cross section from</th>"
            "<th>Made per second</th><th>In the planned run</th></tr></thead><tbody>" + "".join(
                f"<tr><td>{e(nuclide_words(x['residue']))}" + (" (studied)" if x["channel"] == studied else "")
                + f"</td><td>{e(x['channel'])}</td><td>{e(number(x['cross_section_mb']))}</td>"
                f"<td>{e(x['source'])}</td><td>{e(number(x['rate_per_s']))}</td>"
                f"<td>{e(number(x['residues_in_run']))}</td></tr>" for x in top) + "</tbody></table>"
            + (f"<p>The {len(top)} strongest of {len(rows)} channels; all are in appendix B.</p>"
               if len(rows) > len(top) else ""))
    _runs_against_plan(sec, planner, runs)
    sec["result"].append("<p>These are the plan's rates, not a measured quantity: a fusion-evaporation experiment "
                         "has no result in the results ledger yet, so there is nothing to set beside other "
                         "routes.</p>")
    reached = [lim for lim in f.get("limits", []) if lim["flag"]]
    failing = [x for x in f.get("checks", []) if x["status"] in ("known failure", "fail")]
    if reached:
        sec["limits"].append("<p>Of the limits of the fusion-evaporation models, this run reaches these:</p><ul>"
                             + "".join(f"<li><b>{e(lim['title'])}.</b> {e(lim['flag'])}</li>" for lim in reached)
                             + "</ul>")
    if failing:
        who = sorted({x["checker"] for x in failing})
        own = rows and not any(str(x.get("source", "")).startswith("physim") for x in rows)
        sec["limits"].append(
            f"<p>{len(failing)} check{'s' * (len(failing) != 1)} of the models against {e(', '.join(who))} "
            f"{'are' if len(failing) != 1 else 'is a'} known failure{'s' * (len(failing) != 1)}: disagreements "
            "the register documents, each given with its numbers beside the result it bears on in appendix B."
            + (" The channel cross sections here are the setup's own, so the evaporation checks bear on them only "
               "if the table is filled from physim's evaporation." if own else "") + "</p>")


def _elastic_results(sec: dict, planner, runs: list, plan, pl) -> None:
    """What an elastic-scattering experiment gives: the rate in each detector and the beam time."""
    e = _esc
    if pl:
        sec["result"].append(f"<p>Beam time needed: {_time(pl['beam_time_needed_s'])}"
                             + (f" ({e(str(pl['limiting_detector']))}, for {pl['counts_wanted']} "
                                f"{e(_counted(pl))})" if pl["limiting_detector"] else "")
                             + f"; planned: {_time(pl['beam_time_s'])}.</p>")
    rows = list(getattr(plan, "detectors", None) or []) if plan is not None else []
    if rows:
        sec["result"].append(
            "<table><thead><tr><th>Detector</th><th>θ (deg)</th><th>Rate (1/s)</th>"
            "<th>Counts in the planned run</th><th>Statistical precision</th></tr></thead><tbody>" + "".join(
                f"<tr><td>{e(x['detector'])}</td><td>{x['theta_min_deg']:.1f}–{x['theta_max_deg']:.1f}</td>"
                f"<td>{e(number(x['rate_per_s']))}</td><td>{e(number(x['counts_in_run']))}</td>"
                f"<td>{e(percent(100 * x['relative_error']) if x.get('relative_error') else '—')}</td></tr>"
                for x in rows) + "</tbody></table>")
    elif pl:
        sec["result"].append("<table><thead><tr><th>Detector</th><th>Rate (1/s)</th></tr></thead><tbody>" + "".join(
            f"<tr><td>{e(d['detector'])}</td><td>{e(number(d['rate_per_s']))}</td></tr>" for d in pl["detectors"])
            + "</tbody></table>")
    _runs_against_plan(sec, planner, runs)
    sec["result"].append("<p>These are rates and beam times, not a measured quantity of a nucleus: an "
                         "elastic-scattering experiment has no result in the results ledger, so there is nothing "
                         "to set beside other routes.</p>")


def _runs_against_plan(sec: dict, planner, runs: list) -> None:
    """A line on the runs taken, and the current run against the plan, for an experiment without a result."""
    beam_runs = [s for s in runs if s["kind"] != "source"]
    if not beam_runs:
        sec["result"].append("<p>No run has been taken with this setup: these are the plan's numbers.</p>")
        return
    c = {"rows": []}
    if planner.run is not None and planner.run.kind != "source":
        try:
            c = planner.compare()
        except Exception:  # noqa: BLE001 - the report stands without the comparison
            c = {"rows": []}
    judged = [r for r in c["rows"] if r["agree"] is not None]
    text = f"{len(beam_runs)} run{'s' * (len(beam_runs) != 1)} taken."
    if judged:
        differ = [f"{r['quantity'].lower()} ({r['what']})" for r in judged if not r["agree"]]
        text += (f" Run {c['run']} against the plan: " + (
            f"{len(judged)} of {len(judged)} numbers agree" if not differ else
            f"{len(judged) - len(differ)} of {len(judged)} numbers agree; these differ: {', '.join(differ)}")
            + " (appendix C has the table).")
    sec["result"].append(f"<p>{_esc(text)}</p>")


def experiment_html(planner, result=None, scene_png: Optional[str] = None, plan=None, plan_events: int = 50_000,
                    seed: int = 1, journal: str = "physical_review", width: str = "single",
                    validation: bool = False) -> str:
    """The report of the experiment, one self-contained page that leads with the results.

    First the results (:data:`RESULTS`): the result with its statistical and systematic uncertainties, its
    method and its runs, and plainly that it is projected when it comes from a simulated run; the result beside
    the other routes to the same quantity (the independent evidence against the reference, and, apart, the
    projected rows with their closure); the comparisons still to run; what limits the result; and what it
    means. The rows come from the experiment's results ledger (:func:`experiment_rows`). An elastic or
    fusion-evaporation experiment has no result in the ledger: it opens with what it gives (rates, residues).

    Then the appendices (:data:`APPENDICES`), lettered: the setup with its warnings, level schemes and how it
    changed between runs; the plan (what the beam-time report holds: detectors, rates and beam time, kinematics,
    expected peaks, Coulomb excitation, energy loss, predicted spectra); every run with its summary and its
    spectra (the particle detectors, the energy against ring with the gates drawn on it, every γ-ray detector
    as measured and Doppler-corrected to the top of the level scheme); the gates and what each does to the
    γ-ray spectrum; the Doppler correction and its parameters; the analysis of every analysed run, the alignment
    and the multi-step solution; every number explained; and the data, the method and what the simulation leaves
    out.

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
    what = reaction_words(planner)
    beam_runs = [s for s in runs if s["kind"] != "source"]
    keys = {s["number"]: run_key(planner, s) for s in runs if s["number"] in analyses}  # each run's own setup
    try:
        pl = planner.plan()
    except Exception:  # noqa: BLE001 - the report stands without the plan
        pl = None

    # -- the results: what the experiment found, first --------------------------------------------------------
    titles = dict(SECTIONS)
    try:
        if exp.lifetime is not None:
            _lifetime_results(sec, planner, runs, analyses, plan)
        elif exp.excitation is not None:
            _coulex_results(sec, planner, runs, analyses, result, pl)
        elif exp.fusion is not None:
            titles["result"] = "What this setup gives"
            _fusion_results(sec, planner, runs, plan)
        else:
            titles["result"] = "What this setup gives"
            _elastic_results(sec, planner, runs, plan, pl)
    except Exception as err:  # noqa: BLE001 - the appendices stand without the leading sections
        sec["result"].append(f"<p>The results could not be put together: {e(err)}. The appendices hold the "
                             "analysis of every run.</p>")

    # -- A. the setup -----------------------------------------------------------------------------------------
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

    # -- B. the plan ------------------------------------------------------------------------------------------
    sec["plan"].append("<p>What the current setup should give, computed from the setup and not from the runs: "
                       "the runs are compared with it in the next appendix.</p>")
    if pl:
        sec["plan"].append("<h3>Beam time and load</h3>"
                           f"<p>Beam time needed: {_time(pl['beam_time_needed_s'])}"
                           + (f" ({e(str(pl['limiting_detector']))}, for {pl['counts_wanted']} {e(_counted(pl))})"
                              if pl["limiting_detector"] else "") + f"; planned: {_time(pl['beam_time_s'])}.</p>")
        sec["plan"].append("<table><thead><tr><th>Detector</th><th>Rate (1/s)</th><th>Busy</th>"
                           "<th>Coincidences (1/s)</th><th>Safe rings</th></tr></thead><tbody>"
                           + "".join(f"<tr><td>{e(d['detector'])}</td><td>{number(d['rate_per_s'])}</td>"
                                     f"<td>{100 * d['dead_time_fraction']:.2g}%</td>"
                                     f"<td>{number(d['coincidence_per_s'])}</td>"
                                     f"<td>{len(d['safe_rings'])} of {d['rings']}</td></tr>" for d in pl["detectors"])
                           + "</tbody></table>")
        if pl.get("other_peak"):
            op = pl["other_peak"]
            sec["plan"].append(
                f"<p>{e(op['nuclide'])} is excited too (its level scheme is in the setup). Its {e(op['label'])} peak "
                f"at {op['energy_kev']:.1f} keV, in coincidence: "
                + "; ".join(f"{e(d['detector'])} {number(d['other_coincidence_per_s'])} /s" for d in pl["detectors"]
                            if d.get("other_coincidence_per_s") is not None)
                + ". The counts wanted and the beam time are for the nucleus of the reaction.</p>")
    if pp:
        sec["plan"] += ["<h3>Detectors: coverage, cross sections, rates and beam time</h3>", pp["detectors"],
                        "<h3>Kinematics</h3>", pp["kinematics"], "<h3>Expected peaks</h3>", pp["peaks"],
                        pp["excitation_plain"]] + ([pp["fusion_plain"]] if pp.get("fusion_plain") else []) + (
                        [pp["lifetime_plain"]] if pp.get("lifetime_plain") else []) + [
                        "<h3>Energy loss</h3>", pp["energy_loss"],
                        "<h3>Simulated spectra</h3>", pp["spectra"]]
    elif why_no_plan:
        sec["plan"].append(f"<p>The plan's tables and figures could not be built: {e(why_no_plan)}.</p>")

    # -- C. the runs ------------------------------------------------------------------------------------------
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
                               + e(", ".join(f"{d} {c:,} ({number(s['rates_per_s'].get(d, 0))}/s)"
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
            sec["runs"].append(f"<h3>Run {c['run']}: predicted against measured</h3>" + predicted_html(c))

    # -- D. the gates -----------------------------------------------------------------------------------------
    if gates:
        sec["gates"].append("<p>A gate is a condition on the particles; a γ-ray spectrum under a gate holds the "
                            "γ rays in coincidence with the particles that pass it.</p><ul>"
                            + "".join(f"<li><b>{e(g.name)}</b>: {e(g.describe())}"
                                      + (f" ({e(g.note)})" if g.note else "") + ".</li>" for g in gates) + "</ul>")
        sec["gates"] += [pieces[s["number"]]["gates"] for s in beam_runs]
    else:
        sec["gates"].append("<p>No gates: the spectra take every counted particle. A gate is made on the Data "
                            "tab, from the energy-against-ring view or by naming detectors and rings.</p>")

    # -- E. the Doppler correction ----------------------------------------------------------------------------
    try:
        dop = doppler_html(planner, analyses, {n: x["seen"] for n, x in pieces.items() if x["seen"] is not None})
    except Exception as err:  # noqa: BLE001 - the report stands without the section
        dop = f"<p>The Doppler correction could not be described: {e(err)}.</p>"
    sec["doppler"].append(dop or "<p>No excited state or no γ-ray detector: there is no γ ray to correct.</p>")

    # -- F. the analysis --------------------------------------------------------------------------------------
    for s in runs:
        r = analyses.get(s["number"])
        if r is not None:
            sec["analysis"].append(analysis_html(r, f"<h3>Run {s['number']}: from the peak area to B(E2)</h3>",
                                                 keys.get(s["number"])))
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
    lifetime_done = bool(((getattr(plan, "lifetime", None) or {}).get("analysis") or {}).get("available"))
    if not analyses and result is not None:  # the analysis of a simulated sample, without a run
        sec["analysis"].append(analysis_html(result, "<h3>The simulated sample</h3>", run_key(planner, None)))
    elif not analyses and lifetime_done:
        sec["analysis"].append("<p>The lifetime analysis of the runs (the decay curve or the line shapes, its "
                               "budget and its notes) is in appendix B, under “Lifetime measurement”, with the "
                               "plan it belongs to.</p>")
    elif not analyses and exp.excitation is None:
        sec["analysis"].append("<p>No run has been analysed: " + ("a fusion-evaporation" if exp.fusion is not None
                                                                    else "an elastic-scattering")
                               + " run has no analysis step. What it gives is the comparison of the run with the "
                               "plan, in appendix C.</p>")
    elif not analyses:
        sec["analysis"].append("<p>No run has been analysed yet: take a beam run and press Analyse on the "
                               "Analysis step. The analysis is then kept with the run and appears here.</p>")
    m = planner.kept_multistep()
    if m and m.get("available"):
        sec["analysis"].append("<h3>All orders: multi-step excitation and reorientation</h3>" + multistep_html(m))

    # -- G. every number, explained ---------------------------------------------------------------------------
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
            key_now = keys.get(planner.run.number) if planner.run is not None else run_key(planner, None)
            for x in ax:
                shown = weisskopf_shown(x, key_now) if x.key == "weisskopf" else x
                sec["numbers"].append(
                    _explanation_html(shown) if shown is not None else
                    "<p>B(E2) is not given in Weisskopf units here: they are quoted for the decay, B(E2↓), which "
                    "needs the spins of the two levels, and those are not known for this transition.</p>")

    # -- H. the data, the method and what is left out ---------------------------------------------------------
    data_rows = []
    for q, v, src in provenance(exp):
        if exp.excitation is not None and src == "user" and str(q).startswith("State at "):
            src = reaction_b_source(exp)  # the [reaction]'s B: "user" only where the setup holds its own value
        data_rows.append(f"<tr><td>{e(q)}</td><td>{e(str(v))}</td><td>{e(wu_marked(src))}</td></tr>")
    sec["method"].append("<h3>The nuclear data and where it comes from</h3>"
                         "<table><thead><tr><th>Quantity</th><th>Value</th><th>Source</th></tr></thead><tbody>"
                         + "".join(data_rows) + "</tbody></table>"
                         "<p>Sources: <em>ensdf</em> is the evaluated data as read; <em>derived</em> is computed "
                         "from ENSDF's values by physim; <em>assumed</em> is a typical value, not a measurement; "
                         "<em>user</em> was given in the setup. A B(E2) in Weisskopf units is the decay's (↓), as "
                         "ENSDF quotes it; the matrix elements and the state's B(E2↑) are the excitation's.</p>")
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

    # -- the page: the head, the results numbered, the appendices lettered ------------------------------------
    for key in sec:
        sec[key] = [x for x in sec[key] if x]
    shown = [key for key, _ in RESULTS if sec[key]]
    analysed = len(analyses)
    counts = (f"{len(runs)} run{'s' * (len(runs) != 1)}: {len(beam_runs)} with beam, {len(runs) - len(beam_runs)} "
              "with a source"
              + (f", {analysed} analysed" if exp.excitation is not None and exp.lifetime is None and beam_runs
                 else "") if runs else "no run yet: the report holds the setup and the plan")
    body = [f"<h1>{e(planner.name)}</h1>",
            f'<p class="sub">{e(b.nuclide)} at {e(str(b.energy))} on {e(t.material)} ({e(str(t.thickness))}): '
            f"{e(what)} · {counts}" + (f" · {pp['meta']}" if pp else "") + "</p>",
            '<ol class="toc">' + "".join(f'<li><a href="#{key}">{e(titles[key])}</a></li>' for key in shown)
            + '</ol><p style="margin-bottom:0"><a href="#appendices">Appendices</a>: the report as it was, kept '
            'whole.</p><ol class="toc" style="list-style:upper-alpha;margin-top:0.2em">'
            + "".join(f'<li><a href="#{key}">{e(title)}</a></li>' for key, title in APPENDICES) + "</ol>"]
    for i, key in enumerate(shown, 1):
        body.append(f'<h2 id="{key}">{i}. {e(titles[key])}</h2>')
        body += sec[key]
    body.append('<h2 id="appendices" style="font-size:1.45em;border-top:2px solid #999;padding-top:0.6em;'
                'margin-top:2em;break-before:page">Appendices</h2><p>Everything behind the results: the setup, the '
                "plan, every run with its spectra, the gates, the Doppler correction, the analysis step by step, "
                "every number explained, and the method with its validation.</p>")
    for i, (key, title) in enumerate(APPENDICES):
        body.append(f'<h2 id="{key}">{chr(65 + i)}. {e(title)}</h2>')
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


# -- the report of a nucleus, across its experiments ------------------------------------------------------------------

#: The appendices of the report of a nucleus, after its results: (anchor, title), lettered.
NUCLEUS_APPENDICES = (("experiments", "The experiments"), ("history", "The history of the rows"))

#: A route's state in the report of a nucleus.
_ROUTE_STATE = {nucleus.DONE: "done", nucleus.ONLY_PROJECTED: "simulated runs only", nucleus.IN_PHYSIM:
                "physim can plan it", nucleus.LATER: "a kind physim does not have yet",
                nucleus.NO_EXPERIMENT: "no experiment"}
_APPLIES = {True: "yes", False: "no", None: "not tested"}
_MEANS = {"calculation": "a calculation", "evaluation": "the evaluation", "literature": "the literature"}


def _units_written(text) -> str:
    """Text with the ledger's units in it ("884 ± 13 e2fm4↑") as the report writes them ("e²fm⁴")."""
    return re.sub(r"(?<![\w])(e2fm\d|e2b\d?|efm\d?|mu_N2?)(?![\w])", lambda m: unit_words(m.group(1)), str(text))


def _nucleus_quantity_html(name: str, t: dict, ov: dict, by_id: dict) -> list:
    """The parts of the report of a nucleus for one observable of one level or transition: the values side by
    side in their two groups, the routes, what is still to run, the checks and what the value gives."""
    def e(text) -> str:
        return _esc(_units_written(text))

    out = []
    c = ov.get("comparison")
    counts = [f"{n} {k}" for k, n in (("evaluated", t["evaluated"] is not None and 1 or 0),
                                     ("measured", (t["measured"] or {}).get("count", 0)),
                                     ("from an outside program", (t["outside_program"] or {}).get("count", 0)),
                                     ("calculated", (t["theory"] or {}).get("count", 0))) if n]
    projected = (t["projected"] or {}).get("count", 0)
    out.append(f"<p>{e(', '.join(counts) or 'No independent value')}"
               + (f"; {projected} projected from simulated runs, kept apart" if projected else "") + ".</p>")
    out.append("<h3>The values, side by side</h3>")
    if c is not None:
        out.append(comparison_html(c, rule=False, by_id=by_id))
    else:
        out.append("<p>No row holds for it now.</p>")
    # The routes, each with its state for this nucleus and whether it applies.
    rows = []
    for r in ov.get("routes") or []:
        state = _ROUTE_STATE.get(r["status"], r["status"])
        if r["status"] == nucleus.DONE:
            state += f": {r['independent']} value{'s' * (r['independent'] != 1)}"
        if r["projected"]:
            state += f" ({r['projected']} simulated run{'s' * (r['projected'] != 1)})" if r["status"] == \
                nucleus.DONE else f": {r['projected']}"
        if r.get("is_reference"):
            state += "; holds the reference"
        kind = r["kind_label"] or _MEANS.get(r["means"], "")
        rows.append(f"<tr><td>{e(r['label'])}</td><td>{e(kind)}</td><td>{e(state)}</td>"
                    f"<td>{e(_APPLIES[r['applies']])}</td><td>{e(r['reason'] or r['note'])}</td></tr>")
    if rows:
        out.append("<h3>The routes to it</h3><table><thead><tr><th>Route</th><th>Experiment</th><th>State</th>"
                   "<th>Applies here</th><th>Why, or how</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>")
    # What is still to run: the routes physim can plan that nothing of one's own has taken, and the independent
    # evidence that is missing.
    items = []
    for r in ov.get("routes") or []:
        if r["open"] and r["applies"] is not False:
            items.append(f"<li>☐ <b>{e(r['label'])}.</b> physim has this kind of experiment ({e(r['kind_label'])})"
                         + (" and a published value by it is here, but no run of one's own" if r["rows"] else
                            " and no value by it is here") + ".</li>")
    group = [r for r in results.current(list(by_id.values())) if r.observable == ov["observable"]
             and r.key == ov["key"]]
    for x in evidence_missing(group, routes=ov.get("routes") or []):
        if not x["have"]:
            items.append(f"<li>☐ <b>{e(x['label'])}.</b> {e(x['missing'])}</li>")
    out.append("<h3>Still to run</h3>" + (f'<ul style="{_PLAIN_LIST}">' + "".join(items) + "</ul>" if items
                                         else "<p>Nothing: every route physim can plan has a run of one's own, and "
                                              "every kind of independent evidence is entered.</p>"))
    if ov.get("checks"):
        out.append("<h3>Checks between observables</h3><ul>" + "".join(
            f"<li><b>{e(x['label'])}:</b> {e(x['text'])}<br><small>{e(x['how'])}"
            + (f" ({e(x['reference'])})" if x.get("reference") else "") + "</small></li>" for x in ov["checks"])
                   + "</ul>")
    ref = ov.get("reference")
    gives = [g for g in ov.get("gives") or []]
    if ref is not None and gives:
        out.append(f"<h3>What the value gives</h3><p>From the reference, {e(ref['text'])} ({e(ref['source'])}):</p>"
                   "<ul>" + "".join(
                       f"<li>{e(g['to'])} = {e(g['text'])}<br><small>{e(g['how'])}</small></li>" if g["value"]
                       is not None else f"<li>{e(g['to'])}: not worked out ({e(g['missing'])})</li>" for g in gives)
                   + "</ul>")
    if ov.get("used_for"):
        out.append("<h3>What it is used for</h3><ul>" + "".join(
            f"<li>{e(u['text'])}" + (f"<br><small>{e(u['reference'])}</small>" if u.get("reference") else "")
            + "</li>" for u in ov["used_for"]) + "</ul>")
    return out


def nucleus_html(nuclide: str, root=None, home=None, rows: Optional[list] = None) -> str:
    """The report of a nucleus across all its experiments, one self-contained HTML page that leads with the
    results, in the style of the experiment's report (:func:`experiment_html`).

    First, for each observable and level or transition that has more than the evaluation (the headline rows of
    :func:`physim.nuclear.nucleus.view`), a numbered section from :func:`physim.nuclear.nucleus.observable_view`:
    the values side by side in the two groups that are never mixed (the independent evidence against the
    reference with difference, pull and verdict; apart, the simulated runs, each closed against the value its
    simulation was given and never pulled against the reference: the table the nucleus page and the experiment's
    report show, :func:`comparison_html`), the routes with their state and whether they apply to this nucleus,
    what is still to run, the checks between observables, and what the value gives and is used for. Then, as
    appendices, the experiments (name, kind, runs, the latest result) and the history of the rows, every row ever
    written, the newest first.

    ``root`` is the experiments folder and ``home`` the nuclei home (the defaults of
    :func:`physim.nuclear.results.of_nuclide`); ``rows`` the nucleus's rows where the caller has them. Building it
    writes nothing: it neither takes the evaluation nor backfills an experiment."""
    def e(text) -> str:
        return _esc(_units_written(text))

    name = results.nuclide_name(nuclide)
    every = results.of_nuclide(name, root, home) if rows is None else [r for r in rows if r.nuclide == name]
    v = nucleus.view(name, root, home, rows=every)
    by_id = {r.id: r for r in every}
    title = nuclide_words(name)
    quantities = [t for t in v["table"] if t["headline"] and t["rows"]]
    sections = []
    for i, t in enumerate(quantities, 1):
        head = f'<h2 id="q{i}">{i}. {e(quantity_words(t["observable"], t["key"]))}</h2>'
        try:
            ov = nucleus.observable_view(name, t["observable"], t["key"], rows=every)
            sections.append((t, head, _nucleus_quantity_html(name, t, ov, by_id)))
        except (ValueError, TypeError, KeyError) as err:
            sections.append((t, head, [f"<p>Its values could not be set beside each other: {e(err)}.</p>"]))
    studied = [x for x in v["experiments"] if x["role"] == "studied"]
    ev = v["evaluation"]
    body = [f"<h1>{e(title)}: the report of the nucleus</h1>",
            f'<p class="sub">Z = {v["z"]}, N = {v["n"]} · {len(studied)} experiment{"s" * (len(studied) != 1)} '
            f"that stud{'y' if len(studied) != 1 else 'ies'} it, {v['analysed_runs']} analysed "
            f"run{'s' * (v['analysed_runs'] != 1)}, {v['points']} data point{'s' * (v['points'] != 1)} · "
            + (f"the evaluation: {e(ev['reference'])} ({ev['rows']} values)" if ev["rows"]
               else "the evaluation is not taken") + "</p>",
            '<ol class="toc">' + "".join(f'<li><a href="#q{i}">{e(quantity_words(t["observable"], t["key"]))}</a>'
                                         "</li>" for i, (t, _, _) in enumerate(sections, 1)) + "</ol>"
            + '<p style="margin-bottom:0"><a href="#appendices">Appendices</a></p><ol class="toc" '
              'style="list-style:upper-alpha;margin-top:0.2em">'
            + "".join(f'<li><a href="#{k}">{e(x)}</a></li>' for k, x in NUCLEUS_APPENDICES) + "</ol>"]
    if sections:
        body.append(f"<p>Every value of {e(title)} that is entered, across its experiments and its own data "
                    "points, one section per observable and level or transition that has more than the "
                    "evaluation. Each value is shown in the observable's convention and, where it was given in "
                    f"another, as it was given.</p><p><small>{e(RULE)}</small></p>")
    else:
        body.append(f"<p>Nothing beyond the evaluation is entered for {e(title)} yet: no data point, no analysed run "
                    "and no calculation.</p>")
    if v["open"]:
        body.append("<p><b>Still open for the nucleus:</b> " + e("; ".join(o["text"] for o in v["open"])) + ".</p>")
    for _, head, parts in sections:
        body.append(head)
        body += parts
    body.append('<h2 id="appendices" style="font-size:1.45em;border-top:2px solid #999;padding-top:0.6em;'
                'margin-top:2em;break-before:page">Appendices</h2>')
    # A. the experiments
    rows_ = []
    for x in v["experiments"]:
        latest = x["latest"]
        said = ("—" if latest is None else
                f"{quantity_words(latest['observable'], latest['key'])} = {latest['text']}, "
                f"{results.EVIDENCE[latest['evidence']].lower()} ({latest['method_label']})")
        role = {"studied": "studies it", "partner": "the partner"}.get(x["role"], "results only")
        rows_.append(f"<tr><td>{e(x['name'])}<br><small>{e(x['folder'])}</small></td><td>{e(x['kind_label'] or '')}"
                     f"</td><td>{e(role)}</td><td>{x['runs']} ({x['analysed_runs']} analysed)</td>"
                     f"<td>{e(said)}</td></tr>")
    body.append(f'<h2 id="experiments">A. {e(NUCLEUS_APPENDICES[0][1])}</h2>'
                + ("<table><thead><tr><th>Experiment</th><th>Kind</th><th>Role</th><th>Runs</th><th>Latest result"
                   "</th></tr></thead><tbody>" + "".join(rows_) + "</tbody></table>" if rows_ else
                   f"<p>No experiment with {e(title)} yet.</p>"))
    # B. the history of the rows
    status = {"holds": "holds", "superseded": "superseded by a newer row", "retracted": "withdrawn",
              "retraction": "the retraction"}
    past = nucleus.history(every)
    hist = []
    for h in past:
        hist.append(f"<tr><td>{e(h['created'][:10])}</td><td>{e(quantity_words(h['observable'], h['key']))}</td>"
                    f"<td>{e(results.EVIDENCE[h['evidence']])} · {e(h['method_label'])}"
                    + (f"<br><small>{e(h['detail'])}</small>" if h["detail"] else "")
                    + f"</td><td>{e(h['source'])}</td><td>{e('' if h['status'] == 'retraction' else h['text'])}</td>"
                    f"<td>{e(status.get(h['status'], h['status']))}"
                    + (f": {e(h['reason'])}" if h["reason"] else "") + "</td></tr>")
    body.append(f'<h2 id="history">B. {e(NUCLEUS_APPENDICES[1][1])}</h2><p>Every row ever written for '
                f"{e(title)}, the newest first: rows are never rewritten; a correction or a retraction is a new "
                "row.</p>" + ("<table><thead><tr><th>Written</th><th>Of</th><th>Evidence · method</th><th>From</th>"
                              "<th>Value as given</th><th>Status</th></tr></thead><tbody>" + "".join(hist)
                              + "</tbody></table>" if hist else "<p>No row yet.</p>"))
    return ("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>"
            f"{e(name)}: the report of the nucleus</title><style>{record_css()}</style></head><body>"
            + "\n".join(body) + "</body></html>")


def nucleus_pdf(nuclide: str, page: Optional[str] = None, **options) -> Optional[bytes]:
    """The report of a nucleus as a PDF: :func:`nucleus_html` (or ``page``, if already made) printed by a browser
    of this computer (:func:`html_to_pdf`). None when no such browser is found."""
    return html_to_pdf(page if page is not None else nucleus_html(nuclide, **options))


# The names these had before the app called them.
_weisskopf_shown, _closure_words, _known_words = weisskopf_shown, closure_words, known_words

__all__ = ["APPENDICES", "GAMMA_BIN_KEV", "NUCLEUS_APPENDICES", "closure_words", "known_words", "nucleus_html",
           "nucleus_pdf", "weisskopf_shown", "NOT_AVAILABLE", "RESULTS", "RULE", "SECTIONS", "alignment_html",
           "analysis_html", "analysis_png", "analysis_row", "asymmetric", "budget_html", "comparison_html",
           "decay_of", "decay_words", "doppler_html", "evidence_missing", "experiment_html", "experiment_pdf",
           "experiment_rows", "experiment_zip", "gamma_lines", "html_to_pdf", "level_scheme_png", "level_schemes",
           "level_schemes_html", "level_words", "line_png", "matched", "multistep_html", "note_words",
           "nucleus_rows", "number", "observable_rows", "pdf_browser", "percent", "plain_numbers", "plus_minus", "predicted_html",
           "quantity_words", "reaction_b_source", "reaction_words", "row_value", "run_analyses",
           "run_explanations", "run_gamma_png", "run_gates_png", "run_key", "run_rings_png", "run_spectra_png",
           "run_views", "search_png", "setup_history", "source_words", "to_run_html", "unit_words",
           "with_lifetime_rows", "wu_marked"]
