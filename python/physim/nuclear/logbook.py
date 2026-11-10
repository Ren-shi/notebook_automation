"""The record of an experiment and its runs (backlog item 69): the report of the experiment, and the explanations of
the numbers the runs add (their counters, the gates, predicted against measured).

The beam-time report (:mod:`physim.nuclear.report`) and the run record (:mod:`physim.nuclear.record`) describe one
setup; this module describes the experiment: its current setup and how it changed from run to run, the Plan, every
run with its summary and the setup it was taken with, the gates of the Data tab, the analysis of the run the user
chose, and the links to the physics register::

    from physim.nuclear import logbook

    page = logbook.experiment_html(planner, result)          # one self-contained page, every run's spectra in it
    data = logbook.experiment_zip(planner, result)           # it, the run record, the setup, the beam-time report,
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


def run_views(planner, run, result=None) -> dict:
    """What the record draws of a run: the spectra of every particle detector and, for every γ-ray detector, the
    raw spectrum and the Doppler-corrected one (:mod:`physim.nuclear.dataviews`), with the correction and the gate
    the analysis used (``result``) or, without an analysis, the correction for the nucleus that was excited and no
    gate. Returns {"particles": {name: spectrum}, "gammas": {name: {"raw": spectrum, "corrected": spectrum}},
    "correction", "gate", "source": the crystals' spectra of a source run}."""
    from . import dataviews as dv

    if run.kind == "source":
        sr = run.source()
        return {"particles": {}, "gammas": {}, "correction": None, "gate": None,
                "source": {n: {"counts": sr.spectra[n], "edges": sr.edges} for n in sr.names()}}
    exp = run.experiment
    settings = dict(getattr(result, "settings", {}) or {})
    correction = settings.get("correction", "emitter")
    if correction not in ("projectile", "recoil"):
        correction = "recoil" if exp.excitation is not None and exp.excitation.excite == "target" else "projectile"
    gate = _gate_of(planner, settings)
    ev = run.events()
    out = {"particles": {d: dv.particle_spectrum(run, d, gate=gate) for d in ev.detectors}, "gammas": {},
           "correction": correction, "gate": gate, "source": None}
    g = run.gammas()
    if g is not None:
        for n in g.detector_names():
            out["gammas"][n] = {"raw": dv.gamma_spectrum(run, n, correction="off", gate=gate, randoms="none"),
                                "corrected": dv.gamma_spectrum(run, n, correction=correction, gate=gate,
                                                               randoms="none")}
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


def run_spectra_png(planner, run, result=None, dpi: int = 110) -> tuple:
    """The run's spectra as one PNG (bytes) for the record, and its caption: a panel per particle detector, and
    one per γ-ray detector with the raw and the Doppler-corrected spectrum over each other (a source run: a panel
    per crystal). Needs matplotlib."""
    from matplotlib.figure import Figure

    v = run_views(planner, run, result)
    panels = []
    if v["source"]:
        for n, s in v["source"].items():
            panels.append((n, [("counts", s["edges"], s["counts"])], "γ-ray energy (keV)", True))
    else:
        for n, s in v["particles"].items():
            panels.append((n, [("measured", s["edges"] * 1e3, s["counts"])], "energy (keV)", False))
        for n, s in v["gammas"].items():
            panels.append((n, [("as measured", s["raw"]["edges"] * 1e3, s["raw"]["counts"]),
                               (f"corrected for the {v['correction']}", s["corrected"]["edges"] * 1e3,
                                s["corrected"]["counts"])], "γ-ray energy (keV)", False))
    if not panels:
        return b"", ""
    cols = min(3, len(panels))
    rows = -(-len(panels) // cols)
    fig = Figure(figsize=(3.4 * cols, 2.5 * rows), dpi=dpi)
    axes = fig.subplots(rows, cols, squeeze=False)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    for ax, (title, lines, xlabel, log) in zip(axes.flat, panels):
        for i, (label, edges, counts) in enumerate(lines):
            y = np.asarray(counts, dtype=float)
            ax.stairs(np.where(y > 0, y, np.nan) if log else y, np.asarray(edges, dtype=float), label=label,
                      color=("#0072B2", "#D55E00")[i % 2], lw=0.9 if i == 0 else 1.1)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("counts", fontsize=8)
        ax.tick_params(labelsize=7)
        if log:
            ax.set_yscale("log")
        if len(lines) > 1:
            ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    if v["source"]:
        caption = f"Run {run.number}: the counts each crystal recorded from the source, over {_time(run.duration_s)}."
    else:
        caption = (f"Run {run.number}: the real part of the run ({_time(run.real_s)}), counts per bin. γ-ray "
                   f"spectra as measured and Doppler-corrected for the {v['correction']}, in coincidence with "
                   + (f"the particles of gate “{v['gate'].name}” ({v['gate'].describe()})" if v["gate"]
                      else "every counted particle") + "; random coincidences not added; add-back and "
                   "suppression as the setup has them.")
    return buf.getvalue(), caption


def _run_figure_html(planner, s: dict, result) -> str:
    """The run's spectra as an embedded figure, or a line saying why not."""
    try:
        run = planner.run if planner.run is not None and planner.run.number == s["number"] \
            else planner.folder.load_run(s["number"])
        png, caption = run_spectra_png(planner, run, result)
    except ImportError:
        return "<p>The spectra are not drawn: matplotlib is not installed.</p>"
    except Exception as err:  # noqa: BLE001 - the record stands without the figure
        return f"<p>The spectra could not be drawn: {html.escape(str(err), quote=False)}.</p>"
    if not png:
        return ""
    return (f'<figure><img alt="Spectra of run {s["number"]}" style="max-width:100%" '
            f'src="data:image/png;base64,{base64.b64encode(png).decode("ascii")}"><figcaption>'
            f"{html.escape(caption, quote=False)}</figcaption></figure>")


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
                       "nucleus (Physics tab), so no other level is drawn.")
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


def experiment_html(planner, result=None, scene_png: Optional[str] = None) -> str:
    """The report of the experiment as one self-contained page: the setup and how it changed between runs, the
    Plan, every run with its summary and its spectra (every particle detector; every γ-ray detector as measured
    and Doppler-corrected, :func:`run_spectra_png`), the gates, the analysis of the current run (``result``, a
    :class:`~physim.nuclear.analysis.Result`), the run's numbers explained, and the physics register."""
    from .planner import REGISTER

    def e(text) -> str:  # text, not attributes: quotes stay as they are
        return html.escape(str(text), quote=False)

    exp = planner.experiment
    parts = [f"<h1>{e(planner.name)}: the experiment</h1>"]
    if planner.folder is not None:
        parts.append(f"<p>Folder: <code>{e(str(planner.experiment_folder))}</code></p>")
    b, t = exp.beam, exp.target
    parts.append("<h2>1. The setup</h2><ul>"
                 f"<li>Beam: {e(b.nuclide)} at {e(str(b.energy))}, {e(str(b.current))}.</li>"
                 f"<li>Target: {e(t.material)}, {e(str(t.thickness))}.</li>"
                 f"<li>Particle detectors: {e(', '.join(d.name or f'D{i + 1}' for i, d in enumerate(exp.detectors)))}.</li>"
                 f"<li>γ-ray detectors: {e(', '.join(g.name or f'G{i + 1}' for i, g in enumerate(exp.gamma_detectors))) or 'none'}.</li>"
                 f"<li>Beam time planned: {e(str(exp.run.beam_time))}.</li></ul>")
    if scene_png:
        parts.append(f'<p><img src="data:image/png;base64,{scene_png}" alt="the scene" style="max-width:100%"></p>')
    schemes = level_schemes_html(planner)
    if schemes:
        parts.append("<h3>The level schemes</h3>" + schemes)
    analyses = run_analyses(planner, result)
    history = setup_history(planner)
    if history:
        parts.append("<h3>How the setup changed between runs</h3><ul>")
        for h in history:
            what = f"run {h['run']} ({e(h['label'])})" if h["run"] is not None else e(h["label"])
            if h["run"] is not None and h is history[0]:
                parts.append(f"<li>{what}: the setup as it started.</li>")
            else:
                parts.append(f"<li>{what}: " + (e("; ".join(h["changes"])) if h["changes"] else
                                                 "no change to the physics") + "</li>")
        parts.append("</ul>")
    try:
        pl = planner.plan()
    except Exception:  # noqa: BLE001 - the report stands without the plan
        pl = None
    if pl:
        parts.append("<h2>2. The Plan</h2>")
        parts.append(f"<p>Beam time needed: {_time(pl['beam_time_needed_s'])}"
                     + (f" ({e(str(pl['limiting_detector']))}, for {pl['counts_wanted']} {e(pl['measured'])})"
                        if pl["limiting_detector"] else "") + f"; planned: {_time(pl['beam_time_s'])}.</p>")
        parts.append("<table><thead><tr><th>Detector</th><th>Rate (1/s)</th><th>Busy</th><th>Coincidences (1/s)</th>"
                     "<th>Safe rings</th></tr></thead><tbody>"
                     + "".join(f"<tr><td>{e(d['detector'])}</td><td>{_g(d['rate_per_s'])}</td>"
                               f"<td>{100 * d['dead_time_fraction']:.2g}%</td><td>{_g(d['coincidence_per_s'])}</td>"
                               f"<td>{len(d['safe_rings'])} of {d['rings']}</td></tr>" for d in pl["detectors"])
                     + "</tbody></table>")
    runs = planner.runs()
    parts.append("<h2>3. The runs</h2>")
    if not runs:
        parts.append("<p>No run yet.</p>")
    for s in runs:
        parts.append(f"<h3>Run {s['number']}: {e(s['describe'])}</h3><ul>"
                     f"<li>Folder <code>{e(s['label'])}</code>; seed {s.get('seed')}; finished "
                     f"{e(str(s.get('finished', ''))[:19])} (UTC).</li>")
        if s["kind"] == "source":
            src = s.get("source", {})
            parts.append(f"<li>{e(str(src.get('nuclide')))} source of {e(str(src.get('activity')))} at "
                         f"{e(str(src.get('position') or 'the target'))}; counts per crystal: "
                         + e(", ".join(f"{k} {v:,.0f}" for k, v in s.get("gamma_counts", {}).items())) + ".</li>")
        else:
            parts.append(f"<li>{_time(s['real_s'])} simulated event by event of {_time(s['duration_s'])}"
                         + ("; the rest scaled" if s.get("scaled") else "") + f"; {s.get('particles', 0):,} particles "
                         f"and {s.get('gammas', 0):,} γ rays stored.</li>"
                         "<li>Counts (rate): " + e(", ".join(f"{d} {c:,} ({s['rates_per_s'].get(d, 0):.4g}/s)"
                                                             for d, c in s.get("counts", {}).items())) + ".</li>")
            if s.get("live_fraction") is not None:
                parts.append(f"<li>Live fraction {s['live_fraction']:.4f}.</li>")
            if s.get("assumed_offset_mm") is not None:
                parts.append(f"<li>Alignment run: the target assumed {s['assumed_offset_mm']:+g} mm off.</li>")
        parts.append(f"<li>{'Taken with another setup than the current one.' if s['stale'] else 'Taken with the current setup.'}"
                     "</li></ul>")
        parts.append(_run_figure_html(planner, s, analyses.get(s["number"])))
    if planner.run is not None and planner.run.kind != "source":
        c = planner.compare()
        if c["rows"]:
            parts.append(f"<h3>Run {c['run']}: predicted against measured</h3><table><thead><tr><th>Quantity</th>"
                         "<th></th><th>Predicted</th><th>Measured</th><th></th></tr></thead><tbody>"
                         + "".join(f"<tr><td>{e(r['quantity'])}</td><td>{e(r['what'])}</td><td>{_g(r['predicted'])}"
                                   f"</td><td>{_g(r['measured'])} ± {_g(r['uncertainty'])}</td><td>"
                                   f"{ {True: 'agree', False: 'differ', None: ''}[r['agree']] }</td></tr>"
                                   for r in c["rows"]) + "</tbody></table>")
    gates = planner.gates()
    parts.append("<h2>4. The data: gates and conditions</h2>")
    parts.append("<ul>" + "".join(f"<li><b>{e(g.name)}</b>: {e(g.describe())}.</li>" for g in gates) + "</ul>"
                 if gates else "<p>No gates: the spectra take every counted particle.</p>")
    parts.append("<h2>5. The analysis</h2>")
    for s in runs:
        r = analyses.get(s["number"])
        if r is not None:
            parts.append(analysis_html(r, f"<h3>Run {s['number']}: from the peak area to B(E2)</h3>"))
        if s["kind"] != "source":
            try:
                al = planner.kept_alignment(planner.run if planner.run is not None
                                            and planner.run.number == s["number"]
                                            else planner.folder.load_run(s["number"]))
            except Exception:  # noqa: BLE001 - the report stands without it
                al = None
            if al:
                parts.append((f"<h3>Run {s['number']}: the alignment</h3>" if r is None else "") + alignment_html(al))
    if not analyses and result is not None:  # the analysis of a simulated sample, without a run
        parts.append(analysis_html(result, "<h3>The simulated sample</h3>"))
    elif not analyses:
        parts.append("<p>No run has been analysed yet: take a beam run and press Analyse on the Analysis tab. "
                     "The analysis is then kept with the run and appears here.</p>")
    m = planner.kept_multistep()
    if m and m.get("available"):
        parts.append("<h3>All orders: multi-step excitation and reorientation</h3>" + multistep_html(m))
    xs = run_explanations(planner)
    if xs:
        parts.append("<h2>6. The run's numbers, explained</h2>")
        parts += [_explanation_html(x) for x in xs]
    parts.append("<h2>7. The physics behind it</h2><ul>"
                 + "".join(f'<li><a href="{DOCS}{e(path)}.html">{e(tab.replace("_", " "))}</a></li>'
                           for tab, path in REGISTER.items()) + "</ul>")
    return ("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>"
            f"{e(planner.name)}: the experiment</title><style>{record_css()}</style></head><body>"
            + "\n".join(parts) + "</body></html>")


def experiment_zip(planner, result=None, beam_time_report: bool = True, folder: bool = True) -> bytes:
    """The experiment's record as a zip: ``experiment.html`` (this report), ``record.html`` (the run record),
    ``setup.toml``, the beam-time report of the current setup under ``report/`` (``beam_time_report``), and the
    experiment folder with every run under ``experiment/`` (``folder``)."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("experiment.html", experiment_html(planner, result))
        z.writestr("record.html", planner.record_html(result))
        z.writestr("setup.toml", planner.to_toml())
        if beam_time_report:
            import tempfile

            from .report import build

            with tempfile.TemporaryDirectory() as tmp:
                build(planner.experiment, events=50_000, validation=False).write(tmp, figures=False, root=False)
                for p in sorted(Path(tmp).rglob("*")):
                    if p.is_file():
                        z.write(p, "report/" + p.relative_to(tmp).as_posix())
        if folder and planner.folder is not None:
            root = planner.experiment_folder
            for p in sorted(root.rglob("*")):
                if p.is_file():
                    z.write(p, "experiment/" + p.relative_to(root).as_posix())
    return buf.getvalue()


__all__ = ["alignment_html", "analysis_html", "analysis_png", "experiment_html", "experiment_zip",
           "level_scheme_png", "level_schemes", "level_schemes_html", "multistep_html", "run_analyses",
           "run_explanations", "run_spectra_png", "run_views", "setup_history"]
