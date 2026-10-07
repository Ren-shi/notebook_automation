"""The record of an experiment and its runs (backlog item 69): the report of the experiment, and the explanations of
the numbers the runs add (their counters, the gates, predicted against measured).

The beam-time report (:mod:`physim.nuclear.report`) and the run record (:mod:`physim.nuclear.record`) describe one
setup; this module describes the experiment: its current setup and how it changed from run to run, the Plan, every
run with its summary and the setup it was taken with, the gates of the Data tab, the analysis of the run the user
chose, and the links to the physics register::

    from physim.nuclear import logbook

    page = logbook.experiment_html(planner, result)          # one self-contained page
    data = logbook.experiment_zip(planner, result)           # it, the run record, the setup, the beam-time report,
                                                             # and the experiment folder with every run
    for x in logbook.run_explanations(planner, planner.run):  # each with its formula, numbers and meaning
        assert x.recompute() == x.value
"""

from __future__ import annotations

import html
import io
import math
import zipfile
from pathlib import Path
from typing import Optional

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


def experiment_html(planner, result=None, scene_png: Optional[str] = None) -> str:
    """The report of the experiment as one self-contained page: the setup and how it changed between runs, the
    Plan, every run with its summary, the gates, the analysis of the current run (``result``, a
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
    if result is not None:
        parts.append(f"<p>B(E2↑) = {result.b_e2fm4:.4g} e²fm⁴ ± {100 * result.total_unc:.1f}% "
                     f"(statistical {100 * result.statistical:.1f}%, systematic {100 * result.systematic:.1f}%), "
                     f"from {e(result.statistics_from)}; the value put in was {result.truth_e2fm4:.4g} e²fm⁴ "
                     f"(pull {result.pull:+.2f}).</p><ol>"
                     + "".join(f"<li><b>{e(st['step'])}</b>: {e(st['text'])}</li>" for st in result.steps) + "</ol>")
    else:
        parts.append("<p>No analysis has been run on the current run.</p>")
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


__all__ = ["experiment_html", "experiment_zip", "run_explanations", "setup_history"]
