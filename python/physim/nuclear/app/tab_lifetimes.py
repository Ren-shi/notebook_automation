"""Lifetime measurements in the planner app: a plunger (RDDS) or line shapes (DSAM), after Coulomb excitation or
fusion-evaporation, through the steps.

- Plan: the lifetimes panel from :meth:`~physim.nuclear.planner.Planner.lifetimes` (per level its method, τ and
  where τ comes from, the recommended distances, the separation of the shifted and unshifted peaks at each
  detector, F(τ) for DSAM, the beam time per distance; problems and warnings) and the predicted spectra from
  :meth:`~physim.nuclear.planner.Planner.lifetime_spectra`.
- Run: the plunger distance of a run, and the series of runs one per distance.
- Data: the runs of a plunger setup grouped by distance, with the spectra per distance and ring
  (:func:`physim.nuclear.dataviews.distance_spectra`); a fusion run's singles, γγ matrix and gated spectra.
- Analysis: the decay curve and the DDCM τ(x), or the DSAM line-shape fit, from
  :meth:`~physim.nuclear.planner.Planner.analyse_lifetimes`, with τ, its budget and B(Eλ); then τ and the B(E2)
  that follows from it beside what is known about the nucleus, in short (:func:`.tab_analysis.short_comparison`,
  as under a Coulex result). After the app is opened again the analysis kept with the experiment is shown
  (:meth:`~physim.nuclear.planner.Planner.kept_lifetimes`), marked as kept, and is not run again until asked.

Each planner function is read through :func:`call`: a physim without it (or one that fails) gives a plain
``{"available": False, "reason": ...}``, which the panels show as a note. The panels read only the keys of the
shapes frozen for wave 3 of the lifetime work."""

from __future__ import annotations

import math
import os
import time
from typing import Optional

import numpy as np

from .. import workbench as wb
from .tab_plan import _fmt, _time

#: What a method is called on the page.
METHODS = {"rdds": "Plunger (RDDS)", "dsam": "Line shapes (DSAM)", "outside": "Neither"}
ANALYSES = {"ddcm": "Differential decay curve (DDCM)", "decay_curve": "Decay curve fit",
            "ddcm_three_foil": "DDCM, three foils", "lineshape": "Line-shape fit"}


def call(planner, name: str, *args, **kwargs) -> dict:
    """``planner.<name>(...)`` as a dict with ``available``: ``{"available": False, "missing": True, "reason"}``
    when this physim has no such method yet, ``{"available": False, "reason"}`` when it raises."""
    if name == "lifetimes" and not args and not kwargs:
        return wb.lifetime_plan(planner)
    fn = getattr(planner, name, None)
    if fn is None:
        return {"available": False, "missing": True,
                "reason": f"Planner.{name} is not in this version of physim yet: it comes with the planner's "
                          "lifetime methods."}
    try:
        out = fn(*args, **kwargs)
    except Exception as err:  # noqa: BLE001 - the panel says why
        return {"available": False, "reason": f"Planner.{name} failed: {err}"}
    return out if isinstance(out, dict) else {"available": False, "reason": f"Planner.{name} gave {type(out)}"}


def distance_spectra(series, edges_kev=None, doppler: bool = True) -> dict:
    """:func:`physim.nuclear.dataviews.distance_spectra` with ``available``, guarded as :func:`call` is."""
    from .. import dataviews as dv

    fn = getattr(dv, "distance_spectra", None)
    if fn is None:
        return {"available": False, "missing": True,
                "reason": "dataviews.distance_spectra is not in this version of physim yet: the spectra per "
                          "distance and ring come with the planner's lifetime methods."}
    try:
        out = dict(fn(series, edges_kev=edges_kev, doppler=doppler))
    except Exception as err:  # noqa: BLE001
        return {"available": False, "reason": f"The spectra per distance could not be made: {err}"}
    out.setdefault("available", True)
    return out


def _unit(tau_ps) -> tuple:
    """(unit, factor from ps) for a lifetime: fs below 1 ps, ns from 1000 ps."""
    return ("fs", 1e3) if tau_ps < 1 else (("ns", 1e-3) if tau_ps >= 1e3 else ("ps", 1.0))


def _tau(tau_ps, err=None, like=None) -> str:
    """τ in its natural unit (fs, ps or ns), or in that of ``like``, with its uncertainty."""
    if tau_ps is None:
        return "—"
    unit, f = _unit(like if like is not None else tau_ps)
    text = f"{tau_ps * f:.3g}"
    if err is not None and err > 0:
        text += f" ± {err * f:.2g}"
    return f"{text} {unit}"


def _texts(items) -> list:
    out = []
    for x in items or []:
        out.append(x if isinstance(x, str) else (x.get("text") if isinstance(x, dict) else str(x)))
    return [t for t in out if t]


def _note(ui, text: str) -> None:
    ui.label(text).classes("text-sm ps-note rounded p-2 w-full")


def _checks(ui, plan: dict) -> None:
    for level, key in (("error", "problems"), ("warning", "warnings")):
        for text in _texts(plan.get(key)):
            with ui.row().classes("items-start no-wrap gap-2"):
                ui.icon("error" if level == "error" else "warning").classes(
                    "ps-bad" if level == "error" else "ps-warn")
                ui.label(text).classes(("ps-bad" if level == "error" else "ps-warn") + " text-sm")


def _readouts(ui, cells: list) -> None:
    with ui.element("div").classes("ps-readouts"):
        for title, value, sub, colour in cells:
            with ui.element("div").classes("ps-readout"):
                ui.label(title).classes("text-xs ps-muted")
                ui.label(value).classes("ps-readout-value " + colour)
                ui.label(sub).classes("text-xs ps-muted")


def setup_line(exp) -> str:
    """The lifetime measurement as the setup gives it, in a line."""
    lt = exp.lifetime
    if lt is None:
        return ""
    if lt.method == "dsam":
        b = exp.target.backing
        return ("Line shapes (DSAM): the recoil slows down in the "
                + (f"{b.thickness} {b.material} backing." if b is not None else "target."))
    ds = lt.distances_um
    text = f"Plunger (RDDS): {len(ds)} distances, {min(ds):g} to {max(ds):g} µm" if ds else "Plunger (RDDS)"
    if lt.offset_um:
        text += f", + {lt.offset_um:g} µm at contact"
    if lt.stopper is not None:
        text += f"; {lt.stopper.thickness} {lt.stopper.material} {lt.stopper.kind}"
    if lt.degrader2 is not None:
        text += f", then {lt.degrader2.thickness} {lt.degrader2.material} {lt.degrader2.gap} behind it"
    return text + "."


# -- Plan ---------------------------------------------------------------------------------------------------------

def plan_panel(ctx) -> None:
    """The Plan step's lifetimes panel (refreshed with the plan); nothing without a lifetime measurement."""
    ui, p = ctx.ui, ctx.P()
    exp = p.experiment
    if exp.lifetime is None:
        return
    plan = call(p, "lifetimes")
    sub = "τ per level, distances, peak separation, beam time per distance" if exp.lifetime.method == "rdds" else \
        "τ per level, the attenuation F(τ), beam time"
    with ctx.panel(f"Lifetimes · {METHODS[exp.lifetime.method]}", sub):
        ui.label(setup_line(exp)).classes("text-sm")
        if not plan.get("available"):
            _note(ui, plan.get("reason") or "No lifetime plan.")
            try:
                ui.label("The recoil: " + wb._recoil_line(p) + ".").classes("text-sm ps-muted")
            except Exception:  # noqa: BLE001 -- the setup's own numbers, when they can be had
                pass
            return
        render_plan(ctx, plan)
    spectra_block(ctx)


def render_plan(ctx, plan: dict) -> None:
    """The lifetimes plan (the frozen shape of ``Planner.lifetimes``): read-outs for the studied level, the
    levels, the studied level's peak separation per detector, and the checks."""
    ui = ctx.ui
    levels = plan.get("levels") or []
    first = levels[0] if levels else None
    cells = []
    if first is not None:
        cells.append(("Studied level", _tau(first.get("tau_ps")),
                      f"{first.get('level', '?')}, {first.get('gamma_kev', 0):.1f} keV γ ray"
                      + (f"; τ from {first['tau_source']}" if first.get("tau_source") else ""), ""))
        setup_method = plan.get("method")
        fits = setup_method is None or first.get("method") in (setup_method, None)
        short = {"rdds": "RDDS", "dsam": "DSAM", "outside": "neither"}
        cells.append(("Method for it", short.get(first.get("method"), first.get("method") or "—"),
                      ("as the setup has it" if setup_method == first.get("method") else
                       f"the setup has {METHODS.get(setup_method, setup_method)}") if setup_method else
                      "no method in the setup yet", "" if fits else "ps-warn"))
        sep = first.get("separation") or []
        if first.get("method") == "rdds" and sep:
            apart = sum(bool(s.get("resolved")) for s in sep)
            cells.append(("Peaks apart", f"{apart} of {len(sep)}", "detectors where the shifted peak clears the "
                          "unshifted one", "ps-ok" if apart == len(sep) else ("ps-bad" if not apart else "ps-warn")))
        if first.get("attenuation") is not None:
            cells.append(("Attenuation F(τ)", f"{first['attenuation']:.3f}",
                          "share of the full Doppler shift the line keeps", ""))
        advice = first.get("backing_advice") or []
        if advice:
            best = advice[0]
            mine = next((c for c in advice if c.get("setup")), None)
            sub = f"{best['thickness_mg_cm2']:.2g} mg/cm² to stop it, F = {best['attenuation']:.2f}"
            colour = ""
            if mine is not None and mine is not best:
                sub += f"; the setup's {mine['material']}: {100 * mine['relative_sensitivity']:.0f}% as sensitive"
                colour = "ps-ok" if mine.get("good") else "ps-warn"
            cells.append(("Best backing", best["material"], sub, colour))
        if first.get("hours_per_distance") is not None:
            what = "per distance" if first.get("method") == "rdds" else "in all"
            cells.append(("Beam time " + what, _time(3600 * first["hours_per_distance"]),
                          f"for {100 * (plan.get('precision') or 0.1):.0f} % on τ", ""))
        if first.get("beta") is not None:
            ttl = plan.get("time_to_stop_ps")
            cells.append(("Recoil", f"β {first['beta']:.4f}",
                          f"{_tau(ttl)} to stop" if ttl is not None else "", ""))
    if cells:
        _readouts(ui, cells)
    if levels:
        levels_table(ctx, plan)
    if first is not None and first.get("separation"):
        separation_table(ctx, first)
    distances_line(ctx, plan, first)
    advised = next((lev for lev in levels if lev.get("backing_advice")), None)
    if advised is not None:
        backing_table(ctx, advised)
    if plan.get("stopper"):
        stopper_table(ctx, plan)
    _checks(ui, plan)


def backing_table(ctx, lev: dict) -> None:
    """The DSAM backing advice for a level (``backing_advice`` of ``Planner.lifetimes``): the usual backings
    ranked by how sensitive F(τ) is to τ on each."""
    ui = ctx.ui
    rows = []
    for k, c in enumerate(lev["backing_advice"]):
        um = (lambda x: f" ({x:.2g} µm)" if x else "")
        rows.append({
            "k": k, "m": ("● " if c.get("setup") else "") + c["material"],
            "f": f"{c['attenuation']:.3f}", "s": f"{c['sensitivity']:.3f}",
            "r": f"{100 * c['relative_sensitivity']:.0f}%" + (" ✓" if c.get("good") else ""),
            "a": _tau(c.get("slowing_time_ps")), "t": _tau(c.get("time_to_stop_ps")),
            "stop": f"{c['thickness_to_stop_mg_cm2']:.3g}{um(c.get('thickness_to_stop_um'))}",
            "adv": f"{c['thickness_mg_cm2']:.2g}{um(c.get('thickness_um'))}"})
    ui.label(f"DSAM backing for {lev.get('level', '')} (τ = {_tau(lev.get('tau_ps'))})").classes("ps-section mt-2")
    ui.table(columns=ctx.columns((("m", "Backing"), ("f", "F(τ)"), ("s", "|dF/d ln τ|"), ("r", "Of the best"),
                                  ("a", "Slowing time"), ("t", "Time to stop"), ("stop", "To stop (mg/cm²)"),
                                  ("adv", "Advised (mg/cm²)"))),
             rows=rows, row_key="k").props("dense flat").classes("w-full")
    ui.label("Ranked by |dF/d ln τ|: σ_τ/τ = σ_F / |dF/d ln τ|, best (about 1/4) when the slowing-down time matches "
             "τ and F is near ½. ✓ within 80% of the best; ● the setup's backing material. Advised: 1.5 times the "
             "thickness that stops the typical recoil.").classes("text-xs ps-muted")


def stopper_table(ctx, plan: dict) -> None:
    """The beam in the plunger's foils (``stopper`` of ``Planner.lifetimes``): its power and the foils' own
    Coulomb-excitation lines."""
    ui = ctx.ui
    rows = []
    for k, f in enumerate(plan["stopper"]):
        lines = ", ".join(f"{x['nuclide']} {x['gamma_kev']:.0f}" for x in (f.get("lines") or [])[:3])
        if f.get("missing"):
            lines = (lines + "; " if lines else "") + "no level data: " + ", ".join(
                m["nuclide"] for m in f["missing"])
        over = f["power_w"] > f["power_limit_w"]
        rows.append({"k": k, "foil": f"{f['thickness_mg_cm2']:.3g} mg/cm² {f['material']} {f['foil']}",
                     "e": f"{f['beam_in_mev']:.3g} → {f['beam_out_mev']:.3g}",
                     "p": f"{f['power_w']:.2g} W" + (" ⚠" if over else ""), "lim": f"{f['power_limit_w']:g} W",
                     "lines": lines or "—"})
    ui.label("The beam in the plunger's foils").classes("ps-section mt-2")
    ui.table(columns=ctx.columns((("foil", "Foil"), ("e", "Beam in → out (MeV)"), ("p", "Power"),
                                  ("lim", "Limit"), ("lines", "Its Coulex lines (keV)"))),
             rows=rows, row_key="k").props("dense flat").classes("w-full")
    ui.label("The beam that passes the target stops in (or crosses) the plunger's foils: it heats them (a thin "
             "stretched foil cooled by radiation takes about half a watt) and Coulomb-excites their nuclei, whose "
             "lines are warned about when they sit near the studied ones.").classes("text-xs ps-muted")


def levels_table(ctx, plan: dict) -> None:
    ui = ctx.ui
    rows = []
    for k, lev in enumerate(plan["levels"]):
        sep = lev.get("separation") or []
        rec = lev.get("recommended_distances_um") or []
        sens = lev.get("sensitive_um")
        rows.append({
            "k": k, "level": ("● " if k == 0 else "") + str(lev.get("level", "?")),
            "e": f"{lev['gamma_kev']:.1f}" if lev.get("gamma_kev") is not None else "—",
            "tau": _tau(lev.get("tau_ps")), "src": lev.get("tau_source") or "—",
            "beta": f"{lev['beta']:.4f}" if lev.get("beta") is not None else "—",
            "method": METHODS.get(lev.get("method"), lev.get("method") or "—"),
            "rec": (f"{min(rec):.3g}–{max(rec):.3g} ({len(rec)})" if rec else "—"),
            "sens": f"{sens[0]:.3g}–{sens[1]:.3g}" if sens else "—",
            "sep": f"{sum(bool(s.get('resolved')) for s in sep)} of {len(sep)}" if sep else "—",
            "f": f"{lev['attenuation']:.3f}" if lev.get("attenuation") is not None else "—",
            "h": f"{lev['hours_per_distance']:.3g}" if lev.get("hours_per_distance") is not None else "—"})
    ui.label("Levels").classes("ps-section mt-2")
    ui.table(columns=ctx.columns((("level", "Level"), ("e", "Eγ (keV)"), ("tau", "τ"), ("src", "τ from"),
                                  ("beta", "β"), ("method", "Method"), ("rec", "Distances advised (µm)"),
                                  ("sens", "Sensitive (µm)"), ("sep", "Peaks apart"), ("f", "F(τ)"),
                                  ("h", "Hours per distance"))),
             rows=rows, row_key="k").props("dense flat").classes("w-full")
    ui.label("The studied level first (●). Distances advised: log-spaced around vτ; sensitive: where the "
             "unshifted share changes fastest with the distance. A level whose τ is outside what a plunger or "
             "DSAM measures, or not known, is marked Neither.").classes("text-xs ps-muted")


def separation_table(ctx, lev: dict) -> None:
    ui = ctx.ui
    rows = [{"d": s.get("detector", "?"), "t": f"{s['theta_deg']:.1f}" if s.get("theta_deg") is not None else "—",
             "s": f"{s['separation_kev']:.2f}" if s.get("separation_kev") is not None else "—",
             "w": f"{s['fwhm_kev']:.2f}" if s.get("fwhm_kev") is not None else "—",
             "r": "✓ apart" if s.get("resolved") else "✗ overlap"} for s in lev["separation"]]
    ui.label(f"Shifted against unshifted peak, {lev.get('level', '')} → {lev.get('gamma_kev', 0):.1f} keV").classes(
        "ps-section mt-2")
    ui.table(columns=ctx.columns((("d", "Detector"), ("t", "θ (deg)"), ("s", "Separation (keV)"),
                                  ("w", "FWHM (keV)"), ("r", ""))), rows=rows).props("dense flat")


def distances_line(ctx, plan: dict, first) -> None:
    """The setup's distances beside those advised for the studied level."""
    ds = plan.get("distances_um") or []
    if not ds or first is None or first.get("method") != "rdds":
        return
    rec = first.get("recommended_distances_um") or []
    text = f"The setup's distances: {', '.join(f'{d:g}' for d in ds)} µm."
    if rec:
        text += f" Advised for {first.get('level', 'the studied level')}: {', '.join(f'{d:.3g}' for d in rec)} µm."
        sens = first.get("sensitive_um")
        if sens:
            inside = sum(sens[0] <= d <= sens[1] for d in ds)
            text += f" {inside} of the setup's {len(ds)} lie in the sensitive range."
    ctx.ui.label(text).classes("text-xs ps-muted mt-1")


def spectra_of(ctx, p) -> dict:
    """The predicted spectra of the setup (``Planner.lifetime_spectra``), computed once per setup and kept."""
    return ctx.keep(p, ("lifetimes", "spectra"), lambda: call(p, "lifetime_spectra"))


def spectra_block(ctx) -> None:
    """The predicted spectra, at each distance or as line shapes: built when opened, computed on the worker."""
    ui, p = ctx.ui, ctx.P()
    rdds = p.experiment.lifetime.method == "rdds"
    title = "Predicted spectra at each plunger distance" if rdds else "Predicted line shapes"
    exp = ui.expansion(title, icon="ssid_chart").classes("w-full ps-card")
    ctx.when_opened(exp, lambda: render_spectra(ctx, spectra_of(ctx, ctx.P())), "lifetime_spectra")


def render_spectra(ctx, spec: dict) -> None:
    from .figures import figure_lifetime_spectra

    ui = ctx.ui
    if not spec.get("available"):
        _note(ui, spec.get("reason") or "No predicted spectra.")
        return
    distances = sorted({s["distance_um"] for s in spec.get("spectra", []) if s.get("distance_um") is not None})
    box = ui.column().classes("w-full gap-1")
    chosen = {"d": distances[0] if distances else None}

    def draw() -> None:
        box.clear()
        with box:
            ui.plotly(ctx.in_theme(figure_lifetime_spectra(spec, chosen["d"]))).classes("w-full")

    if distances:
        ui.select({d: f"{d:g} µm" for d in distances}, value=chosen["d"], label="Distance",
                  on_change=lambda e: (chosen.update(d=e.value), draw())).props("dense outlined").classes("w-40")
    draw()
    note = ("One line per detector angle: the shifted peak (dashed) holds the decays in flight, the unshifted one "
            "(dotted) those after the stopper." if distances else
            "One line per detector angle: the line shape between the full Doppler shift and rest.")
    ui.label(note + (" " + " ".join(_texts(spec.get("notes"))) if spec.get("notes") else "")).classes(
        "text-xs ps-muted")


# -- Run ----------------------------------------------------------------------------------------------------------

def plunger(exp) -> bool:
    return exp is not None and exp.lifetime is not None and exp.lifetime.method == "rdds"


def taken_distances(p) -> dict:
    """{distance µm: number of runs} of the current setup's runs (not stale), for the run form."""
    out: dict = {}
    for r in p.runs():
        if r.get("stale") or r.get("distance_um") is None or r.get("kind") in ("source", "alignment"):
            continue
        d = float(r["distance_um"])
        out[d] = out.get(d, 0) + 1
    return out


def distance_controls(ctx, f: dict) -> None:
    """The run form's plunger distance, and what the series has so far."""
    ui, p = ctx.ui, ctx.P()
    ds = p.experiment.lifetime.distances_um
    if not ds:
        ui.label("The plunger has no distances: set them in the Lifetime measurement card.").classes("text-sm ps-warn")
        return
    taken = taken_distances(p)
    if f.get("distance") not in ds:  # the first distance without a run, else the first
        f["distance"] = next((d for d in ds if d not in taken), ds[0])
    options = {d: f"{d:g} µm" + (f" · {taken[d]} run{'s' * (taken[d] != 1)}" if d in taken else "") for d in ds}
    with ui.row().classes("items-end gap-2 w-full"):
        ui.select(options, value=f["distance"], label="Plunger distance",
                  on_change=lambda e: f.update(distance=e.value)).props("dense outlined").classes("w-48")
        missing = [d for d in ds if d not in taken]
        if missing and not p.running:
            ui.button(f"Run each distance not taken ({len(missing)})", icon="playlist_play",
                      on_click=lambda: ctx.start_series(dict(f), missing)).props("flat no-caps").tooltip(
                "One run per distance, one after the other, each with the duration and real part above")
    queue = ctx.state.get("distance_queue") or []
    line = f"{len(taken)} of {len(ds)} distances taken"
    if queue:
        line += f"; {len(queue)} more queued: {', '.join(f'{d:g}' for d in queue)} µm"
    ui.label(line + ". The distance is the plunger's dial (without the offset); each run is one point of the decay "
             "curve.").classes("text-xs ps-muted")


def fusion_live(ctx, g: dict) -> None:
    """The live counters of a fusion run: residues made, events the trigger kept, per channel and per crystal."""
    ui = ctx.ui
    trig = g["counts"].get("trigger", 0)
    _readouts(ui, [("Residues made", f"{g['residues']:,}", f"{_fmt(g['residues_per_s'])} /s", ""),
                   ("Kept by the trigger", f"{trig:,}", f"{_fmt(g['rates_per_s'].get('trigger', 0))} /s, "
                    f"γ-ray fold ≥ {g.get('trigger', {}).get('fold', '?')}", ""),
                   ("By the end", _fmt(g["projected"].get("trigger", 0)), "events in the whole run", "")]
             + ([("Plunger distance", f"{g['distance_um']:g} µm", "", "")] if g.get("distance_um") is not None
                else []))
    chans = g.get("channels") or {}
    if chans:
        ui.table(columns=ctx.columns((("c", "Channel"), ("r", "Residues"), ("t", "Kept"))),
                 rows=[{"c": c, "r": f"{v['residues']:,}", "t": f"{v['triggered']:,}"} for c, v in chans.items()]
                 ).props("dense flat")
    if g.get("gamma_singles"):
        ui.table(columns=ctx.columns((("c", "Crystal"), ("s", "γ rays, all residues"), ("k", "In kept events"))),
                 rows=[{"c": c, "s": f"{n:,}", "k": f"{g['gamma'].get(c, 0):,}"}
                       for c, n in g["gamma_singles"].items()]).props("dense flat").classes("max-h-64")


def fusion_sample(ctx, sample: list) -> None:
    ui = ctx.ui
    ui.label("The last events kept").classes("ps-section mt-2")
    ui.table(columns=ctx.columns((("e", "Event"), ("c", "Channel"), ("f", "Fold"), ("g", "γ rays (keV)"))),
             rows=[{"e": x["event"], "c": x["channel"], "f": x["fold"],
                    "g": ", ".join(f"{e:.1f}" for e in x["energies_kev"])} for x in sample]).props("dense flat")


# -- Data ---------------------------------------------------------------------------------------------------------

def series_of(p, run):
    """The distance series the run belongs to (:class:`physim.nuclear.runs.DistanceSeries`), or None."""
    if run is None or run.distance_um is None or p.folder is None:
        return None
    for s in p.folder.distance_series():
        for runs in s.runs.values():
            if any(str(r.folder) == str(run.folder) for r in runs):
                return s
    return None


def _signature(series) -> tuple:
    return tuple(sorted((str(r.folder), str(r.summary.get("finished"))) for runs in series.runs.values()
                        for r in runs))


def series_view(ctx, p) -> tuple:
    """(series, its spectra per distance and ring) for the current run, kept by the series' runs and the
    choice of energies; (None, None) for a run without a plunger distance."""
    run = p.run
    series = ctx.keep(p, ("lifetimes", "series", str(run.folder) if run is not None else None,
                          len(p.runs())), lambda: series_of(p, run))
    if series is None:
        return None, None
    doppler = ctx.data_options.setdefault("lt_doppler", True)
    view = ctx.keep(p, ("lifetimes", "distance spectra", _signature(series), doppler),
                    lambda: distance_spectra(series, doppler=doppler))
    return series, view


def fusion_range(run) -> tuple:
    """The γ-ray energy range of a fusion run's spectra (keV): up to the highest corrected energy recorded."""
    g = run.gammas()
    top = float(np.max(g["corrected"])) * 1e3 if len(g) else 2000.0
    return 0.0, float(min(max(100 * math.ceil(1.05 * top / 100), 500), 6000))


def fusion_views(ctx, p) -> dict:
    """A fusion run's views for the current choices, kept with the run: the singles per γ-ray detector, the γγ
    matrix and the gated spectrum."""
    run = p.run
    o = ctx.data_options
    corrected, channel = o.setdefault("fusion_corrected", True), o.setdefault("fusion_channel", None)
    bins = int(o.get("bins") or 200)
    lo, hi = ctx.keep(p, ("fusion", "range"), lambda: fusion_range(run), of_run=True)
    names = [gd.name or f"γ{i + 1}" for i, gd in enumerate(run.experiment.gamma_detectors)]
    singles = []
    for n in names:
        def make(n=n) -> dict:
            counts, edges = run.gammas().spectrum(n, corrected=corrected, channel=channel, bins=bins,
                                                  range_kev=(lo, hi))
            return dict(counts=counts, edges=edges, label=n, xlabel="γ-ray energy (keV)")
        try:
            singles.append(ctx.keep(p, ("fusion", "singles", n, corrected, channel, bins), make, of_run=True))
        except KeyError:  # a detector with no crystal of its own in the run
            continue
    mbins = int(o.setdefault("fusion_matrix_bins", 400))
    matrix = ctx.keep(p, ("fusion", "matrix", corrected, mbins),
                      lambda: run.matrix(bins=mbins, range_kev=(lo, hi), corrected=corrected), of_run=True)
    gate = o.get("fusion_gate")
    width = float(o.setdefault("fusion_gate_width", 3.0))
    gated = None
    if gate:
        gated = ctx.keep(p, ("fusion", "gated", float(gate), width, corrected, bins),
                         lambda: run.gated(float(gate), width, bins=bins, range_kev=(lo, hi), corrected=corrected),
                         of_run=True)
    return {"singles": singles, "matrix": matrix, "gated": gated, "range": (lo, hi)}


def prepare_data(ctx, p) -> None:
    """What the Data step's lifetime and fusion views show, computed on the worker (``p`` is a snapshot)."""
    run = p.run
    if run is None or run.kind == "source":
        return
    if run.kind == "fusion":
        fusion_views(ctx, p)
    series_view(ctx, p)


def series_block(ctx, run) -> None:
    """The distance series of the current run's setup: its runs by distance, the distances still missing, and the
    spectra per distance and ring."""
    from .figures import figure_distance_spectra

    ui, p = ctx.ui, ctx.P()
    series, view = series_view(ctx, p)
    if series is None:
        return
    o = ctx.data_options
    missing = series.missing()
    with ctx.panel("Distance series", f"{len(series.distances)} distance{'s' * (len(series.distances) != 1)} taken"
                   + (f", {len(missing)} to go" if missing else "")):
        rows = []
        for d in series.distances:
            runs = series.runs[d]
            rows.append({"d": f"{d:g}", "runs": ", ".join(str(r.number) for r in runs),
                         "counts": _fmt(sum(float(r.counts()) if r.kind == "fusion" else
                                            float(sum((r.summary.get("counts") or {}).values())) for r in runs)),
                         "time": _time(sum(r.duration_s for r in runs)),
                         "now": "● current" if any(str(r.folder) == str(run.folder) for r in runs) else ""})
        ui.table(columns=ctx.columns((("d", "Distance (µm)"), ("runs", "Runs"), ("counts", "Events (real part)"),
                                      ("time", "Beam time"), ("now", ""))), rows=rows).props("dense flat")
        if missing:
            ui.label(f"Not taken yet: {', '.join(f'{d:g}' for d in missing)} µm (the Run step takes them).").classes(
                "text-xs ps-muted")
        if not view or not view.get("available", True):
            _note(ui, (view or {}).get("reason") or "No spectra per distance.")
            return
        first = next(iter(view.get("by_distance", {}).values()), {})
        thetas = list(first.get("theta_deg", []))
        ring = min(int(o.setdefault("lt_ring", 0)), max(len(thetas) - 1, 0))
        normalise = o.setdefault("lt_normalise", True)
        with ui.row().classes("items-end gap-2"):
            if thetas:
                ui.select({k: f"ring at θ {t:.1f}°" for k, t in enumerate(thetas)}, value=ring, label="Ring",
                          on_change=lambda e: (o.update(lt_ring=e.value), ctx.refresh_data())).props(
                    "dense outlined").classes("w-44")
            ui.toggle({True: "as measured (shifted)", False: "Doppler corrected"}, value=o["lt_doppler"],
                      on_change=lambda e: (o.update(lt_doppler=e.value), ctx.refresh_data())).props("dense no-caps")
            ui.switch("Each to its own total", value=normalise,
                      on_change=lambda e: (o.update(lt_normalise=e.value), ctx.refresh_data())).props("dense")
        ui.plotly(ctx.in_theme(figure_distance_spectra(view, ring, normalise=normalise))).classes("w-full")
        ui.label("Each line is one distance, from the shortest (dark blue) to the longest (orange): as the "
                 "distance grows, more of the decays happen in flight and the shifted peak takes over from the "
                 "unshifted one. Counts in the runs' real parts.").classes("text-xs ps-muted")


def fusion_data(ctx, run) -> None:
    """A fusion run on the Data step: the singles per γ-ray detector, the γγ matrix, a gated spectrum."""
    from .figures import figure_gg_matrix, figure_histogram

    ui, p = ctx.ui, ctx.P()
    o = ctx.data_options
    views = fusion_views(ctx, p)
    with ui.row().classes("items-end gap-2 w-full"):
        ui.toggle({True: "Doppler corrected (β of the compound nucleus)", False: "as measured"},
                  value=o["fusion_corrected"],
                  on_change=lambda e: (o.update(fusion_corrected=e.value), ctx.refresh_data())).props("dense no-caps")
        ui.select({"": "All channels"} | {c: c for c in run.channels}, value=o["fusion_channel"] or "",
                  label="Channel (as simulated)",
                  on_change=lambda e: (o.update(fusion_channel=e.value or None), ctx.refresh_data())).props(
            "dense outlined").classes("w-44")
        ui.number("Bins", value=o.get("bins") or 200, min=20, max=4000, step=50,
                  on_change=lambda e: o.update(bins=int(e.value or 200))).props("dense outlined").classes(
            "w-24").on("blur", ctx.refresh_data)
    ui.label("The channel is the one each event was simulated in: a perfect selection. The ancillaries' gates "
             "are in the run's trigger.").classes("text-xs ps-muted")
    with ui.element("div").classes("w-full").style("display: flex; flex-wrap: wrap; gap: 8px"):
        for spec in views["singles"]:
            with ui.column().classes("ps-plate gap-0").style("flex: 1 1 300px; min-width: 0"):
                with ui.row().classes("w-full items-center no-wrap"):
                    ui.label(spec["label"]).classes("text-sm font-medium")
                    ui.label(f"{spec['counts'].sum():,.0f} counts").classes("text-xs ps-muted")
                ui.plotly(ctx.keep(p, ("fusion", "small", ctx.state["theme"], spec["label"], o["fusion_corrected"],
                                       o["fusion_channel"], int(o.get("bins") or 200)),
                                   lambda s=spec: ctx.in_theme(figure_histogram([s], height=190)), of_run=True)
                          ).classes("w-full")
    ui.label("γγ coincidences").classes("text-lg font-semibold mt-4")
    both = ui.element("div").classes("w-full grid gap-3").style(
        "grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); align-items: start")
    with both:
        with ui.column().classes("w-full gap-1").style("min-width: 0"):
            counts, edges = views["matrix"]
            ui.plotly(ctx.in_theme(figure_gg_matrix(counts, edges))).classes("w-full")
            ui.label(f"Every pair of γ rays in one kept event, both ways round: {counts.sum() / 2:,.0f} pairs in the "
                     "real part.").classes("text-xs ps-muted")
        with ui.column().classes("w-full gap-1").style("min-width: 0"):
            with ui.row().classes("items-end gap-2"):
                ui.number("Gate on (keV)", value=o.get("fusion_gate"),
                          on_change=lambda e: o.update(fusion_gate=e.value)).props("dense outlined").classes(
                    "w-32").on("blur", ctx.refresh_data)
                ui.number("± (keV)", value=o["fusion_gate_width"], min=0.1, step=0.5,
                          on_change=lambda e: o.update(fusion_gate_width=e.value or 3.0)).props(
                    "dense outlined").classes("w-24").on("blur", ctx.refresh_data)
            if views["gated"] is None:
                ui.label("Give a γ-ray energy to see what is in coincidence with it: gated on a band's "
                         "transition, the band shows.").classes("text-sm ps-muted")
            else:
                counts, edges = views["gated"]
                ui.plotly(ctx.in_theme(figure_histogram(
                    [dict(counts=counts, edges=edges, label="gated", xlabel="γ-ray energy (keV)")],
                    height=320))).classes("w-full")
                ui.label(f"In coincidence with a γ ray at {float(o['fusion_gate']):g} ± {o['fusion_gate_width']:g} "
                         f"keV: {counts.sum():,.0f} γ rays (the gating one left out).").classes("text-xs ps-muted")


# -- Analysis -----------------------------------------------------------------------------------------------------

def analysis_choice(ctx) -> dict:
    """The level and method chosen for the analysis, and whether one was run in this session ("session"); begun
    afresh for another experiment."""
    of = str(getattr(ctx.state.get("planner"), "experiment_folder", None))
    if ctx.state.get("lifetime_analysis", {}).get("of") != of:
        ctx.state["lifetime_analysis"] = {"level": None, "method": None, "session": False, "of": of}
    return ctx.state["lifetime_analysis"]


def kept_analysis(p) -> Optional[dict]:
    """The lifetime analysis kept with the experiment from an earlier session
    (:func:`physim.nuclear.workbench.kept_lifetimes`), with "kept" (when it was written) added; None if there is
    none, the setup's physics changed since, or a run was taken that it did not read."""
    kept = wb.kept_lifetimes(p)
    if kept is None:
        return None
    try:
        if wb.runs_not_read(p, kept):
            return None
    except Exception:  # noqa: BLE001 -- the runs cannot be listed: the kept analysis is not vouched for
        return None
    try:
        when = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(
            os.path.join(str(p.experiment_folder), "lifetime.json"))))
    except OSError:
        when = ""
    return dict(kept, kept=when or "earlier")


def analysis_of(ctx, p) -> dict:
    """The lifetime analysis the step shows for the level and method chosen: the one kept with the experiment
    from an earlier session (:func:`kept_analysis`) until one is run in this session; else
    ``Planner.analyse_lifetimes``, kept by the runs it reads."""
    choice = analysis_choice(ctx)
    if not choice.get("session") and choice["level"] is None and choice["method"] is None:
        kept = kept_analysis(p)
        if kept is not None:
            return kept
    runs = tuple((r["number"], str(r.get("finished"))) for r in p.runs()) if p.folder is not None else ()

    def analyse() -> dict:
        choice["session"] = True
        return call(p, "analyse_lifetimes", level=choice["level"], method=choice["method"])

    return ctx.keep(p, ("lifetimes", "analysis", choice["level"], choice["method"], runs), analyse)


def prepare_analysis(ctx, p) -> None:
    if p.experiment.lifetime is not None and getattr(p, "analyse_lifetimes", None) is not None:
        analysis_of(ctx, p)


def analysis_panel(ctx) -> None:
    """The Analysis step's lifetimes part."""
    ui, p = ctx.ui, ctx.P()
    exp = p.experiment
    if exp.lifetime is None:
        return
    ui.label("Lifetimes").classes("text-lg font-semibold")
    ui.label(setup_line(exp)).classes("text-sm ps-muted")
    plan = call(p, "lifetimes")
    choice = analysis_choice(ctx)
    with ui.row().classes("items-end gap-2"):
        levels = {lev.get("level"): f"{lev.get('level')} ({lev['gamma_kev']:.1f} keV)" if lev.get("gamma_kev")
                  else str(lev.get("level")) for lev in plan.get("levels") or [] if lev.get("level") is not None}
        if levels:
            ui.select({None: "the studied level", **levels}, value=choice["level"], label="Level",
                      on_change=lambda e: rechoose(ctx, level=e.value)).props("dense outlined").classes("w-56")
        if exp.lifetime.method == "rdds":
            methods = {None: "the default", "ddcm": ANALYSES["ddcm"], "decay_curve": ANALYSES["decay_curve"]}
            if exp.lifetime.stopper is not None and exp.lifetime.stopper.kind == "differential":
                methods["ddcm_three_foil"] = ANALYSES["ddcm_three_foil"]
            ui.select(methods, value=choice["method"], label="Analysis",
                      on_change=lambda e: rechoose(ctx, method=e.value)).props("dense outlined").classes("w-64")
    res = analysis_of(ctx, p)
    if not res.get("available"):
        _note(ui, res.get("reason") or "No lifetime analysis.")
        if plunger(exp):
            taken = taken_distances(p)
            ui.label(f"Runs at {len(taken)} of the plunger's {len(exp.lifetime.distances_um)} distances so far.").classes(
                "text-sm ps-muted")
        return
    if res.get("kept"):
        n = len(res.get("runs") or [])
        with ui.row().classes("w-full items-center gap-3 ps-note rounded p-2"):
            ui.icon("history").classes("ps-muted")
            ui.label("Kept from an earlier session"
                     + (f" ({res['kept']})" if res["kept"] != "earlier" else "")
                     + f": the analysis of {res.get('level') or 'the studied level'}"
                     + (f" by {ANALYSES.get(res.get('method'), res.get('method'))}" if res.get("method") else "")
                     + (f", of {n} run{'s' * (n != 1)}" if n else "")
                     + ". It is shown as it was left, and was not run again.").classes("text-sm")
            ui.space()
            ui.button("Analyse again", icon="refresh", on_click=lambda: rechoose(ctx, session=True)).props(
                "dense outline no-caps").classes("px-2")
    render_analysis(ctx, res)


def rechoose(ctx, **choice) -> None:
    """A level or method chosen: the analysis is worked out on the worker, then the step is rebuilt."""
    from nicegui import run

    analysis_choice(ctx).update(choice)

    async def go() -> None:
        await run.io_bound(analysis_of, ctx, ctx.P().snapshot())
        ctx.analysis_view.refresh()
        sync = getattr(ctx, "sync_steps", None)
        if sync is not None:  # the rail's status line shows the τ that is kept
            sync()

    ctx.later(go)


def lifetime_comparisons(ctx, res: dict) -> tuple:
    """τ of the analysed level and the B(E2) that follows from it, each beside what is known about the nucleus
    in short (:func:`.tab_analysis.short_comparison`): ([the comparisons], the sources of the analysis's own
    comparison table that the τ comparison says too: "level scheme" when it closes the result against the value
    the simulation was given, "ENSDF" when it has the evaluated value as its reference). Nothing for an analysis
    that does not name its line, or of a nucleus without rows."""
    from .. import results
    from . import tab_analysis as ta

    line = res.get("transition")
    if not line:
        return [], set()
    p = ctx.P()
    try:
        nuclide = results.nuclide_name(line["nuclide"])
        upper, lower = (results.level_key(nuclide, x["energy_kev"], x.get("spin"), x.get("parity"))
                        for x in (line["initial"], line["final"]))
    except Exception:  # noqa: BLE001 -- the result stands without its comparison
        return [], set()
    here = p.experiment_folder.name if p.experiment_folder is not None else None
    held = any(str(c.get("source", "")).startswith("level scheme") for c in res.get("compare") or [])
    missing = ("Not in the results ledger yet: analyse again to enter it." if held else
               "No row in the results ledger: the setup's level scheme holds no lifetime for this level, so there "
               "is no value the simulation was given to close the result against.")
    out, said = [], set()
    view = ta.comparison_of(p, "tau", upper.text, nuclide)
    if view is not None:
        mine = ta.own_projection(view, here, res.get("tau_ps"), analysis="lifetime")
        out.append(ta.short_comparison(view, mine, missing))
        if mine is not None:
            said.add("level scheme")
        if view["comparison"] is not None and view["comparison"].reference is not None:
            said.add("ENSDF")
    b = res.get("b_down")
    if view is not None and res.get("b_multipole", "E2") == "E2" and b and b[0] is not None:
        try:
            key = results.transition_key(upper, lower, "E2")
        except ValueError:
            key = None
        b_view = ta.comparison_of(p, "B(E2)", key.text, nuclide) if key is not None else None
        if b_view is not None:
            mine = ta.own_projection(b_view, here, b[0], analysis="lifetime")
            out.append(ta.short_comparison(b_view, mine, missing))
    return out, said


def render_analysis(ctx, res: dict) -> None:
    """A lifetime analysis (the frozen shape of ``Planner.analyse_lifetimes``): τ with its budget and B(Eλ),
    τ and the B(E2) that follows from it beside what is known about the nucleus (:func:`lifetime_comparisons`),
    the decay curve and τ(x) (or the line shape), and what the analysis's own comparison with τ from elsewhere
    holds beyond that."""
    from . import tab_analysis as ta
    from .figures import figure_decay_curve, figure_lineshape, figure_tau_distance

    ui = ctx.ui
    budget = {k: float(v) for k, v in (res.get("budget") or {}).items() if v is not None}
    total = math.sqrt(sum(v * v for v in budget.values())) if budget else res.get("stat_ps")
    b = res.get("b_down")
    b_text = "—"
    if b is not None:
        value, unc = (b[0], b[1]) if isinstance(b, (list, tuple)) else (b, None)
        if value is not None:
            b_text = f"{_fmt(value)}" + (f" ± {_fmt(unc, 2)}" if unc else "") + f" {res.get('b_unit') or ''}"
    _readouts(ui, [
        ("τ", _tau(res.get("tau_ps"), res.get("stat_ps")), f"{res.get('level') or ''}, "
         f"{ANALYSES.get(res.get('method'), res.get('method') or '')}; statistical", ""),
        ("Total uncertainty", _tau(total, like=res.get("tau_ps")) if total is not None else "—",
         "all of the budget in quadrature", ""),
        ("B(Eλ)↓", b_text.strip(), "from τ, the branching and the conversion", "")])
    try:
        beside, said = lifetime_comparisons(ctx, res)
    except Exception:  # noqa: BLE001 -- the analysis stands without its comparison
        beside, said = [], set()
    if beside:
        with ui.element("div").classes("w-full grid gap-3").style(
                "grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); align-items: start"):
            for short in beside:
                ta.comparison_block(ctx, short)
    # The analysis's own table of τ from elsewhere: only the lines the comparison above does not say already.
    elsewhere = [c for c in res.get("compare") or []
                 if not any(str(c.get("source", "")).startswith(x) for x in said)]
    if res.get("method") == "lineshape":
        ui.plotly(ctx.in_theme(figure_lineshape(res))).classes("w-full")
    else:
        with ui.element("div").classes("w-full grid gap-3").style(
                "grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); align-items: start"):
            with ui.column().classes("gap-0 min-w-0"):
                ui.label("Decay curve").classes("text-xs ps-muted")
                ui.plotly(ctx.in_theme(figure_decay_curve(res))).classes("w-full")
            if any(pt.get("tau_ps") is not None for pt in res.get("points") or []):
                with ui.column().classes("gap-0 min-w-0"):
                    ui.label("τ at each distance (DDCM)").classes("text-xs ps-muted")
                    ui.plotly(ctx.in_theme(figure_tau_distance(res))).classes("w-full")
    tables = ui.element("div").classes("w-full grid gap-3").style(
        "grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); align-items: start")
    with tables:
        if budget:
            with ui.column().classes("gap-0 min-w-0"):
                ui.label("Uncertainty budget").classes("ps-section")
                ui.table(columns=ctx.columns((("k", "Source"), ("v", "τ (ps)"), ("s", "Share"))),
                         rows=[{"k": k.replace("_", " "), "v": f"{v:.3g}",
                                "s": f"{100 * v * v / (total * total):.0f}%" if total else "—"}
                               for k, v in sorted(budget.items(), key=lambda kv: -abs(kv[1]))]).props("dense flat")
        if elsewhere:
            with ui.column().classes("gap-0 min-w-0"):
                ui.label("Beside τ from elsewhere").classes("ps-section")
                rows = []
                for c in elsewhere:
                    bd = c.get("b_down")
                    if isinstance(bd, (list, tuple)):
                        bd = f"{_fmt(bd[0])}" + (f" ± {_fmt(bd[1], 2)}" if len(bd) > 1 and bd[1] else "") \
                            if bd and bd[0] is not None else "—"
                    rows.append({"s": c.get("source", "?"), "t": _tau(c.get("tau_ps")),
                                 "b": bd if isinstance(bd, str) else (_fmt(bd) if bd is not None else "—")})
                ui.table(columns=ctx.columns((("s", "Source"), ("t", "τ"), ("b", f"B(Eλ)↓ ({res.get('b_unit') or ''})"))),
                         rows=rows).props("dense flat")
    for note in _texts(res.get("notes")):
        ui.label(note).classes("text-xs ps-muted")


__all__ = ["analysis_panel", "call", "distance_controls", "distance_spectra", "fusion_data", "fusion_live",
           "plan_panel", "prepare_analysis", "prepare_data", "render_analysis", "render_plan", "series_block",
           "spectra_block", "taken_distances"]
