"""The Run tab (backlog item 66): one button, a duration, live counters.

Run for a duration (the setup's beam time by default) as a beam run, a source run or an alignment check; the part
simulated event by event and what it costs; progress in experiment time; live counters per detector and crystal,
the particle × γ matrix filling in and the counts-wanted bars; Stop and Extend; the last events; what is not
simulated; and the list of runs, one of which is current for the Data and Analysis tabs."""

from __future__ import annotations

import math

from .. import response as _response
from .. import runs as _runs
from .. import workbench as wb

BUDGETS = ["1 min", "5 min", "10 min", "30 min", "1 h"]


def _time(seconds) -> str:
    if seconds is None or not math.isfinite(seconds):
        return "—"
    if seconds >= 86400:
        return f"{seconds / 86400:.3g} d"
    if seconds >= 3600:
        return f"{seconds / 3600:.3g} h"
    if seconds >= 60:
        return f"{seconds / 60:.3g} min"
    return f"{seconds:.3g} s"


def _clock(seconds: float) -> str:
    s = int(round(seconds))
    return f"{s // 3600:d}:{(s % 3600) // 60:02d}:{s % 60:02d}"


def _n(x) -> str:
    if x is None:
        return "—"
    return f"{x:,.0f}" if abs(x) >= 100 else f"{x:.3g}"


def describe_runs(planner) -> list:
    """The run list as table rows."""
    rows = []
    for s in planner.runs():
        rows.append({"n": s["number"], "label": s["label"], "what": s["describe"],
                     "when": str(s.get("finished", ""))[:16].replace("T", " "),
                     "stale": "setup changed since" if s["stale"] else ""})
    return rows


def form(ctx) -> None:
    """The kind of run, its duration and real part, and the Run (or Stop) button."""
    ui = ctx.ui
    p = ctx.P()
    f = ctx.run_form
    running = p.running
    with ui.row().classes("items-end gap-2 w-full"):
        ui.select({"beam": "Beam run", "source": "Source run (beam off)", "alignment": "Alignment check"},
                  value=f["kind"], label="Kind of run",
                  on_change=lambda e: (f.update(kind=e.value), ctx.refresh_run_form())).props(
            "dense outlined").classes("w-48")
        ui.input("Run for", value=f["duration"] or str(p.experiment.run.beam_time),
                 on_change=lambda e: f.update(duration=e.value)).props("dense outlined").classes("w-28").on(
            "blur", lambda: ctx.refresh_run_form())
        if f["kind"] != "source":
            ui.select(BUDGETS, value=f["budget"], label="Event by event",
                      on_change=lambda e: (f.update(budget=e.value), ctx.refresh_run_form())).props(
                "dense outlined").classes("w-36")
        if running:
            ui.button("Stop", icon="stop", on_click=ctx.stop_run).props("color=negative unelevated")
        else:
            ui.button("Run", icon="play_arrow", on_click=lambda: ctx.start_run(dict(f))).props("unelevated")
    if f["kind"] == "source":
        with ui.row().classes("items-end gap-2"):
            ui.select(_response.source_names(), value=f["source"], label="Source",
                      on_change=lambda e: f.update(source=e.value)).props("dense outlined").classes("w-32")
            ui.input("Activity", value=f["activity"], on_change=lambda e: f.update(activity=e.value)).props(
                "dense outlined").classes("w-32")
            ui.input("Position (x, y, z mm)", value=f["position"],
                     on_change=lambda e: f.update(position=e.value)).props("dense outlined").classes("w-48")
        ui.label("Beam off: the source's lines in every crystal, with the counting statistics of the duration. The "
                 "position is from the target; empty puts the source at the target.").classes("text-xs ps-muted")
    elif f["kind"] == "alignment":
        ui.number("Target assumed off along the beam (mm)", value=f["offset_mm"], step=0.5,
                  on_change=lambda e: f.update(offset_mm=e.value)).props("dense outlined").classes("w-72")
        ui.label("A beam run analysed as if the target stood elsewhere: the Analysis tab shows what the Doppler "
                 "correction then does, and fits the offset back.").classes("text-xs ps-muted")
    if f["kind"] != "source" and not running:
        try:
            e = p.estimate_run(f["duration"] or None, f["budget"])
        except Exception:  # noqa: BLE001 - an estimate must never break the tab
            return
        scaled = f", then {_time(e['scaled_s'])} scaled from it" if e["scaled_s"] > 0 else ", all of it"
        ui.label(f"Simulates the first {_time(e['real_s'])} event by event{scaled}: about {_n(e['particles'])} "
                 f"particles, about {_time(e['cpu_s'])} of CPU and {e['size_mb']:.3g} MB on disk.").classes(
            "text-sm")


def live(ctx) -> None:
    """The progress and the live counters of the run being taken."""
    ui = ctx.ui
    p = ctx.P()
    g = p.run_progress()
    if g is None:
        return
    with ui.column().classes("w-full ps-plate gap-1"):
        with ui.row().classes("items-baseline gap-3"):
            ui.label(_clock(g["real_s"])).classes("ps-readout-value")
            ui.label(f"of {_clock(g['real_target_s'])} simulated event by event (run {g['number']})").classes(
                "text-sm ps-muted")
        ui.linear_progress(value=g["fraction"], show_value=False).classes("w-full")
        if g["scaled_s"] > 0:
            ui.label(f"Then the remaining {_time(g['scaled_s'])} of the {_time(g['duration_s'])} are scaled from "
                     "the real part, with their own Poisson fluctuations.").classes("text-xs ps-muted")
        if g["stopping"]:
            ui.label("Stopping after this step; what is accumulated is kept.").classes("text-sm ps-warn")
        counters(ctx, g)
    sample = p.run_sample(8)
    if sample:
        ui.label("The last particles counted").classes("ps-section mt-2")
        ui.table(columns=ctx.columns((("detector", "Detector"), ("segment", "Ring/strip, sector"),
                                      ("energy", "Energy (MeV)"), ("channel", "Channel"))),
                 rows=[{"detector": x["detector"], "segment": f"{x['segment'][0]}, {x['segment'][1]}",
                        "energy": f"{x['energy_mev']:.3f}", "channel": x["channel"]} for x in sample]).props(
            "dense flat")


def counters(ctx, g: dict) -> None:
    ui = ctx.ui
    wanted = g["counts_wanted"]
    what = {"all": "counts", "excitations": "excitations", "coincidences": "coincidences"}[g["measured"]]
    rows = []
    for d, c in g["counts"].items():
        rows.append({"d": d, "c": _n(c), "r": _n(g["rates_per_s"][d]), "b": f"{100 * g['busy'][d]:.2g}%",
                     "w": (f"{_n(g['projected'][d])} of {wanted}" if wanted else "—")})
    ui.table(columns=ctx.columns((("d", "Particle detector"), ("c", "Counts"), ("r", "Rate (1/s)"),
                                  ("b", "Busy"), ("w", f"{what.capitalize()} by the end, wanted"))),
             rows=rows).props("dense flat")
    if wanted:
        for d, proj in g["projected"].items():
            with ui.row().classes("items-center gap-2 w-full no-wrap"):
                ui.label(d).classes("text-xs w-24")
                ui.linear_progress(value=min(proj / wanted, 1.0), show_value=False).classes("w-full").props(
                    f"color={'positive' if proj >= wanted else 'warning'}")
    ui.label(f"Live fraction {g['live_fraction']:.3f}.").classes("text-xs ps-muted")
    if g["gamma_singles"]:
        crystals = sorted(g["gamma_singles"])
        ui.table(columns=ctx.columns((("c", "Crystal"), ("s", "Singles (from the rates)"),
                                      ("k", "In coincidence, any energy"))),
                 rows=[{"c": c, "s": _n(g["gamma_singles"][c]), "k": _n(g["gamma"].get(c, 0))}
                       for c in crystals]).props("dense flat").classes("max-h-64")
        gdets = list(dict.fromkeys(k.split("|")[1] for k in g["coincidences"])) or []
        if gdets:
            ui.label("Particle × γ coincidences so far: in the full-energy peak / at any energy").classes(
                "ps-section mt-1")
            ui.table(columns=ctx.columns([("p", "")] + [(f"g{i}", n) for i, n in enumerate(gdets)]),
                     rows=[dict({"p": d}, **{f"g{i}": f"{_n(g['peak_coincidences'].get(f'{d}|{n}', 0))} / "
                                                     f"{_n(g['coincidences'].get(f'{d}|{n}', 0))}"
                                             for i, n in enumerate(gdets)}) for d in g["counts"]]).props(
                "dense flat")


def runs_list(ctx) -> None:
    ui = ctx.ui
    p = ctx.P()
    rows = describe_runs(p)
    ui.label("Runs").classes("ps-section mt-3")
    if not rows:
        ui.label("No run yet.").classes("text-sm ps-muted")
        return
    ui.table(columns=ctx.columns((("n", "#"), ("what", "Run"), ("when", "Finished (UTC)"), ("stale", ""))),
             rows=rows).props("dense flat")
    current = p.run.number if p.run is not None else None
    with ui.row().classes("items-end gap-2 w-full"):
        ui.select({r["n"]: f"{r['n']}: {r['what']}" for r in rows}, value=current, label="Current run",
                  on_change=lambda e: ctx.load_run(e.value)).props("dense outlined").classes("w-96")
        if p.run is not None and p.run.kind != "source" and not p.running:
            more = ui.input("Extend by", value="10 min").props("dense outlined").classes("w-28")
            ui.button("Extend", icon="update", on_click=lambda: ctx.extend_run(more.value)).props("flat no-caps")
            ui.button("Watch events", icon="timeline", on_click=ctx.watch_events).props("flat no-caps").tooltip(
                "A sample of this run's events as tracks in the scene")
    if p.run is not None:
        s = p.run.summary
        ui.label(p.run.describe()).classes("text-sm")
        if s.get("kind") != "source":
            ui.label("Counts: " + ", ".join(f"{d} {_n(c)} ({_n(s['rates_per_s'].get(d))}/s)"
                                            for d, c in s.get("counts", {}).items())).classes("text-xs ps-muted")


def not_simulated(ctx) -> None:
    ui = ctx.ui
    with ui.expansion("Not simulated yet", icon="info").classes("w-full ps-sunk mt-3"):
        for item, what, why in wb.NOT_SIMULATED:
            with ui.row().classes("items-start no-wrap gap-2"):
                ui.label(f"{what}:").classes("text-sm font-medium")
                ui.label(f"{why} (backlog item {item}).").classes("text-sm ps-muted")


def render(ctx) -> None:
    ctx.run_view()


def render_view(ctx) -> None:
    """Everything on the tab (refreshed while a run is taken)."""
    form(ctx)
    live(ctx)
    runs_list(ctx)
    not_simulated(ctx)


__all__ = ["describe_runs", "render", "render_view", "_runs"]
