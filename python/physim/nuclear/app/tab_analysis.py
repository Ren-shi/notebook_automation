"""The Analysis tab (backlog item 68): predicted against measured first, then the analysis of the current run with
a gate from the Data tab, the alignment check and the multi-step solution, all reading the run."""

from __future__ import annotations

import math
from typing import Optional


def _value(x, unit: str, time_text) -> str:
    if x is None:
        return "—"
    if unit == "s":
        return time_text(x)
    if unit == "":
        return f"{100 * x:.3g}%"
    if isinstance(x, float) and not math.isfinite(x):
        return "—"
    from ..logbook import number

    return f"{number(float(x))} {unit}".strip()


def compare_block(ctx) -> None:
    """The Plan's numbers kept in the run's summary, beside what the run measured."""
    ui = ctx.ui
    p = ctx.P()
    c = p.compare()
    ui.label("Predicted against measured").classes("text-lg font-semibold")
    if c["run"] is None:
        ui.label("Take a beam run on the Run step: the Plan's numbers are kept with it, and set here beside what it "
                 "measured.").classes("text-sm ps-muted")
        return
    if not c["rows"]:
        ui.label(f"Run {c['run']} has no predictions kept with it.").classes("text-sm ps-muted")
        return
    rows = []
    for r in c["rows"]:
        mark = {True: "✓ agree", False: "✗ differ", None: ""}[r["agree"]]
        unc = r["uncertainty"]
        measured = _value(r["measured"], r["unit"], ctx.time_text)
        if unc is not None and r["measured"] is not None:
            measured += f" ± {_value(unc, r['unit'], ctx.time_text)}"
        rows.append({"q": r["quantity"], "w": r["what"], "p": _value(r["predicted"], r["unit"], ctx.time_text),
                     "m": measured, "a": mark})
    ui.table(columns=ctx.columns((("q", "Quantity"), ("w", ""), ("p", "Predicted"), ("m", "Measured in the run"),
                                  ("a", ""))), rows=rows).props("dense flat").classes("w-full")
    note = (f"Run {c['run']}, its real part; agreement within three standard deviations of the run's statistics. "
            "The efficiency is measured by a source run")
    note += f" (run {c['source_run']})." if c["source_run"] else ": take one on the Run step to measure it."
    ui.label(note).classes("text-xs ps-muted")


# -- the result beside what is known ------------------------------------------------------------------------------

def comparison_of(p, observable: str, key, nuclide: Optional[str] = None) -> Optional[dict]:
    """One observable of one level or transition of the nucleus the experiment studies, by every route
    (:func:`physim.nuclear.nucleus.observable_view`): the rows of the experiments beside this one and of the
    nucleus's own ledger. None for an experiment that studies no nucleus, or where the view cannot be made."""
    from .. import nucleus, results

    try:
        nuclide = nuclide or results.studied_nuclide(p.experiment)
        if nuclide is None:
            return None
        folder = p.experiment_folder
        return nucleus.observable_view(nuclide, observable, key, root=folder.parent if folder is not None else None)
    except Exception:  # noqa: BLE001 -- the result stands without its comparison
        return None


def same_value(a, b) -> bool:
    """Whether two values are the same one, as written and read back."""
    try:
        return abs(float(a) - float(b)) <= 1e-6 * abs(float(b))
    except (TypeError, ValueError):
        return False


def own_projection(view: Optional[dict], here: Optional[str], value: Optional[float], run: Optional[str] = None,
                   analysis: Optional[str] = None):
    """This experiment's own projected row in a comparison (a :class:`physim.nuclear.compare.Projection`) of the
    result shown, which is ``value``: the one of the run named ``run`` (a Coulex analysis: the B(E2) from the
    yield, not a fit of the same run to all orders), or of the analysis named ``analysis`` ("lifetime"); the newest
    if there are several. ``here`` is the experiment's folder name. None without one, and when the row the ledger
    holds is of another analysis than the one shown."""
    c = view.get("comparison") if view else None
    if c is None or here is None:
        return None
    found = []
    for p in c.projected:
        origin = p.row.origin  # (a value converted from another row is no run's: its origin is the conversion)
        if origin.get("kind") != "run" or origin.get("experiment") != here:
            continue
        if (run is not None and origin.get("run") == run and p.method == "coulex-yield") or \
                (analysis is not None and origin.get("analysis") == analysis):
            found.append(p)
    mine = max(found, key=lambda p: p.row.created) if found else None
    return mine if mine is not None and same_value(mine.row.value, value) else None


def short_comparison(view: dict, mine=None, missing: str = "") -> dict:
    """A result beside what is known, in short and in words, from :func:`comparison_of`.

    ``mine`` is this experiment's projected row in the comparison (:func:`own_projection`). The two kinds of
    evidence stay apart: a projected row is closed against the value its simulation was given and its
    uncertainty is set beside how well the value is known; it is never pulled against the reference. Without
    ``mine`` (the result has no row in the results ledger: ``missing`` says why) the reference and the routes
    stand alone.

    Returns {"title" (what the result is a value of, and of which nucleus), "lines": [(label, text)], "link":
    {"text", "href"} (the nucleus's page at this observable)}. The lines: "Projected" (the result in
    the comparison's convention, and plainly that it comes from a simulated run), "Reference" (the value compared against, with its source), "Closure", "Reach"
    ("this setup reaches ±x %; the value is known to ±y %; the best single measurement to ±z %") and "Other
    routes" (how many hold a value, how many are still open). The words of the closure and the reach are the
    report's (:mod:`physim.nuclear.logbook`).
    """
    from .. import logbook, results
    from . import tab_nucleus

    c = view["comparison"]
    nuclide = logbook.nuclide_words(view["nuclide"])
    what = logbook.quantity_words(view["observable"], view["key"])
    lines = []
    if mine is not None:
        several = len(mine.row.origin.get("runs") or []) > 1
        given = f" (given as {logbook.row_value(mine.row, direction=True)})" if mine.converted else ""
        lines.append(("Projected", f"{logbook.row_value(mine.row, canonical=True, direction=True)}{given}, from "
                      + ("simulated runs: the runs were" if several else "a simulated run: the run was")
                      + f" simulated and analysed back, so this is the precision the setup would reach, not a "
                        f"measurement of {nuclide}."))
    elif missing:
        lines.append(("This result", missing))
    ref = c.reference if c is not None else None
    if ref is not None:
        given = f"; given as {logbook.row_value(ref.row, direction=True)}" if ref.converted else ""
        lines.append(("Reference", f"{logbook.row_value(ref.row, canonical=True, direction=True)}, "
                                   f"{results.EVIDENCE[ref.evidence].lower()} ({view['reference']['source']}{given})."))
    else:
        others = len(c.entries) if c is not None else 0
        lines.append(("Reference", f"None: no evaluated value of {what} is entered for {nuclide}"
                      + (f" ({others} independent value{'s' * (others != 1)} without a reference to compare with)"
                         if others else "") + "."))
    if mine is not None:
        lines.append(("Closure", logbook.closure_words(mine)))
        lines.append(("Reach", logbook.known_words(mine, c)))
    lines.append(("Other routes", _routes_words(view["routes"], mine.method if mine is not None else None)))
    return {"title": f"{what} of {nuclide}", "lines": lines,
            "link": {"text": f"Compare on the {nuclide} page",
                     "href": tab_nucleus.url(view["nuclide"], view["observable"], view["key"])}}


def _routes_words(routes: list, own: Optional[str] = None) -> str:
    """How many of the routes to an observable (:func:`physim.nuclear.nucleus.routes`) hold a value and how many
    are still open, in a sentence; ``own`` is the method of the result shown, which is left out. The reference
    itself is not counted as a value a route holds (the route's "is_reference": the evaluation, usually): a route
    holds a value when it has an independent row besides the reference."""
    others = [r for r in routes if r["method"] != own]

    def named(found: list) -> str:
        return ", ".join(r["label"] for r in found)

    hold = [r for r in others if r.get("others", r["independent"])]
    simulated = [r for r in others if r["projected"] and not r["independent"]]
    open_ = [r for r in others if r["open"]]
    parts = [f"{len(hold)} hold{'s' * (len(hold) == 1)} a value ({named(hold)})" if hold else "none holds a value yet"]
    if simulated:
        parts.append(f"{len(simulated)} {'has' if len(simulated) == 1 else 'have'} simulated runs only "
                     f"({named(simulated)})")
    parts.append(f"{len(open_)} {'is' if len(open_) == 1 else 'are'} still open ({named(open_)})" if open_
                 else "none is still open in physim")
    return ("Of the other routes, " if own is not None else "Of the routes, ") + "; ".join(parts) + "."


def comparison_block(ctx, short: Optional[dict]) -> None:
    """The short comparison under a result (:func:`short_comparison`): a panel with its lines and the link to
    the nucleus's page. Nothing without one."""
    if short is None:
        return
    ui = ctx.ui
    with ui.element("div").classes("ps-panel ps-beside-panel"):
        with ui.element("div").classes("ps-panel-head"):
            ui.label("Beside what is known")
            ui.label("· " + short["title"]).classes("ps-muted")
        with ui.element("div").classes("ps-panel-body"):
            with ui.element("div").classes("ps-beside"):
                for label, text in short["lines"]:
                    if label == "Projected":
                        with ui.element("div"):
                            ui.label("PROJECTED").classes("ps-projected")
                    else:
                        ui.label(label).classes("ps-muted")
                    ui.label(text)
            ui.link(short["link"]["text"] + " →", short["link"]["href"]).classes("ps-nucleus-link")


def coulex_comparison(ctx, result=None) -> None:
    """Under the B(E2) of the current run: the comparison in short for the transition of the ``[reaction]``."""
    from .. import logbook, results

    p = ctx.P()
    try:
        exp = p._data_experiment()
        if exp.excitation is None or exp.excitation.multipolarity != "E2":
            return
        key = results.key_of_reaction(exp)
    except Exception:  # noqa: BLE001 -- the result stands without its comparison
        return
    view = comparison_of(p, "B(E2)", key.text, key.nuclide)
    if view is None:
        return
    run = p._beam_run()
    here = p.experiment_folder.name if p.experiment_folder is not None else None
    mine = own_projection(view, here, getattr(result, "b_e2fm4", None), run=run.label) if run is not None else None
    missing = ""
    if run is None:
        missing = ("Analysed on a simulated sample without a run: it is not kept and has no row in the results "
                   "ledger. Take a run and analyse it to set it beside the other routes.")
    elif mine is None:
        why = ""
        if result is not None:
            try:
                summary = next((s for s in p.runs() if s["label"] == run.label), None)
                why = logbook.analysis_row(p, summary, result)[1]
            except Exception:  # noqa: BLE001
                why = ""
        missing = (f"No row in the results ledger: {why}." if why else
                   "Not in the results ledger yet: analyse the run again to enter it.")
    comparison_block(ctx, short_comparison(view, mine, missing))


def render(ctx) -> None:
    ctx.analysis_view()


#: The parts of the step, in order: (class of its anchor, name in the list at the top).
PARTS = (("lifetimes", "Lifetimes"), ("compare", "Predicted against measured"), ("b-e2", "B(E2) from the peak"),
         ("alignment", "Target alignment"), ("all-orders", "All orders"))


def _jump(ctx, anchor: str) -> None:
    ctx.ui.run_javascript("(() => { const el = document.querySelector('.analysis-%s'); "
                          "if (el) el.scrollIntoView({behavior: 'smooth', block: 'start'}); })()" % anchor)


def not_analysed_note(ctx, gap: dict) -> None:
    """At the top of the step for a kind of experiment whose runs physim does not analyse back yet
    (:func:`physim.nuclear.workbench.not_analysed`): what an analysis would need, and where what the experiment
    does have is."""
    ui = ctx.ui
    with ui.column().classes("w-full gap-2 ps-note rounded p-3"):
        ui.label(f"No analysis for {gap['what']} yet").classes("font-semibold")
        ui.label(gap["line"][0].upper() + gap["line"][1:]).classes("text-sm")
        ui.label(f"What this experiment has: {gap['has']} on the Plan step, its runs beside the plan below, and "
                 "the report, which holds them all.").classes("text-sm")
        if gap["kind"] == "fusion-evaporation":
            ui.label("A lifetime measurement is analysed: add a plunger or line shapes (DSAM) from \"Measure the "
                     "lifetimes\" on the Plan step, and this step gets τ and B(E2) back from its runs.").classes(
                "text-sm ps-muted")
        if gap["kind"] == "elastic":
            ui.label("For a B(E2), make the reaction a Coulomb excitation in the Target and reaction card and add "
                     "γ-ray detectors: a run of that kind is analysed.").classes("text-sm ps-muted")
        with ui.row().classes("gap-2"):
            ui.button("The plan", on_click=lambda: ctx.go_to("plan")).props("dense outline no-caps").classes("px-2")
            ui.button("The report", on_click=lambda: ctx.go_to("report")).props("dense outline no-caps").classes(
                "px-2")


def render_view(ctx) -> None:
    """Everything on the step (rebuilt when the current run or the setup changes): a list of its parts, then
    each part in turn."""
    from .. import workbench as wb

    ui = ctx.ui
    p = ctx.P()
    exp = p._data_experiment()
    coulex = exp.excitation is not None and bool(exp.gamma_detectors)
    lifetime = p.experiment.lifetime is not None
    # Fusion-evaporation and elastic scattering: no analysis of a run yet. The step says so, and what it has. (A
    # current run taken while the setup was a Coulomb excitation is still analysed.)
    gap = wb.not_analysed(p.experiment) if not coulex else None
    if gap is not None:
        not_analysed_note(ctx, gap)
        with ui.column().classes("w-full gap-2 analysis-compare ps-panel").style("padding: 10px 14px 12px"):
            compare_block(ctx)
        if p.experiment.levels and p.experiment.fusion is None:  # solved from the level scheme, without a run
            with ui.column().classes("w-full gap-2 analysis-all-orders ps-panel").style("padding: 10px 14px 12px"):
                ctx.multistep_block()
        return
    shown = [(a, t) for a, t in PARTS if (a == "lifetimes" and lifetime)
             or (a != "lifetimes" and (coulex or a in ("compare", "all-orders")))]
    with ui.row().classes("w-full items-center gap-1"):
        ui.label("On this step").classes("text-xs ps-muted")
        for anchor, title in shown:
            ui.button(title, on_click=lambda a=anchor: _jump(ctx, a)).props("dense outline no-caps").classes("px-2")
    if p.run is None:
        ui.label("No run yet: the analysis and the alignment check below work on a fresh simulated sample of the "
                 "setup. Take a run to analyse what it recorded; the result is then kept with the run and goes "
                 "into the report.").classes("text-sm ps-note rounded p-2")
    elif p.run.kind == "source":
        ui.label(f"Run {p.run.number} is a source run: it measures the γ-ray efficiency. Choose a beam run on the "
                 "Run step to analyse it.").classes("text-sm ps-note rounded p-2")
    if lifetime:
        from .tab_lifetimes import analysis_panel

        with ui.column().classes("w-full gap-2 analysis-lifetimes ps-panel").style("padding: 10px 14px 12px"):
            try:
                analysis_panel(ctx)
            except Exception as err:  # noqa: BLE001 -- the other parts still stand
                ui.label(f"The lifetime analysis could not be shown: {err}").classes("text-sm ps-warn")
    # Side by side where the window has the width: the comparison and the choices of the analysis, then the
    # result across the page, then the alignment check and the all-orders solution.
    grid = "grid-template-columns: repeat(auto-fit, minmax(480px, 1fr)); align-items: start"
    pad = "padding: 10px 14px 12px"
    with ui.column().classes("w-full gap-3"):
        # (The comparison's table needs 640 px: under twice that, the two stand one over the other.)
        top = ui.element("div").classes("w-full grid gap-3").style(grid.replace("480px", "640px"))
        results = ui.column().classes("w-full gap-2")
    with top:
        with ui.column().classes("w-full gap-2 analysis-compare ps-panel").style(pad):
            compare_block(ctx)
        if coulex:
            with ui.column().classes("w-full gap-2 analysis-b-e2 ps-panel").style(pad):
                ctx.analysis_block(results)
    if not coulex and not lifetime:
        ui.markdown("The analysis needs Coulomb excitation with γ-ray detectors: set the reaction in the "
                    "**Target and reaction** card and add γ-ray detectors.")
    with ui.element("div").classes("w-full grid gap-3").style(grid):
        if coulex:
            with ui.column().classes("w-full gap-2 analysis-alignment ps-panel").style(pad):
                ctx.alignment_block()
        with ui.column().classes("w-full gap-2 analysis-all-orders ps-panel").style(pad):
            ctx.multistep_block()
