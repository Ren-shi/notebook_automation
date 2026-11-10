"""The Analysis tab (backlog item 68): predicted against measured first, then the analysis of the current run with
a gate from the Data tab, the alignment check and the multi-step solution, all reading the run."""

from __future__ import annotations

import math


def _value(x, unit: str, time_text) -> str:
    if x is None:
        return "—"
    if unit == "s":
        return time_text(x)
    if unit == "":
        return f"{100 * x:.3g}%"
    if isinstance(x, float) and not math.isfinite(x):
        return "—"
    return f"{x:,.4g} {unit}".strip() if abs(x) >= 1000 else f"{x:.4g} {unit}".strip()


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
                                  ("a", ""))), rows=rows).props("dense flat")
    note = (f"Run {c['run']}, its real part; agreement within three standard deviations of the run's statistics. "
            "The efficiency is measured by a source run")
    note += f" (run {c['source_run']})." if c["source_run"] else ": take one on the Run step to measure it."
    ui.label(note).classes("text-xs ps-muted")


def render(ctx) -> None:
    ctx.analysis_view()


#: The parts of the step, in order: (class of its anchor, name in the list at the top).
PARTS = (("compare", "Predicted against measured"), ("b-e2", "B(E2) from the peak"),
         ("alignment", "Target alignment"), ("all-orders", "All orders"))


def _jump(ctx, anchor: str) -> None:
    ctx.ui.run_javascript("(() => { const el = document.querySelector('.analysis-%s'); "
                          "if (el) el.scrollIntoView({behavior: 'smooth', block: 'start'}); })()" % anchor)


def render_view(ctx) -> None:
    """Everything on the step (rebuilt when the current run or the setup changes): a list of its parts, then
    each part in turn."""
    ui = ctx.ui
    p = ctx.P()
    exp = p._data_experiment()
    coulex = exp.excitation is not None and bool(exp.gamma_detectors)
    shown = [(a, t) for a, t in PARTS if coulex or a in ("compare", "all-orders")]
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
    with ui.column().classes("w-full gap-2 analysis-compare"):
        compare_block(ctx)
    if not coulex:
        ui.markdown("The analysis needs Coulomb excitation with γ-ray detectors: set the reaction in the "
                    "**Target and reaction** card and add γ-ray detectors.")
    else:
        with ui.column().classes("w-full gap-2 analysis-b-e2"):
            ctx.analysis_block()
        with ui.column().classes("w-full gap-2 analysis-alignment"):
            ctx.alignment_block()
    with ui.column().classes("w-full gap-2 analysis-all-orders"):
        ctx.multistep_block()
