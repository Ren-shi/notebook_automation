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
        ui.label("Take a beam run on the Run tab: the Plan's numbers are kept with it, and set here beside what it "
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
    note += f" (run {c['source_run']})." if c["source_run"] else ": take one on the Run tab to measure it."
    ui.label(note).classes("text-xs ps-muted")


def render(ctx) -> None:
    ctx.analysis_view()


def render_view(ctx) -> None:
    """Everything on the tab (rebuilt when the current run or the setup changes)."""
    compare_block(ctx)
    exp = ctx.P()._data_experiment()
    if exp.excitation is None or not exp.gamma_detectors:
        ctx.ui.markdown("The analysis needs Coulomb excitation with γ-ray detectors: set the reaction in the "
                        "**Target and reaction** card and add γ-ray detectors.")
    else:
        ctx.analysis_block()
        ctx.alignment_block()
    ctx.multistep_block()
