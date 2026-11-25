"""One observable of one nucleus, every value side by side: the comparison with the literature.

The view of ``/?nucleus=62Ni&observable=...&key=...`` (and of a row of the nucleus's table, :mod:`.tab_nucleus`):
every value of the observable in the two groups that are never mixed (independent evidence against the reference;
simulated runs apart, each against the value its simulation was given); the plot of the values; the routes to the
observable with "Plan it" where physim has the kind of experiment; what the value gives; the form that takes a data
point; and the history of the rows.

Nothing is computed here: every number comes from :mod:`physim.nuclear.nucleus` through :func:`.nucleus_parts.load`.
"""

from __future__ import annotations

from typing import Optional

from .. import nucleus
from ..compare import VERDICTS
from ..results import EVIDENCE
from .figures import COLORS, THEMES, themed
from .nucleus_parts import (add_style, as_given, cell_text, detail_words, display_factor, group_row, label_text,
                            later, load, load_observable, method_words, nuclide_text, open_nucleus, panel, percent,
                            plan, point_form, point_of, scroll_to_point_form, shown, signed, stat_sys, table_head,
                            td, unit_text, uses, waiting, words)

#: A verdict of the comparison as the page says it, and its colour.
VERDICT_TEXT = {"agree": ("✓ " + VERDICTS["agree"].lower(), "ps-ok"),
                "tension": ("△ " + VERDICTS["tension"].lower(), "ps-warn"),
                "differ": ("✗ " + VERDICTS["differ"].lower(), "ps-bad")}
#: What a simulated run's closure is called: it is no test of the physics, only of the analysis on its own input.
CLOSURE_TEXT = {"agree": ("✓ closes", "ps-ok"), "tension": ("△ closes poorly", "ps-warn"),
                "differ": ("✗ does not close", "ps-bad")}
#: A route that is no experiment of the planner, in words.
MEANS_TEXT = {"calculation": "a calculation: run it in the Shell model panel of an experiment",
              "evaluation": "the evaluation: taken from ENSDF",
              "literature": "known from the literature: its value comes in as a data point"}
#: The order the independent evidence is listed in.
INDEPENDENT_ORDER = ("evaluated", "measured", "outside program", "theory")
#: How many rows of history the page lists at once.
HISTORY_ROWS = 25
#: How much of the line beside a method the comparison table shows (the whole of it on hover).
DETAIL_SHOWN = 60


# ---------------------------------------------------------------------------------------------------------------
# The plot


def figure_values(view: dict, theme: str = "light", display: Optional[str] = None):
    """The values of one observable of one level or transition (an :func:`~physim.nuclear.nucleus.observable_view`)
    with their uncertainties, one line each: the independent evidence against the reference, whose uncertainty is
    the band; under a rule the simulated runs, open markers, each with the value its simulation was given."""
    import plotly.graph_objects as go

    c, cells = view["comparison"], view["cells"]
    factor, unit = display_factor(view["observable"], display)
    ink, grid = THEMES[theme]["ink"], THEMES[theme]["grid"]
    fig = go.Figure()
    names = []

    def trace(found: list, name: str, colour: str, symbol: str = "circle", truth: bool = False) -> None:
        if not found:
            return
        y = list(range(-len(names), -len(names) - len(found), -1))
        names.extend(cells[e.id]["source"] for e in found)
        widths = [(e.plus, e.minus) if not truth else (e.total, e.total) for e in found]
        fig.add_trace(go.Scatter(
            x=[e.value * factor for e in found], y=y, mode="markers", name=name,
            marker=dict(size=10, color=colour, symbol=symbol, line=dict(width=2, color=colour)),
            error_x=dict(type="data", symmetric=False, array=[(p or 0.0) * factor for p, _ in widths],
                         arrayminus=[(m or 0.0) * factor for _, m in widths], color=colour, thickness=1.6, width=4),
            customdata=[[method_words(cells[e.id]), cell_text(cells[e.id], display)] for e in found],
            hovertemplate="%{customdata[1]}<br>%{customdata[0]}<extra>" + name + "</extra>"))
        if truth:
            fig.add_trace(go.Scatter(
                x=[e.truth * factor for e in found], y=y, mode="markers", name="given to the simulation",
                marker=dict(size=14, color=colour, symbol="line-ns-open", line=dict(width=2, color=colour)),
                hovertemplate="the simulation was given %{x:.5g}<extra></extra>"))

    ref = c.reference
    if ref is not None:
        if ref.plus is not None:
            fig.add_vrect(x0=(ref.value - ref.minus) * factor, x1=(ref.value + ref.plus) * factor,
                          fillcolor=COLORS[4], opacity=0.22, line_width=0, layer="below")
        fig.add_vline(x=ref.value * factor, line=dict(color=COLORS[4], dash="dash", width=1))
    for evidence, colour in zip(INDEPENDENT_ORDER, (ink, COLORS[0], COLORS[3], COLORS[2])):
        trace([e for e in c.entries if e.evidence == evidence], EVIDENCE[evidence].lower(), colour)
    if c.projected and names:
        fig.add_hline(y=0.5 - len(names), line=dict(color=grid, dash="dot", width=1))
    trace(list(c.projected), "projected (simulated runs)", COLORS[1], "circle-open", truth=True)
    fig.update_layout(height=max(230, 120 + 26 * len(names)), margin=dict(l=10, r=14, t=10, b=10),
                      legend=dict(orientation="h", yanchor="bottom", y=1.0, x=0), hovermode="closest")
    fig.update_xaxes(title_text=f"{label_text(view['label'])} ({unit_text(unit)})")
    fig.update_yaxes(tickmode="array", tickvals=list(range(0, -len(names), -1)), ticktext=names, automargin=True,
                     range=[0.6 - len(names), 0.6], showgrid=False)
    return themed(fig, theme)


# ---------------------------------------------------------------------------------------------------------------
# The comparison of one observable


def render(ctx) -> None:
    """The comparison of the observable ``ctx.state["start"]["observable"]`` (of the level or transition "key", of
    the nucleus "nucleus"), into the current container. ``ctx`` is the start page's (:func:`.start.render`). What it
    shows is read on a worker thread (:func:`.nucleus_parts.load`), with the evaluation taken first if it never was;
    a spinner stands in meanwhile."""
    ui = ctx.ui
    st = ctx.state["start"]
    name = st["nucleus"]
    add_style(ctx)
    held: dict = {}
    with ui.column().classes("w-full gap-3").style("max-width: 1240px; margin: 0 auto"):
        with ui.row().classes("w-full items-center gap-1"):
            ui.button(f"← {nuclide_text(name)}", on_click=lambda: open_nucleus(ctx, name)).props(
                "flat no-caps dense").classes("ps-accent").tooltip("The page of the nucleus")
            ui.label("▸ compared with the literature").classes("text-sm ps-muted")
        body = ui.column().classes("w-full gap-3")
    with body:
        waiting(ui, f"Reading what is known about {nuclide_text(name)}…")

    async def reload() -> None:
        """Read everything again and build the view from it, in place."""
        from nicegui import run

        try:
            data = await run.io_bound(load, name, st.get("observable"), st.get("key"), st.get("reference"),
                                      not held)
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty page
            data = {"error": f"The page of {name} could not be read: {err}"}
        if data is None or body.is_deleted:
            return
        held.clear()
        held.update(data)
        fill()

    async def choose(observable: str, key: str, reference: Optional[str] = None) -> None:
        """Another observable, or another reference for the one shown: only that one is read again."""
        from nicegui import run

        notes: list = []
        try:
            found = await run.io_bound(load_observable, name, held["view"], observable, key, reference, notes,
                                       held.get("rows"))
        except Exception as err:  # noqa: BLE001
            ui.notify(str(err), type="negative", multi_line=True)
            return
        if body.is_deleted:
            return
        held.update(observable=found, notes=notes)
        st["reference"] = reference
        fill()

    def fill() -> None:
        ov = held.get("observable")
        if ov is not None:  # what the page shows now is what a reload shows again
            ref = ov["comparison"].reference if ov["comparison"] is not None else None
            st.update(observable=ov["observable"], key=ov["key"],
                      reference=ref.id if ref is not None and st.get("reference") else None)
        body.clear()
        with body:
            _page(ctx, held, lambda: later(ui, reload), lambda *a: later(ui, lambda: choose(*a)), fill)
            if st.pop("to_point_form", False):
                scroll_to_point_form(ui)

    later(ui, reload)


def _page(ctx, held: dict, reload, choose, redraw) -> None:
    """The view from what was read. ``reload()`` reads everything again, ``choose(observable, key, reference)``
    another observable or reference, ``redraw()`` builds the view again from what is held."""
    ui = ctx.ui
    if held.get("error"):
        ui.label(held["error"]).classes("text-sm ps-warnbox rounded p-3 w-full")
        return
    name, view, ov = held["nuclide"], held["view"], held["observable"]
    for note in held.get("notes") or ():
        ui.label(note).classes("text-sm ps-note p-2 w-full")
    _picker(ctx, view, ov, choose)
    if ov is None:
        ui.label(f"No value of {nuclide_text(name)} is here yet. Add a data point below, or plan an experiment for "
                 "it.").classes("text-sm ps-note p-3 w-full")
        point_form(ctx, name, ov, reload, redraw)
        return
    with ui.element("div").classes("ps-two"):
        _comparison(ctx, held, reload, choose, redraw)
        with ui.column().classes("w-full gap-3 min-w-0"):
            _plot(ctx, ov)
            _projected_note(ctx, ov)
    _routes(ctx, name, ov)
    with ui.element("div").classes("ps-two"):
        with ui.column().classes("w-full gap-3 min-w-0"):
            _gives(ctx, ov)
        with ui.column().classes("w-full gap-3 min-w-0"):
            point_form(ctx, name, ov, reload, redraw)
            _history(ctx, ov)


def _picker(ctx, view: dict, ov: Optional[dict], choose) -> None:
    """The observables of the nucleus to compare: a button per headline row (and for the one shown), and the rest,
    which hold the evaluation only, in a list to search."""
    ui = ctx.ui
    on = (ov["observable"], ov["key"]) if ov is not None else None
    listed = [t for t in view["table"] if t["headline"] or (t["observable"], t["key"]) == on]
    rest = [t for t in view["table"] if not (t["headline"] or (t["observable"], t["key"]) == on)]
    if not view["table"]:
        return
    with ui.row().classes("w-full items-center gap-1"):
        for t in listed:
            this = (t["observable"], t["key"]) == on
            ui.button(label_text(t["label"]), color="primary" if this else None,
                      on_click=lambda t=t: choose(t["observable"], t["key"])).props(
                "no-caps dense " + ("unelevated" if this else "outline")).classes("px-2")
        if rest:
            ui.select({f"{t['observable']}|{t['key']}": label_text(t["label"]) for t in rest}, with_input=True,
                      label=f"{len(rest)} more hold the evaluation only",
                      on_change=lambda e: e.value and choose(*e.value.split("|", 1))).props(
                "dense outlined options-dense").style("min-width: 260px")


def _comparison(ctx, held: dict, reload, choose, redraw) -> None:
    """Every value of the chosen row, side by side, in two tables that are never mixed: the independent evidence
    against the reference (difference, pull, verdict), and the simulated runs, each against the value its
    simulation was given."""
    ui = ctx.ui
    st = ctx.state["start"]
    name, ov = held["nuclide"], held["observable"]
    c, cells = ov["comparison"], ov["cells"]
    display = st.get("unit")
    factor, unit = display_factor(ov["observable"], display)
    judged = {x["id"]: x["verdict_text"] for g in ov["groups"] for x in g["entries"]}
    heading = f"{label_text(ov['label'])}"
    head, body = panel(ui, heading, "every value, side by side")
    if c is None:
        with body:
            ui.label("No value yet: add a data point, or analyse a run of an experiment.").classes(
                "text-sm ps-muted p-3")
        return
    ref = c.reference
    with head:
        ui.space()
        ui.label("against").classes("text-xs ps-muted")
        options = {e.id: f"{cells[e.id]['source']} · {cell_text(cells[e.id], display)}" for e in c.entries}
        ui.select(options, value=ref.id if ref is not None else None, label=None if ref is not None else "pick a row",
                  on_change=lambda e: (st.update(reference=e.value), choose(ov["observable"], ov["key"], e.value))
                  ).props("dense outlined options-dense").style("min-width: 220px; max-width: 340px").tooltip(
            "The row every other value is compared with: the evaluated one unless another is picked")
    with body:
        # -- independent: measured, evaluated, theory, outside program, each against the reference
        value_head = f"{c.label} ({unit_text(unit)})"
        with ui.element("div").classes("ps-nt-box"), table_head(
                ui, ("Source", "Method", value_head, "Known to", "Difference", "Pull", "Verdict", ""),
                numeric=(2, 3, 4, 5)):
            group_row(ui, "INDEPENDENT · " + ("against the reference: " + cells[ref.id]["source"] if ref is not None
                                             else "no reference yet: pick a row to compare against"), 8)
            for evidence in INDEPENDENT_ORDER:
                found = [e for e in c.entries if e.evidence == evidence]
                if not found:
                    continue
                group_row(ui, EVIDENCE[evidence].upper(), 8)
                for e in found:
                    cell = cells[e.id]
                    with ui.element("tr"):
                        td(ui, cell["source"], ["the reference" if e.is_reference else "",
                                                 "◦ not independent of the reference"
                                                 if e.independent_of_reference is False else ""])
                        given = as_given(cell)
                        detail = detail_words(cell)
                        short = detail if len(detail) <= DETAIL_SHOWN else detail[:DETAIL_SHOWN].rstrip() + "…"
                        method_cell = td(ui, method_words(cell), [short, f"given as {given}: {words(cell['converted'])}"
                                                                   if given and not cell.get("via") else "",
                                                                   cell["problem"]])
                        if short != detail:  # a literature point's method in words: whole on hover
                            method_cell.tooltip(detail)
                        td(ui, shown(e.value, e.plus, e.minus, factor=factor), classes="ps-nt-num")
                        td(ui, percent(e.uncertainty_percent) if e.uncertainty_percent is not None
                            else "no uncertainty", classes="ps-nt-num",
                            colour="" if e.uncertainty_percent is not None else "ps-muted")
                        if e.is_reference or ref is None:
                            for _ in range(2):
                                td(ui, "—", classes="ps-nt-num", colour="ps-muted")
                            td(ui, "reference" if e.is_reference else "—", colour="ps-muted")
                        else:
                            td(ui, signed(e.difference_percent, 1, " %"), classes="ps-nt-num")
                            td(ui, signed(e.pull) if e.pull is not None else (
                                f"× {e.ratio:.2f}" if e.ratio is not None else "—"), classes="ps-nt-num")
                            text, colour = VERDICT_TEXT.get(e.verdict, (judged[e.id], "ps-muted"))
                            td(ui, text, colour=colour)
                        with ui.element("td"):
                            if nucleus.is_data_point(e.row):
                                _row_menu(ctx, name, e.row, cell, reload, redraw)
        ui.label(f"Pull = the difference over the two uncertainties in quadrature (the side of each that faces the "
                 f"other value). Below {c.agree_below:g}: agree; {c.agree_below:g} to {c.differ_above:g}: tension; "
                 f"above {c.differ_above:g}: differ. A value without an uncertainty gets its ratio and no verdict. "
                 f"◦ not independent of the reference: the two share a source, so their agreement is no test."
                 ).classes("text-xs ps-muted").style("padding: 8px 12px; line-height: 1.5")
        # -- projected: simulated runs analysed back, apart, each against what its simulation was given
        with ui.element("div").classes("ps-nt-box").style("border-top: 1px solid var(--ps-line)"):
            with table_head(ui, ("Simulated run", "Method", f"{c.label} ({unit_text(unit)})", "Given to the simulation",
                             "Closure", "Reaches"), numeric=(2, 3, 5)):
                group_row(ui, "PROJECTED · simulated runs analysed back. Each is closed against the value its "
                           "simulation was given, never against the reference", 6, "ps-open")
                if not c.projected:
                    with ui.element("tr"):
                        td(ui, "No simulated run of this is analysed yet.", colour="ps-muted").props("colspan=6")
                for p in c.projected:
                    cell = cells[p.id]
                    with ui.element("tr"):
                        td(ui, cell["source"])
                        td(ui, method_words(cell), detail_words(cell))
                        value, _, parts = stat_sys(p.value, p.stat, p.sys, p.total, factor).partition(" ± ")
                        split = "(" in parts  # the statistical and systematic uncertainty on a line of their own
                        td(ui, value if split else f"{value} ± {parts}" if parts else value,
                            [f"± {x}" for x in parts.split(" ± ")] if split else "", classes="ps-nt-num ps-nt-keep")
                        td(ui, f"{p.truth * factor:.4g}", classes="ps-nt-num")
                        text, colour = CLOSURE_TEXT.get(p.closure_verdict, ("no uncertainty", "ps-muted"))
                        td(ui, f"{text} ({signed(p.closure)}σ)" if p.closure is not None else text, colour=colour)
                        td(ui, percent(p.reach_percent), classes="ps-nt-num")


def _row_menu(ctx, name: str, row, cell: dict, reload, redraw) -> None:
    """Correct and retract, on a data point of the ledger. Neither changes the row: each adds one."""
    ui = ctx.ui
    st = ctx.state["start"]

    def correct() -> None:
        st.update(point=point_of(row), point_problems=[])
        client = ui.context.client  # the menu is gone once the page is drawn again
        redraw()
        client.run_javascript("document.querySelector('.ps-point-form')?.scrollIntoView({behavior: 'smooth'})")

    def retract() -> None:
        with ui.dialog() as dialog, ui.card().classes("ps-card").style("min-width: min(480px, 95vw)"):
            ui.label("Retract this data point").classes("text-base font-semibold")
            ui.label(f"{cell['source']} · {method_words(cell)} · {cell_text(cell)}").classes("text-sm ps-muted")
            ui.label("Nothing is deleted: the point leaves what holds now and stays in the history, with the "
                     "reason beside it.").classes("text-xs ps-muted")
            reason = ui.input("Why it is withdrawn").props("dense outlined autofocus").classes("w-full")
            problem = ui.label("").classes("text-sm ps-bad")

            def go() -> None:
                try:
                    nucleus.retract_data_point(name, row.id, reason.value)
                except ValueError as e:
                    problem.set_text(str(e))
                    return
                dialog.close()
                if st.get("reference") == row.id:
                    st["reference"] = None
                ui.notify("The data point is withdrawn; it stays in the history.", type="positive")
                reload()

            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", color=None, on_click=dialog.close).props("flat no-caps").classes("ps-muted")
                ui.button("Retract", on_click=go).props("unelevated no-caps")
        dialog.open()

    with ui.button(icon="more_vert", color=None).props("flat round dense size=sm").classes("ps-muted").tooltip(
            "Correct or retract this data point"):
        with ui.menu():
            ui.menu_item("Correct…", on_click=correct)
            ui.menu_item("Retract…", on_click=retract)


def _plot(ctx, ov: dict) -> None:
    """The values with their uncertainties and the reference band, drawn as the app's other figures are and in the
    page's theme (drawn again when the theme changes)."""
    ui = ctx.ui
    st = ctx.state["start"]
    if ov["comparison"] is None:
        return
    _, body = panel(ui, "The values", "the band is the reference's uncertainty")
    drawn = {"theme": ctx.state.get("theme", "light")}

    @ui.refreshable
    def figure() -> None:
        ui.plotly(themed_json(ov, drawn["theme"], st.get("unit"))).classes("w-full")

    def follow_theme() -> None:
        if ctx.state.get("theme", "light") != drawn["theme"]:
            drawn["theme"] = ctx.state.get("theme", "light")
            figure.refresh()

    with body:
        with ui.element("div").classes("w-full").style("padding: 6px 8px"):
            figure()
        ui.timer(0.5, follow_theme)


def themed_json(ov: dict, theme: str, display: Optional[str]) -> dict:
    """The figure of the values as the browser takes it."""
    return figure_values(ov, theme, display).to_plotly_json()


def _projected_note(ctx, ov: dict) -> None:
    """What the simulated runs say about the setup, beside what is already known: never a test of the value."""
    c = ov["comparison"]
    reaching = [p for p in (c.projected if c is not None else ()) if p.reach_percent is not None]
    if not reaching:
        return
    ui = ctx.ui
    best = min(reaching, key=lambda p: p.reach_percent)
    cell = ov["cells"][best.id]
    known = (f"the value is known to {percent(best.known_percent)} (the reference)" if best.known_percent is not None
             else "there is no reference with an uncertainty to set it beside")
    measured = (f", and to {percent(best.best_measured_percent)} by the best single measurement"
                if best.best_measured_percent is not None else "")
    share = ""
    if best.sys is not None and best.value:
        share = f", of which {percent(100 * best.sys / abs(best.value))[1:]} is systematic"
    with ui.column().classes("w-full ps-note p-3 gap-1"):
        ui.label("What the projected runs say").classes("text-sm font-semibold")
        ui.label(f"This setup reaches {percent(best.reach_percent)}{share}; {known}{measured}. "
                 f"({cell['source']}, {method_words(cell)}: a simulated run analysed back. It says what the setup "
                 f"would reach, not what the value is.)").classes("text-sm").style("line-height: 1.5")


def _routes(ctx, name: str, ov: dict) -> None:
    """The registry's routes to the observable, each with its state for this level or transition, and "Plan it"
    where physim has the kind of experiment: it opens that kind's new-experiment form for this nucleus."""
    ui = ctx.ui
    title = nuclide_text(name)
    ui.label(f"Routes to {label_text(ov['label'])} in {title}").classes("ps-section").style("margin-top: 6px")
    quiet = []
    with ui.element("div").classes("ps-routes"):
        for r in ov["routes"]:
            status, off = r["status"], r.get("applies") is False
            if status == nucleus.NO_EXPERIMENT and not off:
                quiet.append(r)
                continue
            values = f"{r['independent']} value{'s' * (r['independent'] != 1)}"
            runs_ = f"{r['projected']} simulated run{'s' * (r['projected'] != 1)}"
            tag, colour = (("does not apply", "ps-tag-off") if off else {
                nucleus.DONE: (f"done · {values}" + (f", {runs_}" if r["projected"] else ""), "ps-tag-done"),
                nucleus.ONLY_PROJECTED: (f"projected only · {runs_}", "ps-tag-can"),
                nucleus.IN_PHYSIM: ("in physim", "ps-tag-can"),
                nucleus.LATER: ("later", "")}.get(status, (status, "")))
            with ui.element("div").classes("ps-route" + (" ps-route-later" if off or status == nucleus.LATER else "")):
                with ui.row().classes("w-full items-start justify-between gap-x-2 gap-y-1"):
                    ui.label(r["label"] if r["method"] != "other" else (r["kind_label"] or r["label"])).classes(
                        "text-sm font-semibold")
                    ui.label(tag).classes("ps-tag " + colour)
                how = [r["kind_label"] if r["method"] != "other" else "", r["note"]]
                if any(how):
                    ui.label(": ".join(x for x in how if x)).classes("ps-muted")
                if status == nucleus.ONLY_PROJECTED and not off:
                    ui.label("Simulated runs only: nothing is measured by this route yet.").classes("ps-muted")
                if status == nucleus.LATER:
                    ui.label("A kind of experiment physim does not have yet.").classes("ps-muted")
                if status == nucleus.NO_EXPERIMENT:
                    ui.label(MEANS_TEXT.get(r["means"], "")).classes("ps-muted")
                if r.get("reason"):
                    ui.label(r["reason"]).classes("ps-warn" if off else "")
                if r["live"] and r["kind"] and not off:
                    ui.button(f"Plan it for {title}" if r["open"] else f"Plan another for {title}",
                              icon="arrow_forward", on_click=lambda kind=r["kind"]: plan(ctx, name, kind)).props(
                        "flat dense no-caps").classes("self-start").tooltip(
                        f"Opens the new-experiment form for {r['kind_label']}, with {name} filled in")
    if quiet:
        done = [r for r in quiet if r["rows"]]
        words = [f"{r['label']}" + (f" ({r['rows']} value{'s' * (r['rows'] != 1)} here)" if r["rows"] else "")
                 for r in done + [r for r in quiet if not r["rows"]]]
        ui.label("Not experiments of the planner (the evaluation, a calculation, methods known from the literature: "
                 "their values come in as data points): " + " · ".join(words) + ".").classes("text-xs ps-muted")
        for r in quiet:  # the calculation: whether it fits this computer now
            if r["method"] == "shell-model" and r.get("reason"):
                fits = {True: "fits this computer", False: "does not fit this computer now",
                        None: "not tested"}[r.get("applies")]
                ui.label(f"{r['label']}: {fits}. {r['reason']}").classes("text-xs ps-muted ps-route-shell")


def _gives(ctx, ov: dict) -> None:
    """What the reference value gives (the registry's conversions), the checks between observables when there are
    any, and what the value is used for."""
    ui = ctx.ui
    display = ctx.state["start"].get("unit")
    ref = ov["reference"]
    checks = ov.get("checks", []) or []
    used = uses(ov["observable"])
    if ref is None and not checks and not used:
        return
    _, body = panel(ui, "What this value gives",
                     f"from the reference, {cell_text(ref, display)} ({ref['source']})" if ref else "")
    with body.classes("ps-panel-body"):
        with ui.column().classes("w-full gap-2"):
            if ref is None:
                ui.label("No reference value yet: nothing is converted.").classes("text-sm ps-muted")
            for g in ov["gives"]:
                to = g["to"]
                if g.get("key"):
                    try:
                        to = label_text(nucleus.label(to, g["key"]))
                    except ValueError:
                        pass
                with ui.column().classes("gap-0"):
                    if g["value"] is None:
                        ui.label(f"{to}: not worked out").classes("text-sm ps-muted")
                        ui.label(g["missing"]).classes("text-xs ps-muted")
                        continue
                    factor, unit = (display_factor("B(E2)", display) if g["unit"] == "e2fm4" else (1.0, g["unit"]))
                    ui.label(f"{to} = {shown(g['value'], g['plus'], g['minus'], unit, g['direction'], factor)}"
                             ).classes("text-sm ps-num font-medium")
                    ui.label(words(g["how"])).classes("text-xs ps-muted")
            if checks:
                ui.label("Checks against other observables").classes("ps-section").style("margin-top: 4px")
                for check in checks:
                    with ui.column().classes("gap-0"):
                        ui.label(": ".join(str(x) for x in (check.get("label"), check.get("text")) if x)).classes(
                            "text-sm")
                        if check.get("how"):
                            ui.label(str(check["how"])).classes("text-xs ps-muted")
            if used:
                ui.label("What it is used for").classes("ps-section").style("margin-top: 4px")
                for item in used:
                    with ui.column().classes("gap-0"):
                        ui.label(item["text"]).classes("text-sm").style("line-height: 1.5")
                        if item["reference"]:
                            ui.label(item["reference"]).classes("text-xs ps-muted")


def _history(ctx, ov: dict) -> None:
    """Every row of the chosen observable that was ever written, the newest first: rows are never rewritten."""
    ui = ctx.ui
    rows = ov["history"]
    _, body = panel(ui, "History", "every row is kept")
    words = {"superseded": "superseded by a newer row", "retracted": "withdrawn", "retraction": "the retraction"}
    with body:
        if not rows:
            ui.label("No row yet.").classes("text-sm ps-muted p-3")
            return
        with ui.element("table").classes("ps-nt"):
            for h in rows[:HISTORY_ROWS]:
                gone = h["status"] != "holds"
                with ui.element("tr"):
                    td(ui, h["created"][:10], classes="ps-nt-num", colour="ps-muted")
                    td(ui, h["source"], [f"{EVIDENCE[h['evidence']].lower()} · {h['method_label']}"
                                          + (f" · {h['detail']}" if h["detail"] else ""),
                                          words.get(h["status"], "") + (f": {h['reason']}" if h["reason"] else "")],
                        colour="ps-muted" if gone else "")
                    td(ui, "" if h["status"] == "retraction" else h["text"].replace(
                        " " + h["unit"], " " + unit_text(h["unit"])), classes="ps-nt-num",
                        colour="ps-muted" if gone else "")
        if len(rows) > HISTORY_ROWS:
            ui.label(f"… and {len(rows) - HISTORY_ROWS} older rows.").classes("text-xs ps-muted p-2")
