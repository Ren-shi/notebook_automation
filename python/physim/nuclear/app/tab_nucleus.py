"""The page of one nucleus: what it is, what it holds and its experiments.

The view of ``/?nucleus=62Ni`` (and of a card of the nuclei screen, :mod:`.start`): the nucleus with what is known
and still open, its observables as rows, and its experiments. A row of the observables opens the comparison of that
observable (:mod:`.tab_compare`).

Nothing is computed here: every number comes from :mod:`physim.nuclear.nucleus` through :func:`.nucleus_parts.load`.
What the views of a nucleus share is in :mod:`.nucleus_parts`; other modules import its names from here too.
"""

from __future__ import annotations

from .. import observables  # noqa: F401 -- taken from here by the tests
from ..results import METHODS
from .nucleus_parts import (BE2_UNITS, HEADLINES, add_style, as_given, cell_text, experiments_table, import_button,
                            label_text, later, load, method_words, nuclide_text, open_nuclei, open_observable, panel,
                            plan, point_form, projected_text, report_button, scroll_to_point_form, summary_text,
                            table_head, td, unit_text, waiting)
from .nucleus_parts import (  # noqa: F401 -- other modules and the tests take these from here
    add_point, blank_point, display_factor, import_file, import_lines, load_nuclei, load_observable, open_nucleus,
    point_arguments, point_of, report_html, shown, url, uses)
from .tab_compare import figure_values  # noqa: F401 -- taken from here by the tests


# ---------------------------------------------------------------------------------------------------------------
# The page of one nucleus


def render(ctx) -> None:
    """The page of the nucleus ``ctx.state["start"]["nucleus"]``, into the current container. ``ctx`` is the start
    page's (:func:`.start.render`). What it shows is read on a worker thread (:func:`.nucleus_parts.load`), with the
    evaluation taken first if it never was; a spinner stands in meanwhile."""
    ui = ctx.ui
    st = ctx.state["start"]
    name = st["nucleus"]
    add_style(ctx)
    held: dict = {}
    with ui.column().classes("w-full gap-3").style("max-width: 1240px; margin: 0 auto"):
        with ui.row().classes("w-full items-center gap-1"):
            ui.button("Nuclei", on_click=lambda: open_nuclei(ctx)).props("flat no-caps dense").classes("ps-accent")
            ui.label("▸ " + nuclide_text(name)).classes("text-sm ps-muted")
        body = ui.column().classes("w-full gap-3")
    with body:
        waiting(ui, f"Reading what is known about {nuclide_text(name)}…")

    async def reload() -> None:
        """Read everything again and build the page from it, in place."""
        from nicegui import run

        try:
            data = await run.io_bound(lambda: load(name, first=not held, compare=False))
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty page
            data = {"error": f"The page of {name} could not be read: {err}"}
        if data is None or body.is_deleted:
            return
        held.clear()
        held.update(data)
        fill()

    def fill() -> None:
        body.clear()
        with body:
            _page(ctx, held, lambda: later(ui, reload), fill)

    later(ui, reload)


def _page(ctx, held: dict, reload, redraw) -> None:
    """The page from what was read. ``reload()`` reads everything again, ``redraw()`` builds the page again from what
    is held (the unit switch)."""
    ui = ctx.ui
    st = ctx.state["start"]
    if held.get("error"):
        ui.label(held["error"]).classes("text-sm ps-warnbox rounded p-3 w-full")
        return
    name, view = held["nuclide"], held["view"]
    display = st.get("unit")
    title = nuclide_text(name)

    # -- the head: the nucleus, what it holds, and the two things one does from here
    studied = [x for x in view["experiments"] if x["role"] == "studied"]
    with ui.row().classes("w-full items-end gap-3"):
        with ui.column().classes("gap-0"):
            ui.label(title).classes("text-2xl font-semibold")
            ui.label(f"Z = {view['z']}, N = {view['n']} · everything known and projected about it, by method · "
                     f"{len(studied)} experiment{'s' * (len(studied) != 1)}, {view['analysed_runs']} analysed "
                     f"run{'s' * (view['analysed_runs'] != 1)}, {view['points']} data "
                     f"point{'s' * (view['points'] != 1)}").classes("text-sm ps-muted")
        ui.space()
        ui.button("Add a data point", icon="add", on_click=lambda: _to_point_form(ctx, view)).props(
            "outline no-caps").tooltip("The form is with the values of an observable")
        import_button(ctx, name, reload)
        report_button(ctx, name)
        ui.button(f"New experiment for {title}", icon="add", on_click=lambda: plan(ctx, name)).props(
            "unelevated no-caps")
    for note in held.get("notes") or ():
        ui.label(note).classes("text-sm ps-note p-2 w-full")
    evaluation = view["evaluation"]
    if evaluation["taken"]:
        ui.label(f"The evaluation: {evaluation['reference']}, {evaluation['rows']} values, taken "
                 f"{evaluation['taken'][:10]}.").classes("text-xs ps-muted")

    # -- the readouts: the headline rows, and what is still open
    headline = [t for t in view["table"] if t["headline"]]
    with ui.element("div").classes("ps-readouts"):
        for t in headline[:HEADLINES - 1]:
            best = t["evaluated"] or next((t[k]["best"] for k in ("measured", "outside_program", "theory")
                                           if t[k] is not None), None)
            with ui.element("div").classes("ps-readout"):
                ui.label(label_text(t["label"]) + (", evaluated" if t["evaluated"] else "")).classes(
                    "text-xs ps-muted")
                ui.label(cell_text(best, display) if best else "no independent value").classes(
                    "ps-readout-value" + ("" if best else " ps-muted"))
                ui.label(("" if t["evaluated"] or not best else method_words(best) + " · ")
                         + (f"projected: {projected_text(t['projected'])[0]}" if t["projected"] else
                            "no simulated run")).classes("text-xs ps-muted")
        with ui.element("div").classes("ps-readout"):
            ui.label("Still open").classes("text-xs ps-muted")
            ui.label(str(len(view["open"]))).classes("ps-readout-value" + (" ps-warn" if view["open"] else ""))
            kinds = {"route": "routes physim can plan", "runs": "runs to analyse", "evaluation": "the evaluation"}
            ui.label(", ".join(dict.fromkeys(kinds.get(o["what"], o["what"]) for o in view["open"]))
                     or "nothing").classes("text-xs ps-muted")

    _observables(ctx, view, redraw)
    if view["table"]:
        _experiments(ctx, view)
        return
    ui.label(f"No value of {title} is here yet. Add a data point below, or plan an experiment for it."
             ).classes("text-sm ps-note p-3 w-full")
    with ui.element("div").classes("ps-two"):
        with ui.column().classes("w-full gap-3 min-w-0"):
            _experiments(ctx, view)
        with ui.column().classes("w-full gap-3 min-w-0"):
            point_form(ctx, name, None, reload, redraw)


def _to_point_form(ctx, view: dict) -> None:
    """The form of a data point: on this page for a nucleus without a row, otherwise with the values of the first
    observable of its table (the form takes any other too)."""
    if not view["table"]:
        scroll_to_point_form(ctx.ui)
        return
    first = view["table"][0]
    open_observable(ctx, view["nuclide"], first["observable"], first["key"], to_point_form=True)


def _observables(ctx, view: dict, redraw) -> None:
    """The observables of the nucleus, one row per level or transition: evaluated, measured, theory and what the
    simulated runs reach in columns of their own, what is still open at the end. A row opens its comparison
    (:mod:`.tab_compare`), with its routes."""
    ui = ctx.ui
    st = ctx.state["start"]
    display = st.get("unit")
    listed = [t for t in view["table"] if t["headline"]]
    rest = [t for t in view["table"] if not t["headline"]]
    outside = any(t["outside_program"] for t in listed)
    head, body = panel(ui, "Observables", "one row per level or transition; click a row for its routes")
    with head:
        ui.space()
        if any(t["observable"] == "B(E2)" for t in view["table"]):
            def switch(e) -> None:
                st["unit"] = e.value
                redraw()

            ui.toggle({u: f"B(E2)↑ in {unit_text(u)}" for u in BE2_UNITS}, value=display or BE2_UNITS[0],
                      on_change=switch).props("dense no-caps unelevated").tooltip(
                "How a B(E2) is written here. The rows are not changed; a value in W.u. is always the ↓ value.")
    with body:
        if not view["table"]:
            ui.label("No row yet.").classes("text-sm ps-muted p-3")
            return
        columns = ("Observable", "Evaluated", "Measured") + (("Outside program",) if outside else ()) + (
            "Theory", "Projected (simulated runs)", "Still open")
        if listed:
            with table_head(ui, columns):
                for t in listed:
                    row = ui.element("tr").classes("ps-nt-click")
                    row.on("click", lambda t=t: open_observable(ctx, view["nuclide"], t["observable"], t["key"]))
                    with row:
                        td(ui, label_text(t["label"]), colour="font-medium").style("white-space: nowrap")
                        td(ui, cell_text(t["evaluated"], display) if t["evaluated"] else "—",
                            as_given(t["evaluated"]) and "given as " + as_given(t["evaluated"])
                            if t["evaluated"] else "", classes="ps-nt-num")
                        for kind in ("measured",) + (("outside_program",) if outside else ()) + ("theory",):
                            td(ui, *summary_text(t[kind], display), classes="ps-nt-num")
                        td(ui, *projected_text(t["projected"]), classes="ps-nt-num")
                        td(ui, " · ".join(METHODS[m] for m in t["open"]) or "nothing",
                            colour="ps-open" if t["open"] else "ps-muted")
        if rest:
            with ui.expansion(f"{len(rest)} more hold the evaluation only", icon="list").props("dense").classes(
                    "w-full text-sm"):
                rows = [{"id": f"{t['observable']}|{t['key']}", "label": label_text(t["label"]),
                         "value": cell_text(t["evaluated"], display) if t["evaluated"] else "—"} for t in rest]
                more = ui.table(columns=[{"name": "label", "label": "Observable", "field": "label", "align": "left"},
                                         {"name": "value", "label": "Evaluated", "field": "value", "align": "right"}],
                                rows=rows, row_key="id", pagination=12).props("dense flat").classes("w-full")
                more.add_slot("body-cell-label", '<q-td :props="props" style="cursor: pointer; font-family: inherit">'
                                                 '<span class="ps-accent">{{ props.value }}</span></q-td>')
                more.on("rowClick", lambda e: open_observable(ctx, view["nuclide"],
                                                              *e.args[1]["id"].split("|", 1)))


def _experiments(ctx, view: dict) -> None:
    ui = ctx.ui
    found = view["experiments"]
    _, body = panel(ui, "Experiments", f"{len(found)} with this nucleus; a row opens the experiment")
    with body:
        if not found:
            ui.label("No experiment with this nucleus yet.").classes("text-sm ps-muted p-3")
            return
        experiments_table(ctx, found, ctx.state["start"].get("unit"))
