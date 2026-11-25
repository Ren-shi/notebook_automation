"""One observable of one nucleus, every value side by side: the comparison with the literature.

The view of ``/?nucleus=62Ni&observable=...&key=...`` (and of a card or row of the nucleus's page, :mod:`.tab_nucleus`),
from top to bottom: the way back and the actions; the observables of the nucleus as tabs; a row of five tiles (the
reference, the measured values, the best single measurement, theory, what the setup reaches); one table in which every
value is drawn with its uncertainty on one axis, against the reference's band, grouped by method or ordered by year or
by pull; what stands out and what the projected runs say; and under them the routes to the value, what it gives and is
used for, the checks and the history. The form of a data point is a window.

Nothing is computed here: every number, text, place on the axis and count comes from
:func:`physim.nuclear.nucleus.observable_view` through :func:`.nucleus_parts.load`. The page only turns a place
between 0 and 1 into a position.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from .. import nucleus
from ..results import EVIDENCE
from .figures import COLORS, THEMES, themed
from .nucleus_parts import (add_style, cell_text, display_factor, import_button, label_text, later, load,
                            load_observable, method_words, nuclide_text, open_experiment, open_nucleus, panel, plan,
                            point_form, point_of, report_button, shown, td, unit_text, uses, waiting, words)

#: The mark and the colour of a verdict; the words are the data's own (``verdict_text``).
VERDICT_MARK = {"agree": ("✓", "ps-ok"), "tension": ("△", "ps-warn"), "differ": ("✗", "ps-bad")}
#: A route that is no experiment of the planner, in words.
MEANS_TEXT = {"calculation": "a calculation: run it in the Shell model panel of an experiment",
              "evaluation": "the evaluation: taken from ENSDF",
              "literature": "known from the literature: its value comes in as a data point"}
#: The order the independent evidence is listed in (the figure of the values).
INDEPENDENT_ORDER = ("evaluated", "measured", "outside program", "theory")
#: How many rows of history the page lists at once.
HISTORY_ROWS = 25
#: How much of a row's note the table shows (the whole of it on hover).
DETAIL_SHOWN = 60
#: The orders of the table.
ORDERS = {"method": "By method", "year": "By year", "pull": "By pull"}
#: The key to the marks, in the order it is listed: the family and what it is called there.
KEY_ORDER = (("coulex", "Coulomb excitation"), ("electron-scattering", "Electron scattering"),
             ("lifetime", "From the lifetime"), ("other", "Other measured"), ("outside program", "Outside program"),
             ("evaluation", "Other evaluations"), ("theory", "Theory"), ("projected", "Projected"))
#: What a finding of "What stands out" is coloured by.
FINDING_COLOUR = {"differ": "ps-bad", "tension": "ps-warn", "superseded": "ps-warn"}

#: The screen's own style: the head, the tabs, the tiles, the table with its bars and the boxes, in the page's
#: colours (both themes). The marks of the families have colours of their own, as the figures do.
COMPARE_STYLE = """
.ps-cmp { --ps-m1: #0072B2; --ps-m2: #A8477F; --ps-m3: #8A6A00; --ps-m4: #56616C; --ps-band: rgba(143, 84, 0, .16);
  display: flex; flex-direction: column; gap: 20px; width: 100%; max-width: 1376px; margin: 0 auto;
  font-size: 14px; line-height: 1.45; }
body.body--dark .ps-cmp { --ps-m1: #5CC8F0; --ps-m2: #E28BC0; --ps-m3: #F0D060; --ps-m4: #A7B1BA;
  --ps-band: rgba(242, 163, 58, .18); }
.ps-cmp > * { min-width: 0; }
.ps-cmp-where { display: flex; flex-wrap: wrap; align-items: center; gap: 4px 8px; font-size: 13px;
  color: var(--ps-muted); }
.ps-cmp .q-btn.ps-cmp-back { color: var(--ps-accent) !important; padding: 0 4px; min-height: 36px; font-size: 13px; }
.ps-cmp-head { display: flex; flex-wrap: wrap; align-items: baseline; justify-content: space-between; gap: 8px 20px; }
.ps-cmp-title { margin: 0; font-size: 28px; font-weight: 600; line-height: 1.2; }
.ps-cmp-title > div { display: inline; }
.ps-cmp-sym { font-family: var(--ps-mono); font-weight: 500; }
.ps-cmp-actions { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
.ps-cmp-tabs { display: flex; flex-wrap: wrap; align-items: center; gap: 2px;
  border-bottom: 1px solid var(--ps-line); }
.ps-cmp .q-btn.ps-cmp-tab { min-height: 44px; padding: 0 16px; border-radius: 0; font-family: var(--ps-mono);
  font-weight: 400; color: var(--ps-muted) !important; border-bottom: 2px solid transparent; }
.ps-cmp .q-btn.ps-cmp-tab-on { color: var(--ps-ink) !important; font-weight: 600; border-bottom-color: var(--ps-fill); }
.ps-cmp-more { margin-left: 8px; min-width: 240px; padding: 4px 0; }
.ps-cmp-tiles { display: flex; flex-wrap: wrap; gap: 1px; background: var(--ps-line); border: 1px solid var(--ps-line);
  border-radius: 8px; overflow: hidden; }
.ps-cmp-tile { flex: 1 1 190px; background: var(--ps-surface); padding: 14px 18px; min-width: 0; }
.ps-cmp-tile-ref { flex: 1.4 1 240px; }
.ps-cmp-tile-h, .ps-cmp-tile-s { font-size: 12px; color: var(--ps-muted); }
.ps-cmp-tile-s > div { display: inline; }
.ps-cmp-tile-v { font-family: var(--ps-mono); font-size: 26px; font-weight: 500; line-height: 1.3; }
.ps-cmp-tile-v > div { display: inline; }
.ps-cmp-tile-v .ps-cmp-tile-unit { font-size: 13px; color: var(--ps-muted); font-weight: 400; }
.ps-cmp-tile-ref .ps-cmp-tile-v { color: var(--ps-accent); }
.ps-cmp-tile-v.ps-cmp-none { color: var(--ps-muted) !important; font-size: 20px; line-height: 1.6; }
.ps-cmp-box { background: var(--ps-surface); border: 1px solid var(--ps-line); border-radius: 8px; }
.ps-cmp-controls { display: flex; flex-wrap: wrap; align-items: center; gap: 10px 20px; padding: 10px 18px;
  border-bottom: 1px solid var(--ps-line); }
.ps-cmp-control { display: inline-flex; align-items: center; gap: 8px; font-size: 12px; color: var(--ps-muted); }
.ps-cmp-key { display: flex; flex-wrap: wrap; align-items: center; gap: 6px 16px; font-size: 12px;
  color: var(--ps-muted); margin-left: auto; }
.ps-cmp-key > div { display: inline-flex; align-items: center; gap: 6px; }
.ps-cmp-scroll { overflow-x: auto; }
.ps-cmp-grid { min-width: 1180px; }
.ps-cmp-row { display: grid; grid-template-columns: 236px 214px 132px minmax(280px, 1fr) 72px 64px 112px 28px;
  gap: 12px; align-items: stretch; padding: 0 18px; min-height: 46px; border-bottom: 1px solid var(--ps-hair); }
.ps-cmp-row > div { align-self: center; min-width: 0; padding: 5px 0; }
.ps-cmp-row > .ps-cmp-bar { align-self: stretch; position: relative; overflow: hidden; padding: 0; }
.ps-cmp-row.ps-cmp-colhead { align-items: end; min-height: 0; padding-top: 8px; padding-bottom: 8px; font-size: 12px;
  color: var(--ps-muted); border-bottom-color: var(--ps-line); }
.ps-cmp-row.ps-cmp-colhead > div { padding: 0; }
.ps-cmp-ticks { position: relative; height: 18px; font-family: var(--ps-mono); font-size: 11px; }
.ps-cmp-ticks > div { position: absolute; transform: translateX(-50%); white-space: nowrap; }
.ps-cmp-group { display: flex; flex-wrap: wrap; align-items: baseline; gap: 2px 12px; padding: 9px 18px 7px;
  background: var(--ps-sunk); border-bottom: 1px solid var(--ps-hair); }
.ps-cmp-group-label { font-size: 12px; font-weight: 600; letter-spacing: .04em; }
.ps-cmp-group-note { font-size: 12px; color: var(--ps-muted); }
.ps-cmp-src { font-weight: 500; }
.ps-cmp-sub { font-size: 11px; color: var(--ps-muted); }
.ps-cmp-note { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; cursor: help; }
.ps-cmp-num { font-family: var(--ps-mono); font-size: 13px; font-variant-numeric: tabular-nums; }
.ps-cmp-num + .ps-cmp-sub { font-family: var(--ps-mono); }
.ps-cmp-right { text-align: right; }
.ps-cmp-verdict { font-size: 13px; font-weight: 500; }
.ps-cmp-band { position: absolute; top: 0; bottom: 0; background: var(--ps-band); }
.ps-cmp-refline { position: absolute; top: 0; bottom: 0; width: 1px; background: var(--ps-accent); opacity: .7; }
.ps-cmp-err { position: absolute; top: 50%; height: 2px; margin-top: -1px; background: currentColor; }
.ps-cmp-truth { position: absolute; top: 50%; width: 2px; height: 14px; margin: -7px 0 0 -1px;
  background: var(--ps-accent); }
.ps-mk { position: absolute; top: 50%; width: 11px; height: 11px; margin: -5.5px 0 0 -5.5px; box-sizing: border-box;
  border: 2px solid currentColor; background: currentColor; border-radius: 50%; }
.ps-cmp-key .ps-mk, .ps-cmp-key .ps-cmp-truth { position: relative; top: auto; margin: 0; display: inline-block; }
.ps-cmp-key .ps-cmp-truth { height: 12px; }
.ps-cmp-keyband { width: 12px; height: 12px; background: var(--ps-band); border-left: 1px solid var(--ps-accent); }
.ps-fam-reference { color: var(--ps-accent); }
.ps-fam-coulex { color: var(--ps-m1); }
.ps-fam-electron-scattering { color: var(--ps-m2); }
.ps-mk.ps-fam-electron-scattering { width: 10px; height: 10px; margin: -5px 0 0 -5px; border-radius: 1px;
  transform: rotate(45deg); }
.ps-fam-lifetime { color: var(--ps-m3); }
.ps-mk.ps-fam-lifetime { border-radius: 1px; }
.ps-fam-other { color: var(--ps-ink); }
.ps-mk.ps-fam-other { width: 12px; height: 11px; margin: -6px 0 0 -6px; border: 0; border-radius: 0;
  clip-path: polygon(50% 0, 100% 100%, 0 100%); }
.ps-fam-outside-program { color: var(--ps-m2); }
.ps-mk.ps-fam-outside-program { background: transparent; border-radius: 1px; }
.ps-fam-evaluation { color: var(--ps-m4); }
.ps-mk.ps-fam-evaluation { width: 10px; height: 10px; margin: -5px 0 0 -5px; background: transparent;
  border-radius: 1px; transform: rotate(45deg); }
.ps-fam-theory { color: var(--ps-m4); }
.ps-mk.ps-fam-theory { background: transparent; }
.ps-fam-projected { color: var(--ps-accent); }
.ps-mk.ps-fam-projected { background: transparent; border-style: dashed; }
.ps-cmp-foot { padding: 12px 18px; font-size: 12px; color: var(--ps-muted); max-width: 980px; }
.ps-cmp-boxes { display: flex; flex-wrap: wrap; gap: 20px; align-items: stretch; }
.ps-cmp-card { flex: 1 1 420px; min-width: 0; background: var(--ps-surface); border: 1px solid var(--ps-line);
  border-radius: 8px; padding: 14px 18px; display: flex; flex-direction: column; gap: 8px; }
.ps-cmp-h2 { font-size: 16px; font-weight: 600; }
.ps-cmp .q-btn.ps-cmp-link { color: var(--ps-accent) !important; font-weight: 600; min-height: 40px; padding: 0 2px;
  align-self: flex-start; }
.ps-cmp-lower { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1fr); gap: 20px; align-items: start; }
.ps-cmp-lower > div { display: flex; flex-direction: column; gap: 12px; min-width: 0; }
.ps-cmp-window { padding: 0 !important; position: relative; max-height: 92vh; overflow-y: auto; }
.ps-cmp-window .ps-panel { border: 0; }
.ps-cmp-window .ps-panel-head { padding-right: 52px; }
.ps-cmp-window .ps-cmp-close { position: absolute; top: 4px; right: 6px; z-index: 1; }
@media (max-width: 1000px) { .ps-cmp-lower { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 700px) { .ps-cmp-title { font-size: 22px; } .ps-cmp-key { margin-left: 0; } }
"""


# ---------------------------------------------------------------------------------------------------------------
# The plot (no longer on the screen: the bars of the table replaced it)


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


def themed_json(ov: dict, theme: str, display: Optional[str]) -> dict:
    """The figure of the values as the browser takes it."""
    return figure_values(ov, theme, display).to_plotly_json()


# ---------------------------------------------------------------------------------------------------------------
# The comparison of one observable


def render(ctx) -> None:
    """The comparison of the observable ``ctx.state["start"]["observable"]`` (of the level or transition "key", of
    the nucleus "nucleus"), into the current container. ``ctx`` is the start page's (:func:`.start.render`). What it
    shows is read on a worker thread (:func:`.nucleus_parts.load`), with the evaluation taken first if it never was;
    a spinner stands in meanwhile. With ``to_point_form`` in the state the window of a data point opens at once."""
    ui = ctx.ui
    st = ctx.state["start"]
    name = st["nucleus"]
    add_style(ctx)
    if not ctx.state.get("compare_style"):
        ui.add_css(COMPARE_STYLE)
        ctx.state["compare_style"] = True
    held: dict = {}
    with ui.element("div").classes("ps-cmp"):
        with ui.element("div").classes("ps-cmp-where"):
            ui.button(f"← {nuclide_text(name)}", color=None, on_click=lambda: open_nucleus(ctx, name)).props(
                "flat no-caps dense").classes("ps-cmp-back").tooltip("The page of the nucleus")
            ui.label("/")
            ui.label("Compared with the literature")
        body = ui.element("div").classes("ps-cmp").style("max-width: none")
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
            window = _page(ctx, held, lambda: later(ui, reload), lambda *a: later(ui, lambda: choose(*a)), fill)
        opening = st.pop("to_point_form", False)
        if st.pop("point_window", False) or opening:
            if window is not None:
                window.open()

    later(ui, reload)


def _page(ctx, held: dict, reload, choose, redraw):
    """The view from what was read. ``reload()`` reads everything again, ``choose(observable, key, reference)``
    another observable or reference, ``redraw()`` builds the view again from what is held. Returns the window of a
    data point (None when the page could not be read)."""
    ui = ctx.ui
    if held.get("error"):
        ui.label(held["error"]).classes("text-sm ps-warnbox rounded p-3 w-full")
        return None
    name, view, ov = held["nuclide"], held["view"], held["observable"]
    window = _point_window(ctx, name, ov, reload, redraw)
    _head(ctx, name, ov, window, reload)
    for note in held.get("notes") or ():
        ui.label(note).classes("text-sm ps-note p-2 w-full")
    _tabs(ctx, view, ov, choose)
    if ov is None:
        ui.label(f"No value of {nuclide_text(name)} is here yet. Add a data point, or plan an experiment for it."
                 ).classes("text-sm ps-note p-3 w-full")
        return window
    unit = _unit(ctx, ov)
    _tiles(ui, ov, unit)
    _table(ctx, held, unit, reload, choose, redraw)
    with ui.element("div").classes("ps-cmp-boxes"):
        _stands_out(ctx, ov)
        _projected(ctx, name, view, ov)
    with ui.element("div").classes("ps-cmp-lower"):
        with ui.element("div"):
            _routes(ctx, name, ov)
            _history(ctx, ov)
        with ui.element("div"):
            _gives(ctx, ov)
    return window


def _point_window(ctx, name: str, ov: Optional[dict], reload, redraw):
    """The window that holds the form of a data point (:func:`.nucleus_parts.point_form`, unchanged). Closed with a
    correction in it, the correction is dropped, as its Cancel does."""
    ui = ctx.ui
    st = ctx.state["start"]

    def hidden() -> None:
        if (st.get("point") or {}).get("corrects"):
            st.update(point=None, point_problems=[])
            redraw()

    with ui.dialog().on("hide", hidden) as dialog, ui.card().classes("ps-card ps-cmp-window").style(
            "width: min(760px, 96vw); max-width: 96vw"):
        ui.button(icon="close", color=None, on_click=dialog.close).props("flat round dense").classes(
            "ps-muted ps-cmp-close").tooltip("Close")
        point_form(ctx, name, ov, reload, redraw)
    return dialog


def _head(ctx, name: str, ov: Optional[dict], window, reload) -> None:
    """The title and the three actions: a data point (the window), an import, the report of the nucleus."""
    ui = ctx.ui
    with ui.element("div").classes("ps-cmp-head"):
        with ui.element("h1").classes("ps-cmp-title"):
            if ov is None:
                ui.label(f"Nothing of {nuclide_text(name)} to compare yet")
            else:
                ui.label(label_text(ov["label"])).classes("ps-cmp-sym")
                ui.label(f" of {nuclide_text(name)}, every value side by side")
        with ui.element("div").classes("ps-cmp-actions"):
            ui.button("Add a data point", icon="add", on_click=window.open).props("outline no-caps").tooltip(
                "A value with its uncertainty, its method and where it is from")
            import_button(ctx, name, reload)
            report_button(ctx, name)


def _tabs(ctx, view: dict, ov: Optional[dict], choose) -> None:
    """The observables of the nucleus as tabs: a tab per headline row (and for the one shown), the shown one marked;
    the rest, which hold the evaluation only, in a list to search."""
    ui = ctx.ui
    if not view["table"]:
        return
    on = (ov["observable"], ov["key"]) if ov is not None else None
    listed = [t for t in view["table"] if t["headline"] or (t["observable"], t["key"]) == on]
    rest = [t for t in view["table"] if not (t["headline"] or (t["observable"], t["key"]) == on)]
    with ui.element("div").classes("ps-cmp-tabs"):
        for t in listed:
            this = (t["observable"], t["key"]) == on
            ui.button(label_text(t["label"]), color=None, on_click=lambda t=t: choose(t["observable"], t["key"])
                      ).props("flat no-caps" + (" aria-current=page" if this else "")).classes(
                "ps-cmp-tab" + (" ps-cmp-tab-on" if this else ""))
        if rest:
            ui.select({f"{t['observable']}|{t['key']}": label_text(t["label"]) for t in rest}, with_input=True,
                      label=f"{len(rest)} more hold the evaluation only",
                      on_change=lambda e: e.value and choose(*e.value.split("|", 1))).props(
                "dense outlined options-dense").classes("ps-cmp-more")


def _unit(ctx, ov: dict) -> Optional[dict]:
    """The unit the values are shown in (an entry of the view's ``units``), as chosen; the canonical one when none
    is, or the chosen one is not this observable's. A B(E2) keeps the choice of the nucleus's page."""
    st = ctx.state["start"]
    units = ov.get("units") or []
    chosen = st.get("unit") if ov["observable"] == "B(E2)" else (st.get("units") or {}).get(ov["observable"])
    return next((u for u in units if u["unit"] == chosen), units[0] if units else None)


def _set_unit(ctx, ov: dict, unit: str) -> None:
    st = ctx.state["start"]
    if ov["observable"] == "B(E2)":
        st["unit"] = unit
    else:
        st.setdefault("units", {})[ov["observable"]] = unit


def _text(entry: dict, unit: Optional[dict]) -> str:
    """An entry's value and uncertainty in the unit shown."""
    return (entry.get("texts") or {}).get(unit["unit"], entry["text"]) if unit else entry["text"]


def _entries(ov: dict) -> dict:
    """Every entry of the view's groups, by id."""
    return {e["id"]: e for g in ov["groups"] for e in g["entries"]}


def _plural(count: int, word: str) -> str:
    return f"{count} {word}{'s' * (count != 1)}"


# ---------------------------------------------------------------------------------------------------------------
# The summary row


def _tiles(ui, ov: dict, unit: Optional[dict]) -> None:
    """Five tiles from the view's ``counts`` and its reference: the reference, the measured values with their
    verdicts, the best single measurement, theory, and what the setup reaches. A tile with nothing says so."""
    counts = ov["counts"]
    entries = _entries(ov)
    ref = next((e for e in entries.values() if e["is_reference"]), None)

    def tile(head: str, value: str, sub=(), extra: str = "", unit_label: str = "", none: bool = False) -> None:
        with ui.element("div").classes("ps-cmp-tile " + extra):
            ui.label(head).classes("ps-cmp-tile-h")
            with ui.element("div").classes("ps-cmp-tile-v" + (" ps-cmp-none" if none else "")):
                ui.label(value)
                if unit_label:
                    ui.label(" " + unit_label).classes("ps-cmp-tile-unit")
            with ui.element("div").classes("ps-cmp-tile-s"):
                for i, (text, colour) in enumerate([(s, "") if isinstance(s, str) else s for s in sub]):
                    if i:
                        ui.label(" · ")
                    ui.label(text).classes(colour)

    with ui.element("section").classes("ps-cmp-tiles").props('aria-label="Summary"'):
        if ref is not None:
            tile("The reference", _text(ref, unit), [" · ".join(x for x in (ref["source"], ref["known_text"]) if x)],
                 "ps-cmp-tile-ref", unit["label"] if unit else "")
        else:
            tile("The reference", "none", ["nothing evaluated yet: pick a row to compare against"],
                 "ps-cmp-tile-ref", none=True)
        measured = counts["measured"]
        if measured["count"]:
            tile("Measured", _plural(measured["count"], "value"),
                 [(p["text"], VERDICT_MARK.get(p["verdict"], ("", "ps-muted"))[1]) for p in measured["parts"]])
        else:
            tile("Measured", "none", ["no measured value"], none=True)
        best = counts["best_measured"]
        if best is not None:
            shown_best = _text(entries[best["id"]], unit) if best["id"] in entries else best["text"]
            tile("Best single measurement", shown_best,
                 [" · ".join(x for x in (best["source"], best["method"], best["made_from"], best["known_text"]) if x)])
        else:
            tile("Best single measurement", "none", ["no measured value with an uncertainty"], none=True)
        theory = counts["theory"]
        if theory["count"]:
            tile("Theory", _plural(theory["count"], "value"), [theory["text"]] if theory["text"] else [])
        else:
            tile("Theory", "none", ["no theory value"], none=True)
        projected = counts["projected"]
        if projected["count"]:
            tile("Your setup, projected", projected["reach_text"] or "—",
                 [x for x in (_plural(projected["count"], "simulated run"), projected["closing_text"]) if x])
        else:
            tile("Your setup, projected", "none", ["no simulated run"], none=True)


# ---------------------------------------------------------------------------------------------------------------
# The table


def _table(ctx, held: dict, unit: Optional[dict], reload, choose, redraw) -> None:
    """Every value of the observable in one table: the controls (order, unit, the reference, the key to the marks),
    then a row per value with its bar on the axis all rows share, under the group headers when ordered by method;
    under the table the sentence that explains the pull."""
    ui = ctx.ui
    st = ctx.state["start"]
    name, ov = held["nuclide"], held["observable"]
    c, cells = ov["comparison"], ov["cells"]
    entries = _entries(ov)
    order = st.get("order") if st.get("order") in ORDERS else "method"
    axis = ov.get("axis")
    rows = {e.id: e.row for e in c.entries} if c is not None else {}
    with ui.element("section").classes("ps-cmp-box"):
        with ui.element("div").classes("ps-cmp-controls"):
            _controls(ctx, ov, c, unit, order, choose, redraw)
            _key(ui, ov, axis)
        if not entries:
            ui.label("No value yet: add a data point, or analyse a run of an experiment.").classes(
                "text-sm ps-muted").style("padding: 12px 18px")
            return
        with ui.element("div").classes("ps-cmp-scroll"), ui.element("div").classes("ps-cmp-grid"):
            with ui.element("div").classes("ps-cmp-row ps-cmp-colhead"):
                for text in ("Source", "Method", f"Value ({unit['label']})" if unit else "Value"):
                    ui.label(text)
                with ui.element("div").classes("ps-cmp-ticks"):
                    for tick in (unit["ticks"] if unit else (axis or {}).get("ticks") or ()):
                        ui.label(tick["text"]).style(f"left: {100 * tick['x']:.2f}%")
                ui.label("Difference").classes("ps-cmp-right")
                ui.label("Pull").classes("ps-cmp-right")
                ui.label("Verdict")
                ui.label("")
            superseded = {i: s["title"] for s in ov["stands_out"] if s["kind"] == "superseded" for i in s["ids"]}

            def row(e: dict, group: str) -> None:
                _row(ctx, name, e, group, unit, axis, cells, rows, superseded, reload, redraw)

            if order == "method":
                for g in ov["groups"]:
                    with ui.element("div").classes("ps-cmp-group"):
                        ui.label(g["label"]).classes("ps-cmp-group-label")
                        for text in (g["note"], f"converted {g['conversion']}" if g["conversion"] else ""):
                            if text:
                                ui.label(words(text)).classes("ps-cmp-group-note")
                    for e in g["entries"]:
                        row(e, g["id"])
            else:
                group_of = {e["id"]: g["id"] for g in ov["groups"] for e in g["entries"]}
                for i in ov["orders"][order]:
                    if i in entries:
                        row(entries[i], group_of[i])
        if c is not None:
            ui.label(f"Pull = the difference over the two uncertainties in quadrature (the side of each that faces the "
                     f"other value). Below {c.agree_below:g}: agree; {c.agree_below:g} to {c.differ_above:g}: tension; "
                     f"above {c.differ_above:g}: differ. A value without an uncertainty gets its ratio and no verdict. "
                     f"◦ not independent of the reference: the two share a source, so their agreement is no test. "
                     f"A simulated run is closed against the value its simulation was given, never against the "
                     f"reference.").classes("ps-cmp-foot")


def _controls(ctx, ov: dict, c, unit: Optional[dict], order: str, choose, redraw) -> None:
    """Order, unit and the reference every value is compared with."""
    ui = ctx.ui
    st = ctx.state["start"]
    cells = ov["cells"]

    def reorder(e) -> None:
        if e.value and e.value != order:
            st["order"] = e.value
            redraw()

    def re_unit(e) -> None:
        if e.value and unit and e.value != unit["unit"]:
            _set_unit(ctx, ov, e.value)
            redraw()

    with ui.element("div").classes("ps-cmp-control").props('role=group aria-label="Order"'):
        ui.label("Order")
        ui.toggle(dict(ORDERS), value=order, on_change=reorder).props("dense no-caps unelevated").classes(
            "ps-cmp-order")
    units = ov.get("units") or []
    if len(units) > 1:
        with ui.element("div").classes("ps-cmp-control").props('role=group aria-label="Unit"'):
            ui.label("Unit")
            ui.toggle({u["unit"]: u["label"] for u in units}, value=unit["unit"], on_change=re_unit).props(
                "dense no-caps unelevated").classes("ps-cmp-unit").tooltip(
                "How the values are written. The rows are not changed; a value in W.u. is always the ↓ value.")
    if c is not None and c.entries:
        ref = c.reference
        entries = _entries(ov)
        with ui.element("div").classes("ps-cmp-control"):
            ui.label("against")
            options = {e.id: f"{cells[e.id]['source']} · " + (_text(entries[e.id], unit) if e.id in entries
                                                              else cell_text(cells[e.id])) for e in c.entries}
            ui.select(options, value=ref.id if ref is not None else None,
                      label=None if ref is not None else "pick a row",
                      on_change=lambda e: (st.update(reference=e.value), choose(ov["observable"], ov["key"], e.value))
                      ).props("dense outlined options-dense").style("min-width: 220px; max-width: 340px").tooltip(
                "The row every other value is compared with: the evaluated one unless another is picked")


def _key(ui, ov: dict, axis: Optional[dict]) -> None:
    """The key to the marks of the families that are in the table, the reference's band and the value a simulated
    run was given."""
    families = {e["family"] for g in ov["groups"] if g["id"] != "reference" for e in g["entries"]}
    has_ref = any(g["id"] == "reference" for g in ov["groups"])
    with ui.element("div").classes("ps-cmp-key").props('aria-label="Key to the marks"'):
        if has_ref:
            with ui.element("div"):
                ui.element("span").classes("ps-mk ps-fam-reference")
                ui.label("The reference")
        for family, text in KEY_ORDER:
            if family in families:
                with ui.element("div"):
                    ui.element("span").classes(f"ps-mk ps-fam-{_slug(family)}")
                    ui.label(text)
        if "projected" in families:
            with ui.element("div"):
                ui.element("span").classes("ps-cmp-truth")
                ui.label("Given to the simulation")
        if axis and axis.get("band"):
            with ui.element("div"):
                ui.element("span").classes("ps-cmp-keyband")
                ui.label("Reference ± its uncertainty")


def _slug(family: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", str(family).lower())


def _row(ctx, name: str, e: dict, group: str, unit: Optional[dict], axis: Optional[dict], cells: dict, rows: dict,
         superseded: dict, reload, redraw) -> None:
    """One value: source, method, value, its bar, difference, pull and verdict, and the menu of a data point."""
    ui = ctx.ui
    family = "reference" if e["is_reference"] else _slug(e["family"])
    projected = e["family"] == "projected"
    cell = cells.get(e["id"]) or {}
    with ui.element("div").classes("ps-cmp-row"):
        # -- source, with where in the paper and the row's note (whole on hover)
        with ui.element("div"):
            ui.label(e["source"]).classes("ps-cmp-src")
            if projected:
                if e["stat_sys_text"]:
                    ui.label(e["stat_sys_text"]).classes("ps-cmp-sub ps-num")
            else:
                if e["where"] and e["where"] not in e["source"]:
                    ui.label(e["where"]).classes("ps-cmp-sub")
                note = words(e["detail"])
                if note:
                    short = note if len(note) <= DETAIL_SHOWN else note[:DETAIL_SHOWN].rstrip() + "…"
                    ui.label(short).classes("ps-cmp-sub ps-cmp-note").tooltip(note)
            if e["independent_of_reference"] is False and group != "evaluation":
                ui.label("◦ not independent of the reference").classes("ps-cmp-sub")
            if e["id"] in superseded:
                ui.label(superseded[e["id"]]).classes("ps-cmp-sub ps-warn")
        # -- method, with what the value was made from or given as
        with ui.element("div"):
            ui.label(e["method"]).style("font-size: 13px")
            for text, colour in ((words(e["detail"]) if projected else "", ""),
                                 (words(e["made_from"]), ""),
                                 (f"given as {e['given_as']}" if e["given_as"] and not e["made_from"] else "", ""),
                                 (cell.get("problem") or "", "ps-warn")):
                if text:
                    ui.label(text).classes("ps-cmp-sub " + colour)
        # -- value, in the unit shown, and how well it is known
        with ui.element("div"):
            ui.label(_text(e, unit)).classes("ps-cmp-num")
            known = e["known_text"] + (" · the best" if e["is_best"] else "")
            if known:
                ui.label(known).classes("ps-cmp-sub")
        # -- the bar on the shared axis: the reference's band and line, the uncertainty, the mark
        with ui.element("div").classes("ps-cmp-bar").props("aria-hidden=true"):
            band = (axis or {}).get("band")
            if band:
                ui.element("div").classes("ps-cmp-band").style(
                    f"left: {100 * band['x_low']:.2f}%; width: {100 * (band['x_high'] - band['x_low']):.2f}%")
                ui.element("div").classes("ps-cmp-refline").style(f"left: {100 * band['x']:.2f}%")
            if e["x"] is not None:
                if e["x_low"] is not None and e["x_high"] is not None and e["x_high"] > e["x_low"]:
                    ui.element("div").classes(f"ps-cmp-err ps-fam-{family}").style(
                        f"left: {100 * e['x_low']:.2f}%; width: {100 * (e['x_high'] - e['x_low']):.2f}%")
                if e["x_truth"] is not None:
                    ui.element("div").classes("ps-cmp-truth").style(f"left: {100 * e['x_truth']:.2f}%")
                ui.element("div").classes(f"ps-mk ps-fam-{family}").style(f"left: {100 * e['x']:.2f}%")
        # -- against the reference (a simulated run: against what its simulation was given)
        ui.label(e["difference_text"]).classes("ps-cmp-num ps-cmp-right")
        ui.label(e["pull_text"]).classes("ps-cmp-num ps-cmp-right")
        mark, colour = VERDICT_MARK.get(e["verdict"], ("", "ps-muted"))
        ui.label(f"{mark} {e['verdict_text']}".strip()).classes("ps-cmp-verdict " + colour)
        with ui.element("div"):
            row = rows.get(e["id"])
            if e["is_data_point"] and row is not None and nucleus.is_data_point(row):
                _row_menu(ctx, name, row, cell, reload, redraw)


def _row_menu(ctx, name: str, row, cell: dict, reload, redraw) -> None:
    """Correct and retract, on a data point of the ledger. Neither changes the row: each adds one. Correct opens the
    window of a data point with the point in it."""
    ui = ctx.ui
    st = ctx.state["start"]

    def correct() -> None:
        st.update(point=point_of(row), point_problems=[], point_window=True)
        redraw()

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


# ---------------------------------------------------------------------------------------------------------------
# What stands out, and what the projected runs say


def _stands_out(ctx, ov: dict) -> None:
    """The view's findings (``stands_out``): each one's title and sentence, its note shortened (whole on hover)."""
    ui = ctx.ui
    with ui.element("div").classes("ps-cmp-card"):
        ui.label("What stands out").classes("ps-cmp-h2")
        if not ov["stands_out"]:
            ui.label("Nothing stands out: no value differs from the reference or is in tension with it, none is "
                     "superseded, and no theory value lies to one side.").classes("ps-muted")
            return
        for s in ov["stands_out"]:
            with ui.element("div"):
                ui.label(s["title"]).classes("font-medium " + FINDING_COLOUR.get(s["kind"], ""))
                ui.label(s["text"])
                note = words(s["note"])
                if note:
                    short = note if len(note) <= 2 * DETAIL_SHOWN else note[:2 * DETAIL_SHOWN].rstrip() + "…"
                    ui.label(short).classes("text-xs ps-muted").tooltip(note)


def _projected(ctx, name: str, view: dict, ov: dict) -> None:
    """What the simulated runs say about the setup (the sentence of ``counts.projected``): never a test of the value.
    A link opens the experiment of the best run; without a run, the way to plan one."""
    ui = ctx.ui
    projected = ov["counts"]["projected"]
    with ui.element("div").classes("ps-cmp-card"):
        ui.label("What the projected runs say").classes("ps-cmp-h2")
        if not projected["count"]:
            ui.label("No simulated run of this is analysed yet. A simulated run analysed back says what a setup would "
                     "reach, not what the value is.").classes("ps-muted")
            live = next((r for r in ov["routes"] if r["live"] and r["kind"] and r.get("applies") is not False), None)
            ui.button(f"Plan an experiment for {nuclide_text(name)} →", color=None,
                      on_click=lambda: plan(ctx, name, live["kind"] if live else None)).props(
                "flat no-caps").classes("ps-cmp-link")
            return
        if projected["text"]:
            ui.label(projected["text"])
        ui.label("A simulated run analysed back says what the setup would reach, not what the value is.").classes(
            "ps-muted")
        best = _entries(ov).get(projected["best"]) or {}
        folder = best.get("experiment")
        path = next((x["path"] for x in view["experiments"] if folder and Path(x["path"]).name == folder), None)
        if path:
            ui.button(f"Open {folder} →", color=None, on_click=lambda: open_experiment(ctx, path)).props(
                "flat no-caps").classes("ps-cmp-link").tooltip(f"The experiment of the best run, {best['source']}")


# ---------------------------------------------------------------------------------------------------------------
# Under them: the routes, what the value gives, the history


def _routes(ctx, name: str, ov: dict) -> None:
    """The registry's routes to the observable, each with its state for this level or transition, and "Plan it"
    where physim has the kind of experiment: it opens that kind's new-experiment form for this nucleus."""
    ui = ctx.ui
    title = nuclide_text(name)
    ui.label(f"Routes to {label_text(ov['label'])} in {title}").classes("ps-cmp-h2")
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
        listed = [f"{r['label']}" + (f" ({r['rows']} value{'s' * (r['rows'] != 1)} here)" if r["rows"] else "")
                  for r in done + [r for r in quiet if not r["rows"]]]
        ui.label("Not experiments of the planner (the evaluation, a calculation, methods known from the literature: "
                 "their values come in as data points): " + " · ".join(listed) + ".").classes("text-xs ps-muted")
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
    said = {"superseded": "superseded by a newer row", "retracted": "withdrawn", "retraction": "the retraction"}
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
                                          said.get(h["status"], "") + (f": {h['reason']}" if h["reason"] else "")],
                        colour="ps-muted" if gone else "")
                    td(ui, "" if h["status"] == "retraction" else h["text"].replace(
                        " " + h["unit"], " " + unit_text(h["unit"])), classes="ps-nt-num",
                        colour="ps-muted" if gone else "")
        if len(rows) > HISTORY_ROWS:
            ui.label(f"… and {len(rows) - HISTORY_ROWS} older rows.").classes("text-xs ps-muted p-2")
