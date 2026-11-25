"""The landing page of one nucleus: what it is, what is known and projected about it, its experiments and what is
still open.

The view of ``/?nucleus=62Ni`` (and of a card of the nuclei screen, :mod:`.start`), from the top: the nucleus with
its counts and the four things one does from here, beside the card of its lowest levels; one card per headline
observable (the first a large one, with every value on one axis); the experiments on it beside what is still open;
and at the bottom all its observables as a table. A card, a row or an arrow of the levels opens the comparison of
that observable with the literature (:mod:`.tab_compare`).

Nothing is computed here: every number and count comes from :mod:`physim.nuclear.nucleus` through
:func:`.nucleus_parts.load_landing` as text; the page only places the dots of a strip at their ``x`` and the levels
at heights that follow their energies. What the views of a nucleus share is in :mod:`.nucleus_parts`; other modules
import its names from here too.
"""

from __future__ import annotations

import math
import re
from typing import Optional

from .. import observables  # noqa: F401 -- taken from here by the tests
from ..compare import VERDICTS
from ..results import EVIDENCE, METHODS
from ..workbench import STEPS
from .nucleus_parts import (BE2_UNITS, add_style, as_given, cell_text, experiments_table, import_button,
                            label_text, later, load_landing, nuclide_text, open_observable, open_routes, panel,
                            percent, plan, point_form, projected_text, range_text, report_button,
                            scroll_to_point_form, standing, summary_text, table_head, td, unit_text, value_parts,
                            waiting)
from .nucleus_parts import (  # noqa: F401 -- other modules and the tests take these from here
    add_point, blank_point, display_factor, import_file, import_lines, load, load_nuclei, load_observable,
    open_nucleus, point_arguments, point_of, report_html, shown, url, uses)
from .tab_compare import figure_values  # noqa: F401 -- taken from here by the tests

#: The landing view's own style, in the page's colours (both themes): the head, the levels, the cards, the strip.
LANDING_STYLE = """
.ps-l { display: flex; flex-direction: column; gap: 28px; width: 100%; max-width: 1280px; margin: 0 auto;
  box-sizing: border-box; }
.ps-l-top { display: flex; flex-wrap: wrap; gap: 24px 40px; align-items: stretch; }
.ps-l-head { flex: 999 1 520px; min-width: 0; display: flex; flex-direction: column; gap: 14px; }
.ps-l-eyebrow { font-size: 12px; font-weight: 600; letter-spacing: .08em; text-transform: uppercase;
  color: var(--ps-muted); }
.ps-l-name { display: flex; flex-wrap: wrap; align-items: baseline; gap: 6px 18px; }
.ps-l-nuclide { font-family: var(--ps-mono); font-size: 76px; font-weight: 500; line-height: 1;
  letter-spacing: -.02em; white-space: nowrap; }
.ps-l-nuclide > div { display: inline; }
.ps-l-nuclide .ps-l-mass { font-size: .5em; vertical-align: .72em; margin-right: .04em; letter-spacing: 0; }
.ps-l-element { font-size: 20px; color: var(--ps-muted); }
.ps-l-zn { font-family: var(--ps-mono); font-size: 14px; }
.ps-l-said { color: var(--ps-muted); max-width: 640px; line-height: 1.45; }
.ps-l-actions { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 6px; }
.ps-l-card { background: var(--ps-surface); border: 1px solid var(--ps-line); border-radius: 8px;
  box-sizing: border-box; min-width: 0; }
.ps-l-levels { flex: 1 1 360px; padding: 16px 20px; display: flex; flex-direction: column; gap: 10px; }
.ps-l-spread { display: flex; flex-wrap: wrap; justify-content: space-between; align-items: baseline;
  gap: 4px 12px; }
.ps-l-scheme { position: relative; font-family: var(--ps-mono); font-size: 13px; }
.ps-l-line { position: absolute; height: 2px; background: var(--ps-ink); }
.ps-l-lead { position: absolute; height: 1px; background: var(--ps-ink); }
.ps-l-jpi, .ps-l-e { position: absolute; line-height: 16px; height: 16px; white-space: nowrap; }
.ps-l-e { text-align: right; }
.ps-l-arrow { position: absolute; width: 2px; background: var(--ps-accent); }
.ps-l-tip { position: absolute; width: 0; height: 0; border-left: 5px solid transparent;
  border-right: 5px solid transparent; border-top: 9px solid var(--ps-accent); }
.ps-l-arrow-text { position: absolute; font-size: 12px; line-height: 16px; color: var(--ps-accent);
  background: var(--ps-surface); padding: 0 4px; cursor: pointer; }
.ps-l-foot { display: flex; flex-wrap: wrap; justify-content: space-between; gap: 6px 12px; font-size: 13px;
  border-top: 1px solid var(--ps-hair); padding-top: 10px; }
.ps-l a.ps-l-link, .ps-l .ps-l-link { color: var(--ps-accent); text-decoration: none; cursor: pointer; }
.ps-l a.ps-l-link:hover { text-decoration: underline; }
.ps-l-section { display: flex; flex-direction: column; gap: 12px; }
.ps-l-title { display: flex; flex-wrap: wrap; align-items: baseline; gap: 4px 12px; }
.ps-l-h2 { margin: 0; font-size: 18px; font-weight: 600; }
.ps-l-cards { display: flex; flex-wrap: wrap; gap: 12px; align-items: stretch; }
.ps-l-main { flex: 2 1 520px; padding: 18px 20px; display: flex; flex-direction: column; gap: 14px; }
.ps-l-small { flex: 1 1 220px; padding: 18px 20px; display: flex; flex-direction: column; gap: 8px; }
.ps-l-label { font-family: var(--ps-mono); font-weight: 600; }
.ps-l-value { display: flex; flex-wrap: wrap; align-items: baseline; gap: 4px 12px; }
.ps-l-big { font-family: var(--ps-mono); font-size: 40px; font-weight: 500; line-height: 1.1;
  color: var(--ps-accent); }
.ps-l-unit { font-family: var(--ps-mono); font-size: 16px; color: var(--ps-muted); }
.ps-l-mid { font-family: var(--ps-mono); font-size: 24px; font-weight: 500; line-height: 1.2; }
.ps-l-mid-unit { font-family: var(--ps-mono); font-size: 14px; color: var(--ps-muted); }
.ps-l-small-text { font-size: 13px; }
.ps-l-strip { position: relative; height: 34px; border-bottom: 1px solid var(--ps-line); }
.ps-l-band { position: absolute; top: 0; bottom: 0;
  background: color-mix(in srgb, var(--ps-accent) 18%, transparent); }
.ps-l-ref { position: absolute; top: 0; bottom: 0; width: 1px; margin-left: -0.5px; background: var(--ps-accent); }
.ps-l-ref-hit { position: absolute; top: 0; bottom: 0; width: 9px; margin-left: -4.5px; }
.ps-l-dot { position: absolute; top: 12px; width: 9px; height: 9px; margin-left: -4.5px; border-radius: 50%;
  box-sizing: border-box; border: 2px solid var(--ps-muted); background: var(--ps-muted); }
.ps-l-dot-tension { border-color: var(--ps-warn); background: var(--ps-warn); }
.ps-l-dot-differ { border-color: var(--ps-bad); background: var(--ps-bad); }
.ps-l-dot-theory { background: transparent; }
.ps-l-dot-outside { background: transparent; border-radius: 1px; }
.ps-l-dot-projected { border-color: var(--ps-accent); background: transparent; }
.ps-l-dot-evaluation { border-radius: 1px; transform: rotate(45deg); width: 8px; height: 8px; margin-left: -4px;
  border-color: var(--ps-ink); background: var(--ps-surface); }
.ps-l-ticks { position: relative; height: 16px; font-family: var(--ps-mono); font-size: 11px;
  color: var(--ps-muted); }
.ps-l-tick { position: absolute; transform: translateX(-50%); white-space: nowrap; }
.ps-l-axis-unit { position: absolute; right: 0; top: 16px; white-space: nowrap; }
.ps-l-counts { display: flex; flex-wrap: wrap; gap: 6px 20px; font-size: 13px; }
.ps-l-counts b { font-weight: 600; }
.ps-l-end { display: flex; flex-wrap: wrap; justify-content: space-between; align-items: center; gap: 8px 12px;
  border-top: 1px solid var(--ps-hair); padding-top: 12px; }
.ps-l .q-btn.ps-l-go { color: var(--ps-accent) !important; font-weight: 600; }
.ps-l .q-btn.ps-l-go.q-btn--outline::before { border-color: var(--ps-accent) !important; }
.ps-l-row { display: flex; flex-wrap: wrap; gap: 20px; align-items: flex-start; }
.ps-l-exps { flex: 999 1 620px; min-width: 0; }
.ps-l-open { flex: 1 1 340px; min-width: 0; }
.ps-l-item { display: flex; flex-direction: column; align-items: flex-start; gap: 4px; padding: 12px 20px;
  border-bottom: 1px solid var(--ps-hair); }
.ps-l-item:last-child { border-bottom: 0; }
.ps-nt td.ps-nt-open { color: var(--ps-accent); font-weight: 600; white-space: nowrap; text-align: right;
  vertical-align: middle; }
.ps-l-all { scroll-margin-top: 70px; }
@media (max-width: 700px) { .ps-l-nuclide { font-size: 56px; } }
"""

#: The drawing of the levels, in px: its height, the room below the ground state and above the highest level,
#: the least distance between two labels, the widths of the columns of Jπ and of the energies, and of the short
#: slanted pieces that join a label moved apart to its line.
SCHEME = {"height": 200, "bottom": 12, "top": 14, "gap": 16, "jpi": 52, "energy": 92, "lead": 14}

#: How a dot of the strip is drawn, by what it is (the reference is the band and its line, not a dot).
DOT_WORDS = {"measured": "measured", "tension": "measured, in tension", "differ": "measured, differs",
             "evaluation": "another evaluation", "outside": "from an outside program", "theory": "theory",
             "projected": "projected (a simulated run)"}


# ---------------------------------------------------------------------------------------------------------------
# The page of one nucleus


def render(ctx) -> None:
    """The landing page of the nucleus ``ctx.state["start"]["nucleus"]``, into the current container. ``ctx`` is
    the start page's (:func:`.start.render`). What it shows is read on a worker thread
    (:func:`.nucleus_parts.load_landing`), with the evaluation taken first if it never was; a spinner stands in
    meanwhile. The routes still open are tested after the page shows (:func:`.nucleus_parts.open_routes`)."""
    ui = ctx.ui
    st = ctx.state["start"]
    name = st["nucleus"]
    add_style(ctx)
    if not ctx.state.get("landing_style"):
        ui.add_css(LANDING_STYLE)
        ctx.state["landing_style"] = True
    held: dict = {}
    body = ui.column().classes("w-full gap-0")
    with body:
        waiting(ui, f"Reading what is known about {nuclide_text(name)}…")

    async def reload() -> None:
        """Read everything again and build the page from it, in place; then test the routes still open."""
        from nicegui import run

        try:
            data = await run.io_bound(lambda: load_landing(name, first=not held))
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty page
            data = {"error": f"The page of {name} could not be read: {err}"}
        if data is None or body.is_deleted:
            return
        held.clear()
        held.update(data)
        fill()
        if held.get("error"):
            return
        try:
            found = await run.io_bound(open_routes, held["nuclide"], held["view"], held["rows"])
        except Exception as err:  # noqa: BLE001 -- the rest of the page stands
            found = {"error": f"The routes could not be tested: {err}"}
        if body.is_deleted or held.get("view") is not data.get("view"):
            return
        held["routes"] = found
        if held.get("show_open"):
            held["show_open"]()

    def fill() -> None:
        body.clear()
        with body:
            _page(ctx, held, lambda: later(ui, reload), fill)

    later(ui, reload)


def _page(ctx, held: dict, reload, redraw) -> None:
    """The page from what was read. ``reload()`` reads everything again, ``redraw()`` builds the page again from what
    is held (the unit switch)."""
    ui = ctx.ui
    if held.get("error"):
        ui.label(held["error"]).classes("text-sm ps-warnbox rounded p-3 w-full")
        return
    name, view = held["nuclide"], held["view"]
    with ui.element("div").classes("ps-l"):
        with ui.element("section").classes("ps-l-top"):
            _head(ctx, held, reload)
            _levels(ctx, held)
        if view["table"]:
            _cards(ctx, view, redraw)
            with ui.element("section").classes("ps-l-row"):
                with ui.element("div").classes("ps-l-exps"):
                    _experiments(ctx, view)
                with ui.element("div").classes("ps-l-open"):
                    _still_open(ctx, held)
            _observables(ctx, view, redraw, held)
            return
        ui.label(f"No value of {nuclide_text(name)} is here yet. Add a data point below, or plan an experiment "
                 "for it.").classes("text-sm ps-note p-3 w-full")
        with ui.element("div").classes("ps-two"):
            with ui.column().classes("w-full gap-3 min-w-0"):
                _experiments(ctx, view)
                _still_open(ctx, held)
            with ui.column().classes("w-full gap-3 min-w-0"):
                point_form(ctx, name, None, reload, redraw)


def _plural(count: int, word: str) -> str:
    return f"{count} {word}{'s' * (count != 1)}"


def _head(ctx, held: dict, reload) -> None:
    """The nucleus large, what it is, what the page holds about it, and the four things one does from here."""
    ui = ctx.ui
    view = held["view"]
    name = view["nuclide"]
    title = nuclide_text(name)
    studied = [x for x in view["experiments"] if x["role"] == "studied"]
    evaluation = view["evaluation"]
    with ui.element("div").classes("ps-l-head"):
        ui.label("The nucleus you are studying").classes("ps-l-eyebrow")
        with ui.element("div").classes("ps-l-name"):
            # The mass number raised by the style, not as superscript characters: a monospaced font gives those the
            # width of a whole digit each.
            mass, symbol = re.fullmatch(r"(\d*)(.*)", name).groups()
            with ui.element("div").classes("ps-l-nuclide").props(
                    f'role="heading" aria-level="1" aria-label="{title}"'):
                ui.label(mass).classes("ps-l-mass")
                ui.label(symbol).classes("ps-l-symbol")
            if held.get("element"):
                ui.label(held["element"]).classes("ps-l-element")
        ui.label(" · ".join(x for x in (f"Z = {view['z']}", f"N = {view['n']}", held.get("stability")) if x)
                 ).classes("ps-l-zn")
        said = (f"Everything known and projected about it, by method. {_plural(len(studied), 'experiment')}, "
                f"{_plural(view['analysed_runs'], 'analysed run')} and {_plural(view['points'], 'data point')}.")
        if evaluation["taken"]:
            said += (f" The evaluation is {evaluation['reference']}: {_plural(evaluation['rows'], 'value')}, taken "
                     f"{evaluation['taken'][:10]}.")
        ui.label(said).classes("ps-l-said")
        for note in held.get("notes") or ():
            ui.label(note).classes("text-sm ps-note p-2 w-full")
        with ui.element("div").classes("ps-l-actions"):
            ui.button(f"New experiment for {title}", on_click=lambda: plan(ctx, name)).props(
                "unelevated no-caps").tooltip("The kinds of experiment, with this nucleus filled in")
            ui.button("Add a data point", on_click=lambda: _to_point_form(ctx, view)).props(
                "outline no-caps").tooltip("The form is with the values of an observable")
            import_button(ctx, name, reload)
            report_button(ctx, name)


def _to_point_form(ctx, view: dict) -> None:
    """The form of a data point: on this page for a nucleus without a row, otherwise with the values of the first
    observable of its table (the form takes any other too)."""
    if not view["table"]:
        scroll_to_point_form(ctx.ui)
        return
    first = view["table"][0]
    open_observable(ctx, view["nuclide"], first["observable"], first["key"], to_point_form=True)


# ---------------------------------------------------------------------------------------------------------------
# The lowest levels


def label_places(heights: list, gap: float, low: float, high: float) -> list:
    """Where the labels of lines at ``heights`` (px from the bottom, the lowest first) go so that no two are closer
    than ``gap``: each at its line where there is room; where lines are close, the labels are pushed up from below
    one after the other, and if the top one would then pass ``high`` they are pushed down from there. A label that
    moved is joined to its line by a short slanted piece (:func:`_levels`)."""
    out: list = []
    for h in heights:
        out.append(max(h, low, out[-1] + gap) if out else max(h, low))
    if out and out[-1] > high:
        out[-1] = high
        for i in range(len(out) - 2, -1, -1):
            out[i] = min(out[i], out[i + 1] - gap)
    return out


def _levels(ctx, held: dict) -> None:
    """The card of the lowest levels (:func:`physim.nuclear.nucleus.lowest_levels`): each level a line at a height
    that follows its energy, Jπ at the left and the energy at the right; a transition that has a headline
    observable is an arrow in the accent colour with its value, and opens its comparison. Under it E(4⁺₁)/E(2⁺₁)
    when the data give it, and the way to all the evaluated values."""
    ui = ctx.ui
    found = held.get("levels") or {}
    view = held["view"]
    display = ctx.state["start"].get("unit")
    levels = found.get("levels") or []
    with ui.element("div").classes("ps-l-card ps-l-levels"):
        with ui.element("div").classes("ps-l-spread"):
            ui.label("Lowest levels").classes("font-semibold")
            if found.get("reference"):
                ui.label(found["reference"]).classes("text-xs ps-muted")
        if not levels:
            ui.label(found.get("reason") or "No level is here.").classes("text-sm ps-muted ps-l-no-levels")
        else:
            _scheme(ctx, view["nuclide"], levels, display)
        ratio = found.get("ratio")
        if ratio or found.get("evaluated"):
            with ui.element("div").classes("ps-l-foot"):
                if ratio:
                    ui.label(f"{ratio['label']} = {ratio['text'].split(' = ', 1)[0]}").classes(
                        "ps-num ps-l-ratio").tooltip(f"{ratio['text']}. {ratio['how']}")
                else:
                    ui.element("span")
                if found.get("evaluated"):
                    link = ui.link(f"All {found['evaluated']} evaluated values", "#ps-all-observables").classes(
                        "ps-l-link")
                    link.on("click", lambda: _open_all(held))


def _open_all(held: dict) -> None:
    """The link of the levels card: the section of all observables, with its list of the evaluation's values
    open."""
    more = held.get("more")
    if more is not None and not more.is_deleted:
        more.open()


def _scheme(ctx, nuclide: str, levels: list, display: Optional[str]) -> None:
    ui = ctx.ui
    g = SCHEME
    top = levels[-1]["energy"] or 1.0
    span = g["height"] - g["bottom"] - g["top"]
    heights = [g["bottom"] + span * lv["energy"] / top for lv in levels]
    labels = label_places(heights, g["gap"], g["gap"] / 2, g["height"] - g["gap"] / 2)
    keys = [lv["key"] for lv in levels]
    arrows = []
    for i, lv in enumerate(levels):
        for o in lv["observables"]:
            if o["role"] == "from" and o["headline"] and o["other"] and o["other"]["key"] in keys \
                    and all((a[0]["observable"], a[0]["key"]) != (o["observable"], o["key"]) for a in arrows):
                arrows.append((o, keys.index(o["other"]["key"]), i))
    said = ", ".join(f"{lv['jpi']} at {lv['short']}" for lv in levels)
    said += "".join(f"; {label_text(o['label'])}" + (f" = {cell_text(o['value'], display)}" if o["value"] else "")
                    for o, _, _ in arrows)
    scheme = ui.element("div").classes("ps-l-scheme").style(f"height: {g['height']}px")
    scheme.props(f'role="img" aria-label="Level scheme: {said.replace(chr(34), "")}"')
    with scheme:
        for lv, y, at in zip(levels, heights, labels):
            ui.element("div").classes("ps-l-line").style(
                f"left: {g['jpi'] + g['lead']}px; right: {g['energy'] + g['lead']}px; bottom: {y - 1:.1f}px")
            rise = y - at
            length, angle = math.hypot(g["lead"], rise), math.degrees(math.atan2(rise, g["lead"]))
            ui.element("div").classes("ps-l-lead").style(
                f"left: {g['jpi']}px; bottom: {at - 0.5:.1f}px; width: {length:.1f}px; transform-origin: 0 50%; "
                f"transform: rotate({-angle:.2f}deg)")
            ui.element("div").classes("ps-l-lead").style(
                f"right: {g['energy']}px; bottom: {at - 0.5:.1f}px; width: {length:.1f}px; "
                f"transform-origin: 100% 50%; transform: rotate({angle:.2f}deg)")
            ui.label(lv["jpi"]).classes("ps-l-jpi").style(f"left: 0; width: {g['jpi'] - 4}px; bottom: {at - 8:.1f}px")
            ui.label(lv["short"]).classes("ps-l-e").style(
                f"right: 0; width: {g['energy'] - 4}px; bottom: {at - 8:.1f}px").tooltip(lv["text"])
        for n, (o, low, high) in enumerate(arrows):
            x = g["jpi"] + g["lead"] + 26 + 34 * n
            ui.element("div").classes("ps-l-arrow").style(
                f"left: {x}px; bottom: {heights[low] + 9:.1f}px; height: {heights[high] - heights[low] - 10:.1f}px")
            ui.element("div").classes("ps-l-tip").style(f"left: {x - 4}px; bottom: {heights[low] + 1:.1f}px")
            multipolarity = o["key"].rsplit(":", 1)[-1]
            text = f"{multipolarity} · {observables.get(o['observable']).label}"
            if o["value"]:
                text += f" {cell_text(o['value'], display)}"
            mid = (heights[low] + heights[high]) / 2
            tag = ui.label(text).classes("ps-l-arrow-text").style(
                f"left: {x + 8}px; bottom: {mid - 8:.1f}px; max-width: calc(100% - {x + 8}px)")
            tag.tooltip(f"{label_text(o['label'])}: compare with the literature")
            tag.on("click", lambda o=o: open_observable(ctx, nuclide, o["observable"], o["key"]))


# ---------------------------------------------------------------------------------------------------------------
# What is known and projected: one card per headline observable


def _cards(ctx, view: dict, redraw) -> None:
    """One card per headline observable: the first large, with every value on one axis; the others small."""
    ui = ctx.ui
    st = ctx.state["start"]
    display = st.get("unit")
    headline = [t for t in view["table"] if t["headline"]]
    with ui.element("section").classes("ps-l-section"):
        with ui.element("div").classes("ps-l-title w-full"):
            ui.label("What is known and projected").classes("ps-l-h2").props('role="heading" aria-level="2"')
            ui.label("one card per observable; open one to compare every value with the literature").classes(
                "ps-muted")
            ui.element("span").style("flex: 1 1 auto")
            if any(t["observable"] == "B(E2)" for t in view["table"]):
                def switch(e) -> None:
                    st["unit"] = e.value
                    redraw()

                ui.toggle({u: f"B(E2)↑ in {unit_text(u)}" for u in BE2_UNITS}, value=display or BE2_UNITS[0],
                          on_change=switch).props("dense no-caps unelevated").classes("ps-l-units").tooltip(
                    "How a B(E2) is written here. The rows are not changed; a value in W.u. is always the ↓ value.")
        if not headline:
            ui.label("Nothing beyond the evaluation yet: its values are under All observables, at the bottom."
                     ).classes("text-sm ps-muted")
            return
        with ui.element("div").classes("ps-l-cards"):
            _main_card(ctx, view["nuclide"], headline[0], display)
            for t in headline[1:]:
                _small_card(ctx, view["nuclide"], t, display)


def _evidence_line(t: dict, best: Optional[dict]) -> str:
    """What the value of a card is: "evaluated · given as 12.13 W.u.↓", or for a row without an evaluated value
    the evidence and source of its best independent value."""
    if t["evaluated"]:
        given = as_given(t["evaluated"])
        return "evaluated" + (f" · given as {given}" if given else "")
    if best:
        return f"{EVIDENCE[best['evidence']].lower()} · {best['source']} · no evaluated value"
    return "no value yet"


def _open_line(t: dict) -> str:
    return ("Still open: " + ", ".join(METHODS[m] for m in t["open"])) if t["open"] else "Nothing physim can plan"


def _projected_line(t: dict) -> str:
    counts = (t.get("counts") or {}).get("projected") or {}
    if counts.get("count"):
        return "Projected: " + " · ".join(x for x in (counts["reach_text"], counts["closing_text"]) if x)
    if t["projected"]:
        return "Projected: " + projected_text(t["projected"])[0]
    return "No simulated run"


def _main_card(ctx, nuclide: str, t: dict, display: Optional[str]) -> None:
    """The large card: the value large in the accent colour, every value on one axis, the counts in the data's
    own words, what is still open and "Compare with the literature →"."""
    ui = ctx.ui
    best = standing(t)
    with ui.element("article").classes("ps-l-card ps-l-main"):
        with ui.element("div").classes("ps-l-spread"):
            ui.label(label_text(t["label"])).classes("ps-l-label")
            ui.label(_evidence_line(t, best)).classes("text-xs ps-muted")
        with ui.element("div").classes("ps-l-value"):
            if best:
                number, unit = value_parts(best, display)
                ui.label(number).classes("ps-l-big")
                known = percent(best["uncertainty_percent"]) if best["uncertainty_percent"] is not None else ""
                ui.label(unit + (f" · known to {known}" if known else "")).classes("ps-l-unit")
            else:
                ui.label("no independent value yet").classes("ps-l-unit")
        if t.get("strip"):
            _strip(ui, t["strip"])
        counts = t.get("counts")
        with ui.element("div").classes("ps-l-counts"):
            if counts:
                _counts(ui, counts)
            else:
                for kind, word in (("measured", "measured"), ("theory", "theory")):
                    ui.label(range_text(t[kind], display, word))
        with ui.element("div").classes("ps-l-end"):
            ui.label(_open_line(t)).classes("text-sm ps-muted")
            ui.button("Compare with the literature →",
                      on_click=lambda: open_observable(ctx, nuclide, t["observable"], t["key"])).props(
                "outline no-caps").classes("ps-l-go")


def _counts(ui, counts: dict) -> None:
    """The counts of a card in the data's own words, a verdict that differs or is in tension in its colour."""
    measured = counts["measured"]
    with ui.element("span"):
        if measured["count"]:
            ui.label(f"{measured['count']} measured").classes("inline font-semibold")
            colours = {"tension": "ps-warn", "differ": "ps-bad", "not_judged": "ps-muted"}
            for i, part in enumerate(measured["parts"]):
                ui.label((" · " if i == 0 else ", ") + part["text"]).classes(
                    "inline " + colours.get(part["verdict"], ""))
        else:
            ui.label("No measured value").classes("inline")
    if counts["evaluations"]:
        with ui.element("span"):
            ui.label(f"{_plural(counts['evaluations'], 'other evaluation')}").classes("inline font-semibold")
    if counts["outside_program"]:
        with ui.element("span"):
            ui.label(f"{counts['outside_program']} from an outside program").classes("inline font-semibold")
    theory = counts["theory"]
    if theory["count"]:
        with ui.element("span"):
            ui.label(f"{theory['count']} theory").classes("inline font-semibold")
            if theory["text"]:
                ui.label(" · " + theory["text"]).classes("inline")
    projected = counts["projected"]
    with ui.element("span"):
        if projected["count"]:
            ui.label("Projected").classes("inline font-semibold")
            ui.label(" · " + " · ".join(x for x in (projected["reach_text"], projected["closing_text"]) if x)
                     ).classes("inline")
        else:
            ui.label("No simulated run").classes("inline ps-muted")


def _dot_kind(dot: dict) -> str:
    if dot["evidence"] == "projected":
        return "projected"
    if dot["evidence"] == "theory":
        return "theory"
    if dot["evidence"] == "evaluated":
        return "evaluation"
    if dot["evidence"] == "outside program":
        return "outside"
    return dot["verdict"] if dot["verdict"] in ("tension", "differ") else "measured"


def _strip(ui, strip: dict) -> None:
    """Every value on one axis: the reference as its band and line, the others as dots (one that differs in the bad
    colour, one in tension in the warning colour, theory hollow, a simulated run hollow in the accent colour,
    another evaluation a diamond), each with its source and value on hovering; the ticks under it. The axis and the
    values are those of the data, in the canonical unit, which is written under the axis."""
    axis = strip["axis"]
    unit = axis["label"]
    dots = [d for d in strip["dots"] if not d["is_reference"]]
    reference = next((d for d in strip["dots"] if d["is_reference"]), None)
    said = (f"Every value on one axis, {axis['ticks'][0]['text'] if axis['ticks'] else ''} to "
            f"{axis['ticks'][-1]['text'] if axis['ticks'] else ''} {unit}: {len(dots)} besides the reference")
    with ui.element("div").classes("w-full").style("display: flex; flex-direction: column; gap: 4px"):
        line = ui.element("div").classes("ps-l-strip").props(f'role="img" aria-label="{said}"')
        with line:
            band = axis.get("band")
            if band:
                ui.element("div").classes("ps-l-band").style(
                    f"left: {100 * band['x_low']:.3f}%; width: {100 * (band['x_high'] - band['x_low']):.3f}%")
                ui.element("div").classes("ps-l-ref").style(f"left: {100 * band['x']:.3f}%")
                if reference:
                    ui.element("div").classes("ps-l-ref-hit").style(f"left: {100 * band['x']:.3f}%").tooltip(
                        f"The reference: {reference['source']}, {reference['text']} {unit}")
            order = ("evaluation", "theory", "outside", "measured", "tension", "differ", "projected")
            for d in sorted(dots, key=lambda d: order.index(_dot_kind(d))):
                kind = _dot_kind(d)
                ui.element("div").classes(f"ps-l-dot ps-l-dot-{kind}").style(f"left: {100 * d['x']:.3f}%").tooltip(
                    f"{d['source']}: {d['text']} {unit} · {DOT_WORDS[kind]}"
                    + (f" · {VERDICTS[d['verdict']].lower()}" if kind == "measured" and d["verdict"] in VERDICTS
                       else ""))
        with ui.element("div").classes("ps-l-ticks").style("height: 30px"):
            for tick in axis["ticks"]:
                ui.label(tick["text"]).classes("ps-l-tick").style(f"left: {100 * tick['x']:.3f}%")
            ui.label(unit).classes("ps-l-axis-unit")


def _small_card(ctx, nuclide: str, t: dict, display: Optional[str]) -> None:
    """A small card: label, value, what the value is, the measured count and range, the projected line, what is
    still open and "Compare →"."""
    ui = ctx.ui
    best = standing(t)
    with ui.element("article").classes("ps-l-card ps-l-small"):
        ui.label(label_text(t["label"])).classes("ps-l-label")
        if best:
            number, unit = value_parts(best, display)
            with ui.element("div").classes("ps-l-value").style("gap: 4px 8px"):
                ui.label(number).classes("ps-l-mid")
                ui.label(unit).classes("ps-l-mid-unit")
        else:
            ui.label("no value yet").classes("ps-l-mid ps-muted")
        ui.label(_evidence_line(t, best)).classes("text-xs ps-muted")
        ui.label(range_text(t["measured"], display)).classes("ps-l-small-text")
        ui.label(_projected_line(t)).classes("ps-l-small-text")
        ui.label(_open_line(t)).classes("ps-l-small-text ps-muted")
        ui.element("span").style("flex: 1 1 auto")
        ui.button("Compare →", on_click=lambda: open_observable(ctx, nuclide, t["observable"], t["key"])).props(
            "flat no-caps").classes("ps-l-go self-start").style("padding-left: 0")


# ---------------------------------------------------------------------------------------------------------------
# The experiments, what is still open, and all observables


def _experiments(ctx, view: dict) -> None:
    ui = ctx.ui
    found = view["experiments"]
    steps = ", ".join(label for _, label in STEPS)
    _, body = panel(ui, f"Your experiments on {nuclide_text(view['nuclide'])}",
                    f"opening one shows its steps: {steps}")
    with body:
        if not found:
            ui.label("No experiment with this nucleus yet.").classes("text-sm ps-muted p-3")
            return
        experiments_table(ctx, found, ctx.state["start"].get("unit"), opener=True)


def _still_open(ctx, held: dict) -> None:
    """What is still open: the routes physim can plan that have no run of one's own yet, each with physim's test of
    it in the data's words and "Plan it", and the runs not analysed. The routes are tested after the page shows."""
    ui = ctx.ui
    view = held["view"]
    name = view["nuclide"]
    title = nuclide_text(name)
    studied = [x for x in view["experiments"] if x["role"] == "studied"]
    _, body = panel(ui, "Still open")

    @ui.refreshable
    def items() -> None:
        found = held.get("routes")
        others = [o for o in view["open"] if o["what"] != "route"]
        if found is None and any(o["what"] == "route" for o in view["open"]):
            waiting(ui, "Testing the routes physim can plan…")
        elif isinstance(found, dict):
            ui.label(found["error"]).classes("text-sm ps-warn p-3")
        for r in found if isinstance(found, list) else ():
            with ui.element("div").classes("ps-l-item ps-l-route"):
                ui.label(r["label"]).classes("font-medium")
                status = ("In physim, no run of yours yet." if not r["independent"] else
                          f"{_plural(r['independent'], 'published value')} by this method for "
                          f"{label_text(r['for'])}, none of yours.")
                why = f" Physim's test: {r['reason']}." if r["reason"] else ""
                gives = ", ".join(label_text(g["label"]) for g in r["gives"])
                ui.label(f"{status}{why} Gives {gives}.").classes("text-sm " + (
                    "ps-warn" if r["applies"] is False else "ps-muted"))
                if r["kind"] and r["applies"] is not False:
                    ui.button(f"Plan it for {title} →", on_click=lambda kind=r["kind"]: plan(ctx, name, kind)).props(
                        "flat no-caps").classes("ps-l-go").style("padding-left: 0").tooltip(
                        f"Opens the new-experiment form for {r['kind_label'] or r['label']}, with {name} filled in")
        for o in others:
            with ui.element("div").classes("ps-l-item"):
                ui.label(o["text"][0].upper() + o["text"][1:]).classes("font-medium")
                if o["what"] == "runs":
                    waiting_in = [x for x in studied if x["runs"] != x["analysed_runs"]]
                    ui.label(f"In {len(waiting_in)} of the experiments beside it. A run adds a projected value "
                             "once it is analysed.").classes("text-sm ps-muted")
        if isinstance(found, list) and not found and not others:
            ui.label("Nothing: every route physim can plan has a run of yours, and every run is analysed.").classes(
                "text-sm ps-muted p-3")

    with body:
        items()
    held["show_open"] = items.refresh


def _observables(ctx, view: dict, redraw, held: Optional[dict] = None) -> None:
    """All observables of the nucleus, one row per level or transition: evaluated, measured, theory and what the
    simulated runs reach in columns of their own, what is still open at the end; then the values that hold the
    evaluation only. A row opens its comparison (:mod:`.tab_compare`), with its routes. The list of those values is
    kept in ``held`` ("more"), for the link of the levels card that opens it."""
    ui = ctx.ui
    st = ctx.state["start"]
    display = st.get("unit")
    listed = [t for t in view["table"] if t["headline"]]
    rest = [t for t in view["table"] if not t["headline"]]
    outside = any(t["outside_program"] for t in listed)
    with ui.element("section").classes("ps-l-section ps-l-all").props('id="ps-all-observables"'):
        _, body = panel(ui, "All observables", "one row per level or transition; a row opens its comparison with "
                                               "the literature and its routes")
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
            more = ui.expansion(f"{len(rest)} more hold the evaluation only", icon="list").props("dense").classes(
                "w-full text-sm")
            if held is not None:
                held["more"] = more
            with more:
                rows = [{"id": f"{t['observable']}|{t['key']}", "label": label_text(t["label"]),
                         "value": cell_text(t["evaluated"], display) if t["evaluated"] else "—"} for t in rest]
                table = ui.table(columns=[{"name": "label", "label": "Observable", "field": "label", "align": "left"},
                                          {"name": "value", "label": "Evaluated", "field": "value",
                                           "align": "right"}],
                                 rows=rows, row_key="id", pagination=12).props("dense flat").classes("w-full")
                table.add_slot("body-cell-label", '<q-td :props="props" style="cursor: pointer; font-family: inherit">'
                                                  '<span class="ps-accent">{{ props.value }}</span></q-td>')
                table.on("rowClick", lambda e: open_observable(ctx, view["nuclide"],
                                                               *e.args[1]["id"].split("|", 1)))
