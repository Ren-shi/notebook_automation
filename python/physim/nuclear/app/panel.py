"""The setup panel: the checks, always visible, then the setup as a list of decisions (:mod:`physim.nuclear.workbench`),
each a card with its state and its consequence, its fields, and the rarely touched ones behind "more"."""

from __future__ import annotations

from .. import guide
from .. import workbench as wb

LEVEL_STYLE = {"error": "ps-bad font-semibold", "warning": "ps-warn", "note": "ps-muted"}
LEVEL_ICON = {"error": "error", "warning": "warning", "note": "info"}


def render_checks(ctx) -> None:
    ui = ctx.ui
    found = wb.checks(ctx.P())
    if not found:
        with ui.row().classes("items-center gap-1"):
            ui.icon("check_circle").classes("ps-ok")
            ui.label("No warnings").classes("text-sm ps-ok")
        return
    with ui.card().classes("w-full ps-warnbox p-2 gap-1"):
        for c in found:
            with ui.row().classes("items-start no-wrap gap-2"):
                ui.icon(LEVEL_ICON[c["level"]]).classes(LEVEL_STYLE[c["level"]])
                ui.label(c["text"]).classes(LEVEL_STYLE[c["level"]] + " text-sm")


def _group_header(ctx, group: str) -> None:
    ui = ctx.ui
    with ui.row().classes("items-center mt-3 w-full"):
        ui.label(group).classes("ps-section")
        ui.space()
        if group == "Particle detectors":
            with ui.button("Add", icon="add").props("dense flat"):
                with ui.menu():
                    for pl in guide.PLACEMENTS:
                        with ui.menu_item(on_click=ctx.add_detector(pl.key)).classes("max-w-sm"):
                            with ui.column().classes("gap-0"):
                                ui.label(pl.label).classes("text-sm font-medium")
                                ui.label(pl.why).classes("text-xs ps-muted")
        elif group == "γ-ray detectors":
            with ui.button("Add", icon="add").props("dense flat"):
                with ui.menu():
                    for key, label, why, _, _ in guide.GAMMA_PLACEMENTS:
                        with ui.menu_item(on_click=ctx.add_gamma_detector(key)).classes("max-w-sm"):
                            with ui.column().classes("gap-0"):
                                ui.label(label).classes("text-sm font-medium")
                                ui.label(why).classes("text-xs ps-muted")


def _card_body(ctx, d: wb.Decision) -> None:
    ui = ctx.ui
    draft = ctx.P().draft
    if d.key.startswith("detector:"):
        i = int(d.key.split(":")[1])
        sel = ui.select(list(wb.SHAPE_FIELDS), value=d.values.get("shape", "circle"), label="Shape",
                        on_change=lambda e, k=i: ctx.edit(f"detector {k + 1}", "shape", e.value)).props(
            "dense outlined").classes("w-full")
        ctx.help_icon(sel, "detector", "shape")
    ctx.section("", d.fields, d.section, d.values)
    if d.key == "target":
        ctx.reaction_section(draft.get("reaction", {"type": "elastic"}), title="")
    if d.key.startswith("gamma:") and wb._crystals(d.values) == 4:
        i = int(d.key.split(":")[1])
        ctx.clover_switches(i, d.values, numbers=False)
    more = list(d.more)
    if more or d.key == "target":
        with ui.expansion("More", icon="tune").classes("w-full").props("dense"):
            if more:
                ctx.section("", more, d.section, d.values)
            if d.key == "target":
                ctx.section("Backing (optional)", wb.BACKING_FIELDS, "backing", draft["target"].get("backing", {}))
    if d.key.startswith(("detector:", "gamma:")):
        i = int(d.key.split(":")[1])
        gamma = d.key.startswith("gamma:")
        with ui.row():
            ui.button("Duplicate", icon="content_copy",
                      on_click=(ctx.duplicate_gamma_detector if gamma else ctx.duplicate_detector)(i)).props(
                "dense flat")
            ui.button("Remove", icon="delete",
                      on_click=(ctx.remove_gamma_detector if gamma else ctx.remove_detector)(i)).props(
                "dense flat color=negative")


def render(ctx) -> None:
    """The whole panel: checks, then the decision cards in order. Read-only while a run is being taken."""
    ui = ctx.ui
    p = ctx.P()
    with ui.column().classes("w-full gap-1"):
        render_checks(ctx)
        if p.running:
            with ui.row().classes("items-center gap-1 mt-1"):
                ui.icon("lock").classes("ps-warn")
                ui.label("A run is being taken: the setup is read-only until it ends.").classes("text-sm ps-warn")
    coulex = p.draft.get("reaction", {}).get("type") == "coulex"
    shown_groups = set()
    holder = ui.column().classes("w-full gap-1")
    if p.running:
        holder.classes("opacity-60").style("pointer-events: none")
    with holder:
        decisions = wb.decisions(p)
        groups_present = {d.group for d in decisions}
        for d in decisions:
            if d.group not in shown_groups:
                if d.group == "Run conditions" and coulex and "γ-ray detectors" not in groups_present:
                    _group_header(ctx, "γ-ray detectors")
                    ui.label("No γ-ray detectors yet: add one to see the γ ray of the excited state.").classes(
                        "text-xs ps-muted")
                shown_groups.add(d.group)
                if d.group in ("Particle detectors", "γ-ray detectors"):
                    _group_header(ctx, d.group)
                else:
                    ui.element("div").classes("mt-2")
            exp = ui.expansion(value=ctx.state.get("open") == d.key).classes("w-full ps-card").props("dense")
            with exp.add_slot("header"):
                with ui.column().classes("gap-0 w-full"):
                    ui.label(d.title).classes("text-sm font-semibold")
                    ui.label(d.state).classes("text-sm")
                    with ui.row().classes("items-start no-wrap gap-1 w-full"):
                        ui.label(d.consequence).classes("text-xs ps-accent ps-consequence")
                        ui.space()
                        ui.button(icon="science", on_click=lambda k=d.key: ctx.show_physics(wb.physics_for(k))).props(
                            "dense flat round size=xs").classes("ps-muted").on("click.stop", lambda: None).tooltip(
                            "The physics behind this line")
            with exp:
                _card_body(ctx, d)
