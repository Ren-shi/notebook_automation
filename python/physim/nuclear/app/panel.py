"""The setup panel: the checks, always visible, then the setup as a list of decisions (:mod:`physim.nuclear.workbench`),
each a card with its state and its consequence, its fields, and the rarely touched ones behind "more"."""

from __future__ import annotations

from types import SimpleNamespace

from .. import guide
from .. import workbench as wb

LEVEL_STYLE = {"error": "ps-bad font-semibold", "warning": "ps-warn", "note": "ps-muted"}
LEVEL_ICON = {"error": "error", "warning": "warning", "note": "info"}


def render_checks(ctx) -> None:
    ui = ctx.ui
    found = ctx.checks()
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
        elif group == "Ancillary detectors":
            draft = ctx.D().draft
            with ui.button("Add", icon="add").props("dense flat"):
                with ui.menu():
                    for key, (title, _, _, why) in wb.ANCILLARY_CARDS.items():
                        if key in draft:
                            continue
                        with ui.menu_item(on_click=lambda k=key: ctx.add_ancillary(k)).classes("max-w-sm"):
                            with ui.column().classes("gap-0"):
                                ui.label(title).classes("text-sm font-medium")
                                ui.label(why).classes("text-xs ps-muted")
                    with ui.menu_item(on_click=ctx.add_detector("forward")).classes("max-w-sm"):
                        with ui.column().classes("gap-0"):
                            ui.label("Beam monitor").classes("text-sm font-medium")
                            ui.label("a small silicon detector counting the elastically scattered beam").classes(
                                "text-xs ps-muted")
        elif group == "Lifetime measurement" and "lifetime" not in ctx.D().draft:
            with ui.button("Add", icon="add").props("dense flat"):
                with ui.menu():
                    for method, title, why in (
                            ("rdds", "Plunger (RDDS)", "a stopper (or degrader) foil at a few µm to mm behind the "
                                                       "target: lifetimes of about 1 ps to 1 ns"),
                            ("dsam", "Line shapes (DSAM)", "the recoil slows in the target's backing: lifetimes of "
                                                           "about 10 fs to 1 ps")):
                        with ui.menu_item(on_click=lambda m=method: ctx.add_lifetime(m)).classes("max-w-sm"):
                            with ui.column().classes("gap-0"):
                                ui.label(title).classes("text-sm font-medium")
                                ui.label(why).classes("text-xs ps-muted")
        elif group == "γ-ray detectors":
            with ui.button("Add", icon="add").props("dense flat"):
                with ui.menu():
                    for key, label, why, _, _ in guide.GAMMA_PLACEMENTS:
                        with ui.menu_item(on_click=ctx.add_gamma_detector(key)).classes("max-w-sm"):
                            with ui.column().classes("gap-0"):
                                ui.label(label).classes("text-sm font-medium")
                                ui.label(why).classes("text-xs ps-muted")


def _lifetime_fields(ctx, fields: list) -> None:
    """The lifetime card's text fields, two to a row: each applies on blur or Enter (``ctx.edit_lifetime``), is
    outlined while not applied, and follows the draft where it stands (``ctx.field_inputs``)."""
    ui = ctx.ui
    draft = ctx.D()._draft
    with ui.grid(columns=2).classes("w-full gap-1"):
        for section, field, label, placeholder in fields:
            text = "" if ctx.value_in(draft, section, field) is None else str(ctx.value_in(draft, section, field))
            inp = ui.input(label, value=text, placeholder=placeholder).props("dense outlined").classes("w-full")

            def apply(_=None, el=inp, s=section, f=field) -> None:
                ctx.edit_lifetime(s, f, el.value)

            inp.on("blur", apply)
            inp.on("keydown.enter", apply)
            ctx.field_inputs[(section, field)] = inp
            ctx.field_shown[(section, field)] = text
            ctx.mark(inp, section, field)


def deorientation_shown(lt: dict, chosen: str = "none") -> str:
    """Which form of the deorientation the card shows: the one its fields are in, else the one last chosen."""
    de = lt.get("deorientation") or {}
    if any(k in de for k in ("g", "c2", "c4")):
        return "g"
    if any(k in de for k in ("gamma2", "gamma4")):
        return "rates"
    return chosen if de or chosen != "none" else "none"


def _lifetime_body(ctx) -> None:
    """The lifetime card: the method; for a plunger its distances and offset, the foil after the gap (and a
    differential plunger's second foil) and the deorientation in vacuum; for DSAM the backing it reads. Then the
    checks about the lifetime measurement, and Remove."""
    ui = ctx.ui
    draft = ctx.D().draft
    lt = draft.get("lifetime") or {}
    method = lt.get("method", "rdds")
    sel = ui.select(wb.LIFETIME_METHODS, value=method, label="Method",
                    on_change=lambda e: ctx.edit_lifetime("lifetime", "method", e.value)).props(
        "dense outlined").classes("w-full")
    ctx.mark(sel, "lifetime", "method")
    if method == "dsam":
        ui.label("The recoil slows down and stops in the target's backing (or a thick target); the line shape of "
                 "each γ ray holds its lifetime. The backing is the target's own:").classes("text-xs ps-muted mt-1")
        ctx.section("", wb.BACKING_FIELDS, "backing", draft["target"].get("backing", {}))
    else:
        _lifetime_fields(ctx, wb.PLUNGER_FIELDS)
        ui.label("The distances as the plunger's dial gives them; the offset (the distance at electrical contact) "
                 "is added to each. Type them with a unit: 10, 20, 50 um.").classes("text-xs ps-muted")
        ui.label("Foil after the gap").classes("ps-section mt-3")
        kind = (lt.get("stopper") or {}).get("kind", "stopper")
        ks = ui.select(wb.FOIL_KINDS, value=kind, label="Kind",
                       on_change=lambda e: ctx.edit_lifetime("lifetime.stopper", "kind", e.value)).props(
            "dense outlined").classes("w-full")
        ctx.mark(ks, "lifetime.stopper", "kind")
        _lifetime_fields(ctx, wb.FOIL_FIELDS)
        if kind == "differential":
            ui.label("Second foil (stops the recoil)").classes("ps-section mt-3")
            _lifetime_fields(ctx, wb.DEGRADER2_FIELDS)
        ui.label("Deorientation in vacuum").classes("ps-section mt-3")
        form = deorientation_shown(lt, ctx.state.get("deorientation_form", "none"))
        ui.toggle({"none": "None", "rates": "Rates Γ₂, Γ₄", "g": "g factor with C₂, C₄"}, value=form,
                  on_change=lambda e: ctx.deorientation_form(e.value)).props("dense no-caps")
        if form != "none":
            _lifetime_fields(ctx, wb.DEORIENTATION_COMMON + (wb.DEORIENTATION_RATES if form == "rates"
                                                             else wb.DEORIENTATION_G))
            ui.label("G_k(t) = α_k + (1 − α_k) e^(−Γ_k t), with Γ_k = |g| / C_k in the second form (the static model "
                     "of Radeck et al., PRC 85, 014301).").classes("text-xs ps-muted")
        else:
            ui.label("The recoil keeps its alignment in flight. After Coulomb excitation a state living longer than "
                     "about 1 ps loses part of it: give α and Γ (or g and C) to correct the shifted and unshifted "
                     "intensities.").classes("text-xs ps-muted")
        with ui.expansion("More", icon="tune").classes("w-full").props("dense"):
            _lifetime_fields(ctx, wb.FOIL_MORE)
    checks = ui.refreshable(lambda: _lifetime_checks(ctx))
    checks()
    ctx.card_hooks.append(checks.refresh)  # the checks follow the draft as its fields are typed in
    ui.button("Remove", icon="delete", on_click=ctx.remove_lifetime).props("dense flat color=negative")


def lifetime_checks(ctx) -> list:
    """The checks the lifetime card shows: the draft's problems and the applied setup's warnings that are about
    the lifetime measurement, then the lifetime plan's own (once the draft is valid)."""
    found = [c for c in ctx.checks() if wb.lifetime_text(c["text"])]
    plan = wb.lifetime_plan(ctx.P()) if not ctx.D().problems else {}
    for level, key in (("error", "problems"), ("warning", "warnings")):
        for item in plan.get(key) or []:
            text = item if isinstance(item, str) else (item.get("text") if isinstance(item, dict) else str(item))
            if text and all(text != c["text"] for c in found):
                found.append({"level": level, "text": text})
    return found


def _lifetime_checks(ctx) -> None:
    ui = ctx.ui
    found = lifetime_checks(ctx)
    if found:
        with ui.column().classes("w-full ps-warnbox p-2 gap-1 mt-2"):
            for c in found:
                with ui.row().classes("items-start no-wrap gap-2"):
                    ui.icon(LEVEL_ICON[c["level"]]).classes(LEVEL_STYLE[c["level"]])
                    ui.label(c["text"]).classes(LEVEL_STYLE[c["level"]] + " text-sm")


def _card_body(ctx, d: wb.Decision) -> None:
    ui = ctx.ui
    if d.key == "lifetime":
        _lifetime_body(ctx)
        return
    draft = ctx.D().draft
    if d.key.startswith("detector:"):
        i = int(d.key.split(":")[1])
        sel = ui.select(list(wb.SHAPE_FIELDS), value=d.values.get("shape", "circle"), label="Shape",
                        on_change=lambda e, k=i: ctx.edit(f"detector {k + 1}", "shape", e.value)).props(
            "dense outlined").classes("w-full")
        ctx.help_icon(sel, "detector", "shape")
        ctx.mark(sel, f"detector {i + 1}", "shape")
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
    if d.key.startswith("anc:"):
        name = d.key.split(":", 1)[1]
        ui.button("Remove", icon="delete", on_click=lambda n=name: ctx.remove_ancillary(n)).props(
            "dense flat color=negative")
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
        ctx.checks_box()
        if p.running:
            with ui.row().classes("items-center gap-1 mt-1"):
                ui.icon("lock").classes("ps-warn")
                ui.label("A run is being taken: the setup is read-only until it ends.").classes("text-sm ps-warn")
    coulex = ctx.D().draft.get("reaction", {}).get("type") == "coulex"
    fusion = ctx.D().draft.get("reaction", {}).get("type") == "fusion-evaporation"
    shown_groups = set()
    holder = ui.column().classes("w-full gap-1")
    if p.running:
        holder.classes("opacity-60").style("pointer-events: none")
    with holder:
        cards = ctx.cards()
        groups_present = {d.group for d, _ in cards}
        for d, changed in cards:
            if d.group not in shown_groups:
                if d.group == "Run conditions" and (coulex or fusion) and "γ-ray detectors" not in groups_present:
                    _group_header(ctx, "γ-ray detectors")
                    ui.label("No γ-ray detectors yet: add one to see the γ ray of the excited state." if coulex else
                             "No γ-ray detectors yet: add the array to see the residues' γ rays.").classes(
                        "text-xs ps-muted")
                if d.group == "Run conditions" and fusion and "Ancillary detectors" not in groups_present:
                    _group_header(ctx, "Ancillary detectors")
                    ui.label("No ancillary detector: add a neutron wall, a charged-particle ball, a recoil separator "
                             "or a γ calorimeter to pick out the channels, or a beam monitor.").classes(
                        "text-xs ps-muted")
                if d.group == "Run conditions" and (coulex or fusion) and "Lifetime measurement" not in groups_present:
                    _group_header(ctx, "Lifetime measurement")
                    ui.label("No lifetime measurement: add a plunger (recoil distance) or line shapes in a backed "
                             "target (DSAM) to measure how long the states live.").classes("text-xs ps-muted")
                shown_groups.add(d.group)
                if d.group in ("Particle detectors", "γ-ray detectors", "Ancillary detectors", "Lifetime measurement"):
                    _group_header(ctx, d.group)
                else:
                    ui.element("div").classes("mt-2")
            exp = ui.expansion(value=ctx.state.get("open") == d.key).classes("w-full ps-card").props("dense")
            with exp.add_slot("header"):
                with ui.column().classes("gap-0 w-full"):
                    with ui.row().classes("items-center no-wrap gap-2"):
                        title = ui.label(d.title).classes("text-sm font-semibold")
                        mark = ui.label("not applied").classes("ps-pending-mark")
                        mark.set_visibility(changed)
                    state = ui.label(d.state).classes("text-sm")
                    with ui.row().classes("items-start no-wrap gap-1 w-full"):
                        consequence = ui.label(d.consequence).classes("text-xs ps-accent ps-consequence")
                        # The page brings these up to date where they stand after an edit (see sync_cards).
                        ctx.card_parts[d.key] = SimpleNamespace(title=title, state=state, consequence=consequence,
                                                                mark=mark)
                        ui.space()
                        ui.button(icon="science", on_click=lambda k=d.key: ctx.show_physics(wb.physics_for(k))).props(
                            "dense flat round size=xs").classes("ps-muted").on("click.stop", lambda: None).tooltip(
                            "The physics behind this line")
            # The fields of a card are built when it is opened: most cards stay closed, and the panel is rebuilt
            # after every change.
            ctx.when_opened(exp, lambda d=d: _card_body(ctx, d))
