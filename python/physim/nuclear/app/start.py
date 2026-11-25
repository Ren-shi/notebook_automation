"""The start views, shown in place of an experiment's steps and without the rail: the page of one nucleus
(:mod:`.tab_nucleus`), the comparison of one of its observables (:mod:`.tab_compare`), the nuclei, each with its
experiments, what is known about it and what is still open (:func:`nuclei_page`, from
:func:`physim.nuclear.nucleus.nuclei`), the kinds of experiment (:func:`kinds_page`) and the Physics reference
(:func:`physics_head`; the page shows the reference itself under it). Above each of them stands the same bar
(:func:`top_bar`). The app opens on a nucleus (:func:`opening`).

"Start from a kind of experiment" is the start page as it was (:func:`kinds_page`): the kind of experiment
(:data:`physim.nuclear.workbench.KINDS`), then its beam, target and detectors, or an experiment already made.
Coulomb excitation and elastic scattering then go on to the steps as they always were; fusion-evaporation to its own
screens (:mod:`.tab_fusion`); lifetimes to a Coulomb-excitation or fusion-evaporation setup with a plunger or DSAM
(:mod:`.tab_lifetimes`). With no experiment and no nucleus yet, the kinds come first."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Optional

from .. import workbench as wb
from . import nucleus_parts as parts
from . import tab_compare, tab_nucleus


#: How many nuclei the switcher of the bar shows: the most recently changed ones, the others are in "All nuclei".
SWITCHER = 5

#: The bar above the start views: over the whole width of the page, wrapping in a narrow window.
TOP_STYLE = """
.ps-top { display: flex; flex-wrap: wrap; align-items: center; gap: 8px 16px; padding: 8px 24px;
  background: var(--ps-sunk); border-bottom: 1px solid var(--ps-line);
  margin: calc(-1 * var(--nicegui-default-padding, 1rem)); margin-bottom: 8px;
  width: calc(100% + 2 * var(--nicegui-default-padding, 1rem)); }
.ps-top-part { display: flex; flex-wrap: wrap; align-items: center; gap: 6px; min-width: 0; }
.ps-top-btn { border: 1px solid var(--ps-line); background: var(--ps-surface); color: var(--ps-ink);
  min-height: 38px; padding: 0 12px; }
.ps-top-btn.ps-top-dashed { border-style: dashed; background: transparent; color: var(--ps-muted); }
.ps-top-btn.ps-top-on { background: var(--ps-fill); border-color: var(--ps-fill); color: var(--ps-on-fill);
  font-weight: 600; }
.ps-top-nuclide { font-size: 16px; }
.ps-top-back { max-width: 340px; }
.ps-top-back .q-btn__content { flex-wrap: nowrap; }
.ps-top-back .q-btn__content .block { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
"""


def experiment_kind(path: str) -> str:
    """The kind of an experiment folder, read from its setup's reaction without building it ("" if unreadable)."""
    from .. import _toml

    try:
        d = _toml.loads((Path(path) / "experiment.toml").read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001 -- an unreadable folder shows no kind
        return ""
    r = d.get("reaction", {})
    kind = r.get("type", "elastic")
    lt = d.get("lifetime")
    if isinstance(lt, dict) and kind in ("coulex", "fusion-evaporation"):
        how = "plunger" if lt.get("method", "rdds") == "rdds" else "DSAM"
        return f"Lifetimes · {how} after " + ("fusion" if kind == "fusion-evaporation" else "Coulex")
    if kind == "fusion-evaporation":
        return "Fusion-evaporation"
    if kind == "coulex":
        return "Coulex · " + ("beam" if r.get("excite") in ("projectile", "beam") else "target")
    return "Elastic"


def fusion_preview(beam: str, energy: str, target: str) -> str:
    """One line about the compound nucleus a fusion-evaporation setup would make: E*, the highest spin that fuses,
    the beam against the Bass barrier and the xn channel the excitation favours ("" if it cannot be worked out)."""
    try:
        from ..exit_channels import fusion_summary
        from ..experiment import main_nuclide
        from ..rates import beam_ion
        from ..residues import likely_xn

        exp = wb.new_setup(beam, energy, target, "1 mg/cm2", "fusion-evaporation", detectors_from="fusion-blank")
        f = fusion_summary(exp)
        k = likely_xn(beam_ion(exp), main_nuclide(exp.target), exp.beam.energy_mev)
        lo, hi = f["excitation_range_mev"]
        above = exp.beam.energy_mev - f["barrier_lab_mev"]
        return (f"{f['compound']}* at E* ≈ {(lo + hi) / 2:.1f} MeV, spins up to {f['max_l']} ħ. The beam is "
                f"{abs(above):.1f} MeV {'above' if above >= 0 else 'below'} the fusion barrier (Bass, "
                f"{f['barrier_lab_mev']:.1f} MeV lab)." + (f" Most likely {k}n." if k >= 1 else ""))
    except Exception:  # noqa: BLE001 -- not a complete beam and target yet
        return ""


def anything_yet() -> bool:
    """Whether there is an experiment or a nucleus with rows of its own: then the nuclei come first."""
    from .. import results
    from .. import runs as _runs

    if _runs.list_experiments():
        return True
    home = results.nuclei_home()
    return home.is_dir() and any((p / results.POINTS_FILE).is_file() for p in home.iterdir())


def last_changed(n: dict, modified: dict) -> float:
    """When a nucleus of :func:`physim.nuclear.nucleus.nuclei` last changed: its data points or an experiment that
    studies it (``modified``: when each experiment changed, by its path)."""
    from .. import results

    points = results.nucleus_folder(n["nuclide"]) / results.POINTS_FILE
    return max([points.stat().st_mtime if points.is_file() else 0.0]
               + [modified.get(x["path"]) or 0.0 for x in n["experiments"] if x["role"] != "partner"])


def switched_between(found: dict) -> list:
    """The nuclides the bar's switcher offers, from :func:`.nucleus_parts.load_nuclei`: the nuclei that get a card
    (:func:`is_card`), the :data:`SWITCHER` most recently changed ones, in the order of the nuclei screen."""
    cards = [n for n in found["nuclei"] if is_card(n)]
    recent = sorted(cards, key=lambda n: -last_changed(n, found["modified"]))[:SWITCHER]
    return [n["nuclide"] for n in cards if n in recent]


def opening(remembered: Optional[str] = None) -> dict:
    """What the app opens on when its address names nothing, as ``state["start"]`` takes it: {"view", "nucleus"}.
    The nucleus opened last (``remembered``, by the browser); if there is none, or it is not here any more, the one
    studied by the most recently changed experiment; without an experiment that studies one, the most recently
    changed nucleus that holds data points; otherwise the nuclei, or with nothing at all the kinds of experiment
    (:func:`anything_yet`). Reads the experiments and the ledgers: for a worker thread."""
    from .. import results

    found = parts.load_nuclei(fill=False)
    try:
        name = results.nuclide_name(str(remembered)) if remembered else None
    except (ValueError, KeyError, TypeError, IndexError):
        name = None
    if name is None or name not in [n["nuclide"] for n in found["nuclei"]]:
        studied = {x["path"]: n["nuclide"] for n in found["nuclei"] for x in n["experiments"]
                   if x["role"] == "studied"}
        cards = [n for n in found["nuclei"] if is_card(n)]
        if studied:
            name = studied[max(studied, key=lambda path: found["modified"].get(path) or 0.0)]
        elif cards:
            name = max(cards, key=lambda n: last_changed(n, found["modified"]))["nuclide"]
        else:
            name = None
    if name is not None:
        return {"view": "nucleus", "nucleus": name}
    return {"view": "nuclei" if anything_yet() else "kinds"}


def render(ctx) -> None:
    """The start view that ``ctx.state["start"]["view"]`` names, into the current container: one nucleus
    ("nucleus", :func:`.tab_nucleus.render`), one observable of it ("observable", :func:`.tab_compare.render`), the
    nuclei ("nuclei"), the kinds of experiment ("kinds"), the head of the Physics reference ("physics"), or the
    waiting line until it is known what the app opens on ("opening").

    ``ctx`` gives ``ui``, ``state``, ``create_experiment``, ``open_experiment``, ``open_example``, ``import_setup``,
    ``refresh_start`` (show the start page, built again), ``close_start`` (the experiment's steps in its place),
    ``experiment_name`` (the name of the open experiment, None when none is open), ``remember_nucleus`` (in the
    browser, for the next opening), ``open_physics`` and ``theme_switch`` (builds a Light/Dark switch)."""
    st = ctx.state.setdefault("start", {"kind": None})
    if st.get("view") is None:
        st["view"] = "nuclei" if anything_yet() else "kinds"
    if st["view"] in ("nucleus", "observable") and st.get("nucleus"):
        ctx.remember_nucleus(st["nucleus"])
        (tab_nucleus if st["view"] == "nucleus" else tab_compare).render(ctx)
    elif st["view"] == "opening":
        parts.waiting(ctx.ui, "Opening the nucleus you were on…")
    elif st["view"] == "physics":
        physics_head(ctx)
    elif st["view"] == "kinds":
        kinds_page(ctx)
    else:
        nuclei_page(ctx)


# ---------------------------------------------------------------------------------------------------------------
# The bar above every start view, and the Physics reference as one


def top_bar(ctx):
    """The bar above the start views, into the current container: the switcher of the nuclei, the way back to the
    open experiment, the Physics reference and the theme. Returns the function that brings it up to date with the
    view shown (and reads the nuclei again, on a worker thread)."""
    ui = ctx.ui
    held: dict = {"nuclei": []}
    with ui.element("header").classes("ps-top") as bar:
        ui.label("physim").classes("ps-wordmark")
        switcher = ui.element("nav").classes("ps-top-part").props('aria-label="Nuclei"')
        with ui.button("Another nucleus", icon="add", color=None).props("flat no-caps").classes(
                "ps-top-btn ps-top-dashed"):
            with ui.menu() as menu, ui.row().classes("items-center gap-2 p-2"):
                other = ui.input("Another nucleus", placeholder="58Ni").props("dense outlined autofocus").style(
                    "width: 160px")
                ui.button("Open its page", icon="arrow_forward", on_click=lambda: another()).props(
                    "flat no-caps dense").tooltip("The page of any nucleus: ENSDF's adopted values, ready for its "
                                                  "first data point or experiment")
        ui.element("span").style("flex: 1 1 auto")
        right = ui.element("div").classes("ps-top-part")
        ctx.theme_switch()

    def another() -> None:
        menu.close()
        name, other.value = other.value, ""
        parts.open_nucleus(ctx, name)

    other.on("keydown.enter", another)

    def build() -> None:
        st = ctx.state["start"]
        view = st.get("view")
        here = st.get("nucleus") if view in ("nucleus", "observable") else None
        switcher.clear()
        with switcher:
            for name in held["nuclei"] + ([here] if here and here not in held["nuclei"] else []):
                button = ui.button(parts.nuclide_text(name), color=None,
                                   on_click=lambda name=name: parts.open_nucleus(ctx, name)).props("flat no-caps")
                button.classes("ps-top-btn ps-top-nuclide" + (" ps-top-on" if name == here else ""))
                if name == here:
                    button.props('aria-current="page"')
            ui.button("All nuclei", color=None, on_click=lambda: parts.open_nuclei(ctx)).props(
                "flat no-caps").classes("ps-top-btn" + (" ps-top-on" if view == "nuclei" else "")).tooltip(
                "Every nucleus with its experiments, what is known about it and what is still open")
        right.clear()
        with right:
            name = ctx.experiment_name()
            if name is not None:
                ui.button(f"Back to {name}", icon="arrow_back", color=None, on_click=lambda: parts.close(ctx)).props(
                    "flat no-caps").classes("ps-top-btn ps-top-back").tooltip("The steps of the open experiment")
            ui.button("Physics", icon="menu_book", color=None, on_click=ctx.open_physics).props(
                "flat no-caps").classes("ps-top-btn" + (" ps-top-on" if view == "physics" else "")).tooltip(
                "The physics behind the numbers, and where each part is validated")

    async def read() -> None:
        from nicegui import run

        try:
            found = await run.io_bound(lambda: switched_between(parts.load_nuclei(fill=False)))
        except Exception:  # noqa: BLE001 -- the nuclei screen says what is wrong; the bar keeps what it has
            return
        if found is not None and not switcher.is_deleted and found != held["nuclei"]:
            held["nuclei"] = found
            build()

    def refresh() -> None:
        build()
        with bar:  # what asked for it may be gone with the view it was on
            parts.later(ui, read)

    return refresh


def physics_head(ctx) -> None:
    """The head of the Physics reference as a start view: the way back to where one was, and whose numbers the
    reference shows."""
    ui = ctx.ui
    st = ctx.state["start"]
    nucleus = parts.nuclide_text(st.get("nucleus") or "")
    back = {"nucleus": nucleus, "observable": f"the comparison for {nucleus}",
            "kinds": "the kinds of experiment"}.get(st.get("before_physics"), "All nuclei")

    def go_back() -> None:
        st["view"] = st.get("before_physics") or "nuclei"
        parts.show_from_top(ctx)

    with ui.row().classes("w-full items-center gap-3"):
        ui.button(f"← {back}", on_click=go_back).props("flat no-caps dense").classes("ps-accent")
        name = ctx.experiment_name()
        ui.label("Physics").classes("text-xl font-semibold")
        ui.label(f"with the numbers of {name}" if name is not None else
                 "with the numbers of a first example: no experiment is open").classes("text-sm ps-muted")


# ---------------------------------------------------------------------------------------------------------------
# The nuclei


def is_card(n: dict) -> bool:
    """Whether a nucleus of :func:`physim.nuclear.nucleus.nuclei` gets a card of its own: it is studied by an
    experiment or has rows besides the evaluation. One that only appears as the partner of an experiment, or whose
    page was only looked at (which takes its evaluation), goes in the line under the cards."""
    return n["rows"] > n["evaluation"]["rows"] or any(x["role"] != "partner" for x in n["experiments"])


def nuclei_page(ctx) -> None:
    """The nuclei: one card per nucleus that is studied or has rows, with what is known and projected, what is
    still open and its experiments; the partner-only nuclei in a line; the experiments tied to no nucleus in the
    table they always had. Read on a worker thread (:func:`.nucleus_parts.load_nuclei`), which also makes the rows
    of experiments from before the ledger."""
    ui = ctx.ui
    parts.add_style(ctx)
    note = ctx.state["start"].pop("note", None)  # said once: why the address's experiment was not opened
    with ui.column().classes("w-full gap-3").style("max-width: 1100px; margin: 0 auto"):
        if note:
            ui.label(note).classes("text-sm ps-warnbox rounded p-3 w-full ps-start-note")
        with ui.row().classes("w-full items-end gap-2"):
            with ui.column().classes("gap-1"):
                ui.label("What are you studying?").classes("text-xl font-semibold")
                ui.label("Each nucleus holds its experiments and everything known about it, by method. Open a "
                         "nucleus to compare the routes to an observable and to see which experiments are still "
                         "worth planning.").classes("text-sm ps-muted").style("max-width: 78ch")
            ui.space()
            ui.button("Start from a kind of experiment", on_click=lambda: parts.plan(ctx, None)).props(
                "outline no-caps")
        with ui.row().classes("w-full items-center gap-2"):
            other = ui.input("Another nucleus", placeholder="58Ni").props("dense outlined").style("width: 180px")
            other.on("keydown.enter", lambda: parts.open_nucleus(ctx, other.value))
            ui.button("Open its page", icon="arrow_forward",
                      on_click=lambda: parts.open_nucleus(ctx, other.value)).props(
                "flat no-caps dense").tooltip("The page of any nucleus: ENSDF's adopted values, ready for its first "
                                              "data point or experiment")
        body = ui.column().classes("w-full gap-3")
    with body:
        parts.waiting(ui, "Reading the nuclei and their experiments…")

    async def fill() -> None:
        from nicegui import run

        try:
            found = await run.io_bound(parts.load_nuclei)
        except Exception as err:  # noqa: BLE001 -- say why, and keep the experiments reachable
            found = {"error": str(err)}
        if found is None or body.is_deleted:
            return
        body.clear()
        with body:
            if "error" in found:
                ui.label(f"The nuclei could not be read: {found['error']}").classes(
                    "text-sm ps-warnbox rounded p-3 w-full")
                experiments_table(ctx)
                return
            _nuclei(ctx, found)

    parts.later(ui, fill)


def _nuclei(ctx, found: dict) -> None:
    ui = ctx.ui
    display = ctx.state["start"].get("unit")
    cards = [n for n in found["nuclei"] if is_card(n)]
    partners = [n for n in found["nuclei"] if not is_card(n)]
    if not cards and not partners and not found["unassigned"]:
        ui.label("No nucleus yet. Open the page of a nucleus above, or start from a kind of experiment."
                 ).classes("text-sm ps-note p-3 w-full")
    for n in cards:
        name, title = n["nuclide"], parts.nuclide_text(n["nuclide"])
        studied = [x for x in n["experiments"] if x["role"] != "partner"]
        with ui.element("div").classes("ps-panel"):
            with ui.element("div").classes("ps-nuc-head"):
                ui.button(f"{title} →", on_click=lambda name=name: parts.open_nucleus(ctx, name)).props(
                    "flat no-caps dense").classes("text-lg font-semibold ps-open").tooltip(
                    "Everything known and projected about it, by method")
                evaluation = n["evaluation"]
                ui.label(" · ".join((
                    f"{len(studied)} experiment{'s' * (len(studied) != 1)}",
                    f"{n['analysed_runs']} analysed run{'s' * (n['analysed_runs'] != 1)}" if n["analysed_runs"]
                    else "no analysed run",
                    f"{n['points']} data point{'s' * (n['points'] != 1)}" if n["points"] else "no data points yet",
                    f"evaluation: {evaluation['reference']}" if evaluation["taken"] else "evaluation not taken yet",
                ))).classes("text-sm ps-muted")
                ui.space()
                ui.button(f"New experiment for {title}", icon="add",
                          on_click=lambda name=name: parts.plan(ctx, name)).props("flat no-caps dense")
            with ui.element("div").classes("ps-nuc-body"):
                with ui.element("div"):
                    ui.label("WHAT IS KNOWN AND PROJECTED").classes("ps-nuc-h")
                    if not n["table"]:
                        ui.label("The evaluation, on its page." if evaluation["taken"] else "Nothing yet.").classes(
                            "text-sm ps-muted")
                    with ui.element("div").classes("ps-kv"):
                        for t in n["table"][:parts.HEADLINES]:
                            row = ui.label(parts.label_text(t["label"])).style("cursor: pointer")
                            row.on("click", lambda name=name, t=t: parts.open_nucleus(
                                ctx, name, t["observable"], t["key"]))
                            ui.label(parts.known_lines(t, display)).classes("ps-num")
                    if len(n["table"]) > parts.HEADLINES:
                        ui.label(f"… and {len(n['table']) - parts.HEADLINES} more on its page.").classes(
                            "text-xs ps-muted")
                    ui.label("STILL OPEN").classes("ps-nuc-h").style("margin-top: 10px")
                    ui.label(" · ".join(parts.label_text(o["text"]) for o in n["open"]) or "nothing").classes(
                        "text-sm " + ("ps-open" if n["open"] else "ps-muted")).style("line-height: 1.6")
                with ui.element("div"):
                    ui.label("EXPERIMENTS").classes("ps-nuc-h")
                    if n["experiments"]:
                        parts.experiments_table(ctx, n["experiments"], display)
                    else:
                        ui.label("No experiment yet: data points only.").classes("text-sm ps-muted")
    if partners:
        with ui.row().classes("w-full items-center gap-1 ps-note").style("padding: 6px 12px"):
            ui.label("Only as the partner of an experiment, or looked up:" if any(
                not n["experiments"] for n in partners) else "Only as the partner of an experiment:").classes(
                "text-sm ps-muted")
            for n in partners:
                count = len(n["experiments"])
                ui.button(f"{parts.nuclide_text(n['nuclide'])} · " + (
                    f"in {count} experiment{'s' * (count != 1)}" if count else "the evaluation only"),
                          on_click=lambda name=n["nuclide"]: parts.open_nucleus(ctx, name)).props(
                    "flat no-caps dense").classes("text-sm")
    if found["unassigned"]:
        modified = found.get("modified", {})
        listed = [{"name": u["name"], "path": u["path"], "runs": u["runs"], "kind": u["kind_label"],
                   "modified": modified.get(u["path"]), "reason": u["reason"]} for u in found["unassigned"]]
        experiments_table(ctx, sorted(listed, key=lambda x: -(x["modified"] or 0)),
                          "Experiments tied to no nucleus")


# ---------------------------------------------------------------------------------------------------------------
# The kinds of experiment: the start page as it was


def kinds_page(ctx) -> None:
    """The kinds of experiment, the form of the kind picked, and the experiments already made. When a nucleus
    asked for the experiment (``state["start"]["for_nucleus"]``, from "Plan it" or "New experiment for" a
    nucleus), the form has it filled in where it has a field for it."""
    ui = ctx.ui
    st = ctx.state.setdefault("start", {"kind": None})
    nuclide = st.get("for_nucleus")
    with ui.column().classes("w-full gap-3").style("max-width: 1060px; margin: 0 auto"):
        with ui.row().classes("w-full items-end"):
            with ui.column().classes("gap-1"):
                ui.label("What experiment are you planning?").classes("text-xl font-semibold")
                ui.label("Pick the kind of experiment first: it decides which detectors the setup offers, what the "
                         "plan computes and which spectra you see. The beam, target and detectors can change "
                         "afterwards.").classes("text-sm ps-muted").style("max-width: 70ch")
            ui.space()
            ui.button("The nuclei", icon="arrow_back", on_click=lambda: leave(None)).props("outline no-caps").tooltip(
                "The nuclei, each with its experiments and what is known about it")
        if nuclide:
            with ui.row().classes("w-full items-center gap-2 ps-note").style("padding: 6px 12px"):
                ui.label(f"An experiment for {parts.nuclide_text(nuclide)}: pick its kind. Where the form has a field "
                         f"for the nucleus that is excited, {nuclide} is filled in.").classes("text-sm")
                ui.space()
                ui.button(f"Back to {parts.nuclide_text(nuclide)}", on_click=lambda: leave(nuclide)).props(
                    "flat no-caps dense")
        with ui.element("div").classes("w-full").style(
                "display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 10px"):
            for k in wb.KINDS:
                if not k.live:
                    continue
                on = st["kind"] == k.key
                card = ui.element("div").classes("ps-card ps-kind" + (" ps-kind-on" if on else "")).style(
                    "padding: 14px 16px; cursor: pointer; display: flex; flex-direction: column; gap: 8px")
                card.on("click", lambda key=k.key: pick(key))
                with card:
                    ui.label(k.label).classes("text-base font-semibold")
                    ui.label(k.what).classes("text-sm ps-muted")
                    ui.label("→ " + k.gets).classes("text-xs ps-num ps-muted").style(
                        "margin-top: auto; border-top: 1px solid var(--ps-hair); padding-top: 8px")
        with ui.element("div").classes("w-full").style(
                "display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 10px"):
            for k in wb.KINDS:
                if k.live:
                    continue
                with ui.element("div").classes("ps-card").style(
                        "padding: 10px 16px; border-style: dashed; opacity: 0.6"):
                    ui.label(k.label).classes("text-sm font-semibold")
                    ui.label(k.what + " Later.").classes("text-xs ps-muted")
        form = ui.column().classes("w-full")
        with form:
            if st["kind"]:
                kind_form(ctx, next(k for k in wb.KINDS if k.key == st["kind"]))
        experiments_table(ctx)

    def pick(key: str) -> None:
        st["kind"] = key
        st.pop("values", None)
        ctx.refresh_start()

    def leave(to) -> None:
        st.update(kind=None, for_nucleus=None)
        st.pop("values", None)
        parts.open_nucleus(ctx, to) if to else parts.open_nuclei(ctx)


def filled_for(k, nuclide) -> dict:
    """What the form of the kind ``k`` starts with for an experiment that a nucleus asked for: the nucleus as the
    one that is excited, where the form has that field (as the beam, or as the target where it is the kind's own
    target). Everything else is the kind's own default; nothing is computed. Empty for a kind whose form has no
    such field, and without a nucleus."""
    if not nuclide or k.key not in ("coulex", "lifetimes"):
        return {}
    if nuclide == k.defaults[2]:
        return {"target": nuclide, "excite": "target"}
    return {"beam": nuclide, "excite": "beam"}


def kind_form(ctx, k) -> None:
    ui = ctx.ui
    st = ctx.state["start"]
    nuclide = st.get("for_nucleus")
    filled = filled_for(k, nuclide)
    values = st.setdefault("values", {**dict(zip(("beam", "energy", "target", "thickness"), k.defaults),
                                             template=k.templates[0], name="",
                                             excite="beam" if k.key == "lifetimes" else "target", fill=True),
                                      **filled})
    with ui.column().classes("w-full ps-card gap-0"):
        with ui.row().classes("w-full items-center ps-panel-head"):
            ui.label(k.label)
            ui.label("the beam, the target and where the detectors come from").classes("ps-muted text-sm")
            ui.space()
            ui.button("Change the kind", on_click=lambda: (st.update(kind=None), ctx.refresh_start())).props(
                "flat no-caps dense")
        if nuclide:
            as_, other = ("target", "beam") if filled.get("excite") == "target" else ("beam", "target")
            ui.label((f"{nuclide} is filled in as the {as_}, the nucleus that is excited. The {other}, the energy and "
                      f"the target's thickness are this kind's usual ones: set them for {nuclide} (the Setup step "
                      f"says whether the energy is safe). To excite it as the {other} instead, swap the two."
                      if filled else
                      f"This form has no field for the nucleus studied: choose the beam and target that reach "
                      f"{nuclide}.")).classes("text-sm ps-note p-2").style("margin: 10px 12px 0")
        with ui.element("div").classes("w-full p-3").style(
                "display: grid; grid-template-columns: minmax(0, 1.25fr) minmax(0, 1fr); gap: 14px"):
            with ui.column().classes("gap-2 min-w-0"):
                ui.label("Beam and target").classes("ps-section")
                with ui.grid(columns=2).classes("w-full gap-1"):
                    for key, label in (("beam", "Beam"), ("energy", "Energy"), ("target", "Target"),
                                       ("thickness", "Thickness")):
                        inp = ui.input(label, value=values[key]).props("dense outlined").classes("w-full")
                        inp.on_value_change(lambda e, key=key: values.update({key: e.value}))
                        if k.key == "fusion-evaporation":
                            inp.on("blur", lambda: preview.refresh())
                if k.key in ("coulex", "lifetimes"):
                    ui.select({"target": "The target", "beam": "The beam"}, value=values["excite"],
                              label="Which nucleus is excited" + (" (after Coulomb excitation)"
                                                                  if k.key == "lifetimes" else ""),
                              on_change=lambda e: values.update(excite=e.value)).props("dense outlined").classes(
                        "w-full")
                if k.key == "lifetimes":
                    ui.toggle({True: "Fill the exit channels from physim", False: "Enter them by hand"},
                              value=values["fill"], on_change=lambda e: values.update(fill=e.value)).props(
                        "dense no-caps")
                    ui.label("The excited nucleus is for the Coulomb-excitation template, the exit channels for "
                             "the fusion-evaporation ones. The plunger's stopper is made thick enough to stop this "
                             "beam and target's recoil; DSAM takes the template's backing. The worked examples "
                             "are the published measurements they are named after.").classes("text-xs ps-muted")
                if k.key == "fusion-evaporation":
                    @ui.refreshable
                    def preview() -> None:
                        text = fusion_preview(values["beam"], values["energy"], values["target"])
                        if text:
                            ui.label(text).classes("text-sm ps-note p-2 w-full")

                    preview()
                    ui.label("Exit channels").classes("ps-section")
                    ui.toggle({True: "Fill from physim", False: "Enter by hand"}, value=values["fill"],
                              on_change=lambda e: values.update(fill=e.value)).props("dense no-caps")
                    ui.label("Fill from physim runs physim's own fusion and evaporation (an estimate; PACE4 or a "
                             "measurement is the better source). Importing a PACE4 or GEMINI++ run comes later. The "
                             "table stays editable in the setup.").classes("text-xs ps-muted")
                ui.input("Name (made from the beam and target if empty)", value=values["name"],
                         on_change=lambda e: values.update(name=e.value)).props("dense outlined").classes("w-full")
            with ui.column().classes("gap-2 min-w-0"):
                ui.label("Detectors from").classes("ps-section")
                ui.radio({t: f"{wb.template_label(t)[0]} — {wb.template_label(t)[1]}" for t in k.templates},
                         value=values["template"], on_change=lambda e: values.update(template=e.value)).props(
                    "dense").classes("text-sm")
                if k.examples:
                    ui.label("Or start from a worked example").classes("ps-section")
                    with ui.row().classes("gap-1").style("flex-wrap: wrap"):
                        for ex in k.examples:
                            ui.button(ex, on_click=lambda ex=ex: ctx.open_example(ex)).props(
                                "dense outline no-caps rounded").classes("text-xs")
        with ui.row().classes("w-full items-center p-3").style("border-top: 1px solid var(--ps-line)"):
            ui.label("Next: the Setup step" + (", with the fusion-evaporation scene." if k.key == "fusion-evaporation"
                                                else ", with the lifetime measurement's card." if k.key == "lifetimes"
                                                else ".")).classes("text-sm ps-muted")
            ui.space()
            ui.button("Make the experiment", icon="arrow_forward", on_click=lambda: make()).props(
                "unelevated no-caps")

    async def make() -> None:
        measure = k.measure
        coulex = k.key == "coulex"
        if k.key == "lifetimes":  # the template says which reaction the lifetime measurement follows
            reaction = next(t[4] for t in wb.LIFETIME_TEMPLATES if t[0] == values["template"])
            coulex = reaction == "coulex"
            measure = "fusion-evaporation" if not coulex else measure
        if coulex:
            measure = "coulex-beam" if values["excite"] == "beam" else "coulex-target"
        await ctx.create_experiment(values["beam"], values["energy"], values["target"], values["thickness"], measure,
                                    values["template"], values["name"],
                                    fill=measure == "fusion-evaporation" and values["fill"])


def experiments_table(ctx, found=None, title: str = "Or open an experiment") -> None:
    """The experiments as a table: name, kind, runs, when last changed; a row opens the experiment. ``found``
    lists them ({"name", "path", "runs", "modified"}, and "kind" and "reason" where known: an experiment tied to
    no nucleus says why); every experiment when left out."""
    ui = ctx.ui
    from .. import runs as _runs

    if found is None:
        found = _runs.list_experiments()
    reasons = any(x.get("reason") for x in found)
    with ui.column().classes("w-full ps-card gap-0"):
        with ui.row().classes("w-full items-center ps-panel-head"):
            ui.label(title)
            ui.label(f"in {_runs.home()}").classes("ps-muted text-sm")
            ui.space()
            ui.button("Import a setup file…", on_click=ctx.import_setup).props("flat no-caps dense")
        if not found:
            ui.label("No experiment yet.").classes("text-sm p-3 ps-muted")
            return
        rows = [{"name": x["name"], "kind": x.get("kind") or experiment_kind(x["path"]), "runs": x["runs"],
                 "changed": time.strftime("%Y-%m-%d %H:%M", time.localtime(x["modified"])) if x.get("modified")
                 else "", "reason": x.get("reason", ""), "path": x["path"]}
                for x in found[:30]]
        columns = (("name", "Experiment"), ("kind", "Kind"), ("runs", "Runs"), ("changed", "Changed")) + (
            (("reason", "Why it belongs to no nucleus"),) if reasons else ())
        table = ui.table(columns=[{"name": c, "label": lab, "field": c, "align": "left"} for c, lab in columns],
                         rows=rows, row_key="path").props("dense flat" + (" wrap-cells" if reasons else "")).classes(
            "w-full")
        table.add_slot("body-cell-name", '<q-td :props="props" style="cursor: pointer"><span class="ps-accent">'
                                         '{{ props.value }}</span></q-td>')
        table.on("rowClick", lambda e: parts.open_experiment(ctx, e.args[1]["path"]))
