"""The start page: the nuclei first, each with its experiments, what is known about it and what is still open
(:func:`nuclei_page`, from :func:`physim.nuclear.nucleus.nuclei`); a nucleus opens its own page (:mod:`.tab_nucleus`).

"Start from a kind of experiment" is the start page as it was (:func:`kinds_page`): the kind of experiment
(:data:`physim.nuclear.workbench.KINDS`), then its beam, target and detectors, or an experiment already made.
Coulomb excitation and elastic scattering then go on to the steps as they always were; fusion-evaporation to its own
screens (:mod:`.tab_fusion`); lifetimes to a Coulomb-excitation or fusion-evaporation setup with a plunger or DSAM
(:mod:`.tab_lifetimes`). With no experiment and no nucleus yet, the kinds come first."""

from __future__ import annotations

import time
from pathlib import Path

from .. import workbench as wb
from . import tab_nucleus as tn


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


def render(ctx) -> None:
    """The start page, into the current container: the nuclei (:func:`nuclei_page`), one nucleus
    (:func:`.tab_nucleus.render`) or the kinds of experiment (:func:`kinds_page`), as ``ctx.state["start"]["view"]``
    says. ``ctx`` gives ``ui``, ``create_experiment``, ``open_experiment``, ``open_example``, ``import_setup``,
    ``refresh_start``, ``close_start`` (None when no experiment is open) and ``state``."""
    st = ctx.state.setdefault("start", {"kind": None})
    if st.get("view") is None:
        st["view"] = "nuclei" if anything_yet() else "kinds"
    if st["view"] == "nucleus" and st.get("nucleus"):
        tn.render(ctx)
    elif st["view"] == "kinds":
        kinds_page(ctx)
    else:
        nuclei_page(ctx)


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
    table they always had. Read on a worker thread (:func:`.tab_nucleus.load_nuclei`), which also makes the rows
    of experiments from before the ledger."""
    ui = ctx.ui
    tn.add_style(ctx)
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
            if ctx.close_start is not None:
                ui.button("Back to the experiment", icon="arrow_back", on_click=lambda: tn.close(ctx)).props(
                    "flat no-caps").classes("ps-muted")
            ui.button("Start from a kind of experiment", on_click=lambda: tn.plan(ctx, None)).props(
                "outline no-caps")
        with ui.row().classes("w-full items-center gap-2"):
            other = ui.input("Another nucleus", placeholder="58Ni").props("dense outlined").style("width: 180px")
            other.on("keydown.enter", lambda: tn.open_nucleus(ctx, other.value))
            ui.button("Open its page", icon="arrow_forward", on_click=lambda: tn.open_nucleus(ctx, other.value)).props(
                "flat no-caps dense").tooltip("The page of any nucleus: ENSDF's adopted values, ready for its first "
                                              "data point or experiment")
        body = ui.column().classes("w-full gap-3")
    with body:
        tn.waiting(ui, "Reading the nuclei and their experiments…")

    async def fill() -> None:
        from nicegui import run

        try:
            found = await run.io_bound(tn.load_nuclei)
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

    tn.later(ui, fill)


def _nuclei(ctx, found: dict) -> None:
    ui = ctx.ui
    display = ctx.state["start"].get("unit")
    cards = [n for n in found["nuclei"] if is_card(n)]
    partners = [n for n in found["nuclei"] if not is_card(n)]
    if not cards and not partners and not found["unassigned"]:
        ui.label("No nucleus yet. Open the page of a nucleus above, or start from a kind of experiment."
                 ).classes("text-sm ps-note p-3 w-full")
    for n in cards:
        name, title = n["nuclide"], tn.nuclide_text(n["nuclide"])
        studied = [x for x in n["experiments"] if x["role"] != "partner"]
        with ui.element("div").classes("ps-panel"):
            with ui.element("div").classes("ps-nuc-head"):
                ui.button(f"{title} →", on_click=lambda name=name: tn.open_nucleus(ctx, name)).props(
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
                          on_click=lambda name=name: tn.plan(ctx, name)).props("flat no-caps dense")
            with ui.element("div").classes("ps-nuc-body"):
                with ui.element("div"):
                    ui.label("WHAT IS KNOWN AND PROJECTED").classes("ps-nuc-h")
                    if not n["table"]:
                        ui.label("The evaluation, on its page." if evaluation["taken"] else "Nothing yet.").classes(
                            "text-sm ps-muted")
                    with ui.element("div").classes("ps-kv"):
                        for t in n["table"][:tn.HEADLINES]:
                            row = ui.label(tn.label_text(t["label"])).style("cursor: pointer")
                            row.on("click", lambda name=name, t=t: tn.open_nucleus(
                                ctx, name, t["observable"], t["key"]))
                            ui.label(tn.known_lines(t, display)).classes("ps-num")
                    if len(n["table"]) > tn.HEADLINES:
                        ui.label(f"… and {len(n['table']) - tn.HEADLINES} more on its page.").classes(
                            "text-xs ps-muted")
                    ui.label("STILL OPEN").classes("ps-nuc-h").style("margin-top: 10px")
                    ui.label(" · ".join(tn.label_text(o["text"]) for o in n["open"]) or "nothing").classes(
                        "text-sm " + ("ps-open" if n["open"] else "ps-muted")).style("line-height: 1.6")
                with ui.element("div"):
                    ui.label("EXPERIMENTS").classes("ps-nuc-h")
                    if n["experiments"]:
                        tn.experiments_table(ctx, n["experiments"], display)
                    else:
                        ui.label("No experiment yet: data points only.").classes("text-sm ps-muted")
    if partners:
        with ui.row().classes("w-full items-center gap-1 ps-note").style("padding: 6px 12px"):
            ui.label("Only as the partner of an experiment, or looked up:" if any(
                not n["experiments"] for n in partners) else "Only as the partner of an experiment:").classes(
                "text-sm ps-muted")
            for n in partners:
                count = len(n["experiments"])
                ui.button(f"{tn.nuclide_text(n['nuclide'])} · " + (
                    f"in {count} experiment{'s' * (count != 1)}" if count else "the evaluation only"),
                          on_click=lambda name=n["nuclide"]: tn.open_nucleus(ctx, name)).props(
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
            if ctx.close_start is not None:
                ui.button("Back to the experiment", icon="arrow_back", on_click=lambda: tn.close(ctx)).props(
                    "flat no-caps").classes("ps-muted")
            ui.button("The nuclei", icon="arrow_back", on_click=lambda: leave(None)).props("outline no-caps").tooltip(
                "The nuclei, each with its experiments and what is known about it")
        if nuclide:
            with ui.row().classes("w-full items-center gap-2 ps-note").style("padding: 6px 12px"):
                ui.label(f"An experiment for {tn.nuclide_text(nuclide)}: pick its kind. Where the form has a field "
                         f"for the nucleus that is excited, {nuclide} is filled in.").classes("text-sm")
                ui.space()
                ui.button(f"Back to {tn.nuclide_text(nuclide)}", on_click=lambda: leave(nuclide)).props(
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
        tn.open_nucleus(ctx, to) if to else tn.open_nuclei(ctx)


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
        table.on("rowClick", lambda e: tn.open_experiment(ctx, e.args[1]["path"]))
