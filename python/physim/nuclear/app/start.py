"""The start page: the kind of experiment first (:data:`physim.nuclear.workbench.KINDS`), then its beam, target and
detectors, or an experiment already made. Coulomb excitation and elastic scattering then go on to the steps as
they always were; fusion-evaporation to its own screens (:mod:`.tab_fusion`); lifetimes to a Coulomb-excitation or
fusion-evaporation setup with a plunger or DSAM (:mod:`.tab_lifetimes`)."""

from __future__ import annotations

import time
from pathlib import Path

from .. import workbench as wb


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


def render(ctx) -> None:
    """The start page, into the current container. ``ctx`` gives ``ui``, ``create_experiment``,
    ``open_experiment``, ``open_example``, ``import_setup``, ``close_start`` (None when no experiment is open) and
    ``state``."""
    ui = ctx.ui
    st = ctx.state.setdefault("start", {"kind": None})
    with ui.column().classes("w-full gap-3").style("max-width: 1060px; margin: 0 auto"):
        with ui.row().classes("w-full items-end"):
            with ui.column().classes("gap-1"):
                ui.label("What experiment are you planning?").classes("text-xl font-semibold")
                ui.label("Pick the kind of experiment first: it decides which detectors the setup offers, what the "
                         "plan computes and which spectra you see. The beam, target and detectors can change "
                         "afterwards.").classes("text-sm ps-muted").style("max-width: 70ch")
            ui.space()
            if ctx.close_start is not None:
                ui.button("Back to the experiment", icon="arrow_back", on_click=ctx.close_start).props(
                    "flat no-caps").classes("ps-muted")
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


def kind_form(ctx, k) -> None:
    ui = ctx.ui
    st = ctx.state["start"]
    values = st.setdefault("values", dict(zip(("beam", "energy", "target", "thickness"), k.defaults),
                                          template=k.templates[0], name="",
                                          excite="beam" if k.key == "lifetimes" else "target", fill=True))
    with ui.column().classes("w-full ps-card gap-0"):
        with ui.row().classes("w-full items-center ps-panel-head"):
            ui.label(k.label)
            ui.label("the beam, the target and where the detectors come from").classes("ps-muted text-sm")
            ui.space()
            ui.button("Change the kind", on_click=lambda: (st.update(kind=None), ctx.refresh_start())).props(
                "flat no-caps dense")
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


def experiments_table(ctx) -> None:
    ui = ctx.ui
    from .. import runs as _runs

    found = _runs.list_experiments()
    with ui.column().classes("w-full ps-card gap-0"):
        with ui.row().classes("w-full items-center ps-panel-head"):
            ui.label("Or open an experiment")
            ui.label(f"in {_runs.home()}").classes("ps-muted text-sm")
            ui.space()
            ui.button("Import a setup file…", on_click=ctx.import_setup).props("flat no-caps dense")
        if not found:
            ui.label("No experiment yet.").classes("text-sm p-3 ps-muted")
            return
        rows = [{"name": x["name"], "kind": experiment_kind(x["path"]), "runs": x["runs"],
                 "changed": time.strftime("%Y-%m-%d %H:%M", time.localtime(x["modified"])), "path": x["path"]}
                for x in found[:30]]
        table = ui.table(columns=[{"name": c, "label": lab, "field": c, "align": "left"} for c, lab in
                                  (("name", "Experiment"), ("kind", "Kind"), ("runs", "Runs"), ("changed", "Changed"))],
                         rows=rows, row_key="path").props("dense flat").classes("w-full")
        table.add_slot("body-cell-name", '<q-td :props="props" style="cursor: pointer"><span class="ps-accent">'
                                         '{{ props.value }}</span></q-td>')
        table.on("rowClick", lambda e: ctx.open_experiment(e.args[1]["path"]))
