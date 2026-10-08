"""The planner page: the setup panel, the status strip and the stage tabs, assembled for one browser client."""

from __future__ import annotations

import base64
import math
import time
from typing import Optional


from .. import guide
from ..analysis import SHIFT_H
from ..record import record_css
from ..planner import Planner
from ..quantity import Quantity

from .. import dataviews as dv
from .figures import (CLOVER_FIELDS, REACTION_FIELDS, REACTION_TYPES, STYLE, THEMES, _shown, _value,
                      figure_alignment, figure_gamma_crystal, figure_detector_spectrum, figure_efficiency, figure_energy_loss,
                      figure_excitation, figure_gamma_spectra, figure_kinematics, figure_levels, figure_overlay,
                      figure_source_spectrum, figure_spectra, figure_strips, figure_sweep, figure_trajectories,
                      report_zip, themed)

# ---------------------------------------------------------------------------------------------------------------
# The page


def _fmt(x, digits=4) -> str:
    if x is None:
        return "—"
    if isinstance(x, float) and not math.isfinite(x):
        return "never" if x > 0 else "—"
    return f"{x:.{digits}g}"


def _time(seconds) -> str:
    if seconds is None:
        return "—"
    if not math.isfinite(seconds):
        return "never"
    if seconds >= 3600:
        return f"{seconds / 3600:.3g} h"
    return f"{seconds / 60:.3g} min" if seconds >= 60 else f"{seconds:.3g} s"


def _prepared(example: str) -> Planner:
    """An example's planner with its rates already computed (run in a worker thread while the page waits)."""
    p = Planner.example(example)
    p.rates()
    return p


def _starting_planner(example: Optional[str], template: Optional[str], experiment: Optional[str]) -> tuple:
    """The planner a page opens with, and whether to offer the experiment list first: the experiment asked for,
    else a template or example asked for, else the most recent experiment, else the first template."""
    from .. import runs as _runs
    from .. import workbench as wb

    if experiment:
        return Planner.open_experiment(experiment), False
    if template:
        return Planner(wb.template(template)), False
    if example:
        return Planner.example(example), False
    found = _runs.list_experiments()
    if found:
        try:
            return Planner.open_experiment(found[0]["path"]), False
        except Exception:  # noqa: BLE001 - a broken folder: start from a template instead
            pass
    return Planner(wb.template(wb.TEMPLATES[0][0])), True


def build_page(example: Optional[str] = None, events: int = 100_000, theme: Optional[str] = None,
               template: Optional[str] = None, experiment: Optional[str] = None) -> None:
    """Build the planner page for the current client (call inside a NiceGUI page function).

    It opens ``experiment`` (a folder), or a new setup from ``template`` (:data:`physim.nuclear.workbench.TEMPLATES`)
    or ``example``, or else the most recent experiment; with no experiment yet it offers the experiment list and
    "New experiment" first. ``theme`` is "light" or "dark"; ``None`` takes the one this browser used last, the
    system's setting the first time."""
    from types import SimpleNamespace

    from nicegui import run, ui

    from .. import paper
    from .. import runs as _runs
    from .. import workbench as wb
    from ..experiment import Experiment, SetupError
    from . import panel as _panel
    from . import strip as _strip
    from . import tab_analysis, tab_data, tab_physics, tab_plan, tab_report, tab_run, tab_setup

    first, offer_list = _starting_planner(example, template, experiment)
    state = {"planner": first, "events": events, "theme": theme or "light", "open": None, "stage": "setup",
             "export": None, "journal": "physical_review", "width": "single", "format": "pdf"}

    def P() -> Planner:  # noqa: N802
        return state["planner"]

    # -- actions ----------------------------------------------------------------------------------------------
    def changed(ok: bool) -> None:
        if ok:
            refresh_results()
        else:
            ui.notify("; ".join(P().problems), type="negative", multi_line=True)
        setup_panel.refresh()
        status_strip.refresh()

    def card_of(section: str) -> Optional[str]:
        if section in ("beam", "run"):
            return section
        if section in ("target", "backing", "reaction"):
            return "target"
        if section.startswith("gamma detector"):
            return f"gamma:{int(section.split()[2]) - 1}"
        if section.startswith("detector"):
            return f"detector:{int(section.split()[1]) - 1}"
        return None

    def open_card(key: Optional[str]) -> None:
        state["open"] = key if key != "checks" else None
        drawer.show()
        setup_panel.refresh()

    def edit(section: str, field: str, text) -> None:
        value = _value(field, text)
        d = P().draft
        if section == "title":
            old = d.get("title")
        elif section == "backing":
            old = d["target"].get("backing", {}).get(field)
        elif section.startswith("gamma detector"):
            old = d["gamma_detectors"][int(section.split()[2]) - 1].get(field)
        elif section.startswith("detector"):
            old = d["detectors"][int(section.split()[1]) - 1].get(field)
        else:
            old = d.get(section, {}).get(field)
        if (None if old is None else str(old)) == (None if value is None else str(value)):
            return
        if P().running:
            ui.notify("A run is being taken: the setup is read-only until it ends.", type="warning")
            return
        state["open"] = card_of(section) or state["open"]
        changed(P().set(section, field, value))

    def add_detector(key: str = "side"):
        def run_() -> None:
            names = [det.get("name") for det in P().draft["detectors"]]
            changed(P().add_detector(**guide.placement(key, names)))
            setup_panel.refresh()
        return run_

    def duplicate_detector(k: int):
        def run() -> None:
            changed(P().duplicate_detector(k))
            setup_panel.refresh()
        return run

    def remove_detector(k: int):
        def run() -> None:
            if len(P().draft["detectors"]) == 1:
                ui.notify("A setup needs at least one particle detector.", type="warning")
                return
            changed(P().remove_detector(k))
            setup_panel.refresh()
        return run

    def add_gamma_detector(key: str = "disc"):
        def run() -> None:
            names = [g.get("name") for g in P().draft.get("gamma_detectors", [])]
            changed(P().add_gamma_detector(**guide.gamma_placement(key, names)))
            setup_panel.refresh()
        return run

    def duplicate_gamma_detector(k: int):
        def run() -> None:
            changed(P().duplicate_gamma_detector(k))
            setup_panel.refresh()
        return run

    def remove_gamma_detector(k: int):
        def run() -> None:
            changed(P().remove_gamma_detector(k))
            setup_panel.refresh()
        return run

    def use(planner: Planner) -> None:
        state["planner"] = planner
        state["open"] = None
        state["compound"] = False
        setup_panel.refresh()
        status_strip.refresh()
        refresh_results()
        main_area.refresh()

    async def open_experiment(path: str) -> None:
        experiments_dialog.close()
        note = ui.notification("Opening the experiment…", spinner=True, timeout=None)
        try:
            planner = await run.io_bound(Planner.open_experiment, path)
        except Exception as err:  # noqa: BLE001 - say why instead of failing
            ui.notify(f"Could not open {path}: {err}", type="negative", multi_line=True)
            return
        finally:
            note.dismiss()
        use(planner)

    async def create_experiment(beam: str, energy: str, target: str, thickness: str, measure: str,
                                template_key: str, name: str) -> None:
        def make() -> tuple:
            exp = wb.new_setup(beam, energy, target, thickness or "0.5 mg/cm2", measure,
                               detectors_from=template_key)
            notes = wb.new_setup_notes(exp)
            planner = Planner(exp)
            planner.create_experiment(name or None)
            planner.rates()
            return planner, notes

        note = ui.notification("Making the experiment…", spinner=True, timeout=None)
        try:
            planner, notes = await run.io_bound(make)
        except (SetupError, ValueError, KeyError) as err:
            ui.notify(f"Could not make the experiment: {err}", type="negative", multi_line=True)
            return
        finally:
            note.dismiss()
        new_dialog.close()
        experiments_dialog.close()
        use(planner)
        for n in notes:
            ui.notify(n, type="warning", multi_line=True, close_button=True)

    async def import_file(e) -> None:
        upload_dialog.close()
        if P().running:
            ui.notify("A run is being taken: the setup is read-only until it ends.", type="warning")
            return
        try:
            exp = Experiment.from_toml(await e.file.text())
        except (SetupError, ValueError) as err:
            ui.notify(f"Could not read the setup: {err}", type="negative", multi_line=True)
            return
        folder, current = P().folder, P().run
        planner = Planner(exp)
        if folder is not None:  # the imported setup becomes the experiment's current setup
            planner.folder, planner.run = folder, current
            folder.save_setup(exp)
        use(planner)

    def export_setup() -> None:
        if P().problems:
            ui.notify("Fix the setup problems first.", type="warning")
            return
        ui.download.content(P().to_toml(), "setup.toml")

    # -- setup panel ------------------------------------------------------------------------------------------
    def bind(el, sec_name: str, field: str) -> None:
        """Apply an input's text when it loses focus or Enter is pressed (not on every keystroke)."""

        def apply() -> None:
            edit(sec_name, field, el.value)

        el.on("blur", apply)
        el.on("keydown.enter", apply)

    def help_icon(element, sec_name: str, field: str) -> None:
        """A help icon inside an input: hover for the help, or click it (on a touch screen)."""
        h = guide.help_for(sec_name, field)
        if h is None:
            return
        with element.add_slot("append"):
            icon = ui.icon("help_outline", size="xs").classes("cursor-help ps-muted")
            with icon, ui.tooltip().classes("bg-slate-800 text-white max-w-xs p-2"):
                ui.label(h.what).classes("text-sm")
                ui.label(f"Typical: {h.typical}").classes("text-xs text-slate-300 mt-1")
                if h.effect != "—":
                    ui.label(f"Raise it: {h.effect}").classes("text-xs text-slate-300")
            icon.on("click.stop", lambda: ui.notify(h.text(), multi_line=True, close_button=True,
                                                    classes="whitespace-pre-line"))

    def section(title, fields, sec_name, values):
        if title:
            ui.label(title).classes("ps-section mt-3")
        # The guided steps have few fields each: one per line leaves room for the labels.
        with ui.grid(columns=2).classes("w-full gap-1"):
            for field, label, placeholder in fields:
                v = values.get(field)
                if (sec_name, field) in (("beam", "nuclide"), ("target", "material")):
                    nuclide_picker(sec_name, field, v, placeholder)
                    continue
                inp = ui.input(label, value=_shown(field, v), placeholder=placeholder).props(
                    "dense outlined").classes("w-full")
                bind(inp, sec_name, field)
                help_icon(inp, sec_name, field)

    def nuclide_picker(sec_name: str, field: str, value, placeholder: str) -> None:
        """The beam's nuclide or the target's material as two choices, the element and its mass number (a target may
        also be the element's natural mix, or a compound typed in), in place of one text field."""
        target = sec_name == "target"
        parts = wb.split_nuclide(value)
        options = dict(wb.elements())
        if target:
            options = {wb.COMPOUND: "Compound or other material", **options}
        compound = target and (state.get("compound") or (value and parts is None))
        symbol = wb.COMPOUND if compound else (parts[0] if parts else None)

        def set_element(new: str) -> None:
            if new == wb.COMPOUND:
                state["compound"] = True
                setup_panel.refresh()
                return
            state["compound"] = False
            a = wb.most_abundant(new)
            edit(sec_name, field, (new if target else f"{a}{new}") if a else f"{next(iter(wb.mass_numbers(new)))}{new}")

        el = ui.select(options, value=symbol if symbol in options else None, with_input=True,
                       label="Target element" if target else "Beam element",
                       on_change=lambda e: set_element(e.value)).props("dense outlined options-dense").classes("w-full")
        help_icon(el, sec_name, field)
        if compound:
            inp = ui.input("Material", value=_shown(field, value), placeholder=placeholder).props(
                "dense outlined").classes("w-full")
            bind(inp, sec_name, field)
            return
        if symbol is None:
            ui.label("")
            return
        masses = {str(a): label for a, label in wb.mass_numbers(symbol).items()}
        if target:
            masses = {wb.NATURAL: "natural mix", **masses}
        current = str(parts[1]) if parts and parts[1] is not None else (wb.NATURAL if target else None)
        ui.select(masses, value=current if current in masses else None, label="Mass number",
                  on_change=lambda e: edit(sec_name, field, symbol if e.value == wb.NATURAL else f"{e.value}{symbol}")
                  ).props("dense outlined options-dense").classes("w-full")

    @ui.refreshable
    def setup_panel():
        _panel.render(ctx)

    def clover_switches(i: int, gd: dict, numbers: bool = True) -> None:
        """Add-back and Compton suppression for a clover: two switches, and (``numbers``) the numbers behind
        them."""
        sec = f"gamma detector {i + 1}"
        with ui.row().classes("items-center gap-4"):
            ui.switch("Add-back", value=bool(gd.get("addback")),
                      on_change=lambda e: edit(sec, "addback", "yes" if e.value else None)).props("dense").tooltip(
                guide.help_for("gamma", "addback").text())
            ui.switch("Compton suppression (BGO shield)", value=gd.get("shield") is not None,
                      on_change=lambda e: edit(sec, "shield", "BGO" if e.value else None)).props("dense").tooltip(
                guide.help_for("gamma", "shield").text())
        if numbers:
            section("", CLOVER_FIELDS, sec, gd)

    def reaction_section(reaction: dict, title: str = "Reaction") -> None:
        if title:
            ui.label(title).classes("ps-section mt-3")
        kind = reaction.get("type", "elastic")
        sel = ui.select(REACTION_TYPES, value=kind, label="What happens in the target",
                        on_change=lambda e: edit("reaction", "type", e.value)).props("dense outlined").classes(
            "w-full")
        help_icon(sel, "reaction", "type")
        if kind != "coulex":
            return
        with ui.grid(columns=2).classes("w-full gap-1"):
            exc = ui.select({"target": "Target nucleus", "projectile": "Beam nucleus"},
                            value=reaction.get("excite", "target"), label="Excited nucleus",
                            on_change=lambda e: edit("reaction", "excite", e.value)).props("dense outlined")
            help_icon(exc, "reaction", "excite")
            mul = ui.select(["E1", "E2", "E3"], value=reaction.get("multipolarity", "E2"), label="Multipolarity",
                            on_change=lambda e: edit("reaction", "multipolarity", e.value)).props("dense outlined")
            help_icon(mul, "reaction", "multipolarity")
        state_picker(reaction)
        section("", REACTION_FIELDS, "reaction", reaction)

    def state_picker(reaction: dict) -> None:
        """The excited state, chosen from the excited nucleus's level scheme (read from ENSDF); choosing one fills in
        its energy and B(Eλ↑) below, which can still be typed over."""
        role = "target" if reaction.get("excite", "target") == "target" else "beam"
        multipolarity = reaction.get("multipolarity", "E2")
        scheme = P().experiment.levels.get(role) if not P().problems else None
        choices = wb.states(scheme, multipolarity)
        nuclide = P().level_nuclides()[role] or role
        if not choices:
            ui.label(f"No state of {nuclide} with a known B({multipolarity}↑) in the local ENSDF copy: enter the "
                     "state's energy and B(Eλ↑) below.").classes("text-xs ps-muted")
            return
        try:
            kev = Quantity.parse(str(reaction.get("energy"))).to("keV")
        except Exception:  # noqa: BLE001 - no energy yet, or one being typed
            kev = None
        current = next((n for n in choices if kev is not None and abs(scheme.levels[n].energy.value - kev) < 0.5),
                       None)

        def pick(level) -> None:
            if level is None or level == current:
                return
            try:
                changed(P().use_state(role, level, multipolarity))
            except ValueError as err:
                ui.notify(str(err), type="warning")

        ui.select(choices, value=current, label=f"State of {nuclide}",
                  on_change=lambda e: pick(e.value)).props("dense outlined options-dense").classes("w-full")

    # -- figures and their export ---------------------------------------------------------------------------------
    def plot(fig, name: str, **options) -> None:
        """A figure on its plate, with the button that exports it in a journal's style."""
        with ui.column().classes("w-full ps-plate gap-0"):
            with ui.row().classes("w-full items-center"):
                ui.space()
                ui.button("Paper figure", icon="article", on_click=lambda: open_export(name, options)).props(
                    "dense flat no-caps").tooltip("Export this figure in a journal's style")
            ui.plotly(themed(fig, state["theme"])).classes("w-full")

    def open_export(name: str, options: dict) -> None:
        state["export"] = (name, options)
        export_controls.refresh()
        export_dialog.open()

    async def render_preview() -> None:
        name, options = state["export"]
        j, w = state["journal"], state["width"]
        preview_box.clear()
        with preview_box:
            ui.spinner(size="lg")
        try:
            png = await run.io_bound(lambda: paper.preview(P(), name, j, w, **options))
        except Exception as err:  # noqa: BLE001 -- show the reason instead of an empty dialog
            preview_box.clear()
            with preview_box:
                ui.label(f"Could not draw the figure: {err}").classes("ps-bad")
            return
        if (j, w) != (state["journal"], state["width"]):  # the choice changed while this one was drawn
            return
        preview_box.clear()
        mm = paper.width_mm(j, w)
        style = paper.journal(j)
        with preview_box:
            with ui.element("div").classes("ps-paper"):
                ui.html(f'<img alt="{paper.FIGURES[name]} in the {style.label} style" style="width: {mm}mm" '
                        f'src="data:image/png;base64,{base64.b64encode(png).decode()}">', sanitize=False)
            ui.label(f"{mm:g} mm wide · {style.size:g} pt {style.family} lettering").classes("text-xs ps-muted")

    async def choose(key: str, value) -> None:
        if state[key] == value:
            return
        state[key] = value
        if key == "journal" and state["width"] not in paper.JOURNALS[value].widths:
            state["width"] = "single"
        if key != "format":
            export_controls.refresh()
            await render_preview()

    async def download_figure() -> None:
        name, options = state["export"]
        j, w, fmt = state["journal"], state["width"], state["format"]
        note = ui.notification("Drawing the figure…", spinner=True, timeout=None)
        try:
            data = await run.io_bound(lambda: paper.export(P(), name, fmt, j, w, **options))
        finally:
            note.dismiss()
        ui.download.content(data, f"{name}-{j}-{w}.{fmt}")

    def download_data() -> None:
        name, options = state["export"]
        ui.download.content(paper.data_csv(P(), name, **options), f"{name}.csv")

    @ui.refreshable
    def export_controls() -> None:
        if state["export"] is None:
            return
        style = paper.JOURNALS[state["journal"]]
        ui.label(f"Paper figure: {paper.FIGURES[state['export'][0]]}").classes("text-lg font-semibold")
        ui.select({k: j.label for k, j in paper.JOURNALS.items()}, value=state["journal"], label="Journal style",
                  on_change=lambda e: choose("journal", e.value)).props("dense outlined").classes("w-full")
        ui.label(style.note + " Check the journal's current guidelines before you submit.").classes(
            "text-xs ps-muted")
        ui.label("Width").classes("ps-section")
        ui.toggle({k: f"{paper.WIDTHS[k]} · {mm:g} mm" for k, mm in style.widths.items()}, value=state["width"],
                  on_change=lambda e: choose("width", e.value)).props("dense no-caps unelevated")
        ui.label("Format").classes("ps-section")
        ui.toggle({k: label for k, (label, _) in paper.FORMATS.items()}, value=state["format"],
                  on_change=lambda e: choose("format", e.value)).props("dense no-caps unelevated")
        ui.button("Download figure", icon="download", on_click=download_figure).props("no-caps unelevated")
        ui.button("Download the plotted data (CSV)", icon="table_view", on_click=download_data).props(
            "no-caps flat")
        ui.label("White background in either theme. Text stays editable text in PDF and SVG, and every curve has "
                 "its own line style, so the figure also reads in greyscale.").classes("text-xs ps-muted")

    # -- key results ------------------------------------------------------------------------------------------
    @ui.refreshable
    def readouts() -> None:
        try:
            r, t, ws = P().rates(), P().trajectories(), P().warnings()
        except Exception:  # noqa: BLE001 -- the strip must never break the page
            return
        rows = r["rows"]
        what = {"all": "counts", "excitations": "excitations", "coincidences": "coincidences"}[r["measured"]]
        fast = max(rows, key=lambda x: x["rate_per_s"])
        rate = fast["rate_per_s"]
        cells = [("Highest rate", (f"{rate:,.0f}" if rate >= 1000 else _fmt(rate, 3)) + " /s",
                  f"particles in {fast['detector']}", "")]
        timed = [x for x in rows if x["beam_time_s"] is not None]
        if r["counts_wanted"] and timed:
            slow = max(timed, key=lambda x: x["beam_time_s"])
            cells.append(("Longest beam time", _time(slow["beam_time_s"]),
                          f"{slow['detector']}, for {r['counts_wanted']} {what}", ""))
        else:
            few = min(rows, key=lambda x: x["counts_in_run"])
            cells.append(("Fewest counts in the run", _fmt(few["counts_in_run"], 3), f"{what} in {few['detector']}",
                          ""))
        cells.append(("Closest approach, head-on", f"{t['d0_fm']:.1f} fm",
                      f"nuclear range {t['interaction_radius_fm']:.1f} fm", ""))
        n = {level: sum(w.level == level for w in ws) for level in ("error", "warning", "note")}
        if n["error"]:
            checks = (f"{n['error']} error" + "s" * (n["error"] > 1), "ps-bad")
        elif n["warning"]:
            checks = (f"{n['warning']} warning" + "s" * (n["warning"] > 1), "ps-warn")
        else:
            checks = ("No warnings", "ps-ok")
        cells.append(("Setup checks", checks[0], f"{n['note']} note" + "s" * (n["note"] != 1), checks[1]))
        with ui.element("div").classes("ps-readouts"):
            for title, value, sub, colour in cells:
                with ui.element("div").classes("ps-readout"):
                    ui.label(title).classes("text-xs ps-muted")
                    ui.label(value).classes("ps-readout-value " + colour)
                    ui.label(sub).classes("text-xs ps-muted")

    # -- result tabs ------------------------------------------------------------------------------------------
    def columns(spec):
        return [{"name": k, "label": lab, "field": k, "align": "left"} for k, lab in spec]

    # -- the scene and its side panel ---------------------------------------------------------------------------
    scene_state = {"view": None, "key": None, "element": None, "live": {}, "track": None, "full": False, "extras": [],
                   "tracks": {"n": 30, "select": "all", "weighted": True, "playing": True, "speed": 1.0},
                   "source": {"nuclide": "152Eu", "activity": "37 kBq", "time": "1 h", "run": False}}

    def scene_selected(key, element) -> None:
        scene_state["key"], scene_state["element"] = key, element
        scene_state["track"] = None
        selection_panel.refresh()

    def scene_track(track) -> None:
        scene_state["track"] = track
        selection_panel.refresh()

    def scene_live(values: dict) -> None:
        """Follow a detector while it is dragged: the quick numbers go straight into the side panel."""
        labels = scene_state["live"]
        text = {}
        if "theta" in values:
            text["place"] = (f"θ {values['theta']:.1f}°, φ {values['phi']:.1f}°, "
                             f"{values['distance_mm']:.1f} mm")
        if "theta_range" in values:
            text["covers"] = f"{values['theta_range'][0]:.1f}–{values['theta_range'][1]:.1f}°"
            text["omega"] = f"{_fmt(values['solid_angle_msr'])} msr"
            text["hidden"] = f"{100 * values['hidden']:.0f} %"
            text["rate"] = f"{_fmt(values['rate_per_s'])} /s"
        if "half_angle_deg" in values:
            text["half"] = f"{values['half_angle_deg']:.1f}°"
            text["geometric"] = f"{100 * values['geometric_efficiency']:.2f} %"
        text["status"] = f"Not allowed here: {values['problem']}" if values.get("problem") else "Moving…"
        for name, value in text.items():
            if name in labels:
                labels[name].set_text(value)

    def scene_moved(key: str, position, mode: str) -> None:
        """A detector was dropped at a new place: the scene and its side panel follow at once, the other results
        (and the Monte Carlo) a moment later."""
        ok = P().move(key, position, mode)
        selection_panel.refresh()

        def rest() -> None:
            changed(ok)
            setup_panel.refresh()

        ui.timer(0.05, rest, once=True)

    def scene_current() -> bool:
        """Redraw the scene in place (the camera stays), if it is on the page and in the right theme."""
        view = scene_state["view"]
        if view is None or not view.alive or view.theme != state["theme"]:
            return False
        if not view.busy and view.drawn is not P().experiment:
            view.draw()
        selection_panel.refresh()
        geometry_tables.refresh()
        return True

    def typed_place(key: str, field: str, unit: str):
        def run_(e) -> None:
            text = str(e.sender.value or "").strip()
            if not text:
                return
            try:
                float(text)
                text = f"{text} {unit}"
            except ValueError:
                pass
            kind, i = key.split(":")
            entry = P().draft["detectors" if kind == "detector" else "gamma_detectors"][int(i)]
            if str(entry.get(field)) == text:
                return
            changed(P().place(key, **{field: text}))
            setup_panel.refresh()
        return run_

    def open_explanation(key: str, detector: Optional[str] = None, gamma_detector: Optional[str] = None,
                         result=None) -> None:
        """A dialog with one number's explanation: formula, the numbers substituted, meaning, assumptions."""
        try:
            found = [x for x in P().explanations(detector, gamma_detector, result) if x["key"] == key]
        except Exception as err:  # noqa: BLE001 -- an explanation must never break the page
            ui.notify(f"No explanation: {err}", type="warning")
            return
        if not found:
            ui.notify("No explanation for this number yet.", type="warning")
            return
        x = found[0]
        with ui.dialog() as dialog, ui.card().classes("ps-card").style("max-width: min(760px, 95vw)"):
            ui.label(x["title"]).classes("text-base font-medium")
            ui.label(f"{_fmt(x['value'])} {x['unit']}".strip()).classes("ps-readout-value")
            ui.label("Formula").classes("ps-section mt-2")
            ui.label(x["formula"]).classes("text-sm ps-num")
            ui.label("With this run's numbers").classes("ps-section mt-2")
            ui.label(x["substituted"]).classes("text-sm ps-num")
            ui.label("What it means").classes("ps-section mt-2")
            ui.label(x["meaning"]).classes("text-sm")
            ui.label("Assumptions").classes("ps-section mt-2")
            ui.label(x["assumptions"]).classes("text-sm ps-muted")
            ui.link("The theory page", x["reference"], new_tab=True).classes("text-sm")
            with ui.row().classes("gap-2"):
                section = wb.physics_for(x["key"])
                ui.button("See the physics", icon="science",
                          on_click=lambda s=section: (dialog.close(), show_physics(s))).props("flat no-caps").tooltip(
                    "The panel of the Physics tab that shows it: " + dict((k, t) for k, t, _, _ in wb.PHYSICS)[section])
                ui.button("Close", on_click=dialog.close).props("flat no-caps")
        dialog.open()

    def track_panel(t) -> None:
        """One event's numbers, for a selected track."""
        with ui.column().classes("w-full ps-card p-3 gap-1"):
            ui.label("One event").classes("ps-section")
            ui.label(f"Event {t.event}: {t.channel}").classes("text-base font-medium")
            ui.label(f"Stands for {t.weight:.3g} events per second; "
                     + ("a particle–γ coincidence." if t.coincidence else "no γ ray recorded.")).classes(
                "text-sm ps-muted")
            for pr in t.particles:
                ui.label(f"{'Recoil' if pr['recoil'] else 'Scattered beam'} in {P().experiment.detectors[pr['detector']].name} "
                         f"segment ({pr['segment_i'] + 1}, {pr['segment_j'] + 1}): θ {pr['theta']:.1f}°, "
                         f"φ {pr['phi']:.1f}°, {pr['energy']:.2f} MeV at the reaction, {pr['measured']:.2f} MeV "
                         f"measured" + ("" if pr["counted"] else " (below threshold)")).classes("text-sm")
            ui.label(f"Reaction at {t.particles[0]['depth']:.3f} mg/cm² into the target, beam at "
                     f"{t.particles[0]['beam_energy']:.2f} MeV, θ_CM {t.particles[0]['theta_cm']:.1f}°").classes(
                "text-sm")
            for g in t.gammas:
                ui.label(f"γ ray of {1e3 * g['energy0']:.1f} keV in {t.crystal_hits[t.gammas.index(g)]}: "
                         f"{1e3 * g['energy_lab']:.1f} keV in the laboratory, {1e3 * g['measured']:.1f} keV measured, "
                         f"corrected {1e3 * g['corrected_recoil']:.1f} (recoil) / {1e3 * g['corrected_projectile']:.1f} "
                         f"(projectile) keV; β = {g['beta']:.4f}").classes("text-sm")
            ui.button("Back to the detectors", on_click=lambda: scene_state["view"].select_track(None)).props(
                "dense flat no-caps")

    @ui.refreshable
    def selection_panel() -> None:
        key, element = scene_state["key"], scene_state["element"]
        if scene_state["track"] is not None:
            track_panel(scene_state["track"])
            return
        try:
            s = P().selection(key, element)
        except (IndexError, KeyError, ValueError):
            key, s = None, P().selection()
        live = scene_state["live"] = {}

        def row(label: str, value: str, name: Optional[str] = None, classes: str = "",
                explain: Optional[str] = None) -> None:
            with ui.row().classes("w-full justify-between no-wrap gap-3 items-center"):
                with ui.row().classes("items-center no-wrap gap-1"):
                    ui.label(label).classes("text-sm ps-muted")
                    if explain:
                        ui.button(icon="help_outline", on_click=lambda e_=explain: open_explanation(
                            e_, s.get("name") if s["kind"] == "detector" else None,
                            s.get("name") if s["kind"] == "gamma" else None)).props(
                            "dense flat round size=xs").tooltip("How this number is obtained")
                lab = ui.label(value).classes("text-sm ps-num text-right " + classes)
            if name:
                live[name] = lab

        def place_inputs(fields) -> None:
            kind, i = key.split(":")
            entry = P().draft["detectors" if kind == "detector" else "gamma_detectors"][int(i)]
            if "position" in entry:
                ui.label("Placed by position: edit it in the setup panel.").classes("text-xs ps-muted")
                return
            ui.label("Exact values").classes("ps-section mt-2")
            with ui.row().classes("w-full no-wrap gap-2"):
                for field, label, unit in fields:
                    box = ui.input(label, value=str(entry.get(field, "") or "")).props("dense outlined").classes(
                        "min-w-0")
                    box.on("blur", typed_place(key, field, unit))
                    box.on("keydown.enter", typed_place(key, field, unit))

        with ui.column().classes("w-full ps-card p-3 gap-1"):
            if s["kind"] == "experiment":
                ui.label("The whole experiment").classes("ps-section")
                ui.label(s["title"]).classes("text-base font-medium")
                row("Particle detectors", str(s["detectors"]))
                row("γ-ray detectors", str(s["gamma_detectors"]))
                row("Solid angle, all particle detectors", f"{_fmt(s['solid_angle_msr'])} msr")
                row("Particles counted", f"{_fmt(s['rate_per_s'])} /s")
                if s["measured"] != "all":
                    row(f"Of these, {s['measured']}", f"{_fmt(s['measured_rate_per_s'])} /s")
                if s["gamma_efficiency"] is not None:
                    row("γ-ray efficiency, all detectors", f"{100 * s['gamma_efficiency']:.2f} %")
                t = s["target"]
                row("Target", f"{t['material']}, {_fmt(t['thickness_um'], 3)} µm"
                    + (f", tilted {t['tilt_deg']:g}°" if t["tilt_deg"] else ""))
                for line in s["advice"]:
                    with ui.row().classes("items-start no-wrap gap-2 mt-1"):
                        ui.icon("lightbulb").classes("ps-accent mt-0.5")
                        ui.label(line).classes("text-sm")
                ui.label("Click a detector in the scene for its numbers.").classes("text-xs ps-muted mt-1")
                return
            ui.label("Particle detector" if s["kind"] == "detector" else "γ-ray detector").classes("ps-section")
            ui.label(s["name"] + (f" · {s['model']}" if s.get("model") else "")).classes("text-base font-medium")
            row("Centre", f"θ {s['theta']:.1f}°, φ {s['phi']:.1f}°, {s['distance_mm']:.1f} mm", "place")
            if s["kind"] == "detector":
                row("Covers θ", f"{s['theta_range'][0]:.1f}–{s['theta_range'][1]:.1f}°", "covers")
                row("Solid angle", f"{_fmt(s['solid_angle_msr'])} msr", "omega", explain="solid_angle")
                row("Hidden by others", f"{100 * s['hidden']:.0f} %", "hidden")
                row("Particles reaching it", f"{_fmt(s['rate_per_s'])} /s", "rate", explain="rate")
                if s["measured"] != "all":
                    row(f"Of these, {s['measured']}", f"{_fmt(s['measured_rate_per_s'])} /s")
                row("Counts in the run", _fmt(s["counts_in_run"]), explain="counts")
                if s["rate_per_s"] == 0:
                    ui.label("Nothing reaches this detector: the kinematics send no particle here, or another "
                             "detector hides it.").classes("text-sm ps-warn")
                if s["unsafe"]:
                    word = "strips" if s["shape"] == "rectangle" else "rings"
                    ui.label(f"{len(s['unsafe'])} of its {word} (violet in the scene) see collisions closer than "
                             "Cline's safe distance.").classes("text-sm ps-warn")
            else:
                row("Half-angle", f"{s['half_angle_deg']:.1f}°", "half")
                row("Geometric coverage", f"{100 * s['geometric_efficiency']:.2f} %", "geometric")
                row(f"Full-energy-peak efficiency at {s['gamma_energy_kev']:.0f} keV",
                    f"{100 * s['peak_efficiency']:.3f} %", explain="gamma_efficiency")
                if s["coincidence_rate_per_s"] is not None:
                    row("Particle–γ coincidences", f"{_fmt(s['coincidence_rate_per_s'])} /s", explain="coincidences")
                if len(s["crystals"]) > 1:
                    row("Crystals", ", ".join(s["crystals"]))
            status = ui.label("").classes("text-xs ps-warn")
            live["status"] = status
            el = s.get("element")
            if el:
                ui.separator().classes("my-1")
                ui.label(el["label"]).classes("ps-section")
                if "theta_range" in el:
                    row("Covers θ", f"{el['theta_range'][0]:.1f}–{el['theta_range'][1]:.1f}°")
                    row("Solid angle", f"{_fmt(el['solid_angle_msr'])} msr")
                    row("Particles reaching it", f"{_fmt(el['rate_per_s'])} /s")
                    if s["measured"] != "all":
                        row(f"Of these, {s['measured']}", f"{_fmt(el['measured_rate_per_s'])} /s")
                    row("Cline's safe distance", "kept" if el["safe"] else "not kept",
                        classes="" if el["safe"] else "ps-warn")
                else:
                    row("Centre at θ", f"{el['theta']:.1f}°")
                    row("Half-angle", f"{el['half_angle_deg']:.1f}°")
            if s["kind"] == "detector" and s["on_axis"]:
                place_inputs((("distance", "Distance", "mm"),))
            else:
                place_inputs((("theta", "θ", "deg"), ("phi", "φ", "deg"), ("distance", "Distance", "mm")))
            if s["kind"] == "gamma":
                gamma_response(key, s)
            if s["kind"] == "detector":
                spectrum_box = ui.column().classes("w-full")

                async def spectrum() -> None:
                    spectrum_box.clear()
                    with spectrum_box:
                        ui.spinner(size="md")
                    fig = await run.io_bound(figure_detector_spectrum, P(), s["name"], state["events"])
                    spectrum_box.clear()
                    with spectrum_box:
                        ui.plotly(themed(fig, state["theme"])).classes("w-full")

                ui.button("Simulate its spectrum", icon="play_arrow", on_click=spectrum).props(
                    "dense flat no-caps").classes("mt-1")

    def gamma_response(key: str, s: dict) -> None:
        """The efficiency curve of the selected γ-ray detector, what it is built from, and a calibration-source
        run that measures it."""
        from .. import response as _response

        src = scene_state["source"]
        resp = s["response"]
        ui.label("Efficiency against energy").classes("ps-section mt-2")
        where = {"model": f"Typical response of a {resp['material']} crystal, not a calibration of this detector.",
                 "fixed": "The single efficiency of the setup, used at every energy.",
                 "curve": "The measured efficiency curve of the setup."}[resp["source"]]
        through = ", ".join(f"{a['material']} {a['g_cm2']:.3g} g/cm² ({a['origin']})" for a in resp["absorbers"])
        ui.label(where + (" The crystal is taken as long as it is wide." if resp["assumed_length"] else "")
                 + f" On the way: {through}.").classes("text-xs ps-muted")
        curve_box = ui.column().classes("w-full")
        spectrum_box = ui.column().classes("w-full")

        def draw(run_=None) -> None:
            points = run_.efficiency_points(s["name"]) if run_ is not None else None
            curve_box.clear()
            with curve_box:
                ui.plotly(themed(figure_efficiency(P(), key, points), state["theme"])).classes("w-full")
            spectrum_box.clear()
            if run_ is not None:
                with spectrum_box:
                    ui.plotly(themed(figure_source_spectrum(run_, s["name"]), state["theme"])).classes("w-full")
                    ui.label(f"{run_.source.nuclide} at the target position, {run_.decays:.3g} decays. The points "
                             "above are the peak areas of its strong lines divided by the γ rays emitted.").classes(
                        "text-xs ps-muted")

        async def run_source() -> None:
            src.update(nuclide=nuclide.value, activity=activity.value, time=duration.value)
            try:
                run_ = await run.io_bound(P().source_run, src["nuclide"], src["activity"], src["time"])
            except ValueError as err:
                ui.notify(str(err), type="negative", multi_line=True)
                return
            src["run"] = True
            draw(run_)

        try:
            draw(P().source_run(src["nuclide"], src["activity"], src["time"]) if src["run"] else None)
        except ValueError:
            draw()
        ui.label("Calibration source, in place of the beam").classes("ps-section mt-2")
        with ui.row().classes("w-full no-wrap gap-2 items-end"):
            nuclide = ui.select(_response.source_names(), value=src["nuclide"], label="Source").props(
                "dense outlined").classes("min-w-0")
            activity = ui.input("Activity", value=src["activity"]).props("dense outlined").classes("min-w-0")
            duration = ui.input("Time", value=src["time"]).props("dense outlined").classes("min-w-0")
        ui.button("Run the source", icon="play_arrow", on_click=run_source).props("dense flat no-caps")

    def tracks_toolbar() -> None:
        """Animated tracks of a sample of events: how many, which, and the play controls."""
        ts = scene_state["tracks"]
        note = ui.label("").classes("text-xs ps-muted")

        async def show() -> None:
            view = scene_state["view"]
            scene_state["show_tracks"] = show
            if view is None:
                return
            try:
                r = await run.io_bound(P().tracks, int(count.value), ts["select"], ts["weighted"], 1,
                                       max(state["events"], 200_000))
            except ValueError as err:
                ui.notify(str(err), type="warning")
                return
            view.show_tracks(r["tracks"], ts["speed"])
            view.control_tracks(playing=ts["playing"])
            note.set_text(r["description"])
            channels = {"all": "all events", "coincidences": "particle–γ coincidences"}
            channels.update({ch: ch for ch in r["channels"]})
            which.set_options(channels, value=ts["select"] if ts["select"] in channels else "all")

        def toggle_play() -> None:
            ts["playing"] = not ts["playing"]
            play.set_text("Pause" if ts["playing"] else "Play")
            play.props(f"icon={'pause' if ts['playing'] else 'play_arrow'}")
            if scene_state["view"] is not None:
                scene_state["view"].control_tracks(playing=ts["playing"])

        def set_speed(e) -> None:
            ts["speed"] = float(e.value)
            if scene_state["view"] is not None:
                scene_state["view"].control_tracks(speed=ts["speed"])

        with ui.row().classes("w-full items-center gap-2"):
            count = ui.number("Tracks", value=ts["n"], min=1, max=200, step=5).props("dense outlined").classes(
                "w-24").tooltip("Up to about 60 stay smooth on a laptop")
            which = ui.select({"all": "all events", "coincidences": "particle–γ coincidences"}, value=ts["select"],
                              label="Which events", on_change=lambda e: ts.update(select=e.value)).props(
                "dense outlined").classes("w-52")
            ui.select({True: "as in a run (by rate)", False: "as generated (rare ones show)"}, value=ts["weighted"],
                      label="Sample", on_change=lambda e: ts.update(weighted=e.value)).props(
                "dense outlined").classes("w-56")
            ui.button("Show tracks", icon="timeline", on_click=show).props("dense flat no-caps")
            scene_state["show_tracks"] = show
            play = ui.button("Pause", icon="pause", on_click=toggle_play).props("dense flat no-caps")
            ui.slider(min=0.2, max=4.0, step=0.2, value=ts["speed"], on_change=set_speed).classes("w-32").tooltip(
                "Speed")
            ui.button("Clear", on_click=lambda: (scene_state["view"].clear_tracks(), note.set_text(""))).props(
                "dense flat no-caps")

    @ui.refreshable
    def scene_full(on: bool) -> None:
        """The scene on the whole screen: the setup panel, the status strip, the tabs and the panels beside and
        below the scene are hidden, and the scene fills the window; off brings them back."""
        scene_state["full"] = on
        drawer = state.get("drawer")
        if drawer is not None:
            (drawer.hide if on else drawer.show)()
        for el in (state.get("strip_box"), state.get("tabs"), *scene_state.get("extras", [])):
            if el is not None:
                el.set_visibility(not on)
        geometry_panel.refresh()
        # The scene follows its element's size once the window says it changed: once right away, and once more
        # after the drawer has slid out of the way.
        for delay in (0.15, 0.6):
            ui.timer(delay, lambda: ui.run_javascript("window.dispatchEvent(new Event('resize'))"), once=True)

    @ui.refreshable
    def geometry_panel():
        from ..scene_view import SceneView

        full = scene_state["full"]
        with ui.row().classes("w-full items-start gap-3").style("flex-wrap: wrap"):
            with ui.column().classes("ps-plate gap-1").style(
                    "flex: 1 1 100%; min-width: 0" if full else "flex: 1 1 520px; min-width: 0"):
                with ui.row().classes("w-full items-center gap-1"):
                    ui.toggle({"angle": "Drag changes the angle", "distance": "the distance"}, value="angle",
                              on_change=lambda e: scene_state["view"].set_mode(e.value)).props(
                        "dense no-caps unelevated").tooltip(
                        "What dragging a selected detector changes. A detector around the beam always slides "
                        "along it.")
                    ui.space()
                    for name, label in (("default", "3D"), ("side", "Side"), ("top", "Top"), ("beam", "Along beam")):
                        ui.button(label, on_click=lambda n=name: scene_state["view"].look(n)).props(
                            "dense flat no-caps")
                    ui.button("Paper figure", icon="article", on_click=lambda: open_export("geometry", {})).props(
                        "dense flat no-caps").tooltip("Export the layout as a figure in a journal's style")
                    ui.button(icon="fullscreen_exit" if full else "fullscreen",
                              on_click=lambda: scene_full(not scene_state["full"])).props(
                        "dense flat round").tooltip("Back to the page (Esc)" if full
                                                    else "The scene on the whole screen")
                view = SceneView(P, state["theme"], scene_selected, scene_live, scene_moved,
                                 height="calc(100vh - 175px)" if full else None)
                view.on_track = scene_track
                scene_state["view"] = view
                tracks_toolbar()
                if scene_state["key"] is not None:
                    view.select(scene_state["key"], scene_state["element"], notify=False)
                ui.label(("Full screen: the setup panel, the tabs and the numbers beside the scene are hidden; "
                          "Esc or the button brings them back. " if full else "")
                         + "To scale; the numbers along the beam are mm from the target. Click a detector, a ring "
                         "or a crystal to select it; drag a selected detector to move it. Drag the background to "
                         "turn the view, scroll to zoom.").classes(
                    "text-xs ps-muted")
            if not full:
                with ui.column().classes("gap-2").style("flex: 0 1 360px; min-width: 300px"):
                    selection_panel()
        if full:
            ui.keyboard(on_key=lambda e: scene_full(False) if e.key.escape and e.action.keydown else None)
        else:
            geometry_tables()

    @ui.refreshable
    def geometry_tables():
        rows = [{"detector": d["name"], "omega": _fmt(d["solid_angle_msr"]),
                 "theta": f"{d['theta_range'][0]:.1f}–{d['theta_range'][1]:.1f}",
                 "phi": f"{d['phi_range'][0]:.1f}–{d['phi_range'][1]:.1f}", "segments": d["segments"]}
                for d in P().geometry()["detectors"]]
        ui.table(columns=columns((("detector", "Particle detector"), ("omega", "Ω (msr)"), ("theta", "θ (deg)"),
                                  ("phi", "φ (deg)"), ("segments", "Segments"))), rows=rows).props("dense flat")
        gammas = P().geometry()["gamma_detectors"]
        if gammas:
            ui.table(columns=columns((("detector", "γ-ray detector"), ("theta", "θ (deg)"),
                                      ("distance", "Distance (mm)"), ("half", "Half-angle (deg)"))),
                     rows=[{"detector": g["name"], "theta": f"{g['theta']:.1f}", "distance": f"{g['distance_mm']:g}",
                            "half": f"{g['half_angle_deg']:.1f}"} for g in gammas]).props("dense flat")

    @ui.refreshable
    def kinematics_panel():
        plot(figure_kinematics(P()), "kinematics")

    @ui.refreshable
    def rates_panel():
        r = P().rates()
        rows = [{"detector": x["detector"], "omega": _fmt(x["solid_angle_msr"]),
                 "theta": f"{x['theta_range'][0]:.1f}–{x['theta_range'][1]:.1f}", "rate": _fmt(x["rate_per_s"]),
                 "counts": _fmt(x["counts_in_run"]),
                 "error": f"{100 * x['relative_error']:.2g}%" if math.isfinite(x["relative_error"]) else "—",
                 "time": _time(x["beam_time_s"])} for x in r["rows"]]
        what = {"all": "counts", "excitations": "excitations", "coincidences": "particle–γ coincidences"}[
            r["measured"]]
        wanted = f" for {r['counts_wanted']} {what}" if r["counts_wanted"] else ""
        cols = [("detector", "Detector"), ("omega", "Ω (msr)"), ("theta", "θ (deg)"), ("rate", "Rate (1/s)")]
        if r["measured"] != "all":
            for row, x in zip(rows, r["rows"]):
                row["exc"] = _fmt(x["excitation_per_s"])
                row["coinc"] = _fmt(x["coincidence_per_s"]) if x["coincidence_per_s"] is not None else "—"
            cols += [("exc", "Excitations (1/s)"), ("coinc", "With γ ray (1/s)")]
        cols += [("counts", f"{what.capitalize()} in run"), ("error", "Stat. error"), ("time", f"Beam time{wanted}")]
        ui.table(columns=columns(cols), rows=rows).props("dense flat")
        names = [x["detector"] for x in r["rows"]]

        @ui.refreshable
        def strips(name):
            plot(figure_strips(P(), name), "strips", detector=name)

        ui.select(names, value=names[0], label="Rates per strip of",
                  on_change=lambda e: strips.refresh(e.value)).props("dense outlined").classes("w-48")
        strips(names[0])
        with ui.card().classes("w-full mt-2"):
            ui.label("Sweep one parameter").classes("font-semibold")
            with ui.row().classes("items-end gap-2"):
                par = ui.select(["beam energy", "target thickness", "detector angle"], value="beam energy",
                                label="Vary").props("dense outlined").classes("w-44")
                vals = ui.input("Values (comma-separated, with units)",
                                value="4 MeV, 5 MeV, 6 MeV, 7 MeV").props("dense outlined").classes("w-72")
                qty = ui.select(["rate", "beam time", "peak energy", "peak width"], value="rate",
                                label="Show").props("dense outlined").classes("w-36")
                det = ui.select(names, value=names[0], label="Detector").props("dense outlined").classes("w-36")

                async def do_sweep():
                    values = [v.strip() for v in vals.value.split(",") if v.strip()]
                    try:
                        s = await run.io_bound(P().sweep, par.value, values, qty.value, det.value)
                    except (SetupError, ValueError, KeyError) as err:
                        ui.notify(str(err), type="negative", multi_line=True)
                        return
                    sweep_out.clear()
                    with sweep_out:
                        plot(figure_sweep(s), "sweep", sweep=s)

                ui.button("Run sweep", on_click=do_sweep)
            sweep_out = ui.column().classes("w-full")

    @ui.refreshable
    def energy_loss_panel():
        el = P().energy_loss()
        ui.table(columns=columns((("layer", "Layer"), ("material", "Material"), ("t", "Thickness (mg/cm²)"),
                                  ("ein", "Beam in (MeV)"), ("eout", "Beam out (MeV)"), ("loss", "Loss (keV)"),
                                  ("strag", "Straggling FWHM (keV)"))),
                 rows=[{"layer": x["layer"], "material": x["material"], "t": _fmt(x["thickness_mg_cm2"]),
                        "ein": f"{x['energy_in_mev']:.4f}", "eout": f"{x['energy_out_mev']:.4f}",
                        "loss": f"{x['loss_mev'] * 1e3:.1f}", "strag": f"{x['straggling_fwhm_mev'] * 1e3:.1f}"}
                       for x in el["layers"]]).props("dense flat")
        plot(figure_energy_loss(P()), "energy_loss")
        ui.table(columns=columns((("detector", "Detector"), ("e", "Scattered beam (MeV)"),
                                  ("dl", "After dead layer (MeV)"), ("pt", "Punch-through above (MeV)"))),
                 rows=[{"detector": x["detector"], "e": f"{x['ejectile_energy_mev']:.4f}",
                        "dl": f"{x['after_dead_layer_mev']:.4f}", "pt": f"{x['punch_through_mev']:.4g}"}
                       for x in el["detectors"]]).props("dense flat")

    @ui.refreshable
    def spectra_panel():
        with ui.row().classes("items-end gap-2"):
            n = ui.number("Events", value=state["events"], format="%d", min=1000, step=100_000).props(
                "dense outlined").classes("w-40")
            seed = ui.number("Seed", value=1, format="%d", min=0).props("dense outlined").classes("w-28")

            async def simulate():
                state["events"] = int(n.value)
                fig = await run.io_bound(figure_spectra, P(), int(n.value), int(seed.value))
                plot_area.clear()
                with plot_area:
                    plot(fig, "spectra", events=int(n.value), seed=int(seed.value))

            ui.button("Simulate", icon="play_arrow", on_click=simulate)
        plot_area = ui.column().classes("w-full")
        with plot_area:
            plot(figure_spectra(P(), state["events"], 1), "spectra", events=state["events"], seed=1)
        gamma_block()

    def gamma_block() -> None:
        """Particle–γ coincidences: the matrix, and each γ-ray detector's raw and corrected spectrum."""
        if P().experiment.excitation is None or not P().experiment.gamma_detectors:
            return
        ui.label("γ rays in coincidence with the particles").classes("font-semibold mt-3")
        ui.label("Every excited event whose particle reached a detector emits its γ ray, with the angular "
                 "correlation, from the moving nucleus; the crystals record it with their response. The Doppler "
                 "correction uses the centre of the segment and of the crystal that fired. Random coincidences "
                 "come from the singles rates and the coincidence window, and the room background from the "
                 "[run] section. A clover with add-back or a BGO shield (its switches are in the setup panel) "
                 "also shows its spectrum without them.").classes("text-xs ps-muted")
        holder = ui.column().classes("w-full")

        async def build() -> None:
            holder.clear()
            with holder:
                ui.spinner(size="md")
            try:
                g = await run.io_bound(P().gamma_spectra, max(state["events"], 400_000), 1)
            except ValueError as err:
                holder.clear()
                with holder:
                    ui.label(str(err)).classes("text-sm ps-warn")
                return
            holder.clear()
            with holder:
                if not g["available"]:
                    ui.label(g["reason"]).classes("text-sm")
                    return
                co = g["coincidences"]
                cols = [("p", "Particle detector")] + [(f"g{i}", f"{n} (true / random)")
                                                        for i, n in enumerate(g["gamma_detectors"])]
                rows = [dict({"p": d}, **{f"g{i}": f"{_fmt(co['true'][(d, n)], 3)} / "
                                                   f"{_fmt(co['random'][(d, n)], 2)}"
                                          for i, n in enumerate(g["gamma_detectors"])})
                        for d in g["particle_detectors"]]
                ui.table(columns=columns(cols), rows=rows).props("dense flat")
                ui.label(f"Counts in the run, after a live fraction of {g['live_fraction']:.3f}; coincidence "
                         f"window {1e9 * g['window_s']:.0f} ns. {g['notes'][0]}").classes("text-xs ps-muted")
                sel = ui.select(g["gamma_detectors"], value=g["gamma_detectors"][0], label="γ-ray detector").props(
                    "dense outlined").classes("w-48")
                fig_box = ui.column().classes("w-full")

                def draw() -> None:
                    fig_box.clear()
                    with fig_box:
                        ui.plotly(themed(figure_gamma_spectra(P(), sel.value, max(state["events"], 400_000), 1),
                                         state["theme"])).classes("w-full")

                sel.on_value_change(lambda e: draw())
                draw()

        ui.button("Simulate the γ rays", icon="play_arrow", on_click=build).props("dense flat no-caps")

    def analysis_block() -> None:
        """From the γ-ray peak to B(E2): the steps, each adjustable, and what the beam time determines."""
        ui.label("Analysis: from the peak area to B(E2)").classes("font-semibold mt-3")
        ui.label("The simulated events are analysed as an experimentalist would: a gate on the inelastic particle "
                 "group, the Doppler correction, a peak fit over a line with the random coincidences subtracted, "
                 "the yield over the efficiency and the angular correlation, a normalisation, and B(E2) from the "
                 "first-order proportionality. Each systematic is found by changing its input by one standard "
                 "deviation and running again.").classes("text-xs ps-muted")
        with ui.row().classes("w-full items-end gap-2"):
            norm = ui.select({"rutherford": "to the elastic particles in the same rings",
                              "target": "to a known transition"}, value="rutherford",
                             label="Normalisation").props("dense outlined").classes("w-72")
            gate = ui.select({"inelastic": "inelastic group", "all": "every particle"}, value="inelastic",
                             label="Particle gate").props("dense outlined").classes("w-44")
            corr = ui.select({"emitter": "the excited nucleus", "projectile": "the projectile",
                              "recoil": "the recoil"}, value="emitter", label="Doppler correction").props(
                "dense outlined").classes("w-48")
            width = ui.number("Fit window (FWHM)", value=4.0, min=1.5, max=10, step=0.5).props(
                "dense outlined").classes("w-36")
            rings = ui.input("Rings or strips", placeholder="all, or 4-15").props("dense outlined").classes("w-36")
            named = ui.select({"": "none (the choices here)"} | {g.name: g.name for g in P().gates()}, value="",
                              label="Gate from the Data tab").props("dense outlined").classes("w-56").tooltip(
                "A named gate replaces the detectors, rings and particle gate here")
        with ui.row().classes("w-full items-end gap-2"):
            ref_e = ui.input("Reference transition", placeholder="328 keV").props("dense outlined").classes("w-40")
            ref_b = ui.input("Its B(E2↑)", placeholder="1.65 e2b2").props("dense outlined").classes("w-40")
            ref_u = ui.number("Its uncertainty (%)", value=3.0, min=0, max=100).props("dense outlined").classes(
                "w-40")
            eff_u = ui.number("Efficiency unc. (%)", value=5.0, min=0, max=100).props("dense outlined").classes(
                "w-40")
            prec = ui.number("Wanted precision (%)", value=5.0, min=0.1, max=100).props("dense outlined").classes(
                "w-40")
        out_box = ui.column().classes("w-full")

        def parse_rings(text: str):
            text = (text or "").strip()
            if not text or text == "all":
                return None
            picked = []
            for part in text.split(","):
                a, _, b = part.strip().partition("-")
                picked += list(range(int(a), int(b or a) + 1))
            return picked

        async def analyse() -> None:
            settings = {"normalisation": norm.value, "particle_gate": gate.value, "correction": corr.value,
                        "fit_half_width": float(width.value), "efficiency_unc": float(eff_u.value) / 100,
                        "wanted_precision": float(prec.value) / 100}
            try:
                settings["rings"] = parse_rings(rings.value)
            except ValueError:
                ui.notify("Rings: write 'all', a number, or a range such as 4-15", type="warning")
                return
            if norm.value == "target":
                if not ref_e.value or not ref_b.value:
                    ui.notify("Give the reference transition's energy and B(E2↑)", type="warning")
                    return
                settings["reference"] = {"energy": ref_e.value, "b_up": ref_b.value, "unc": float(ref_u.value) / 100}
            out_box.clear()
            with out_box:
                ui.spinner(size="md")
            r = await run.io_bound(P().analyse, settings, max(state["events"], 400_000), 1, named.value or None)
            out_box.clear()
            with out_box:
                show(r)

        def show(r: dict) -> None:
            if not r["available"]:
                ui.label(r["reason"]).classes("text-sm ps-warn")
                return
            for step in r["steps"]:
                with ui.row().classes("items-start no-wrap gap-2"):
                    ui.label(step["step"]).classes("text-xs ps-section w-40")
                    ui.label(step["text"]).classes("text-sm")
            b = r["b_e2fm4"]
            if b != b:
                ui.label("No result: " + "; ".join(r["notes"])).classes("text-sm ps-warn")
                return
            last = P().analysis(max(state["events"], 400_000), 1).history[-1]
            with ui.element("div").classes("ps-readouts mt-2"):
                for title, value, sub, key in (
                        ("B(E2↑)", f"{r['b_e2b2']:.4g} e²b²", f"{b:.4g} e²fm⁴ · {r['b_wu']:.3g} W.u.", "b_e2"),
                        ("Uncertainty", f"± {100 * r['total_unc']:.1f} %",
                         f"statistical {100 * r['statistical']:.1f} %, systematic {100 * r['systematic']:.1f} %; "
                         + (f"the Monte Carlo sample itself adds {100 * r['monte_carlo']:.1f} %"
                            if r.get("whole_run_statistical") is None else
                            f"statistics from {r['statistics_from']}; the whole run: "
                            f"{100 * r['whole_run_statistical']:.1f} %"), "uncertainty"),
                        ("Put in", f"{1e-4 * r['truth_e2fm4']:.4g} e²b²", f"pull {r['pull']:+.2f} σ", None),
                        ("Beam time", f"{r['hours_for_precision']:.3g} h",
                         f"for {100 * r['settings']['wanted_precision']:.3g} % statistics; "
                         f"{r['counts_per_shift']:.3g} counts in the peak per {SHIFT_H:g} h shift", None)):
                    with ui.element("div").classes("ps-readout"):
                        with ui.row().classes("items-center no-wrap gap-1"):
                            ui.label(title).classes("text-xs ps-muted")
                            if key:
                                ui.button(icon="help_outline", on_click=lambda k=key: open_explanation(
                                    k, result=last)).props("dense flat round size=xs").tooltip(
                                    "How this number is obtained")
                        ui.label(value).classes("ps-readout-value")
                        ui.label(sub).classes("text-xs ps-muted")
            with ui.row().classes("gap-1 mt-1 items-center"):
                ui.label("The chain, step by step:").classes("text-xs ps-muted")
                for key, label in (("area", "area"), ("yield", "yield"), ("mean_probability", "⟨P⟩"),
                                   ("b_e2", "B(E2)"), ("weisskopf", "W.u."), ("beta2", "β₂"), ("q0", "Q₀"),
                                   ("lifetime", "lifetime")):
                    ui.button(label, on_click=lambda k=key: open_explanation(k, result=last)).props(
                        "dense flat no-caps size=sm")
            with ui.row().classes("items-start gap-4 mt-2"):
                ui.table(columns=columns((("source", "Systematic"), ("value", "Relative (%)"))),
                         rows=[{"source": k, "value": f"{100 * v:.2f}"} for k, v in r["budget"].items()]).props(
                    "dense flat")
                sh = r["shape"]
                if sh:
                    rows = [{"q": "β₂", "v": f"{sh['beta2']:.3f}"}, {"q": "B(E2↑) in W.u.", "v": f"{r['b_wu']:.3g}"},
                            {"q": "Q₀ of a rigid rotor", "v": f"{sh['q0_efm2']:.3g} e fm²"},
                            {"q": "Q_s(2⁺) of that rotor", "v": f"{sh['qs_2plus_efm2']:.3g} e fm²"}]
                    if "e4_over_e2" in sh:
                        rows.append({"q": "E(4⁺)/E(2⁺)", "v": f"{sh['e4_over_e2']:.3f}"})
                    ui.table(columns=columns((("q", "Shape"), ("v", "Value"))), rows=rows).props("dense flat")
            if r["shape"]:
                ui.label(r["shape"]["note"]).classes("text-xs ps-muted")
            for note in r["notes"]:
                ui.label(note).classes("text-xs ps-muted")

        ui.button("Analyse", icon="functions", on_click=analyse).props("dense flat no-caps")

    def alignment_block() -> None:
        """A misplaced target: what the analysis sees, and how it finds the offset."""
        ui.label("A misplaced target").classes("font-semibold mt-3")
        ui.label("The analysis may assume the target somewhere else than it is. The Doppler correction then uses "
                 "wrong angles: each corrected peak shifts and broadens, ring by ring and crystal by crystal. The "
                 "diagnostic plot below, the corrected centroid against ring, is what one uses on real data to "
                 "find a misplaced target; the fit gives the offset back.").classes("text-xs ps-muted")
        current = P().run
        default_offset = float(current.summary.get("assumed_offset_mm", 2.0)) if current is not None and \
            current.kind == "alignment" else 2.0
        if current is not None and current.kind == "alignment":
            ui.label(f"Run {current.number} is an alignment run: the target assumed {default_offset:+g} mm off.").classes(
                "text-sm ps-accent")
        with ui.row().classes("items-end gap-2"):
            off = ui.number("Assumed offset along the beam (mm)", value=default_offset, min=-20, max=20,
                            step=0.5).props(
                "dense outlined").classes("w-64")
            fit_too = ui.checkbox("Fit the offset back (about 15 s)", value=True)
        box = ui.column().classes("w-full")

        async def check() -> None:
            box.clear()
            with box:
                ui.spinner(size="md")
            r = await run.io_bound(P().alignment, float(off.value), max(state["events"], 400_000), 1,
                                   bool(fit_too.value))
            if scene_state["view"] is not None and scene_state["view"].alive:
                scene_state["view"].show_ghost(float(off.value))
            box.clear()
            with box:
                if not r["available"]:
                    ui.label(r["reason"]).classes("text-sm ps-warn")
                    return
                o = r["overlay"]
                ui.label(f"With the target assumed {r['offset_mm']:+g} mm off, the corrected peak shifts by "
                         f"{o['shift_kev']:+.2f} keV and its width changes by {o['broadening_kev']:+.2f} keV "
                         f"(FWHM {o['fwhm_true_kev']:.1f} → {o['fwhm_assumed_kev']:.1f} keV).").classes("text-sm")
                ui.plotly(themed(figure_overlay(r), state["theme"])).classes("w-full")
                ui.plotly(themed(figure_alignment(P(), r), state["theme"])).classes("w-full")
                if P().run is not None and P().run.kind != "source" and P().run.gammas() is not None:
                    exc = P()._data_experiment().excitation
                    corr_key = "recoil" if exc.excite == "target" else "projectile"
                    view = dv.gamma_vs_crystal(P().run, corr_key, offset_mm=float(off.value))
                    ui.label("The run's γ rays, corrected with the assumed geometry, crystal by crystal (as on "
                             "the Data tab)").classes("ps-section mt-2")
                    ui.plotly(themed(figure_gamma_crystal(view), state["theme"])).classes("w-full")
                if "fit" in r:
                    f = r["fit"]
                    ui.label(f"Fitted: the target is {f['offset_mm']:+.2f} ± {f['uncertainty_mm']:.2f} mm from "
                             f"where the analysis assumed it (truth {f['truth_mm']:+g} mm; χ² {f['chi2']:.0f} for "
                             f"{f['degrees_of_freedom']} points, judged against a simulation with the assumed "
                             "geometry). The scene shows the assumed target as a faint outline.").classes("text-sm")

        ui.button("Check the alignment", icon="straighten", on_click=check).props("dense flat no-caps")

    @ui.refreshable
    def trajectories_panel():
        t = P().trajectories()
        ui.label(f"{P().experiment.beam.nuclide} on {t['target']}: head-on distance d₀ = {t['d0_fm']:.2f} fm, "
                 f"nuclear range {t['interaction_radius_fm']:.1f} fm, grazing angle "
                 f"{t['grazing_angle_deg']:.1f}° (CM).").classes("text-sm")
        plot(figure_trajectories(P()), "trajectories")

    def lookup_levels(role: str):
        def run() -> None:
            try:
                ok = P().lookup_levels(role)
            except (LookupError, ValueError) as e:
                ui.notify(str(e), type="warning", multi_line=True)
                return
            changed(ok)
        return run

    def remove_levels(role: str):
        return lambda: changed(P().remove_levels(role))

    def use_state(role: str, level: int, multipolarity: str):
        def run() -> None:
            changed(P().use_state(role, level, multipolarity))
            setup_panel.refresh()
        return run

    def matrix_element_input(role: str, m: dict) -> None:
        """A matrix element as an editable field: a new value, with its unit, becomes the user's."""
        lam = int(m["multipolarity"][1])
        unit = ("efm", "efm2", "efm3")[lam - 1]
        shown = f"{_fmt(m['value_efm'])} {unit}"
        box = ui.input(value=shown).props("dense outlined").classes("w-36")
        if m["unc_efm"]:
            box.tooltip(f"± {_fmt(m['unc_efm'], 2)} {unit}")

        def apply() -> None:
            text = (box.value or "").strip()
            if text == shown:
                return
            try:
                ok = P().set_matrix_element(role, m["from"], m["to"], m["multipolarity"], text)
            except ValueError as e:
                ui.notify(str(e), type="negative", multi_line=True)
                box.value = shown
                return
            changed(ok)

        box.on("blur", apply)
        box.on("keydown.enter", apply)

    def levels_block() -> None:
        """The level schemes of the target and the beam: look-up, diagram, and matrix elements with their source."""
        info = P().levels()
        ui.label("Level schemes").classes("font-semibold")
        if not info["ensdf"]:
            ui.label(f"There is no local copy of ENSDF in {info['ensdf_folder']}. Download one with "
                     "scripts/fetch_ensdf.py, or type the level scheme in the setup file.").classes(
                "text-xs ps-muted")
        for role, title in (("target", "Target"), ("beam", "Beam")):
            nuclide, table = info["nuclides"][role], info[role]
            with ui.row().classes("items-center gap-2"):
                ui.label(f"{title}: {nuclide or 'not set'}").classes("text-sm")
                look = ui.button("Look up in ENSDF", icon="search", on_click=lookup_levels(role)).props(
                    "dense flat no-caps")
                if not info["ensdf"] or nuclide is None:
                    look.disable()
                if table:
                    ui.button("Remove", icon="delete", on_click=remove_levels(role)).props(
                        "dense flat no-caps color=negative")
            if not table:
                continue
            ui.label(table["reference"]).classes("text-xs ps-muted")
            for note in table["notes"]:
                ui.label(note).classes("text-xs ps-muted")
            ui.plotly(themed(figure_levels(table), state["theme"])).classes("w-full")
            names = [f"{lv['jpi']} {lv['energy_kev']:g} keV" for lv in table["levels"]]
            ui.label("Matrix elements").classes("text-sm font-semibold")
            ui.label("Sizes come from ENSDF's transition strengths or half-lives; ENSDF gives no signs, so each "
                     "sign is assumed. Type a new value with its unit (efm2 or eb for E2) to use your own, with a minus "
                     "sign where you want one. \"Plan\" sets the reaction to excite that state.").classes(
                "text-xs ps-muted")
            with ui.grid(columns=6).classes("items-center gap-x-4 gap-y-0 text-sm"):
                for head in ("Between", "", "Matrix element", "B↑ or Q", "Source", ""):
                    ui.label(head).classes("text-xs ps-muted")
                for m in table["matrix_elements"]:
                    lam = int(m["multipolarity"][1])
                    ui.label(names[m["from"]] if m["from"] == m["to"] else f"{names[m['from']]} ↔ {names[m['to']]}")
                    ui.label(m["multipolarity"])
                    matrix_element_input(role, m)
                    if m["q_efm2"] is not None:
                        ui.label(f"Q = {_fmt(m['q_efm2'] / 100)} e b")
                    elif m["b_up_wu"] is not None:
                        ui.label(f"{_fmt(m['b_up_e2fm'])} e²fm^{2 * lam} ({_fmt(m['b_up_wu'], 3)} W.u.)")
                    else:
                        ui.label("")
                    ui.label(m["source"]).tooltip(m["note"] or m["source"])
                    if m["from"] == 0 and m["to"] != 0:
                        ui.button("Plan", on_click=use_state(role, m["to"], m["multipolarity"])).props(
                            "dense flat no-caps")
                    else:
                        ui.label("")
        ui.separator()

    @ui.refreshable
    def gamma_panel(with_levels: bool = True):
        if with_levels:
            levels_block()
        g = P().gamma()
        if not g["available"]:
            ui.markdown("This setup has no excited state. To plan Coulomb excitation, set **Reaction** to "
                    "*Coulomb excitation* in the setup panel, enter the state's energy and B(Eλ↑), and add "
                    "γ-ray detectors.")
            return
        s = g["state"]
        ui.label(f"{s['excite'].capitalize()} excited to {s['energy_kev']:g} keV ({s['multipolarity']}, "
                 f"B↑ = {s['b_up_e2fm']:.4g} e²fm^{2 * int(s['multipolarity'][1])}): ξ = {g['xi']:.2f}, "
                 f"η = {g['eta']:.1f}, total {g['total_mb']:.3g} mb; safe up to "
                 f"{g['max_safe_angle']:.0f}° CM.").classes("text-sm")
        plot(figure_excitation(P()), "gamma")
        ui.table(columns=columns((("detector", "Particle detector"), ("rate", "Excitation events (1/s)"))),
                 rows=[{"detector": k, "rate": _fmt(v)} for k, v in g["rates"].items()]).props("dense flat")
        if g["particles"]:
            ui.label("Particle energies: elastic and after exciting the state").classes("font-semibold mt-2")
            ui.label("At each particle detector's smallest, central and largest angle, at the reaction point "
                     "(before energy loss in the target; the Spectra tab includes it). β is the speed of the "
                     "excited nucleus in that event, which the Doppler correction needs.").classes(
                "text-xs ps-muted")
            ui.table(columns=columns((("detector", "Detector"), ("particle", "Particle"), ("theta", "θ lab (deg)"),
                                      ("el", "Elastic (MeV)"), ("ex", "Excited (MeV)"), ("diff", "Difference (MeV)"),
                                      ("beta", "β excited"))),
                     rows=[{"detector": r["detector"], "particle": f"{r['nuclide']} ({r['particle']})",
                            "theta": f"{r['theta_lab']:.1f}", "el": f"{r['elastic_mev']:.2f}",
                            "ex": f"{r['excited_mev']:.2f}", "diff": f"{r['difference_mev']:.3f}",
                            "beta": f"{r['beta_excited']:.4f}"} for r in g["particles"]]).props("dense flat")
        populations_block()
        if g["correlation"]:
            ui.label("γ rays in coincidence: the angular correlation").classes("font-semibold mt-2")
            ui.label("How many γ rays each crystal sees when the particle is in each particle detector, relative "
                     "to γ rays sent evenly in all directions (1 = even). The excited nucleus is left oriented by "
                     "the collision, so its γ rays favour some directions; they are also thrown forward by its "
                     "motion. The coincidence rates use these factors.").classes("text-xs ps-muted")
            ui.toggle({"correlated": "Correlated with the particle", "isotropic": "Isotropic"}, value=g["emission"],
                      on_change=lambda e: edit("reaction", "emission", e.value)).props("dense no-caps unelevated")
            names = list(dict.fromkeys((r["gamma_detector"], r["crystal"]) for r in g["correlation"]))
            cols = [("p", "Particle detector")] + [(f"c{i}", f"{n} {c}".strip()) for i, (n, c) in enumerate(names)]
            by = {}
            for r in g["correlation"]:
                by.setdefault(r["particle_detector"], {})[(r["gamma_detector"], r["crystal"])] = r["factor"]
            ui.table(columns=columns(cols),
                     rows=[dict({"p": p_}, **{f"c{i}": f"{f.get(n, 1.0):.2f}" for i, n in enumerate(names)})
                           for p_, f in by.items()]).props("dense flat")
        if g["doppler"]:
            ui.label("γ rays: Doppler-shifted energy and width").classes("font-semibold mt-2")
            ui.table(columns=columns((("p", "Particle detector"), ("g", "γ detector"), ("mean", "E_γ (keV)"),
                                      ("shift", "Shift (keV)"), ("fwhm", "FWHM (keV)"))),
                     rows=[{"p": r["particle_detector"], "g": r["gamma_detector"], "mean": f"{r['mean_kev']:.2f}",
                            "shift": f"{r['shift_kev']:+.2f}", "fwhm": f"{r['fwhm_kev']:.2f}"}
                           for r in g["doppler"]]).props("dense flat")
        else:
            ui.label("Add γ-ray detectors in the setup panel for Doppler shifts.").classes("text-sm")

    def multistep_block() -> None:
        """All orders: the coupled equations on the level scheme, the reorientation comparison, GOSIA input."""
        if not P().experiment.levels:
            return
        ui.label("All orders: multi-step excitation and reorientation").classes("font-semibold mt-2")
        ui.label("The coupled equations follow every substate of every level along the orbit, so excitation in "
                 "several steps, the reorientation effect of the quadrupole moments and interference between paths "
                 "are included. The same 2⁺ state is then run with its quadrupole moment at the prolate rotor "
                 "value, zero and the oblate value, to see whether the planned run can tell them apart.").classes(
            "text-xs ps-muted")
        holder = ui.column().classes("w-full")

        async def solve() -> None:
            holder.clear()
            with holder:
                ui.spinner(size="md")
                ui.label("Solving on a grid of angles and energies; a few levels take seconds, twenty take a "
                         "minute.").classes("text-xs ps-muted")
            m = await run.io_bound(P().multistep)
            holder.clear()
            with holder:
                if not m["available"]:
                    ui.label(m["reason"]).classes("text-sm ps-warn")
                    return
                t = m["beam_time_s"]
                dets = m["detectors"]
                cols = [("t", "γ ray"), ("e", "keV")] + [(f"d{i}", f"{d} (all / first order)")
                                                         for i, d in enumerate(dets)]
                rows = [dict({"t": x["label"], "e": f"{x['energy_kev']:.1f}"},
                             **{f"d{i}": f"{_fmt(x['all_orders'][d] * t, 3)} / {_fmt(x['first_order'][d] * t, 3)}"
                                for i, d in enumerate(dets)}) for x in m["transitions"]]
                ui.label("γ rays in the run, with all orders and with first order").classes("ps-section mt-1")
                ui.table(columns=columns(cols), rows=rows).props("dense flat")
                if "shapes" in m:
                    sh = m["shapes"]
                    q = sh["q_efm2"]
                    ui.label(f"Prolate or oblate: the {m['nuclide']} 2⁺ state with Q(2⁺) = {q['prolate']:+.1f} "
                             f"(prolate rotor), 0 and {q['oblate']:+.1f} e fm² (oblate)").classes("ps-section mt-2")
                    ui.label(("The planned run can tell prolate from oblate: " if sh["separable"] else
                              "The planned run cannot tell prolate from oblate: ")
                             + f"over all rings the two differ by {sh['total_difference_sigma']:.1f} standard "
                             "deviations of the counts.").classes("text-sm " + ("ps-ok" if sh["separable"] else
                                                                                   "ps-warn"))
                    ui.table(columns=columns((("d", "Detector"), ("r", "Ring / strip"), ("p", "Prolate"),
                                              ("s", "Q = 0"), ("o", "Oblate"), ("u", "± counts"), ("z", "Δ / σ"))),
                             rows=[{"d": r_["detector"], "r": r_["ring"] + 1, "p": _fmt(r_["counts"]["prolate"], 3),
                                    "s": _fmt(r_["counts"]["spherical"], 3), "o": _fmt(r_["counts"]["oblate"], 3),
                                    "u": _fmt(r_["uncertainty"], 2), "z": f"{r_['difference_sigma']:.1f}"}
                                   for r_ in sh["rings"]]).props("dense flat").classes("max-h-80")
                for note in m["notes"][:4]:
                    ui.label(note).classes("text-xs ps-muted")
                ui.button("GOSIA input file for this setup", icon="download",
                          on_click=lambda: ui.download.content(m["gosia_input"].encode("utf-8"),
                                                               "gosia.inp")).props("dense flat no-caps").tooltip(
                    "Written from the GOSIA manual's format and not checked against a run: GOSIA is not on this "
                    "machine.")

        ui.button("Solve with all orders", icon="calculate", on_click=solve).props("dense flat no-caps")
        fit_box = ui.column().classes("w-full")

        async def fit_run() -> None:
            fit_box.clear()
            with fit_box:
                ui.spinner(size="md")
            try:
                f = await run.io_bound(P().fit_matrix_elements)
            except (ValueError, KeyError) as err:
                fit_box.clear()
                with fit_box:
                    ui.label(str(err)).classes("text-sm ps-warn")
                return
            fit_box.clear()
            with fit_box:
                for k, v in f["values"].items():
                    unc = f["uncertainties"].get(k)
                    ui.label(f"⟨{k[1]}‖E2‖{k[0]}⟩ fitted to the run's yields: {v:.4g}"
                             + (f" ± {unc:.2g}" if unc is not None else "") + f" (χ² {f['chi2']:.2g})").classes(
                        "text-sm")

        if P().run is not None and P().run.kind != "source":
            ui.button("Fit the matrix element to this run's yields", icon="tune", on_click=fit_run).props(
                "dense flat no-caps").tooltip("The yields are the run's peaks per particle detector, over the "
                                              "efficiency and the correlation")

    def populations_block() -> None:
        """Every level of the scheme that first-order excitation reaches, and the γ rays that follow."""
        pop = P().populations()
        if not pop["available"] or not pop["levels"]:
            return
        ui.label(f"{pop['nuclide']}: all levels excited from the ground state").classes("font-semibold mt-2")
        ui.label(f"First-order cross sections over all scattering angles at {pop['beam_energy_mev']:.2f} MeV, "
                 "from the matrix elements of the level scheme. A level is also fed by the decay of the levels "
                 "above it; internal conversion takes its share of each transition. Excitation in two steps is "
                 "not included.").classes("text-xs ps-muted")
        with ui.row().classes("items-start gap-4"):
            ui.table(columns=columns((("lev", "Level"), ("e", "Energy (keV)"), ("d", "Excited directly (mb)"),
                                      ("p", "With feeding (mb)"))),
                     rows=[{"lev": x["label"], "e": f"{x['energy_kev']:.1f}", "d": _fmt(x["direct_mb"]),
                            "p": _fmt(x["populated_mb"])} for x in pop["levels"]]).props("dense flat")
            ui.table(columns=columns((("t", "γ ray"), ("e", "Energy (keV)"), ("g", "γ rays (mb)"))),
                     rows=[{"t": x["label"], "e": f"{x['energy_kev']:.1f}", "g": _fmt(x["gamma_mb"])}
                           for x in pop["gammas"][:12]]).props("dense flat")
        for note in pop["notes"][:4]:
            ui.label(note).classes("text-xs ps-muted")

    @ui.refreshable
    def report_panel():
        ui.markdown("The report collects the setup, every warning, the detector table, kinematics, peaks, energy "
                    "loss, spectra and each model's validation status. The zip holds `report.html` (print it to "
                    "PDF from the browser), CSV tables, figures (PNG and PDF), the setup file and, with uproot "
                    "installed, `events.root`.")

        async def download():
            note = ui.notification("Building the report…", spinner=True, timeout=None)
            try:
                data = await run.io_bound(report_zip, P(), 1, 200_000, state["journal"], state["width"])
            except ImportError as err:
                ui.notify(f"The report export is not available in this version: {err}", type="warning")
                return
            finally:
                note.dismiss()
            ui.download.content(data, "physim-report.zip")

        def pick(key: str, value) -> None:
            state[key] = value
            if state["width"] not in paper.JOURNALS[state["journal"]].widths:
                state["width"] = "single"
            report_panel.refresh()

        ui.label("The report's figures are drawn in a journal's style, at its column width.").classes(
            "text-sm ps-muted")
        with ui.row().classes("items-end gap-2"):
            ui.select({k: j.label for k, j in paper.JOURNALS.items()}, value=state["journal"], label="Journal style",
                      on_change=lambda e: pick("journal", e.value)).props("dense outlined").classes("w-80")
            ui.select({k: f"{paper.WIDTHS[k]} · {mm:g} mm"
                       for k, mm in paper.JOURNALS[state["journal"]].widths.items()}, value=state["width"],
                      label="Figure width", on_change=lambda e: pick("width", e.value)).props(
                "dense outlined").classes("w-56")
        ui.button("Build and download the report", icon="description", on_click=download)
        record_block()

    def record_block() -> None:
        """The run record: every number with its explanation, viewed here and saved as a page to print."""
        ui.label("The run record").classes("font-semibold mt-4")
        ui.label("Every number of the plan, and of the analysis once one has been run, with its formula, the "
                 "formula with this run's numbers in it, what it means and what it assumes; the nuclear data and "
                 "its provenance; the method step by step; and what the simulation leaves out. The same page is "
                 "in the report's zip as record.html, to print to PDF.").classes("text-xs ps-muted")
        holder = ui.column().classes("w-full")

        def scene_png() -> Optional[str]:
            try:
                return base64.b64encode(paper.preview(P(), "geometry", state["journal"], state["width"])).decode()
            except Exception:  # noqa: BLE001 -- the record stands without the picture
                return None

        async def show() -> None:
            holder.clear()
            with holder:
                ui.spinner(size="md")
            page = await run.io_bound(lambda: P().record_html(scene_png=scene_png()))
            holder.clear()
            with holder:
                ui.add_css(record_css(".ps-record"))
                ui.html('<div class="ps-record">' + page.split("<body>", 1)[1].rsplit("</body>", 1)[0] + "</div>",
                        sanitize=False).classes("w-full ps-plate")
                ui.button("Save the record (HTML; print it to PDF)", icon="download",
                          on_click=lambda: ui.download.content(page.encode("utf-8"), "run-record.html")).props(
                    "dense flat no-caps")

        ui.button("Show the run record", icon="menu_book", on_click=show).props("dense flat no-caps")

    @ui.refreshable
    def efficiency_panel() -> None:
        """The full-energy and total efficiency of each γ-ray detector against energy."""
        gds = P().experiment.gamma_detectors
        if not gds:
            ui.label("No γ-ray detectors.").classes("text-sm ps-muted")
            return
        with ui.row().classes("w-full gap-2").style("flex-wrap: wrap"):
            for i, gd in enumerate(gds):
                with ui.column().classes("ps-plate gap-0").style("flex: 1 1 360px; min-width: 0"):
                    ui.label(gd.name or f"γ{i + 1}").classes("text-sm font-medium")
                    ui.plotly(themed(figure_efficiency(P(), f"gamma:{i}"), state["theme"])).classes("w-full")

    @ui.refreshable
    def plan_answers() -> None:
        tab_plan.render_answers(ctx)

    panels = {"geometry": geometry_panel, "kinematics": kinematics_panel, "rates": rates_panel,
              "energy_loss": energy_loss_panel, "spectra": spectra_panel, "trajectories": trajectories_panel,
              "gamma": gamma_panel, "report": report_panel, "efficiency": efficiency_panel,
              "plan": plan_answers}

    def refresh_results() -> None:
        readouts.refresh()
        reading_box.refresh()
        run_panel.refresh()
        data_view.refresh()
        analysis_view.refresh()
        report_view.refresh()
        for name, p in panels.items():
            if name == "geometry" and scene_current():
                continue  # the scene is redrawn in place, so the camera stays where it is
            p.refresh()

    @ui.refreshable
    def reading_box(tab: str) -> None:
        try:
            lines = guide.reading(P(), tab)
        except Exception:  # noqa: BLE001 -- a reading must never break the page
            lines = []
        if not lines:
            return
        with ui.row().classes("w-full items-start no-wrap gap-2 ps-note rounded p-2"):
            ui.icon("lightbulb").classes("ps-accent mt-0.5")
            with ui.column().classes("gap-1"):
                ui.label("How to read this").classes("text-xs font-semibold ps-accent uppercase")
                ui.label(" ".join(lines)).classes("text-sm")

    def explain(tab: str):
        e = P().explain(tab)
        with ui.expansion("Explain: " + e["title"], icon="school").classes("w-full ps-sunk mt-2"):
            ui.markdown(f"**Formula.** {e['formula']}\n\n**Assumptions.** {e['assumptions']}\n\n"
                        f"**Where it stops being valid.** {e['limits']}\n\n"
                        f"**Validation:** see the physics register page `{e['register']}` in physim's docs.")

    def show_physics(section: str) -> None:
        """Open the Physics tab at a section (:data:`physim.nuclear.workbench.PHYSICS`)."""
        if state.get("tabs") is not None:
            state["tabs"].set_value("physics")
        ui.run_javascript(f"setTimeout(() => document.querySelector('.physics-{section}')"
                          "?.scrollIntoView({behavior: 'smooth', block: 'start'}), 400)")

    def block(tab: str, title: Optional[str] = None) -> None:
        """One of the planner's results: its heading, how to read it, the panel, and its explanation."""
        if title:
            ui.label(title).classes("text-lg font-semibold mt-2")
        if tab == "geometry":  # the scene's full-screen mode hides what stands above and below it
            with ui.element("div").classes("w-full") as above:
                reading_box(tab)
            panels[tab]()
            with ui.element("div").classes("w-full") as below:
                explain(tab)
            scene_state["extras"] = [above, below]
            return
        reading_box(tab)
        panels[tab]()
        explain(tab)

    # -- the data tab ---------------------------------------------------------------------------------------------
    data_options = {"correction": "recoil" if P().experiment.excitation is not None
                    and P().experiment.excitation.excite == "target" else "projectile",
                    "mode": "coincidence", "randoms": "shown", "addback": True, "crystals": False, "gate": None,
                    "bins": 200, "lo_kev": None, "hi_kev": None, "scaled": False, "view_detector": None,
                    "view_correction": "off", "other": None, "compare_spectrum": None}
    gate_form = {"name": "", "detector": None, "rings": "", "group": "excited", "energy": ""}

    @ui.refreshable
    def data_view() -> None:
        try:
            tab_data.render_view(ctx)
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty tab
            ui.label(f"The data could not be shown: {err}").classes("ps-bad")

    @ui.refreshable
    def report_view() -> None:
        try:
            tab_report.render_view(ctx)
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty tab
            ui.label(f"The report could not be shown: {err}").classes("text-sm ps-warn")

    @ui.refreshable
    def analysis_view() -> None:
        try:
            tab_analysis.render_view(ctx)
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty tab
            ui.label(f"The analysis could not be shown: {err}").classes("text-sm ps-warn")

    # -- the run tab ----------------------------------------------------------------------------------------------
    run_state = {"watching": False}
    run_form = {"kind": "beam", "duration": None, "budget": "10 min", "source": "152Eu", "activity": "37 kBq",
                "position": "", "offset_mm": 2.0}

    @ui.refreshable
    def run_panel() -> None:
        tab_run.render_view(ctx)

    def load_run(n) -> None:
        if n is None or (P().run is not None and P().run.number == n):
            return
        P().load_run(n)
        status_strip.refresh()
        setup_panel.refresh()
        refresh_results()

    def start_run(f: dict) -> None:
        kind = f["kind"]
        options = {}
        if kind == "source":
            options = {"source": f["source"], "activity": f["activity"]}
            if str(f.get("position") or "").strip():
                try:
                    options["position"] = [float(x) for x in str(f["position"]).replace(";", ",").split(",")]
                except ValueError:
                    ui.notify("The position is three numbers in mm: x, y, z.", type="warning")
                    return
        elif kind == "alignment":
            options = {"offset_mm": float(f["offset_mm"] or 0.0)}
        try:
            out = P().start_run(kind, f["duration"] or None, f["budget"], background=kind != "source", **options)
        except (SetupError, ValueError, _runs.RunInProgress) as err:
            ui.notify(str(err), type="negative", multi_line=True)
            return
        after_start(out)

    def after_start(out) -> None:
        setup_panel.refresh()
        run_panel.refresh()
        if isinstance(out, _runs.RunData):
            finished()
            return
        run_state["watching"] = True  # the page's timer follows the run (see tick)

    def extend_run(more: str) -> None:
        try:
            out = P().extend_run(more, background=True)
        except (ValueError, _runs.RunInProgress) as err:
            ui.notify(str(err), type="negative", multi_line=True)
            return
        after_start(out)

    def tick() -> None:
        """Every half second: the live counters while a run is taken, and the views once it ends."""
        if not run_state["watching"]:
            return
        if P().running:
            run_panel.refresh()
            return
        run_state["watching"] = False
        bg = P()._background
        if bg is not None and bg.error is not None:
            ui.notify(f"The run failed: {bg.error}", type="negative", multi_line=True)
        finished()

    def finished() -> None:
        setup_panel.refresh()
        status_strip.refresh()
        refresh_results()
        if P().run is not None:
            ui.notify(f"Run {P().run.number} stored: {P().run.describe()}", type="positive", multi_line=True)

    def stop_run() -> None:
        bg = P()._background
        if bg is not None:
            bg.stop()
            run_panel.refresh()

    async def watch_events() -> None:
        """A sample of the current run's events as tracks in the scene, on the Setup tab."""
        if state.get("tabs") is not None:
            state["tabs"].set_value("setup")
        show = scene_state.get("show_tracks")
        if show is not None:
            await show()

    # -- experiments ----------------------------------------------------------------------------------------------
    @ui.refreshable
    def experiments_list() -> None:
        found = _runs.list_experiments()
        ui.label("Experiments").classes("text-lg font-semibold")
        ui.label(f"In {_runs.home()}").classes("text-xs ps-muted")
        if not found:
            ui.label("No experiment yet: make one with New experiment.").classes("text-sm")
        for x in found[:30]:
            with ui.row().classes("w-full items-center no-wrap gap-2"):
                with ui.column().classes("gap-0"):
                    ui.label(x["name"]).classes("text-sm font-medium")
                    ui.label(f"{x['runs']} run{'s' * (x['runs'] != 1)} · "
                             f"{time.strftime('%Y-%m-%d %H:%M', time.localtime(x['modified']))}").classes(
                        "text-xs ps-muted")
                ui.space()
                ui.button("Open", on_click=lambda path=x["path"]: open_experiment(path)).props("dense flat no-caps")
        ui.button("New experiment", icon="add", on_click=lambda: (experiments_dialog.close(), new_dialog.open())).props(
            "unelevated no-caps")

    def new_experiment_form() -> None:
        ui.label("New experiment").classes("text-lg font-semibold")
        ui.label("The beam, the target and what to measure; the detectors come after, from the Add menus, or "
                 "start from a template.").classes("text-xs ps-muted")
        with ui.grid(columns=2).classes("w-full gap-1"):
            beam = ui.input("Beam", value="16O", placeholder="16O").props("dense outlined")
            energy = ui.input("Energy", value="64 MeV", placeholder="64 MeV or 4 MeV/u").props("dense outlined")
            target = ui.input("Target", value="208Pb", placeholder="208Pb, Au, CD2").props("dense outlined")
            thickness = ui.input("Thickness", value="0.5 mg/cm2").props("dense outlined")
        measure = ui.select(dict(wb.MEASUREMENTS), value="coulex-target", label="What to measure").props(
            "dense outlined").classes("w-full")
        tmpl = ui.select({k: label for k, label, _ in wb.TEMPLATES}, value="blank", label="Detectors from").props(
            "dense outlined").classes("w-full")
        name = ui.input("Name (made from the beam and target if empty)").props("dense outlined").classes("w-full")
        ui.label("A template brings only its detectors; the beam, target and measurement are the ones above.").classes(
            "text-xs ps-muted")
        ui.button("Make the experiment", icon="science",
                  on_click=lambda: create_experiment(beam.value, energy.value, target.value, thickness.value,
                                                     measure.value, tmpl.value, name.value)).props(
            "unelevated no-caps")

    # -- the main area --------------------------------------------------------------------------------------------
    stage_modules = {"setup": tab_setup, "plan": tab_plan, "run": tab_run, "data": tab_data,
                     "analysis": tab_analysis, "report": tab_report, "physics": tab_physics}

    @ui.refreshable
    def status_strip() -> None:
        _strip.render(ctx)

    @ui.refreshable
    def main_area() -> None:
        with ui.tabs(value=state["stage"], on_change=lambda e: state.update(stage=e.value)).classes("w-full").props(
                "dense no-caps align=left") as tabs:
            for key, label in wb.STAGES:
                ui.tab(key, label=label)
        state["tabs"] = tabs
        with ui.tab_panels(tabs, value=state["stage"]).classes("w-full"):
            for key, _ in wb.STAGES:
                with ui.tab_panel(key):
                    stage_modules[key].render(ctx)

    def set_theme(name: str) -> None:
        if name == state["theme"]:
            return
        state["theme"] = name
        dark.value = name == "dark"
        ui.run_javascript(f"try {{ localStorage.setItem('physim-planner-theme', '{name}') }} catch (e) {{}}")
        refresh_results()  # the figures are drawn in the theme's colours

    async def restore_choices() -> None:
        """The theme this browser used last (the system's light or dark setting the first time), unless the
        address gives it (?theme=...)."""
        try:
            saved_theme, system_dark = await ui.run_javascript(
                "(() => { let t = null; try { t = localStorage.getItem('physim-planner-theme') } catch (e) {} "
                "return [t, window.matchMedia('(prefers-color-scheme: dark)').matches] })()", timeout=3)
        except Exception:  # noqa: BLE001 -- no answer from the browser: keep the defaults
            return
        if theme is None:
            wanted = saved_theme if saved_theme in THEMES else ("dark" if system_dark else "light")
            if wanted != state["theme"]:
                theme_toggle.value = wanted

    ctx = SimpleNamespace(
        ui=ui, P=P, state=state, edit=edit, section=section, help_icon=help_icon, bind=bind,
        reaction_section=reaction_section, clover_switches=clover_switches, add_detector=add_detector,
        duplicate_detector=duplicate_detector, remove_detector=remove_detector,
        add_gamma_detector=add_gamma_detector, duplicate_gamma_detector=duplicate_gamma_detector,
        remove_gamma_detector=remove_gamma_detector, open_card=open_card, block=block, readouts=readouts, analysis_block=analysis_block, alignment_block=alignment_block,
        multistep_block=multistep_block, open_explanation=open_explanation, efficiency_panel=efficiency_panel,
        plan_answers=plan_answers, run_view=run_panel, run_form=run_form, refresh_run_form=run_panel.refresh,
        start_run=start_run, stop_run=stop_run, extend_run=extend_run, load_run=load_run, watch_events=watch_events,
        columns=columns, data_view=data_view, refresh_data=data_view.refresh, data_options=data_options,
        gate_form=gate_form, fmt=_fmt, time_text=_time,
        analysis_view=analysis_view, report_view=report_view, show_physics=show_physics, levels_block=levels_block,
        geometry_tables=geometry_tables, explain=explain, reading_box=reading_box, gamma_panel=gamma_panel)

    # -- layout -----------------------------------------------------------------------------------------------
    ui.add_css(STYLE)
    ui.add_css(record_css(".ps-record"))
    dark = ui.dark_mode(state["theme"] == "dark")
    with ui.header(elevated=False).classes("items-center ps-header py-1"):
        ui.label("physim").classes("text-lg font-semibold")
        ui.label("experiment workbench").classes("ps-muted")
        ui.space()
        ui.button("Experiments", icon="folder_open",
                  on_click=lambda: (experiments_list.refresh(), experiments_dialog.open())).props(
            "flat no-caps")
        ui.button("New experiment", icon="add", on_click=lambda: new_dialog.open()).props("flat no-caps")
        ui.button("Import a setup file", icon="upload", on_click=lambda: upload_dialog.open()).props("flat no-caps")
        ui.button("Export the setup", icon="download", on_click=export_setup).props("flat no-caps")
        theme_toggle = ui.toggle({"light": "Light", "dark": "Dark"}, value=state["theme"],
                                 on_change=lambda e: set_theme(e.value)).props("dense no-caps unelevated")
    with ui.dialog() as upload_dialog, ui.card():
        ui.label("Import a setup file (.toml) as this experiment's setup")
        ui.upload(auto_upload=True, on_upload=import_file).props("accept=.toml max-files=1")
    with ui.dialog() as experiments_dialog, ui.card().classes("ps-card").style("min-width: min(560px, 95vw)"):
        experiments_list()
    with ui.dialog() as new_dialog, ui.card().classes("ps-card").style("min-width: min(560px, 95vw)"):
        new_experiment_form()
    with ui.dialog() as export_dialog, ui.card().classes("ps-card").style("max-width: min(1100px, 95vw)"):
        with ui.row().classes("items-start gap-6"):
            with ui.column().classes("gap-2").style("width: 320px; max-width: 100%"):
                export_controls()
            preview_box = ui.column().classes("items-center gap-2")
    export_dialog.on("show", render_preview)
    with ui.left_drawer(value=True).classes("ps-rail").props("width=440 bordered behavior=desktop") as drawer:
        setup_panel()
    state["drawer"] = drawer
    with ui.column().classes("w-full gap-2"):
        with ui.element("div").classes("w-full") as strip_box:
            status_strip()
        state["strip_box"] = strip_box
        main_area()
    if offer_list:
        experiments_dialog.open()
    ui.timer(0.2, restore_choices, once=True)
    ui.timer(0.5, tick)
