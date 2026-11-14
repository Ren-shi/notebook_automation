"""The planner page: the rail on the left (the experiment, the steps, the Physics reference and the theme), the
stages, and the setup panel on the right of the Setup step, assembled for one browser client."""

from __future__ import annotations

import asyncio
import base64
import math
import threading
import time
from typing import Optional


from .. import guide
from ..analysis import SHIFT_H
from ..record import record_css
from ..planner import Planner
from ..quantity import Quantity

from .. import dataviews as dv
from .figures import (CLOVER_FIELDS, REACTION_FIELDS, REACTION_TYPES, STYLE, THEMES, _shown, _value,
                      figure_alignment, figure_gamma_crystal, figure_peak, figure_search, figure_detector_spectrum, figure_efficiency, figure_energy_loss,
                      figure_excitation, figure_gamma_spectra, figure_kinematics, figure_levels, figure_overlay,
                      figure_source_spectrum, figure_spectra, figure_strips, figure_sweep, figure_trajectories,
                      themed)

#: How a change that is not applied yet looks: an outline on its field, a mark on its card, the bar that holds
#: Apply and Discard at the top of the setup panel, and the note above the tabs.
PENDING_STYLE = """
.ps-pending.q-field--outlined .q-field__control:before { border: 2px solid var(--ps-warn) !important; }
.ps-pending .q-field__label { color: var(--ps-warn) !important; }
.ps-pending-bar { padding: 0 2px; }
.ps-pending-bar.ps-pending-on { position: sticky; top: 0; z-index: 3; background: var(--ps-warnbox);
  border: 1px solid var(--ps-warn); border-radius: 8px; padding: 8px 10px; margin: 0 0 6px 0; }
.ps-pending-note { background: var(--ps-warnbox); border: 1px solid var(--ps-warn); margin-top: 6px; }
.ps-pending-mark { color: var(--ps-warn); border: 1px solid var(--ps-warn); border-radius: 8px; padding: 0 6px;
  font-size: 11px; line-height: 16px; white-space: nowrap; }
"""

#: The rail on the left (the experiment, then the stages as numbered steps, each with its state and a line of
#: status) and the foot of a step.
STEP_STYLE = """
.ps-nav-col { min-height: calc(100vh - 32px); padding: 4px 2px; }
.ps-wordmark { font-family: var(--ps-mono); font-size: 16px; font-weight: 600; padding-left: 6px; }
.ps-experiment { border: 1px solid var(--ps-line); background: var(--ps-surface); color: var(--ps-ink);
  min-height: 40px; }
.ps-experiment .q-btn__content { justify-content: space-between; flex-wrap: nowrap; text-align: left; }
.ps-steps { counter-reset: ps-step; width: 100%; }
.ps-steps .q-tab { counter-increment: ps-step; justify-content: flex-start; min-height: 54px; padding: 6px 8px;
  margin: 1px 0; border-radius: 6px; overflow: hidden; }
.ps-steps .q-tab__content { align-items: flex-start; width: 100%; min-width: 0; padding: 0 0 0 36px;
  position: relative; }
.ps-steps .q-tab__content::before { content: counter(ps-step); position: absolute; left: 0; top: 50%;
  transform: translateY(-50%); width: 26px; height: 26px; border-radius: 4px; background: var(--ps-raise);
  color: var(--ps-ink); font-family: var(--ps-mono); font-size: 12px; font-weight: 600; line-height: 26px;
  text-align: center; }
.ps-steps .ps-step-done .q-tab__content::before { content: "\\2713"; color: var(--ps-ok);
  background: color-mix(in srgb, var(--ps-ok) 18%, transparent); }
.ps-steps .q-tab--active { background: var(--ps-raise); }
.ps-steps .q-tab--active .q-tab__content::before { content: counter(ps-step); background: var(--ps-fill);
  color: var(--ps-on-fill); }
.ps-steps .q-tab__label { font-size: 14px; font-weight: 600; line-height: 1.35; }
.ps-steps .ps-step-waiting .q-tab__label { color: var(--ps-muted); }
.ps-steps .q-tab__indicator { display: none; }
.ps-step-status { font-family: var(--ps-mono); font-size: 11px; font-weight: 400; color: var(--ps-muted);
  width: 100%; text-align: left; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.ps-step-status.ps-warn { color: var(--ps-warn); }
.ps-step-hidden { display: none !important; }
.q-tab-panels { background: transparent !important; }
.q-tab-panel { padding: 0; }
.ps-foot { border-top: 1px solid var(--ps-line); margin-top: 18px; padding: 10px 2px 4px; }
"""

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
               template: Optional[str] = None, experiment: Optional[str] = None, prebuild: bool = False) -> None:
    """Build the planner page for the current client (call inside a NiceGUI page function).

    It opens ``experiment`` (a folder), or a new setup from ``template`` (:data:`physim.nuclear.workbench.TEMPLATES`)
    or ``example``, or else the most recent experiment; with no experiment yet it offers the experiment list and
    "New experiment" first. ``theme`` is "light" or "dark"; ``None`` takes the one this browser used last, the
    system's setting the first time."""
    from types import SimpleNamespace

    from nicegui import Client, background_tasks, run, ui

    from .. import paper
    from .. import runs as _runs
    from .. import workbench as wb
    from ..experiment import Experiment, SetupError
    from . import panel as _panel
    from . import tab_analysis, tab_data, tab_physics, tab_plan, tab_report, tab_run, tab_setup

    first, offer_list = _starting_planner(example, template, experiment)
    state = {"planner": first, "events": events, "theme": theme or "light", "open": None, "stage": "setup",
             "export": None, "journal": "physical_review", "width": "single", "format": "pdf",
             # Whether the setup's numbers are there (the panel and the step bar show them); a page built whole
             # (``prebuild``) computes them as it is built.
             "ready": prebuild,
             # Whether a change of the setup is applied as it is made, or waits for Apply (the default).
             "auto_apply": False}
    client = ui.context.client

    def P() -> Planner:  # noqa: N802
        """The planner of the setup that is applied: every number on the page is its."""
        return state["planner"]

    def draft_of(planner: Planner) -> Planner:
        """A planner of its own for the setup as it is being edited. Edits go to it and compute nothing; Apply
        makes the page's planner adopt it, Discard makes it adopt the page's planner."""
        draft = Planner(planner.experiment)
        draft.adopt(planner)
        return draft

    state["draft"] = draft_of(first)

    def D() -> Planner:  # noqa: N802
        """The planner of the draft: the setup with the changes not applied yet."""
        return state["draft"]

    # -- actions ----------------------------------------------------------------------------------------------
    def changed(ok: bool, panel: bool = False) -> None:
        """After an edit of the setup. ``panel``: the cards themselves changed (a detector added, removed or
        moved), so the setup panel is rebuilt even under a field being typed in."""
        if ok:
            state["panel_due"] = panel or state.get("panel_due", False)
            refresh_results()  # rebuilds the setup panel and the step bar once the numbers are ready
            return
        ui.notify("; ".join(P().problems), type="negative", multi_line=True)
        setup_panel.refresh()
        sync_steps()

    def pending() -> list:
        """The changes not applied, one line each (:func:`physim.nuclear.workbench.changes`)."""
        return wb.changes(P()._draft, D()._draft)

    def value_in(setup: dict, section: str, field: str):
        """A field of a setup (in the structure of a setup file) by the planner's section name; None if the setup
        has no such field (or no such detector)."""
        if section == "title":
            return setup.get("title")
        if section == "backing":
            return setup["target"].get("backing", {}).get(field)
        if section.startswith(("gamma detector", "detector")):
            gamma = section.startswith("gamma")
            items = setup.get("gamma_detectors" if gamma else "detectors", [])
            i = int(section.split()[2 if gamma else 1]) - 1
            return items[i].get(field) if i < len(items) else None
        return setup.get(section, {}).get(field)

    def is_pending(section: str, field: str) -> bool:
        """Whether a field of the draft differs from the applied setup (a detector is found there by its name)."""
        now = value_in(D()._draft, section, field)
        if section.startswith(("gamma detector", "detector")):
            key = card_of(section)
            mine = wb.card_setup(D()._draft, key) or {}
            theirs = wb.card_setup(P()._draft, key, name=mine.get("name"))
            return theirs is None or theirs.get(field) != now
        return value_in(P()._draft, section, field) != now

    #: The setup panel's parts as they stand on the page, so an edit can update them where they are instead of
    #: rebuilding the panel under the next field being typed in: each card's labels by its key, each field's
    #: element by (section, field), and the text inputs among them.
    card_parts: dict = {}
    field_marks: list = []
    field_inputs: dict = {}
    field_shown: dict = {}

    def scene_shows() -> tuple:
        """What a drawing of the scene depends on: the draft's setup, and the applied one's numbers once they are
        there (the rings that see unsafe collisions are painted from them)."""
        return D().experiment, P().experiment if numbers_kept() else None

    def draw_scene() -> None:
        """Draw the scene again if what it shows has changed (a drawing takes the tracks off)."""
        view = scene_state["view"]
        if view is None or not view.alive or view.busy:
            return
        now = scene_shows()
        was = scene_state.get("shows")
        if was is None or any(a is not b for a, b in zip(was, now)):
            view.draw()
            scene_state["shows"] = now

    def mark(element, section: str, field: str, register: bool = True) -> None:
        """Outline a field's element if its value is not applied yet."""
        if register:
            field_marks.append((section, field, element))
        if is_pending(section, field):
            element.classes(add="ps-pending")
        else:
            element.classes(remove="ps-pending")

    def cards() -> list:
        """The setup panel's cards for the draft: (card, whether it holds a change not applied). A card as it is
        applied has its consequence line from the applied setup's numbers; a changed one says to press Apply; and
        while the numbers are on their way the lines say so."""
        d, p = D(), P()
        applied = p._cache.get(("page", "decisions", repr(p._draft), tuple(p.problems)))
        if applied is None and (prebuild or numbers_kept()):
            applied = decisions(p)
        if applied is not None and d._draft == p._draft and not d.problems:
            return [(c, False) for c in applied]
        by_title = {(c.group, c.title): c for c in applied or []}
        out = []
        for c in wb.draft_decisions(d):
            name = c.values.get("name") if c.key.startswith(("detector:", "gamma:")) else None
            same = wb.card_setup(d._draft, c.key) == wb.card_setup(p._draft, c.key, name=name)
            was = by_title.get((c.group, c.title))
            if not same:
                c.consequence = "Not applied yet: Apply works out its numbers"
            elif was is not None:
                c.consequence = was.consequence
            else:
                c.consequence = "…"
            out.append((c, not same))
        return out

    def checks() -> list:
        """What stands at the top of the setup panel: the draft's problems, then the applied setup's warnings."""
        found = [{"level": "error", "text": text} for text in D().problems]
        return found + (wb.checks(P()) if prebuild or numbers_kept() else [])

    def sync_cards() -> None:
        """Bring the setup panel up to date with the draft where it stands: the cards' lines, the fields' values
        and their marks. The panel is rebuilt only if its cards are no longer the draft's."""
        now = cards()
        if [c.key for c, _ in now] != list(card_parts) or any(part.title.is_deleted for part in card_parts.values()):
            setup_panel.refresh()
            return
        for c, changed_ in now:
            part = card_parts[c.key]
            part.title.set_text(c.title)
            part.state.set_text(c.state)
            part.consequence.set_text(c.consequence)
            part.mark.set_visibility(changed_)
        draft = D()._draft
        for (section, field), box in list(field_inputs.items()):
            if box.is_deleted:
                field_inputs.pop((section, field))
                continue
            # A field follows the draft when the draft's value changed, and is left alone otherwise: it may hold
            # what is being typed into it.
            text = _shown(field, value_in(draft, section, field)) or ""
            if text != field_shown.get((section, field)):
                field_shown[(section, field)] = text
                if (box.value or "") != text:
                    box.value = text
        field_marks[:] = [m for m in field_marks if not m[2].is_deleted]
        for section, field, element in field_marks:
            mark(element, section, field, register=False)
        checks_box.refresh()
        checks_line.refresh()

    def draft_views() -> None:
        """What follows the draft without any number being worked out: the scene, the quick numbers beside it and
        the geometry tables."""
        draw_scene()
        if "setup" in built and "setup" not in dirty:
            selection_panel.refresh()
            if detectors_of(D()._draft) != scene_state.get("tables"):
                geometry_tables.refresh()

    def detectors_of(setup: dict) -> str:
        return repr((setup.get("detectors"), setup.get("gamma_detectors")))

    def locked() -> bool:
        if P().running:
            ui.notify("A run is being taken: the setup is read-only until it ends.", type="warning")
        return P().running

    def staged(ok: bool, rebuild: bool = False) -> None:
        """After an edit of the draft. Nothing is computed: the cards, the bar and the scene follow the draft, and
        the numbers wait for Apply (unless every change is applied as it is made). ``rebuild``: the cards or their
        fields changed, so the panel is built again."""
        if not ok:
            ui.notify("; ".join(D().problems), type="negative", multi_line=True)
        elif state["auto_apply"] and not P().running:
            state["panel_due"] = rebuild or state.get("panel_due", False)
            apply_changes()
            if not rebuild:
                sync_cards()
            return
        if rebuild:
            setup_panel.refresh()
        else:
            sync_cards()
        pending_bar.refresh()
        pending_note.refresh()
        sync_foot()
        draft_views()
        if ok:
            compute_ahead()

    def compute_ahead() -> None:
        """Once the draft has stood for a second: its rates, warnings and Doppler table on the worker, on the
        draft's own planner. Apply adopts them, so its numbers are usually there at once. A later edit makes an
        earlier request stand down before it starts."""
        state["draft_edit"] = n = state.get("draft_edit", 0) + 1
        d = D()

        def numbers() -> None:
            with work:
                if state.get("draft_edit") != n or D() is not d or d.problems:
                    return
                snapshot = d.snapshot()
                try:
                    snapshot.rates()
                    snapshot.warnings()
                    if snapshot.experiment.excitation is not None:
                        snapshot.gamma()
                except Exception:  # noqa: BLE001 -- Apply says what is wrong
                    pass

        async def go() -> None:
            await asyncio.sleep(1.0)
            if gone() or state.get("draft_edit") != n or not pending():
                return
            await run.io_bound(numbers)

        later(go)

    def apply_changes() -> None:
        """Apply the draft: the page's planner takes it (with the numbers already worked out for it, if any), the
        experiment folder's setup is saved, and everything on the page is brought up to date. The setup panel is
        updated where it stands if a field is being typed in (Ctrl+Enter from a field keeps the cursor there)."""
        d, p = D(), P()
        if d._draft == p._draft:
            return
        if p.running:
            ui.notify("A run is being taken: the changes wait until it ends.", type="warning")
            return
        if d.problems:
            ui.notify("Not applied: " + "; ".join(d.problems), type="negative", multi_line=True)
            return
        state["draft_edit"] = state.get("draft_edit", 0) + 1
        p.adopt(d)
        sync_cards()
        pending_bar.refresh()
        pending_note.refresh()
        sync_foot()
        refresh_results()

    def discard_changes() -> None:
        """Throw the draft away: the fields, the cards and the scene go back to the applied setup."""
        state["draft_edit"] = state.get("draft_edit", 0) + 1
        D().adopt(P())
        state["compound"] = False
        setup_panel.refresh()
        pending_bar.refresh()
        pending_note.refresh()
        sync_foot()
        draft_views()

    def direct(action) -> None:
        """A change made with a button of its own on a result (a level scheme looked up, a matrix element, the
        emission pattern): ``action(planner)`` is applied at once, and made in the draft as well, where the
        changes not applied stay as they are. An error of the action is the caller's to report."""
        if locked():
            return
        ok = action(P())
        try:
            action(D())
        except Exception:  # noqa: BLE001 -- the draft may have moved away from where the action makes sense
            pass
        if D()._draft == P()._draft:
            D().adopt(P())  # nothing pending: the draft has the applied setup's results again
        changed(ok, panel=True)
        pending_bar.refresh()
        pending_note.refresh()
        sync_foot()

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
        go_to("setup")  # the setup is edited from the Setup step
        setup_dialog.open()
        if numbers_kept():
            setup_panel.refresh()
        else:  # the numbers of a change are on their way: the panel is rebuilt when they are there
            state["panel_due"] = True

    def edit(section: str, field: str, text, rebuild: bool = False) -> None:
        """A field of the setup was changed: the draft takes the value. ``rebuild``: the choice changes which
        fields or cards there are (a select does; a text field is updated where it stands)."""
        value = _value(field, text)
        old = value_in(D()._draft, section, field)
        if (None if old is None else str(old)) == (None if value is None else str(value)):
            return
        if locked():
            return
        state["open"] = card_of(section) or state["open"]
        staged(D().set(section, field, value), rebuild or field in ("model", "shape", "crystals", "nuclide", "material"))

    def add_detector(key: str = "side"):
        def run_() -> None:
            names = [det.get("name") for det in D().draft["detectors"]]
            if not locked():
                staged(D().add_detector(**guide.placement(key, names)), rebuild=True)
        return run_

    def duplicate_detector(k: int):
        def run() -> None:
            if not locked():
                staged(D().duplicate_detector(k), rebuild=True)
        return run

    def remove_detector(k: int):
        def run() -> None:
            if len(D().draft["detectors"]) == 1:
                ui.notify("A setup needs at least one particle detector.", type="warning")
                return
            if not locked():
                staged(D().remove_detector(k), rebuild=True)
        return run

    def add_gamma_detector(key: str = "disc"):
        def run() -> None:
            names = [g.get("name") for g in D().draft.get("gamma_detectors", [])]
            if not locked():
                staged(D().add_gamma_detector(**guide.gamma_placement(key, names)), rebuild=True)
        return run

    def duplicate_gamma_detector(k: int):
        def run() -> None:
            if not locked():
                staged(D().duplicate_gamma_detector(k), rebuild=True)
        return run

    def remove_gamma_detector(k: int):
        def run() -> None:
            if not locked():
                staged(D().remove_gamma_detector(k), rebuild=True)
        return run

    def use(planner: Planner) -> None:
        lost = len(pending())
        state["planner"] = planner
        if lost:
            ui.notify(f"{lost} change{'s' * (lost != 1)} of the setup you left had not been applied, and "
                      f"{'were' if lost != 1 else 'was'} dropped.", type="warning", multi_line=True,
                      close_button=True)
        state["draft"] = draft_of(planner)
        state["draft_edit"] = state.get("draft_edit", 0) + 1
        state["open"] = None
        state["compound"] = False
        state["ready"] = prebuild
        setup_panel.refresh()
        pending_bar.refresh()
        pending_note.refresh()
        sync_foot()
        sync_steps()
        main_area.refresh()
        refresh_results()

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
        if pending():
            ui.notify("The file holds the setup as it is applied, without the changes not applied yet.",
                      type="warning")
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
                field_inputs[(sec_name, field)] = inp
                field_shown[(sec_name, field)] = _shown(field, v) or ""
                mark(inp, sec_name, field)

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
            edit(sec_name, field, (new if target else f"{a}{new}") if a else f"{next(iter(wb.mass_numbers(new)))}{new}",
                 rebuild=True)

        el = ui.select(options, value=symbol if symbol in options else None, with_input=True,
                       label="Target element" if target else "Beam element",
                       on_change=lambda e: set_element(e.value)).props("dense outlined options-dense").classes("w-full")
        help_icon(el, sec_name, field)
        mark(el, sec_name, field)
        if compound:
            inp = ui.input("Material", value=_shown(field, value), placeholder=placeholder).props(
                "dense outlined").classes("w-full")
            bind(inp, sec_name, field)
            mark(inp, sec_name, field)
            return
        if symbol is None:
            ui.label("")
            return
        masses = {str(a): label for a, label in wb.mass_numbers(symbol).items()}
        if target:
            masses = {wb.NATURAL: "natural mix", **masses}
        current = str(parts[1]) if parts and parts[1] is not None else (wb.NATURAL if target else None)
        mass = ui.select(masses, value=current if current in masses else None, label="Mass number",
                         on_change=lambda e: edit(sec_name, field,
                                                  symbol if e.value == wb.NATURAL else f"{e.value}{symbol}",
                                                  rebuild=True)).props("dense outlined options-dense").classes("w-full")
        mark(mass, sec_name, field)

    @ui.refreshable
    def setup_panel():
        card_parts.clear()
        field_marks.clear()
        field_inputs.clear()
        field_shown.clear()
        if not state["ready"]:  # the rates are being computed on the worker (refresh_results)
            with ui.row().classes("items-center gap-2 ps-muted p-2"):
                ui.spinner(size="md")
                ui.label("Computing the setup's numbers…").classes("text-sm")
            return
        _panel.render(ctx)

    @ui.refreshable
    def checks_box() -> None:
        _panel.render_checks(ctx)

    def edit_setup() -> None:
        """Open the window the setup is edited in."""
        setup_dialog.open()

    def done_editing() -> None:
        """The window's main button: apply what was changed, if anything, and close. It stays open when the
        changes cannot be applied (a value is wrong, or a run is being taken), and says why."""
        d, p = D(), P()
        if d._draft != p._draft:
            if p.running or d.problems:
                apply_changes()  # says why not
                return
            apply_changes()
        setup_dialog.close()

    def cancel_editing() -> None:
        """Close the window without the changes made in it."""
        if pending():
            discard_changes()
        setup_dialog.close()

    @ui.refreshable
    def checks_line() -> None:
        """On the Setup step: errors and warnings in full, since they need something done; the notes as a count
        that opens the window, where all the checks are."""
        found = checks()
        loud = [c for c in found if c["level"] in ("error", "warning")]
        notes = len(found) - len(loud)
        with ui.row().classes("items-center gap-2"):
            if not loud and (prebuild or numbers_kept()):
                ui.icon("check_circle").classes("ps-ok")
                ui.label("No warnings").classes("text-sm ps-ok")
            if notes:
                ui.button(f"{notes} note{'s' * (notes != 1)} on the setup", color=None, on_click=edit_setup).props(
                    "dense flat no-caps").classes("ps-muted")
        if loud:
            with ui.card().classes("w-full ps-warnbox p-2 gap-1"):
                for c in loud:
                    with ui.row().classes("items-start no-wrap gap-2"):
                        ui.icon(_panel.LEVEL_ICON[c["level"]]).classes(_panel.LEVEL_STYLE[c["level"]])
                        ui.label(c["text"]).classes(_panel.LEVEL_STYLE[c["level"]] + " text-sm")

    def changes_and_buttons(found: list, lines: int) -> None:
        for text in found[:lines]:
            ui.label(text).classes("text-xs ps-num")
        if len(found) > lines:
            ui.label(f"… and {len(found) - lines} more").classes("text-xs ps-muted")
        with ui.row().classes("items-center gap-2"):
            go = ui.button("Apply", icon="check", on_click=lambda: apply_changes()).props(
                "dense unelevated no-caps").tooltip("Work out the numbers of the setup as you have changed it "
                                                    "(Ctrl+Enter)")
            if D().problems or P().running:
                go.disable()
            ui.button("Discard", icon="undo", on_click=discard_changes).props("dense flat no-caps").tooltip(
                "Back to the setup as it is applied")

    @ui.refreshable
    def pending_bar() -> None:
        """At the top of the setup panel: the changes not applied, with Apply and Discard, and the switch that
        applies every change as it is made."""
        found = pending()
        with ui.column().classes("w-full gap-1 ps-pending-bar" + (" ps-pending-on" if found else "")):
            if found:
                with ui.row().classes("items-center gap-1 no-wrap"):
                    ui.icon("edit_note").classes("ps-warn")
                    ui.label(f"{len(found)} change{'s' * (len(found) != 1)} not applied").classes(
                        "text-sm font-semibold ps-warn")
                changes_and_buttons(found, 5)
            ui.switch("Apply each change as I make it", value=state["auto_apply"],
                      on_change=lambda e: set_auto_apply(e.value)).props("dense size=xs").classes(
                "text-xs ps-muted").tooltip("Off: change several fields, then press Apply. On: the numbers are "
                                            "worked out after every change, as before.")

    @ui.refreshable
    def pending_note() -> None:
        """Above the tabs while changes wait: whose numbers the page shows."""
        found = pending()
        if not found:
            return
        with ui.row().classes("w-full items-center gap-3 ps-pending-note rounded px-3 py-1"):
            ui.icon("edit_note").classes("ps-warn")
            ui.label(f"The numbers on this page are for the setup before your {len(found)} "
                     f"change{'s' * (len(found) != 1)}. The scene shows the changed setup.").classes("text-sm")
            ui.space()
            ui.button("Apply", icon="check", on_click=lambda: apply_changes()).props("dense unelevated no-caps")
            ui.button("Discard", icon="undo", on_click=discard_changes).props("dense flat no-caps")

    def set_auto_apply(on: bool) -> None:
        state["auto_apply"] = bool(on)
        ui.run_javascript("try { localStorage.setItem('physim-planner-auto-apply', '%s') } catch (e) {}"
                          % ("on" if on else "off"))
        if on:
            apply_changes()

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
                        on_change=lambda e: edit("reaction", "type", e.value, rebuild=True)).props(
            "dense outlined").classes("w-full")
        help_icon(sel, "reaction", "type")
        mark(sel, "reaction", "type")
        if kind != "coulex":
            return
        with ui.grid(columns=2).classes("w-full gap-1"):
            exc = ui.select({"target": "Target nucleus", "projectile": "Beam nucleus"},
                            value=reaction.get("excite", "target"), label="Excited nucleus",
                            on_change=lambda e: edit("reaction", "excite", e.value, rebuild=True)).props(
                "dense outlined")
            help_icon(exc, "reaction", "excite")
            mark(exc, "reaction", "excite")
            mul = ui.select(["E1", "E2", "E3"], value=reaction.get("multipolarity", "E2"), label="Multipolarity",
                            on_change=lambda e: edit("reaction", "multipolarity", e.value, rebuild=True)).props(
                "dense outlined")
            help_icon(mul, "reaction", "multipolarity")
            mark(mul, "reaction", "multipolarity")
        state_picker(reaction)
        section("", REACTION_FIELDS, "reaction", reaction)

    def state_picker(reaction: dict) -> None:
        """The excited state, chosen from the excited nucleus's level scheme (read from ENSDF); choosing one fills in
        its energy and B(Eλ↑) below, which can still be typed over."""
        role = "target" if reaction.get("excite", "target") == "target" else "beam"
        multipolarity = reaction.get("multipolarity", "E2")
        scheme = D().experiment.levels.get(role) if not D().problems else None
        choices = wb.states(scheme, multipolarity)
        nuclide = D().level_nuclides()[role] or role
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
            if level is None or level == current or locked():
                return
            try:
                staged(D().use_state(role, level, multipolarity), rebuild=True)
            except ValueError as err:
                ui.notify(str(err), type="warning")

        ui.select(choices, value=current, label=f"State of {nuclide}",
                  on_change=lambda e: pick(e.value)).props("dense outlined options-dense").classes("w-full")

    # -- results kept for the setup on screen, computed off the event loop ------------------------------------------
    #: One computation at a time for this page: a second request for the same numbers waits for the first, then
    #: finds them kept. Loading a run's events has its own turn, so the setup's numbers do not wait behind it.
    work = threading.Lock()
    run_work = threading.Lock()

    def numbers_kept() -> bool:
        """Whether the current setup's rates are there. While the worker computes them, nothing on the event loop
        asks for them (it would compute them a second time, with the page waiting)."""
        return "rates" in P()._cache

    # Whatever the page computes on a worker thread, a button's computation included, it computes on a snapshot
    # of the planner (``P().snapshot()``): a setup applied meanwhile gets none of its results, and what the
    # planner writes beside a run or the experiment (analysis.json, alignment.json, multistep.json) names the
    # setup it was computed for.
    def keep(p: Planner, key: tuple, make, of_run: bool = False):
        """A figure or a table of the page, made once and kept: with the planner's own results for the setup
        ``p`` holds, so it goes when the setup changes; or (``of_run``) with the current run, for what is read
        from the run's events alone, so a change of the setup leaves it. A tab shown again, a second tab with the
        same figure and the worker that prepares a tab all find it here."""
        if of_run and p.run is not None:
            kept = p.run.__dict__.setdefault("_page_kept", {})
            if key not in kept:
                kept[key] = make()
            return kept[key]
        cache = p._cache
        if ("page",) + key not in cache:
            cache[("page",) + key] = make()
        return cache[("page",) + key]

    def sendable(fig) -> dict:
        """A Plotly figure as the browser takes it. ``ui.plotly`` turns a figure into this every time the element
        is built or styled (5 ms a time, three times per element); a figure kept in this form is turned once."""
        return fig.to_plotly_json()

    def in_theme(fig):
        """The figure in the page's theme, ready to send (:func:`themed` takes 20 to 60 ms: the kept figures are
        themed when they are made, and come through here as they are)."""
        return fig if isinstance(fig, dict) else sendable(themed(fig, state["theme"]))

    #: The figures the panels show, by name: each takes the planner and what the panel chose.
    makers = {"kinematics": figure_kinematics, "energy_loss": figure_energy_loss,
              "trajectories": figure_trajectories, "gamma": figure_excitation, "strips": figure_strips,
              "spectra": figure_spectra, "efficiency": figure_efficiency,
              "levels": lambda p, role: figure_levels(p.levels()[role])}

    def shown(name: str, *chosen, p: Optional[Planner] = None):
        """The figure ``name`` of the setup, in the page's theme: drawn once per setup and theme."""
        p = p or P()
        theme = state["theme"]
        return keep(p, ("figure", theme, name) + chosen, lambda: sendable(themed(makers[name](p, *chosen), theme)))

    def decisions(p: Optional[Planner] = None) -> list:
        """The setup panel's cards (:func:`physim.nuclear.workbench.decisions`), worked out once per setup and
        per state of the fields: each card's consequence line reads the rates, the stopping tables or the
        kinematics."""
        p = p or P()
        return keep(p, ("decisions", repr(p._draft), tuple(p.problems)), lambda: wb.decisions(p))

    def reading(tab: str, p: Optional[Planner] = None) -> list:
        """The "How to read this" lines of a panel (:func:`physim.nuclear.guide.reading`), once per setup."""
        p = p or P()
        return keep(p, ("reading", tab), lambda: guide.reading(p, tab))

    def needs(tab: str, p: Planner) -> None:
        """Compute and keep what the panel ``tab`` shows (on a worker thread, before the panel is built)."""
        if tab in guide._READINGS:
            try:
                reading(tab, p)
            except Exception:  # noqa: BLE001 -- a reading must never break the page
                pass
        if tab == "geometry":
            p.geometry()
            p.selection()
        elif tab in ("kinematics", "energy_loss", "trajectories"):
            shown(tab, p=p)
        elif tab == "gamma":
            for role in ("target", "beam"):
                if p.levels()[role]:
                    shown("levels", role, p=p)
            if p.gamma()["available"]:
                shown("gamma", p=p)
                p.populations()
        elif tab == "efficiency":
            for i in range(len(p.experiment.gamma_detectors)):
                shown("efficiency", f"gamma:{i}", p=p)
        elif tab == "rates":
            shown("strips", p.rates()["rows"][0]["detector"], p=p)
        elif tab == "spectra":
            shown("spectra", state["events"], 1, p=p)

    #: The panels each stage shows at once, and those it shows when asked (the Plan tab's checks, which open).
    stage_needs = {"setup": (("geometry",), ()),
                   "plan": ((), ("kinematics", "energy_loss", "trajectories", "gamma", "efficiency", "rates",
                                 "spectra")),
                   "physics": (("gamma", "kinematics", "energy_loss", "trajectories", "geometry", "efficiency"), ()),
                   "data": ((), ()), "run": ((), ()), "analysis": ((), ()), "report": ((), ())}

    def warm(p: Planner) -> None:
        """The numbers every tab starts from: the rates, the warnings, the Doppler table."""
        p.rates()
        p.warnings()
        if p.experiment.excitation is not None:
            p.gamma()
        decisions(p)

    def prepare(key: str, p: Planner, everything: bool = False) -> None:
        """Compute and keep what the stage ``key`` shows when it is built (``everything``: also what it shows only
        when asked), for the snapshot ``p``. Runs on a worker thread; the page builds the stage afterwards from
        what is kept, which takes the event loop a few tens of milliseconds instead of seconds."""
        with work:
            try:
                warm(p)
            except Exception:  # noqa: BLE001 -- the panels say what is wrong
                return
            first, later = stage_needs[key]
            if key == "plan":  # the checks left open are built with the tab
                first = tuple(tab for tab in later if tab in state.get("checks_open", ()))
            for tab in first + (later if everything else ()):
                try:
                    needs(tab, p)
                except Exception:  # noqa: BLE001 -- the panel says what is wrong when it is built
                    pass
        if key == "data" and p.run is not None:
            with run_work:
                try:
                    tab_data.prepare(ctx, p)
                except Exception:  # noqa: BLE001 -- the tab says what is wrong when it is built
                    pass

    def with_work(fn, *args) -> None:
        with work:
            try:
                fn(*args)
            except Exception:  # noqa: BLE001 -- the panel says what is wrong when it is built
                pass

    def when_opened(expansion, build, tab: Optional[str] = None) -> None:
        """Build what an expansion holds when it is first opened, not before (most stay closed): ``build()`` fills
        it, after the worker has computed what the panel ``tab`` shows (a spinner meanwhile; no wait if it is
        kept already). An expansion that starts open, or a page built whole, is filled at once."""
        if prebuild or expansion.value:
            with expansion:
                build()
            return
        with expansion:
            box = ui.column().classes("w-full gap-2")
        filled = []

        async def fill(e) -> None:
            if not e.value or filled:
                return
            filled.append(True)
            if tab is not None:
                with box:
                    ui.spinner(size="md")
                await run.io_bound(with_work, needs, tab, P().snapshot())
                if box.is_deleted:
                    return
                box.clear()
            with box:
                build()

        expansion.on_value_change(fill)

    def later(job) -> None:
        """Start a coroutine of the page on its own, not inside the element whose event asked for it: that element
        may be rebuilt before the coroutine ends, and NiceGUI cancels what a deleted element started."""
        async def run_() -> None:
            with client:
                await job()

        background_tasks.create(run_(), name="planner page")

    def gone() -> bool:
        """Whether this page's browser tab is closed (nothing is left to update)."""
        return client.id not in Client.instances

    # -- figures and their export ---------------------------------------------------------------------------------
    def plot(fig, name: str, **options) -> None:
        """A figure on its plate, with the button that exports it in a journal's style."""
        with ui.column().classes("w-full ps-plate gap-0"):
            with ui.row().classes("w-full items-center"):
                ui.space()
                ui.button("Paper figure", icon="article", on_click=lambda: open_export(name, options)).props(
                    "dense flat no-caps").tooltip("Export this figure in a journal's style")
            ui.plotly(in_theme(fig)).classes("w-full")

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
            png = await run.io_bound(lambda p=P().snapshot(): paper.preview(p, name, j, w, **options))
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
            data = await run.io_bound(lambda p=P().snapshot(): paper.export(p, name, fmt, j, w, **options))
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

    @ui.refreshable
    def setup_strip() -> None:
        """Above the scene: the whole setup in a few numbers (those of the selection panel with nothing
        selected)."""
        if not numbers_kept():
            return
        try:
            s = P().selection()
        except Exception:  # noqa: BLE001 -- the strip must never break the page
            return
        cells = [("Particles counted", _fmt(s["rate_per_s"]), "/s", "")]
        if s["measured"] != "all":
            cells.append((f"Of these, {s['measured']}", _fmt(s["measured_rate_per_s"]), "/s", "ps-ok"))
        cells.append(("Solid angle, particle detectors", _fmt(s["solid_angle_msr"]), "msr", ""))
        if s["gamma_efficiency"] is not None:
            cells.append(("γ-ray efficiency, all detectors", f"{100 * s['gamma_efficiency']:.2f}", "%", ""))
        with ui.element("div").classes("ps-readouts"):
            for title, value, unit, colour in cells:
                with ui.element("div").classes("ps-readout"):
                    ui.label(title).classes("text-xs ps-muted")
                    with ui.row().classes("items-baseline gap-1 no-wrap"):
                        ui.label(value).classes("ps-readout-value " + colour)
                        ui.label(unit).classes("text-xs ps-muted ps-num")

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
        """A detector was dropped at a new place: the draft takes it. The scene and the quick numbers beside it
        follow at once; the other numbers wait for Apply."""
        if locked():
            return
        staged(D().move(key, position, mode))

    class SceneSetup:
        """What the scene asks of a planner, answered from the draft: the setup to draw, the quick numbers of a
        drag, and the rings that see unsafe collisions once they are known (they need the rates)."""

        @property
        def experiment(self):
            return D().experiment

        def live(self, key: str, fields=None) -> dict:
            return D().live(key, fields)

        def safety(self) -> dict:
            if D()._draft == P()._draft:
                return P().safety() if numbers_kept() else {}
            return D().safety() if "rates" in D()._cache else {}

    scene_setup = SceneSetup()

    def scene_current() -> bool:
        """Redraw the scene in place (the camera stays), if it is on the page and in the right theme."""
        view = scene_state["view"]
        if view is None or not view.alive or view.theme != state["theme"]:
            return False
        draw_scene()
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
            entry = D().draft["detectors" if kind == "detector" else "gamma_detectors"][int(i)]
            if str(entry.get(field)) == text or locked():
                return
            staged(D().place(key, **{field: text}))
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
        if D()._draft != P()._draft:  # changes wait for Apply: the quick numbers of the draft
            draft_selection(key)
            return
        if not numbers_kept():  # a change is being computed: the panel is rebuilt when its numbers are there
            with ui.row().classes("items-center gap-2 ps-muted ps-card p-3 w-full"):
                ui.spinner(size="sm")
                ui.label("Computing…").classes("text-sm")
            scene_state["live"] = {}
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
            exact_values(key, fields)

        with ui.column().classes("w-full ps-card p-3 gap-1"):
            if s["kind"] == "experiment":
                ui.label("The whole experiment").classes("ps-section")
                ui.label(s["title"]).classes("text-base font-medium")
                row("Particle detectors", str(s["detectors"]))
                row("γ-ray detectors", str(s["gamma_detectors"]))
                # (The rates, the solid angle and the γ-ray efficiency are in the strip above the scene.)
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
                    fig = await run.io_bound(figure_detector_spectrum, P().snapshot(), s["name"], state["events"])
                    spectrum_box.clear()
                    with spectrum_box:
                        ui.plotly(themed(fig, state["theme"])).classes("w-full")

                ui.button("Simulate its spectrum", icon="play_arrow", on_click=spectrum).props(
                    "dense flat no-caps").classes("mt-1")

    def exact_values(key: str, fields) -> None:
        """The place of the selected detector as fields to type in (they change the draft, as a drag does)."""
        kind, i = key.split(":")
        entry = D().draft["detectors" if kind == "detector" else "gamma_detectors"][int(i)]
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

    def draft_selection(key: Optional[str]) -> None:
        """Beside the scene while changes wait for Apply: what follows from the draft's geometry alone. For a
        detector, where it stands, what it covers and a coarse estimate of its rate (as during a drag)."""
        from .. import scene as _scene

        live = scene_state["live"] = {}
        d = D()
        values = None
        if key is not None and key != "target":
            try:
                values = d.live(key)
            except (IndexError, KeyError, ValueError):
                values = None

        def row(label: str, name: str) -> None:
            with ui.row().classes("w-full justify-between no-wrap gap-3 items-center"):
                ui.label(label).classes("text-sm ps-muted")
                live[name] = ui.label("").classes("text-sm ps-num text-right")

        with ui.column().classes("w-full ps-card p-3 gap-1"):
            if values is None or "name" not in values:
                draft = d._draft
                ui.label("The whole experiment").classes("ps-section")
                ui.label(f"{len(draft['detectors'])} particle detectors, "
                         f"{len(draft.get('gamma_detectors', []))} γ-ray detectors").classes("text-sm")
            else:
                gamma = key.startswith("gamma")
                ui.label("γ-ray detector" if gamma else "Particle detector").classes("ps-section")
                ui.label(values["name"]).classes("text-base font-medium")
                row("Centre", "place")
                if gamma:
                    row("Half-angle", "half")
                    row("Geometric coverage", "geometric")
                else:
                    row("Covers θ", "covers")
                    row("Solid angle", "omega")
                    row("Hidden by others", "hidden")
                    row("Particles reaching it (estimate)", "rate")
                live["status"] = ui.label("").classes("text-xs ps-warn")
                scene_live(values)
                live["status"].set_text(f"Not allowed here: {values['problem']}" if values.get("problem") else "")
                if not gamma and _scene.on_axis(d.experiment, key):
                    exact_values(key, (("distance", "Distance", "mm"),))
                else:
                    exact_values(key, (("theta", "θ", "deg"), ("phi", "φ", "deg"), ("distance", "Distance", "mm")))
            ui.label("These follow from where things stand in the setup as you have changed it. Apply works out "
                     "the rates, the coincidences and the rest.").classes("text-xs ps-accent mt-1")

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
                run_ = await run.io_bound(P().snapshot().source_run, src["nuclide"], src["activity"], src["time"])
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
                r = await run.io_bound(P().snapshot().tracks, int(count.value), ts["select"], ts["weighted"], 1,
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
        nav = state.get("nav")
        if nav is not None:
            (nav.hide if on else nav.show)()
        for el in (state.get("strip_box"), footer_boxes.get("setup"),
                   *scene_state.get("extras", [])):
            if el is not None:
                el.set_visibility(not on)
        geometry_panel.refresh()
        # The scene follows its element's size once the window says it changed: once right away, and once more
        # after the rail has slid out of the way.
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
                view = SceneView(lambda: scene_setup, state["theme"], scene_selected, scene_live, scene_moved,
                                 height="calc(100vh - 175px)" if full else None)
                view.on_track = scene_track
                scene_state["view"] = view
                scene_state["shows"] = scene_shows()
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
                # Beside the scene: what is selected, then the tables of angles, in the height the scene has.
                with ui.column().classes("gap-2").style("flex: 0 1 460px; min-width: 320px"):
                    selection_panel()
                    geometry_tables()
        if full:
            ui.keyboard(on_key=lambda e: scene_full(False) if e.key.escape and e.action.keydown else None)

    @ui.refreshable
    def geometry_tables():
        scene_state["tables"] = detectors_of(D()._draft)
        rows = [{"detector": d["name"], "omega": _fmt(d["solid_angle_msr"]),
                 "theta": f"{d['theta_range'][0]:.1f}–{d['theta_range'][1]:.1f}",
                 "phi": f"{d['phi_range'][0]:.1f}–{d['phi_range'][1]:.1f}", "segments": d["segments"]}
                for d in D().geometry()["detectors"]]
        ui.table(columns=columns((("detector", "Particle detector"), ("omega", "Ω (msr)"), ("theta", "θ (deg)"),
                                  ("phi", "φ (deg)"), ("segments", "Segments"))), rows=rows).props("dense flat").classes(
            "w-full")
        gammas = D().geometry()["gamma_detectors"]
        if gammas:
            ui.table(columns=columns((("detector", "γ-ray detector"), ("theta", "θ (deg)"),
                                      ("distance", "Distance (mm)"), ("half", "Half-angle (deg)"))),
                     rows=[{"detector": g["name"], "theta": f"{g['theta']:.1f}", "distance": f"{g['distance_mm']:g}",
                            "half": f"{g['half_angle_deg']:.1f}"} for g in gammas]).props("dense flat").classes("w-full")

    @ui.refreshable
    def kinematics_panel():
        plot(shown("kinematics"), "kinematics")

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
            plot(shown("strips", name), "strips", detector=name)

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
                        s = await run.io_bound(P().snapshot().sweep, par.value, values, qty.value, det.value)
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
        plot(shown("energy_loss"), "energy_loss")
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
                fig = await run.io_bound(figure_spectra, P().snapshot(), int(n.value), int(seed.value))
                plot_area.clear()
                with plot_area:
                    plot(fig, "spectra", events=int(n.value), seed=int(seed.value))

            ui.button("Simulate", icon="play_arrow", on_click=simulate)
        plot_area = ui.column().classes("w-full")
        with plot_area:
            plot(shown("spectra", state["events"], 1), "spectra", events=state["events"], seed=1)
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
                g = await run.io_bound(P().snapshot().gamma_spectra, max(state["events"], 400_000), 1)
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

    def panel(title: str, sub: str = ""):
        """A bordered panel with a head line (its name, and after a dot what it shows); its body is returned to
        build in."""
        with ui.element("div").classes("ps-panel"):
            with ui.element("div").classes("ps-panel-head"):
                ui.label(title)
                if sub:
                    ui.label("· " + sub).classes("ps-muted")
            return ui.element("div").classes("ps-panel-body")

    def analysis_block(results=None) -> None:
        """From the γ-ray peak to B(E2): the steps, each adjustable, and what the beam time determines."""
        ui.label("Analysis: from the peak area to B(E2)").classes("font-semibold mt-3")
        ui.label("The events are analysed as an experimentalist would: a gate on the particles in the silicon, the "
                 "Doppler correction for each nucleus, a search for the peaks in each corrected spectrum, a fit of "
                 "each peak over a line and the area under it with the random coincidences subtracted, the yield "
                 "over the efficiency and the angular correlation, a normalisation, and B(E2) from the first-order "
                 "proportionality. Each systematic is found by changing its input by one standard deviation and "
                 "running again.").classes("text-xs ps-muted")
        with ui.row().classes("w-full items-end gap-2"):
            norm = ui.select({"rutherford": "to the particles in the same rings",
                              "target": "to the other nucleus's peak (a known transition)"}, value="rutherford",
                             label="Normalisation").props("dense outlined").classes("w-80")
            gate = ui.select({"inelastic": "inelastic group", "all": "every particle"}, value="inelastic",
                             label="Particle gate").props("dense outlined").classes("w-44")
            corr = ui.select({"emitter": "the excited nucleus", "projectile": "the projectile",
                              "recoil": "the recoil"}, value="emitter", label="Doppler correction").props(
                "dense outlined").classes("w-48")
            width = ui.number("Fit window (FWHM)", value=4.0, min=1.5, max=10, step=0.5).props(
                "dense outlined").classes("w-36")
            rings = ui.input("Rings or strips", placeholder="all, or 4-15").props("dense outlined").classes("w-36")
            named = ui.select({"": "none (the choices here)"} | {g.name: g.name for g in P().gates()}, value="",
                              label="Gate from the Data step").props("dense outlined").classes("w-56").tooltip(
                "A named gate replaces the detectors, rings and particle gate here")
        with ui.row().classes("w-full items-end gap-2"):
            ref_e = ui.input("Reference transition", placeholder="328 keV").props("dense outlined").classes(
                "w-40").tooltip("Only without a level scheme of the other nucleus: with one, its peak in the run "
                                "is the reference, and these two fields stay empty")
            ref_b = ui.input("Its B(E2↑)", placeholder="1.65 e2b2").props("dense outlined").classes("w-40")
            ref_u = ui.number("Its uncertainty (%)", value=3.0, min=0, max=100).props("dense outlined").classes(
                "w-40")
            eff_u = ui.number("Efficiency unc. (%)", value=5.0, min=0, max=100).props("dense outlined").classes(
                "w-40")
            prec = ui.number("Wanted precision (%)", value=5.0, min=0.1, max=100).props("dense outlined").classes(
                "w-40")
        ui.button("Analyse", icon="functions", on_click=lambda: analyse()).props("dense unelevated no-caps")
        # The result goes where the caller says (the Analysis step shows it across the page, under the choices).
        out_box = results if results is not None else ui.column().classes("w-full")

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
            if norm.value == "target" and (ref_e.value or ref_b.value):
                if not ref_e.value or not ref_b.value:
                    ui.notify("Give the reference transition's energy and B(E2↑), or leave both empty to use the "
                              "other nucleus's peak in the run", type="warning")
                    return
                settings["reference"] = {"energy": ref_e.value, "b_up": ref_b.value, "unc": float(ref_u.value) / 100}
            out_box.clear()
            with out_box:
                ui.spinner(size="md")
            r = await run.io_bound(P().snapshot().analyse, settings, max(state["events"], 400_000), 1, named.value or None)
            out_box.clear()
            with out_box:
                show(r, P().analysis(max(state["events"], 400_000), 1).history[-1] if r["available"] else None)
            sync_steps()

        def part(number: int, title: str, text: str) -> None:
            ui.label(f"{number} · {title}").classes("font-semibold mt-2")
            ui.label(text).classes("text-xs ps-muted")

        def sci(x: float) -> str:
            """A large number as a mantissa and a power of ten."""
            if not x or abs(x) < 1e4:
                return f"{x:,.4g}"
            e = int(math.floor(math.log10(abs(x))))
            return f"{x / 10 ** e:.3f} × 10" + str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))

        def facts(rows) -> None:
            """Quiet lines of (what, value)."""
            for what, value in rows:
                with ui.element("div").classes("ps-fact"):
                    ui.label(what)
                    ui.label(value)

        def show_search(r: dict) -> None:
            """The peaks the search found in each corrected spectrum, drawn and listed."""
            search = r.get("search")
            if not search:
                return
            part(1, "Find the peaks", "Each Doppler-corrected spectrum is searched over its whole range. A nucleus's "
                 "lines are sharp in the spectrum corrected for it and smeared in the other, so each nucleus shows "
                 "its peaks in its own spectrum. Every peak found is fitted; one that sits on a transition of the "
                 "level scheme is named.")
            with ui.element("div").classes("w-full grid gap-3").style(
                    "grid-template-columns: repeat(auto-fit, minmax(360px, 1fr))"):
                for which in ("projectile", "recoil"):
                    sp = search[which]
                    with panel(f"Corrected for the {which}", f"the lines of {sp['nuclide']}" if sp.get("nuclide")
                               else ""):
                        ui.label(f"In coincidence with {sp['gate']}; random coincidences subtracted.").classes(
                            "text-xs ps-muted")
                        ui.plotly(themed(figure_search(sp, which), state["theme"])).classes("w-full")
                        if not sp["peaks"]:
                            ui.label("No peak found.").classes("text-sm ps-muted")
            found = [(which, search[which], pk) for which in ("projectile", "recoil") for pk in search[which]["peaks"]]
            if found:
                counts = ", ".join(f"{len(search[which]['peaks'])} corrected for the {which}"
                                   for which in ("projectile", "recoil"))
                with panel("Peaks found", counts):
                    ui.table(columns=columns((("n", "Corrected for"), ("e", "Peak at (keV)"), ("f", "FWHM (keV)"),
                                              ("a", "Area (counts)"), ("s", "Significance"), ("t", "Transition"))),
                             rows=[{"n": f"the {which}" + (f" ({sp['nuclide']})" if sp.get("nuclide") else ""),
                                    "e": f"{pk['energy_kev']:.1f}", "f": f"{pk['fwhm_kev']:.1f}",
                                    "a": f"{pk['area']:,.0f} ± {pk['area_unc']:,.0f}",
                                    "s": f"{pk['significance']:.0f} σ",
                                    "t": (f"{pk['label']} ({pk['line_kev']:.1f} keV)" if pk["label"] else
                                          (f"none (smeared {pk['near']})" if pk.get("near") else "none"))}
                                   for which, sp, pk in found]).props("dense flat").classes("w-full")
                    if any(pk.get("near") for _, _, pk in found):
                        ui.label("Smeared: a line of the other nucleus, which this correction shifts and "
                                 "broadens.").classes("text-xs ps-muted")

        def show_lines(r: dict, beside=None) -> None:
            """Each integrated peak: the fit drawn over the spectrum, its four numbers that matter, and behind a
            fold every number from the area to the yield."""
            lines = r.get("lines")
            if not lines:
                return
            part(2, "Integrate the peaks", "The peak of each nucleus in its own corrected spectrum: a Gaussian on a "
                 "straight line is fitted, and the area is the counts within ±3σ of its centre over that line, "
                 "less the random coincidences expected there. The area over the efficiency and the angular "
                 "correlation is the number of γ rays emitted.")
            with ui.element("div").classes("w-full grid gap-3").style(
                    "grid-template-columns: repeat(auto-fit, minmax(360px, 1fr))"):
                for ln in lines:
                    fit = ln["fit"]
                    with panel(ln["label"], f"{ln['energy_kev']:.1f} keV"):
                        ui.label(f"Corrected for the {ln['correction']}; in coincidence with {ln['gate']}.").classes(
                            "text-xs ps-muted")
                        ui.plotly(themed(figure_peak(ln, float(r["settings"].get("fit_half_width", 4.0))),
                                         state["theme"])).classes("w-full")
                        with ui.element("div").classes("ps-readouts").style(
                                "grid-template-columns: repeat(2, minmax(0, 1fr))"):
                            for title, value in (("Area", f"{ln['area']:,.0f} ± {ln['area_unc']:,.0f}"),
                                                 ("FWHM", f"{1e3 * 2.3548 * fit['sigma']:.1f} keV"),
                                                 ("Efficiency", f"{100 * ln['efficiency']:.3f} %"),
                                                 ("Per particle", f"{ln['per_particle']:.4g}")):
                                with ui.element("div").classes("ps-readout").style("padding: 8px 12px"):
                                    ui.label(title).classes("text-xs ps-muted")
                                    ui.label(value).classes("ps-num")
                        found = ln.get("found")
                        rows = [
                            ("Found by the search", f"at {found['energy_kev']:.1f} keV, {found['significance']:.0f} σ"
                             if found else "no: fitted at the transition energy"),
                            ("Fitted centre and width", f"{1e3 * fit['centroid']:.1f} keV, FWHM "
                                                        f"{1e3 * 2.3548 * fit['sigma']:.1f} keV"),
                            ("Counts over the line, ±3σ", f"{fit['area']:,.1f} ± {fit['area_unc']:,.1f}"),
                            ("Random coincidences in the window", f"{ln['randoms']:,.1f}, subtracted"),
                            ("Area", f"{ln['area']:,.1f} ± {ln['area_unc']:,.1f} counts"),
                            ("Full-energy efficiency", f"{100 * ln['efficiency']:.3f} %"),
                            ("Angular correlation", f"{ln['correlation']:.3f}")]
                        if ln["gate_efficiency"] != 1.0:
                            rows.append(("Share of the excitations the gate keeps", f"{100 * ln['gate_efficiency']:.1f} %"))
                        rows.append(("γ rays emitted (the yield)", f"{ln['yield']:,.4g}"))
                        rows.append(("Per particle in the same rings", f"{ln['per_particle']:.4g}"))
                        if ln.get("expected_per_particle") is not None:
                            rows.append(("Expected from its level scheme, all orders",
                                         f"{ln['expected_per_particle']:.4g} (measured / expected {ln['ratio']:.3f})"))
                        if ln.get("blends"):
                            rows.append(("Also inside the window", ", ".join(ln["blends"])))
                        with ui.expansion("Every number of this peak").props("dense").classes("w-full"):
                            ui.table(columns=columns((("q", "Quantity"), ("v", "Value"))),
                                     rows=[{"q": q, "v": v} for q, v in rows]).props("dense flat")
                if beside is not None:  # what stands beside the fitted peaks, in the same row
                    with ui.column().classes("gap-3").style("min-width: 0"):
                        beside()

        def show_result(r: dict, last=None) -> None:
            """The strip at the top: B(E2), its uncertainty, how it compares with what was put in, the beam time."""
            b = r["b_e2fm4"]
            with ui.element("div").classes("ps-readouts"):
                for title, value, sub, key, main in (
                        ("B(E2↑)", f"{b:.4g} e²fm⁴", f"{r['b_e2b2']:.4g} e²b² · {r['b_wu']:.3g} W.u.", "b_e2", True),
                        ("Uncertainty", f"± {100 * r['total_unc']:.1f} %",
                         f"statistical {100 * r['statistical']:.1f} %, systematic {100 * r['systematic']:.1f} %; "
                         + (f"the Monte Carlo sample itself adds {100 * r['monte_carlo']:.1f} %"
                            if r.get("whole_run_statistical") is None else
                            f"statistics from {r['statistics_from']}; the whole run: "
                            f"{100 * r['whole_run_statistical']:.1f} %"), "uncertainty", False),
                        ("Put in", f"{r['truth_e2fm4']:.4g} e²fm⁴", f"pull {r['pull']:+.2f} σ", None, False),
                        ("Beam time", f"{r['hours_for_precision']:.3g} h",
                         f"for {100 * r['settings']['wanted_precision']:.3g} % statistics; "
                         f"{r['counts_per_shift']:.3g} counts in the peak per {SHIFT_H:g} h shift", None, False)):
                    with ui.element("div").classes("ps-readout" + (" ps-readout-main" if main else "")):
                        with ui.row().classes("items-center no-wrap gap-1"):
                            ui.label(title).classes("text-xs ps-muted")
                            if key:
                                ui.button(icon="help_outline", color=None, on_click=lambda k=key: open_explanation(
                                    k, result=last)).props("dense flat round size=xs").classes("ps-muted").tooltip(
                                    "How this number is obtained")
                        ui.label(value).classes("ps-readout-value")
                        ui.label(sub).classes("text-xs ps-muted")

        def show_chain(r: dict, last=None) -> None:
            """From the area to B(E2): the few numbers of the chain, large, with what is done between them; the
            smaller facts under them; every step in words behind a fold."""
            lines = r.get("lines") or []
            own = next((ln for ln in lines if ln.get("nucleus", 0) == 0), None)
            all_orders = any(step["step"] == "all orders" for step in r["steps"])
            box = ui.element("div").classes("ps-panel")
            with box:
                with ui.element("div").classes("ps-panel-head"):
                    ui.label("3 · From the area to B(E2)")
                if own is not None:
                    by = "÷ efficiency %.4g %% ÷ correlation %.3f" % (100 * own["efficiency"], own["correlation"])
                    if own["gate_efficiency"] != 1.0:
                        by += " ÷ gate %.3f" % own["gate_efficiency"]
                    norm = ("normalised to the other nucleus's peak" if r["settings"].get("normalisation") == "target"
                            else "normalised to the particles")
                    chain = [("Area of the peak", f"{own['area']:,.0f} ± {own['area_unc']:,.0f}", "counts", by),
                             ("γ rays emitted", sci(own["yield"]), "",
                              "÷ %s particles in the same rings" % sci(own["yield"] / own["per_particle"])
                              if own["per_particle"] else ""),
                             ("γ rays per particle", f"{own['per_particle']:.4g}", "",
                              norm + ("; coupled equations of the level scheme" if all_orders
                                      else "; first-order proportionality")),
                             ("B(E2↑)", f"{r['b_e2fm4']:.4g}", f"e²fm⁴ ± {100 * r['total_unc']:.1f} %", "")]
                    with ui.element("div").classes("ps-chain"):
                        for k, (title, value, unit, op) in enumerate(chain):
                            last_one = k == len(chain) - 1
                            with ui.element("div").classes("ps-chain-step" + (" ps-chain-last" if last_one else "")):
                                ui.label(title).classes("text-xs ps-muted")
                                with ui.row().classes("items-baseline gap-2"):
                                    ui.label(value).classes("ps-chain-value")
                                    if unit:
                                        ui.label(unit).classes("text-xs ps-muted ps-num")
                            if op:
                                ui.label(op).classes("ps-chain-op")
                    fit = own["fit"]
                    rows = [("Fitted at", f"{1e3 * fit['centroid']:.1f} keV, FWHM {1e3 * 2.3548 * fit['sigma']:.1f} keV"),
                            ("Randoms subtracted", f"{own['randoms']:,.1f} counts")]
                    first = r.get("b_first_order")
                    if all_orders and first and first == first:
                        rows.append(("First order would read", f"{first:.4g} e²fm⁴"))
                        rows.append(("All orders against first order", f"{100 * (r['b_e2fm4'] / first - 1):+.1f} %"))
                    if own.get("blends"):
                        rows.append(("Also inside the window", ", ".join(own["blends"])))
                    with ui.element("div").classes("ps-facts"):
                        ui.label("Behind the numbers").classes("text-xs ps-muted")
                        facts(rows)
                with ui.expansion("Every step in words").props("dense").classes("w-full"):
                    for step in r["steps"]:
                        with ui.column().classes("gap-0 mb-2"):
                            ui.label(step["step"]).classes("text-xs ps-section")
                            ui.label(step["text"]).classes("text-sm")
                with ui.row().classes("gap-1 items-center").style("padding: 4px 10px 8px"):
                    ui.label("How each is obtained:").classes("text-xs ps-muted")
                    for key, label in (("area", "area"), ("yield", "yield"), ("mean_probability", "⟨P⟩"),
                                       ("b_e2", "B(E2)"), ("weisskopf", "W.u."), ("beta2", "β₂"), ("q0", "Q₀"),
                                       ("lifetime", "lifetime")):
                        ui.button(label, color=None, on_click=lambda k=key: open_explanation(k, result=last)).props(
                            "dense flat no-caps size=sm").classes("ps-accent")

        def show_budget(r: dict) -> None:
            """The systematic uncertainties as bars, the largest first, and the shape the B(E2) stands for."""
            budget = sorted(r["budget"].items(), key=lambda kv: -kv[1])
            top = max([v for _, v in budget] + [1e-12])
            with panel("Systematics", f"{100 * r['systematic']:.1f} % in all"):
                for name, value in budget:
                    with ui.element("div").style("padding: 4px 0"):
                        with ui.element("div").classes("ps-fact").style("padding: 0"):
                            ui.label(name)
                            ui.label(f"{100 * value:.2f} %")
                        with ui.element("div").classes("ps-bar"):
                            ui.element("div").style(f"width: {100 * value / top:.0f}%")
            sh = r["shape"]
            if sh:
                rows = [("β₂", f"{sh['beta2']:.3f}"), ("B(E2↑) in W.u.", f"{r['b_wu']:.3g}"),
                        ("Q₀ of a rigid rotor", f"{sh['q0_efm2']:.3g} e fm²"),
                        ("Q_s(2⁺) of that rotor", f"{sh['qs_2plus_efm2']:.3g} e fm²")]
                if "e4_over_e2" in sh:
                    rows.append(("E(4⁺)/E(2⁺)", f"{sh['e4_over_e2']:.3f}"))
                with panel("Shape"):
                    facts(rows)
                    ui.label(sh["note"]).classes("text-xs ps-muted")

        def show(r: dict, last=None) -> None:
            if not r["available"]:
                ui.label(r["reason"]).classes("text-sm ps-warn")
                return
            b = r["b_e2fm4"]
            if b != b:
                show_search(r)
                show_lines(r)
                for step in r["steps"]:
                    with ui.row().classes("items-start no-wrap gap-2"):
                        ui.label(step["step"]).classes("text-xs ps-section w-40")
                        ui.label(step["text"]).classes("text-sm")
                ui.label("No result: " + "; ".join(r["notes"])).classes("text-sm ps-warn")
                return
            show_result(r, last)
            with ui.row().classes("w-full items-start gap-3").style("flex-wrap: wrap"):
                with ui.column().classes("gap-2").style("flex: 3 1 560px; min-width: 0"):
                    show_search(r)
                with ui.column().classes("gap-3 mt-2").style("flex: 1 1 320px; min-width: 280px"):
                    show_chain(r, last)
            if r.get("lines"):
                show_lines(r, lambda: show_budget(r))
            else:
                show_budget(r)
            for note in r["notes"]:
                ui.label(note).classes("text-xs ps-muted")

        # The analysis kept with the current run is shown at once: the result that the step bar and the report give.
        kept = P().run_analysis() if P().run is not None else None
        if kept is not None:
            from dataclasses import asdict

            shown = dict(asdict(kept), available=True, total_unc=kept.total_unc)
            with out_box:
                ui.label(f"The analysis kept with run {P().run.number}"
                         + ("" if kept.lines else " (made before the peaks were searched for: press Analyse to see "
                                                  "them found and integrated)") + ":").classes("text-xs ps-muted")
                show(shown, kept)

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
            fit_too = ui.checkbox("Fit the offset back (15 s on a sample, a few minutes on a long run)",
                                  value=True)
        box = ui.column().classes("w-full")

        async def check() -> None:
            box.clear()
            with box:
                ui.spinner(size="md")
            r = await run.io_bound(P().snapshot().alignment, float(off.value), max(state["events"], 400_000), 1,
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
                             "the Data step)").classes("ps-section mt-2")
                    ui.plotly(themed(figure_gamma_crystal(view), state["theme"])).classes("w-full")
                if "fit" in r:
                    f = r["fit"]
                    if f.get("against") == "flat":
                        judged = ("flat lines: a simulation with the assumed geometry gave too few γ rays "
                                  "for a pattern")
                    elif f.get("reference_at_mm"):
                        judged = "a simulation with the target where the fit put it"
                    else:
                        judged = "a simulation with the assumed geometry"
                    ui.label(f"Fitted: the target is {f['offset_mm']:+.2f} ± {f['uncertainty_mm']:.2f} mm from "
                             f"where the analysis assumed it (truth {f['truth_mm']:+g} mm; χ² {f['chi2']:.0f} for "
                             f"{f['degrees_of_freedom']} degrees of freedom, judged against {judged}). The scene "
                             "shows the assumed target as a faint outline.").classes("text-sm")

        ui.button("Check the alignment", icon="straighten", on_click=check).props("dense flat no-caps")

    @ui.refreshable
    def trajectories_panel():
        t = P().trajectories()
        ui.label(f"{P().experiment.beam.nuclide} on {t['target']}: head-on distance d₀ = {t['d0_fm']:.2f} fm, "
                 f"nuclear range {t['interaction_radius_fm']:.1f} fm, grazing angle "
                 f"{t['grazing_angle_deg']:.1f}° (CM).").classes("text-sm")
        plot(shown("trajectories"), "trajectories")

    def lookup_levels(role: str):
        def run() -> None:
            try:
                direct(lambda planner: planner.lookup_levels(role))
            except (LookupError, ValueError) as e:
                ui.notify(str(e), type="warning", multi_line=True)
        return run

    def remove_levels(role: str):
        return lambda: direct(lambda planner: planner.remove_levels(role))

    def use_state(role: str, level: int, multipolarity: str):
        return lambda: direct(lambda planner: planner.use_state(role, level, multipolarity))

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
                direct(lambda planner: planner.set_matrix_element(role, m["from"], m["to"], m["multipolarity"], text))
            except ValueError as e:
                ui.notify(str(e), type="negative", multi_line=True)
                box.value = shown

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
            ui.plotly(shown("levels", role)).classes("w-full")
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
        plot(shown("gamma"), "gamma")
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
                      on_change=lambda e: direct(lambda planner: planner.set("reaction", "emission", e.value))
                      ).props("dense no-caps unelevated")
            if g["emission"] == "correlated":
                ui.label("A nucleus that decays after leaving the target is deoriented by the hyperfine field of "
                         "its highly charged ion (the two-state model, with GOSIA's parameters): its γ rays come "
                         "out more evenly, the more so the longer its level lives.").classes("text-xs ps-muted")
                ui.toggle({"none": "No deorientation", "vacuum": "Deorientation in vacuum"},
                          value=g.get("deorientation", "none"),
                          on_change=lambda e: direct(lambda planner: planner.set("reaction", "deorientation",
                                                                                 e.value))
                          ).props("dense no-caps unelevated")
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
            m = await run.io_bound(P().snapshot().multistep)
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
                f = await run.io_bound(P().snapshot().fit_matrix_elements)
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

    def report_panel():
        """The report of the experiment: one document (:mod:`physim.nuclear.app.tab_report`)."""
        tab_report.experiment_block(ctx)

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
                    ui.plotly(shown("efficiency", f"gamma:{i}")).classes("w-full")

    @ui.refreshable
    def plan_answers() -> None:
        tab_plan.render_answers(ctx)

    panels = {"geometry": geometry_panel, "kinematics": kinematics_panel, "rates": rates_panel,
              "energy_loss": energy_loss_panel, "spectra": spectra_panel, "trajectories": trajectories_panel,
              "gamma": gamma_panel, "report": report_panel, "efficiency": efficiency_panel,
              "plan": plan_answers}

    #: The container of each stage's content; the stages whose content is stale (rebuilt when shown, or by the
    #: pass that follows a change); and the stages built at least once since the tabs were laid out.
    stage_boxes: dict = {}
    dirty: set = set()
    built: set = set()
    #: How many times each stage was built (a build ahead of time gives way to a later one).
    builds: dict = {}
    #: After the stage on screen, the others in the order one is likely to look at them.
    look_order = ("plan", "setup", "run", "data", "analysis", "report", "physics")

    def stage_parts(key: str) -> list:
        """A stage as the steps that build it (a module may give several, see ``tab_physics.parts``)."""
        module = stage_modules[key]
        return module.parts(ctx) if hasattr(module, "parts") else [lambda: module.render(ctx)]

    def build_part(part) -> None:
        try:
            part()
        except Exception as err:  # noqa: BLE001 -- say why instead of an empty tab
            ui.label(f"This tab could not be shown: {err}").classes("text-sm ps-warn")

    def render_stage(key: str) -> None:
        """Build a stage's content into its container (emptying what was there)."""
        box = stage_boxes.get(key)
        if box is None:
            return
        box.clear()
        builds[key] = builds.get(key, 0) + 1
        with box:
            for part in stage_parts(key):
                build_part(part)
        dirty.discard(key)
        built.add(key)

    async def render_stage_ahead(key: str, superseded) -> None:
        """Build a stage that is not on screen, a part per turn of the event loop, so a click on the page is
        answered between the parts. It gives way if the stage is built in one go meanwhile (its tab was chosen)
        or the setup changed."""
        box = stage_boxes.get(key)
        if box is None:
            return
        box.clear()
        builds[key] = n = builds.get(key, 0) + 1
        for part in stage_parts(key):
            if superseded() or builds[key] != n or box.is_deleted:
                return
            with box:
                build_part(part)
            await asyncio.sleep(0.02)
        if builds[key] == n and not superseded():
            dirty.discard(key)
            built.add(key)

    def waiting(key: str) -> None:
        """A stage's container while its numbers are computed."""
        box = stage_boxes.get(key)
        if box is None:
            return
        box.clear()
        with box:
            with ui.row().classes("items-center gap-2 ps-muted"):
                ui.spinner(size="md")
                ui.label("Computing…").classes("text-sm")

    async def build_stage(key: str) -> None:
        """Bring a stale stage up to date: its numbers on the worker (a spinner in the tab if that takes more than
        a moment), then its content from what is kept. Nothing is computed on the event loop, so the page keeps
        answering. A change of the setup meanwhile supersedes this: its own refresh rebuilds the stage."""
        token = state.get("generation", 0)
        task = asyncio.ensure_future(run.io_bound(prepare, key, P().snapshot()))
        done, _ = await asyncio.wait({task}, timeout=0.15)
        if not done:
            if key in dirty:
                waiting(key)
            await task
        if gone() or state.get("generation", 0) != token or key not in dirty:
            return
        render_stage(key)

    async def show_stage(key: str) -> None:
        """A tab was chosen: it is shown as it is if it was built for the current setup, and built now if not."""
        state["stage"] = key
        if key in dirty:
            await build_stage(key)

    def go_to(key: str) -> None:
        """Show a stage, as choosing it in the step bar does."""
        if state.get("tabs") is not None and state["stage"] != key:
            state["tabs"].set_value(key)

    def advance(key: str) -> None:
        """The button at the foot of a step: on to the next one, from its top."""
        go_to(key)
        ui.run_javascript("window.scrollTo({top: 0})")

    def working(on: bool) -> None:
        """The spinner at the top of the rail: on while anything asked for is being computed."""
        state["working"] = max(0, state.get("working", 0) + (1 if on else -1))
        if state.get("busy") is not None:
            state["busy"].set_visibility(state["working"] > 0)

    def refresh_current() -> None:
        """Rebuild the stage on screen from what is kept (the scene in place, so the camera stays where it is)."""
        key = state["stage"]
        if key == "setup" and "setup" in built:
            if not scene_current():
                geometry_panel.refresh()
            reading_box.refresh()
            setup_strip.refresh()
            checks_line.refresh()
            dirty.discard("setup")
            return
        render_stage(key)

    async def build_ahead(p: Planner, superseded) -> None:
        """The stages not on screen, one by one in the order one is likely to look at them: the numbers on the
        worker, then the content, so each is there when its tab is chosen. Then what the Plan tab shows only when
        a check is opened, and last a run's events (seconds, and memory, for a long run)."""
        last = ["data"] if p.run is not None else []
        for key in [k for k in look_order if k not in last] + ["plan+"] + last:
            if key != "plan+" and key not in dirty:
                continue
            await asyncio.sleep(0.05)  # the page's own events go first
            await run.io_bound(prepare, key.rstrip("+"), p, key == "plan+")
            if superseded():
                return
            if key in dirty:
                await render_stage_ahead(key, superseded)

    def refresh_results() -> None:
        """After the setup, the run or the theme changed. Every stage is now stale, and those not on screen are
        emptied. The numbers are computed on the worker (a spinner in the header while they are); then the strip,
        the setup panel and the stage on screen are rebuilt from them, and after that the other stages one by
        one, so each is there when its tab is chosen. Another change in the meantime supersedes this one."""
        state["generation"] = token = state.get("generation", 0) + 1
        current = state["stage"]
        for key, box in stage_boxes.items():
            dirty.add(key)
            if key != current:
                box.clear()
        working(True)

        def superseded() -> bool:
            return gone() or state.get("generation") != token

        async def numbers(p: Planner) -> None:
            await run.io_bound(with_work, warm, p)
            if superseded():
                return
            state["ready"] = True
            typing = False
            if client.has_socket_connection:
                try:  # the cards are not rebuilt under a field being typed in; the next blur brings them up to date
                    typing = await ui.run_javascript(
                        "!!document.activeElement && ['INPUT', 'TEXTAREA'].includes(document.activeElement.tagName)",
                        timeout=1.0)
                except Exception:  # noqa: BLE001
                    typing = False
                if superseded():
                    return
            if state.pop("panel_due", False) or not typing:
                setup_panel.refresh()
            else:
                sync_cards()
            pending_bar.refresh()
            pending_note.refresh()
            sync_foot()
            sync_steps()
            while True:  # the stage on screen (the user may choose another tab while its numbers are computed)
                key = state["stage"]
                await run.io_bound(prepare, key, p)
                if superseded():
                    return
                if key == state["stage"]:
                    break
            if key in dirty:
                refresh_current()

        async def go() -> None:
            p = P().snapshot()
            try:
                await numbers(p)
            finally:
                working(False)
            if not superseded():
                await build_ahead(p, superseded)

        later(go)

    def refresh_data() -> None:
        """A choice on the Data tab changed: what it changes is computed on the worker (the spinner in the header
        if that takes more than a moment), then the view is rebuilt from what is kept."""
        state["data_choice"] = n = state.get("data_choice", 0) + 1
        token = state.get("generation", 0)

        async def go() -> None:
            if P().run is not None:
                task = asyncio.ensure_future(run.io_bound(prepare, "data", P().snapshot()))
                done, _ = await asyncio.wait({task}, timeout=0.15)
                if not done:
                    working(True)
                    try:
                        await task
                    finally:
                        working(False)
                if gone() or state.get("data_choice") != n or state.get("generation", 0) != token:
                    return
            data_view.refresh()
            sync_steps()

        later(go)

    @ui.refreshable
    def reading_box(tab: str) -> None:
        try:
            lines = reading(tab)
        except Exception:  # noqa: BLE001 -- a reading must never break the page
            lines = []
        if not lines:
            return
        with ui.row().classes("w-full items-start no-wrap gap-2 ps-note rounded p-2"):
            ui.icon("lightbulb").classes("ps-accent mt-0.5")
            with ui.column().classes("gap-1"):
                ui.label("How to read this").classes("text-xs font-semibold ps-muted")
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
        # The section may still be on its way (the tab is built when first shown): look for it for a few seconds.
        ui.run_javascript(
            "(() => { let n = 0; const go = () => { const el = document.querySelector('.physics-%s'); "
            "if (el) el.scrollIntoView({behavior: 'smooth', block: 'start'}); "
            "else if (n++ < 60) setTimeout(go, 100); }; setTimeout(go, 300); })()" % section)

    def block(tab: str, title: Optional[str] = None) -> None:
        """One of the planner's results: its heading, how to read it, the panel, and its explanation."""
        if title:
            ui.label(title).classes("text-lg font-semibold mt-2")
        if tab == "geometry":  # the scene's full-screen mode hides what stands above and below it
            with ui.element("div").classes("w-full flex flex-col gap-2") as above:
                with ui.row().classes("w-full items-stretch gap-2").style("flex-wrap: wrap"):
                    with ui.element("div").style("flex: 1 1 420px; min-width: 0"):
                        setup_strip()
                    ui.button("Edit the setup", icon="tune", on_click=edit_setup).props(
                        "unelevated no-caps").tooltip("The beam, the target and reaction, the detectors and the "
                                                      "run conditions, in a window of their own")
                checks_line()
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
        sync_steps()
        setup_panel.refresh()
        refresh_results()

    def start_run(f: dict) -> None:
        if pending():
            ui.notify("Changes of the setup are not applied yet: Apply or Discard them, then start the run.",
                      type="warning", multi_line=True)
            return
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

    def retake_run() -> None:
        try:
            out = P().retake_run(background=True)
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
        sync_steps()
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

    #: The step bar's elements, {step: (its tab, the line under its name)}, and the foot of each stage.
    step_els: dict = {}
    footer_boxes: dict = {}

    def sync_foot() -> None:
        """The foot of each step: what there is to do next, and the button that goes there
        (:func:`physim.nuclear.workbench.next_step`)."""
        waiting_changes = len(pending())
        for key, box in footer_boxes.items():
            if box.is_deleted:
                continue
            box.clear()
            try:
                nxt = wb.next_step(P(), key, waiting_changes)
            except Exception:  # noqa: BLE001 -- the foot must never break the page
                nxt = None
            if nxt is None:
                continue
            with box:
                with ui.row().classes("w-full items-center gap-3 ps-foot"):
                    ui.label(nxt["hint"]).classes("text-sm ps-muted")
                    ui.space()
                    if nxt["label"]:
                        ui.button(nxt["label"], on_click=lambda to=nxt["to"]: advance(to)).props(
                            "unelevated no-caps icon-right=arrow_forward")

    def sync_steps() -> None:
        """The step bar, the foot of each step and the experiment's name in the header, from the experiment as
        it stands (:func:`physim.nuclear.workbench.steps`): after the setup's numbers are there, a run is taken or
        chosen, a gate is saved, a run is analysed. While the numbers of a change are on their way the bar keeps
        what it showed; the pass that brings them calls this again."""
        experiment_button.refresh()
        sync_foot()
        if not state["ready"] or not numbers_kept():
            return
        try:
            found = wb.steps(P())
        except Exception:  # noqa: BLE001 -- the step bar must never break the page
            return
        for step in found:
            tab, status = step_els.get(step["key"], (None, None))
            if tab is None or tab.is_deleted:
                continue
            tab.classes(remove="ps-step-done ps-step-ready ps-step-waiting", add=f"ps-step-{step['state']}")
            status.set_text(step["status"])
            status.classes(remove="ps-warn", add="ps-warn" if step["warn"] else "")

    @ui.refreshable
    def main_area() -> None:
        step_els.clear()
        footer_boxes.clear()
        # The steps of the experiment, in order down the rail on the left, each with its state and a line of
        # status; the Physics reference is a stage without a step (the rail's button and the ? beside a number
        # show it).
        rail = state["steps_box"]
        rail.clear()
        with rail:
            with ui.tabs(value=state["stage"], on_change=lambda e: show_stage(e.value)).classes("ps-steps").props(
                    "vertical dense no-caps") as tabs:
                for key, label in wb.STAGES:
                    tab = ui.tab(key, label=label)
                    if key not in dict(wb.STEPS):
                        tab.classes("ps-step-hidden")
                        continue
                    with tab:
                        step_els[key] = (tab, ui.label("").classes("ps-step-status"))
        state["tabs"] = tabs
        # The stages are built by refresh_results, once their numbers are computed on the worker: the one on
        # screen first, then the others. ``prebuild`` builds every stage now (the served-page test).
        with ui.tab_panels(tabs, value=state["stage"]).classes("w-full"):
            for key, _ in wb.STAGES:
                with ui.tab_panel(key):
                    stage_boxes[key] = ui.column().classes("w-full gap-2")
                    footer_boxes[key] = ui.element("div").classes("w-full")
                    dirty.add(key)
        built.clear()
        if prebuild:
            for key, _ in wb.STAGES:
                render_stage(key)
        else:
            waiting(state["stage"])
        sync_steps()

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
            saved_theme, system_dark, auto = await ui.run_javascript(
                "(() => { let t = null, a = null; try { t = localStorage.getItem('physim-planner-theme'); "
                "a = localStorage.getItem('physim-planner-auto-apply') } catch (e) {} "
                "return [t, window.matchMedia('(prefers-color-scheme: dark)').matches, a] })()", timeout=3)
        except Exception:  # noqa: BLE001 -- no answer from the browser: keep the defaults
            return
        if auto == "on" and not state["auto_apply"]:
            state["auto_apply"] = True
            pending_bar.refresh()
        if theme is None:
            wanted = saved_theme if saved_theme in THEMES else ("dark" if system_dark else "light")
            if wanted != state["theme"]:
                theme_toggle.value = wanted

    ctx = SimpleNamespace(
        ui=ui, P=P, state=state, edit=edit, section=section, help_icon=help_icon, bind=bind,
        reaction_section=reaction_section, clover_switches=clover_switches, add_detector=add_detector,
        D=D, cards=cards, checks=checks, checks_box=checks_box, card_parts=card_parts,
        mark=mark,
        duplicate_detector=duplicate_detector, remove_detector=remove_detector,
        add_gamma_detector=add_gamma_detector, duplicate_gamma_detector=duplicate_gamma_detector,
        remove_gamma_detector=remove_gamma_detector, open_card=open_card, block=block, readouts=readouts, panel=panel, analysis_block=analysis_block, alignment_block=alignment_block,
        multistep_block=multistep_block, open_explanation=open_explanation, efficiency_panel=efficiency_panel,
        plan_answers=plan_answers, run_view=run_panel, run_form=run_form, refresh_run_form=run_panel.refresh,
        start_run=start_run, stop_run=stop_run, extend_run=extend_run, retake_run=retake_run, load_run=load_run, watch_events=watch_events,
        columns=columns, data_view=data_view, refresh_data=refresh_data, data_options=data_options,
        keep=keep, in_theme=in_theme, when_opened=when_opened,
        gate_form=gate_form, fmt=_fmt, time_text=_time,
        analysis_view=analysis_view, report_view=report_view, show_physics=show_physics, levels_block=levels_block,
        geometry_tables=geometry_tables, explain=explain, reading_box=reading_box, gamma_panel=gamma_panel,
        go_to=go_to)

    @ui.refreshable
    def experiment_button() -> None:
        """At the top of the rail: the experiment's name, with what is done to an experiment as its menu, and
        under it the beam and the target."""
        p = P()
        d = p.draft
        with ui.column().classes("w-full gap-1"):
            with ui.button(p.name, color=None).props("flat no-caps icon-right=expand_more").classes(
                    "w-full ps-experiment"):
                with ui.menu():
                    ui.menu_item("Open an experiment…",
                                 on_click=lambda: (experiments_list.refresh(), experiments_dialog.open()))
                    ui.menu_item("New experiment…", on_click=lambda: new_dialog.open())
                    ui.separator()
                    ui.menu_item("Import a setup file…", on_click=lambda: upload_dialog.open())
                    ui.menu_item("Export the setup…", on_click=export_setup)
            ui.label(f"{d['beam'].get('nuclide', '?')} at {d['beam'].get('energy', '?')} on "
                     f"{d['target'].get('material', '?')} {d['target'].get('thickness', '')}".strip()).classes(
                "text-xs ps-muted").style("padding-left: 6px")

    # -- layout -----------------------------------------------------------------------------------------------
    ui.add_css(STYLE)
    ui.add_css(PENDING_STYLE)
    ui.add_css(STEP_STYLE)
    ui.add_css(record_css(".ps-record"))
    dark = ui.dark_mode(state["theme"] == "dark")
    # The rail on the left: the experiment, the steps (built by main_area), and at its foot the Physics reference
    # and the theme.
    with ui.left_drawer(value=True).classes("ps-nav").props("width=232 bordered behavior=desktop") as nav:
        with ui.column().classes("w-full gap-3 ps-nav-col"):
            with ui.row().classes("w-full items-center gap-2 no-wrap"):
                ui.label("physim").classes("ps-wordmark")
                state["busy"] = ui.spinner("dots", size="sm").classes("ps-muted").tooltip("Computing the numbers")
                state["busy"].set_visibility(False)
            experiment_button()
            state["steps_box"] = ui.element("div").classes("w-full")
            ui.space()
            with ui.row().classes("w-full items-center justify-between gap-1 no-wrap"):
                ui.button("Physics", icon="menu_book", color=None, on_click=lambda: go_to("physics")).props(
                    "flat no-caps dense").classes("ps-muted").tooltip(
                    "The physics behind the numbers, and where each part is validated")
                theme_toggle = ui.toggle({"light": "Light", "dark": "Dark"}, value=state["theme"],
                                         on_change=lambda e: set_theme(e.value)).props("dense no-caps unelevated")
    state["nav"] = nav
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
    # The setup is edited in a window of its own, at the right of the page so the scene stays in view beside it:
    # "Edit the setup" on the Setup step opens it, "Done editing" applies what was changed and closes it.
    with ui.dialog().props("position=right full-height") as setup_dialog:
        with ui.card().classes("ps-setup-window"):
            with ui.row().classes("w-full items-center no-wrap ps-panel-head"):
                ui.label("Edit the setup")
                ui.space()
                ui.button(icon="close", color=None, on_click=setup_dialog.close).props(
                    "flat round dense").classes("ps-muted").tooltip("Close. Changes not applied wait: the page "
                                                                    "says so, with Apply and Discard.")
            with ui.element("div").classes("ps-setup-body"):
                pending_bar()
                setup_panel()
            with ui.row().classes("w-full items-center gap-2 ps-setup-foot"):
                ui.button("Done editing", icon="check", on_click=done_editing).props("unelevated no-caps").tooltip(
                    "Apply what was changed and close (Ctrl+Enter applies without closing)")
                ui.button("Cancel", color=None, on_click=cancel_editing).props("flat no-caps").classes(
                    "ps-muted").tooltip("Close without the changes made here")
    with ui.column().classes("w-full gap-2"):
        with ui.element("div").classes("w-full") as strip_box:
            pending_note()
        state["strip_box"] = strip_box
        main_area()
    if offer_list:
        experiments_dialog.open()
    if not prebuild:
        refresh_results()  # the numbers on the worker, then the panel, the strip and the stages
    # Ctrl+Enter applies the changes from anywhere on the page, a field being typed in included (its own Enter
    # reaches the server first and puts its text in the draft).
    ui.on("physim_apply", lambda: apply_changes())
    ui.add_head_html("<script>document.addEventListener('keydown', (e) => { if (e.key === 'Enter' && "
                     "(e.ctrlKey || e.metaKey) && !e.repeat) emitEvent('physim_apply'); });</script>")
    ui.timer(0.2, restore_choices, once=True)
    ui.timer(0.5, tick)
