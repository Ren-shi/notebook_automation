"""The Physics tab (backlog item 70): the physics behind the numbers, as a reference, in the order of the physics
register: nuclear data, kinematics, energy loss, the Rutherford orbit, the detectors and the γ-ray response, rates
and events, Coulomb excitation. Each section links to its register page and the tests behind it; the ? of a
number and the consequence line of a setup card open the section that shows its physics. Nothing here needs a
run: it is computed from the setup."""

from __future__ import annotations

from .. import workbench as wb
from ..record import DOCS

#: The planner panels in each section (planner tab names, or "efficiency" and "levels").
PANELS = {"data": ("levels",), "kinematics": ("kinematics",), "stopping": ("energy_loss",),
          "rutherford": ("trajectories",), "detectors": ("geometry_tables", "efficiency"), "rates": ("rates_explain",),
          "coulex": ("gamma",)}


def anchor(key: str) -> str:
    return f"physics-{key}"


def render(ctx) -> None:
    for part in parts(ctx):
        part()


def parts(ctx) -> list:
    """The tab as steps, the head and then one per section, each building its share into the current container.
    The page takes them in one go when the tab is on screen, and one per turn of the event loop when it builds
    the tab ahead of time, so the page's own events go between them."""
    return [lambda: head(ctx)] + [lambda entry=entry: section(ctx, *entry) for entry in wb.PHYSICS]


def head(ctx) -> None:
    ui = ctx.ui
    ui.label("The physics behind the numbers, in the order of the physics register. Each section says where it "
             "is validated; the ? beside a number on a step opens the section that shows its physics.").classes(
        "text-sm ps-muted")
    with ui.row().classes("gap-1"):
        for key, title, _, _ in wb.PHYSICS:
            ui.button(title.split(":")[0], on_click=lambda k=key: ctx.show_physics(k)).props(
                "dense flat no-caps size=sm")


def section(ctx, key: str, title: str, page: str, tests) -> None:
    ui = ctx.ui
    with ui.column().classes(f"w-full gap-1 mt-3 {anchor(key)}"):
        with ui.row().classes("w-full items-baseline gap-3"):
            ui.label(title).classes("text-lg font-semibold")
            ui.link("Register page", DOCS + page + ".html", new_tab=True).classes("text-xs")
            ui.label("Validated by " + ", ".join(tests)).classes("text-xs ps-muted")
        if key == "coulex" and ctx.P().experiment.excitation is None:
            ui.label("The setup has no excited state: set Coulomb excitation in the Target and reaction card.").classes(
                "text-sm ps-muted")
            return
        for panel in PANELS[key]:
            if panel == "levels":
                ctx.levels_block()
            elif panel == "efficiency":
                if ctx.P().experiment.gamma_detectors:
                    ui.label("The γ-ray detectors' efficiency against energy").classes("ps-section mt-1")
                    ctx.efficiency_panel()
            elif panel == "geometry_tables":
                ctx.geometry_tables()
                ctx.explain("geometry")
            elif panel == "rates_explain":
                ctx.reading_box("rates")
                ctx.explain("rates")
                ctx.explain("spectra")
            elif panel == "gamma":
                ctx.reading_box("gamma")
                ctx.gamma_panel(False)
                ctx.explain("gamma")
            else:
                ctx.block(panel)
