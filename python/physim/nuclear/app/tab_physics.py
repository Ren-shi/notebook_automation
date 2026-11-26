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


#: The deorientation in vacuum, as the Physics tab explains it (docs/theory/coulex.md, "Deorientation in vacuum").
#: Markdown: the asterisk of Λ* and the underscores of τ_c and G_k are escaped.
DEORIENTATION_TEXT = (
    "**The model.** GOSIA's two-state model (Bosch and Spehl 1977; Brenn et al. 1977; the GOSIA manual). A nucleus "
    "that has left the target is an ion whose electron shells are still decaying: the hyperfine field at the nucleus "
    "fluctuates, turned at random every τ\\_c, until at the rate Λ\\* the ion reaches its ground state and the field "
    "becomes static. The field turns the nuclear spin, so the orientation the collision gave fades: the tensors of "
    "rank k of a level that has been in vacuum for a time t are G\\_k(t) times what the collision (or the decay that "
    "fed it) left. Averaged over the decay of a level of mean life τ this is GOSIA's G\\_k. 1 is no loss.\n\n"
    "**Its parameters**, GOSIA's defaults (its input VAC): atomic spin J = 3, Γ = 0.02 /ps, Λ\\* = 0.0345 /ps, "
    "τ\\_c = 3.5 ps, field H = K Z β^x with K = 600 T and x = 0.6. The g factor is Z/A, as GOSIA has it, unless the "
    "setup gives one (`[reaction] deorientation_g`, for every level of the excited nucleus); the other nucleus of the "
    "collision has Z/A. Only the g factor can be set.\n\n"
    "**Where it acts.** Only in vacuum: a level is deoriented while it lives outside the target and what is behind "
    "it, from where its nucleus leaves the target (each simulated event with its own depth, direction and speed), or "
    "from its feeding if it is fed outside. A nucleus that decays in the target, or stops there, is not deoriented; "
    "nor is a level without a half-life in the level scheme, which decays at once. Down a cascade each level is "
    "deoriented before it emits and before it feeds the next one, as GOSIA feeds with the deoriented tensors. The "
    "plan's angular correlation, so the coincidence rates and the analysis, the simulated events, the stored runs "
    "and the GOSIA input of the all-orders solution have it. A run stored before physim had it is loaded without "
    "it, as its events were made.\n\n"
    "**Checked against GOSIA, blind.** 174 coefficients G₂, G₄, G₆ of levels of ⁶²Ni, ¹⁹⁴Pt, ¹⁵²Sm, ⁴⁸Ti and ⁹⁸Ru "
    "at v/c from 0.02 to 0.08, with GOSIA's default parameters and with one set where all of them were changed, "
    "the app's numbers stored before GOSIA was run: none missed the band (0.002 + half the last digit GOSIA prints). "
    "The largest difference, 0.00058, is GOSIA's g = Z/mass against physim's Z/A.\n\n"
    "**Not included.** The fields in solids (a nucleus in the target, its backing or a stopper), the transient field "
    "of a nucleus crossing a magnetised foil, and the deorientation of fusion-evaporation residues. A lifetime "
    "measurement has its own: a plunger only the static model of its `lifetime.deorientation` (or none), DSAM none.")


def deorientation(ctx) -> None:
    """The deorientation in vacuum of the setup: whether it is in effect and why, G₂ and G₄ of every level with a
    half-life of both nuclei at a few speeds, and the model's explanation."""
    ui = ctx.ui
    try:
        d = wb.deorientation_levels(ctx.P())
    except Exception:  # noqa: BLE001 - this part must never break the tab
        return
    with ui.column().classes("w-full gap-1 ps-deorientation"):
        _deorientation_body(ctx, d)


def _deorientation_body(ctx, d: dict) -> None:
    ui = ctx.ui
    ui.label("Deorientation in vacuum").classes("font-semibold mt-2")
    if d["in_effect"] == "vacuum" and d["g"] is not None:
        head = (f"On for {d['nuclide'] or 'the excited nucleus'}, g = {d['g']:.3g}"
                + (" (Z/A)" if d["g_source"] == "Z/A" else " (set)") + "; the other nucleus has Z/A. ")
    else:
        head = ""
    ui.label(head + d["why"].replace(" The Physics tab gives each level's G₂ and G₄.", "")).classes("text-sm")
    if d["levels"] and d["in_effect"] == "vacuum":
        ui.label("G₂ / G₄ of each level with a half-life, averaged over its decay, for a nucleus that is in vacuum "
                 "from the reaction on (one that decays in the target keeps more of its orientation); 1 is no loss. "
                 "The speed of the excited nucleus at each particle detector is in the table of particle energies "
                 "above.").classes("text-xs ps-muted")
        cols = (("nucleus", "Nucleus"), ("level", "Level"), ("energy", "Energy (keV)"), ("t", "Half-life (ps)"))
        cols += tuple((f"b{i}", f"β = {b:g}") for i, b in enumerate(wb.DEORIENTATION_BETAS))

        def cell(pair) -> str:
            g2, g4 = pair
            return f"{g2:.3f}" + ("" if g4 is None else f" / {g4:.3f}")

        ui.table(columns=ctx.columns(cols),
                 rows=[dict({"nucleus": r["nucleus"], "level": r["level"], "energy": f"{r['energy_kev']:g}",
                             "t": ctx.fmt(r["half_life_ps"], 3)},
                            **{f"b{i}": cell(r["g"][b]) for i, b in enumerate(wb.DEORIENTATION_BETAS)})
                       for r in d["levels"] if "g" in r]).props("dense flat")
    with ui.expansion("Explain: Deorientation in vacuum", icon="school").classes("w-full ps-sunk mt-2"):
        ui.markdown(DEORIENTATION_TEXT)


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
                deorientation(ctx)
            else:
                ctx.block(panel)
