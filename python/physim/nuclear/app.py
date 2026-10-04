"""The experiment planner web app: build a setup by pointing and clicking, see every result, export the report.

Start it with ``physim app`` (or ``python -m physim app``); it opens in the browser at http://localhost:8080. It runs
on your computer only: no accounts, nothing uploaded. Needs the ``app`` extra: ``pip install physim-engine[app]``.

The app is a thin layer over :class:`physim.nuclear.planner.Planner`: every number it shows can be had from Python
too. The figure functions here (``figure_geometry``, ``figure_kinematics``, ...) return Plotly figures and work in
a notebook as well.
"""

from __future__ import annotations

import io
import math
import os
import time
import zipfile
from typing import Optional

import numpy as np

from . import guide
from .planner import TABS, Planner

#: Fields of each setup section shown in the setup panel, in order: (field, label, placeholder).
BEAM_FIELDS = [("nuclide", "Nuclide", "4He"), ("energy", "Energy", "5.5 MeV or 4 MeV/u"),
               ("current", "Current", "1 pnA or 10 enA"), ("charge_state", "Charge state", "6"),
               ("energy_spread", "Energy spread (FWHM)", "0.1 %"), ("spot_size", "Spot size (FWHM)", "2 mm")]
TARGET_FIELDS = [("material", "Material", "Au, 208Pb, CD2"), ("thickness", "Thickness", "0.5 mg/cm2 or 1 um"),
                 ("tilt", "Tilt", "0 deg"), ("density", "Density", "19.3 g/cm3")]
BACKING_FIELDS = [("material", "Backing material", "C"), ("thickness", "Backing thickness", "20 ug/cm2")]
RUN_FIELDS = [("beam_time", "Beam time", "12 h"), ("counts_wanted", "Counts wanted", "5000")]
DETECTOR_FIELDS = [("name", "Name", "D1"), ("theta", "θ", "45 deg"), ("phi", "φ", "0 deg"),
                   ("distance", "Distance", "100 mm"), ("width", "Width", "50 mm"), ("height", "Height", "50 mm"),
                   ("strips_x", "Strips x", "16"), ("strips_y", "Strips y", "16"), ("radius", "Radius", "5 mm"),
                   ("inner_radius", "Inner radius", "9 mm"), ("outer_radius", "Outer radius", "41 mm"),
                   ("rings", "Rings", "16"), ("sectors", "Sectors", "24"), ("thickness", "Thickness", "300 um"),
                   ("dead_layer", "Dead layer", "0.5 um"), ("resolution", "Resolution (FWHM)", "20 keV"),
                   ("threshold", "Threshold", "200 keV"), ("material", "Material", "Si")]
#: The excited state of a Coulomb-excitation setup ([reaction] in a setup file); type, excite and multipolarity are
#: drop-downs.
REACTION_FIELDS = [("energy", "State energy", "1.454 MeV"), ("b_up", "B(Eλ↑)", "0.0695 e2b2 or 695 e2fm4")]
REACTION_TYPES = {"elastic": "Elastic (Rutherford) scattering", "coulex": "Coulomb excitation"}
GAMMA_FIELDS = [("name", "Name", "Ge1"), ("theta", "θ", "90 deg"), ("phi", "φ", "90 deg"),
                ("distance", "Distance", "120 mm"), ("radius", "Crystal radius", "35 mm"),
                ("resolution", "Resolution (FWHM)", "2.5 keV")]
#: Which size fields each shape uses.
SHAPE_FIELDS = {"rectangle": {"width", "height", "strips_x", "strips_y"}, "circle": {"radius"},
                "annular": {"inner_radius", "outer_radius", "rings", "sectors"}}
INT_FIELDS = {"charge_state", "counts_wanted", "strips_x", "strips_y", "rings", "sectors"}

COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf", "#8c564b", "#e377c2", "#7f7f7f",
          "#bcbd22"]


def _go():
    import plotly.graph_objects as go

    return go


def _value(field: str, text):
    """A typed field value from the text the user entered ("" removes an optional field)."""
    if text is None:
        return None
    text = str(text).strip()
    if text == "":
        return None
    if field in INT_FIELDS:
        try:
            return int(text)
        except ValueError:
            return text  # the setup checks report it
    return text


# ---------------------------------------------------------------------------------------------------------------
# Figures (Plotly), one per tab


def figure_geometry(planner: Planner):
    """Beam, target and detectors in 3D (mm)."""
    go = _go()
    g = planner.geometry()
    fig = go.Figure()
    b = g["beam"]
    fig.add_trace(go.Scatter3d(x=b[:, 0], y=b[:, 1], z=b[:, 2], mode="lines", name="beam",
                               line=dict(color="#888", dash="dash")))
    t = g["target_size_mm"]
    fig.add_trace(go.Scatter3d(x=[-t, t, t, -t, -t], y=[-t, -t, t, t, -t], z=[0] * 5, mode="lines", name="target",
                               line=dict(color="#c9a227", width=6)))
    for i, d in enumerate(g["detectors"]):
        o = d["outline"]
        fig.add_trace(go.Scatter3d(x=o[:, 0], y=o[:, 1], z=o[:, 2], mode="lines", name=d["name"],
                                   line=dict(color=COLORS[i % len(COLORS)], width=5)))
    for d in g.get("gamma_detectors", []):
        o = d["outline"]
        fig.add_trace(go.Scatter3d(x=o[:, 0], y=o[:, 1], z=o[:, 2], mode="lines", name=f"{d['name']} (γ)",
                                   line=dict(color="#7c3aed", width=4, dash="dash")))
    e = g["extent"]
    axis = dict(range=[-e, e])
    fig.update_layout(scene=dict(xaxis=dict(title="x (mm)", **axis), yaxis=dict(title="y (mm)", **axis),
                                 zaxis=dict(title="z, beam (mm)", **axis), aspectmode="cube"),
                      margin=dict(l=0, r=0, t=10, b=0), height=520, legend=dict(itemsizing="constant"))
    return fig


def figure_kinematics(planner: Planner):
    """Lab energy against lab angle for ejectiles and recoils, with each detector's coverage shaded."""
    go = _go()
    k = planner.kinematics()
    fig = go.Figure()
    for i, d in enumerate(k["detectors"]):
        lo, hi = d["theta_range"]
        fig.add_vrect(x0=lo, x1=hi, fillcolor=COLORS[i % len(COLORS)], opacity=0.12, line_width=0,
                      annotation_text=d["name"], annotation_position="top")
    for c in k["curves"]:
        th, en = np.asarray(c["theta"]), np.asarray(c["energy"])
        ok = th <= c["max_angle"] + 1e-9
        fig.add_trace(go.Scatter(x=th[ok], y=en[ok], mode="lines", name=c["label"],
                                 line=dict(dash=("dot" if c.get("inelastic") else "solid")
                                           if c["particle"] == "ejectile" else "dash")))
    fig.update_layout(xaxis=dict(title="lab angle θ (deg)", range=[0, 180]), yaxis=dict(title="energy (MeV)"),
                      margin=dict(l=50, r=10, t=30, b=40), height=440)
    return fig


def figure_strips(planner: Planner, detector: str):
    """Counts per second in each strip pair or ring–sector of one detector."""
    go = _go()
    strips = planner.rates()["strips"][detector]
    ni = max(s[0] for s in strips) + 1
    nj = max(s[1] for s in strips) + 1
    z = np.zeros((nj, ni))
    for (i, j), r in strips.items():
        z[j, i] = r
    fig = go.Figure(go.Heatmap(z=z, colorscale="Viridis", colorbar=dict(title="1/s")))
    fig.update_layout(xaxis=dict(title="strip / ring"), yaxis=dict(title="strip / sector"),
                      margin=dict(l=50, r=10, t=10, b=40), height=360)
    return fig


def figure_energy_loss(planner: Planner):
    """Beam energy through the target and backing."""
    go = _go()
    el = planner.energy_loss()
    fig = go.Figure(go.Scatter(x=el["depth_mg_cm2"], y=el["beam_energy_mev"], mode="lines", name="beam"))
    fig.update_layout(xaxis=dict(title="depth (mg/cm²)"), yaxis=dict(title="beam energy (MeV)"),
                      margin=dict(l=50, r=10, t=10, b=40), height=320)
    return fig


def figure_spectra(planner: Planner, events: int = 100_000, seed: int = 1):
    """Simulated measured-energy spectra, counts per bin in the planned beam time (log scale)."""
    go = _go()
    s = planner.spectra(events=events, seed=seed)["spectra"]
    fig = go.Figure()
    for i, (name, h) in enumerate(s.items()):
        e = h["edges"]
        x = np.repeat(e, 2)[1:-1]
        y = np.repeat(h["counts"], 2)
        fig.add_trace(go.Scatter(x=x, y=np.where(y > 0, y, np.nan), mode="lines", name=name,
                                 line=dict(color=COLORS[i % len(COLORS)], width=1.2)))
    fig.update_layout(xaxis=dict(title="measured energy (MeV)"), yaxis=dict(title="counts / bin", type="log"),
                      margin=dict(l=50, r=10, t=10, b=40), height=420)
    return fig


def figure_trajectories(planner: Planner):
    """Coulomb orbits in the CM frame, with the interaction radius drawn as a circle."""
    go = _go()
    t = planner.trajectories()
    fig = go.Figure()
    for o in t["orbits"]:
        xy = o["xy"]
        fig.add_trace(go.Scatter(x=xy[:, 0], y=xy[:, 1], mode="lines",
                                 name=f"b = {o['b_fm']:.1f} fm → {o['deflection_deg']:.0f}°"))
    a = np.linspace(0, 2 * math.pi, 200)
    r = t["interaction_radius_fm"]
    fig.add_trace(go.Scatter(x=r * np.cos(a), y=r * np.sin(a), mode="lines", name="nuclear range",
                             line=dict(color="#888", dash="dot")))
    lim = 12 * t["d0_fm"]
    fig.update_layout(xaxis=dict(title="x (fm)", range=[-lim, lim]),
                      yaxis=dict(title="y (fm)", range=[-lim, lim], scaleanchor="x"),
                      margin=dict(l=50, r=10, t=10, b=40), height=480)
    return fig


def figure_sweep(sweep: dict):
    go = _go()
    fig = go.Figure(go.Scatter(x=sweep["x"], y=sweep["y"], mode="lines+markers", name=sweep["quantity"]))
    fig.update_layout(xaxis=dict(title=sweep["parameter"]), yaxis=dict(title=f"{sweep['quantity']} "
                                                                             f"({sweep['detector']})"),
                      margin=dict(l=60, r=10, t=10, b=40), height=320)
    return fig


def figure_excitation(planner: Planner):
    """Coulomb-excitation probability against CM angle (when the setup has an excited state)."""
    go = _go()
    g = planner.gamma()
    fig = go.Figure()
    if g["available"]:
        fig.add_trace(go.Scatter(x=g["theta_cm"], y=g["probability"], mode="lines", name="P(θ)"))
        if g["max_safe_angle"] < 180:
            fig.add_vrect(x0=g["max_safe_angle"], x1=180, fillcolor="#d62728", opacity=0.1, line_width=0,
                          annotation_text="not safe", annotation_position="top")
    fig.update_layout(xaxis=dict(title="CM angle θ (deg)", range=[0, 180]),
                      yaxis=dict(title="excitation probability", type="log"), margin=dict(l=60, r=10, t=10, b=40),
                      height=320)
    return fig


FIGURES = {"geometry": figure_geometry, "kinematics": figure_kinematics, "energy_loss": figure_energy_loss,
           "spectra": figure_spectra, "trajectories": figure_trajectories, "gamma": figure_excitation}


def report_zip(planner: Planner, seed: int = 1, events: int = 200_000) -> bytes:
    """The report folder (HTML, CSV, figures, setup, ROOT file if uproot is installed) as a zip archive."""
    import tempfile
    from pathlib import Path

    from .report import build

    with tempfile.TemporaryDirectory() as tmp:
        build(planner.experiment, seed=seed, events=events).write(tmp)
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for p in sorted(Path(tmp).rglob("*")):
                if p.is_file():
                    z.write(p, p.relative_to(tmp).as_posix())
        return buf.getvalue()


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


def build_page(example: str = "alpha_on_gold", events: int = 100_000, mode: Optional[str] = None) -> None:
    """Build the planner page for the current client (call inside a NiceGUI page function).

    ``mode`` is "guided" (a step-by-step workflow) or "expert" (every input and result at once); ``None`` takes the
    one this browser used last, guided the first time."""
    from nicegui import run, ui

    from .experiment import Experiment, SetupError

    state = {"planner": Planner.example(example), "events": events, "mode": mode or "guided",
             "step": guide.STEPS[0].key, "example": example}

    def P() -> Planner:  # noqa: N802
        return state["planner"]

    # -- actions ----------------------------------------------------------------------------------------------
    def changed(ok: bool) -> None:
        step_status.refresh()
        if ok:
            refresh_results()
        else:
            ui.notify("; ".join(P().problems), type="negative", multi_line=True)
            warnings_banner.refresh()

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
        changed(P().set(section, field, value))
        if (section.startswith(("detector", "gamma detector")) and field in ("name", "shape")) or (
                section == "reaction" and field == "type"):
            setup_panel.refresh()

    def add_detector() -> None:
        n = len(P().draft["detectors"])
        changed(P().add_detector(name=f"D{n + 1}", shape="circle", theta="45 deg", distance="100 mm",
                                 radius="5 mm", thickness="300 um"))
        setup_panel.refresh()

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

    def add_gamma_detector() -> None:
        n = len(P().draft.get("gamma_detectors", []))
        changed(P().add_gamma_detector(name=f"Ge{n + 1}", theta="90 deg", phi="90 deg", distance="120 mm",
                                       radius="35 mm", resolution="2.5 keV"))
        setup_panel.refresh()

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

    async def load_example(name: str, then: Optional[str] = None) -> None:
        note = ui.notification(f"Loading {name}…", spinner=True, timeout=None)
        try:
            planner = await run.io_bound(_prepared, name)
        finally:
            note.dismiss()
        state["planner"] = planner
        state["example"] = name
        example_select.value = name  # its handler sees the name is already loaded
        if then:
            state["step"] = then
        setup_panel.refresh()
        refresh_results()
        main_area.refresh()

    async def load_file(e) -> None:
        upload_dialog.close()
        try:
            state["planner"] = Planner(Experiment.from_toml(await e.file.text()))
        except (SetupError, ValueError) as err:
            ui.notify(f"Could not read the setup: {err}", type="negative", multi_line=True)
            return
        setup_panel.refresh()
        refresh_results()

    def save_setup() -> None:
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
            icon = ui.icon("help_outline", size="xs").classes("cursor-help text-slate-400")
            with icon, ui.tooltip().classes("bg-slate-800 text-white max-w-xs p-2"):
                ui.label(h.what).classes("text-sm")
                ui.label(f"Typical: {h.typical}").classes("text-xs text-slate-300 mt-1")
                if h.effect != "—":
                    ui.label(f"Raise it: {h.effect}").classes("text-xs text-slate-300")
            icon.on("click.stop", lambda: ui.notify(h.text(), multi_line=True, close_button=True,
                                                    classes="whitespace-pre-line"))

    def section(title, fields, sec_name, values):
        if title:
            ui.label(title).classes("text-sm font-semibold mt-3")
        # The guided steps have few fields each: one per line leaves room for the labels.
        with ui.grid(columns=1 if state["mode"] == "guided" else 2).classes("w-full gap-1"):
            for field, label, placeholder in fields:
                v = values.get(field)
                inp = ui.input(label, value="" if v is None else str(v), placeholder=placeholder).props(
                    "dense outlined").classes("w-full")
                bind(inp, sec_name, field)
                help_icon(inp, sec_name, field)

    @ui.refreshable
    def setup_panel():
        if state["mode"] == "guided":
            stepper()
            return
        d = P().draft
        title = ui.input("Title", value=d.get("title", "")).props("dense outlined").classes("w-full")
        bind(title, "title", "")
        help_icon(title, "title", "")
        section("Beam", BEAM_FIELDS, "beam", d["beam"])
        target_inputs(d)
        section("Run", RUN_FIELDS, "run", d["run"])
        reaction_section(d.get("reaction", {"type": "elastic"}))
        particle_detectors(d)
        if d.get("reaction", {}).get("type") == "coulex":
            gamma_detectors(d)

    def target_inputs(d: dict) -> None:
        section("Target", TARGET_FIELDS, "target", d["target"])
        section("Backing (optional)", BACKING_FIELDS, "backing", d["target"].get("backing", {}))

    def particle_detectors(d: dict) -> None:
        with ui.row().classes("items-center mt-3 w-full"):
            ui.label("Particle detectors").classes("text-sm font-semibold")
            ui.space()
            ui.button("Add", icon="add", on_click=add_detector).props("dense flat")
        ui.label("Silicon detectors: they measure the energy of the scattered beam particles and recoils.").classes(
            "text-xs text-slate-500")
        size_fields = set().union(*SHAPE_FIELDS.values())
        for i, det in enumerate(d["detectors"]):
            shape = det.get("shape", "circle")
            with ui.expansion(f"{det.get('name') or f'D{i + 1}'} · {shape} silicon").classes("w-full bg-white"):
                sel = ui.select(list(SHAPE_FIELDS), value=shape, label="Shape",
                                on_change=lambda e, k=i: edit(f"detector {k + 1}", "shape", e.value)).props(
                    "dense outlined").classes("w-full")
                help_icon(sel, "detector", "shape")
                wanted = [f for f in DETECTOR_FIELDS if f[0] not in size_fields or f[0] in SHAPE_FIELDS[shape]]
                section("", wanted, f"detector {i + 1}", det)
                with ui.row():
                    ui.button("Duplicate", icon="content_copy", on_click=duplicate_detector(i)).props("dense flat")
                    ui.button("Remove", icon="delete", on_click=remove_detector(i)).props(
                        "dense flat color=negative")

    def gamma_detectors(d: dict) -> None:
        with ui.row().classes("items-center mt-3 w-full"):
            ui.label("γ-ray detectors").classes("text-sm font-semibold")
            ui.space()
            ui.button("Add", icon="add", on_click=add_gamma_detector).props("dense flat")
        ui.label("Germanium detectors: they see the γ ray emitted when the excited state decays.").classes(
            "text-xs text-slate-500")
        for i, gd in enumerate(d.get("gamma_detectors", [])):
            with ui.expansion(f"{gd.get('name') or f'γ{i + 1}'} · germanium").classes("w-full bg-white"):
                section("", GAMMA_FIELDS, f"gamma detector {i + 1}", gd)
                with ui.row():
                    ui.button("Duplicate", icon="content_copy",
                              on_click=duplicate_gamma_detector(i)).props("dense flat")
                    ui.button("Remove", icon="delete", on_click=remove_gamma_detector(i)).props(
                        "dense flat color=negative")

    def reaction_section(reaction: dict, title: str = "Reaction") -> None:
        if title:
            ui.label(title).classes("text-sm font-semibold mt-3")
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
        section("", REACTION_FIELDS, "reaction", reaction)

    # -- guided workflow --------------------------------------------------------------------------------------
    def go_to(key: str):
        def run_() -> None:
            state["step"] = key
            setup_panel.refresh()
            main_area.refresh()
        return run_

    def step_inputs(key: str) -> None:
        d = P().draft
        if key == "goal":
            for goal, (name, text, ex) in guide.GOALS.items():
                with ui.card().classes("w-full p-3 gap-1 bg-white"):
                    ui.label(name).classes("font-semibold")
                    ui.label(text).classes("text-sm text-slate-600")
                    ui.button("Start here", icon="arrow_forward",
                              on_click=lambda ex=ex: load_example(ex, then="beam")).props(
                        "dense flat no-caps").classes("self-start")
                    ui.label(f"Starts from the example {ex}.").classes("text-xs text-slate-500")
            ui.label("Or keep the current setup and set the reaction here:").classes("text-sm mt-2")
            reaction_section(d.get("reaction", {"type": "elastic"}), title="")
        elif key == "beam":
            section("", BEAM_FIELDS, "beam", d["beam"])
        elif key == "target":
            target_inputs(d)
        elif key == "detectors":
            particle_detectors(d)
        elif key == "gamma":
            gamma_detectors(d)
        elif key == "rates":
            section("", RUN_FIELDS, "run", d["run"])

    def stepper() -> None:
        steps = guide.steps_for(P())
        keys = [s["step"].key for s in steps]
        if state["step"] not in keys:
            state["step"] = keys[0]

        def picked(e) -> None:
            if e.value != state["step"]:
                state["step"] = e.value
                main_area.refresh()
                step_status.refresh()

        with ui.stepper(value=state["step"], on_value_change=picked).props(
                "vertical flat header-nav animated=false").classes("w-full bg-transparent"):
            for i, s in enumerate(steps):
                step = s["step"]
                with ui.step(step.key, title=f"{i + 1}. {step.title}"):
                    ui.label(step.intro).classes("text-sm text-slate-600")
                    step_inputs(step.key)
                    step_status(step.key, keys[i - 1] if i else None, keys[i + 1] if i + 1 < len(keys) else None)

    @ui.refreshable
    def step_status(key: str = "", back: Optional[str] = None, nxt: Optional[str] = None) -> None:
        if not key:
            return
        problems = next((s["problems"] for s in guide.steps_for(P()) if s["step"].key == key), [])
        for p in problems:
            with ui.row().classes("items-start no-wrap gap-1 mt-1"):
                ui.icon("error").classes("text-red-700")
                ui.label(p).classes("text-sm text-red-700")
        with ui.row().classes("mt-2 gap-2"):
            if nxt:
                btn = ui.button("Next", icon="arrow_downward", on_click=go_to(nxt))
                if problems:
                    btn.disable()
                    btn.tooltip("Fix the problem above first")
            if back:
                ui.button("Back", on_click=go_to(back)).props("flat")

    # -- warnings banner --------------------------------------------------------------------------------------
    @ui.refreshable
    def warnings_banner():
        ws = P().warnings()
        if not ws:
            ui.label("No warnings.").classes("text-green-700")
            return
        style = {"error": "text-red-700 font-semibold", "warning": "text-amber-800", "note": "text-slate-600"}
        icon = {"error": "error", "warning": "warning", "note": "info"}
        with ui.card().classes("w-full bg-amber-50 p-2 gap-1"):
            for w in ws:
                with ui.row().classes("items-start no-wrap gap-2"):
                    ui.icon(icon[w.level]).classes(style[w.level])
                    ui.label(w.text).classes(style[w.level] + " text-sm")

    # -- result tabs ------------------------------------------------------------------------------------------
    def columns(spec):
        return [{"name": k, "label": lab, "field": k, "align": "left"} for k, lab in spec]

    @ui.refreshable
    def geometry_panel():
        ui.plotly(figure_geometry(P())).classes("w-full")
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
        ui.plotly(figure_kinematics(P())).classes("w-full")

    @ui.refreshable
    def rates_panel():
        r = P().rates()
        rows = [{"detector": x["detector"], "omega": _fmt(x["solid_angle_msr"]),
                 "theta": f"{x['theta_range'][0]:.1f}–{x['theta_range'][1]:.1f}", "rate": _fmt(x["rate_per_s"]),
                 "counts": _fmt(x["counts_in_run"]),
                 "error": f"{100 * x['relative_error']:.2g}%" if math.isfinite(x["relative_error"]) else "—",
                 "time": _time(x["beam_time_s"])} for x in r["rows"]]
        wanted = f" for {r['counts_wanted']} counts" if r["counts_wanted"] else ""
        ui.table(columns=columns((("detector", "Detector"), ("omega", "Ω (msr)"), ("theta", "θ (deg)"),
                                  ("rate", "Rate (1/s)"), ("counts", "Counts in run"), ("error", "Stat. error"),
                                  ("time", f"Beam time{wanted}"))), rows=rows).props("dense flat")
        names = [x["detector"] for x in r["rows"]]

        @ui.refreshable
        def strips(name):
            ui.plotly(figure_strips(P(), name)).classes("w-full")

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
                        ui.plotly(figure_sweep(s)).classes("w-full")

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
        ui.plotly(figure_energy_loss(P())).classes("w-full")
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
                    ui.plotly(fig).classes("w-full")

            ui.button("Simulate", icon="play_arrow", on_click=simulate)
        plot_area = ui.column().classes("w-full")
        with plot_area:
            ui.plotly(figure_spectra(P(), state["events"], 1)).classes("w-full")

    @ui.refreshable
    def trajectories_panel():
        t = P().trajectories()
        ui.label(f"{P().experiment.beam.nuclide} on {t['target']}: head-on distance d₀ = {t['d0_fm']:.2f} fm, "
                 f"nuclear range {t['interaction_radius_fm']:.1f} fm, grazing angle "
                 f"{t['grazing_angle_deg']:.1f}° (CM).").classes("text-sm")
        ui.plotly(figure_trajectories(P())).classes("w-full")

    @ui.refreshable
    def gamma_panel():
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
        ui.plotly(figure_excitation(P())).classes("w-full")
        ui.table(columns=columns((("detector", "Particle detector"), ("rate", "Excitation events (1/s)"))),
                 rows=[{"detector": k, "rate": _fmt(v)} for k, v in g["rates"].items()]).props("dense flat")
        if g["particles"]:
            ui.label("Particle energies: elastic and after exciting the state").classes("font-semibold mt-2")
            ui.label("At each particle detector's smallest, central and largest angle, at the reaction point "
                     "(before energy loss in the target; the Spectra tab includes it). β is the speed of the "
                     "excited nucleus in that event, which the Doppler correction needs.").classes(
                "text-xs text-slate-500")
            ui.table(columns=columns((("detector", "Detector"), ("particle", "Particle"), ("theta", "θ lab (deg)"),
                                      ("el", "Elastic (MeV)"), ("ex", "Excited (MeV)"), ("diff", "Difference (MeV)"),
                                      ("beta", "β excited"))),
                     rows=[{"detector": r["detector"], "particle": f"{r['nuclide']} ({r['particle']})",
                            "theta": f"{r['theta_lab']:.1f}", "el": f"{r['elastic_mev']:.2f}",
                            "ex": f"{r['excited_mev']:.2f}", "diff": f"{r['difference_mev']:.3f}",
                            "beta": f"{r['beta_excited']:.4f}"} for r in g["particles"]]).props("dense flat")
        if g["doppler"]:
            ui.label("γ rays: Doppler-shifted energy and width").classes("font-semibold mt-2")
            ui.table(columns=columns((("p", "Particle detector"), ("g", "γ detector"), ("mean", "E_γ (keV)"),
                                      ("shift", "Shift (keV)"), ("fwhm", "FWHM (keV)"))),
                     rows=[{"p": r["particle_detector"], "g": r["gamma_detector"], "mean": f"{r['mean_kev']:.2f}",
                            "shift": f"{r['shift_kev']:+.2f}", "fwhm": f"{r['fwhm_kev']:.2f}"}
                           for r in g["doppler"]]).props("dense flat")
        else:
            ui.label("Add γ-ray detectors in the setup panel for Doppler shifts.").classes("text-sm")

    @ui.refreshable
    def report_panel():
        ui.markdown("The report collects the setup, every warning, the detector table, kinematics, peaks, energy "
                    "loss, spectra and each model's validation status. The zip holds `report.html` (print it to "
                    "PDF from the browser), CSV tables, figures (PNG and PDF), the setup file and, with uproot "
                    "installed, `events.root`.")

        async def download():
            note = ui.notification("Building the report…", spinner=True, timeout=None)
            try:
                data = await run.io_bound(report_zip, P(), 1, 200_000)
            except ImportError as err:
                ui.notify(f"The report export is not available in this version: {err}", type="warning")
                return
            finally:
                note.dismiss()
            ui.download.content(data, "physim-report.zip")

        ui.button("Build and download the report", icon="description", on_click=download)

    panels = {"geometry": geometry_panel, "kinematics": kinematics_panel, "rates": rates_panel,
              "energy_loss": energy_loss_panel, "spectra": spectra_panel, "trajectories": trajectories_panel,
              "gamma": gamma_panel, "report": report_panel}

    def refresh_results() -> None:
        warnings_banner.refresh()
        reading_box.refresh()
        for p in panels.values():
            p.refresh()

    @ui.refreshable
    def reading_box(tab: str) -> None:
        try:
            lines = guide.reading(P(), tab)
        except Exception:  # noqa: BLE001 -- a reading must never break the page
            lines = []
        if not lines:
            return
        with ui.row().classes("w-full items-start no-wrap gap-2 bg-sky-50 rounded p-2"):
            ui.icon("lightbulb").classes("text-sky-700 mt-0.5")
            with ui.column().classes("gap-1"):
                ui.label("How to read this").classes("text-xs font-semibold text-sky-800 uppercase")
                ui.label(" ".join(lines)).classes("text-sm text-slate-800")

    def explain(tab: str):
        e = P().explain(tab)
        with ui.expansion("Explain: " + e["title"], icon="school").classes("w-full bg-slate-50 mt-2"):
            ui.markdown(f"**Formula.** {e['formula']}\n\n**Assumptions.** {e['assumptions']}\n\n"
                        f"**Where it stops being valid.** {e['limits']}\n\n"
                        f"**Validation:** see the physics register page `{e['register']}` in physim's docs.")

    labels = {"geometry": "Geometry", "kinematics": "Kinematics", "rates": "Rates and beam time",
              "energy_loss": "Energy loss", "spectra": "Spectra", "trajectories": "Trajectories",
              "gamma": "Excitation and γ rays", "report": "Report"}

    @ui.refreshable
    def main_area() -> None:
        warnings_banner()
        if state["mode"] == "guided":
            step = next(s for s in guide.STEPS if s.key == state["step"])
            if not step.tabs:
                welcome()
            for t in step.tabs:
                ui.label(labels[t]).classes("text-lg font-semibold mt-2")
                reading_box(t)
                panels[t]()
                explain(t)
            return
        with ui.tabs().classes("w-full") as tabs:
            tab = {t: ui.tab(labels[t]) for t in TABS}
        with ui.tab_panels(tabs, value=tab["geometry"]).classes("w-full"):
            for t in TABS:
                with ui.tab_panel(tab[t]):
                    reading_box(t)
                    panels[t]()
                    explain(t)

    def welcome() -> None:
        with ui.column().classes("max-w-3xl gap-2 mt-2"):
            ui.label("Plan an experiment, step by step").classes("text-xl font-semibold")
            ui.markdown(
                "The steps on the left take you through a plan in the order you would decide it: **what to "
                "measure**, the **beam**, the **target**, the **detectors**, then the **rates and beam time**, the "
                "**spectra** you will see, and the **report**.\n\n"
                "- Every step starts filled in with values that work, so the results on this side are always "
                "complete. Change one value at a time and watch what it does.\n"
                "- Hover over (or tap) the **?** in any field for what it means, a typical value, and what raising "
                "it does.\n"
                "- Each result opens with **How to read this**: what to look for, with your numbers. The "
                "**Explain** panel underneath has the formula and its limits.\n"
                "- Problems show in red inside the step; **Next** waits until they are fixed.\n\n"
                "**Expert view** (top right) shows every input and result at once.")

    def set_mode(mode: str) -> None:
        if mode == state["mode"]:
            return
        state["mode"] = mode
        ui.run_javascript(f"try {{ localStorage.setItem('physim-planner-mode', '{mode}') }} catch (e) {{}}")
        setup_panel.refresh()
        main_area.refresh()

    async def restore_mode() -> None:
        if mode is not None:  # given in the address (?mode=...)
            return
        try:
            saved = await ui.run_javascript("localStorage.getItem('physim-planner-mode')", timeout=3)
        except Exception:  # noqa: BLE001 -- no answer from the browser: keep the default
            return
        if saved in ("guided", "expert") and saved != state["mode"]:
            mode_toggle.value = saved

    # -- layout -----------------------------------------------------------------------------------------------
    ui.query("body").style("background-color: #ffffff; color: #1b1b1b")  # light page in a dark-mode browser too
    with ui.header().classes("items-center bg-slate-800 py-1"):
        ui.label("physim · experiment planner").classes("text-lg font-medium")
        ui.space()
        mode_toggle = ui.toggle({"guided": "Guided", "expert": "Expert view"}, value=state["mode"],
                                on_change=lambda e: set_mode(e.value)).props(
            "dense no-caps toggle-color=white toggle-text-color=slate-800 text-color=white")
        example_select = ui.select(
            Planner.examples(), value=example, label="Start from example",
            on_change=lambda e: None if e.value == state["example"] else load_example(e.value)).props(
            "dense dark outlined").classes("w-60")
        ui.button("Load setup", icon="upload", on_click=lambda: upload_dialog.open()).props("flat color=white")
        ui.button("Save setup", icon="download", on_click=save_setup).props("flat color=white")
    with ui.dialog() as upload_dialog, ui.card():
        ui.label("Load a setup file (.toml)")
        ui.upload(auto_upload=True, on_upload=load_file).props("accept=.toml max-files=1")
    with ui.left_drawer(value=True).classes("bg-slate-50").props("width=420 bordered behavior=desktop"):
        setup_panel()
    with ui.column().classes("w-full gap-2"):
        main_area()
    ui.timer(0.2, restore_mode, once=True)


# -- starting the server --------------------------------------------------------------------------------------------

#: What ``/physim-planner`` answers, so a second launch can recognise a planner already running on the port.
MARKER = "physim experiment planner"


def log_path():
    """Where ``physim app --desktop`` writes its log (there is no console to write to)."""
    from pathlib import Path

    base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), ".local", "state")
    return Path(base) / "physim" / "planner.log"


def planner_at(port: int, host: str = "127.0.0.1") -> bool:
    """Whether a physim planner is already serving at ``host:port``."""
    import json
    import urllib.request

    try:
        with urllib.request.urlopen(f"http://{host}:{port}/physim-planner", timeout=2) as r:
            return json.loads(r.read().decode()).get("app") == MARKER
    except (OSError, ValueError):
        return False


def port_free(port: int, host: str = "127.0.0.1") -> bool:
    import socket

    with socket.socket() as s:
        try:
            s.bind((host, port))
        except OSError:
            return False
    return True


def pick_port(port: int, host: str = "127.0.0.1") -> int:
    """``port`` if it is free, otherwise a free port chosen by the system."""
    import socket

    if port_free(port, host):
        return port
    with socket.socket() as s:
        s.bind((host, 0))
        return s.getsockname()[1]


def main(argv: Optional[list] = None) -> None:
    """``physim app``: start the planner in the browser."""
    import argparse

    ap = argparse.ArgumentParser(prog="physim app", description="Start the experiment planner in the browser.")
    ap.add_argument("--example", default="alpha_on_gold", help="example setup to start from")
    ap.add_argument("--port", type=int, default=8080)
    ap.add_argument("--host", default="127.0.0.1",
                    help="address to listen on (default: this computer only; 0.0.0.0 serves the local network)")
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    ap.add_argument("--desktop", action="store_true",
                    help="started from a shortcut: log to a file, reuse a planner that is already running, use "
                         "another port if this one is taken, and stop when the last browser tab has been closed")
    ap.add_argument("--idle-exit", type=float, default=None, metavar="SECONDS",
                    help="stop when no browser tab has been open for this long (default with --desktop: 60)")
    args = ap.parse_args(argv)
    idle_exit = args.idle_exit if args.idle_exit is not None else (60.0 if args.desktop else None)
    local = "127.0.0.1" if args.host in ("0.0.0.0", "::") else args.host

    if args.desktop:
        import sys

        path = log_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        # Started with pythonw.exe there is no console: sys.stdout is None, and printing would fail.
        log = open(path, "a", encoding="utf-8", buffering=1)  # noqa: SIM115 -- open for the life of the server
        sys.stdout = sys.stderr = log
        print(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')} physim app {' '.join(argv or [])}")
        if planner_at(args.port, local):
            print(f"a planner is already running on port {args.port}; opening it")
            if not args.no_browser:
                import webbrowser

                webbrowser.open(f"http://{local}:{args.port}/")
            return
        args.port = pick_port(args.port, args.host)

    try:
        import matplotlib

        matplotlib.use("Agg")  # the report draws its figures in a worker thread, without a screen
    except ImportError:
        pass
    try:
        from nicegui import app, ui
    except ImportError:
        raise SystemExit("The planner app needs NiceGUI: pip install physim-engine[app]") from None

    @ui.page("/")
    def index(example: str = args.example, mode: Optional[str] = None):
        build_page(example, mode=mode if mode in ("guided", "expert") else None)

    @app.get("/physim-planner")
    def marker():
        from .. import __version__

        return {"app": MARKER, "version": __version__}

    if idle_exit is not None:
        from nicegui import Client

        idle = {"since": time.monotonic()}

        def check_idle():
            if any(c.has_socket_connection for c in list(Client.instances.values())):
                idle["since"] = time.monotonic()
            elif time.monotonic() - idle["since"] > idle_exit:
                print(f"no browser tab open for {idle_exit:g} s; stopping")
                app.shutdown()

        def start_clock():
            idle["since"] = time.monotonic()
            _start_idle_timer(check_idle, min(5.0, idle_exit / 4))

        app.on_startup(start_clock)

    ui.run(title="physim experiment planner", host=args.host, port=args.port, show=not args.no_browser,
           reload=False, show_welcome_message=True)


def _start_idle_timer(check, every: float) -> None:
    import asyncio

    async def loop():
        while True:
            await asyncio.sleep(every)
            check()

    asyncio.get_running_loop().create_task(loop())


__all__ = ["FIGURES", "MARKER", "build_page", "log_path", "main", "pick_port", "planner_at", "figure_energy_loss", "figure_geometry", "figure_kinematics",
           "figure_excitation", "figure_spectra", "figure_strips", "figure_sweep", "figure_trajectories",
           "report_zip"]
