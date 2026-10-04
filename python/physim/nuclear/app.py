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
import zipfile
from typing import Optional

import numpy as np

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


def build_page(example: str = "alpha_on_gold", events: int = 100_000) -> None:
    """Build the planner page for the current client (call inside a NiceGUI page function)."""
    from nicegui import run, ui

    from .experiment import Experiment, SetupError

    state = {"planner": Planner.example(example), "events": events}

    def P() -> Planner:  # noqa: N802
        return state["planner"]

    # -- actions ----------------------------------------------------------------------------------------------
    def changed(ok: bool) -> None:
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
        elif section.startswith("detector"):
            old = d["detectors"][int(section.split()[1]) - 1].get(field)
        else:
            old = d.get(section, {}).get(field)
        if (None if old is None else str(old)) == (None if value is None else str(value)):
            return
        changed(P().set(section, field, value))
        if section.startswith("detector") and field in ("name", "shape"):
            setup_panel.refresh()

    def add_detector() -> None:
        n = len(P().draft["detectors"])
        changed(P().add_detector(name=f"D{n + 1}", shape="circle", theta="45 deg", distance="100 mm",
                                 radius="5 mm", thickness="300 um"))
        setup_panel.refresh()

    def duplicate_detector(k: int) -> None:
        changed(P().duplicate_detector(k))
        setup_panel.refresh()

    def remove_detector(k: int) -> None:
        if len(P().draft["detectors"]) == 1:
            ui.notify("A setup needs at least one detector.", type="warning")
            return
        changed(P().remove_detector(k))
        setup_panel.refresh()

    def load_example(name: str) -> None:
        state["planner"] = Planner.example(name)
        setup_panel.refresh()
        refresh_results()

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

    def section(title, fields, sec_name, values):
        if title:
            ui.label(title).classes("text-sm font-semibold mt-3")
        with ui.grid(columns=2).classes("w-full gap-1"):
            for field, label, placeholder in fields:
                v = values.get(field)
                inp = ui.input(label, value="" if v is None else str(v), placeholder=placeholder).props(
                    "dense outlined").classes("w-full")
                bind(inp, sec_name, field)

    @ui.refreshable
    def setup_panel():
        d = P().draft
        title = ui.input("Title", value=d.get("title", "")).props("dense outlined").classes("w-full")
        bind(title, "title", "")
        section("Beam", BEAM_FIELDS, "beam", d["beam"])
        section("Target", TARGET_FIELDS, "target", d["target"])
        section("Backing (optional)", BACKING_FIELDS, "backing", d["target"].get("backing", {}))
        section("Run", RUN_FIELDS, "run", d["run"])
        with ui.row().classes("items-center mt-3 w-full"):
            ui.label("Detectors").classes("text-sm font-semibold")
            ui.space()
            ui.button("Add", icon="add", on_click=add_detector).props("dense flat")
        size_fields = set().union(*SHAPE_FIELDS.values())
        for i, det in enumerate(d["detectors"]):
            shape = det.get("shape", "circle")
            with ui.expansion(f"{det.get('name') or f'D{i + 1}'} · {shape}").classes("w-full bg-white"):
                ui.select(list(SHAPE_FIELDS), value=shape, label="Shape",
                          on_change=lambda e, k=i: edit(f"detector {k + 1}", "shape", e.value)).props(
                    "dense outlined").classes("w-full")
                wanted = [f for f in DETECTOR_FIELDS if f[0] not in size_fields or f[0] in SHAPE_FIELDS[shape]]
                section("", wanted, f"detector {i + 1}", det)
                with ui.row():
                    ui.button("Duplicate", icon="content_copy",
                              on_click=lambda k=i: duplicate_detector(k)).props("dense flat")
                    ui.button("Remove", icon="delete",
                              on_click=lambda k=i: remove_detector(k)).props("dense flat color=negative")

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
        ui.table(columns=columns((("detector", "Detector"), ("omega", "Ω (msr)"), ("theta", "θ (deg)"),
                                  ("phi", "φ (deg)"), ("segments", "Segments"))), rows=rows).props("dense flat")

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
            ui.markdown(g["reason"] + " To plan Coulomb excitation, set the reaction to `coulex` in the setup "
                        "file (see the guide) and add `[[gamma_detectors]]`.")
            return
        s = g["state"]
        ui.label(f"{s['excite'].capitalize()} excited to {s['energy_kev']:g} keV ({s['multipolarity']}, "
                 f"B↑ = {s['b_up_e2fm']:.4g} e²fm^{2 * int(s['multipolarity'][1])}): ξ = {g['xi']:.2f}, "
                 f"η = {g['eta']:.1f}, total {g['total_mb']:.3g} mb; safe up to "
                 f"{g['max_safe_angle']:.0f}° CM.").classes("text-sm")
        ui.plotly(figure_excitation(P())).classes("w-full")
        ui.table(columns=columns((("detector", "Detector"), ("rate", "Excitation events (1/s)"))),
                 rows=[{"detector": k, "rate": _fmt(v)} for k, v in g["rates"].items()]).props("dense flat")
        if g["doppler"]:
            ui.label("γ rays: Doppler-shifted energy and width").classes("font-semibold mt-2")
            ui.table(columns=columns((("p", "Particle detector"), ("g", "γ detector"), ("mean", "E_γ (keV)"),
                                      ("shift", "Shift (keV)"), ("fwhm", "FWHM (keV)"))),
                     rows=[{"p": r["particle_detector"], "g": r["gamma_detector"], "mean": f"{r['mean_kev']:.2f}",
                            "shift": f"{r['shift_kev']:+.2f}", "fwhm": f"{r['fwhm_kev']:.2f}"}
                           for r in g["doppler"]]).props("dense flat")
        else:
            ui.label("Add [[gamma_detectors]] to the setup file for Doppler shifts.").classes("text-sm")

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
        for p in panels.values():
            p.refresh()

    def explain(tab: str):
        e = P().explain(tab)
        with ui.expansion("Explain: " + e["title"], icon="school").classes("w-full bg-slate-50 mt-2"):
            ui.markdown(f"**Formula.** {e['formula']}\n\n**Assumptions.** {e['assumptions']}\n\n"
                        f"**Where it stops being valid.** {e['limits']}\n\n"
                        f"**Validation:** see the physics register page `{e['register']}` in physim's docs.")

    labels = {"geometry": "Geometry", "kinematics": "Kinematics", "rates": "Rates and beam time",
              "energy_loss": "Energy loss", "spectra": "Spectra", "trajectories": "Trajectories",
              "gamma": "Excitation and γ rays", "report": "Report"}

    # -- layout -----------------------------------------------------------------------------------------------
    ui.query("body").style("background-color: #ffffff; color: #1b1b1b")  # light page in a dark-mode browser too
    with ui.header().classes("items-center bg-slate-800 py-1"):
        ui.label("physim · experiment planner").classes("text-lg font-medium")
        ui.space()
        ui.select(Planner.examples(), value=example, label="Start from example",
                  on_change=lambda e: load_example(e.value)).props("dense dark outlined").classes("w-60")
        ui.button("Load setup", icon="upload", on_click=lambda: upload_dialog.open()).props("flat color=white")
        ui.button("Save setup", icon="download", on_click=save_setup).props("flat color=white")
    with ui.dialog() as upload_dialog, ui.card():
        ui.label("Load a setup file (.toml)")
        ui.upload(auto_upload=True, on_upload=load_file).props("accept=.toml max-files=1")
    with ui.left_drawer(value=True).classes("bg-slate-50").props("width=400 bordered behavior=desktop"):
        setup_panel()
    with ui.column().classes("w-full gap-2"):
        warnings_banner()
        with ui.tabs().classes("w-full") as tabs:
            tab = {t: ui.tab(labels[t]) for t in TABS}
        with ui.tab_panels(tabs, value=tab["geometry"]).classes("w-full"):
            for t in TABS:
                with ui.tab_panel(tab[t]):
                    panels[t]()
                    explain(t)


def main(argv: Optional[list] = None) -> None:
    """``physim app``: start the planner in the browser."""
    import argparse

    ap = argparse.ArgumentParser(prog="physim app", description="Start the experiment planner in the browser.")
    ap.add_argument("--example", default="alpha_on_gold", help="example setup to start from")
    ap.add_argument("--port", type=int, default=8080)
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    args = ap.parse_args(argv)
    try:
        from nicegui import ui
    except ImportError:
        raise SystemExit("The planner app needs NiceGUI: pip install physim-engine[app]") from None

    @ui.page("/")
    def index(example: str = args.example):
        build_page(example)

    ui.run(title="physim experiment planner", port=args.port, show=not args.no_browser, reload=False,
           show_welcome_message=True)


__all__ = ["FIGURES", "build_page", "main", "figure_energy_loss", "figure_geometry", "figure_kinematics",
           "figure_excitation", "figure_spectra", "figure_strips", "figure_sweep", "figure_trajectories",
           "report_zip"]
