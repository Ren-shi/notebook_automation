"""The experiment planner web app: build a setup by pointing and clicking, see every result, export the report.

Start it with ``physim app`` (or ``python -m physim app``); it opens in the browser at http://localhost:8080. It runs
on your computer only: no accounts, nothing uploaded. Needs the ``app`` extra: ``pip install physim-engine[app]``.

The app is a thin layer over :class:`physim.nuclear.planner.Planner`: every number it shows can be had from Python
too. The figure functions here (``figure_geometry``, ``figure_kinematics``, ...) return Plotly figures and work in
a notebook as well.
"""

from __future__ import annotations

import base64
import io
import math
import os
import time
import zipfile
from typing import Optional

import numpy as np

from . import guide
from .analysis import SHIFT_H
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
                ("resolution", "Resolution (FWHM)", "2.5 keV"), ("efficiency", "Efficiency (full peak)", "2 %"),
                ("material", "Crystal material", "Ge"), ("absorbers", "Absorbers", "Pb 1 mm, Cu 0.5 mm")]
#: Shown in place of the radius for a γ-ray detector with real crystals (a model of the catalogue).
CRYSTAL_FIELDS = [("crystals", "Crystals", "4"), ("crystal_diameter", "Crystal diameter", "50 mm"),
                  ("crystal_length", "Crystal length", "70 mm"), ("crystal_pitch", "Crystal pitch", "45 mm"),
                  ("housing_side", "Housing side", "101 mm"), ("window_gap", "Window to crystal", "5 mm")]
#: Which size fields each shape uses.
SHAPE_FIELDS = {"rectangle": {"width", "height", "strips_x", "strips_y"}, "circle": {"radius"},
                "annular": {"inner_radius", "outer_radius", "rings", "sectors"}}
INT_FIELDS = {"charge_state", "counts_wanted", "strips_x", "strips_y", "rings", "sectors", "crystals"}

#: Detector and curve colours (Okabe–Ito first): distinguishable for colour-blind readers, legible on the light and
#: the dark page.
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F", "#882255", "#44AA99",
          "#DDCC77"]

#: Figure colours of each page theme.
THEMES = {"light": {"paper": "#FFFFFF", "ink": "#16191C", "grid": "#E4E7E4"},
          "dark": {"paper": "#171B1E", "ink": "#E9EBE9", "grid": "#2A3035"}}

FONT = '"IBM Plex Sans", "Segoe UI", system-ui, -apple-system, sans-serif'
MONO = '"IBM Plex Mono", ui-monospace, "Cascadia Mono", Consolas, monospace'

#: The page's style sheet. Colours are variables, set once for the light page and once for the dark one
#: (``body--dark`` is Quasar's dark-mode class). No fonts or scripts are fetched: the app works offline.
STYLE = """
body { --ps-ground: #F1F2F0; --ps-surface: #FFFFFF; --ps-sunk: #F7F8F6; --ps-ink: #16191C; --ps-muted: #555C63;
  --ps-line: #D5D9D6; --ps-accent: #0B5FA5; --ps-ok: #1E6B45; --ps-warn: #8A4B00; --ps-bad: #B3261E;
  --ps-warnbox: #FBF3E4; --ps-note: #EAF2FA; --q-primary: #0B5FA5 !important;
  background: var(--ps-ground) !important; color: var(--ps-ink); font-family: %(font)s; }
body.body--dark { --ps-ground: #0F1214; --ps-surface: #171B1E; --ps-sunk: #1D2226; --ps-ink: #E9EBE9;
  --ps-muted: #A0A8AE; --ps-line: #2C3338; --ps-accent: #7CBDF2; --ps-ok: #7FD0A4; --ps-warn: #F0B45A;
  --ps-bad: #FF8A80; --ps-warnbox: #2A2316; --ps-note: #16222E; --q-primary: #3A8AD0 !important; }
.ps-header { background: var(--ps-surface) !important; color: var(--ps-ink) !important;
  border-bottom: 1px solid var(--ps-line); }
.ps-rail { background: var(--ps-sunk) !important; }
.ps-card, .ps-plate, .q-table__card { background: var(--ps-surface) !important; border: 1px solid var(--ps-line);
  box-shadow: none !important; }
.ps-card { border-radius: 8px; }
.ps-plate { border-radius: 10px; padding: 6px 12px 10px; }
.ps-sunk { background: var(--ps-sunk); }
.ps-muted { color: var(--ps-muted); } .ps-accent { color: var(--ps-accent); } .ps-ok { color: var(--ps-ok); }
.ps-warn { color: var(--ps-warn); } .ps-bad { color: var(--ps-bad); }
.ps-warnbox { background: var(--ps-warnbox) !important; border: 1px solid var(--ps-line);
  box-shadow: none !important; }
.ps-note { background: var(--ps-note); }
.ps-section { font-size: 12px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase;
  color: var(--ps-muted); }
.q-table th { color: var(--ps-muted); font-weight: 500; }
.q-table td, .ps-num { font-family: %(mono)s; font-variant-numeric: tabular-nums; }
.ps-readouts { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1px; width: 100%%;
  background: var(--ps-line); border: 1px solid var(--ps-line); border-radius: 10px; overflow: hidden; }
.ps-readout { background: var(--ps-surface); padding: 12px 16px; }
.ps-readout-value { font-family: %(mono)s; font-size: 22px; font-weight: 500; line-height: 1.3; }
.ps-paper { background: #FFFFFF; padding: 12px; max-width: 100%%;
  box-shadow: 0 1px 0 rgba(0, 0, 0, 0.25), 0 8px 28px rgba(0, 0, 0, 0.18); }
.ps-paper img { display: block; max-width: 100%%; }
""" % {"font": FONT, "mono": MONO}


def _go():
    import plotly.graph_objects as go

    return go


def _shown(field: str, value) -> str:
    """A setup value as the text of its input."""
    if value is None:
        return ""
    if field == "absorbers" and isinstance(value, list):
        return ", ".join(f"{m} {t}" for m, t in value)
    return str(value)


def _value(field: str, text):
    """A typed field value from the text the user entered ("" removes an optional field)."""
    if text is None:
        return None
    text = str(text).strip()
    if text == "":
        return None
    if field == "absorbers":  # "Pb 1 mm, Cu 0.5 mm"
        return [part.strip().split(None, 1) if " " in part.strip() else [part.strip(), ""]
                for part in text.split(",") if part.strip()]
    if field in INT_FIELDS:
        try:
            return int(text)
        except ValueError:
            return text  # the setup checks report it
    return text


# ---------------------------------------------------------------------------------------------------------------
# Figures (Plotly), one per tab


def themed(fig, theme: str = "light"):
    """Give a figure the page's look in the "light" or "dark" theme: the page's background and ink, boxed axes
    with inward ticks, a quiet grid. Returns the figure."""
    t = THEMES[theme]
    fig.update_layout(template="plotly_dark" if theme == "dark" else "plotly_white", paper_bgcolor=t["paper"],
                      plot_bgcolor=t["paper"], font=dict(family=FONT, color=t["ink"], size=13), colorway=COLORS,
                      legend=dict(bgcolor="rgba(0,0,0,0)"))
    axis = dict(showline=True, linewidth=1, linecolor=t["ink"], mirror="ticks", ticks="inside", tickcolor=t["ink"],
                gridcolor=t["grid"], zeroline=False)
    fig.update_xaxes(**axis)
    fig.update_yaxes(**axis)
    wall = dict(backgroundcolor=t["paper"], gridcolor=t["grid"], color=t["ink"], showbackground=True)
    fig.update_scenes(xaxis=wall, yaxis=wall, zaxis=wall)
    return fig


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
    for d in g.get("blocking", []):
        o = d["outline"]
        fig.add_trace(go.Scatter3d(x=o[:, 0], y=o[:, 1], z=o[:, 2], mode="lines", name=d["name"],
                                   line=dict(color="#888", width=2), showlegend=False))
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


def figure_efficiency(planner: Planner, key: str, points: Optional[list] = None):
    """A γ-ray detector's efficiency against energy (log–log): the full-energy-peak and total efficiency, the
    energy of the excited state's γ ray, and the ``points`` a simulated source run gave
    (:meth:`physim.nuclear.response.SourceRun.efficiency_points`)."""
    go = _go()
    c = planner.efficiency(key)
    e = 1e3 * c["energy_mev"]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=e, y=100 * c["peak"], mode="lines", name="full-energy peak",
                             line=dict(color=COLORS[0], width=2)))
    fig.add_trace(go.Scatter(x=e, y=100 * c["total"], mode="lines", name="any energy",
                             line=dict(color=COLORS[1], width=1.2, dash="dash")))
    if points:
        fig.add_trace(go.Scatter(
            x=[1e3 * p["energy_mev"] for p in points], y=[100 * p["efficiency"] for p in points],
            error_y=dict(type="data", array=[100 * p["uncertainty"] for p in points], width=0, thickness=1),
            mode="markers", name="source run", marker=dict(color=COLORS[3], size=6)))
    if planner.experiment.excitation is not None:
        fig.add_vline(x=1e3 * planner.gamma_energy_mev(), line=dict(color="#888", width=1, dash="dot"))
    fig.update_layout(xaxis=dict(title="γ-ray energy (keV)", type="log"),
                      yaxis=dict(title="efficiency (%)", type="log"), height=330,
                      margin=dict(l=55, r=10, t=10, b=40), legend=dict(orientation="h", y=-0.32, x=0))
    return fig


def figure_source_spectrum(run, name: str):
    """A γ-ray detector's spectrum from a simulated calibration-source run (counts per bin, log scale)."""
    go = _go()
    h = run.detector(name)
    x = 1e3 * np.repeat(run.edges, 2)[1:-1]
    y = np.repeat(h, 2)
    fig = go.Figure(go.Scatter(x=x, y=np.where(y > 0, y, np.nan), mode="lines",
                               line=dict(color=COLORS[0], width=1)))
    width = 1e3 * (run.edges[1] - run.edges[0])
    fig.update_layout(xaxis=dict(title="energy (keV)"), yaxis=dict(title=f"counts / {width:.2g} keV", type="log"),
                      height=260, margin=dict(l=55, r=10, t=10, b=40), showlegend=False)
    return fig


def figure_detector_spectrum(planner: Planner, detector: str, events: int = 100_000, seed: int = 1):
    """One detector's simulated measured-energy spectrum, for the scene's side panel."""
    go = _go()
    h = planner.spectra(events=events, seed=seed)["spectra"].get(detector)
    fig = go.Figure()
    if h is None:
        fig.add_annotation(text="no events in this detector", showarrow=False, x=0.5, y=0.5, xref="paper",
                           yref="paper")
    else:
        y = np.repeat(h["counts"], 2)
        fig.add_trace(go.Scatter(x=np.repeat(h["edges"], 2)[1:-1], y=np.where(y > 0, y, np.nan), mode="lines",
                                 name=detector, line=dict(color=COLORS[0], width=1.2)))
    fig.update_layout(xaxis=dict(title="measured energy (MeV)"), yaxis=dict(title="counts / bin", type="log"),
                      margin=dict(l=50, r=10, t=10, b=40), height=260, showlegend=False)
    return fig


def figure_gamma_spectra(planner: Planner, detector: str, events: int = 400_000, seed: int = 1):
    """A γ-ray detector's coincidence spectrum: as measured, and Doppler-corrected for the emitting nucleus."""
    go = _go()
    g = planner.gamma_spectra(events, seed)
    fig = go.Figure()
    if not g["available"]:
        fig.add_annotation(text=g["reason"], showarrow=False, x=0.5, y=0.5, xref="paper", yref="paper")
        return fig
    s = g["spectra"][detector]
    x = 1e3 * np.repeat(s["edges"], 2)[1:-1]
    for i, (key, label) in enumerate((("measured", "measured"), (g["emitter"], f"corrected for the {g['emitter']}"),
                                      ("projectile" if g["emitter"] == "recoil" else "recoil",
                                       "corrected for the wrong nucleus"))):
        y = np.repeat(s[key], 2)
        fig.add_trace(go.Scatter(x=x, y=np.where(y > 0, y, np.nan), mode="lines", name=label,
                                 line=dict(color=COLORS[i % len(COLORS)], width=1.2 if i < 2 else 0.8),
                                 visible=True if i < 2 else "legendonly"))
    fig.add_vline(x=g["energy_kev"], line=dict(color="#888", width=1, dash="dot"))
    fig.update_layout(xaxis=dict(title="γ-ray energy (keV)"), yaxis=dict(title="counts / bin in the run"),
                      margin=dict(l=50, r=10, t=10, b=40), height=320, legend=dict(orientation="h", y=-0.3, x=0))
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


def figure_levels(table: dict):
    """A level scheme: a line per level with its spin, parity and energy, and an arrow per γ-ray transition."""
    go = _go()
    fig = go.Figure()
    levels, transitions = table["levels"], table["transitions"]
    n = max(len(transitions), 1)
    for lv in levels:
        e = lv["energy_kev"]
        fig.add_trace(go.Scatter(x=[0, n + 1], y=[e, e], mode="lines", line=dict(color=COLORS[0], width=2),
                                 hovertemplate=f"{lv['jpi']}  {e:g} keV<extra></extra>", showlegend=False))
        fig.add_annotation(x=0, y=e, text=lv["jpi"], xanchor="right", showarrow=False, xshift=-4)
        fig.add_annotation(x=n + 1, y=e, text=f"{e:g}", xanchor="left", showarrow=False, xshift=4)
    for k, t in enumerate(transitions, start=1):
        top, bottom = levels[t["from"]]["energy_kev"], levels[t["to"]]["energy_kev"]
        share = "" if t["branching"] is None else f", {100 * t['branching']:.3g}% of the γ rays"
        fig.add_annotation(x=k, y=bottom, ax=k, ay=top, xref="x", yref="y", axref="x", ayref="y", showarrow=True,
                           arrowhead=2, arrowwidth=1.5, arrowcolor=COLORS[1],
                           hovertext=f"{t['energy_kev']:g} keV {t['multipolarity']}{share}")
    fig.update_layout(xaxis=dict(visible=False, range=[-1, n + 2]), yaxis=dict(title="level energy (keV)"),
                      margin=dict(l=60, r=10, t=10, b=10), height=360)
    return fig


FIGURES = {"geometry": figure_geometry, "kinematics": figure_kinematics, "energy_loss": figure_energy_loss,
           "spectra": figure_spectra, "trajectories": figure_trajectories, "gamma": figure_excitation}


def report_zip(planner: Planner, seed: int = 1, events: int = 200_000, journal: str = "physical_review",
               width: str = "single") -> bytes:
    """The report folder (HTML, CSV, figures in a journal's style, setup, ROOT file if uproot is installed) as a
    zip archive."""
    import tempfile
    from pathlib import Path

    from .report import build

    with tempfile.TemporaryDirectory() as tmp:
        build(planner.experiment, seed=seed, events=events, journal=journal, width=width).write(tmp)
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


def build_page(example: str = "alpha_on_gold", events: int = 100_000, mode: Optional[str] = None,
               theme: Optional[str] = None) -> None:
    """Build the planner page for the current client (call inside a NiceGUI page function).

    ``mode`` is "guided" (a step-by-step workflow) or "expert" (every input and result at once); ``None`` takes the
    one this browser used last, guided the first time. ``theme`` is "light" or "dark"; ``None`` takes the one this
    browser used last, the system's setting the first time."""
    from nicegui import run, ui

    from . import paper
    from .experiment import Experiment, SetupError

    state = {"planner": Planner.example(example), "events": events, "mode": mode or "guided",
             "step": guide.STEPS[0].key, "example": example, "theme": theme or "light",
             "export": None, "journal": "physical_review", "width": "single", "format": "pdf"}

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
        with ui.grid(columns=1 if state["mode"] == "guided" else 2).classes("w-full gap-1"):
            for field, label, placeholder in fields:
                v = values.get(field)
                inp = ui.input(label, value=_shown(field, v), placeholder=placeholder).props(
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
            ui.label("Particle detectors").classes("ps-section")
            ui.space()
            with ui.button("Add", icon="add").props("dense flat"):
                with ui.menu():
                    for pl in guide.PLACEMENTS:
                        with ui.menu_item(on_click=add_detector(pl.key)).classes("max-w-sm"):
                            with ui.column().classes("gap-0"):
                                ui.label(pl.label).classes("text-sm font-medium")
                                ui.label(pl.why).classes("text-xs ps-muted")
        ui.label("Silicon detectors: they measure the energy of the scattered beam particles and recoils. Set θ "
                 "below 90° for forward angles, above 90° for backward ones, or use Add for a ready-made detector "
                 "in each region.").classes("text-xs ps-muted")
        size_fields = set().union(*SHAPE_FIELDS.values())
        for i, det in enumerate(d["detectors"]):
            shape = det.get("shape", "circle")
            with ui.expansion(f"{det.get('name') or f'D{i + 1}'} · {shape} silicon").classes("w-full ps-card"):
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
            ui.label("γ-ray detectors").classes("ps-section")
            ui.space()
            with ui.button("Add", icon="add").props("dense flat"):
                with ui.menu():
                    for key, label, why, _, _ in guide.GAMMA_PLACEMENTS:
                        with ui.menu_item(on_click=add_gamma_detector(key)).classes("max-w-sm"):
                            with ui.column().classes("gap-0"):
                                ui.label(label).classes("text-sm font-medium")
                                ui.label(why).classes("text-xs ps-muted")
        ui.label("They see the γ ray emitted when the excited state decays. Add offers a single crystal, a "
                 "germanium clover and a LaBr₃ scintillator with their real dimensions.").classes(
            "text-xs ps-muted")
        for i, gd in enumerate(d.get("gamma_detectors", [])):
            real = "crystal_diameter" in gd or "model" in gd
            kind = gd.get("model") or ("crystals" if real else "disc")
            with ui.expansion(f"{gd.get('name') or f'γ{i + 1}'} · {kind}").classes("w-full ps-card"):
                fields = [f for f in GAMMA_FIELDS if not (real and f[0] == "radius")]
                section("", fields + (CRYSTAL_FIELDS if real else []), f"gamma detector {i + 1}", gd)
                with ui.row():
                    ui.button("Duplicate", icon="content_copy",
                              on_click=duplicate_gamma_detector(i)).props("dense flat")
                    ui.button("Remove", icon="delete", on_click=remove_gamma_detector(i)).props(
                        "dense flat color=negative")

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
                with ui.card().classes("w-full p-3 gap-1 ps-card"):
                    ui.label(name).classes("font-semibold")
                    ui.label(text).classes("text-sm ps-muted")
                    ui.button("Start here", icon="arrow_forward",
                              on_click=lambda ex=ex: load_example(ex, then="beam")).props(
                        "dense flat no-caps").classes("self-start")
                    ui.label(f"Starts from the example {ex}.").classes("text-xs ps-muted")
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
                    ui.label(step.intro).classes("text-sm ps-muted")
                    step_inputs(step.key)
                    step_status(step.key, keys[i - 1] if i else None, keys[i + 1] if i + 1 < len(keys) else None)

    @ui.refreshable
    def step_status(key: str = "", back: Optional[str] = None, nxt: Optional[str] = None) -> None:
        if not key:
            return
        problems = next((s["problems"] for s in guide.steps_for(P()) if s["step"].key == key), [])
        for p in problems:
            with ui.row().classes("items-start no-wrap gap-1 mt-1"):
                ui.icon("error").classes("ps-bad")
                ui.label(p).classes("text-sm ps-bad")
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
        if not ws:  # the key results above say so
            return
        style = {"error": "ps-bad font-semibold", "warning": "ps-warn", "note": "ps-muted"}
        icon = {"error": "error", "warning": "warning", "note": "info"}
        with ui.card().classes("w-full ps-warnbox p-2 gap-1"):
            for w in ws:
                with ui.row().classes("items-start no-wrap gap-2"):
                    ui.icon(icon[w.level]).classes(style[w.level])
                    ui.label(w.text).classes(style[w.level] + " text-sm")

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
    scene_state = {"view": None, "key": None, "element": None, "live": {},
                   "source": {"nuclide": "152Eu", "activity": "37 kBq", "time": "1 h", "run": False}}

    def scene_selected(key, element) -> None:
        scene_state["key"], scene_state["element"] = key, element
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

    @ui.refreshable
    def selection_panel() -> None:
        key, element = scene_state["key"], scene_state["element"]
        try:
            s = P().selection(key, element)
        except (IndexError, KeyError, ValueError):
            key, s = None, P().selection()
        live = scene_state["live"] = {}

        def row(label: str, value: str, name: Optional[str] = None, classes: str = "") -> None:
            with ui.row().classes("w-full justify-between no-wrap gap-3"):
                ui.label(label).classes("text-sm ps-muted")
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
                row("Solid angle", f"{_fmt(s['solid_angle_msr'])} msr", "omega")
                row("Hidden by others", f"{100 * s['hidden']:.0f} %", "hidden")
                row("Particles reaching it", f"{_fmt(s['rate_per_s'])} /s", "rate")
                if s["measured"] != "all":
                    row(f"Of these, {s['measured']}", f"{_fmt(s['measured_rate_per_s'])} /s")
                row("Counts in the run", _fmt(s["counts_in_run"]))
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
                    f"{100 * s['peak_efficiency']:.3f} %")
                if s["coincidence_rate_per_s"] is not None:
                    row("Particle–γ coincidences", f"{_fmt(s['coincidence_rate_per_s'])} /s")
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
        from . import response as _response

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

    @ui.refreshable
    def geometry_panel():
        from .scene_view import SceneView

        with ui.row().classes("w-full items-start gap-3").style("flex-wrap: wrap"):
            with ui.column().classes("ps-plate gap-1").style("flex: 1 1 520px; min-width: 0"):
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
                view = SceneView(P, state["theme"], scene_selected, scene_live, scene_moved)
                scene_state["view"] = view
                if scene_state["key"] is not None:
                    view.select(scene_state["key"], scene_state["element"], notify=False)
                ui.label("To scale; the numbers along the beam are mm from the target. Click a detector, a ring "
                         "or a crystal to select it; drag a selected detector to move it. Drag the background to "
                         "turn the view, scroll to zoom.").classes(
                    "text-xs ps-muted")
            with ui.column().classes("gap-2").style("flex: 0 1 360px; min-width: 300px"):
                selection_panel()
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
                 "[run] section.").classes("text-xs ps-muted")
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
        analysis_block()

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
            r = await run.io_bound(P().analyse, settings, max(state["events"], 400_000), 1)
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
            with ui.element("div").classes("ps-readouts mt-2"):
                for title, value, sub in (
                        ("B(E2↑)", f"{r['b_e2b2']:.4g} e²b²", f"{b:.4g} e²fm⁴ · {r['b_wu']:.3g} W.u."),
                        ("Uncertainty", f"± {100 * r['total_unc']:.1f} %",
                         f"statistical {100 * r['statistical']:.1f} %, systematic {100 * r['systematic']:.1f} %; "
                         f"the Monte Carlo sample itself adds {100 * r['monte_carlo']:.1f} %"),
                        ("Put in", f"{1e-4 * r['truth_e2fm4']:.4g} e²b²", f"pull {r['pull']:+.2f} σ"),
                        ("Beam time", f"{r['hours_for_precision']:.3g} h",
                         f"for {100 * r['settings']['wanted_precision']:.3g} % statistics; "
                         f"{r['counts_per_shift']:.3g} counts in the peak per {SHIFT_H:g} h shift")):
                    with ui.element("div").classes("ps-readout"):
                        ui.label(title).classes("text-xs ps-muted")
                        ui.label(value).classes("ps-readout-value")
                        ui.label(sub).classes("text-xs ps-muted")
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
    def gamma_panel():
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

    panels = {"geometry": geometry_panel, "kinematics": kinematics_panel, "rates": rates_panel,
              "energy_loss": energy_loss_panel, "spectra": spectra_panel, "trajectories": trajectories_panel,
              "gamma": gamma_panel, "report": report_panel}

    def refresh_results() -> None:
        readouts.refresh()
        warnings_banner.refresh()
        reading_box.refresh()
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
        with ui.tabs().classes("w-full").props("dense no-caps align=left") as tabs:
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

    def set_theme(name: str) -> None:
        if name == state["theme"]:
            return
        state["theme"] = name
        dark.value = name == "dark"
        ui.run_javascript(f"try {{ localStorage.setItem('physim-planner-theme', '{name}') }} catch (e) {{}}")
        refresh_results()  # the figures are drawn in the theme's colours

    async def restore_choices() -> None:
        """The mode and theme this browser used last (the system's light or dark setting the first time), unless
        the address gives them (?mode=...&theme=...)."""
        try:
            saved_mode, saved_theme, system_dark = await ui.run_javascript(
                "(() => { let m = null, t = null; try { m = localStorage.getItem('physim-planner-mode'); "
                "t = localStorage.getItem('physim-planner-theme') } catch (e) {} "
                "return [m, t, window.matchMedia('(prefers-color-scheme: dark)').matches] })()", timeout=3)
        except Exception:  # noqa: BLE001 -- no answer from the browser: keep the defaults
            return
        if theme is None:
            wanted = saved_theme if saved_theme in THEMES else ("dark" if system_dark else "light")
            if wanted != state["theme"]:
                theme_toggle.value = wanted
        if mode is None and saved_mode in ("guided", "expert") and saved_mode != state["mode"]:
            mode_toggle.value = saved_mode

    # -- layout -----------------------------------------------------------------------------------------------
    ui.add_css(STYLE)
    dark = ui.dark_mode(state["theme"] == "dark")
    with ui.header(elevated=False).classes("items-center ps-header py-1"):
        ui.label("physim").classes("text-lg font-semibold")
        ui.label("experiment planner").classes("ps-muted")
        ui.space()
        mode_toggle = ui.toggle({"guided": "Guided", "expert": "Expert view"}, value=state["mode"],
                                on_change=lambda e: set_mode(e.value)).props("dense no-caps unelevated")
        theme_toggle = ui.toggle({"light": "Light", "dark": "Dark"}, value=state["theme"],
                                 on_change=lambda e: set_theme(e.value)).props("dense no-caps unelevated")
        example_select = ui.select(
            Planner.examples(), value=example, label="Start from example",
            on_change=lambda e: None if e.value == state["example"] else load_example(e.value)).props(
            "dense outlined").classes("w-60")
        ui.button("Load setup", icon="upload", on_click=lambda: upload_dialog.open()).props("flat no-caps")
        ui.button("Save setup", icon="download", on_click=save_setup).props("flat no-caps")
    with ui.dialog() as upload_dialog, ui.card():
        ui.label("Load a setup file (.toml)")
        ui.upload(auto_upload=True, on_upload=load_file).props("accept=.toml max-files=1")
    with ui.dialog() as export_dialog, ui.card().classes("ps-card").style("max-width: min(1100px, 95vw)"):
        with ui.row().classes("items-start gap-6"):
            with ui.column().classes("gap-2").style("width: 320px; max-width: 100%"):
                export_controls()
            preview_box = ui.column().classes("items-center gap-2")
    export_dialog.on("show", render_preview)
    with ui.left_drawer(value=True).classes("ps-rail").props("width=420 bordered behavior=desktop"):
        setup_panel()
    with ui.column().classes("w-full gap-2"):
        readouts()
        main_area()
    ui.timer(0.2, restore_choices, once=True)


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
    def index(example: str = args.example, mode: Optional[str] = None, theme: Optional[str] = None):
        build_page(example, mode=mode if mode in ("guided", "expert") else None,
                   theme=theme if theme in THEMES else None)

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
           "report_zip", "themed"]
