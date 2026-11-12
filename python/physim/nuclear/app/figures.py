"""The planner app's figures (Plotly), its page style and colours, and the report download: module-level, so
they work in a notebook as well."""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from .. import workbench as _wb
from ..planner import Planner

#: Fields of each setup section, (field, label, placeholder): a decision card's own fields and those behind its
#: "more" fold (:mod:`physim.nuclear.workbench`).
BEAM_FIELDS = _wb.BEAM_FIELDS + _wb.BEAM_MORE
TARGET_FIELDS = _wb.TARGET_FIELDS + _wb.TARGET_MORE
BACKING_FIELDS = _wb.BACKING_FIELDS
RUN_FIELDS = _wb.RUN_FIELDS + _wb.RUN_MORE
DETECTOR_FIELDS = _wb.DETECTOR_FIELDS + _wb.DETECTOR_MORE
REACTION_FIELDS = _wb.REACTION_FIELDS
REACTION_TYPES = _wb.REACTION_TYPES
GAMMA_FIELDS = _wb.GAMMA_FIELDS + _wb.GAMMA_MORE
CRYSTAL_FIELDS = _wb.CRYSTAL_FIELDS
CLOVER_FIELDS = _wb.CLOVER_FIELDS
SHAPE_FIELDS = _wb.SHAPE_FIELDS
INT_FIELDS = {"charge_state", "counts_wanted", "strips_x", "strips_y", "rings", "sectors", "crystals"}

#: Detector and curve colours (Okabe–Ito first): distinguishable for colour-blind readers, legible on the light and
#: the dark page.
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F", "#882255", "#44AA99",
          "#DDCC77"]

#: Figure colours of each page theme.
THEMES = {"light": {"paper": "#FFFFFF", "ink": "#12171C", "grid": "#E7EBEE"},
          "dark": {"paper": "#12171C", "ink": "#E6EAED", "grid": "#1D242B"}}

#: The curve colours on the dark page, in the order of :data:`COLORS`: the same roles, light enough to read on
#: the dark panels (the first is the data's cyan, the second the accent's amber).
DARK_COLORS = ["#5CC8F0", "#F2A33A", "#5BD08A", "#E28BC0", "#F0D060", "#8FA8FF", "#A7B1BA", "#D77A9B", "#6FD6C4",
               "#E8DC9A"]

FONT = '"IBM Plex Sans", "Segoe UI", system-ui, -apple-system, sans-serif'
MONO = '"IBM Plex Mono", ui-monospace, "Cascadia Mono", Consolas, monospace'

#: The page's style sheet. Colours are variables, set once for the light page and once for the dark one
#: (``body--dark`` is Quasar's dark-mode class). No fonts or scripts are fetched: the app works offline.
STYLE = """
body { --ps-ground: #EEF1F3; --ps-surface: #FFFFFF; --ps-sunk: #F6F8F9; --ps-raise: #E7EBEE; --ps-ink: #12171C;
  --ps-muted: #56616C; --ps-line: #D5DBE0; --ps-hair: #E7EBEE; --ps-accent: #8F5400; --ps-fill: #F2A33A;
  --ps-on-fill: #1A1204; --ps-data: #0072B2; --ps-ok: #17663A; --ps-warn: #B5470B; --ps-bad: #B3261E;
  --ps-warnbox: #FBEFE6; --ps-note: #F6F8F9; --q-primary: #8F5400 !important;
  --ps-mono: %(mono)s;
  background: var(--ps-ground) !important; color: var(--ps-ink); font-family: %(font)s; font-size: 14px; }
body.body--dark { --ps-ground: #0B0E11; --ps-surface: #12171C; --ps-sunk: #0F1317; --ps-raise: #1D242B;
  --ps-ink: #E6EAED; --ps-muted: #9AA6B1; --ps-line: #242C34; --ps-hair: #1A2026; --ps-accent: #F2A33A;
  --ps-fill: #F2A33A; --ps-on-fill: #1A1204; --ps-data: #5CC8F0; --ps-ok: #5BD08A; --ps-warn: #FF8F5A;
  --ps-bad: #FF7B7B; --ps-warnbox: #2A1A12; --ps-note: #171D23; --q-primary: #F2A33A !important; }
.ps-header { background: var(--ps-surface) !important; color: var(--ps-ink) !important;
  border-bottom: 1px solid var(--ps-line); }
.ps-rail, .ps-nav { background: var(--ps-sunk) !important; color: var(--ps-ink); }
.q-drawer { background: var(--ps-sunk) !important; }
.q-drawer--bordered { border-color: var(--ps-line) !important; }
.ps-card, .ps-plate, .q-table__card { background: var(--ps-surface) !important; border: 1px solid var(--ps-line);
  box-shadow: none !important; }
.ps-card, .q-table__card { border-radius: 6px; }
.ps-plate { border-radius: 6px; padding: 6px 12px 10px; }
.ps-sunk { background: var(--ps-sunk); }
.ps-muted { color: var(--ps-muted) !important; } .ps-accent { color: var(--ps-data) !important; }
.ps-ok { color: var(--ps-ok); }
.ps-warn { color: var(--ps-warn); } .ps-bad { color: var(--ps-bad); }
.ps-warnbox { background: var(--ps-warnbox) !important; border: 1px solid var(--ps-line);
  box-shadow: none !important; }
.ps-note { background: var(--ps-note); border: 1px solid var(--ps-line); border-radius: 6px; }
.ps-section { font-size: 13px; font-weight: 600; color: var(--ps-muted); }
/* Controls: square corners of 6 px, filled buttons in the accent with dark text, quiet menus and dialogs. */
.q-btn { border-radius: 6px; font-weight: 500; text-transform: none; }
.q-btn--round, .q-btn--fab { border-radius: 50%%; }
/* (Quasar's own !important rules live in a cascade layer, where they beat any !important from outside it: these two
   are put in that layer, and win there by their selectors.) */
@layer quasar_importants {
  body .bg-primary { background: var(--ps-fill) !important; }
  body .q-btn.bg-primary, body .q-btn.bg-primary .q-icon, body .bg-primary.text-white {
    color: var(--ps-on-fill) !important; }
}
.q-btn--outline { color: var(--ps-ink) !important; }
.q-btn--outline::before { border-color: var(--ps-line) !important; }
.q-btn-group { border-radius: 6px; box-shadow: none; border: 1px solid var(--ps-line); }
.q-btn-group > .q-btn-item { border-radius: 5px !important; }
.q-field--outlined .q-field__control { border-radius: 6px; }
.q-field--outlined .q-field__control:before { border-color: var(--ps-line); }
.q-menu, .q-dialog .q-card { background: var(--ps-surface); color: var(--ps-ink); border: 1px solid var(--ps-line);
  border-radius: 6px; box-shadow: 0 12px 40px rgba(0, 0, 0, 0.35) !important; }
.q-expansion-item { border-radius: 6px; }
.q-separator { background: var(--ps-line); }
.q-table th { color: var(--ps-muted); font-weight: 500; }
.q-table td, .ps-num { font-family: %(mono)s; font-variant-numeric: tabular-nums; }
.q-table tbody td { border-color: var(--ps-hair); }
.q-table thead th { border-color: var(--ps-line); }
/* The window the setup is edited in: the whole height at the right of the page, its cards scrolling between
   its head and its foot. */
.ps-setup-window { width: min(560px, 96vw); max-width: 96vw !important; height: 100vh;
  max-height: 100vh !important; padding: 0; gap: 0; display: flex; flex-direction: column; flex-wrap: nowrap;
  border-radius: 0 !important; background: var(--ps-sunk) !important; border-width: 0 0 0 1px !important; }
.q-dialog__inner:has(> .ps-setup-window) { padding: 0 !important; }
.ps-setup-body { flex: 1 1 auto; min-height: 0; overflow-y: auto; padding: 10px 14px; width: 100%%; }
.ps-setup-foot { border-top: 1px solid var(--ps-line); padding: 10px 14px; }
/* A panel: a bordered surface with a head line (its name and what it shows) over its content. */
.ps-panel { background: var(--ps-surface); border: 1px solid var(--ps-line); border-radius: 6px; width: 100%%;
  min-width: 0; }
.ps-panel-head { display: flex; flex-wrap: wrap; align-items: center; gap: 4px 8px; padding: 9px 14px;
  border-bottom: 1px solid var(--ps-line); font-weight: 600; }
.ps-panel-head .ps-muted { font-weight: 400; }
.ps-panel-body { padding: 10px 14px 12px; }
.ps-readouts { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 1px; width: 100%%;
  background: var(--ps-line); border: 1px solid var(--ps-line); border-radius: 6px; overflow: hidden; }
.ps-readout { background: var(--ps-surface); padding: 12px 16px; }
.ps-readout-value { font-family: %(mono)s; font-size: 24px; font-weight: 500; line-height: 1.3; }
.ps-readout-main .ps-readout-value { font-size: 40px; line-height: 1.15; color: var(--ps-accent); }
.ps-readout-main { grid-column: span 2; }
/* The chain of an analysis: a few large numbers, what is done between them, the last one in the accent. */
.ps-chain { display: flex; flex-direction: column; gap: 8px; padding: 14px; }
.ps-chain-step { border: 1px solid var(--ps-line); border-radius: 6px; padding: 10px 14px; }
.ps-chain-step .ps-chain-value { font-family: %(mono)s; font-size: 20px; font-weight: 500; }
.ps-chain-step.ps-chain-last { background: var(--ps-raise); border-color: var(--ps-accent); }
.ps-chain-step.ps-chain-last .ps-chain-value { font-size: 30px; color: var(--ps-accent); }
.ps-chain-op { font-family: %(mono)s; font-size: 12px; color: var(--ps-muted); padding: 0 14px; }
.ps-chain-op::before { content: "\\2193"; margin-right: 8px; }
.ps-facts { padding: 6px 14px 10px; border-top: 1px solid var(--ps-line); }
.ps-fact { display: flex; justify-content: space-between; gap: 12px; padding: 4px 0; font-size: 13px; }
.ps-fact > :first-child { color: var(--ps-muted); }
.ps-fact > :last-child { font-family: %(mono)s; text-align: right; }
.ps-bar { height: 4px; margin-top: 6px; background: var(--ps-raise); }
.ps-bar > div { height: 4px; background: var(--ps-fill); }
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
    if field in ("absorbers", "contaminants") and isinstance(value, list):
        return ", ".join(f"{m} {t}" for m, t in value)
    return str(value)


def _crystals_of(gd: dict) -> Optional[int]:
    """How many crystals a γ-ray detector of the draft has, from its fields or its model."""
    if gd.get("crystals") is not None:
        return gd["crystals"] if isinstance(gd["crystals"], int) else None
    if gd.get("model"):
        from .. import catalogue

        try:
            return catalogue.model(gd["model"], "gamma").fields.get("crystals")
        except ValueError:
            return None
    return None


def _value(field: str, text):
    """A typed field value from the text the user entered ("" removes an optional field)."""
    if text is None:
        return None
    text = str(text).strip()
    if text == "":
        return None
    if field in ("absorbers", "contaminants"):  # "Pb 1 mm, Cu 0.5 mm"
        return [part.strip().split(None, 1) if " " in part.strip() else [part.strip(), ""]
                for part in text.split(",") if part.strip()]
    if field in INT_FIELDS:
        try:
            return int(text)
        except ValueError:
            return text  # the setup checks report it
    if field == "addback":
        return text.lower() in ("true", "yes", "on", "1")
    if field in ("addback_factor", "suppression_factor"):
        try:
            return float(text)
        except ValueError:
            return text
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
    if theme == "dark":
        _recolour(fig, dict(zip(COLORS, DARK_COLORS)))
        fig.update_layout(colorway=DARK_COLORS)
    return fig


def _recolour(fig, to: dict) -> None:
    """Replace the curve colours of a figure by their counterparts in ``to`` ({colour: colour}): the lines, markers
    and fills of its traces, and its shapes and annotations. A colour that is not a key stays."""
    def swap(obj, *names) -> None:
        for name in names:
            try:
                value = getattr(obj, name, None)
                if isinstance(value, str) and value.upper() in to:
                    setattr(obj, name, to[value.upper()])
            except Exception:  # noqa: BLE001 -- a trace type without this property
                pass

    for trace in fig.data:
        swap(trace, "fillcolor")
        for part in ("line", "marker", "textfont", "error_y", "error_x"):
            child = getattr(trace, part, None)
            if child is not None:
                swap(child, "color")
                inner = getattr(child, "line", None)
                if inner is not None:
                    swap(inner, "color")
    for shape in fig.layout.shapes or ():
        swap(shape, "fillcolor")
        swap(shape.line, "color")
    for note in fig.layout.annotations or ():
        swap(note, "arrowcolor", "bordercolor")
        swap(note.font, "color")


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
    if "plain" in s:
        # The same γ rays without add-back or suppression, for the comparison.
        m = g["modes"][detector]
        what = " and ".join(w for w, on in (("add-back", m["addback"]), ("suppression", m["shield"])) if on)
        for i, (key, label) in enumerate(((g["emitter"], f"corrected, without {what}"),
                                          ("measured", f"measured, without {what}"))):
            y = np.repeat(s["plain"][key], 2)
            fig.add_trace(go.Scatter(x=x, y=np.where(y > 0, y, np.nan), mode="lines", name=label,
                                     line=dict(color=COLORS[(3 + i) % len(COLORS)], width=1, dash="dash"),
                                     visible=True if i == 0 else "legendonly"))
    fig.add_vline(x=g["energy_kev"], line=dict(color="#888", width=1, dash="dot"))
    fig.update_layout(xaxis=dict(title="γ-ray energy (keV)"), yaxis=dict(title="counts / bin in the run"),
                      margin=dict(l=50, r=10, t=10, b=40), height=320, legend=dict(orientation="h", y=-0.3, x=0))
    return fig


def figure_alignment(planner: Planner, result: dict):
    """The diagnostic plot of a misplaced target: the corrected peak's centroid against ring, one line per
    crystal, with the transition energy as a dashed line. Flat when the geometry assumed is right."""
    go = _go()
    fig = go.Figure()
    rows = result["diagnostic"]
    crystals = list(dict.fromkeys(r["crystal"] for r in rows))
    dets = list(dict.fromkeys(r["detector"] for r in rows))
    for i, cr in enumerate(crystals):
        for j, d in enumerate(dets):
            mine = [r for r in rows if r["crystal"] == cr and r["detector"] == d]
            if not mine:
                continue
            fig.add_trace(go.Scatter(x=[r["ring"] + 1 for r in mine], y=[r["centroid_kev"] for r in mine],
                                     error_y=dict(type="data", array=[r["error_kev"] for r in mine], width=0,
                                                  thickness=1),
                                     mode="lines+markers", name=f"{cr} with {d}",
                                     line=dict(color=COLORS[i % len(COLORS)], width=1, dash=("solid", "dash", "dot")[j % 3]),
                                     marker=dict(size=5)))
    fig.add_hline(y=result["energy_kev"], line=dict(color="#888", width=1, dash="dot"))
    fig.update_layout(xaxis=dict(title="ring or strip"), yaxis=dict(title="corrected centroid (keV)"), height=360,
                      margin=dict(l=55, r=10, t=10, b=40), legend=dict(font=dict(size=10)))
    return fig


def figure_overlay(result: dict):
    """The corrected peak with the true geometry and with the assumed one."""
    go = _go()
    o = result["overlay"]
    x = 1e3 * np.repeat(o["edges"], 2)[1:-1]
    fig = go.Figure()
    for i, (key, label) in enumerate((("true", "true geometry"), ("assumed", f"target assumed {result['offset_mm']:+g} mm off"))):
        y = np.repeat(o[key], 2)
        fig.add_trace(go.Scatter(x=x, y=y, mode="lines", name=label, line=dict(color=COLORS[i], width=1.2)))
    fig.update_layout(xaxis=dict(title="corrected γ-ray energy (keV)"), yaxis=dict(title="counts / bin in the run"),
                      height=300, margin=dict(l=55, r=10, t=10, b=40), legend=dict(orientation="h", y=-0.3, x=0))
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
    """The report as a zip archive (:func:`physim.nuclear.logbook.experiment_zip`): ``report.html`` and, where a
    browser can print it, ``report.pdf``; the setup file; the plan's tables (``tables/``, CSV) and figures in a
    journal's style (``figures/``, PNG and PDF); and the experiment folder with every run, if there is one."""
    from .. import logbook

    return logbook.experiment_zip(planner, None, figures=True, plan_events=events, seed=seed, journal=journal,
                                  width=width, validation=True)


# -- the Data tab's figures (backlog item 67) -------------------------------------------------------------------

def figure_histogram(spectra: list, height: int = 200, log: bool = False, title: Optional[str] = None):
    """One or more spectra (dicts with "counts", "edges", "label", and optionally "scale" for the edges,
    "xlabel" and "view", a (low, high) window of the axis to show) as step lines; small for the grid."""
    go = _go()
    fig = go.Figure()
    for i, s in enumerate(spectra):
        e = np.asarray(s["edges"], dtype=float) * s.get("scale", 1.0)
        fig.add_trace(go.Scatter(x=np.repeat(e, 2)[1:-1], y=np.repeat(np.asarray(s["counts"]), 2),
                                 mode="lines", name=s.get("label", ""), line=dict(color=COLORS[i % len(COLORS)],
                                                                                   width=1.2)))
    fig.update_layout(height=height, margin=dict(l=45, r=8, t=24 if title else 6, b=30),
                      showlegend=len(spectra) > 1, legend=dict(orientation="h", y=-0.3, x=0),
                      title=dict(text=title, font=dict(size=12)) if title else None,
                      xaxis=dict(title=spectra[0].get("xlabel", "energy (MeV)") if spectra else ""),
                      yaxis=dict(title="counts", type="log" if log else "linear"))
    view = next((s["view"] for s in spectra if s.get("view")), None)
    if view is not None:
        fig.update_xaxes(range=[float(view[0]), float(view[1])])
    return fig


def figure_search(spectrum: dict, correction: str, height: int = 300):
    """A whole Doppler-corrected spectrum with the peaks the search found (``spectrum`` is one entry of an
    analysis result's ``search``): a named peak carries its transition, one that is not named a grey mark."""
    go = _go()
    fig = go.Figure()
    edges = 1e3 * np.asarray(spectrum["edges"], dtype=float)
    counts = np.asarray(spectrum["counts"], dtype=float)
    # Two bins to one: the search works on bins a quarter of the resolution wide, finer than the eye needs.
    m = len(counts) // 2 * 2
    counts, edges = counts[:m].reshape(-1, 2).sum(axis=1), edges[:m + 1:2]
    y = np.repeat(counts, 2)
    fig.add_trace(go.Scatter(x=np.repeat(edges, 2)[1:-1], y=np.where(y > 0, y, np.nan), mode="lines",
                             name=f"corrected for the {correction}", line=dict(color=COLORS[0], width=1)))
    centres = (edges[:-1] + edges[1:]) / 2
    for named, color, label in ((True, COLORS[1], "peak found, on a transition"),
                                (False, "#888888", "peak found, not on a transition")):
        peaks = [pk for pk in spectrum["peaks"] if bool(pk["label"]) == named]
        if not peaks:
            continue
        tops = []
        for pk in peaks:
            near = np.abs(centres - pk["energy_kev"]) <= max(pk["fwhm_kev"], edges[1] - edges[0])
            tops.append(float(np.nanmax(counts[near])) if near.any() else 1.0)
        fig.add_trace(go.Scatter(
            x=[pk["energy_kev"] for pk in peaks], y=[max(t, 1.0) * 1.6 for t in tops], mode="markers",
            name=label, marker=dict(symbol="triangle-down", size=9, color=color),
            hovertext=[(pk["label"] + ": " if pk["label"] else "") + f"{pk['energy_kev']:.1f} keV, FWHM "
                       f"{pk['fwhm_kev']:.1f} keV, area {pk['area']:.0f} ± {pk['area_unc']:.0f}"
                       + (f"; near {pk['near']}" if pk.get("near") else "") for pk in peaks], hoverinfo="text"))
        for pk, t in zip(peaks, tops):
            if pk["label"]:  # slanted, so the names of neighbouring lines do not run into each other
                name = pk["label"].split(" ", 1)[1] if pk.get("nuclide") and " " in pk["label"] else pk["label"]
                fig.add_annotation(x=pk["energy_kev"], y=math.log10(max(t, 1.0) * 2.2), text=name, showarrow=False,
                                   textangle=-50, xanchor="left", yanchor="bottom", font=dict(size=10))
    top = float(np.nanmax(counts)) if len(counts) and np.nanmax(counts) > 0 else 1.0
    fig.update_layout(height=height, margin=dict(l=55, r=10, t=10, b=40),
                      xaxis=dict(title=f"γ-ray energy corrected for the {correction} (keV)"),
                      yaxis=dict(title="counts / bin", type="log", range=[-0.3, math.log10(top) + 1.7]),
                      legend=dict(orientation="h", y=-0.28, x=0))
    return fig


def figure_peak(line: dict, fit_half_width: float = 4.0, height: int = 300):
    """One integrated peak (an entry of an analysis result's ``lines``): the gated, corrected spectrum round the
    transition, the random coincidences expected, the fitted Gaussian on its line, and the ±3σ window whose
    counts over the line are the area."""
    from ..response import FWHM_PER_SIGMA

    go = _go()
    fig = go.Figure()
    sp, fit = line["spectrum"], line["fit"]
    edges = 1e3 * np.asarray(sp["edges"], dtype=float)
    counts = np.asarray(sp["counts"], dtype=float)
    rnd = np.asarray(sp.get("randoms") or [], dtype=float)
    k = 4 if len(counts) >= 40 else 1   # the fit's bins are a quarter of the resolution: four to one for the eye
    m = len(counts) // k * k
    counts, edges = counts[:m].reshape(-1, k).sum(axis=1), edges[:m + 1:k]
    x = np.repeat(edges, 2)[1:-1]
    fig.add_trace(go.Scatter(x=x, y=np.repeat(counts, 2), mode="lines", name=f"corrected for the {line['correction']}",
                             line=dict(color=COLORS[0], width=1.2)))
    if rnd.size >= m and rnd.any():
        fig.add_trace(go.Scatter(x=x, y=np.repeat(rnd[:m].reshape(-1, k).sum(axis=1), 2), mode="lines",
                                 name="random coincidences expected", line=dict(color="#888888", width=1, dash="dash")))
    sigma, centre = float(fit.get("sigma") or 0.0), float(fit.get("centroid") or 0.0)
    if sigma > 0 and (fit.get("height", 0) > 0 or fit.get("level", 0) > 0):
        half = fit_half_width * FWHM_PER_SIGMA * sigma
        e = np.linspace(centre - half, centre + half, 300)
        curve = (fit["height"] * np.exp(-0.5 * ((e - centre) / sigma) ** 2) + fit["level"]
                 + fit["slope"] * (e - fit["origin"]))
        base = fit["level"] + fit["slope"] * (e - fit["origin"])
        fig.add_trace(go.Scatter(x=1e3 * e, y=k * base, mode="lines", name="the line under the peak",
                                 line=dict(color=COLORS[1], width=1, dash="dot")))
        fig.add_trace(go.Scatter(x=1e3 * e, y=k * curve, mode="lines", name="fit: Gaussian on a line",
                                 line=dict(color=COLORS[1], width=1.6)))
        fig.add_vrect(x0=1e3 * (centre - 3 * sigma), x1=1e3 * (centre + 3 * sigma), fillcolor=COLORS[1], opacity=0.10,
                      line_width=0, annotation_text=f"area {line['area']:.0f} ± {line['area_unc']:.0f}",
                      annotation_position="top left", annotation_font_size=11)
    fig.add_vline(x=line["energy_kev"], line=dict(color="#888", width=1, dash="dot"))
    # The legend stands above the plot: below it, it ran into the axis title when it wrapped to several lines.
    fig.update_layout(height=height, margin=dict(l=55, r=10, t=10, b=40),
                      xaxis=dict(title=f"γ-ray energy corrected for the {line['correction']} (keV)"),
                      yaxis=dict(title=f"counts / {edges[1] - edges[0]:.2g} keV"),
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, font=dict(size=11)))
    return fig


def figure_energy_ring(view: dict, height: int = 380):
    """Measured energy against ring, with the kinematic line of each group."""
    go = _go()
    e = view["edges"]
    centres = (e[:-1] + e[1:]) / 2
    fig = go.Figure(go.Heatmap(x=view["rings"] + 1, y=centres, z=np.log10(np.maximum(view["counts"].T, 0.5)),
                               colorscale="Viridis", colorbar=dict(title="log₁₀ counts"), hoverinfo="skip"))
    for i, line in enumerate(view["lines"]):
        fig.add_trace(go.Scatter(x=np.asarray(line["rings"]) + 1, y=line["energy"], mode="lines+markers",
                                 name=line["label"], line=dict(color=COLORS[(i + 1) % len(COLORS)], width=1.5,
                                                               dash="dot" if line["group"] == "elastic" else "solid"),
                                 marker=dict(size=4)))
    fig.update_layout(height=height, margin=dict(l=55, r=10, t=10, b=40), dragmode="select",
                      xaxis=dict(title=f"{view['detector']}: ring or strip"), yaxis=dict(title="energy (MeV)"),
                      legend=dict(orientation="h", y=-0.22, x=0))
    return fig


def figure_gamma_crystal(view: dict, height: int = 360):
    """γ-ray energy against crystal, raw or corrected."""
    go = _go()
    e = view["edges"]
    centres = 1e3 * (e[:-1] + e[1:]) / 2
    fig = go.Figure(go.Heatmap(x=view["crystals"], y=centres, z=view["counts"].T, colorscale="Viridis",
                               colorbar=dict(title="counts")))
    fig.update_layout(height=height, margin=dict(l=60, r=10, t=10, b=60),
                      xaxis=dict(title="crystal"), yaxis=dict(title=f"γ-ray energy, {view['correction']} (keV)"))
    return fig
