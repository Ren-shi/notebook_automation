"""Paper-ready figures in a journal's style (matplotlib, imported on first use).

Each of the planner's figures can be drawn at a journal's column width, with its typeface and label size, and saved
as a vector PDF or SVG or a 600 dpi PNG::

    from physim.nuclear import paper
    from physim.nuclear.planner import Planner

    p = Planner.example("alpha_on_gold")
    fig = paper.figure(p, "kinematics", journal="physical_review", width="single")
    fig.savefig("kinematics.pdf")

    pdf = paper.export(p, "spectra", fmt="pdf", journal="nature", width="double", events=200_000)
    csv = paper.data_csv(p, "spectra", events=200_000)       # the plotted numbers

The figures are white with black, boxed axes and inward ticks; every curve has its own line style as well as its
own colour (the Okabe–Ito palette), so they survive greyscale printing and colour-blind readers. Text stays text in
the PDF and SVG files, so it can be edited afterwards.

The widths and type sizes in :data:`JOURNALS` follow each publisher's figure guidelines at the time of writing;
check the journal's current instructions before submitting.
"""

from __future__ import annotations

import io
import math
from dataclasses import dataclass

import numpy as np

#: Okabe–Ito colours, darkest first: distinguishable with every common colour-vision deficiency.
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
#: Line styles cycled together with the colours, so curves differ without colour too.
DASHES = ["-", "--", "-.", ":", (0, (5, 1, 1, 1, 1, 1)), (0, (8, 2))]

MM = 1 / 25.4  # inches per millimetre


@dataclass(frozen=True)
class Journal:
    """A journal family's figure style: column widths in mm, typeface and label size in points."""

    label: str
    widths: dict
    family: str  # "serif" or "sans-serif"
    size: float
    note: str


#: Figure styles by key. Widths are the publishers' column widths.
JOURNALS = {
    "physical_review": Journal("Physical Review (APS)", {"single": 86.0, "double": 178.0}, "serif", 8.0,
                               "Phys. Rev. C, Phys. Rev. Lett.: 8.6 cm and 17.8 cm columns."),
    "nature": Journal("Nature", {"single": 89.0, "double": 183.0}, "sans-serif", 7.0,
                      "Nature and Nature Physics: 89 mm and 183 mm, sans-serif lettering of 5–7 pt."),
    "science": Journal("Science", {"single": 57.0, "middle": 121.0, "double": 184.0}, "sans-serif", 7.0,
                       "Science: 5.7 cm, 12.1 cm and 18.4 cm, Helvetica lettering."),
    "elsevier": Journal("Elsevier (Nucl. Phys. A, Phys. Lett. B, NIM)", {"single": 90.0, "middle": 140.0,
                                                                         "double": 190.0}, "sans-serif", 8.0,
                        "Elsevier journals: 90 mm, 140 mm and 190 mm."),
    "springer": Journal("Springer (Eur. Phys. J. A)", {"single": 84.0, "double": 174.0}, "sans-serif", 8.0,
                        "Springer journals: 84 mm and 174 mm, sans-serif lettering of 8–12 pt."),
}

#: Width names in the order they are offered, with their labels.
WIDTHS = {"single": "1 column", "middle": "1.5 columns", "double": "2 columns"}

#: The figures :func:`figure` can draw, with their titles.
FIGURES = {"geometry": "Geometry", "kinematics": "Kinematics", "strips": "Rates per strip",
           "energy_loss": "Energy loss", "spectra": "Spectra", "trajectories": "Trajectories",
           "gamma": "Excitation probability", "sweep": "Sweep"}

#: File formats of :func:`export`: label and savefig options.
FORMATS = {"pdf": ("PDF (vector)", {}), "svg": ("SVG (vector)", {}), "png": ("PNG, 600 dpi", {"dpi": 600})}

#: Height over width of each figure.
_ASPECT = {"geometry": 0.9, "trajectories": 0.72, "strips": 0.8}


def journal(key) -> Journal:
    """The :class:`Journal` for a key of :data:`JOURNALS` (or the journal itself)."""
    if isinstance(key, Journal):
        return key
    try:
        return JOURNALS[key]
    except KeyError:
        raise ValueError(f"unknown journal {key!r}; choose one of {', '.join(JOURNALS)}") from None


def width_mm(journal_key="physical_review", width="single") -> float:
    """A journal's column width in mm; ``width`` is "single", "middle", "double" or a number of mm."""
    if isinstance(width, (int, float)):
        return float(width)
    widths = journal(journal_key).widths
    if width not in widths:
        raise ValueError(f"{journal(journal_key).label} has no {width!r} width; choose one of {', '.join(widths)}")
    return widths[width]


def rc(journal_key="physical_review") -> dict:
    """The matplotlib settings of a journal's style, for ``matplotlib.rc_context``."""
    j = journal(journal_key)
    return {
        "font.family": j.family,
        "font.serif": ["Times New Roman", "Times", "STIXGeneral", "DejaVu Serif"],
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "stix" if j.family == "serif" else "dejavusans",
        "font.size": j.size, "axes.labelsize": j.size, "axes.titlesize": j.size,
        "xtick.labelsize": j.size, "ytick.labelsize": j.size, "legend.fontsize": j.size - 1,
        "axes.linewidth": 0.6, "lines.linewidth": 1.0, "patch.linewidth": 0.6,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
        "axes.grid": False, "legend.frameon": False, "legend.handlelength": 2.6,
        "axes.edgecolor": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
        "text.color": "black", "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",  # text stays editable text
    }


def _style(i: int) -> dict:
    return {"color": COLORS[i % len(COLORS)], "linestyle": DASHES[i % len(DASHES)]}


# -- the plotted numbers ----------------------------------------------------------------------------------------------


def series(planner, name: str, *, detector=None, events: int = 100_000, seed: int = 1, sweep=None) -> dict:
    """What a figure plots: ``{"x": label, "y": label, "series": [(label, x, y), ...]}`` (``"z"`` and a third
    array for the 3D geometry). The strips figure has one row per strip instead: x and y are the strip indices and
    the third array is the rate."""
    if name == "kinematics":
        out = []
        for c in planner.kinematics()["curves"]:
            th, en = np.asarray(c["theta"]), np.asarray(c["energy"])
            ok = th <= c["max_angle"] + 1e-9  # a recoil left at rest has no direction
            out.append((c["label"], th[ok], en[ok]))
        return {"x": "lab angle θ (deg)", "y": "energy (MeV)", "series": out}
    if name == "energy_loss":
        el = planner.energy_loss()
        return {"x": "depth (mg/cm²)", "y": "beam energy (MeV)",
                "series": [("beam", np.asarray(el["depth_mg_cm2"]), np.asarray(el["beam_energy_mev"]))]}
    if name == "spectra":
        out = []
        for det, h in planner.spectra(events=events, seed=seed)["spectra"].items():
            e = np.asarray(h["edges"])
            out.append((det, 0.5 * (e[1:] + e[:-1]), np.asarray(h["counts"])))
        return {"x": "measured energy (MeV)", "y": "counts / bin", "series": out}
    if name == "trajectories":
        t = planner.trajectories()
        out = [(f"b = {o['b_fm']:.1f} fm, {o['deflection_deg']:.0f}°", o["xy"][:, 0], o["xy"][:, 1])
               for o in t["orbits"]]
        return {"x": "x (fm)", "y": "y (fm)", "series": out}
    if name == "gamma":
        g = planner.gamma()
        out = [("P(θ)", np.asarray(g["theta_cm"]), np.asarray(g["probability"]))] if g["available"] else []
        return {"x": "CM angle θ (deg)", "y": "excitation probability", "series": out}
    if name == "sweep":
        if sweep is None:
            raise ValueError("the sweep figure needs sweep=Planner.sweep(...)")
        return {"x": sweep["parameter"], "y": f"{sweep['quantity']} ({sweep['detector']})",
                "series": [(sweep["quantity"], np.asarray(sweep["x"], float), np.asarray(sweep["y"], float))]}
    if name == "strips":
        strips = planner.rates()["strips"]
        detector = detector or next(iter(strips))
        if detector not in strips:
            raise ValueError(f"no detector {detector!r}; choose one of {', '.join(strips)}")
        ij = np.array(sorted(strips[detector]))
        return {"x": "strip / ring", "y": "strip / sector", "z": "rate (1/s)",
                "series": [(detector, ij[:, 0], ij[:, 1],
                            np.array([strips[detector][tuple(k)] for k in ij], float))]}
    if name == "geometry":
        g = planner.geometry()
        t = g["target_size_mm"]
        out = [("beam", *g["beam"].T),
               ("target", *np.array([[-t, -t, 0], [t, -t, 0], [t, t, 0], [-t, t, 0], [-t, -t, 0]], float).T)]
        out += [(d["name"], *d["outline"].T) for d in g["detectors"]]
        out += [(f"{d['name']} (γ)", *d["outline"].T) for d in g.get("gamma_detectors", [])]
        return {"x": "x (mm)", "y": "y (mm)", "z": "z, beam (mm)", "series": out}
    raise ValueError(f"unknown figure {name!r}; choose one of {', '.join(FIGURES)}")


def data_csv(planner, name: str, **options) -> str:
    """The numbers a figure plots, as CSV text: one row per point, with the curve's label in the first column."""
    s = series(planner, name, **options)
    cols = [s["x"], s["y"]] + ([s["z"]] if "z" in s else [])
    lines = [",".join(["series"] + [f'"{c}"' for c in cols])]
    for label, *arrays in s["series"]:
        for row in zip(*arrays):
            lines.append(",".join([f'"{label}"'] + [f"{float(v):.8g}" for v in row]))
    return "\n".join(lines) + "\n"


# -- drawing ----------------------------------------------------------------------------------------------------------


def _draw(planner, name: str, fig, options: dict) -> None:
    s = series(planner, name, **options)
    if name == "geometry":
        ax = fig.add_subplot(projection="3d")
        e = planner.geometry()["extent"]
        n_particle = len(planner.geometry()["detectors"])
        for k, (label, x, y, z) in enumerate(s["series"]):
            if label == "beam":
                ax.plot(x, y, z, color="0.45", lw=0.8, ls="--")
            elif label == "target":
                ax.plot(x, y, z, color="black", lw=1.2)
            else:
                i = k - 2
                st = _style(i) if i < n_particle else {"color": "0.25", "linestyle": ":"}
                ax.plot(x, y, z, lw=1.0, **st)
                c = np.array([x.mean(), y.mean(), z.mean()]) * 1.12
                ax.text(*c, label, color=st["color"], ha="center", va="center")
        ax.set(xlim=(-e, e), ylim=(-e, e), zlim=(-e, e), xlabel=s["x"], ylabel=s["y"], zlabel=s["z"])
        ax.set_box_aspect((1, 1, 1))
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor("white")
            axis.pane.set_edgecolor("0.6")
            axis._axinfo["grid"].update(color="0.88", linewidth=0.4)  # noqa: SLF001 -- no public setting
        ax.tick_params(pad=0)
        ax.xaxis.labelpad = ax.yaxis.labelpad = ax.zaxis.labelpad = -4
        return
    ax = fig.add_subplot()
    ax.set(xlabel=s["x"], ylabel=s["y"])
    if name == "strips":
        _, i, j, rate = s["series"][0]
        z = np.zeros((int(j.max()) + 1, int(i.max()) + 1))
        z[j.astype(int), i.astype(int)] = rate
        mesh = ax.pcolormesh(np.arange(z.shape[1] + 1) - 0.5, np.arange(z.shape[0] + 1) - 0.5, z, cmap="viridis",
                             rasterized=True)
        fig.colorbar(mesh, ax=ax, label=s["z"])
        return
    for k, (label, x, y) in enumerate(s["series"]):
        if name == "spectra":
            ax.stairs(np.where(y > 0, y, np.nan), _edges(x), lw=0.8, label=label, **_style(k))
        elif name == "sweep":
            ax.plot(x, y, marker="o", ms=3, color="black", label=label)
        else:
            ax.plot(x, y, label=label, **_style(k))
    if name == "kinematics":
        for i, d in enumerate(planner.kinematics()["detectors"]):
            lo, hi = d["theta_range"]
            c = COLORS[i % len(COLORS)]
            ax.axvspan(lo, hi, color=c, alpha=0.15, lw=0)
            ax.text(0.5 * (lo + hi), 1.01, d["name"], transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    color=c, fontsize="small")
        ax.set(xlim=(0, 180), xticks=range(0, 181, 30))
        ax.set_ylim(bottom=0)
        ax.legend()
    elif name == "spectra":
        ax.set_yscale("log")
        ax.legend()
    elif name == "trajectories":
        t = planner.trajectories()
        a = np.linspace(0, 2 * math.pi, 200)
        r = t["interaction_radius_fm"]
        ax.plot(r * np.cos(a), r * np.sin(a), color="0.4", ls=":", lw=0.8, label="nuclear range")
        lim = 12 * t["d0_fm"]
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal")
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
    elif name == "gamma":
        g = planner.gamma()
        ax.set(xlim=(0, 180), xticks=range(0, 181, 30))
        if g["available"]:
            ax.set_yscale("log")
            if g["max_safe_angle"] < 180:
                ax.axvspan(g["max_safe_angle"], 180, facecolor="none", edgecolor="0.4", hatch="////", lw=0)
                ax.text(0.5 * (g["max_safe_angle"] + 180), 1.01, "not safe", transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize="small")
    if name != "trajectories" and ax.get_yscale() != "log":
        ax.minorticks_on()


def _edges(centres) -> np.ndarray:
    half = 0.5 * (centres[1] - centres[0]) if len(centres) > 1 else 0.5
    return np.concatenate([centres - half, [centres[-1] + half]])


def figure(planner, name: str, journal="physical_review", width="single", **options):
    """One of the planner's figures (a key of :data:`FIGURES`) in a journal's style, as a matplotlib figure of the
    journal's column width. ``options``: ``detector`` for "strips", ``events`` and ``seed`` for "spectra",
    ``sweep`` (what :meth:`Planner.sweep` returned) for "sweep"."""
    import matplotlib
    from matplotlib.figure import Figure

    if name not in FIGURES:
        raise ValueError(f"unknown figure {name!r}; choose one of {', '.join(FIGURES)}")
    w = width_mm(journal, width) * MM
    with matplotlib.rc_context(rc(journal)):
        fig = Figure(figsize=(w, w * _ASPECT.get(name, 0.68)), layout="constrained")
        _draw(planner, name, fig, options)
    return fig


def export(planner, name: str, fmt: str = "pdf", journal="physical_review", width="single", **options) -> bytes:
    """A figure's file contents: ``fmt`` is "pdf", "svg" (both vector, with editable text) or "png" (600 dpi)."""
    import matplotlib

    if fmt not in FORMATS:
        raise ValueError(f"unknown format {fmt!r}; choose one of {', '.join(FORMATS)}")
    fig = figure(planner, name, journal, width, **options)
    buf = io.BytesIO()
    with matplotlib.rc_context(rc(journal)):
        fig.savefig(buf, format=fmt, **FORMATS[fmt][1])
    return buf.getvalue()


def preview(planner, name: str, journal="physical_review", width="single", dpi: int = 200, **options) -> bytes:
    """The figure as a PNG for the screen."""
    import matplotlib

    fig = figure(planner, name, journal, width, **options)
    buf = io.BytesIO()
    with matplotlib.rc_context(rc(journal)):
        fig.savefig(buf, format="png", dpi=dpi)
    return buf.getvalue()


__all__ = ["COLORS", "FIGURES", "FORMATS", "JOURNALS", "WIDTHS", "Journal", "data_csv", "export", "figure",
           "journal", "preview", "rc", "series", "width_mm"]
