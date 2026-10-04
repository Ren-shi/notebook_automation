"""Pictures of a nuclear experiment setup (matplotlib, imported on first use).

::

    from physim.nuclear import Experiment, plot

    exp = Experiment.example("oxygen_on_lead_array")
    plot.setup_3d(exp)       # beam, target and detectors in 3D
    plot.coverage(exp)       # where each detector sits in (θ, φ)

    from physim.nuclear.events import simulate
    ev = simulate(exp)
    plot.spectra(ev)         # measured-energy spectrum of each detector
    plot.theta_energy(ev)    # measured energy against lab angle
    plot.kinematics(Planner(exp).kinematics())   # E against θ for every ejectile and recoil
"""

from __future__ import annotations

import numpy as np

from ..plot import _plt, color, tidy
from .detectors import Array, angles


def setup_3d(experiment, *, ax=None, figsize=(7, 6), labels=True, strips=True):
    """Draw the beam axis, the target and every detector face in 3D (mm). Returns the axes."""
    plt = _plt()
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(projection="3d")
    array = Array.from_experiment(experiment)
    extent = max(np.linalg.norm(g.centre) for g in array) * 1.25
    ax.plot([0, 0], [0, 0], [-extent, extent], color="#888888", lw=1, ls="--")
    ax.quiver(0, 0, -extent, 0, 0, 0.35 * extent, color="#888888", arrow_length_ratio=0.2)
    t = 0.04 * extent
    ax.plot([-t, t, t, -t, -t], [-t, -t, t, t, -t], [0, 0, 0, 0, 0], color="#c9a227", lw=2)
    for i, g in enumerate(array):
        edge = g.outline()
        ax.plot(edge[:, 0], edge[:, 1], edge[:, 2], color=color(i), lw=1.6)
        if strips and g.shape == "rectangle":
            for k in range(1, g.strips_x):
                a = -g.width / 2 + k * g.width / g.strips_x
                p = g.centre + np.outer([-g.height / 2, g.height / 2], g.v) + a * g.u
                ax.plot(p[:, 0], p[:, 1], p[:, 2], color=color(i), lw=0.3, alpha=0.5)
        if labels:
            c = g.centre * 1.08
            ax.text(c[0], c[1], c[2], g.name, color=color(i), fontsize=8)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("z, beam (mm)")
    lim = (-extent, extent)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_zlim(lim)
    ax.set_title(experiment.title)
    return ax


def coverage(experiment, *, ax=None, figsize=(7, 3.8)):
    """Each detector's outline in (θ, φ): what it covers as seen from the target. Returns the axes."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    array = Array.from_experiment(experiment)
    for i, g in enumerate(array):
        if g.phi_range() == (-180.0, 180.0):  # a ring around the beam: a band between its θ edges
            lo, hi = g.theta_range()
            ax.plot([-180, 180], [lo, lo], color=color(i), lw=1.6, label=g.name)
            ax.plot([-180, 180], [hi, hi], color=color(i), lw=1.6)
            continue
        th, ph = angles(g.outline(360))
        _, ph0 = angles(g.centre)
        # Unwrap φ around the detector's own azimuth, and draw copies shifted by ±360° so a detector straddling
        # ±180° appears whole at both edges of the plot.
        ph = ph0 + (ph - ph0 + 180.0) % 360.0 - 180.0
        for k, shift in enumerate((0.0, -360.0, 360.0)):
            ax.plot(ph + shift, th, color=color(i), lw=1.6, label=g.name if k == 0 else None)
    ax.set_xlim(-180, 180)
    ax.set_ylim(180, 0)
    ax.set_xlabel("φ (deg)")
    ax.set_ylabel("θ (deg)")
    tidy(ax, grid=True, legend=True)
    return ax


def spectra(events, detectors=None, *, bins=300, log=True, by_particle=True, figsize=None):
    """Measured-energy spectrum of each detector (counts in the planned beam time per bin), one panel per
    detector; with ``by_particle`` the ejectiles and recoils of each channel are drawn separately. Returns the
    figure."""
    plt = _plt()
    names = list(events.detectors) if detectors is None else list(detectors)
    fig, axes = plt.subplots(len(names), 1, figsize=figsize or (7, 2.2 * len(names)), squeeze=False)
    for ax, name in zip(axes[:, 0], names):
        m = events.select(name)
        x = events["measured"][m]
        if not len(x):
            ax.text(0.5, 0.5, f"{name}: no counts", transform=ax.transAxes, ha="center")
            continue
        rng = (0.0, float(x.max()) * 1.03)
        h, edges = events.spectrum(name, bins=bins, range=rng)
        ax.stairs(h, edges, color="#444444", lw=0.8, label="all")
        if by_particle:
            k = 0
            for ch in events.channels:
                for particle in ("ejectile", "recoil"):
                    hp, _ = events.spectrum(name, bins=bins, range=rng, channel=ch, particle=particle)
                    if hp.sum() > 1e-3 * h.sum():
                        ax.stairs(hp, edges, color=color(k), lw=1.0, label=f"{ch} {particle}")
                        k += 1
        if log:
            ax.set_yscale("log")
        ax.set_ylabel("counts / bin")
        ax.set_title(name, fontsize=9, loc="left")
        tidy(ax, legend=by_particle)
    axes[-1, 0].set_xlabel("measured energy (MeV)")
    fig.tight_layout()
    return fig


def theta_energy(events, detector=None, *, bins=(180, 200), ax=None, figsize=(7, 4)):
    """Measured energy against lab angle for all counted particles (or one detector): the kinematic curves of
    each channel. Returns the axes."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    m = events.select(detector)
    th, e = events["theta"][m], events["measured"][m]
    w = events["weight"][m] * events.beam_time_s
    if len(th):
        h, xe, ye = np.histogram2d(th, e, bins=bins, weights=w)
        from matplotlib.colors import LogNorm

        h = np.where(h > 0, h, np.nan)
        mesh = ax.pcolormesh(xe, ye, h.T, norm=LogNorm(), cmap="viridis")
        plt.colorbar(mesh, ax=ax, label="counts / bin")
    ax.set_xlabel("lab angle θ (deg)")
    ax.set_ylabel("measured energy (MeV)")
    ax.set_title(detector or "all detectors", fontsize=9, loc="left")
    tidy(ax)
    return ax


def kinematics(kin, *, ax=None, figsize=(7, 4)):
    """Lab energy against lab angle for every ejectile and recoil, with each detector's angular coverage shaded.
    ``kin`` is :meth:`physim.nuclear.planner.Planner.kinematics`. Returns the axes."""
    plt = _plt()
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    for i, d in enumerate(kin["detectors"]):
        lo, hi = d["theta_range"]
        ax.axvspan(lo, hi, color=color(i), alpha=0.12, lw=0)
        ax.text(0.5 * (lo + hi), 1.0, d["name"], transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=7, color=color(i))
    for c in kin["curves"]:
        ok = np.asarray(c["theta"]) <= c["max_angle"] + 1e-9  # a recoil left at rest has no direction
        ax.plot(np.asarray(c["theta"])[ok], np.asarray(c["energy"])[ok], lw=1.4, ls="-" if c["particle"] == "ejectile" else "--", label=c["label"])
    ax.set_xlim(0, 180)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("lab angle θ (deg)")
    ax.set_ylabel("energy (MeV)")
    tidy(ax, legend=True)
    return ax
