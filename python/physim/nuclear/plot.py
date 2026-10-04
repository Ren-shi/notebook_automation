"""Pictures of a nuclear experiment setup (matplotlib, imported on first use).

::

    from physim.nuclear import Experiment, plot

    exp = Experiment.example("oxygen_on_lead_array")
    plot.setup_3d(exp)       # beam, target and detectors in 3D
    plot.coverage(exp)       # where each detector sits in (θ, φ)
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
