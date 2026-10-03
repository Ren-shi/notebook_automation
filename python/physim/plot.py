"""Plotting and animation helpers (matplotlib, imported on first use).

Every function takes a recorded trajectory (:class:`physim.Trajectory` or
:class:`physim.RigidTrajectory`), draws on ``ax`` (a new figure if ``None``) and returns
the axes, so the result can be styled further::

    ax = ps.plot.orbits(traj, labels=["Sun", "Earth"])
    ps.plot.energy_error({"verlet": t1, "rk4": t2})
    anim = ps.plot.animate(traj, trail=40)
    anim.save("orbit.gif", writer="pillow")      # or "orbit.mp4" with ffmpeg installed

:func:`view3d` makes an interactive plotly figure; it needs the optional ``plotly`` package
(``pip install physim[plot3d]``).
"""

from __future__ import annotations

import numpy as np

from .analysis import relative_energy_error

#: Default colour cycle: distinguishable, and legible on white and dark backgrounds.
COLORS = [
    "#2a78d6",
    "#c4482b",
    "#52514e",
    "#d68a2a",
    "#8a5ad6",
    "#2a9d8f",
    "#c2408f",
    "#7a8b2a",
]

_AXES = {"x": 0, "y": 1, "z": 2}


def _plt():
    import matplotlib.pyplot as plt

    return plt


def _axes(ax, figsize):
    if ax is None:
        _, ax = _plt().subplots(figsize=figsize)
    return ax


def _plane(plane):
    if len(plane) != 2 or any(c not in _AXES for c in plane) or plane[0] == plane[1]:
        raise ValueError(f"plane must be two of 'x', 'y', 'z', e.g. 'xy'; got {plane!r}")
    return _AXES[plane[0]], _AXES[plane[1]]


def _particles(traj, particles):
    n = np.asarray(traj.pos).shape[1]
    if particles is None:
        return list(range(n))
    particles = [int(p) for p in np.atleast_1d(particles)]
    for p in particles:
        if not 0 <= p < n:
            raise ValueError(f"particle {p} out of range ({n} particles)")
    return particles


def color(i):
    """The ``i``-th colour of :data:`COLORS` (cycling)."""
    return COLORS[i % len(COLORS)]


def tidy(ax, grid=False, legend=None):
    """House style: no top/right spines, optional light grid, and a frameless legend when
    an artist is labelled and the axes have no legend yet (or when ``legend=True``).
    Returns ``ax``."""
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(alpha=0.25)
    labelled = ax.get_legend_handles_labels()[1]
    if legend or (legend is None and labelled and ax.get_legend() is None):
        ax.legend(frameon=False)
    return ax


def orbits(traj, particles=None, *, ax=None, plane="xy", labels=None, colors=None, start=True,
           lw=1.0, equal=True, figsize=(5, 5)):
    """Paths of ``particles`` (default all) projected on ``plane``, with a dot at each
    starting point when ``start``. ``labels`` and ``colors`` are per-particle lists."""
    ax = _axes(ax, figsize)
    a, b = _plane(plane)
    pos = np.asarray(traj.pos)
    for k, p in enumerate(_particles(traj, particles)):
        c = colors[k] if colors is not None else color(k)
        label = labels[k] if labels is not None else None
        ax.plot(pos[:, p, a], pos[:, p, b], color=c, lw=lw, label=label)
        if start:
            ax.plot(pos[0, p, a], pos[0, p, b], "o", color=c, ms=4)
    if equal:
        ax.set_aspect("equal")
    ax.set_xlabel(plane[0])
    ax.set_ylabel(plane[1])
    return tidy(ax)


def energy_error(trajs, *, ax=None, labels=None, floor=1e-17, time_unit=1.0, xlabel="t",
                 figsize=(7, 3.5)):
    """``|ΔE / E₀|`` against time on a log scale. ``trajs`` is one trajectory, a list, or a
    ``{label: trajectory}`` dict (the labels become the legend). Times are divided by
    ``time_unit`` (e.g. 365.25 to plot a run in days against years)."""
    ax = _axes(ax, figsize)
    if isinstance(trajs, dict):
        labels, trajs = list(trajs), list(trajs.values())
    elif not isinstance(trajs, (list, tuple)):
        trajs = [trajs]
    for k, tr in enumerate(trajs):
        if getattr(tr, "energy", None) is None:
            raise ValueError("the trajectory has no energies (run with energies=True)")
        err = np.abs(relative_energy_error(tr)) + floor
        t = np.asarray(tr.t) / time_unit
        ax.semilogy(t, err, color=color(k), lw=1.2, label=labels[k] if labels else None)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$|\Delta E / E_0|$")
    return tidy(ax)


def phase_space(traj, particle=0, axis="x", *, ax=None, momentum=False, mass=None, lw=0.8,
                figsize=(5, 4), **kwargs):
    """Velocity (or momentum ``m v`` with ``momentum=True`` and the particle's ``mass``)
    against position along ``axis`` for one particle."""
    ax = _axes(ax, figsize)
    if axis not in _AXES:
        raise ValueError(f"axis must be 'x', 'y' or 'z', got {axis!r}")
    (p,) = _particles(traj, particle)
    k = _AXES[axis]
    q = np.asarray(traj.pos)[:, p, k]
    v = np.asarray(traj.vel)[:, p, k]
    if momentum:
        if mass is None:
            raise ValueError("momentum=True needs the particle's mass")
        v = mass * v
    kwargs.setdefault("color", color(0))
    ax.plot(q, v, lw=lw, **kwargs)
    ax.set_xlabel(axis)
    ax.set_ylabel(f"p_{axis}" if momentum else f"v_{axis}")
    return tidy(ax)


def animate(traj, particles=None, *, plane="xy", trail=30, every=1, interval=40, sizes=None,
            colors=None, labels=None, ax=None, limits=None, time_label=True, figsize=(5, 5)):
    """A matplotlib animation of the particles moving on ``plane``, each with a fading tail
    of the last ``trail`` frames (0 for none). ``every`` skips frames; ``interval`` is the
    delay between frames in ms. Save with ``anim.save("out.gif", writer="pillow")`` or
    ``"out.mp4"`` (ffmpeg); in a notebook show it with ``IPython.display.HTML(anim.to_jshtml())``.

    ``limits`` is ``((x0, x1), (y0, y1))`` (default: fits the whole trajectory).
    """
    plt = _plt()
    from matplotlib.animation import FuncAnimation

    a, b = _plane(plane)
    ps_ = _particles(traj, particles)
    pos = np.asarray(traj.pos)[::every][:, ps_]
    t = np.asarray(traj.t)[::every]
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    if limits is None:
        lo = pos[:, :, [a, b]].min(axis=(0, 1))
        hi = pos[:, :, [a, b]].max(axis=(0, 1))
        pad = 0.05 * np.maximum(hi - lo, 1e-12)
        limits = ((lo[0] - pad[0], hi[0] + pad[0]), (lo[1] - pad[1], hi[1] + pad[1]))
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal")
    ax.set_xlabel(plane[0])
    ax.set_ylabel(plane[1])
    tidy(ax, legend=False)
    cols = [colors[k] if colors is not None else color(k) for k in range(len(ps_))]
    size = np.full(len(ps_), 30.0) if sizes is None else np.broadcast_to(sizes, len(ps_)).astype(float)
    tails = [ax.plot([], [], color=c, lw=1, alpha=0.5)[0] for c in cols] if trail else []
    dots = ax.scatter(pos[0, :, a], pos[0, :, b], s=size, c=cols, zorder=3)
    if labels is not None:
        for k, lab in enumerate(labels):
            ax.plot([], [], "o", color=cols[k], label=lab)
        ax.legend(frameon=False, loc="upper right")
    text = ax.text(0.02, 0.97, "", transform=ax.transAxes, va="top") if time_label else None

    def update(f):
        dots.set_offsets(pos[f][:, [a, b]])
        for k, line in enumerate(tails):
            seg = pos[max(0, f - trail):f + 1, k]
            line.set_data(seg[:, a], seg[:, b])
        if text is not None:
            text.set_text(f"t = {t[f]:.3g}")
        return [dots, *tails] + ([text] if text is not None else [])

    anim = FuncAnimation(fig, update, frames=len(t), interval=interval, blit=True)
    plt.close(fig)  # show it via the animation, not as a static figure
    return anim


def view3d(traj, particles=None, *, labels=None, trail=True, marker_size=4, every=1):
    """An interactive 3D plotly figure: each particle's path (``trail``) and a slider over
    frames moving the markers. Needs ``plotly`` (``pip install physim[plot3d]``)."""
    try:
        import plotly.graph_objects as go
    except ImportError as err:  # pragma: no cover - depends on the environment
        raise ImportError("view3d needs plotly: pip install plotly (or physim[plot3d])") from err
    ps_ = _particles(traj, particles)
    pos = np.asarray(traj.pos)[::every][:, ps_]
    t = np.asarray(traj.t)[::every]
    names = labels if labels is not None else [f"particle {p}" for p in ps_]
    data = []
    if trail:
        for k in range(len(ps_)):
            data.append(go.Scatter3d(x=pos[:, k, 0], y=pos[:, k, 1], z=pos[:, k, 2], mode="lines",
                                     line=dict(color=color(k), width=2), name=names[k]))
    marker = dict(size=marker_size, color=[color(k) for k in range(len(ps_))])

    def points(f):
        return go.Scatter3d(x=pos[f, :, 0], y=pos[f, :, 1], z=pos[f, :, 2], mode="markers",
                            marker=marker, name="positions", showlegend=False)

    data.append(points(0))
    m = len(data) - 1
    frames = [go.Frame(data=[points(f)], traces=[m], name=str(f)) for f in range(len(t))]
    steps = [dict(method="animate", label=f"{t[f]:.3g}",
                  args=[[str(f)], dict(mode="immediate", frame=dict(duration=0, redraw=True))])
             for f in range(len(t))]
    fig = go.Figure(data=data, frames=frames)
    fig.update_layout(scene=dict(aspectmode="data"), sliders=[dict(steps=steps, currentvalue=dict(prefix="t = "))],
                      margin=dict(l=0, r=0, t=30, b=0))
    return fig
