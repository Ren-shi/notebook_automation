"""Diagnostics computed from recorded trajectories."""

import numpy as np


def relative_energy_error(traj):
    """``(E(t) - E(0)) / |E(0)|`` for each recorded frame."""
    e = np.asarray(traj.energy)
    return (e - e[0]) / abs(e[0])


def poincare_section(world, dt, steps, event, *, coordinates=None):
    """Runs ``world`` for ``steps`` steps of ``dt`` and returns the states at every crossing of
    ``event`` (a non-terminal :class:`physim.Event`, e.g. ``Event.coordinate(0, "x",
    direction=1)``) as ``(t, pos, vel)`` arrays of shapes (K,), (K, N, 3), (K, N, 3).

    Crossings are located to ~1e-12 of a step, so the points lie on the section to that
    accuracy. ``coordinates``, if given, is a function ``(pos, vel) -> array`` applied to
    each crossing instead, e.g. ``lambda p, v: (p[:, 0, 1], v[:, 0, 1])`` for (y, p_y).
    """
    traj = world.run(dt, steps, record_every=max(steps, 1), events=[event], energies=False)
    if coordinates is not None:
        return coordinates(traj.event_pos, traj.event_vel)
    return traj.event_t, traj.event_pos, traj.event_vel
