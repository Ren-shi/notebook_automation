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


def radial_distribution(positions, box, r_max=None, bins=100):
    """Radial distribution function g(r) of particles in a periodic box (minimum image).

    ``positions`` is ``(N, 3)`` or ``(F, N, 3)`` (frames are averaged); ``box`` the side
    lengths (a number or 3 values). Returns ``(r, g)`` at the bin centres up to ``r_max``
    (default half the smallest side). O(N²) per frame.
    """
    pos = np.asarray(positions, dtype=float)
    if pos.ndim == 2:
        pos = pos[None]
    box = np.broadcast_to(np.asarray(box, dtype=float), (3,))
    r_max = 0.5 * box.min() if r_max is None else r_max
    n = pos.shape[1]
    edges = np.linspace(0.0, r_max, bins + 1)
    counts = np.zeros(bins)
    for frame in pos:
        for i in range(n - 1):
            d = frame[i + 1 :] - frame[i]
            d -= box * np.round(d / box)
            counts += np.histogram(np.linalg.norm(d, axis=1), bins=edges)[0]
    density = n / np.prod(box)
    shells = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = 2.0 * counts / (len(pos) * n * density * shells)
    return 0.5 * (edges[1:] + edges[:-1]), g


def tabulate_pair(potential, r_min, r_max, n=2001):
    """Tabulates a vectorised pair potential ``potential(r)`` for :class:`physim.TabulatedPair`:
    returns ``(r_min, dr, V, dV)`` with the derivative from fourth-order finite differences."""
    r = np.linspace(r_min, r_max, n)
    dr = r[1] - r[0]
    h = 1e-4 * dr
    v = np.asarray(potential(r), dtype=float)
    dv = (
        -np.asarray(potential(r + 2 * h)) + 8 * np.asarray(potential(r + h))
        - 8 * np.asarray(potential(r - h)) + np.asarray(potential(r - 2 * h))
    ) / (12 * h)
    return float(r_min), float(dr), v.tolist(), dv.tolist()
