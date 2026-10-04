"""Reproducible random numbers and random initial conditions.

Everything here draws from physim's counter-based generator (the one the Langevin thermostat uses): each value is
a pure function of ``(seed, counter, index)``. The same seed always gives the same numbers, on any machine and for
any thread count. Setup helpers use the counters from :data:`SETUP` upwards (``SETUP + draw``), so they never
reuse the numbers a thermostat draws during a run::

    w = ps.World()
    w.seed = 42
    pos = ps.random.positions_in_sphere(500, radius=1.0, seed=w.seed)
    for p in pos:
        w.add_particle(p, mass=1.0)
    ps.random.maxwell_boltzmann(w, temperature=0.5)   # sets velocities, total momentum zero

Pass a different ``draw`` to get a fresh, independent set from the same seed.
"""

from __future__ import annotations

import numpy as np

from ._core import RANDOM_SETUP_COUNTER, random_normal, random_uniform

#: First counter used by setup helpers; runs use counters 0, 1, 2, ... (one per step).
SETUP = RANDOM_SETUP_COUNTER


def uniform(n: int, seed: int = 1, counter: int = 0, start: int = 0) -> np.ndarray:
    """``n`` uniform deviates in (0, 1): values ``start, ..., start + n - 1`` of draw ``counter``."""
    return random_uniform(seed, counter, n, start)


def normal(n: int, seed: int = 1, counter: int = 0, start: int = 0) -> np.ndarray:
    """``n`` standard normal deviates, as :func:`uniform`."""
    return random_normal(seed, counter, n, start)


def _setup_counter(draw: int) -> int:
    if draw < 0:
        raise ValueError("draw must be a non-negative integer")
    return SETUP + int(draw)


def positions_in_box(n: int, low, high, seed: int = 1, draw: int = 0) -> np.ndarray:
    """``n`` points uniformly in the box between corners ``low`` and ``high``, an (n, 3) array."""
    low, high = np.asarray(low, dtype=float), np.asarray(high, dtype=float)
    u = uniform(3 * n, seed, _setup_counter(draw)).reshape(n, 3)
    return low + u * (high - low)


def positions_in_sphere(n: int, radius: float, centre=(0.0, 0.0, 0.0), seed: int = 1, draw: int = 0) -> np.ndarray:
    """``n`` points uniformly inside a ball: isotropic directions, radius R u^(1/3)."""
    counter = _setup_counter(draw)
    g = normal(3 * n, seed, counter).reshape(n, 3)
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    r = radius * uniform(n, seed, counter, start=6 * n) ** (1.0 / 3.0)  # indices beyond those the normals used
    return np.asarray(centre, dtype=float) + g * r[:, None]


def maxwell_boltzmann(world, temperature: float, *, seed=None, draw: int = 0, zero_momentum: bool = True,
                      k_b: float = 1.0) -> np.ndarray:
    """Set the velocities of ``world`` from the Maxwell-Boltzmann distribution at ``temperature`` and return them.

    Each component is normal with variance k_B T / m. Pinned particles keep zero velocity. With
    ``zero_momentum`` the centre-of-mass velocity is removed (which lowers the temperature by a factor
    (N - 1)/N). Draws from ``seed`` (the world's :attr:`~physim.World.seed` if not given).
    """
    if temperature < 0:
        raise ValueError("temperature must not be negative")
    seed = world.seed if seed is None else seed
    m = np.asarray(world.masses, dtype=float)
    pinned = np.asarray(world.pinned, dtype=bool)
    n = len(m)
    xi = normal(3 * n, seed, _setup_counter(draw)).reshape(n, 3)
    with np.errstate(divide="ignore"):
        sigma = np.where(pinned | (m <= 0), 0.0, np.sqrt(k_b * temperature / np.where(m > 0, m, 1.0)))
    v = xi * sigma[:, None]
    free = ~pinned
    if zero_momentum and free.sum() > 1:
        p = (m[free, None] * v[free]).sum(axis=0)
        v[free] -= p / m[free].sum()
    world.velocities = v
    return v
