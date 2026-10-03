"""Ready-made systems. Each function returns a configured :class:`physim.World`.

- :func:`two_body`: a binary or planet from orbital elements (barycentric frame).
- :func:`solar_system`: the Sun and planets from JPL's J2000 mean elements.
- :func:`plummer_sphere`: an equilibrium star cluster.
- :func:`figure_eight`: the Chenciner-Montgomery three-body choreography.
- :func:`pendulum`: a rigid pendulum (a rod to a fixed pivot under gravity).
- :func:`spring_lattice`: masses on a 1D, 2D or 3D grid joined by springs.

Orbital helpers: :func:`elements_to_state`, :func:`state_to_elements`. Angles are radians
unless a name says otherwise.
"""

from __future__ import annotations

import math

import numpy as np

from . import units as _units
from ._core import NewtonianGravity, SpringNetwork, UniformField, World

#: Period of the figure-eight orbit with G = m = 1 (Simó's value).
FIGURE_EIGHT_PERIOD = 6.3259139829

# JPL "Approximate positions of the planets" (Standish), Table 1 (1800-2050), J2000 ecliptic:
# a [AU], e, I [deg], L [deg], longitude of perihelion [deg], longitude of the node [deg],
# and the Sun/planet mass ratio. "Earth" is the Earth-Moon barycentre.
PLANETS = {
    "Mercury": (0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593, 6023600.0),
    "Venus": (0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255, 408523.71),
    "Earth": (1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0, 328900.56),
    "Mars": (1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891, 3098708.0),
    "Jupiter": (5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909, 1047.3486),
    "Saturn": (9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448, 3497.898),
    "Uranus": (19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503, 22902.98),
    "Neptune": (30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574, 19412.24),
    "Pluto": (39.48211675, 0.24882730, 17.14001206, 238.92903833, 224.06891629, 110.30393684, 1.35e8),
}


def _solve_kepler(M, e):
    """Eccentric anomaly from the mean anomaly (Newton's method)."""
    M = math.remainder(M, 2 * math.pi)
    E = M + e * math.sin(M) if e < 0.8 else math.pi * math.copysign(1.0, M)
    for _ in range(50):
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < 1e-15:
            break
    return E


def elements_to_state(a, e, i=0.0, Omega=0.0, omega=0.0, *, nu=None, M=None, mu=1.0):
    """Position and velocity of an elliptic orbit (``0 <= e < 1``) about a centre with
    gravitational parameter ``mu``: semi-major axis ``a``, inclination ``i``, longitude of
    the ascending node ``Omega``, argument of periapsis ``omega``, and either the true
    anomaly ``nu`` or the mean anomaly ``M`` (default: ``nu = 0``, at periapsis)."""
    if not (0 <= e < 1 and a > 0 and mu > 0):
        raise ValueError(f"need an ellipse: a > 0, 0 <= e < 1 and mu > 0 (got a={a}, e={e}, mu={mu})")
    if nu is not None and M is not None:
        raise ValueError("give the true anomaly nu or the mean anomaly M, not both")
    if M is not None:
        E = _solve_kepler(M, e)
        nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    nu = 0.0 if nu is None else nu
    p = a * (1 - e * e)
    r = p / (1 + e * math.cos(nu))
    pos = np.array([r * math.cos(nu), r * math.sin(nu), 0.0])
    s = math.sqrt(mu / p)
    vel = np.array([-s * math.sin(nu), s * (e + math.cos(nu)), 0.0])
    co, so = math.cos(omega), math.sin(omega)
    cO, sO = math.cos(Omega), math.sin(Omega)
    ci, si = math.cos(i), math.sin(i)
    R = np.array(
        [
            [cO * co - sO * so * ci, -cO * so - sO * co * ci, sO * si],
            [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
            [so * si, co * si, ci],
        ]
    )
    return R @ pos, R @ vel


def state_to_elements(pos, vel, mu=1.0):
    """Orbital elements of a bound orbit from a relative position and velocity: a dict with
    ``a, e, i, Omega, omega, nu, M, period``."""
    r = np.asarray(pos, float)
    v = np.asarray(vel, float)
    rn = np.linalg.norm(r)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    energy = 0.5 * v @ v - mu / rn
    if energy >= 0:
        raise ValueError("the orbit is not bound (energy >= 0)")
    a = -mu / (2 * energy)
    evec = np.cross(v, h) / mu - r / rn
    e = np.linalg.norm(evec)
    i = math.acos(np.clip(h[2] / hn, -1, 1))
    n = np.array([-h[1], h[0], 0.0])
    nn = np.linalg.norm(n)
    Omega = math.atan2(n[1], n[0]) if nn > 1e-14 * hn else 0.0
    # Reference directions in the orbital plane when the node or periapsis is undefined.
    node_dir = n / nn if nn > 1e-14 * hn else np.array([1.0, 0.0, 0.0])
    q = np.cross(h / hn, node_dir)
    if e > 1e-14:
        omega = math.atan2(evec @ q, evec @ node_dir)
        peri = evec / e
        nu = math.atan2(np.cross(peri, r) @ h / hn, peri @ r)
    else:
        omega = 0.0
        nu = math.atan2(r @ q, r @ node_dir)
    E = 2 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))
    return {
        "a": a,
        "e": e,
        "i": i,
        "Omega": Omega % (2 * math.pi),
        "omega": omega % (2 * math.pi),
        "nu": nu % (2 * math.pi),
        "M": (E - e * math.sin(E)) % (2 * math.pi),
        "period": 2 * math.pi * math.sqrt(a**3 / mu),
    }


def _to_barycentre(world):
    m = world.masses
    total = m.sum()
    com = (m[:, None] * world.positions).sum(0) / total
    vcom = (m[:, None] * world.velocities).sum(0) / total
    world.positions = world.positions - com
    world.velocities = world.velocities - vcom


def two_body(m1=1.0, m2=1e-3, a=1.0, e=0.0, i=0.0, Omega=0.0, omega=0.0, *, nu=None, M=None,
             G=1.0, integrator="yoshida4"):
    """Two bodies on a Keplerian orbit with the given elements (of body 2 relative to body 1),
    in the barycentric frame. The relative orbit has period ``2π sqrt(a³ / (G (m1 + m2)))``."""
    mu = G * (m1 + m2)
    r, v = elements_to_state(a, e, i, Omega, omega, nu=nu, M=M, mu=mu)
    w = World(integrator=integrator)
    w.add_particle([0, 0, 0], mass=m1)
    w.add_particle(r, vel=v, mass=m2)
    _to_barycentre(w)
    w.add_force(NewtonianGravity(G=G))
    return w


def solar_system(planets=None, *, units=None, integrator="wisdom_holman", barycentric=True):
    """The Sun (particle 0) and ``planets`` (names from :data:`PLANETS`, default the eight
    planets in order) at J2000, from JPL's mean elements (accurate to roughly 1e-3 in
    position over 1800-2050; "Earth" is the Earth-Moon barycentre).

    ``units`` defaults to :data:`physim.units.ASTRONOMICAL` (AU, solar masses, years,
    G = 4π²). With the default ``"wisdom_holman"`` integrator, take a step of about 1/20 of
    the innermost period (e.g. ``dt = 0.012`` with Mercury, ``0.6`` from Jupiter out).
    """
    units = _units.ASTRONOMICAL if units is None else units
    names = list(PLANETS)[:8] if planets is None else list(planets)
    for n in names:
        if n not in PLANETS:
            raise ValueError(f"unknown planet {n!r}; choose from {list(PLANETS)}")
    au = _units.ASTRONOMICAL
    G = au.G
    w = World(integrator=integrator)
    w.add_particle([0, 0, 0], mass=1.0)
    deg = math.pi / 180
    for n in names:
        a, e, inc, L, varpi, node, ratio = PLANETS[n]
        m = 1.0 / ratio
        r, v = elements_to_state(
            a, e, inc * deg, node * deg, (varpi - node) * deg, M=(L - varpi) * deg, mu=G * (1 + m)
        )
        w.add_particle(r, vel=v, mass=m)
    if barycentric:
        _to_barycentre(w)
    if units is not au:
        w.positions = _units.convert(w.positions, "length", au, units)
        w.velocities = _units.convert(w.velocities, "velocity", au, units)
        w.masses = _units.convert(w.masses, "mass", au, units)
        G = units.G
    w.add_force(NewtonianGravity(G=G))
    return w


def plummer_sphere(n, total_mass=1.0, scale_radius=1.0, *, G=1.0, seed=0, softening=0.0,
                   max_radius=None, integrator="verlet"):
    """``n`` equal-mass stars sampled from a Plummer sphere in equilibrium (Aarseth, Hénon &
    Wielen 1974), centred and at rest. The virial ratio ``2K / |U|`` is 1 up to sampling
    noise of order ``1 / sqrt(n)``; the total energy is ``-3π G M² / (64 a)``.

    ``max_radius`` (default ``10 a``) truncates the rare far-out stars; ``seed`` makes the
    sample reproducible; ``softening`` is passed to :class:`NewtonianGravity`.
    """
    if n < 2:
        raise ValueError("a Plummer sphere needs at least 2 stars")
    rng = np.random.default_rng(seed)
    a = float(scale_radius)
    r_max = 10 * a if max_radius is None else float(max_radius)
    # Radii from the cumulative mass m(r) = r³ / (r² + a²)^{3/2}, below r_max.
    x_max = (r_max**2 / (r_max**2 + a**2)) ** 1.5
    x = rng.uniform(0, x_max, n)
    r = a / np.sqrt(x ** (-2 / 3) - 1)
    pos = r[:, None] * _isotropic(rng, n)
    # Speeds q v_esc with q from g(q) = q² (1 - q²)^{7/2} by rejection (max g < 0.1).
    q = np.empty(n)
    filled = 0
    while filled < n:
        cand = rng.uniform(0, 1, 2 * (n - filled))
        y = rng.uniform(0, 0.1, cand.size)
        ok = cand[y < cand**2 * (1 - cand**2) ** 3.5]
        take = min(ok.size, n - filled)
        q[filled:filled + take] = ok[:take]
        filled += take
    v_esc = np.sqrt(2 * G * total_mass / np.sqrt(r**2 + a**2))
    vel = (q * v_esc)[:, None] * _isotropic(rng, n)
    w = World(integrator=integrator)
    m = total_mass / n
    for p, v in zip(pos, vel):
        w.add_particle(p, vel=v, mass=m)
    _to_barycentre(w)
    w.add_force(NewtonianGravity(G=G, softening=softening))
    return w


def _isotropic(rng, n):
    z = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * math.pi, n)
    s = np.sqrt(1 - z * z)
    return np.stack([s * np.cos(phi), s * np.sin(phi), z], axis=1)


def figure_eight(mass=1.0, *, G=1.0, integrator="yoshida4"):
    """Three equal masses chasing each other on a figure-eight (Chenciner & Montgomery 2000,
    initial conditions of Simó). With ``G = mass = 1`` the period is
    :data:`FIGURE_EIGHT_PERIOD`; in general it is that divided by ``sqrt(G mass)``."""
    s = math.sqrt(G * mass)
    x1 = np.array([0.97000436, -0.24308753, 0.0])
    v3 = np.array([-0.93240737, -0.86473146, 0.0]) * s
    w = World(integrator=integrator)
    w.add_particle(x1, vel=-v3 / 2, mass=mass)
    w.add_particle(-x1, vel=-v3 / 2, mass=mass)
    w.add_particle([0, 0, 0], vel=v3, mass=mass)
    w.add_force(NewtonianGravity(G=G))
    return w


def pendulum(length=1.0, theta0=0.5, *, omega0=0.0, g=9.81, mass=1.0, pivot=(0.0, 0.0, 0.0),
             integrator="verlet"):
    """A rigid pendulum in the x-y plane: a bob of ``mass`` on a massless rod of ``length``
    from the fixed ``pivot``, released at angle ``theta0`` from straight down with angular
    velocity ``omega0``, under gravity ``g`` along -y. The bob is particle 0. The small-angle
    period is ``2π sqrt(length / g)``."""
    pivot = np.asarray(pivot, float)
    d = np.array([math.sin(theta0), -math.cos(theta0), 0.0])
    t = np.array([math.cos(theta0), math.sin(theta0), 0.0])
    w = World(integrator=integrator)
    w.add_particle(pivot + length * d, vel=length * omega0 * t, mass=mass)
    w.add_rod(0, pivot, length)
    w.add_force(UniformField([0, -g, 0]))
    return w


def spring_lattice(shape, *, spacing=1.0, k=1.0, mass=1.0, boundary="free", diagonals=False,
                   integrator="verlet"):
    """Masses on a regular grid of ``shape`` (``(n,)``, ``(nx, ny)`` or ``(nx, ny, nz)``)
    joined to their nearest neighbours by springs of stiffness ``k`` at rest at ``spacing``,
    with particle index ``(i * ny + j) * nz + l``.

    ``boundary="fixed"`` pins the outermost layer (for a chain: both end masses), so the
    lowest longitudinal mode of a chain of ``n`` masses vibrates at
    ``2 sqrt(k / m) sin(π / (2 (n - 1)))``. ``diagonals=True`` adds face-diagonal springs
    (rest length ``spacing·√2``) so 2D and 3D lattices resist shear.
    """
    shape = tuple(int(s) for s in np.atleast_1d(shape))
    if not 1 <= len(shape) <= 3 or min(shape) < 1:
        raise ValueError(f"shape must have 1 to 3 positive sizes, got {shape}")
    if boundary not in ("free", "fixed"):
        raise ValueError(f"boundary must be 'free' or 'fixed', got {boundary!r}")
    dims = shape + (1,) * (3 - len(shape))
    idx = np.arange(np.prod(dims)).reshape(dims)
    grid = np.stack(np.meshgrid(*[np.arange(d) for d in dims], indexing="ij"), -1).reshape(-1, 3)
    w = World(integrator=integrator)
    for p in grid:
        w.add_particle(p * spacing, mass=mass)
    offsets = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    if diagonals:
        offsets += [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
    I, J, L = [], [], []
    for off in offsets:
        if any(o != 0 and d == 1 for o, d in zip(off, dims)):
            continue
        src = tuple(slice(max(0, -o), d - max(0, o)) for o, d in zip(off, dims))
        dst = tuple(slice(max(0, o), d + min(0, o)) for o, d in zip(off, dims))
        I.extend(idx[src].ravel())
        J.extend(idx[dst].ravel())
        L.extend([spacing * math.sqrt(sum(o * o for o in off))] * idx[src].size)
    if I:
        w.add_force(SpringNetwork(I, J, k, L))
    if boundary == "fixed":
        for p, n in zip(grid, idx.ravel()):
            if any((c == 0 or c == d - 1) and d > 1 for c, d in zip(p, dims)):
                w.pin(int(n))
    return w


__all__ = [
    "FIGURE_EIGHT_PERIOD",
    "PLANETS",
    "elements_to_state",
    "figure_eight",
    "pendulum",
    "plummer_sphere",
    "solar_system",
    "spring_lattice",
    "state_to_elements",
    "two_body",
]
