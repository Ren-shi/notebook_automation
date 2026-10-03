import math

import numpy as np
import pytest

import physim as ps
from physim import scenarios as sc
from physim import units as u


def test_units_round_trip_and_astronomical_g():
    au = u.ASTRONOMICAL
    assert au.G == pytest.approx(4 * math.pi**2, rel=1e-15)
    # G = 4π² in AU³/(M_sun yr²) is the SI G converted.
    assert u.convert(u.G_SI, "G", u.SI, au) == pytest.approx(4 * math.pi**2, rel=1e-12)
    assert au.to_si(1.0, "time") / u.DAY == pytest.approx(365.2568983, rel=1e-9)
    assert au.mass == pytest.approx(1.989e30, rel=1e-3)
    v = au.from_si(29.78e3, "velocity")
    assert v == pytest.approx(2 * math.pi, rel=2e-3)
    x = np.array([1.0, 2.0])
    assert np.allclose(au.from_si(au.to_si(x, (2, 1, -2)), "energy"), x)
    custom = u.scaled("km-kg-s", 1e3, 1.0, 1.0)
    assert custom.G == pytest.approx(u.G_SI * 1e-9)
    assert u.convert(1.0, "length", custom, u.SI) == 1e3
    with pytest.raises(ValueError, match="no physical scale"):
        u.DIMENSIONLESS.to_si(1.0, "length")
    with pytest.raises(ValueError, match="unknown quantity"):
        u.SI.scale("speed")


def test_elements_round_trip():
    for el in [
        dict(a=1.3, e=0.4, i=0.3, Omega=1.0, omega=2.0, nu=0.7),
        dict(a=5.0, e=0.05, i=2.9, Omega=4.0, omega=0.1, nu=5.5),
    ]:
        r, v = sc.elements_to_state(el["a"], el["e"], el["i"], el["Omega"], el["omega"], nu=el["nu"], mu=2.0)
        back = sc.state_to_elements(r, v, mu=2.0)
        for key, val in el.items():
            assert back[key] == pytest.approx(val, abs=1e-10), key
        r2, _ = sc.elements_to_state(el["a"], el["e"], el["i"], el["Omega"], el["omega"], M=back["M"], mu=2.0)
        assert np.allclose(r, r2, atol=1e-12)
    with pytest.raises(ValueError):
        sc.elements_to_state(1.0, 1.2)
    with pytest.raises(ValueError):
        sc.elements_to_state(1.0, 0.1, nu=0.0, M=0.0)


def test_two_body_period_and_energy():
    m1, m2, a, e = 1.0, 0.3, 2.0, 0.6
    w = sc.two_body(m1, m2, a, e, i=0.4, Omega=0.2, omega=1.0, M=0.0)
    T = 2 * math.pi * math.sqrt(a**3 / (m1 + m2))
    mu = m1 * m2 / (m1 + m2)
    assert w.total_energy() == pytest.approx(-m1 * m2 / (2 * a), rel=1e-12)
    assert np.allclose(w.momentum(), 0, atol=1e-15)
    r0 = w.positions[1] - w.positions[0]
    n = 20_000
    w.run(T / n, n, record_every=n)
    assert np.allclose(w.positions[1] - w.positions[0], r0, atol=1e-6)
    assert mu > 0


def test_solar_system():
    w = sc.solar_system()
    assert w.n_particles == 9
    assert w.masses[5] == pytest.approx(1 / 1047.3486)
    assert np.allclose(w.momentum(), 0, atol=1e-15)
    # Earth on 2000-01-01 is ~0.983 AU from the Sun (just before perihelion).
    assert np.linalg.norm(w.positions[3] - w.positions[0]) == pytest.approx(0.9833, abs=1e-3)
    # Jupiter's period from its elements: 11.86 years.
    rel = sc.state_to_elements(w.positions[5] - w.positions[0], w.velocities[5] - w.velocities[0],
                               mu=4 * math.pi**2 * (1 + w.masses[5]))
    assert rel["period"] == pytest.approx(11.86, abs=0.02)
    # Outer planets for 1000 years with Wisdom-Holman at 0.5 yr steps.
    w = sc.solar_system(["Jupiter", "Saturn", "Uranus", "Neptune"])
    e0 = w.total_energy()
    tr = w.run(0.5, 2000, record_every=100)
    assert np.max(np.abs(tr.energy / e0 - 1)) < 1e-5
    # SI units give the same physics.
    w_si = sc.solar_system(["Earth"], units=u.SI, integrator="yoshida4")
    assert np.linalg.norm(w_si.positions[1] - w_si.positions[0]) == pytest.approx(0.9833 * u.AU, rel=1e-3)
    year = u.ASTRONOMICAL.time
    w_si.run(year / 2000, 2000, record_every=2000)
    w_au = sc.solar_system(["Earth"], integrator="yoshida4")
    w_au.run(1 / 2000, 2000, record_every=2000)
    assert np.allclose(w_si.positions / u.AU, w_au.positions, atol=1e-9)
    with pytest.raises(ValueError, match="unknown planet"):
        sc.solar_system(["Vulcan"])


def test_plummer_sphere_is_virialised():
    n = 1000
    w = sc.plummer_sphere(n, seed=3, softening=0.01)
    assert w.n_particles == n
    assert np.allclose(w.center_of_mass(), 0, atol=1e-12)
    assert np.allclose(w.momentum(), 0, atol=1e-12)
    # Total energy -3π/64 for G = M = a = 1 (up to sampling noise and the r < 10a cut).
    assert w.total_energy() == pytest.approx(-3 * math.pi / 64, rel=0.08)
    # Virial ratio 2K/|U| ≈ 1, averaged over a crossing time.
    tr = w.run(0.01, 200, record_every=10)
    q = np.mean(2 * tr.kinetic / np.abs(tr.potential))
    assert q == pytest.approx(1.0, abs=0.06)
    # Reproducible from the seed.
    a, b = sc.plummer_sphere(50, seed=7), sc.plummer_sphere(50, seed=7)
    assert np.array_equal(a.positions, b.positions)


def test_figure_eight_returns_after_one_period():
    w = sc.figure_eight()
    x0, v0 = w.positions.copy(), w.velocities.copy()
    assert np.allclose(w.momentum(), 0, atol=1e-8)
    n = 8000
    w.run(sc.FIGURE_EIGHT_PERIOD / n, n, record_every=n)
    assert np.max(np.abs(w.positions - x0)) < 1e-6
    assert np.max(np.abs(w.velocities - v0)) < 1e-6
    # Heavier masses: faster by sqrt(G m).
    w = sc.figure_eight(mass=4.0)
    w.run(sc.FIGURE_EIGHT_PERIOD / 2 / n, n, record_every=n)
    assert np.max(np.abs(w.positions - x0)) < 1e-6


def test_pendulum_period():
    theta0, L, g = 2.0, 1.5, 9.81
    w = sc.pendulum(L, theta0, g=g)
    # Exact period 4 sqrt(L/g) K(sin²(θ0/2)), via the AGM.
    a, b = 1.0, math.cos(theta0 / 2)
    while abs(a - b) > 1e-15:
        a, b = (a + b) / 2, math.sqrt(a * b)
    T = 2 * math.pi * math.sqrt(L / g) / a
    ev = ps.Event.coordinate(0, "x", direction=1)
    tr = w.run(1e-4, int(1.6 * T / 1e-4), record_every=1000, events=[ev])
    assert len(tr.event_t) >= 1
    # Crossings of x = 0 moving right: one per period, first at 3T/4.
    assert tr.event_t[0] == pytest.approx(0.75 * T, rel=1e-6)
    assert np.allclose(np.linalg.norm(tr.pos[:, 0], axis=1), L, atol=1e-9)


def test_spring_chain_mode_frequency():
    n, k, m = 12, 4.0, 0.5
    w = sc.spring_lattice((n,), k=k, mass=m, boundary="fixed")
    assert w.pinned[0] and w.pinned[-1] and not w.pinned[1:-1].any()
    # Excite the lowest longitudinal mode and time a full oscillation.
    j = np.arange(n)
    shape = np.sin(math.pi * j / (n - 1))
    pos = w.positions.copy()
    pos[:, 0] += 0.01 * shape
    w.positions = pos
    omega = 2 * math.sqrt(k / m) * math.sin(math.pi / (2 * (n - 1)))
    T = 2 * math.pi / omega
    steps = 20_000
    w.run(T / steps, steps, record_every=steps)
    assert np.allclose(w.positions[:, 0] - j, 0.01 * shape, atol=1e-6)


def test_spring_lattices():
    w = sc.spring_lattice((3, 4), boundary="fixed", diagonals=True, spacing=2.0)
    assert w.n_particles == 12
    assert w.potential_energy() == 0.0
    # Edges are pinned: only the two interior nodes move.
    assert np.count_nonzero(~w.pinned) == 2
    # 17 nearest-neighbour bonds and 12 diagonals of length 2√2.
    rest = np.asarray(w.checkpoint()["forces"][0]["rest_length"])
    assert rest.size == 29 and np.count_nonzero(np.isclose(rest, 2 * math.sqrt(2))) == 12
    # 3D cube of side 3: 3 n² (n - 1) = 54 nearest-neighbour bonds.
    w = sc.spring_lattice((3, 3, 3))
    assert list(w.forces.values()) == ["SpringNetwork(54 bonds)"]
    with pytest.raises(ValueError):
        sc.spring_lattice((2, 2), boundary="periodic")
