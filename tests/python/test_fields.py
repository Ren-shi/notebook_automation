import math

import numpy as np
import pytest

import physim as ps


def test_wave_equation_standing_mode():
    n, c = 65, 1.5
    h = 1.0 / (n - 1)
    w = ps.WaveEquation([n], h, c=c, boundary="dirichlet")
    (x,) = w.axes
    w.u = np.sin(math.pi * x)
    e0 = w.energy()
    steps = 400
    dt = 1.0 / steps
    assert dt < w.max_stable_dt
    w.step(dt, steps)
    assert w.t == pytest.approx(1.0)
    assert np.abs(w.u - np.sin(math.pi * x) * math.cos(math.pi * c)).max() < 1e-3
    assert w.energy() == pytest.approx(e0, rel=1e-3)
    with pytest.raises(ValueError, match="stability"):
        w.step(2 * w.max_stable_dt)
    with pytest.raises(ValueError, match="shape"):
        w.u = np.zeros(10)

    # 2D periodic membrane keeps its shape and frequency.
    n = 32
    m = ps.WaveEquation([n, n], 1.0 / n, boundary="periodic")
    X, Y = np.meshgrid(*m.axes, indexing="ij")
    m.u = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)
    assert m.u.shape == (n, n)
    period = 1 / math.sqrt(2)
    m.step(period / 200, 200)
    assert np.abs(m.u - np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)).max() < 0.02


def test_heat_equation_methods():
    n, d = 32, 0.5
    for method, dt in [("crank_nicolson", 0.01), ("explicit", None)]:
        s = ps.HeatEquation([n, n], 1.0 / n, diffusivity=d, boundary="neumann", method=method)
        X, Y = np.meshgrid(*s.axes, indexing="ij")
        s.u = np.cos(np.pi * X) * np.cos(np.pi * Y) + 2.0
        total = s.total()
        step = dt or 0.9 * s.max_explicit_dt
        steps = int(round(0.1 / step))
        s.step(0.1 / steps, steps)
        exact = np.cos(np.pi * X) * np.cos(np.pi * Y) * math.exp(-2 * np.pi**2 * d * 0.1) + 2.0
        assert np.abs(s.u - exact).max() < 2e-3, method
        assert s.total() == pytest.approx(total, rel=1e-10)
    assert s.max_explicit_dt == pytest.approx((1 / n) ** 2 / (4 * d))
    with pytest.raises(ValueError, match="stability"):
        s.step(10 * s.max_explicit_dt)
    with pytest.raises(ValueError, match="method"):
        ps.HeatEquation([8], 0.1, method="implicit")


def test_poisson():
    n = 32
    h = 1.0 / n
    x = np.arange(n) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    exact = np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y) * np.sin(4 * np.pi * Z)
    f = -(4 + 4 + 16) * np.pi**2 * exact
    phi = ps.solve_poisson(f, h, boundary="periodic")
    assert phi.shape == f.shape
    assert np.abs(phi - exact).max() < 0.03

    # Dirichlet with boundary values: a harmonic function is reproduced.
    n = 33
    x = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    harmonic = X**2 - Y**2 + 3 * X * Y
    phi = ps.solve_poisson(np.zeros((n, n)), x[1], boundary="dirichlet", phi=harmonic)
    assert np.abs(phi - harmonic).max() < 1e-9  # quadratics are exact for the 5-point stencil

    with pytest.raises(ValueError, match="powers of two"):
        ps.solve_poisson(np.zeros((12, 12)), 0.1)
    with pytest.raises(ValueError, match="boundary"):
        ps.solve_poisson(np.zeros(8), 0.1, boundary="open")


def test_particle_mesh_gravity():
    cluster = ps.scenarios.plummer_sphere(3000, seed=2, softening=0.0, max_radius=4.0)
    pm = ps.ParticleMesh(box_size=10.0, cells=64)
    direct = ps.World()
    for p, m in zip(cluster.positions, cluster.masses):
        direct.add_particle(p, mass=m)
    direct.add_force(ps.NewtonianGravity(G=1.0, softening=pm.spacing))
    cluster.clear_forces()
    cluster.add_force(pm)
    a_pm, a_dir = cluster.accelerations(), direct.accelerations()
    err = np.linalg.norm(a_pm - a_dir, axis=1) / np.linalg.norm(a_dir, axis=1)
    assert np.median(err) < 0.02
    assert np.abs(cluster.masses @ a_pm).max() < 1e-12
    assert cluster.potential_energy() == pytest.approx(direct.potential_energy(), rel=0.01)

    # A short run (on a coarser grid) conserves energy to the grid accuracy, and checkpoints
    # round-trip.
    cluster.clear_forces()
    cluster.add_force(ps.ParticleMesh(box_size=10.0, cells=32))
    tr = cluster.run(0.01, 50, record_every=10)
    assert np.abs(ps.relative_energy_error(tr)).max() < 0.01
    restored = ps.World.from_checkpoint(cluster.checkpoint())
    assert list(restored.forces.values()) == ["ParticleMesh"]
    assert np.array_equal(restored.accelerations(), cluster.accelerations())
    assert "ParticleMesh(box_size=10.0, cells=64" in repr(pm)

    with pytest.raises(ValueError, match="power of two"):
        ps.ParticleMesh(box_size=1.0, cells=30)
    w = ps.World()
    w.add_particle([20, 0, 0])
    w.add_force(ps.ParticleMesh(box_size=10.0, cells=16))
    with pytest.raises(ValueError, match="outside"):
        w.step(0.1)
