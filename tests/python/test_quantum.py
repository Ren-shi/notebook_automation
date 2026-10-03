import numpy as np
import pytest

import physim as ps


def test_construction_and_defaults():
    s = ps.Schrodinger([256], 0.1)
    assert s.method == "split_step" and s.boundary == "periodic" and s.shape == [256]
    (x,) = s.axes
    assert x[0] == pytest.approx(-12.8) and x[128] == 0.0  # [-L/2, L/2) through zero
    assert s.psi.dtype == np.complex128 and s.psi.shape == (256,)
    assert s.order == 2 and s.hbar == 1.0 and s.mass == 1.0
    assert "Schrodinger(shape=[256]" in repr(s)

    d = ps.Schrodinger([33, 17], 0.25, boundary="dirichlet")
    assert d.method == "crank_nicolson"
    assert d.axes[0][0] == -d.axes[0][-1] == -4.0  # edge nodes symmetric about zero
    # Non-power-of-two periodic grids fall back to Crank-Nicolson.
    assert ps.Schrodinger([100], 0.1).method == "crank_nicolson"
    o = ps.Schrodinger([64], 0.1, origin=0.0)
    assert o.axes[0][0] == 0.0

    with pytest.raises(ValueError, match="periodic"):
        ps.Schrodinger([64], 0.1, boundary="dirichlet", method="split_step")
    with pytest.raises(ValueError, match="method"):
        ps.Schrodinger([64], 0.1, method="euler")
    with pytest.raises(ValueError, match="order"):
        ps.Schrodinger([64], 0.1, order=3)
    with pytest.raises(ValueError, match="shape"):
        s.psi = np.zeros(10)
    with pytest.raises(ValueError, match="origin"):
        ps.Schrodinger([64], 0.1, origin=[0.0, 1.0])


def test_packet_observables_and_helpers():
    s = ps.Schrodinger([1024], 0.08)
    s.set_gaussian(center=-10, width=1.5, momentum=2.0)
    assert s.norm() == pytest.approx(1.0, abs=1e-12)
    assert s.position()[0] == pytest.approx(-10.0)
    assert s.position_variance()[0] == pytest.approx(1.5**2)
    assert s.momentum()[0] == pytest.approx(2.0)
    assert s.kinetic_energy() == pytest.approx((4 + 0.25 / 1.5**2) / 2)

    p_axes, phi = ps.quantum.momentum_density(s)
    dp = p_axes[0][1] - p_axes[0][0]
    assert np.sum(phi) * dp == pytest.approx(1.0)
    assert p_axes[0][np.argmax(phi)] == pytest.approx(2.0, abs=dp)
    # ∫ j dx = ⟨p⟩ / m for a normalised state.
    j = ps.quantum.probability_current(s)
    assert j.shape == (1, 1024)
    assert np.sum(j) * s.spacing == pytest.approx(2.0, rel=1e-9)
    # The packet is symmetric about the grid point at its centre.
    below = ps.quantum.probability(s, lambda x: x < -10)
    assert below + ps.quantum.probability(s, lambda x: x <= -10) == pytest.approx(1.0)
    assert ps.quantum.probability(s, np.ones(1024, bool)) == pytest.approx(1.0)

    t, frames = s.run(0.01, 100, record_every=30)
    assert np.allclose(t, [0, 0.3, 0.6, 0.9, 1.0])
    assert frames.shape == (5, 1024)
    assert np.allclose(frames[-1], s.psi)
    assert s.position()[0] == pytest.approx(-8.0)
    with pytest.raises(ValueError, match="record_every"):
        s.run(0.01, 10, record_every=0)

    s.psi = 2 * s.psi
    assert s.norm() == pytest.approx(4.0)
    s.normalize()
    assert s.norm() == pytest.approx(1.0)


def test_eigenstates_and_analytic_levels():
    s = ps.Schrodinger([32, 32], 0.5, potential=lambda t, x, y: 0.5 * (x**2 + y**2))
    energies, states = s.eigenstates(6)
    assert states.shape == (6, 32, 32)
    assert np.allclose(energies, ps.quantum.harmonic_levels(6, dims=2), atol=1e-6)
    assert np.sum(states[0] ** 2) * s.cell_volume == pytest.approx(1.0)

    box = ps.Schrodinger([201], 0.01, boundary="dirichlet")
    energies, _ = box.eigenstates(3)
    assert np.allclose(energies, ps.quantum.box_levels(3, 2.0), rtol=1e-3)
    with pytest.raises(ValueError, match="count"):
        box.eigenstates(0)


def test_tunnelling_matches_the_analytic_transmission():
    height, a = 1.0, 1.0
    s = ps.Schrodinger([8192], 600 / 8192)
    h = s.spacing
    s.potential = np.where(np.abs(s.axes[0]) < 0.5 * a - 0.25 * h, height, 0.0)
    a_grid = np.count_nonzero(s.potential) * h
    for p0 in [1.1, 1.6]:
        s.t = 0.0
        s.set_gaussian(center=-80, width=8, momentum=p0)
        s.step(0.05, int(130 / p0 / 0.05))
        transmitted = ps.quantum.probability(s, lambda x: x > 0.5 * a)
        expected = ps.quantum.packet_transmission(p0, 8, height, a_grid)
        assert transmitted == pytest.approx(expected, rel=3e-3), p0
    # The analytic coefficient itself: continuous at E = V₀, below 1 except at resonances.
    e = np.array([0.999999, 1.0, 1.000001])
    assert np.ptp(ps.quantum.barrier_transmission(e, 1.0, 1.0)) < 1e-5
    k_res = np.pi  # q a = π above the barrier: full transmission
    assert ps.quantum.barrier_transmission(1.0 + k_res**2 / 2, 1.0, 1.0) == pytest.approx(1.0)
    with pytest.raises(ValueError, match="positive"):
        ps.quantum.barrier_transmission(0.0, 1.0, 1.0)


def test_ehrenfest_against_a_classical_world():
    # Double well V = x⁴/4 - x²: with a small ħ, a narrow packet's ⟨x⟩ follows the classical
    # particle at first and departs once the packet spreads over the anharmonic potential.
    def v(t, x):
        return 0.25 * x**4 - x**2

    s = ps.Schrodinger([1024], 0.025, hbar=0.05, potential=v, order=4)
    s.set_gaussian(center=2.2, width=0.15, momentum=0.0)
    times, frames = s.run(0.002, 5000, record_every=50)
    density = np.abs(frames) ** 2
    x_q = density @ s.axes[0] / density.sum(axis=1)

    w = ps.World(integrator="yoshida4")
    w.add_particle([2.2, 0, 0])
    w.add_force(ps.CustomForce(lambda t, pos, vel, m: np.c_[-(pos[:, 0] ** 3) + 2 * pos[:, 0],
                                                             np.zeros((len(pos), 2))],
                               velocity_dependent=False))
    tr = w.run(0.002, 5000, record_every=50)
    x_c = tr.pos[:, 0, 0]
    assert np.allclose(tr.t, times)
    early = times <= 1.0
    assert np.abs(x_q[early] - x_c[early]).max() < 0.03
    assert np.abs(x_q - x_c).max() > 1.0


def test_time_dependent_potential_and_errors():
    calls = []

    def trap(t, x):
        calls.append(t)
        return 0.5 * (x - np.sin(0.6 * t)) ** 2

    s = ps.Schrodinger([256], 24 / 256, potential=trap)
    assert s.potential_function is trap
    s.set_gaussian(width=np.sqrt(0.5))
    s.step(0.01, 100)
    assert calls[-1] == pytest.approx(1.0)
    assert np.allclose(s.potential, 0.5 * (s.axes[0] - np.sin(0.6)) ** 2)
    s.t = 0.0  # the potential follows the clock
    assert np.allclose(s.potential, 0.5 * s.axes[0] ** 2)

    # A constant from the function broadcasts; None clears the potential.
    s.potential = lambda t, x: 3.0
    assert np.all(s.potential == 3.0)
    s.potential = None
    assert s.potential_function is None and not s.potential.any()

    class Boom(RuntimeError):
        pass

    def failing(t, x):
        if t > 0.25:
            raise Boom("no")
        return 0 * x

    s.potential = failing
    s.t = 0.0
    with pytest.raises(Boom):
        s.step(0.1, 10)
    assert s.t == pytest.approx(0.2)
    with pytest.raises(ValueError, match="shape"):
        s.potential = np.zeros(3)
    with pytest.raises(ValueError, match="non-negative"):
        s.absorber = -np.ones(256)


def test_absorbing_layer():
    s = ps.Schrodinger([512], 60 / 512)
    s.absorbing_layer(10.0, 2.0)
    assert s.absorber[256] == 0.0 and s.absorber[0] > 1.9
    s.set_gaussian(momentum=3.0)
    s.step(0.01, 2000)
    assert s.norm() < 1e-3


def test_crank_nicolson_2d():
    s = ps.Schrodinger([64, 64], 0.2, boundary="dirichlet",
                       potential=lambda t, x, y: 0.5 * (x**2 + y**2))
    s.set_gaussian(center=[1.0, 0.0], width=np.sqrt(0.5), momentum=[0.0, 1.0])
    e0 = s.energy()
    s.step(0.02, 50)
    assert s.last_iterations > 0
    assert s.norm() == pytest.approx(1.0, abs=1e-9)
    assert s.energy() == pytest.approx(e0, rel=1e-9)
    # Coherent state: ⟨x⟩ = cos t, ⟨y⟩ = sin t (up to the O(h²) grid error).
    assert np.allclose(s.position(), [np.cos(1.0), np.sin(1.0)], atol=0.015)
