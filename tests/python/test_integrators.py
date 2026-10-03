import math

import numpy as np
import pytest

import physim as ps


def oscillator(integrator="verlet"):
    w = ps.World(integrator)
    w.add_particle([1, 0, 0])
    w.add_force(ps.AnchorSpring(0, [0, 0, 0], k=1.0))
    return w


def error(w, dt, t_end=2.0):
    w.step(dt, n=round(t_end / dt))
    return math.hypot(w.positions[0, 0] - math.cos(t_end), w.velocities[0, 0] + math.sin(t_end))


def test_new_integrators_are_listed_and_converge():
    for name in ["yoshida6", "yoshida8", "pefrl", "blanes_moan4", "gauss2", "gauss4", "gauss6"]:
        assert name in ps.INTEGRATORS
        order = oscillator(name).integrator_info["order"]
        dt = 0.1 if order >= 5 else 0.02
        measured = math.log2(error(oscillator(name), dt) / error(oscillator(name), dt / 2))
        assert measured == pytest.approx(order, abs=0.25), name
    assert ps.World("forest_ruth").integrator == "yoshida4"
    assert ps.World("implicit_midpoint").integrator == "gauss2"


def test_user_defined_composition_and_splitting():
    c = 2 ** (1 / 3)
    w = oscillator()
    w.use_composition([1 / (2 - c), -c / (2 - c), 1 / (2 - c)], order=4, name="triple_jump")
    assert w.integrator_info == {"name": "triple_jump", "order": 4, "symplectic": True}
    ref = oscillator("yoshida4")
    w.run(0.01, 300, record_every=300)
    ref.run(0.01, 300, record_every=300)
    # The same scheme up to the last bit of the weights: yoshida4 computes 2^(1/3) with cbrt,
    # and `2 ** (1 / 3)` (pow) is not correctly rounded on every platform (it differs on Windows).
    np.testing.assert_allclose(w.positions, ref.positions, rtol=0, atol=1e-14)

    w = oscillator()
    w.use_splitting([("drift", 0.5), ("kick", 1.0), ("drift", 0.5)], order=2, name="dkd")
    restored = ps.World.from_checkpoint(w.checkpoint())
    assert restored.integrator == "dkd"
    assert w.checkpoint()["integrator_scheme"]["kind"] == "splitting"
    w.run(0.01, 100, record_every=100)
    restored.run(0.01, 100, record_every=100)
    assert np.array_equal(w.positions, restored.positions)

    with pytest.raises(ValueError, match="sum to 1"):
        w.use_composition([0.5, 0.4], order=2)
    with pytest.raises(ValueError, match="kick"):
        w.use_splitting([("push", 1.0)], order=1)
    with pytest.raises(ValueError, match="built-in"):
        w.use_composition([1.0], order=2, name="verlet")
    assert w.integrator == "dkd"  # unchanged by the failed calls


def test_wisdom_holman_from_python(tmp_path):
    w = ps.World("wisdom_holman")
    w.add_particle([0, 0, 0], mass=1.0)
    w.add_particle([1, 0, 0], [0, 1, 0], mass=1e-3)
    w.add_particle([0, 2.1, 0], [-0.7, 0, 0.02], mass=5e-4)
    w.add_force(ps.NewtonianGravity(G=1.0))
    traj = w.run(0.05, 4000, record_every=40)
    assert np.abs(ps.relative_energy_error(traj)).max() < 1e-5
    ps.save_checkpoint(w, tmp_path / "wh.npz")
    assert ps.load_checkpoint(tmp_path / "wh.npz").integrator == "wisdom_holman"

    w.pin(0)
    with pytest.raises(ValueError, match="pinned"):
        w.step(0.05)
