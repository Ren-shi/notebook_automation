import numpy as np
import pytest

import physim as ps


def test_force_ids_and_parameter_sweep():
    w = ps.World()
    w.add_particle([1, 0, 0])
    w.add_particle([0, 0, 0])
    g = w.add_force(ps.NewtonianGravity(G=1.0))
    f = w.add_force(ps.UniformField([0, -1, 0]))
    assert w.forces == {g: "NewtonianGravity", f: "UniformField"}

    accs = []
    for G in [1.0, 2.0, 3.0]:
        w.set_force_params(g, G=G)
        accs.append(w.accelerations()[0, 0])
    np.testing.assert_allclose(accs, [-1.0, -2.0, -3.0])

    assert w.force_params(g) == {"G": 3.0, "softening": 0.0}
    np.testing.assert_array_equal(w.force_params(f)["g"], [0, -1, 0])
    w.set_force_params(f, g=[0, 0, -2])
    np.testing.assert_array_equal(w.force_params(f)["g"], [0, 0, -2])

    with pytest.raises(ValueError, match="non-negative"):
        w.set_force_params(g, G=9.0, softening=-1.0)
    assert w.force_params(g)["G"] == 3.0  # rolled back
    with pytest.raises(ValueError, match="no parameter"):
        w.set_force_params(g, mass=1.0)

    w.replace_force(f, ps.UniformField([0, 0, 0]))
    w.remove_force(g)
    assert w.forces == {f: "UniformField"}
    np.testing.assert_array_equal(w.accelerations(), 0.0)
    with pytest.raises(ValueError, match="no force with id"):
        w.remove_force(g)


def test_remove_particle_and_masses():
    w = ps.World()
    for x in range(3):
        w.add_particle([x, 0, 0], mass=x + 1.0)
    s = w.add_force(ps.Spring(1, 2, k=1.0, rest_length=0.5))
    with pytest.raises(ValueError, match="still used by Spring"):
        w.remove_particle(2)
    w.remove_particle(0)
    assert w.n_particles == 2
    np.testing.assert_array_equal(w.masses, [2.0, 3.0])
    assert w.forces == {s: "Spring(0, 1)"}

    w.masses = [4, 5]
    np.testing.assert_array_equal(w.masses, [4.0, 5.0])
    with pytest.raises(ValueError, match="expected 2 masses"):
        w.masses = [1.0]
    with pytest.raises(ValueError, match="non-negative"):
        w.masses = [1.0, -1.0]


def test_tracers_and_pinning():
    w = ps.World(integrator="yoshida4")
    sun = w.add_particle([0, 0, 0], mass=1.0)
    w.add_particle([1, 0, 0], vel=[0, 1, 0], mass=0.0)
    w.add_force(ps.NewtonianGravity())
    w.run(dt=1e-3, steps=5000)
    np.testing.assert_array_equal(w.velocities[sun], 0.0)
    np.testing.assert_allclose(np.linalg.norm(w.positions[1]), 1.0, rtol=1e-9)

    w.pin(1)
    assert list(w.pinned) == [False, True]
    np.testing.assert_array_equal(w.velocities[1], 0.0)
    with pytest.raises(ValueError, match="pinned"):
        w.velocities = [[0, 0, 0], [1, 0, 0]]
    before = w.positions[1].copy()
    w.step(1e-3, n=100)
    np.testing.assert_array_equal(w.positions[1], before)


def test_failed_run_returns_partial_trajectory():
    def fails_late(t, pos, vel, mass):
        if t > 0.55:
            raise RuntimeError("diverged")
        return np.zeros_like(pos)

    w = ps.World()
    w.add_particle([0, 0, 0], vel=[1, 0, 0])
    w.add_force(ps.CustomForce(fails_late))
    with pytest.raises(RuntimeError, match="diverged") as info:
        w.run(dt=0.1, steps=100, record_every=2)
    traj = info.value.trajectory
    np.testing.assert_allclose(traj.t, [0.0, 0.2, 0.4, 0.5])
    np.testing.assert_allclose(w.t, 0.5)
    np.testing.assert_allclose(w.positions[0, 0], 0.5)


def test_integrator_info():
    w = ps.World(integrator="yoshida4")
    assert w.integrator_info == {"name": "yoshida4", "order": 4, "symplectic": True}
    w.integrator = "rk4"
    assert w.integrator_info["symplectic"] is False
