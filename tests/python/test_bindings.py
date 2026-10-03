import numpy as np
import pytest

import physim as ps


def two_body(integrator="verlet"):
    w = ps.World(integrator=integrator)
    w.add_particle([0.5, 0, 0], vel=[0, 0.4, 0], mass=0.5)
    w.add_particle([-0.5, 0, 0], vel=[0, -0.4, 0], mass=0.5)
    return w


def py_gravity(t, pos, vel, mass):
    d = pos[None, :, :] - pos[:, None, :]  # d[i, j] = r_j - r_i
    r3 = np.linalg.norm(d, axis=-1) ** 3
    np.fill_diagonal(r3, np.inf)
    return np.einsum("ijk,j->ik", d / r3[..., None], mass)


def py_gravity_potential(t, pos, mass):
    i, j = np.triu_indices(len(mass), k=1)
    return -np.sum(mass[i] * mass[j] / np.linalg.norm(pos[i] - pos[j], axis=-1))


def test_trajectory_shapes():
    traj = two_body().run(dt=0.01, steps=100, record_every=10)
    assert len(traj) == traj.n_frames == 11
    assert traj.pos.shape == traj.vel.shape == (11, 2, 3)
    assert traj.t.shape == traj.energy.shape == (11,)
    np.testing.assert_allclose(traj.t[-1], 1.0)


def test_custom_force_matches_builtin():
    builtin = two_body()
    builtin.add_force(ps.NewtonianGravity(G=1.0))
    custom = two_body()
    custom.add_force(
        ps.CustomForce(py_gravity, potential=py_gravity_potential, velocity_dependent=False)
    )
    np.testing.assert_allclose(custom.accelerations(), builtin.accelerations(), rtol=1e-14)
    a = builtin.run(dt=1e-3, steps=2000, record_every=100)
    b = custom.run(dt=1e-3, steps=2000, record_every=100)
    np.testing.assert_allclose(b.pos, a.pos, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(b.energy, a.energy, rtol=1e-10)
    assert custom.forces == ["py_gravity"]


def test_energy_conservation_diagnostic():
    w = two_body("yoshida4")
    w.add_force(ps.NewtonianGravity())
    traj = w.run(dt=1e-3, steps=20_000, record_every=100)
    assert np.max(np.abs(ps.relative_energy_error(traj))) < 1e-8


def test_python_exceptions_keep_their_type():
    def broken(t, pos, vel, mass):
        raise ZeroDivisionError("boom")

    w = two_body()
    w.add_force(ps.CustomForce(broken))
    with pytest.raises(ZeroDivisionError, match="boom"):
        w.step(0.01)


def test_bad_custom_return_shape():
    w = two_body()
    w.add_force(ps.CustomForce(lambda t, pos, vel, mass: np.zeros(3)))
    with pytest.raises(ValueError, match=r"expected .*\(2, 3\)"):
        w.step(0.01)


def test_invalid_inputs():
    with pytest.raises(ValueError, match="unknown integrator"):
        ps.World(integrator="leapfrogg")
    w = two_body()
    with pytest.raises(ValueError):
        w.add_force(lambda *a: 0)
    with pytest.raises(ValueError, match="mass"):
        w.add_particle([0, 0, 0], mass=-1.0)
    with pytest.raises(ValueError, match=r"\(2, 3\)"):
        w.positions = np.zeros((3, 3))
    with pytest.raises(ValueError, match="3 components"):
        w.add_particle([0, 0])


def test_state_setters_and_integrator_switch():
    w = two_body()
    w.positions = [[1, 0, 0], [-1, 0, 0]]
    w.velocities = np.zeros((2, 3))
    np.testing.assert_array_equal(w.positions[:, 0], [1, -1])
    w.integrator = "rk4"
    assert w.integrator == "rk4"
    w.t = 5.0
    w.add_force(ps.UniformField([0, -1, 0]))
    w.step(0.5, n=2)
    np.testing.assert_allclose(w.t, 6.0)
    np.testing.assert_allclose(w.positions[:, 1], [-0.5, -0.5])
