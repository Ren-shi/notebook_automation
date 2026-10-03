import numpy as np
import pytest

import physim as ps


def test_head_on_collision_and_restitution():
    w = ps.World()
    w.add_particle([0, 0, 0], [2, 0, 0], mass=1.0, radius=0.5)
    w.add_particle([3, 0, 0], [-1, 0, 0], mass=3.0, radius=0.5)
    w.set_collisions(restitution=0.5)
    p0 = w.velocities.T @ w.masses
    w.run(0.1, 30, record_every=30)
    assert w.collision_count == 1
    assert np.allclose(w.velocities.T @ w.masses, p0, atol=1e-12)
    # Relative speed reversed and scaled by e.
    assert w.velocities[1, 0] - w.velocities[0, 0] == pytest.approx(0.5 * 3.0, abs=1e-12)
    assert np.array_equal(w.radii, [0.5, 0.5])


def test_box_gas_conserves_energy_and_stays_inside(tmp_path):
    rng = np.random.default_rng(3)
    w = ps.World()
    for k in range(64):
        x = 1.0 + 2.0 * np.array([k % 4, (k // 4) % 4, k // 16])
        w.add_particle(x, rng.normal(size=3), radius=0.3)
    w.set_collisions(walls=ps.box_walls([0, 0, 0], [8, 8, 8]))
    e0 = w.kinetic_energy()
    traj = w.run(0.05, 400, record_every=40)
    assert w.kinetic_energy() == pytest.approx(e0, rel=1e-12)
    assert w.collision_count > 50
    assert traj.pos.min() >= 0.3 - 1e-9 and traj.pos.max() <= 7.7 + 1e-9
    assert w.collisions["restitution"] == 1.0 and len(w.collisions["walls"]) == 6

    ps.save_checkpoint(w, tmp_path / "gas.npz")
    back = ps.load_checkpoint(tmp_path / "gas.npz")
    assert back.collisions == w.collisions
    w.run(0.05, 100, record_every=100)
    back.run(0.05, 100, record_every=100)
    assert np.array_equal(w.positions, back.positions)


def test_soft_contact():
    w = ps.World()
    w.add_particle([-1, 0, 0], [1, 0, 0], radius=0.5)
    w.add_particle([1, 0, 0], [-1, 0, 0], radius=0.5)
    w.add_force(ps.SoftContact(k=1e4, law="hertz"))
    assert repr(ps.SoftContact(10.0)) == "SoftContact(k=10.0, damping=0.0, law='linear')"
    traj = w.run(1e-4, 20_000, record_every=100)
    assert w.velocities[1, 0] == pytest.approx(1.0, abs=1e-4)  # elastic rebound
    assert np.abs(ps.relative_energy_error(traj)).max() < 1e-4
    with pytest.raises(ValueError, match="law"):
        ps.SoftContact(1.0, law="spring")


def test_errors():
    w = ps.World()
    w.add_particle([0, 0, 0], radius=0.1)
    with pytest.raises(ValueError, match="restitution"):
        w.set_collisions(restitution=2.0)
    with pytest.raises(ValueError, match="radius"):
        w.add_particle([0, 0, 0], radius=-1.0)
    w.set_collisions()
    with pytest.raises(ValueError, match="collisions"):
        w.run_adaptive(1.0)
    w.clear_collisions()
    assert w.collisions is None
