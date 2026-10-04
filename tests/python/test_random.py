import math

import numpy as np
import pytest

import physim as ps


def test_same_key_same_numbers():
    a = ps.random.uniform(1000, seed=3, counter=7)
    assert np.array_equal(a, ps.random.uniform(1000, seed=3, counter=7))
    assert not np.array_equal(a, ps.random.uniform(1000, seed=4, counter=7))
    assert not np.array_equal(a, ps.random.uniform(1000, seed=3, counter=8))
    # Index ranges are slices of one stream.
    assert np.array_equal(ps.random.uniform(10, 3, 7, start=990), a[990:])
    assert np.all((a > 0) & (a < 1))
    z = ps.random.normal(200_000, seed=11)
    assert abs(z.mean()) < 0.01 and abs(z.var() - 1) < 0.01


def test_langevin_draws_the_same_numbers():
    """One BAOAB step of a free particle: v₁ = c₁ v₀ + c₂ √(T/m) ξ, with ξ the generator's normals for
    (seed, counter 0, indices 0, 1, 2) — the numbers ps.random.normal gives."""
    w = ps.World()
    w.seed = 21
    w.add_particle([0, 0, 0], vel=[1.0, -0.5, 0.25], mass=2.0)
    w.use_langevin(temperature=0.7, friction=3.0)
    dt = 0.01
    w.step(dt)
    c1 = math.exp(-3.0 * dt)
    c2 = math.sqrt(1 - c1 * c1)
    xi = ps.random.normal(3, seed=21, counter=0)
    expect = np.array([1.0, -0.5, 0.25]) * c1 + c2 * math.sqrt(0.7 / 2.0) * xi
    np.testing.assert_allclose(w.velocities[0], expect, rtol=1e-12)


def test_langevin_uses_the_world_seed_and_restarts_exactly(tmp_path):
    def world(seed):
        w = ps.World()
        w.seed = seed
        for p in ps.random.positions_in_box(50, [0, 0, 0], [1, 1, 1], seed=seed):
            w.add_particle(p, mass=1.0)
        w.add_force(ps.LennardJones(epsilon=1.0, sigma=0.1, cutoff=0.3))
        w.use_langevin(temperature=1.0, friction=1.0)
        return w

    a, b, c = world(5), world(5), world(6)
    for w in (a, b, c):
        w.run(1e-4, 20)
    assert np.array_equal(a.positions, b.positions)
    assert not np.array_equal(a.positions, c.positions)
    # Checkpoint mid-run, restore, continue: bit for bit the same as running straight through.
    straight, part = world(5), world(5)
    straight.run(1e-4, 30)
    part.run(1e-4, 12)
    ps.save_checkpoint(part, tmp_path / "half.npz")
    resumed = ps.load_checkpoint(tmp_path / "half.npz")
    assert resumed.seed == 5
    resumed.run(1e-4, 18)
    assert np.array_equal(resumed.positions, straight.positions)
    assert np.array_equal(resumed.velocities, straight.velocities)


def test_old_checkpoints_get_the_default_seed():
    w = ps.World()
    w.add_particle([0, 0, 0])
    cp = w.checkpoint()
    del cp["seed"]
    assert ps.World.from_checkpoint(cp).seed == 1


def test_positions_in_box_and_sphere():
    box = ps.random.positions_in_box(20_000, [-1, 0, 2], [1, 3, 4], seed=2)
    assert box.shape == (20_000, 3)
    assert np.all(box >= [-1, 0, 2]) and np.all(box <= [1, 3, 4])
    np.testing.assert_allclose(box.mean(axis=0), [0, 1.5, 3], atol=0.03)
    ball = ps.random.positions_in_sphere(50_000, radius=2.0, centre=[1, 1, 1], seed=2)
    r = np.linalg.norm(ball - 1.0, axis=1)
    assert r.max() <= 2.0
    # Uniform in volume: P(r < x) = (x/R)³.
    for x in (0.5, 1.0, 1.5):
        assert np.mean(r < x) == pytest.approx((x / 2.0) ** 3, abs=0.01)
    np.testing.assert_allclose((ball - 1.0).mean(axis=0), 0, atol=0.02)
    # A new draw gives new points; the same draw the same points.
    assert np.array_equal(ps.random.positions_in_box(5, 0, 1, seed=2, draw=1),
                          ps.random.positions_in_box(5, 0, 1, seed=2, draw=1))
    assert not np.array_equal(ps.random.positions_in_box(5, 0, 1, seed=2, draw=1), box[:5])


def test_maxwell_boltzmann():
    w = ps.World()
    w.seed = 9
    masses = np.where(np.arange(4000) % 2 == 0, 1.0, 4.0)
    for i, m in enumerate(masses):
        w.add_particle([i, 0, 0], mass=m)
    w.pin(0)
    v = ps.random.maxwell_boltzmann(w, temperature=2.0)
    assert np.array_equal(w.velocities, v)
    assert np.all(v[0] == 0)
    np.testing.assert_allclose((masses[:, None] * v).sum(axis=0), 0, atol=1e-9)
    for m in (1.0, 4.0):
        sel = (masses == m) & (np.arange(4000) != 0)
        assert v[sel].var() == pytest.approx(2.0 / m, rel=0.05)
    # Kinetic temperature (k_B = 1): ⟨m v²⟩/3 per particle.
    assert np.mean((masses[1:, None] * v[1:] ** 2).sum(axis=1)) / 3 == pytest.approx(2.0, rel=0.03)
    again = ps.random.maxwell_boltzmann(w, temperature=2.0)
    assert np.array_equal(again, v)
    assert not np.array_equal(ps.random.maxwell_boltzmann(w, temperature=2.0, draw=1), v)
    with pytest.raises(ValueError):
        ps.random.maxwell_boltzmann(w, temperature=-1.0)


def test_setup_draws_never_reuse_run_draws():
    assert ps.random.SETUP == 2**62
    run = ps.random.normal(10, seed=1, counter=0)
    setup = ps.random.normal(10, seed=1, counter=ps.random.SETUP)
    assert not np.array_equal(run, setup)
