import math

import numpy as np
import pytest

import physim as ps

G = 9.81


def ellipk(k):
    """Complete elliptic integral of the first kind, by the arithmetic-geometric mean."""
    a, b = 1.0, math.sqrt(1.0 - k * k)
    for _ in range(40):
        a, b = 0.5 * (a + b), math.sqrt(a * b)
    return math.pi / (2.0 * a)


def pendulum(theta0, length=1.0, integrator="yoshida4"):
    w = ps.World(integrator)
    w.add_particle([length * math.sin(theta0), -length * math.cos(theta0), 0], mass=1.0)
    w.add_force(ps.UniformField([0, -G, 0]))
    rod = w.add_rod(0, [0, 0, 0])
    return w, rod


def test_large_amplitude_period_matches_elliptic_integral():
    theta0, length = 2.5, 0.6
    w, _ = pendulum(theta0, length)
    exact = 4 * math.sqrt(length / G) * ellipk(math.sin(theta0 / 2))
    crossing = ps.Event.coordinate(0, "x", direction=-1)
    traj = w.run(1e-3, int(6 * exact / 1e-3), record_every=1000, events=[crossing])
    period = np.diff(traj.event_t).mean()
    assert period == pytest.approx(exact, rel=1e-8)


def test_rods_hold_and_report_tension():
    w, rod = pendulum(0.0)
    traj = w.run(1e-3, 200, record_every=50)
    assert traj.tension.shape == (traj.n_frames, 1)
    assert traj.tension[0, 0] == 0.0  # no step taken yet
    assert traj.tension[-1, 0] == pytest.approx(G, rel=1e-9)  # m g, hanging at rest
    assert w.constraint_tensions() == pytest.approx([G], rel=1e-9)
    assert w.constraints == {rod: "Rod(0, [0, 0, 0])"}

    # A chain: rods between particles, lengths taken from the current positions.
    w = ps.World("verlet")
    for k in range(5):
        w.add_particle([k + 1.0, 0, 0])
    w.add_force(ps.UniformField([0, -G, 0]))
    w.add_rod(0, [0, 0, 0])
    for k in range(1, 5):
        w.add_rod(k, k - 1)
    traj = w.run(1e-3, 3000, record_every=100)
    seg = np.linalg.norm(np.diff(traj.pos, axis=1, prepend=0.0), axis=2)
    assert np.abs(seg - 1.0).max() < 1e-9
    assert traj.tension.shape == (traj.n_frames, 5)


def test_rod_errors():
    w = ps.World("rk4")
    w.add_particle([1, 0, 0])
    w.add_particle([0, 0, 0], mass=0.0)
    with pytest.raises(ValueError, match="massless"):
        w.add_rod(1, [0, 0, 0])
    with pytest.raises(ValueError, match="place them first"):
        w.add_rod(0, [0, 0, 0], length=2.0)
    rod = w.add_rod(0, [0, 0, 0], length=1.0)
    with pytest.raises(ValueError, match="verlet or yoshida4"):
        w.step(1e-3)
    w.integrator = "verlet"
    w.step(1e-3)
    w.constraint_tolerance = 1e-6
    assert w.constraint_tolerance == 1e-6
    with pytest.raises(ValueError):
        w.constraint_tolerance = 0.0
    w.remove_constraint(rod)
    assert w.constraints == {}


def test_constraints_survive_checkpoint_and_trajectory_files(tmp_path):
    w = ps.World("yoshida4")
    w.add_particle([1, 0, 0])
    w.add_particle([1, -1, 0], mass=0.5)
    w.add_force(ps.UniformField([0, -G, 0]))
    w.add_rod(0, [0, 0, 0])
    w.add_rod(1, 0)
    w.constraint_tolerance = 1e-11
    straight = ps.World.from_checkpoint(w.checkpoint())
    straight.run(1e-3, 900, record_every=900)

    w.run(1e-3, 400, record_every=400)
    ps.save_checkpoint(w, tmp_path / "c.npz")
    resumed = ps.load_checkpoint(tmp_path / "c.npz")
    assert resumed.constraints == w.constraints
    assert resumed.constraint_tolerance == 1e-11
    resumed.run(1e-3, 500, record_every=500)
    assert np.array_equal(resumed.positions, straight.positions)
    assert np.array_equal(resumed.velocities, straight.velocities)

    for ext in ["npz", "h5"]:
        if ext == "h5":
            pytest.importorskip("h5py")
        traj = resumed.run(1e-3, 100, record_every=10)
        assert traj.metadata["constraints"][1] == {"id": 1, "i": 1, "j": 0, "length": 1.0}
        ps.save_trajectory(traj, tmp_path / f"t.{ext}")
        back = ps.load_trajectory(tmp_path / f"t.{ext}")
        assert np.array_equal(back.tension, traj.tension)
