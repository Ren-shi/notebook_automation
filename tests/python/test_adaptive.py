import math

import numpy as np
import pytest

import physim as ps


def orbit(e):
    """Pinned unit-mass sun, light planet at apoapsis of an a = 1 orbit (period 2π)."""
    w = ps.World()
    w.add_particle([0, 0, 0], mass=1.0)
    w.pin(0)
    w.add_particle([1 + e, 0, 0], [0, math.sqrt((1 - e) / (1 + e)), 0], mass=1e-6)
    w.add_force(ps.NewtonianGravity(G=1.0))
    return w


def test_steps_adapt_to_an_eccentric_orbit():
    w = orbit(0.99)
    traj = w.run_adaptive(2 * math.pi, rtol=1e-10, atol=1e-13)
    stats = traj.metadata["adaptive"]
    assert traj.metadata["integrator"] == "dopri5"
    assert traj.n_frames == stats["accepted"] + 1
    assert stats["evaluations"] == 6 * (stats["accepted"] + stats["rejected"]) + 2
    steps = np.diff(traj.t)
    assert steps.max() / steps.min() > 1000  # tiny at periapsis, large at apoapsis
    energy = traj.kinetic + traj.potential
    assert np.abs(energy / energy[0] - 1).max() < 1e-7
    assert w.t == 2 * math.pi
    assert np.allclose(w.positions[1], [1.99, 0, 0], atol=1e-7)
    assert w.integrator == "verlet"  # the world's own integrator is untouched


def test_output_times_and_events():
    w = orbit(0.9)
    times = np.linspace(0, 4 * math.pi, 9)
    periapsis = ps.Event.radial_velocity(1, 0, direction=1)
    traj = w.run_adaptive(4 * math.pi, times=times, events=[periapsis], rtol=1e-12, atol=1e-15)
    assert np.array_equal(traj.t, times)
    assert traj.event_t == pytest.approx([math.pi, 3 * math.pi], abs=1e-8)
    assert np.linalg.norm(traj.event_pos[:, 1], axis=1) == pytest.approx([0.1, 0.1], abs=1e-9)
    # Half periods alternate between apoapsis and periapsis.
    assert np.linalg.norm(traj.pos[::2, 1], axis=1) == pytest.approx([1.9, 0.1, 1.9, 0.1, 1.9], abs=1e-8)

    stop = ps.Event.coordinate(1, "x", direction=-1, terminal=True)
    traj = w.run_adaptive(0.0, events=[stop])  # backwards in time
    assert traj.terminated_by == 0
    assert w.positions[1, 0] == pytest.approx(0.0, abs=1e-12)
    assert traj.t[-1] == w.t
    assert w.t < 4 * math.pi


def test_streams_to_a_writer(tmp_path):
    w = orbit(0.5)
    with ps.TrajectoryWriter(tmp_path / "t.npz") as writer:
        summary = w.run_adaptive(20.0, sink=writer, chunk_size=7)
    assert summary.n_frames == 0
    back = ps.load_trajectory(tmp_path / "t.npz")
    assert back.n_frames == summary.metadata["adaptive"]["accepted"] + 1
    assert back.metadata["adaptive"]["rtol"] == 1e-9
    assert back.t[-1] == 20.0


def test_errors():
    w = orbit(0.5)
    with pytest.raises(ValueError, match="outside the run"):
        w.run_adaptive(1.0, times=[0.0, 2.0])
    with pytest.raises(ValueError, match="rtol"):
        w.run_adaptive(1.0, rtol=-1.0)
    with pytest.raises(ValueError, match="max_steps") as info:
        w.run_adaptive(10.0, max_steps=4)
    assert info.value.trajectory.n_frames == 5
    assert info.value.trajectory.t[-1] == w.t

    w = ps.World()
    w.add_particle([1, 0, 0])
    w.add_rod(0, [0, 0, 0])
    with pytest.raises(ValueError, match="constraints"):
        w.run_adaptive(1.0)


def test_dopri5_fixed_step_integrator():
    w = orbit(0.5)
    w.integrator = "dopri5"
    assert w.integrator_info["order"] == 5
    traj = w.run(1e-3, 6284, record_every=6284)
    energy = traj.kinetic + traj.potential
    assert abs(energy[-1] / energy[0] - 1) < 1e-9
