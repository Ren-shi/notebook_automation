import numpy as np
import pytest

import physim as ps


def projectile():
    w = ps.World(integrator="verlet")
    w.add_particle([0, 1, 0], vel=[2, 4, 0])
    w.add_force(ps.UniformField([0, -9.81, 0]))
    return w


def oscillator():
    w = ps.World(integrator="rk4")
    w.add_particle([1, 0, 0])
    w.add_force(ps.AnchorSpring(0, anchor=[0, 0, 0], k=1.0))
    return w


def test_terminal_ground_hit():
    w = projectile()
    ground = ps.Event.coordinate(0, "y", direction=-1, terminal=True)
    traj = w.run(dt=0.01, steps=10_000, record_every=10, events=[ground])
    t_hit = (4 + np.sqrt(16 + 2 * 9.81)) / 9.81
    np.testing.assert_allclose(w.t, t_hit, rtol=0, atol=1e-12)
    assert traj.terminated_by == 0
    assert traj.event_t.shape == (1,) and traj.event_index.tolist() == [0]
    assert traj.event_pos.shape == traj.event_vel.shape == (1, 1, 3)
    assert abs(traj.event_pos[0, 0, 1]) < 1e-10
    assert traj.t[-1] == w.t


def test_python_event_matches_builtin():
    builtin = oscillator().run(1e-3, 10_000, 10_000, events=[ps.Event.coordinate(0, 0)])
    custom = oscillator().run(1e-3, 10_000, 10_000, events=[ps.Event(lambda t, pos, vel, mass: pos[0, 0])])
    np.testing.assert_array_equal(custom.event_t, builtin.event_t)
    np.testing.assert_allclose(builtin.event_t, np.pi / 2 + np.pi * np.arange(3), atol=1e-10)
    assert custom.terminated_by is None


def test_periapsis_and_apoapsis_events():
    w = ps.World(integrator="yoshida4")
    w.add_particle([0.5, 0, 0], vel=[0, 0.4, 0], mass=0.5)
    w.add_particle([-0.5, 0, 0], vel=[0, -0.4, 0], mass=0.5)
    w.add_force(ps.NewtonianGravity())
    a = 1 / (2 - 0.64)
    period = 2 * np.pi * a**1.5
    events = [ps.Event.radial_velocity(0, 1, direction=+1), ps.Event.radial_velocity(0, 1, direction=-1)]
    traj = w.run(dt=1e-3, steps=int(3 * period / 1e-3), record_every=1000, events=events)
    peri = traj.event_t[traj.event_index == 0]
    apo = traj.event_t[traj.event_index == 1]
    np.testing.assert_allclose(peri, (np.arange(3) + 0.5) * period, atol=1e-9)
    np.testing.assert_allclose(apo, (np.arange(1, 3)) * period, atol=1e-9)


def test_event_errors_propagate():
    def broken(t, pos, vel, mass):
        if t > 0.5:
            raise ZeroDivisionError("bad event")
        return 1.0

    w = oscillator()
    with pytest.raises(ZeroDivisionError, match="bad event") as info:
        w.run(0.1, 100, events=[ps.Event(broken)])
    assert info.value.trajectory.t[-1] == pytest.approx(0.5)
    assert w.t == pytest.approx(0.5)


def test_event_validation_and_repr():
    with pytest.raises(ValueError, match="direction"):
        ps.Event.coordinate(0, "x", direction=2)
    with pytest.raises(ValueError, match="axis"):
        ps.Event.coordinate(0, "w")
    with pytest.raises(ValueError, match="callable"):
        ps.Event(3.0)
    e = ps.Event.separation(0, 1, 2.0, direction=-1, terminal=True)
    assert (e.direction, e.terminal) == (-1, True)
    assert repr(e) == "Event(Separation(0, 1), direction=-1, terminal=True)"
    assert ps.Event(oscillator, name="custom").name == "custom"
