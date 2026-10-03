import math

import numpy as np
import pytest

import physim as ps


def busy_world():
    """Every new built-in force at once, on a small system."""
    w = ps.World("rk4")
    w.add_particle([0, 0, 0], mass=5.0)
    for k in range(1, 6):
        a = 0.9 * k
        w.add_particle([a, 0.1 * k, 0.05], vel=[0, 1.0 / math.sqrt(a), 0.02], mass=0.1 * k)
    w.add_force(ps.NewtonianGravity(G=1.0, softening=0.05))
    w.add_force(ps.PostNewtonian(0, c=50.0))
    w.add_force(ps.J2Oblateness(0, J2=0.01, radius=0.3, axis=[0.1, 0, 1]))
    w.add_force(ps.DampedSpring(1, 2, k=0.5, rest_length=1.0, c=0.05))
    w.add_force(ps.ModulatedSpring(3, [0, 3, 0], k=0.1, depth=0.2, omega=1.3, phase=0.4))
    w.add_force(ps.ModulatedSpring(4, 5, k=0.2, depth=0.1, omega=0.7, rest_length=0.5))
    w.add_force(ps.SpringNetwork([1, 2], [4, 5], k=[0.1, 0.2], rest_length=2.0))
    w.add_force(ps.PowerLaw(k=0.01, n=1.5, center=[0.2, 0, 0]))
    w.add_force(ps.Yukawa(k=0.1, length=2.0))
    w.add_force(ps.PlummerPotential(GM=0.2, a=1.0))
    w.add_force(ps.HernquistPotential(GM=0.2, a=0.5, center=[0, 0.1, 0]))
    w.add_force(ps.HarmonicTrap([0.01, 0.02, 0.03]))
    w.add_force(ps.PeriodicForce(2, amplitude=[0, 0, 0.01], omega=2.0, phase=1.0))
    return w


def test_new_forces_survive_checkpoints_bit_for_bit(tmp_path):
    w = busy_world()
    w.run(1e-3, 200, record_every=200)
    straight = ps.World.from_checkpoint(w.checkpoint())
    assert straight.forces == w.forces
    ps.save_checkpoint(w, tmp_path / "c.npz")
    resumed = ps.load_checkpoint(tmp_path / "c.npz")
    a = straight.run(1e-3, 300, record_every=300)
    b = resumed.run(1e-3, 300, record_every=300)
    assert np.array_equal(a.pos, b.pos)
    assert np.array_equal(a.energy, b.energy)
    kinds = [f["type"] for f in a.metadata["forces"]]
    assert "SpringNetwork" in kinds and "J2Oblateness" in kinds


def test_parameters_can_be_changed():
    w = ps.World()
    w.add_particle([1, 0, 0])
    drive = w.add_force(ps.PeriodicForce(0, [1, 0, 0], omega=1.0))
    trap = w.add_force(ps.HarmonicTrap(1.0))
    w.set_force_params(drive, omega=2.5, phase=0.3)
    w.set_force_params(trap, omega=[1.0, 2.0, 3.0])
    assert w.force_params(drive)["omega"] == 2.5
    assert np.array_equal(w.force_params(trap)["omega"], [1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        w.set_force_params(trap, omega=1.0)  # a vector parameter
    pn = w.add_force(ps.PostNewtonian(0, c=10.0))
    with pytest.raises(ValueError, match="positive"):
        w.set_force_params(pn, c=-1.0)


def test_invalid_construction():
    with pytest.raises(ValueError, match="equal lengths"):
        ps.SpringNetwork([0, 1], [1], 1.0, 1.0)
    with pytest.raises(ValueError, match="itself"):
        ps.SpringNetwork([0], [0], 1.0, 1.0)
    with pytest.raises(ValueError, match="non-negative"):
        ps.DampedSpring(0, 1, 1.0, 1.0, c=-0.5)
    with pytest.raises(ValueError, match="positive"):
        ps.Yukawa(1.0, length=0.0)
    with pytest.raises(ValueError, match="axis"):
        ps.J2Oblateness(0, 1e-3, 1.0, axis=[0, 0, 0])
    with pytest.raises(ValueError, match="differ"):
        ps.ModulatedSpring(2, 2, 1.0, 0.1, 1.0)
    assert len(ps.SpringNetwork(np.arange(4), np.arange(1, 5), 1.0, np.ones(4))) == 4


def test_spring_network_equals_springs():
    def chain(network):
        w = ps.World()
        for k in range(6):
            w.add_particle([k * 1.1, 0.05 * (k % 2), 0])
        if network:
            w.add_force(ps.SpringNetwork(range(5), range(1, 6), 10.0, 1.0))
        else:
            for k in range(5):
                w.add_force(ps.Spring(k, k + 1, 10.0, 1.0))
        return w

    a = chain(True).run(1e-3, 500, record_every=50)
    b = chain(False).run(1e-3, 500, record_every=50)
    assert np.array_equal(a.pos, b.pos)


def test_driven_oscillator_amplitude_at_resonance():
    # x'' = -x - 0.2 x' + 0.1 cos t: steady amplitude F0 / (γ ω0) = 0.5.
    w = ps.World()
    w.add_particle([0, 0, 0])
    w.pin(0)
    w.add_particle([1, 0, 0])
    w.add_force(ps.DampedSpring(0, 1, k=1.0, rest_length=1.0, c=0.2))
    w.add_force(ps.PeriodicForce(1, [0.1, 0, 0], omega=1.0))
    t0 = 64 * 2 * math.pi
    times = np.linspace(t0, t0 + 2 * math.pi, 401)
    traj = w.run_adaptive(times[-1], times=times, rtol=1e-11, atol=1e-14)
    x = traj.pos[:, 1, 0] - 1.0
    assert (x.max() - x.min()) / 2 == pytest.approx(0.5, rel=1e-6)
