import json
import math

import numpy as np
import pytest

import physim as ps


def gyrating(integrator="boris"):
    w = ps.World(integrator)
    w.add_particle([0, 1, 0], [1, 0, 0], mass=1.0, charge=1.0)
    w.add_force(ps.MagneticField([0, 0, 1]))
    return w


def test_charges_api():
    w = ps.World()
    w.add_particle([0, 0, 0], charge=2.0)
    w.add_particle([1, 0, 0])
    assert np.array_equal(w.charges, [2.0, 0.0])
    w.charges = [1.0, -1.0]
    assert np.array_equal(w.charges, [1.0, -1.0])
    with pytest.raises(ValueError):
        w.charges = [1.0]
    w.remove_particle(0)
    assert np.array_equal(w.charges, [-1.0])


def test_boris_gyration_keeps_speed():
    w = gyrating()
    traj = w.run(0.01, 100_000, record_every=100, energies=False)
    speed = np.linalg.norm(traj.vel[:, 0], axis=1)
    assert np.abs(speed - 1).max() < 1e-13
    assert np.abs(np.linalg.norm(traj.pos[:, 0], axis=1) - 1).max() < 1e-10
    assert traj.metadata["charges"] == [1.0]


def test_e_cross_b_drift_and_field_functions():
    # Built-in uniform fields and the same fields as Python functions agree exactly.
    def build(use_functions):
        w = ps.World("boris")
        w.add_particle([0, 0, 0], [0.3, 0, 0], mass=2.0, charge=4.0)
        if use_functions:
            w.add_force(ps.FieldForce(
                E=lambda t, pos: np.tile([0, 0.1, 0], (len(pos), 1)),
                B=lambda t, pos: np.tile([0, 0, 2.0], (len(pos), 1)),
            ))
        else:
            w.add_force(ps.ElectricField([0, 0.1, 0]))
            w.add_force(ps.MagneticField([0, 0, 2.0]))
        return w

    a, b = build(False), build(True)
    a.run(0.01, 5000, record_every=5000)
    b.run(0.01, 5000, record_every=5000)
    assert np.array_equal(a.positions, b.positions)
    # Mean velocity over a long time approaches E × B / B² = (0.05, 0, 0).
    assert a.positions[0, 0] / a.t == pytest.approx(0.05, abs=2e-3)
    assert abs(a.positions[0, 1] / a.t) < 2e-3


def test_coulomb_and_checkpoints(tmp_path):
    w = ps.World("yoshida4")
    w.add_particle([-1, 0, 0], mass=1.0, charge=1.0)
    w.add_particle([1, 0, 0], mass=1.0, charge=1.0)
    w.add_force(ps.Coulomb(k=2.0))
    assert w.potential_energy() == pytest.approx(1.0)  # k q q / r
    traj = w.run(1e-3, 2000, record_every=100)
    assert np.abs(ps.relative_energy_error(traj)).max() < 1e-10
    assert w.positions[1, 0] > 1.0  # like charges repel
    assert np.linalg.norm(w.velocities.sum(axis=0)) < 1e-14  # momentum conserved

    ps.save_checkpoint(w, tmp_path / "c.npz")
    back = ps.load_checkpoint(tmp_path / "c.npz")
    assert np.array_equal(back.charges, [1.0, 1.0])
    assert "Coulomb" in back.forces.values()

    # Checkpoints written before charges existed still load, with neutral particles.
    cp = w.checkpoint()
    del cp["charges"]
    assert np.array_equal(ps.World.from_checkpoint(cp).charges, [0.0, 0.0])


def test_errors():
    w = gyrating()
    w.add_force(ps.LinearDrag(0.1))
    with pytest.raises(ValueError, match="LinearDrag"):
        w.step(0.01)
    with pytest.raises(ValueError, match="E, B"):
        ps.FieldForce()
    with pytest.raises(ValueError, match="charge must be finite"):
        ps.World().add_particle([0, 0, 0], charge=math.inf)
