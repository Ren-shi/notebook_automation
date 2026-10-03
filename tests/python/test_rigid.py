import math

import numpy as np
import pytest

import physim as ps


def test_free_top_flips_and_conserves():
    s = ps.RigidSystem()
    s.add_body(1.0, [1, 2, 3], angular_velocity=[0.01, 1.0, 0.02])
    e0 = s.total_energy()
    tr = s.run(dt=1e-3, steps=60_000, record_every=10)
    assert len(tr) == tr.n_frames == 6001 and tr.n_bodies == 1
    # Body-frame angular velocity: R^T omega.
    w_body = np.einsum("fbji,fbj->fbi", tr.rotation, tr.angular_velocity)[:, 0]
    signs = np.sign(w_body[:, 1])
    assert np.count_nonzero(np.diff(signs)) >= 3
    assert np.allclose(tr.angular_momentum, tr.angular_momentum[0], atol=1e-14)
    assert np.max(np.abs(tr.energy / e0 - 1)) < 1e-6
    assert np.allclose(np.linalg.norm(tr.orientation, axis=-1), 1.0)


def test_heavy_top_precession():
    m, g, l, i3, spin = 1.0, 1.0, 1.0, 0.5, 50.0
    s = ps.RigidSystem(order=4)
    q = ps.quaternion_from_axis_angle([1, 0, 0], 0.5)
    rot = np.array(
        [[1, 0, 0], [0, math.cos(0.5), -math.sin(0.5)], [0, math.sin(0.5), math.cos(0.5)]]
    )
    s.add_body(
        m,
        [0.5, 0.5, i3],
        pivot=[0, 0, 0],
        center_of_mass=[0, 0, l],
        orientation=q,
        angular_velocity=rot @ [0, 0, spin],
    )
    s.add_force(ps.BodyGravity([0, 0, -g]))
    assert np.allclose(s.inertia[0], [1.5, 1.5, 0.5])
    assert np.allclose(s.centers_of_mass[0], rot @ [0, 0, l])
    tr = s.run(dt=2e-3, steps=25_000, record_every=5)
    axis = tr.rotation[:, 0, :, 2]
    phi = np.unwrap(np.arctan2(axis[:, 1], axis[:, 0]))
    rate = (phi[-1] - phi[0]) / (tr.t[-1] - tr.t[0])
    assert rate == pytest.approx(m * g * l / (i3 * spin), rel=0.01)
    assert np.max(np.abs(tr.energy / tr.energy[0] - 1)) < 1e-6
    assert np.all(tr.pos == 0)


def test_springs_custom_forces_and_conservation():
    s = ps.RigidSystem(order=4)
    a = s.add_body(1.0, ps.inertia_ellipsoid(1.0, 0.5, 0.3, 0.2), angular_velocity=[0.5, 2, -1])
    b = s.add_body(
        3.0,
        ps.inertia_cylinder(3.0, 0.2, 1.0),
        position=[1.5, 0.2, 0],
        velocity=[0, 0.1, 0],
        orientation=ps.quaternion_from_axis_angle([0, 1, 0], 0.4),
    )
    s.add_force(ps.BodySpring((a, [0.4, 0.1, 0]), (b, [0, 0.1, -0.4]), k=20.0, rest_length=0.6))
    p0, l0, e0 = s.momentum(), s.angular_momentum(), s.total_energy()
    tr = s.run(2e-3, 5000, 50)
    assert np.allclose(s.momentum(), p0, atol=1e-12)
    assert np.allclose(s.angular_momentum(), l0, atol=1e-11)
    assert np.max(np.abs(tr.energy / e0 - 1)) < 1e-5

    # A Python torque: a constant twist about z spins a sphere up linearly.
    s = ps.RigidSystem()
    s.add_body(2.0, ps.inertia_sphere(2.0, 0.5))
    calls = []

    def twist(t, pos, quat):
        calls.append(quat.shape)
        return np.zeros_like(pos), np.tile([0.0, 0.0, 0.3], (len(pos), 1))

    s.add_force(twist)
    s.run(0.01, 100)
    assert calls[0] == (1, 4)
    assert s.angular_momenta[0] == pytest.approx([0, 0, 0.3])
    assert s.angular_velocities[0] == pytest.approx([0, 0, 0.3 / 0.2])
    forces, torques = s.wrenches()
    assert torques[0] == pytest.approx([0, 0, 0.3])

    # Fixed anchors, state edits and the box inertia helper.
    s = ps.RigidSystem()
    s.add_body(1.0, ps.inertia_box(1.0, 1, 2, 3))
    s.add_force(ps.BodySpring([0, 0, 2], (0, [0, 0, 1]), k=1.0))
    s.set_body_state(0, position=[0, 0, 0.5], angular_velocity=[1, 0, 0])
    assert s.angular_velocities[0] == pytest.approx([1, 0, 0])
    s.set_body_state(0, orientation=ps.quaternion_from_axis_angle([0, 0, 1], 1.0))
    assert s.angular_velocities[0] == pytest.approx([1, 0, 0])
    assert s.potential_energy() == pytest.approx(0.5 * 0.5**2)
    assert s.rotation_matrices.shape == (1, 3, 3)
    assert "RigidSystem(bodies=1" in repr(s)


def test_errors():
    s = ps.RigidSystem()
    with pytest.raises(ValueError):
        ps.RigidSystem(order=3)
    with pytest.raises(ValueError, match="triangle"):
        s.add_body(1.0, [1, 1, 3])
    with pytest.raises(ValueError, match="center_of_mass"):
        s.add_body(1.0, [1, 1, 1], center_of_mass=[0, 0, 1])
    with pytest.raises(ValueError, match="principal axis"):
        s.add_body(1.0, [1, 1, 1], pivot=[0, 0, 0], center_of_mass=[1, 1, 0])
    with pytest.raises(ValueError, match="pivot"):
        s.add_body(1.0, [1, 1, 1], pivot=[0, 0, 0], position=[1, 0, 0])
    with pytest.raises(ValueError, match="out of range"):
        s.add_force(ps.BodySpring((0, [0, 0, 0]), [0, 0, 0], k=1.0))
    with pytest.raises(ValueError):
        s.add_force(ps.UniformField([0, 0, -1]))
    s.add_body(1.0, [1, 1, 1])

    def bad(t, pos, quat):
        raise RuntimeError("boom")

    s.add_force(bad)
    with pytest.raises(RuntimeError, match="boom"):
        s.step(0.1)
