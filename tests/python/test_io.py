import json

import numpy as np
import pytest

import physim as ps


def pull(t, pos, vel, mass):
    return -0.01 * pos


def busy_world(integrator="yoshida4"):
    w = ps.World(integrator)
    w.add_particle([1, 0, 0], vel=[0, 1, 0], mass=1e-3)
    w.add_particle([0, 0, 0], mass=1.0)
    w.add_particle([0, 2, 0], vel=[0.5, 0, 0.1], mass=1e-3)
    w.add_particle([3, 0, 0], mass=0.5)
    w.pin(3)
    doomed = w.add_force(ps.LinearDrag(0.5))
    w.add_force(ps.NewtonianGravity(G=1.0, softening=0.01))
    w.add_force(ps.Spring(0, 2, k=0.1, rest_length=2.0))
    w.add_force(ps.AnchorSpring(2, [0, 0, 1], k=0.05))
    w.add_force(ps.UniformField([0, 0, -0.01]))
    w.add_force(ps.QuadraticDrag(1e-4))
    w.add_force(ps.CustomForce(pull, velocity_dependent=False, name="pull"))
    w.remove_force(doomed)
    w.set_force_params(1, G=1.5)
    return w


# wisdom_holman needs pure gravity without pins and boris rejects drag; their restarts are
# tested in test_integrators.py and test_em.py
@pytest.mark.parametrize(
    "integrator", [i for i in ps.INTEGRATORS if i not in ("wisdom_holman", "boris")]
)
def test_checkpoint_restart_is_bit_exact(tmp_path, integrator):
    straight = busy_world(integrator)
    interrupted = busy_world(integrator)
    straight.run(1e-3, 800, record_every=100)
    interrupted.run(1e-3, 300, record_every=100)

    ps.save_checkpoint(interrupted, tmp_path / "c.npz")
    resumed = ps.load_checkpoint(
        tmp_path / "c.npz", custom_forces={"pull": ps.CustomForce(pull, velocity_dependent=False)}
    )
    resumed.run(1e-3, 500, record_every=100)

    assert resumed.t == straight.t
    assert np.array_equal(resumed.positions, straight.positions)
    assert np.array_equal(resumed.velocities, straight.velocities)
    assert resumed.forces == straight.forces
    assert resumed.integrator == straight.integrator
    assert resumed.force_params(1) == {"G": 1.5, "softening": 0.01}
    assert resumed.pinned.tolist() == [False, False, False, True]


def test_checkpoint_needs_custom_forces_back(tmp_path):
    w = busy_world()
    ps.save_checkpoint(w, tmp_path / "c.npz")
    with pytest.raises(ValueError, match="pull"):
        ps.load_checkpoint(tmp_path / "c.npz")
    # Forces can also be given by id.
    r = ps.load_checkpoint(tmp_path / "c.npz", {6: ps.CustomForce(pull)})
    assert r.forces[6] == "pull"
    with pytest.raises(ValueError, match="npz"):
        ps.save_checkpoint(w, tmp_path / "c.bin")


def test_checkpoint_dict_is_plain_data():
    data = busy_world().checkpoint()
    arrays = {k: data.pop(k) for k in ("positions", "velocities", "masses", "charges", "radii", "pinned")}
    assert arrays["positions"].shape == (4, 3)
    restored = json.loads(json.dumps(data))  # everything else is JSON
    assert restored["integrator"] == "yoshida4"
    assert restored["next_force_id"] == 7
    assert [f["type"] for f in restored["forces"]][-1] == "external"


@pytest.mark.parametrize("ext", ["npz", "h5"])
def test_trajectory_round_trip(tmp_path, ext):
    if ext == "h5":
        pytest.importorskip("h5py")
    w = busy_world()
    traj = w.run(1e-3, 3000, record_every=7, events=[ps.Event.radial_velocity(0, 1)])
    assert len(traj.event_t) > 0
    ps.save_trajectory(traj, tmp_path / f"t.{ext}")
    back = ps.load_trajectory(tmp_path / f"t.{ext}")
    for name in ["t", "pos", "vel", "kinetic", "potential", "energy", "event_t", "event_index", "event_pos", "event_vel"]:
        assert np.array_equal(getattr(back, name), getattr(traj, name)), name
    assert back.metadata == traj.metadata
    assert back.metadata["dt"] == 1e-3 and back.metadata["integrator"] == "yoshida4"
    assert back.metadata["forces"][0] == {"id": 1, "type": "NewtonianGravity", "G": 1.5, "softening": 0.01}
    assert back.metadata["engine_version"] == ps.__version__


@pytest.mark.parametrize("ext", ["npz", "h5"])
def test_streaming_matches_in_memory_run(tmp_path, ext):
    if ext == "h5":
        pytest.importorskip("h5py")
    events = [ps.Event.coordinate(0, "x", direction=-1, terminal=False)]
    expected = busy_world().run(1e-3, 5000, events=events, energies=False)
    with ps.TrajectoryWriter(tmp_path / f"s.{ext}") as out:
        summary = busy_world().run(1e-3, 5000, events=events, energies=False, sink=out, chunk_size=333)
    assert summary.n_frames == 0 and len(summary.event_t) == 0
    assert out.n_frames == 5001
    got = ps.load_trajectory(tmp_path / f"s.{ext}")
    assert got.kinetic is None and got.energy is None
    for name in ["t", "pos", "vel", "event_t", "event_index", "event_pos"]:
        assert np.array_equal(getattr(got, name), getattr(expected, name)), name


def test_sink_receives_chunks_and_terminal_event():
    chunks = []
    w = ps.World("verlet")
    w.add_particle([0, 0, 0], vel=[0, 1, 0])
    w.add_force(ps.UniformField([0, -1, 0]))
    landed = ps.Event.coordinate(0, "y", direction=-1, terminal=True)
    summary = w.run(1e-3, 10_000, events=[landed], sink=chunks.append, chunk_size=500)
    assert summary.terminated_by == 0
    assert [len(c) for c in chunks][:-1] == [500] * (len(chunks) - 1)
    assert chunks[-1].terminated_by == 0
    assert sum(len(c.event_t) for c in chunks) == 1
    assert sum(len(c) for c in chunks) == 2002  # t = 0 ... 2.0 plus the landing frame
    assert all(c.kinetic is not None for c in chunks)


def test_sink_errors_stop_the_run():
    w = ps.World()
    w.add_particle([0, 0, 0])

    def full(chunk):
        raise OSError("disk full")

    with pytest.raises(OSError, match="disk full"):
        w.run(0.1, 100, sink=full, chunk_size=10)
    assert w.t == pytest.approx(0.9)  # stopped after the 10th frame (t = 0 .. 0.9)


def test_energies_off_skips_energy_arrays():
    traj = busy_world().run(1e-3, 10, energies=False)
    assert traj.kinetic is None and traj.potential is None and traj.energy is None
    assert traj.pos.shape == (11, 4, 3)
    assert traj.metadata["energies"] is False


def test_trajectory_constructor_validates_shapes():
    t = np.arange(3.0)
    pos = np.zeros((3, 2, 3))
    traj = ps.Trajectory(t, pos, pos)
    assert traj.n_frames == 3 and traj.n_particles == 2 and traj.kinetic is None
    with pytest.raises(ValueError, match="shape"):
        ps.Trajectory(t, np.zeros((2, 2, 3)), pos)
    with pytest.raises(ValueError, match="both"):
        ps.Trajectory(t, pos, pos, kinetic=t)
