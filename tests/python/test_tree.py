import numpy as np
import pytest

import physim as ps


def cluster(force, n=800, seed=4):
    rng = np.random.default_rng(seed)
    w = ps.World()
    for x in rng.normal(size=(n, 3)):
        w.add_particle(x, mass=1.0 / n)
    w.add_force(force)
    return w


def test_tree_matches_direct_and_survives_checkpoints(tmp_path):
    direct = cluster(ps.NewtonianGravity(softening=0.01)).accelerations()
    for theta, quad, tol in [(0.5, False, 3e-3), (0.5, True, 1e-3), (0.3, True, 2e-4)]:
        a = cluster(ps.TreeGravity(softening=0.01, theta=theta, quadrupole=quad)).accelerations()
        err = np.sqrt(np.mean(np.sum((a - direct) ** 2, axis=1) / np.sum(direct**2, axis=1)))
        assert err < tol, (theta, quad, err)

    w = cluster(ps.TreeGravity(theta=0.6, quadrupole=True))
    assert repr(ps.TreeGravity()) == "TreeGravity(G=1.0, softening=0.0, theta=0.5, quadrupole=False)"
    fid = next(iter(w.forces))
    w.set_force_params(fid, theta=0.4)
    ps.save_checkpoint(w, tmp_path / "c.npz")
    back = ps.load_checkpoint(tmp_path / "c.npz")
    assert back.force_params(fid)["theta"] == 0.4
    w.run(1e-3, 10, record_every=10)
    back.run(1e-3, 10, record_every=10)
    assert np.array_equal(w.positions, back.positions)
    with pytest.raises(ValueError, match="theta"):
        ps.TreeGravity(theta=1.5)
