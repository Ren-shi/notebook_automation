import math

import numpy as np
import pytest

import physim as ps


def henon_heiles(y, py, energy, integrator="yoshida4"):
    potential = 0.5 * y * y - y**3 / 3
    px = math.sqrt(2 * (energy - potential) - py * py)
    w = ps.World(integrator)
    w.add_particle([0, y, 0], [px, py, 0])
    w.add_force(ps.HenonHeiles())
    return w


def test_lyapunov_and_megno_classify_orbits():
    chaotic = henon_heiles(0.0, 0.0, 1 / 6).lyapunov(0.01, 100_000, record_every=5000)
    regular = henon_heiles(0.0, 0.0, 1 / 12).lyapunov(0.01, 100_000, record_every=5000)
    assert chaotic["exponents"][0] > 0.05
    assert regular["exponents"][0] < 0.01
    assert regular["mean_megno"][-1] == pytest.approx(2.0, abs=0.1)
    assert chaotic["mean_megno"][-1] > 10
    assert chaotic["running"].shape == (20, 1)
    assert chaotic["t"][-1] == pytest.approx(1000.0)


def test_full_spectrum_pairs():
    out = henon_heiles(-0.1, 0.0, 1 / 6).lyapunov(0.01, 50_000, n=6, record_every=50_000)
    lam = out["exponents"]
    assert lam.shape == (6,)
    assert abs(lam.sum()) < 1e-12
    assert lam[0] + lam[5] == pytest.approx(0.0, abs=0.03 * lam[0])


def test_poincare_section_points_lie_on_the_section():
    w = henon_heiles(0.1, 0.0, 1 / 8)
    energy0 = w.total_energy()
    section = ps.Event.coordinate(0, "x", direction=1)
    t, pos, vel = ps.poincare_section(w, 0.01, 50_000, section)
    assert len(t) > 50
    assert np.abs(pos[:, 0, 0]).max() < 1e-12
    assert (vel[:, 0, 0] > 0).all()
    # Energy is conserved at every crossing, so (y, py) stays inside the energy curve.
    y, py = ps.poincare_section(
        henon_heiles(0.1, 0.0, 1 / 8), 0.01, 50_000, section,
        coordinates=lambda p, v: (p[:, 0, 1], v[:, 0, 1]),
    )
    assert np.all(0.5 * py**2 + 0.5 * y**2 - y**3 / 3 <= energy0 + 1e-6)


def test_henon_heiles_force_and_errors(tmp_path):
    f = ps.HenonHeiles(lam=0.5, center=[1, 0, 0])
    assert repr(f) == "HenonHeiles(lam=0.5, center=[1.0, 0.0, 0.0])"
    w = henon_heiles(0.0, 0.0, 0.1)
    fid = next(iter(w.forces))
    w.set_force_params(fid, lam=0.8)
    assert w.force_params(fid)["lam"] == 0.8
    restored = ps.World.from_checkpoint(w.checkpoint())
    assert restored.force_params(fid)["lam"] == 0.8
    with pytest.raises(ValueError, match="n_exponents"):
        w.lyapunov(0.01, 10, n=7)
    w.add_rod(0, [1, 0, 0])
    with pytest.raises(ValueError, match="constraints"):
        w.lyapunov(0.01, 10)
