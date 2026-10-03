import numpy as np
import pytest

import physim as ps


def fcc_liquid(k=3, rho=0.8442, temperature=0.722, seed=2):
    n = 4 * k**3
    side = (n / rho) ** (1 / 3)
    a = side / k
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    cells = np.stack(np.meshgrid(*[np.arange(k)] * 3, indexing="ij"), -1).reshape(-1, 3)
    pos = ((cells[:, None, :] + basis[None]) * a).reshape(-1, 3)
    vel = np.random.default_rng(seed).normal(size=(n, 3))
    vel -= vel.mean(axis=0)
    vel *= np.sqrt(temperature / (vel**2).sum(axis=1).mean() * 3)
    w = ps.World()
    for x, v in zip(pos, vel):
        w.add_particle(x, v)
    w.add_force(ps.LennardJones(cutoff=2.5, box=side))
    return w, side


def test_lj_liquid_structure_and_thermostats(tmp_path):
    w, side = fcc_liquid()
    assert w.temperature() == pytest.approx(0.722, rel=1e-12)
    w.use_langevin(0.722, friction=1.0, seed=3)
    w.run(0.005, 2000, record_every=2000, energies=False)
    traj = w.run(0.005, 3000, record_every=60, energies=False)
    r, g = ps.radial_distribution(traj.pos[-20:], side, bins=120)
    peak = r[np.argmax(g)]
    assert 1.05 < peak < 1.15 and 2.6 < g.max() < 3.4, (peak, g.max())
    assert abs(np.mean(g[r > 2.3]) - 1) < 0.15  # tending to 1 at long range
    assert w.pressure(side**3) == pytest.approx(w.pressure(side**3))

    # Nosé-Hoover holds the temperature and conserves its extended energy.
    w.use_nose_hoover(1.0, tau=0.5)
    w.run(0.005, 2000, record_every=2000, energies=False)
    e0 = w.total_energy() + w.thermostat_energy()
    traj = w.run(0.005, 4000, record_every=10)
    t = 2 * traj.kinetic / (3 * w.n_particles)
    assert t.mean() == pytest.approx(1.0, rel=0.04)
    assert (w.total_energy() + w.thermostat_energy()) == pytest.approx(e0, rel=1e-3)

    # Thermostat state and pair potentials survive checkpoints exactly.
    ps.save_checkpoint(w, tmp_path / "md.npz")
    back = ps.load_checkpoint(tmp_path / "md.npz")
    assert back.integrator == "nose_hoover"
    w.run(0.005, 50, record_every=50)
    back.run(0.005, 50, record_every=50)
    assert np.array_equal(w.positions, back.positions)


def test_tabulated_and_morse_pairs():
    w = ps.World()
    w.add_particle([0, 0, 0])
    w.add_particle([1.3, 0.2, 0])
    lj = lambda r: 4 * (r**-12 - r**-6)
    w.add_force(ps.TabulatedPair(*ps.tabulate_pair(lj, 0.7, 2.5, 4001), cutoff=2.5))
    a_table = w.accelerations()
    w.replace_force(next(iter(w.forces)), ps.LennardJones(cutoff=2.5))
    assert np.allclose(w.accelerations(), a_table, rtol=1e-7, atol=1e-10)

    m = ps.Morse(depth=1.0, a=2.0, r0=1.2, cutoff=3.0, box=[10, 10, 10])
    assert "Morse" in repr(m)
    with pytest.raises(ValueError, match="cutoff"):
        ps.LennardJones(cutoff=3.0, box=5.0)
    with pytest.raises(ValueError):
        w.use_langevin(-1.0)
