"""Recoil transport for lifetimes (``physim.nuclear.transport``, wave 1C): straight slowing exactly as item 73 did it,
vacuum gaps and stoppers against hand-calculated numbers, time dilation, the ZBL scattering Monte Carlo
(reproducible, its energy bookkeeping, its spread against TRIM) and the DSAM kernel."""

import csv
import inspect
import math
import time
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import Experiment, data
from physim.nuclear import transport as tr
from physim.nuclear.gamma_events import _clock as ge_clock
from physim.nuclear.gamma_events import _slowing as ge_slowing
from physim.nuclear.planner import Planner
from physim.nuclear.rates import _after, _exit_paths, stack, stopping
from physim.nuclear.stopping import Stopping

LIFETIMES = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "lifetimes"
C_UM_PS = 299.792458

RU, SM = "98Ru", "143Sm"
M_RU, M_SM = data.nuclide(RU).nuclear_mass_mev, data.nuclide(SM).nuclear_mass_mev
CU, AU, MG = data.material("Cu"), data.material("Au"), data.material("Mg")
Z = np.array([0.0, 0.0, 1.0])


def energy_at(mass, beta):
    return mass * (1 / math.sqrt(1 - beta * beta) - 1)


E_RU = energy_at(M_RU, 0.045)  # ⁹⁸Ru after the Mg target of Radeck et al. (2012)
E_SM = energy_at(M_SM, 0.016)  # ¹⁴³Sm in the Au backing of Rajbanshi et al. (2018)


def plunger(gap_um=100.0):
    return [tr.Segment(MG, 0.7, name="target"), tr.Segment(None, gap_um=gap_um, name="gap"),
            tr.Segment(CU, 15.7, name="stopper")]


def stop_time_by_quadrature(st: Stopping, mass, energy, density, proper=False):
    """∫ dE / (S ρ v) (ps) from 0 to each energy up to ``energy``: S from Stopping.stopping_power on a fine grid
    from the lowest energy of Stopping's table, and below it the residual range crossed at that speed (the clock's
    convention: the last few keV, at β below 10⁻³, take a time that depends on how the end is cut)."""
    e0 = Stopping.GRID[0] * st.mass_u
    e = np.geomspace(e0, energy, 20001)
    s = np.asarray(st.stopping_power(e))  # MeV/(mg/cm²)
    beta = np.sqrt(e * (e + 2 * mass)) / (e + mass)
    speed = beta * (1 + e / mass) if proper else beta
    f = 1.0 / (s * 1e3 * density * 2.99792458e-2 * speed)
    t = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(e))])
    t += st.range(e0) / (1e3 * density * 2.99792458e-2 * speed[0])
    return e, t


# -- the frozen interface ------------------------------------------------------------------------------------------


def test_signatures_match_the_frozen_interface():
    def names(f):
        return list(inspect.signature(f).parameters)

    assert names(tr.slow) == ["segments", "species", "mass_mev", "energy", "direction", "start_depth_mg_cm2",
                              "t_decay_ps", "scatter", "rng"]
    assert names(tr.time_to_stop) == ["segments", "species", "mass_mev", "energy", "start_depth_mg_cm2"]
    assert names(tr.kernel) == ["segments", "species", "mass_mev", "energy", "direction", "start_depth_mg_cm2",
                                "times_ps", "scatter", "rng", "n"]
    for f, kw in ((tr.slow, ("scatter", "rng")), (tr.kernel, ("scatter", "rng", "n"))):
        params = inspect.signature(f).parameters
        assert all(params[k].kind is inspect.Parameter.KEYWORD_ONLY for k in kw)
    assert inspect.signature(tr.slow).parameters["scatter"].default is False
    assert inspect.signature(tr.kernel).parameters["scatter"].default is True
    assert inspect.signature(tr.kernel).parameters["n"].default == 20000
    assert [f.name for f in tr.Segment.__dataclass_fields__.values()] == ["material", "thickness_mg_cm2", "gap_um",
                                                                          "name"]
    assert list(tr.Slowed.__dataclass_fields__)[:5] == ["energy", "direction", "segment", "stopped", "beta"]
    with pytest.raises(Exception):
        tr.Segment(CU, 1.0).thickness_mg_cm2 = 2.0  # frozen


def test_a_stack_needs_matter_and_densities():
    with pytest.raises(ValueError, match="matter"):
        tr.slow([tr.Segment(None, gap_um=10.0)], RU, M_RU, E_RU, Z, 0.0, 1.0)
    with pytest.raises(ValueError, match="density"):
        tr.slow([tr.Segment(data.Material(CU.name, CU.atoms, None, CU.I_eV), 1.0)], RU, M_RU, E_RU, Z, 0.0, 1.0)


# -- straight slowing: item 73 exactly ----------------------------------------------------------------------------


_OLD_CLOCKS: dict = {}


def _old_clock(species, material, density, mass_mev):
    """``gamma_events._clock`` before the move, verbatim (seconds)."""
    key = (species, material.name, density, mass_mev)
    if key in _OLD_CLOCKS:
        return _OLD_CLOCKS[key]
    st = stopping(species, material)
    e, r = np.asarray(st._e, dtype=float), np.asarray(st._range, dtype=float)
    beta = np.sqrt(np.maximum(e * (e + 2 * mass_mev), 0.0)) / (e + mass_mev)
    mid = 0.5 * (beta[1:] + beta[:-1])
    dt = np.diff(r) / (1000.0 * density * 2.99792458e10 * np.maximum(mid, 1e-12))
    t = np.r_[r[0] / (1000.0 * density * 2.99792458e10 * max(beta[0], 1e-12)), 0.0]
    t = np.cumsum(np.r_[t[0], dt])
    _OLD_CLOCKS[key] = (np.r_[0.0, e], np.r_[0.0, t])
    return _OLD_CLOCKS[key]


def _old_slowing(experiment, layers, depth, dirs, species, energy, t_decay, mass_mev):
    """``gamma_events._slowing`` before the move, verbatim."""
    fwd, steps = _exit_paths(experiment, layers, 0, depth, dirs)
    e = np.asarray(energy, dtype=float).copy()
    elapsed = np.zeros(len(e))
    done = np.zeros(len(e), dtype=bool)
    stopped = np.zeros(len(e), dtype=bool)
    out = np.full(len(e), np.nan)
    for j, path, forward in steps:
        m = (fwd == forward) & ~done & (np.asarray(path) > 0)
        if not m.any():
            continue
        mat = layers[j].material
        grid_e, grid_t = _old_clock(species, mat, float(mat.density_g_cm3), float(mass_mev))
        idx = np.flatnonzero(m)
        e_in = e[idx]
        e_out = _after(stopping(species, mat), e_in, np.asarray(path)[idx])
        t_in, t_out = np.interp(e_in, grid_e, grid_t), np.interp(e_out, grid_e, grid_t)
        here = t_decay[idx] < elapsed[idx] + (t_in - t_out)
        out[idx[here]] = np.interp(t_in[here] - (t_decay[idx[here]] - elapsed[idx[here]]), grid_t, grid_e)
        done[idx[here]] = True
        rest = (e_out <= 0) & ~here
        out[idx[rest]] = 0.0
        stopped[idx[rest]] = True
        done[idx[rest]] = True
        go = ~here & ~rest
        elapsed[idx[go]] += (t_in - t_out)[go]
        e[idx[go]] = e_out[go]
    out[~done] = e[~done]
    return out, stopped


def ni58(tilt=None) -> Experiment:
    d = Planner.example("coulex_ni58").draft
    d["target"]["thickness"] = "1 mg/cm2"
    d["target"]["backing"] = {"material": "Au", "thickness": "20 mg/cm2"}
    if tilt is not None:
        d["target"]["tilt"] = tilt
    return Experiment.from_dict(d)


@pytest.mark.parametrize("tilt", [None, "30 deg"])
def test_straight_slowing_is_item_73s(tilt):
    exp = ni58(tilt)
    layers = stack(exp)
    rng = np.random.default_rng(3)
    n = 50_000
    th = np.arccos(rng.uniform(-0.4, 1, n))  # forward and backward, some nearly in the plane of the target
    ph = rng.uniform(0, 2 * np.pi, n)
    dirs = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    depth = rng.uniform(0, layers[0].thickness, n)
    e = rng.uniform(5, 40, n)
    t = rng.exponential(1e-12, n)  # s
    m = data.nuclide("58Ni").nuclear_mass_mev + 1.4542
    want, want_stopped = _old_slowing(exp, layers, depth, dirs, "58Ni", e, t, m)
    got, got_stopped = ge_slowing(exp, layers, depth, dirs, "58Ni", e, t, m)
    assert np.array_equal(got_stopped, want_stopped)
    assert 0.1 < want_stopped.mean() < 0.9
    # The same to rounding: the clock is now kept in ps rather than s.
    np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-12)
    # Through transport.slow directly (no tilt), with the segments of the stack.
    if tilt is None:
        segs = [tr.Segment(lay.material, lay.thickness, name=lay.name) for lay in layers]
        s = tr.slow(segs, "58Ni", m, e, dirs, depth, t * 1e12)
        np.testing.assert_allclose(s.energy, want, rtol=1e-9, atol=1e-12)
        assert np.array_equal(s.direction, dirs / np.linalg.norm(dirs, axis=1)[:, None])
        np.testing.assert_allclose(s.beta, np.sqrt(1 - 1 / (1 + s.energy / m) ** 2), rtol=1e-12)
        # Segment of the decay: target, backing, or out (the back, or the front for backward recoils).
        assert set(np.unique(s.segment)) <= {-1, 0, 1, 2}
        assert np.all(s.energy[s.segment == -1] > 0) and np.all(dirs[s.segment == -1, 2] < 0)


def test_gamma_events_clock_is_the_transport_clock_in_seconds():
    m = data.nuclide("58Ni").nuclear_mass_mev
    e, t = ge_clock("58Ni", AU, float(AU.density_g_cm3), m)
    e0, t0 = _old_clock("58Ni", AU, float(AU.density_g_cm3), m)
    assert np.array_equal(e, e0)
    np.testing.assert_allclose(t, t0, rtol=1e-12)


# -- vacuum gaps and stoppers by hand ------------------------------------------------------------------------------


def test_a_vacuum_gap_is_crossed_at_constant_speed():
    v = 0.045 * C_UM_PS  # 13.4907 µm/ps
    assert v == pytest.approx(13.49066, rel=1e-6)
    stack_ = plunger(100.0)
    # From the back face of the target (no matter left to cross), along the beam: 100 µm take 7.41253 ps.
    s = tr.slow(stack_, RU, M_RU, E_RU, Z, 0.7, [1.0, 5.0, 7.40, 7.42])
    assert list(s.segment) == [1, 1, 1, 2]
    np.testing.assert_allclose(s.energy[:3], E_RU, rtol=1e-12)
    np.testing.assert_allclose(s.beta[:3], 0.045, rtol=1e-9)
    np.testing.assert_allclose(s.position_um[:3, 2], v * np.array([1.0, 5.0, 7.40]), rtol=1e-9)
    assert s.energy[3] < E_RU and not s.stopped.any()
    # At 30° the path is 100/cos 30° = 115.470 µm: 8.55934 ps.
    d30 = np.array([math.sin(math.radians(30)), 0.0, math.cos(math.radians(30))])
    s = tr.slow(stack_, RU, M_RU, E_RU, d30, 0.7, [8.55, 8.57])
    assert list(s.segment) == [1, 2]
    np.testing.assert_allclose(s.position_um[0], d30 * v * 8.55, rtol=1e-9)
    # Backward through a gap upstream of the target: the recoil leaves by the front face of the stack (−1).
    up = [tr.Segment(None, gap_um=50.0), tr.Segment(MG, 0.7), tr.Segment(CU, 15.7)]
    s = tr.slow(up, RU, M_RU, E_RU, -Z, 0.0, [1.0, 10.0])
    assert list(s.segment) == [0, -1]
    np.testing.assert_allclose(s.energy, E_RU, rtol=1e-12)


def test_the_stopper_stops_the_recoil_on_time():
    st = stopping(RU, CU)
    rho = float(CU.density_g_cm3)
    e, t = stop_time_by_quadrature(st, M_RU, E_RU, rho)
    t_cu = t[-1]
    assert 0.5 < t_cu < 3.0  # about a picosecond in copper
    gap = 100.0 / (0.045 * C_UM_PS)
    # Time to stop: the gap, then the copper (the same integral, done independently).
    got = tr.time_to_stop(plunger(100.0), RU, M_RU, E_RU, 0.7)
    assert got[0] == pytest.approx(gap + t_cu, rel=2e-3)
    # From inside the target the recoil first slows in the Mg.
    inside = tr.time_to_stop(plunger(100.0), RU, M_RU, E_RU, 0.0)
    e_out = float(_after(stopping(RU, MG), E_RU, 0.7))
    _, t_from = stop_time_by_quadrature(st, M_RU, e_out, rho)
    assert inside[0] > got[0] - t_cu + t_from[-1] - 1e-9  # slower across the gap, and the Mg takes time
    # Too thin a stopper: the recoil leaves, and never stops.
    thin = [tr.Segment(MG, 0.7), tr.Segment(None, gap_um=10.0), tr.Segment(CU, 2.0)]
    assert np.isinf(tr.time_to_stop(thin, RU, M_RU, E_RU, 0.7)[0])
    # The energy at a time in the stopper: the quadrature inverted.
    for tc in (0.05, 0.2, 0.5):
        s = tr.slow(plunger(100.0), RU, M_RU, E_RU, Z, 0.7, gap + tc)
        want = np.interp(t_cu - tc, t, e)
        assert s.segment[0] == 2 and s.energy[0] == pytest.approx(want, rel=5e-3)
    s = tr.slow(plunger(100.0), RU, M_RU, E_RU, Z, 0.7, gap + t_cu * 1.01)
    assert s.stopped[0] and s.energy[0] == 0 and s.segment[0] == 2
    assert s.position_um[0, 2] == pytest.approx(100.0 + st.range(E_RU) * 10 / rho, rel=1e-9)


class _VelocityStopping:
    """A stopping table S = k √E (stopping proportional to the speed), on Stopping's own grid."""

    def __init__(self, k, mass_u):
        self._e = Stopping.GRID * mass_u
        self.k = k
        self._range = 2 * np.sqrt(self._e) / k
        self._se = k * np.sqrt(self._e) * 1e3
        self._sn = np.zeros_like(self._e)

    def range(self, e):
        e = np.asarray(e, dtype=float)
        return np.where(e > 0, 2 * np.sqrt(np.maximum(e, 0)) / self.k, 0.0)


def test_stopping_proportional_to_speed_gives_the_analytic_attenuation(monkeypatch):
    # v(t) = v₀ e^(−t/α): with S = k √E, dE/dt = −k ρ c √(2/m) E, so α = 2 / (k ρ c √(2/m)) (to O(β²)).
    alpha = 0.2  # ps, as in the analytic fixture
    rho = float(CU.density_g_cm3)
    k = 2 / (alpha * 1e3 * rho * 2.99792458e-2 * math.sqrt(2 / M_RU))
    fake = _VelocityStopping(k, M_RU / (data.U_KEV * 1e-3))
    monkeypatch.setattr(tr, "stopping", lambda ion, material: fake)
    monkeypatch.setattr(tr, "_CLOCKS", {})
    segs = [tr.Segment(CU, 100.0)]
    times = np.array([0.0, 0.05, 0.1, 0.2, 0.4, 0.8])
    s = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, times)
    beta0 = s.beta[0]
    np.testing.assert_allclose(s.beta / beta0, np.exp(-times / alpha), rtol=3e-3)
    # F(τ) = <v/v₀> over the decay times: the fixture's α/(α + τ).
    rows = [r for r in csv.DictReader(ln for ln in open(LIFETIMES / "analytic" / "dsam_exponential.csv")
                                      if not ln.startswith("#")) if float(r["alpha_ps"]) == alpha]
    assert rows
    q = (np.arange(4000) + 0.5) / 4000
    for r in rows:
        tau = float(r["tau_ps"])
        t = -tau * np.log1p(-q)
        f = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, t).beta.mean() / beta0
        assert f == pytest.approx(float(r["attenuation_F"]), rel=4e-3, abs=1e-3), tau
    # The kernel (no scattering) gives the same mean speeds.
    kern = tr.kernel(segs, RU, M_RU, E_RU, Z, 0.0, times, scatter=False, n=50)
    np.testing.assert_allclose(kern.beta.mean(axis=1) / beta0, np.exp(-times / alpha), rtol=3e-3)


# -- time dilation -------------------------------------------------------------------------------------------------


def test_time_dilation_enters_through_lab_time():
    gamma = 1 / math.sqrt(1 - 0.045**2)
    assert gamma == pytest.approx(1.0010140, abs=1e-7)
    stack_ = plunger(1000.0)  # 74 ps of flight
    # In the gap the nucleus's clock runs slow by γ: a proper decay time of 10 ps is 10.0101 ps in the lab.
    lab = tr.lab_time(stack_, RU, M_RU, E_RU, Z, 0.7, [10.0, 50.0])
    np.testing.assert_allclose(lab, gamma * np.array([10.0, 50.0]), rtol=1e-12)
    s = tr.slow(stack_, RU, M_RU, E_RU, Z, 0.7, lab)
    np.testing.assert_allclose(s.proper_ps, [10.0, 50.0], rtol=1e-12)
    # Slowing in the stopper: lab − proper grows by ∫ (1 − 1/γ) dt, then stays (γ = 1 at rest).
    st = stopping(RU, CU)
    rho = float(CU.density_g_cm3)
    _, t_lab = stop_time_by_quadrature(st, M_RU, E_RU, rho)
    _, t_prop = stop_time_by_quadrature(st, M_RU, E_RU, rho, proper=True)
    flight = 1000.0 / (0.045 * C_UM_PS)
    late = flight + t_lab[-1] + 5.0  # 5 ps after stopping
    s = tr.slow(stack_, RU, M_RU, E_RU, Z, 0.7, late)
    assert s.stopped[0]
    want = flight * (1 - 1 / gamma) + (t_lab[-1] - t_prop[-1])
    assert late - s.proper_ps[0] == pytest.approx(want, rel=2e-3)
    assert tr.lab_time(stack_, RU, M_RU, E_RU, Z, 0.7, s.proper_ps)[0] == pytest.approx(late, rel=1e-9)
    # Round trip in the middle of the slowing.
    proper = np.array([flight / gamma + 0.1, flight / gamma + 0.4])
    lab = tr.lab_time(stack_, RU, M_RU, E_RU, Z, 0.7, proper)
    np.testing.assert_allclose(tr.slow(stack_, RU, M_RU, E_RU, Z, 0.7, lab).proper_ps, proper, rtol=1e-6)
    assert np.all(lab > proper)


# -- ZBL scattering ------------------------------------------------------------------------------------------------


def test_the_scattering_angle_table():
    eps = np.array([1e-3, 0.1, 1.0, 10.0, 1000.0])
    b = np.array([0.01, 0.3, 1.0, 3.0, 30.0])
    e2, b2 = np.meshgrid(eps, b, indexing="ij")
    coulomb = tr._scattering_angle(e2, b2, phi=lambda x: np.ones_like(np.asarray(x, dtype=float)))
    np.testing.assert_allclose(coulomb, 2 * np.arctan(1 / (2 * e2 * b2)), rtol=1e-6)  # tan(Θ/2) = 1/(2εb)
    # The reduced nuclear stopping of the universal potential, 2ε ∫ sin²(Θ/2) b db, against ZBL's fit to it.
    inner = tr._tables()["inner"]
    for e in (1e-3, 1e-2, 0.1, 1.0, 3.0, 10.0, 100.0):
        table = 2 * e * tr._bilinear(inner, math.log(e), math.log(100.0))
        fit = math.log1p(1.1383 * e) / (2 * (e + 0.01321 * e**0.21226 + 0.19593 * math.sqrt(e))) if e <= 30 else (
            math.log(e) / (2 * e))
        assert table == pytest.approx(fit, rel=0.015 if e <= 3 else 0.06), e


def mc(segments, species, mass, energy, n, seed, times=np.inf, direction=Z, depth=0.0):
    t = np.broadcast_to(np.asarray(times, dtype=float), (n,))
    return tr._monte_carlo(segments, species, mass, np.full(n, energy), np.tile(direction, (n, 1)),
                           np.full(n, depth), t[:, None], np.random.default_rng(seed))


def test_the_monte_carlo_is_reproducible():
    segs = [tr.Segment(CU, 2.0)]
    a = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, np.full(3000, 0.3), scatter=True, rng=7)
    b = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, np.full(3000, 0.3), scatter=True, rng=np.random.default_rng(7))
    c = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, np.full(3000, 0.3), scatter=True, rng=8)
    d = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, np.full(3000, 0.3), scatter=True)
    e = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, np.full(3000, 0.3), scatter=True, rng=0)
    for key in ("energy", "direction", "segment", "stopped", "beta", "proper_ps", "position_um"):
        assert np.array_equal(getattr(a, key), getattr(b, key)), key
        assert np.array_equal(getattr(d, key), getattr(e, key)), key  # rng=None is seed 0
    assert not np.array_equal(a.direction, c.direction)
    np.testing.assert_allclose(np.linalg.norm(a.direction, axis=1), 1.0, atol=1e-9)
    # Turned, but on average still forward and symmetric about the beam.
    assert np.all(a.direction[:, 2] < 1) and a.direction[:, 2].mean() > 0.9
    assert abs(a.direction[:, 0].mean()) < 0.01 and abs(a.direction[:, 1].mean()) < 0.01


@pytest.mark.parametrize("species,mass,energy,material,thickness", [
    (RU, M_RU, E_RU, CU, 30.0),  # stops
    (RU, M_RU, E_RU, CU, 2.0),  # leaves
    (SM, M_SM, E_SM, AU, 10.0),
    ("16O", data.nuclide("16O").nuclear_mass_mev, 20.0, data.material("Mylar"), 1.0),  # a compound
])
def test_the_energy_bookkeeping(species, mass, energy, material, thickness):
    n = 4000
    r = mc([tr.Segment(material, thickness)], species, mass, energy, n, 3)
    led = r["ledger"]
    seg = r["segment"][:, 0]
    left = np.where(seg != 0, r["energy"][:, 0], 0.0)  # out of the back (1) or, backscattered, the front (−1)
    total = led["electronic"] + led["nuclear_continuous"] + led["nuclear_collisions"] + led["below_cut"] + left
    np.testing.assert_allclose(total, energy, rtol=1e-11)  # every MeV accounted for, per recoil
    assert np.all(led["electronic"] > 0) and np.all(led["nuclear_continuous"] >= -1e-9)
    st = stopping(species, material)
    stopped = seg == 0
    assert np.array_equal(r["stopped"][:, 0], stopped)
    if stopped.mean() > 0.99:
        # On average the recoils lose S_e + S_n as Stopping says: the path length is the CSDA range, and the
        # nuclear share of the loss is ∫ S_n/S dE / E (the sampled collisions are not counted twice).
        assert led["path_mg_cm2"][stopped].mean() == pytest.approx(st.range(energy), rel=0.015)
        e = np.geomspace(1e-6 * energy, energy, 4000)
        share = np.trapezoid(st.stopping_power(e, "nuclear") / st.stopping_power(e), e) / energy
        nuclear = (led["nuclear_continuous"] + led["nuclear_collisions"])[stopped].mean() / energy
        assert nuclear == pytest.approx(share, rel=0.04)
    else:
        assert (seg == 1).mean() > 0.99  # a rare recoil turns hard and stops inside
        # Leaving a layer: on average the energy after it, with a spread from the nuclear collisions.
        assert r["energy"][seg == 1, 0].mean() == pytest.approx(float(_after(st, energy, thickness)), rel=0.01)


def test_without_collisions_the_monte_carlo_is_the_straight_slowing():
    # At times short against the step's, decays in the Monte Carlo are where the straight slowing has them.
    segs = [tr.Segment(CU, 30.0)]
    t = np.array([1e-4, 2e-4, 5e-4])
    a = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, t)
    b = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, t, scatter=True, rng=1)
    np.testing.assert_allclose(b.energy, a.energy, rtol=2e-3)
    # And the mean energy at a later time is close too (collisions spread it, the mean follows the tables).
    a = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, 0.2)
    b = tr.slow(segs, RU, M_RU, E_RU, Z, 0.0, np.full(4000, 0.2), scatter=True, rng=2)
    assert b.energy.mean() == pytest.approx(a.energy[0], rel=0.02)


def test_a_plunger_with_scattering():
    n = 20000
    stack_ = plunger(30.0)
    gap = 30.0 / (0.045 * C_UM_PS)
    t = np.random.default_rng(5).exponential(2.0, n)
    a = tr.slow(stack_, RU, M_RU, E_RU, Z, 0.35, t)
    b = tr.slow(stack_, RU, M_RU, E_RU, Z, 0.35, t, scatter=True, rng=5)
    for seg in (0, 1, 2):
        assert (b.segment == seg).mean() == pytest.approx((a.segment == seg).mean(), abs=0.02), seg
    in_gap = b.segment == 1
    # In the gap: flying straight at the speed it left the target with, not stopped.
    assert not b.stopped[in_gap].any()
    assert np.all(b.position_um[in_gap, 2] > 0) and np.all(b.position_um[in_gap, 2] < 30.0 + 2.5)
    assert np.all(b.stopped[t > 2 * gap + 5.0])
    assert b.stopped.mean() == pytest.approx(a.stopped.mean(), abs=0.02)


def test_speed():
    n = 20000
    t0 = time.perf_counter()
    r = mc([tr.Segment(CU, 30.0)], RU, M_RU, E_RU, n, 1)
    took = time.perf_counter() - t0
    assert r["stopped"].all()
    assert took < 30.0, took  # about 3 s on a laptop (10⁵ recoils, followed to rest, in about 15 s)


# -- against TRIM --------------------------------------------------------------------------------------------------


def _trim() -> dict:
    path = LIFETIMES / "trim" / "transport_trim.csv"
    if not path.exists():
        pytest.skip("no TRIM reference (scripts/transport_reference.py)")
    out: dict = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(ln for ln in f if not ln.startswith("#")):
            out.setdefault((r["case"], r["run"]), {"row": r})[r["quantity"]] = float(r["value"])
    return out


CASES = {"ru98_cu": (RU, M_RU, E_RU, CU), "sm143_au": (SM, M_SM, E_SM, AU)}


def srim_stopping(species, material) -> Stopping:
    """Stopping with SRIM's electronic and nuclear stopping (``srim_stopping.csv``) in place of physim's, on
    Stopping's own grid, with its range integrated as Stopping does: what TRIM slowed the recoils with."""
    path = LIFETIMES / "trim" / "srim_stopping.csv"
    if not path.exists():
        pytest.skip("no SRIM stopping table (scripts/transport_reference.py --stopping-only)")
    with open(path, encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(ln for ln in f if not ln.startswith("#")) if r["ion"] == species]
    e = np.array([float(r["energy_mev"]) for r in rows])
    se = np.array([float(r["electronic_mev_mg_cm2"]) for r in rows])
    sn = np.array([float(r["nuclear_mev_mg_cm2"]) for r in rows])
    st = Stopping(species, material)
    ln = np.log(st._e)
    se_g = np.exp(np.interp(ln, np.log(e), np.log(se)))
    sn_g = np.exp(np.interp(ln, np.log(e), np.log(sn)))
    st._se, st._sn = se_g * 1e3, sn_g * 1e3  # MeV cm²/g
    f = st._e / (se_g + sn_g)
    st._range = np.r_[0.0, np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(ln))] + 2 * st._e[0] / (se_g[0] + sn_g[0])
    return st


@pytest.fixture
def with_srim_stopping(monkeypatch):
    """Make the transport slow recoils with SRIM's stopping (the clocks are recomputed for it)."""
    def use(species, material):
        st = srim_stopping(species, material)
        monkeypatch.setattr(tr, "stopping", lambda ion, mat: st)
        monkeypatch.setattr(tr, "_CLOCKS", {})
        return st
    return use


#: Tolerances against TRIM with the same stopping: TRIM adds electronic energy-loss straggling, which physim's Monte
#: Carlo leaves out, so its range straggling is a little larger; and TRIM drops the distant collisions (p > p_max)
#: that physim adds as a Gaussian turn, so its smallest angles (the median) are a few per cent smaller.
TOL = {"range": 0.03, "straggling": 0.15, "lateral": 0.06, "energy": 0.02, "angle": 0.06, "median": 0.10}


@pytest.mark.parametrize("case", sorted(CASES))
def test_range_and_lateral_spread_against_trim(case, with_srim_stopping):
    ref = _trim()[(case, "range")]
    species, mass, energy, material = CASES[case]
    with_srim_stopping(species, material)
    thickness = float(ref["row"]["thickness_mg_cm2"])
    r = mc([tr.Segment(material, thickness)], species, mass, energy, 20000, 11)
    pos = r["position_um"][:, 0]
    stopped = r["stopped"][:, 0]
    assert stopped.mean() > 0.99  # (a few backscattered out of the front)
    pos = pos[stopped]
    proj, strag = pos[:, 2].mean(), pos[:, 2].std()
    lateral = math.sqrt(0.5 * (pos[:, 0] ** 2 + pos[:, 1] ** 2).mean())
    assert proj == pytest.approx(ref["projected_range_um"], rel=TOL["range"])
    assert strag == pytest.approx(ref["range_straggling_um"], rel=TOL["straggling"])
    assert lateral == pytest.approx(ref["lateral_sigma_um"], rel=TOL["lateral"])


@pytest.mark.parametrize("case,run", [(c, r) for c in sorted(CASES) for r in ("quarter", "half")])
def test_angular_spread_against_trim(case, run, with_srim_stopping):
    ref = _trim()[(case, run)]
    species, mass, energy, material = CASES[case]
    with_srim_stopping(species, material)
    thickness = float(ref["row"]["thickness_mg_cm2"])
    r = mc([tr.Segment(material, thickness)], species, mass, energy, 20000, 12)
    out = r["segment"][:, 0] == 1
    angle = np.degrees(np.arccos(np.clip(r["direction"][out, 0, 2], -1, 1)))
    assert out.mean() == pytest.approx(ref["transmitted"], abs=0.02)
    assert r["energy"][out, 0].mean() == pytest.approx(ref["mean_energy_mev"], rel=TOL["energy"])
    assert math.sqrt((angle**2).mean()) == pytest.approx(ref["rms_angle_deg"], rel=TOL["angle"])
    assert np.median(angle) == pytest.approx(ref["median_angle_deg"], rel=TOL["median"])
    assert np.percentile(angle, 90) == pytest.approx(ref["angle_90_deg"], rel=TOL["angle"])


def test_the_trim_inputs_are_what_the_script_writes():
    import importlib.util

    script = Path(__file__).resolve().parents[2] / "scripts" / "transport_reference.py"
    spec = importlib.util.spec_from_file_location("transport_reference", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for *_, name, text in module.files(5000):
        path = LIFETIMES / "trim" / name
        assert path.read_bytes() == text.encode("latin-1"), name


def test_physim_stopping_against_srim_for_the_recoils():
    """Recorded, not a transport check: physim's electronic stopping of these slow heavy recoils is above SRIM's at
    low speed, where SRIM is itself below measured stopping (stopping_measured_heavy.csv: SRIM's median is 7–16%
    low for ions from Cu up below 0.1 MeV/u), so with physim's own stopping the ranges are shorter than TRIM's. This
    test pins the size of the difference, so that a change to the stopping shows up here."""
    for species, material, energy, low, high in ((RU, CU, E_RU, 0.89, 0.97), (SM, AU, E_SM, 0.76, 0.85)):
        ratio = stopping(species, material).range(energy) / srim_stopping(species, material).range(energy)
        assert low < ratio < high, (species, ratio)


# -- the DSAM kernel -----------------------------------------------------------------------------------------------


def test_the_kernel():
    segs = [tr.Segment(data.material("Sn"), 0.4, name="target"), tr.Segment(AU, 13.0, name="backing")]
    times = np.array([0.5, 0.0, 0.1, 0.2, 1.0, 3.0])
    k = tr.kernel(segs, SM, M_SM, E_SM, Z, 0.2, times, rng=1, n=4000)
    assert np.array_equal(k.times_ps, np.sort(times))
    assert k.beta.shape == (6, 4000) and k.direction.shape == (6, 4000, 3) and k.n == 4000
    np.testing.assert_allclose(k.beta[0], k.beta0, rtol=1e-6)  # at t = 0: the start
    mean_fwd = k.mean_beta_cos(0.0)
    assert np.all(np.diff(mean_fwd) < 0) and mean_fwd[-1] < 0.05 * mean_fwd[0]
    assert np.all(np.diff(k.stopped.mean(axis=1)) >= 0) and k.stopped[-1].mean() > 0.9
    # Backward detectors see the opposite shift; at 90° none on average.
    np.testing.assert_allclose(k.mean_beta_cos(180.0), -mean_fwd, rtol=1e-5)
    assert abs(k.mean_beta_cos(90.0)[2]) < 0.02 * mean_fwd[2]
    # Proper time runs slow while the nucleus moves, and at the lab rate once it has stopped.
    assert np.all(k.proper_ps <= k.times_ps[:, None] + 1e-6)
    rest = k.stopped[-2]
    assert rest.mean() > 0.5
    np.testing.assert_allclose(k.times_ps[-1] - k.proper_ps[-1, rest], k.times_ps[-2] - k.proper_ps[-2, rest],
                               atol=1e-5)
    # Without scattering every sample keeps its direction, and the mean speed is the straight slowing's.
    k0 = tr.kernel(segs, SM, M_SM, E_SM, Z, 0.2, times, scatter=False, n=10)
    np.testing.assert_allclose(k0.direction, np.broadcast_to(Z, k0.direction.shape))
    s = tr.slow(segs, SM, M_SM, E_SM, Z, 0.2, np.sort(times))
    np.testing.assert_allclose(k0.beta[:, 0], s.beta, rtol=1e-6)
    # Scattering slows the forward component faster than the straight path does (the recoils turn).
    assert np.all(mean_fwd[2:4] < k0.beta[2:4, 0])
    # Starting states drawn from arrays.
    e = np.array([E_SM, 0.8 * E_SM])
    k2 = tr.kernel(segs, SM, M_SM, e, Z, [0.1, 0.3], [0.0], rng=2, n=1000)
    assert len(np.unique(k2.beta0)) == 2 and len(np.unique(k2.direction0, axis=0)) == 1


# -- the kernel in lifetimes.py's DSAM functions -------------------------------------------------------------------


def test_the_transport_kernel_gives_the_analytic_attenuation_through_lifetimes(monkeypatch):
    from physim.nuclear import lifetimes

    alpha = 0.2  # ps: stopping proportional to the speed, as in the analytic fixture
    rho = float(CU.density_g_cm3)
    k = 2 / (alpha * 1e3 * rho * 2.99792458e-2 * math.sqrt(2 / M_RU))
    fake = _VelocityStopping(k, M_RU / (data.U_KEV * 1e-3))
    monkeypatch.setattr(tr, "stopping", lambda ion, material: fake)
    monkeypatch.setattr(tr, "_CLOCKS", {})
    times = np.linspace(0.0, 40 * alpha, 8001)
    kern = tr.kernel([tr.Segment(CU, 100.0)], RU, M_RU, E_RU, Z, 0.0, times, scatter=False, n=1)
    assert isinstance(kern, lifetimes.Kernel)
    rows = [r for r in csv.DictReader(ln for ln in open(LIFETIMES / "analytic" / "dsam_exponential.csv")
                                      if not ln.startswith("#")) if float(r["alpha_ps"]) == alpha]
    for r in rows:
        tau = float(r["tau_ps"])
        f = lifetimes.dsam_attenuation(tau, kern)
        assert f == pytest.approx(float(r["attenuation_F"]), rel=4e-3, abs=1e-3), tau


def _centroid_and_width(row, edges):
    mid = 0.5 * (edges[1:] + edges[:-1])
    mean = float(row @ mid / row.sum())
    return mean, float(np.sqrt(row @ (mid - mean) ** 2 / row.sum()))


def test_a_scattered_kernel_s_line_shape():
    from physim.nuclear import lifetimes

    segs = [tr.Segment(data.material("Sn"), 0.4, name="target"), tr.Segment(AU, 13.0, name="backing")]
    times = np.linspace(0.0, 3.0, 601)
    kern = tr.kernel(segs, SM, M_SM, E_SM, Z, 0.2, times, rng=3, n=4000)
    tau, e0 = 0.3, 1000.0
    edges = np.arange(960.0, 1040.0001, 0.05)
    decay = np.exp(-times / tau) / tau
    for theta in (30.0, 150.0, 90.0):
        row = lifetimes.lineshape(e0, [tau], None, [1.0], kern, theta, edges)[0]
        assert row.sum() == pytest.approx(1.0, abs=1e-9)
        # The centroid's shift is the mean β·cos(angle to the detector) at emission (to first order in β).
        m = kern.mean_beta_cos(theta).astype(float)
        expected = np.trapezoid(m * decay, times) + m[-1] * math.exp(-times[-1] / tau)
        centroid, width = _centroid_and_width(row, edges)
        assert (centroid - e0) / e0 == pytest.approx(expected, abs=0.05 * abs(expected) + 4e-4), theta
        if theta == 90.0:
            # Sideways motion counts: at 90° the line is broadened, not a δ at E₀ √(1 − β_z²) as it would be for
            # recoils along the beam.
            beam_only = lifetimes.SampledKernel(kern.times_ps, kern.beta_cos(0.0))
            row_z = lifetimes.lineshape(e0, [tau], None, [1.0], beam_only, theta, edges)[0]
            _, width_z = _centroid_and_width(row_z, edges)
            assert width > 1.0 and width > 5 * width_z, (width, width_z)
    # The cached binning is the detector's: two angles give two different lines.
    a = lifetimes.lineshape(e0, [tau], None, [1.0], kern, 30.0, edges)[0]
    b = lifetimes.lineshape(e0, [tau], None, [1.0], kern, -30.0, edges)[0]
    assert _centroid_and_width(a, edges)[0] == pytest.approx(_centroid_and_width(b, edges)[0], abs=0.3)
    assert not np.array_equal(a, b)  # the scattered recoils are not exactly symmetric about the beam


def test_a_straight_kernel_is_the_equivalent_sampled_kernel():
    from physim.nuclear import lifetimes

    segs = [tr.Segment(data.material("Sn"), 0.4, name="target"), tr.Segment(AU, 13.0, name="backing")]
    times = np.linspace(0.0, 3.0, 301)
    kern = tr.kernel(segs, SM, M_SM, E_SM, Z, 0.2, times, scatter=False, n=5)
    sampled = lifetimes.SampledKernel(kern.times_ps, kern.beta_cos(0.0))
    edges = np.arange(960.0, 1040.0001, 0.1)
    for theta in (0.0, 30.0, 90.0, 150.0):
        np.testing.assert_allclose(lifetimes.lineshape(1000.0, [0.3, 1.0], None, [1.0, 0.5], kern, theta, edges),
                                   lifetimes.lineshape(1000.0, [0.3, 1.0], None, [1.0, 0.5], sampled, theta, edges),
                                   rtol=1e-12, atol=1e-15)
    assert lifetimes.dsam_attenuation(0.3, kern) == pytest.approx(lifetimes.dsam_attenuation(0.3, sampled), rel=1e-12)
