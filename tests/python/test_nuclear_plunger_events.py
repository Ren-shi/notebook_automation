"""Coulomb-excitation γ rays through a plunger and a DSAM backing (lifetimes wave 2A): the emitter follows the recoil
stack (``rates.recoil_stack``), decays after a time drawn in its own frame from the level's lifetime, and its γ ray is
shifted, degraded or unshifted; with deorientation in vacuum the angular distribution is attenuated by G_k(t).

The yardsticks are :mod:`physim.nuclear.lifetimes`: the unshifted share against distance (``unshifted``), the peak
separation (``peak_separation``) and the attenuations of the two components (``deorientation``,
``deorientation_shifted``)."""

import copy
import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import lifetimes as lt
from physim.nuclear.gamma import excitation_of
from physim.nuclear.gamma_events import simulate_gammas

E0_KEV = 652.46
DISTANCES = (6.0, 40.0, 160.0)


def _p(x, k):
    return {2: (3 * x * x - 1) / 2, 4: (35 * x**4 - 30 * x * x + 3) / 8}[k]


@pytest.fixture(scope="module")
def ru98():
    return Experiment.example("coulex_ru98_plunger")


@pytest.fixture(scope="module")
def tau(ru98):
    """The 2⁺'s mean life as the events use it: the setup's B(E2) through lifetimes.level_taus."""
    taus = lt.level_taus(excitation_of(ru98).scheme)
    (tau_ps, source), = taus.values()
    assert source.startswith("B(E2)") and 8.0 < tau_ps < 9.5  # Radeck et al.: 8.36(29) ps
    return tau_ps


@pytest.fixture(scope="module")
def series(ru98):
    """A small distance series: the γ rays at each distance."""
    return {d: simulate_gammas(ru98, events=40_000, seed=7, distance_um=d) for d in DISTANCES}


def test_each_gamma_ray_says_where_it_was_emitted(ru98, series):
    g = series[40.0]
    assert g.distance_um == 40.0 and g.segments == ["target", "gap", "stopper"]
    assert np.all(g["distance_um"] == 40.0)
    for name in ("stopped", "segment", "time_ps", "vacuum_ps", "stopper_ps"):
        assert len(g[name]) == len(g), name
    # Decays in the target and the gap fly, those at rest are in the stopper (the Cu stops the recoil).
    seg, stopped = g["segment"], g["stopped"]
    assert set(np.unique(seg)) <= {0, 1, 2} and np.all(seg[stopped] == 2)
    assert np.all(g["beta"][stopped] == 0) and np.all(g["beta"][seg == 1] > 0.03)
    assert np.all(np.isclose(g["energy_lab"][stopped], g["energy0"][stopped]))
    # Time in vacuum: none in the target, the whole flight once in the stopper.
    assert np.all(g["vacuum_ps"][seg == 0] == 0)
    in_stopper = seg == 2
    assert np.all(g["vacuum_ps"][in_stopper] < g["stopper_ps"][in_stopper])
    assert np.all(g["time_ps"] >= 0)
    assert any("Plunger at 40 µm" in n for n in g.notes)


def test_the_unshifted_share_follows_the_decay_curve(series, tau):
    """Per distance, the share of the γ rays emitted after the recoil reached the stopper is the mean of
    lifetimes.unshifted over the recoils' own flight times, within statistics."""
    shares = []
    for d in DISTANCES:
        g = series[d]
        after = g["segment"] >= 2
        t = g["stopper_ps"]
        assert np.all(np.isfinite(t))
        expected = float(np.mean(lt.unshifted([tau], None, [1.0], t)[0]))
        n = len(g)
        p = float(after.mean())
        sigma = math.sqrt(expected * (1 - expected) / n)
        assert abs(p - expected) < 4 * sigma + 0.003, (d, p, expected, sigma)
        shares.append(p)
        # The decay curve of one level with the setup's typical recoil, as the planner would draw it.
        beta = float(np.median(g["beta"][g["segment"] == 1]))
        typical = lt.unshifted([tau], None, [1.0], lt.flight_times(d + 8.0, beta, proper=True))[0]
        assert abs(p - typical) < 0.06, (d, p, typical)
    assert shares[0] > shares[1] > shares[2] > 0.0


def test_the_shifted_peak_sits_at_the_separation_of_each_ring(series):
    """The in-flight (gap) full-energy peak lies above (forward clovers, 41.5°) or below (backward, 138.5°) the
    unshifted one by lifetimes.peak_separation at the ring's angle, for the recoils' own speeds (within 6%: the
    emitters fly up to 14° off the beam and the crystals span ±7°)."""
    g = series[40.0]
    ring = np.array([g.crystals[c].name[0] for c in g["crystal"]])
    full = g["counted"] & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-9)
    for name, theta in (("F", 41.5), ("B", 138.5)):
        shifted = full & (ring == name) & (g["segment"] == 1)
        unshifted = full & (ring == name) & g["stopped"]
        assert shifted.sum() > 100 and unshifted.sum() > 300
        seen = 1e3 * (g["measured"][shifted].mean() - g["measured"][unshifted].mean())
        expected = float(np.mean(lt.peak_separation(E0_KEV, g["beta"][shifted], theta)))
        assert seen == pytest.approx(expected, rel=0.06), (name, seen, expected)
        assert abs(1e3 * g["measured"][unshifted].mean() - E0_KEV) < 0.3
    # The shifted line is far enough from the unshifted one for the 2.1 keV clovers to separate them.
    assert lt.peak_separation(E0_KEV, 0.04, 41.5) > 5 * 2.1


def _angular_setups(ru98, deorientation=True):
    """⁹⁸Ru on Mg with the Mg recoils gated at 18.8°–20.8° (head-on collisions: an alignment along the beam) and
    seven single crystals over θ that do not shade one another (alternate sides); with and without the paper's
    deorientation."""
    d = copy.deepcopy(ru98.to_dict())
    d["detectors"][0]["inner_radius"] = "8.5 mm"
    d["detectors"][0]["outer_radius"] = "9.5 mm"
    d["gamma_detectors"] = [{"name": f"G{a}", "theta": f"{a} deg", "phi": f"{180 * (i % 2)} deg", "distance": "80 mm",
                             "radius": "25 mm", "resolution": "2.1 keV"}
                            for i, a in enumerate((25, 45, 65, 90, 115, 135, 155))]
    plain = copy.deepcopy(d)
    del plain["lifetime"]["deorientation"]
    return Experiment.from_dict(d), Experiment.from_dict(plain)


def _attenuation(a, b, component, angles):
    """G₂, G₄ (and errors) of the deoriented run ``a`` against the plain run ``b`` (same events), from the counts per
    crystal of one component: with the plain run's W(θ) = 1 + A₂P₂ + A₄P₄ (fitted), the ratio r = W'/W gives
    rW − 1 = G₂A₂P₂ + G₄A₄P₄, linear in G."""
    x = np.cos(np.radians(angles))
    na = np.array([np.sum((a["crystal"] == k) & component(a)) for k in range(len(angles))], dtype=float)
    nb = np.array([np.sum((b["crystal"] == k) & component(b)) for k in range(len(angles))], dtype=float)
    m = np.c_[np.ones_like(x), _p(x, 2), _p(x, 4)]
    c, *_ = np.linalg.lstsq(m / np.sqrt(nb)[:, None], nb / np.sqrt(nb), rcond=None)
    a2, a4 = c[1] / c[0], c[2] / c[0]
    w = 1 + a2 * _p(x, 2) + a4 * _p(x, 4)
    r = na / nb
    y, sy = r * w - 1, r * np.sqrt(1 / na + 1 / nb) * w
    design = np.c_[a2 * _p(x, 2), a4 * _p(x, 4)] / sy[:, None]
    g, *_ = np.linalg.lstsq(design, y / sy, rcond=None)
    err = np.sqrt(np.diag(np.linalg.inv(design.T @ design)))
    return {2: g[0], 4: g[1]}, {2: err[0], 4: err[1]}, a2


def test_deorientation_attenuates_the_shifted_and_unshifted_components_differently(ru98, tau):
    """At 100 µm (a flight of about 7 ps) the unshifted γ rays carry G_k at the flight time, the shifted ones G̃_k,
    the mean over the decays in flight (Radeck et al. Eq. 7): the angular distributions of the two components
    are attenuated as lifetimes.deorientation and lifetimes.deorientation_shifted say."""
    deoriented, plain = _angular_setups(ru98)
    angles = (25, 45, 65, 90, 115, 135, 155)
    a = simulate_gammas(deoriented, events=100_000, seed=11, distance_um=100.0)
    b = simulate_gammas(plain, events=100_000, seed=11, distance_um=100.0)
    assert np.all(a.events.columns["recoil"][a["particle"]])
    dor = deoriented.lifetime.deorientation
    flight = float(np.median(a["vacuum_ps"][a["segment"] >= 2]))
    unshifted = lt.deorientation(flight, dor.alpha, dor.gamma_per_ps)
    shifted = lt.deorientation_shifted(flight, tau, dor.alpha, dor.gamma_per_ps)
    # Each γ ray's own G_k: at the flight time once in the stopper, averaged over the flight for those in the gap.
    g = lt.deorientation(a["vacuum_ps"], dor.alpha, dor.gamma_per_ps)
    for k in (2, 4):
        assert np.mean(g[k][a["segment"] >= 2]) == pytest.approx(unshifted[k], abs=0.01)
        assert np.mean(g[k][a["segment"] == 1]) == pytest.approx(shifted[k], abs=0.015)
    seen_u, err_u, a2 = _attenuation(a, b, lambda e: e["segment"] >= 2, angles)
    seen_s, err_s, _ = _attenuation(a, b, lambda e: e["segment"] == 1, angles)
    assert a2 > 0.3  # a strong alignment to see it on
    for k in (2, 4):
        assert abs(seen_u[k] - unshifted[k]) < 3.5 * err_u[k], (k, seen_u, unshifted, err_u)
        assert abs(seen_s[k] - shifted[k]) < 3.5 * err_s[k], (k, seen_s, shifted, err_s)
    # The two components differ: the bias the correction W_US/W_SH is for.
    assert seen_s[2] - seen_u[2] > 3 * math.hypot(err_s[2], err_u[2])
    # Without deorientation nothing is attenuated: the plain run's own G is 1.
    assert "Deorientation" not in " ".join(b.notes) and "Deorientation" in " ".join(a.notes)


def test_a_differential_plunger_gives_three_components(ru98):
    """A degrader then a stopper (``kind = "differential"``): γ rays in flight before the degrader (fully shifted),
    after it (degraded: less shifted) and at rest."""
    d = copy.deepcopy(ru98.to_dict())
    d["lifetime"]["stopper"] = {"kind": "differential", "material": "Mg", "thickness": "2 mg/cm2"}
    d["lifetime"]["degrader2"] = {"material": "Cu", "thickness": "15.7 mg/cm2", "gap": "60 um"}
    del d["lifetime"]["deorientation"]
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, events=30_000, seed=5, distance_um=40.0)
    assert g.segments == ["target", "gap", "degrader", "gap2", "degrader2"]
    seg = g["segment"]
    first, second = seg == 1, seg == 3
    assert first.sum() > 100 and second.sum() > 100 and g["stopped"].sum() > 100
    assert np.all(g["stopped"] <= (seg == 4))
    # The degrader slows the recoil: slower after it, and the γ rays less shifted.
    assert g["beta"][second].mean() < g["beta"][first].mean() - 0.003
    fwd = np.array([g.crystals[c].name[0] == "F" for c in g["crystal"]])
    shift = g["energy_lab"] * 1e3 - E0_KEV
    assert 0 < shift[fwd & second].mean() < shift[fwd & first].mean()
    # Time in vacuum: both gaps add up for the γ rays at rest.
    assert g["vacuum_ps"][g["stopped"]].mean() > g["vacuum_ps"][second].mean()


def test_dsam_slows_and_scatters_in_the_backing(ru98):
    """DSAM (``method = "dsam"``): the recoil slows in a gold backing, with the ZBL scattering, so the γ rays range
    from the full shift to none. (⁹⁸Ru's 8.6 ps 2⁺ is far too long for DSAM, and the forward-scattered ⁹⁸Ru are too
    fast for this backing: most leave it and decay in flight; the test checks the bookkeeping, not a measurement.)"""
    d = copy.deepcopy(ru98.to_dict())
    d["lifetime"] = {"method": "dsam"}
    d["target"]["backing"] = {"material": "Au", "thickness": "15 mg/cm2"}
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, events=6_000, seed=2, gammas_per_event=2)
    # No plunger: the distance is NaN (0 would be a real 0 µm), and so is the time to a stopper.
    assert g.distance_um is None and np.all(np.isnan(g["distance_um"])) and np.all(np.isnan(g["stopper_ps"]))
    assert g.segments == ["target", "backing"]
    # The fastest ⁹⁸Ru (scattered forward at small CM angles) can leave by the back: segment 2, beyond the backing.
    assert set(np.unique(g["segment"])) <= {0, 1, 2}
    stopped, seg = g["stopped"], g["segment"]
    assert stopped.any() and np.all(seg[stopped] == 1) and np.all(g["beta"][stopped] == 0)
    beyond = seg == 2
    assert beyond.sum() > 0.5 * len(g)
    # Those that left decay at their exit speed, below the speed they were made with; time in vacuum after the
    # backing only.
    vacuum = g["vacuum_ps"]
    assert np.all(g["beta"][beyond] < 0.085) and np.mean(vacuum[beyond] > 0) > 0.98
    assert np.all(vacuum[seg < 2] == 0) and np.all(vacuum <= g["time_ps"] + 1e-9)
    assert any(n.startswith("DSAM") for n in g.notes)
    # Without the lifetime section the backing still slows the recoil only through half-lives (item 73): none here.
    del d["lifetime"]
    g0 = simulate_gammas(Experiment.from_dict(d), events=6_000, seed=2, gammas_per_event=2)
    assert not g0["stopped"].any() and "time_ps" in g0.columns and np.all(g0["time_ps"] == 0)
    # Without a lifetime section the plunger's columns are not stored, but read as NaN.
    assert "distance_um" not in g0.columns and np.all(np.isnan(g0["distance_um"]))


def test_deorientation_scales_the_rank_k_parts_of_a_density_matrix():
    """In a cascade the emitting state's density matrix is attenuated: its rank-k statistical tensors × G_k (trace
    and hermiticity kept), so the γ ray's distribution is W_iso + G₂ W₂ + G₄ W₄."""
    from physim.nuclear.correlation import Emission
    from physim.nuclear.gamma_events import _attenuate

    rho = np.diag([0.05, 0.15, 0.6, 0.15, 0.05]).astype(complex)[None]   # J = 2, aligned along z
    rng = np.random.default_rng(3)
    dirs = rng.normal(size=(200, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    em = Emission(2.0, 0.0, 2)

    def w(g2, g4):
        r = _attenuate(np.repeat(rho, len(dirs), axis=0), {2: np.full(len(dirs), g2), 4: np.full(len(dirs), g4)})
        assert np.allclose(np.einsum("naa->n", r).real, 1.0) and np.allclose(r, np.conj(np.transpose(r, (0, 2, 1))))
        return em.weights(r, dirs)

    iso, full = w(0.0, 0.0), w(1.0, 1.0)
    assert np.allclose(iso, iso[0]) and np.allclose(full, em.weights(np.repeat(rho, len(dirs), axis=0), dirs))
    w2, w4 = w(1.0, 0.0) - iso, w(0.0, 1.0) - iso
    assert np.allclose(w(0.5, 0.2), iso + 0.5 * w2 + 0.2 * w4)
    # The rank-2 part of an alignment along z goes as P₂(cos θ), the rank-4 part as P₄.
    x = dirs[:, 2]
    assert np.allclose(w2 / _p(x, 2), (w2 / _p(x, 2))[0], rtol=1e-6, atol=1e-9)
    assert np.allclose(w4 / _p(x, 4), (w4 / _p(x, 4))[0], rtol=1e-6, atol=1e-9)


def test_a_cascade_through_a_plunger_adds_up_its_decay_times():
    """⁶²Ni's level scheme (ENSDF half-lives) behind a plunger, with deorientation: each cascade's γ rays leave in
    order, and a level fed from above decays later than the same level made directly."""
    d = copy.deepcopy(Experiment.example("coulex_ni62_pt194").to_dict())
    d["lifetime"] = {"method": "rdds", "distances": ["50 um"],
                     "deorientation": {"alpha2": 0.2, "gamma2": "0.2 /ps", "alpha4": 0.3, "gamma4": "0.3 /ps"},
                     "stopper": {"kind": "stopper", "material": "Au", "thickness": "20 mg/cm2"}}
    # The sample the test was written with: its 13 fed 2⁺ states are few, and other samples of this size (4 to 10 of them,
    # a mean time of 3 ± 1 ps against 2.2) are not always enough for the last comparison.
    d["reaction"]["e1_polarisation"] = 0
    exp = Experiment.from_dict(d)
    g = simulate_gammas(exp, events=40_000, seed=3, distance_um=50.0)
    own = g["nucleus"] == 0
    # Backscattered ⁶²Ni leave by the target's front face (segment −1) and decay in flight there.
    assert own.sum() > 200 and set(np.unique(g["segment"])) <= {-1, 0, 1, 2}
    assert not g["stopped"][g["segment"] == -1].any()
    # Within a cascade the times grow step by step.
    cas, t = g["cascade"][own], g["time_ps"][own]
    order = np.lexsort((t, cas))
    same = cas[order][1:] == cas[order][:-1]
    assert same.any()
    lv = g["initial"][own][order]
    assert np.all(lv[1:][same] < lv[:-1][same])   # later γ rays come from lower levels
    # The 2⁺₁ fed from above decays later than the 2⁺₁ made directly.
    first_2 = (g["initial"][own] == 1)
    fed = np.isin(g["cascade"][own], cas[order][1:][same])
    assert t[first_2 & fed].mean() > t[first_2 & ~fed].mean()
    assert any("Deorientation" in n for n in g.notes)


def test_a_plunger_without_a_distance_takes_its_first(ru98):
    g = simulate_gammas(ru98, events=5_000, seed=1)
    assert g.distance_um == 3.0 and any("the setup's first" in n for n in g.notes)
    with pytest.raises(ValueError, match="negative"):
        simulate_gammas(ru98, events=5_000, seed=1, distance_um=-1.0)
