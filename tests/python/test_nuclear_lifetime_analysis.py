"""Lifetime analysis (planning/lifetimes.md, wave 2C): peak areas, decay-curve fits, the DDCM (singles, gated,
three-foil), DSAM line-shape fits, B from τ and deorientation, on spectra built from lifetimes.py (and
transport.kernel with fake stopping) with Poisson noise, and against wave 0's references in
tests/reference/nuclear/lifetimes/ (closed forms, NaPyTau, AhKin, ⁹⁸Ru's Table I)."""

import csv
import inspect
import math
import re
import tomllib
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import data
from physim.nuclear import lifetime_analysis as la
from physim.nuclear import lifetimes as lt
from physim.nuclear import transport as tr
from physim.nuclear.levels import Level, LevelScheme, Transition, Value
from physim.nuclear.stopping import Stopping

REF = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "lifetimes"
AHKIN = REF / "ahkin" / "y83_in_ta"


def rows(name: str) -> list:
    with open(REF / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader([line for line in f if not line.startswith("#")]))


def sqrt_counts(x):
    return np.sqrt(np.maximum(np.asarray(x, dtype=float), 1.0))


def within(result, truth, k=2.0):
    """τ back within k of its quoted (statistical) uncertainty."""
    return abs(result.tau_ps - truth) <= k * result.stat_ps


# ⁹⁸Ru-like plunger (Radeck et al. 2012): distances, speed, lifetime of the 2⁺₁.
D_RU = np.array([3, 8, 15, 25, 40, 60, 90, 130, 200, 300, 450], dtype=float)
BETA_RU, TAU_RU = 0.045, 8.36


# ---------------------------------------------------------------------------------------------------------------
# The interface


def test_signatures_match_the_frozen_interface():
    def params(f):
        return [p for p, v in inspect.signature(f).parameters.items() if v.kind is not v.KEYWORD_ONLY]

    assert params(la.peak_areas) == ["spectra_by_distance", "e0_kev", "beta", "theta_deg", "fwhm_kev"]
    assert params(la.decay_curve) == ["distances_um", "unshifted", "errors", "beta", "feeding"]
    assert params(la.ddcm) == ["distances_um", "unshifted", "shifted", "errors", "beta", "feeder_unshifted",
                               "feeder_shifted", "alpha", "polynomial"]
    assert inspect.signature(la.ddcm).parameters["polynomial"].default == "piecewise"
    assert inspect.signature(la.ddcm).parameters["alpha"].default == 1.0
    assert params(la.fit_lineshape) == ["counts", "edges_kev", "e0_kev", "kernel", "theta_deg", "feeding",
                                        "background"]
    assert inspect.signature(la.fit_lineshape).parameters["background"].default == "linear"
    assert params(la.b_from_result) == ["result", "scheme", "level"]
    assert params(la.deorientation_correction) == ["theta_deg", "a_k", "g_unshifted", "g_shifted"]
    assert params(la.fit_deorientation) == ["distances_um", "angular_unshifted", "angular_shifted", "a_k0",
                                            "tau_ps", "beta"]
    fields = [f for f in la.TauResult.__dataclass_fields__]
    assert fields[:6] == ["tau_ps", "stat_ps", "budget", "method", "per_distance", "notes"]


def test_the_budget_names_every_source():
    u = np.exp(-lt.flight_times(D_RU, BETA_RU, proper=True) / TAU_RU)
    r = la.decay_curve(D_RU, u, np.full_like(u, 0.01), BETA_RU, offset_error_um=1.0, beta_error=0.0006)
    assert set(la.BUDGET_KEYS) <= set(r.budget)
    assert r.budget["statistics"] == r.stat_ps
    assert r.budget["stopping"] == 0.0
    # τ ∝ 1/β for a plunger: the velocity entry is τ σβ/β.
    assert r.budget["velocity"] == pytest.approx(TAU_RU * 0.0006 / BETA_RU, rel=0.02)
    # An offset of 1 µm moves τ by a fraction of a ps: of the order of 1 µm / v (the share is pinned to 1 at d = 0,
    # so it is not exactly that).
    v = BETA_RU * lt.C_UM_PER_PS / lt.lorentz_gamma(BETA_RU)
    assert 0.3 / v < r.budget["distance"] < 3.0 / v
    assert r.total_ps == pytest.approx(math.sqrt(sum(x * x for x in r.budget.values())))


# ---------------------------------------------------------------------------------------------------------------
# Peak areas


def _plunger_spectra(rng, e0, beta, theta, fwhm, counts=40000.0, background=20.0):
    t = lt.flight_times(D_RU, beta, proper=True)
    u = np.exp(-t / TAU_RU)
    edges = np.arange(e0 - 40.0, e0 + 60.0, 0.5)
    centres = 0.5 * (edges[1:] + edges[:-1])
    e_s = e0 + lt.peak_separation(e0, beta, theta)
    sig_u, sig_s = fwhm / 2.3548, 1.3 * fwhm / 2.3548
    spectra = []
    for share in u:
        mu = background + 0.5 * counts * (
            share * np.exp(-0.5 * ((centres - e0) / sig_u) ** 2) / (sig_u * math.sqrt(2 * math.pi))
            + (1 - share) * np.exp(-0.5 * ((centres - e_s) / sig_s) ** 2) / (sig_s * math.sqrt(2 * math.pi)))
        spectra.append((rng.poisson(mu).astype(float), edges))
    return spectra, counts * u, counts * (1 - u)


@pytest.mark.parametrize("theta", [36.5, 143.5, 85.0])
def test_peak_areas_recover_the_components(theta):
    rng = np.random.default_rng(11)
    e0, fwhm = 652.0, 2.2
    spectra, true_u, true_s = _plunger_spectra(rng, e0, BETA_RU, theta, fwhm)
    areas = la.peak_areas(spectra, e0, BETA_RU, theta, fwhm)
    assert len(areas) == len(D_RU)
    s = np.array([a[0] for a in areas])
    u = np.array([a[1] for a in areas])
    err = np.array([a[2] for a in areas])
    pulls = np.concatenate([(s - true_s) / err[:, 0], (u - true_u) / err[:, 1]])
    assert np.all(np.abs(pulls) < 4), pulls
    assert abs(pulls.mean()) < 1.0
    # From the areas, the DDCM gives the lifetime back.
    r = la.ddcm(D_RU, u, s, (err[:, 1], err[:, 0]), BETA_RU)
    assert within(r, TAU_RU, 2.5), (r.tau_ps, r.stat_ps)


def test_peak_areas_take_shared_edges():
    rng = np.random.default_rng(12)
    spectra, _, _ = _plunger_spectra(rng, 652.0, BETA_RU, 36.5, 2.2)
    edges = spectra[0][1]
    a = la.peak_areas(spectra, 652.0, BETA_RU, 36.5, 2.2)
    b = la.peak_areas([c for c, _ in spectra], 652.0, BETA_RU, 36.5, 2.2, edges_kev=edges)
    assert a == b
    with pytest.raises(ValueError):
        la.peak_areas([c for c, _ in spectra], 652.0, BETA_RU, 36.5, 2.2)


# ---------------------------------------------------------------------------------------------------------------
# Decay curves: the closed forms and closure


def test_decay_curve_recovers_the_one_level_closed_form():
    table = rows("analytic/rdds_one_level.csv")
    cases = sorted({(float(r["tau_ps"]), float(r["beta"])) for r in table})
    for tau, beta in cases:
        sel = [r for r in table if float(r["tau_ps"]) == tau and float(r["beta"]) == beta
               and 0.02 * tau <= float(r["flight_time_ps"]) <= 5 * tau]
        d = np.array([float(r["distance_um"]) for r in sel])
        for column, proper in (("unshifted", True), ("unshifted_no_dilation", False)):
            u = np.array([float(r[column]) for r in sel])
            res = la.decay_curve(d, u, 0.01 * u + 1e-4, beta, proper=proper)
            assert res.tau_ps == pytest.approx(tau, rel=1e-6), (tau, beta, column)
            assert res.method == "decay curve"


def test_decay_curve_recovers_the_two_level_closed_form():
    table = rows("analytic/rdds_two_levels.csv")
    cases = sorted({tuple(float(r[k]) for k in ("tau1_ps", "tau2_ps", "p1", "p2", "beta")) for r in table})
    assert len(cases) >= 4
    for tau1, tau2, p1, p2, beta in cases:
        sel = [r for r in table if (float(r["tau1_ps"]), float(r["tau2_ps"]), float(r["p1"]), float(r["p2"]),
                                    float(r["beta"])) == (tau1, tau2, p1, p2, beta)]
        d = np.array([float(r["distance_um"]) for r in sel])
        lower = np.array([float(r["unshifted_lower"]) for r in sel])
        upper = np.array([float(r["unshifted_upper"]) for r in sel])
        keep = lower > 1e-8
        # The fed level, its feeder's τ known.
        feed = la.Feeding([tau1, 1.0], [[0, 1], [0, 0]], [p1, p2], level=1)
        res = la.decay_curve(d[keep], lower[keep], 0.01 * lower[keep] + 1e-6, beta, feed, proper=False)
        assert res.tau_ps == pytest.approx(tau2, rel=1e-5), (tau1, tau2, p1, p2)
        # The feeder on its own.
        keep = upper > 1e-8
        res = la.decay_curve(d[keep], upper[keep], 0.01 * upper[keep] + 1e-6, beta, proper=False)
        assert res.tau_ps == pytest.approx(tau1, rel=1e-5)


def test_decay_curve_closure_with_noise_and_a_free_side_feeding_time():
    rng = np.random.default_rng(21)
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    for seed_tau, tau_sf in ((TAU_RU, 2.0),):
        share = lt.unshifted([tau_sf, seed_tau], [[0, 1], [0, 0]], [0.3, 0.7], t)[1]
        n = rng.poisson(30000 * share)
        total = 30000
        r_meas = n / total
        err = np.sqrt(np.maximum(share * (1 - share), 1e-4) / total)
        fixed = la.decay_curve(D_RU, r_meas, err, BETA_RU, la.Feeding.side(tau_sf, 0.3, side_error_ps=0.5))
        assert within(fixed, seed_tau), (fixed.tau_ps, fixed.stat_ps)
        assert fixed.budget["feeding"] > 0
        free = la.decay_curve(D_RU, r_meas, err, BETA_RU, la.Feeding.side(1.0, 0.3, free=True))
        assert within(free, seed_tau, 2.5), (free.tau_ps, free.stat_ps)
        assert free.params["tau_0_ps"] == pytest.approx(tau_sf, rel=0.5)
        assert free.stat_ps > fixed.stat_ps  # a free side feeding costs precision


def test_decay_curve_with_a_free_normalisation():
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    u = 0.93 * np.exp(-t / TAU_RU)
    r = la.decay_curve(D_RU, u, np.full_like(u, 0.005), BETA_RU, free_norm=True)
    assert r.tau_ps == pytest.approx(TAU_RU, rel=1e-6)
    assert r.params["norm"] == pytest.approx(0.93, rel=1e-6)


# ---------------------------------------------------------------------------------------------------------------
# DDCM


def test_ddcm_agrees_with_napytau():
    """The same input as NaPyTau's (analytic/rdds_one_level.csv, flight times 0.05–3τ, 1% errors), one polynomial
    over all distances fitted unweighted to the shifted intensity: the τ(x) agree with NaPyTau's to its six
    figures. Its final τ is not a usable reference (its uncertainty propagation is broken, so its weighted mean is
    the first τ(x)); ours, the error-weighted mean, is at least as close to the true τ."""
    one = rows("analytic/rdds_one_level.csv")
    refs = rows("napytau/ddcm_one_level.csv")
    assert len(refs) == 9
    for r in refs:
        tau, beta, deg = float(r["tau_ps"]), float(r["beta"]), int(r["polynomial_degree"])
        sel = [x for x in one if float(x["tau_ps"]) == tau and float(x["beta"]) == beta
               and 0.05 * tau <= float(x["flight_time_ps"]) <= 3 * tau]
        d = np.array([float(x["distance_um"]) for x in sel])
        u = np.array([float(x["unshifted_no_dilation"]) for x in sel])
        s = 1 - u
        res = la.ddcm(d, u, s, (0.01 * u + 1e-4, 0.01 * s + 1e-4), beta, polynomial=deg, fit="shifted",
                      weighted=False, proper=False)
        mine = np.array([p["tau_ps"] for p in res.per_distance])
        theirs = np.array([float(x) for x in r["tau_i_ps"].split(";")])
        np.testing.assert_allclose(mine, theirs, rtol=2e-5)
        assert abs(res.tau_ps - tau) <= abs(float(r["tau_final_ps"]) - tau) + 1e-3 * tau
        if deg == 4:  # the README's central values: within 0.3% of the true τ (NaPyTau's degree 4)
            assert res.tau_ps == pytest.approx(tau, rel=3e-3)
            assert float(r["tau_final_ps"]) == pytest.approx(tau, rel=5e-3)


@pytest.mark.parametrize("polynomial", ["piecewise", "smooth", 4])
def test_ddcm_on_the_closed_form(polynomial):
    """Exact decay curves (no noise): the joint fit follows them to a fraction of a percent with each curve."""
    for tau, beta in ((10.0, 0.05), (100.0, 0.05), (10.0, 0.10)):
        t = np.geomspace(0.03 * tau, 4 * tau, 11)
        d = t * beta * lt.C_UM_PER_PS
        u = np.exp(-t / tau)
        s = 1 - u
        res = la.ddcm(d, u, s, (0.01 * u + 1e-4, 0.01 * s + 1e-4), beta, polynomial=polynomial, proper=False)
        assert res.tau_ps == pytest.approx(tau, rel=0.01), (polynomial, tau, beta, res.tau_ps)
        assert res.method == "DDCM"
        assert len(res.per_distance) == len(d)


def _singles(rng, n=30000):
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    u = np.exp(-t / TAU_RU)
    return rng.poisson(n * u).astype(float), rng.poisson(n * (1 - u)).astype(float)


@pytest.mark.parametrize("polynomial", ["piecewise", "smooth"])
def test_ddcm_singles_closure(polynomial):
    rng = np.random.default_rng(31)
    for _ in range(3):
        nu, ns = _singles(rng)
        r = la.ddcm(D_RU, nu, ns, (sqrt_counts(nu), sqrt_counts(ns)), BETA_RU, polynomial=polynomial)
        assert within(r, TAU_RU), (r.tau_ps, r.stat_ps)
        assert 0.003 < r.stat_ps / TAU_RU < 0.03


def test_ddcm_with_feeding_correction_and_gated():
    rng = np.random.default_rng(32)
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    tau_b = 3.0
    share = lt.unshifted([tau_b, TAU_RU], [[0, 1], [0, 0]], [0.4, 0.6], t)
    n = 30000
    nu, ns = rng.poisson(n * share[1]).astype(float), rng.poisson(n * (1 - share[1])).astype(float)
    bu = rng.poisson(0.4 * n * share[0]).astype(float)
    bs = rng.poisson(0.4 * n * (1 - share[0])).astype(float)
    errs = (sqrt_counts(nu), sqrt_counts(ns), sqrt_counts(bu), sqrt_counts(bs))
    corrected = la.ddcm(D_RU, nu, ns, errs, BETA_RU, bu, bs, alpha=1.0, alpha_error=0.05)
    assert within(corrected, TAU_RU), (corrected.tau_ps, corrected.stat_ps)
    assert corrected.budget["feeding"] > 0
    # Without the correction the feeding makes the level look slower.
    uncorrected = la.ddcm(D_RU, nu, ns, errs[:2], BETA_RU)
    assert uncorrected.tau_ps > TAU_RU + 5 * uncorrected.stat_ps
    # Gated on the feeder's shifted component: B decayed before t, A after (or before) it.
    lb, la_ = 1 / tau_b, 1 / TAU_RU
    after = lb / (lb - la_) * (np.exp(-la_ * t) - np.exp(-lb * t))
    before = -np.expm1(-lb * t) - after
    m = 40000
    gu, gs = rng.poisson(m * after).astype(float), rng.poisson(m * before).astype(float)
    for polynomial in ("piecewise", "smooth"):
        gated = la.ddcm(D_RU, gu, gs, (sqrt_counts(gu), sqrt_counts(gs)), BETA_RU, gated=True,
                        norm=np.full(len(D_RU), float(m)), polynomial=polynomial)
        assert gated.method == "DDCM gated"
        assert within(gated, TAU_RU), (polynomial, gated.tau_ps, gated.stat_ps)
    with pytest.raises(ValueError):
        la.ddcm(D_RU, gu, gs, errs, BETA_RU, bu, bs, gated=True)


def test_ddcm_shifted_only_mode_and_region():
    rng = np.random.default_rng(33)
    nu, ns = _singles(rng, 100000)
    r = la.ddcm(D_RU, nu, ns, (sqrt_counts(nu), sqrt_counts(ns)), BETA_RU, polynomial="smooth", fit="shifted")
    assert within(r, TAU_RU, 2.5), (r.tau_ps, r.stat_ps)
    region = la.ddcm(D_RU, nu, ns, (sqrt_counts(nu), sqrt_counts(ns)), BETA_RU, region_um=(10, 250))
    assert within(region, TAU_RU, 2.5)
    assert [p["used"] for p in region.per_distance].count(True) <= 7


def test_three_foil_ddcm_cr50():
    """⁵⁰Cr-like (Beckers et al. 2022): β 3.9% after the target, 1.5% after the degrader; pairs of distances."""
    rng = np.random.default_rng(41)
    b1, b2, tau = 0.039, 0.015, 12.2
    x1 = np.array([17, 30, 42, 62, 90, 122, 200, 350, 750.0])
    x2 = np.array([40, 45, 40, 67, 80, 127, 150, 250, 335.0])
    t1 = lt.flight_times(x1, b1, proper=True)
    t2 = t1 + lt.flight_times(x2, b2, proper=True)
    u1, u2 = np.exp(-t1 / tau), np.exp(-t2 / tau)
    shares = (1 - u1, u1 - u2, u2)
    pairs = np.column_stack([x1, x2])
    exact = la.ddcm_three_foil(pairs, *[30000 * x for x in shares], tuple(sqrt_counts(30000 * x) for x in shares),
                               b1, b2)
    assert exact.tau_ps == pytest.approx(tau, rel=0.005)
    assert exact.method == "DDCM three-foil"
    assert len(exact.per_distance) == 2 * len(x1)
    assert exact.per_distance[len(x1)]["distance_um"] == pytest.approx(x1[0] + x2[0] * b1 / b2)
    for _ in range(3):
        n = [rng.poisson(30000 * x).astype(float) for x in shares]
        r = la.ddcm_three_foil(pairs, *n, tuple(sqrt_counts(x) for x in n), b1, b2, beta_reduced_error=0.001)
        assert within(r, tau), (r.tau_ps, r.stat_ps)
        assert r.budget["velocity"] > 0
    # Ignoring the slower second gap (x_eff = x₁ + x₂) puts the stopper points on a different curve.
    wrong = la.ddcm_three_foil(pairs, *[30000 * x for x in shares], tuple(sqrt_counts(30000 * x) for x in shares),
                               b1, b1 * 0.999)
    assert abs(wrong.tau_ps - tau) > 10 * exact.stat_ps


def test_three_foil_with_a_feeder():
    b1, b2, tau, tau_b = 0.039, 0.015, 12.2, 2.4
    x1 = np.array([17, 30, 42, 62, 90, 122, 200, 350, 750.0])
    x2 = np.array([40, 45, 40, 67, 80, 127, 150, 250, 335.0])
    t1 = lt.flight_times(x1, b1, proper=True)
    t2 = t1 + lt.flight_times(x2, b2, proper=True)
    f, pops = [[0, 1], [0, 0]], [0.5, 0.5]
    a1, a2 = lt.unshifted([tau_b, tau], f, pops, t1), lt.unshifted([tau_b, tau], f, pops, t2)
    n, nb = 40000.0, 20000.0  # the feeder carries half of the level's population
    level = [n * (1 - a1[1]), n * (a1[1] - a2[1]), n * a2[1]]
    feeder = [nb * (1 - a1[0]), nb * (a1[0] - a2[0]), nb * a2[0]]
    errs = tuple(sqrt_counts(x) for x in level + feeder)
    r = la.ddcm_three_foil(np.column_stack([x1, x2]), *level, errs, b1, b2, *feeder, alpha=1.0)
    assert r.tau_ps == pytest.approx(tau, rel=0.01)


# ---------------------------------------------------------------------------------------------------------------
# DSAM line shapes


def _fold(row, fine, edges, fwhm):
    c = 0.5 * (fine[1:] + fine[:-1])
    s = fwhm / 2.3548
    return (la._phi((edges[1:, None] - c[None]) / s) - la._phi((edges[:-1, None] - c[None]) / s)) @ row


def _synthetic_line(kernel, taus, feeding, pops, level, rng, e0=1000.0, theta=30.0, fwhm=2.5, area=20000.0):
    edges = np.arange(e0 - 25.0, e0 + 50.0, 1.0)
    fine = np.arange(e0 - 45.0, e0 + 70.0, 0.1)
    row = lt.lineshape(e0, taus, feeding, pops, kernel, theta, fine)[level]
    centres = 0.5 * (edges[1:] + edges[:-1])
    mu = area * _fold(row, fine, edges, fwhm) + 30.0 + 0.05 * (centres - e0)
    return rng.poisson(mu).astype(float), edges


@pytest.fixture(scope="module")
def spread_kernel():
    """Stopping ∝ speed (α = 0.5 ps) with a 6% spread of initial speeds, recoils along the beam."""
    rng = np.random.default_rng(3)
    times = np.linspace(0.0, 4.0, 401)
    b0 = rng.normal(0.035, 0.002, 300)
    return lt.SampledKernel(times, b0[None, :] * np.exp(-times[:, None] / 0.5))


def test_lineshape_closure_one_level(spread_kernel):
    rng = np.random.default_rng(51)
    for tau in (0.15, 0.6):
        counts, edges = _synthetic_line(spread_kernel, [tau], None, [1.0], 0, rng)
        r = la.fit_lineshape(counts, edges, 1000.0, spread_kernel, 30.0, fwhm_kev=2.5)
        assert within(r, tau), (tau, r.tau_ps, r.stat_ps)
        assert r.method == "line shape"
        # The stopping budget: a uniform scaling of the time axis scales τ by the same share.
        assert r.budget["stopping"] == pytest.approx(la.STOPPING_UNCERTAINTY * r.tau_ps, rel=0.15)
        assert r.budget["distance"] == 0.0


def test_lineshape_with_side_feeding_fixed_and_free(spread_kernel):
    rng = np.random.default_rng(52)
    tau, tau_sf = 0.3, 0.1
    feed = la.Feeding.side(tau_sf, 0.4, side_error_ps=0.03)
    counts, edges = _synthetic_line(spread_kernel, [tau_sf, tau], feed.feeding, feed.populations, 1, rng,
                                    area=60000.0)
    fixed = la.fit_lineshape(counts, edges, 1000.0, spread_kernel, 30.0, feed, fwhm_kev=2.5)
    assert within(fixed, tau), (fixed.tau_ps, fixed.stat_ps)
    assert fixed.budget["feeding"] > 0
    free = la.fit_lineshape(counts, edges, 1000.0, spread_kernel, 30.0, la.Feeding.side(0.3, 0.4, free=True),
                            fwhm_kev=2.5, stopping_uncertainty=0.0)
    assert within(free, tau, 2.5), (free.tau_ps, free.stat_ps)
    assert "tau_0_ps" in free.params
    assert free.stat_ps >= fixed.stat_ps * 0.95


def test_lineshape_with_a_transport_kernel_and_fake_stopping(monkeypatch):
    """transport.kernel through a backing whose stopping is ∝ speed (the fake stopping of the transport tests):
    the line shape fit gives τ back, and its stopping budget is the documented scaling."""
    sm = "143Sm"
    mass = data.nuclide(sm).nuclear_mass_mev
    au = data.material("Au")
    alpha = 0.4
    rho = float(au.density_g_cm3)
    k = 2 / (alpha * 1e3 * rho * 2.99792458e-2 * math.sqrt(2 / mass))

    class VelocityStopping:
        def __init__(self, k, mass_u):
            self._e = Stopping.GRID * mass_u
            self.k = k
            self._range = 2 * np.sqrt(self._e) / k
            self._se = k * np.sqrt(self._e) * 1e3
            self._sn = np.zeros_like(self._e)

        def range(self, e):
            e = np.asarray(e, dtype=float)
            return np.where(e > 0, 2 * np.sqrt(np.maximum(e, 0)) / self.k, 0.0)

    monkeypatch.setattr(tr, "stopping", lambda ion, material: VelocityStopping(k, mass / (data.U_KEV * 1e-3)))
    monkeypatch.setattr(tr, "_CLOCKS", {})
    beta0 = 0.02
    energy = mass * (1 / math.sqrt(1 - beta0 ** 2) - 1)
    rng = np.random.default_rng(61)
    energies = energy * rng.normal(1.0, 0.04, 200)
    times = np.linspace(0.0, 10 * alpha, 401)
    kern = tr.kernel([tr.Segment(au, 50.0, name="backing")], sm, mass, energies, np.array([0.0, 0.0, 1.0]), 0.0,
                     times, scatter=False, n=200, rng=1)
    tau = 0.5
    counts, edges = _synthetic_line(kern, [tau], None, [1.0], 0, rng, e0=800.0, theta=35.0, area=30000.0)
    r = la.fit_lineshape(counts, edges, 800.0, kern, 35.0, fwhm_kev=2.0)
    assert within(r, tau), (r.tau_ps, r.stat_ps)
    assert r.budget["stopping"] == pytest.approx(la.STOPPING_UNCERTAINTY * r.tau_ps, rel=0.2)


def _ahkin_kernel():
    """AhKin's kernel file: the density of β_k (β·cos to the detector) at each time, turned into sample recoils by
    quantiles (the line shape needs only each time's distribution). AhKin's Doppler shift is first order,
    E = E₀(1 + β_k); the recoils' β is set so that the relativistic E₀√(1 − β²)/(1 − β) equals it."""
    lines = (AHKIN / "betak_runge.dat").read_text().splitlines()
    head = [line for line in lines if line.startswith("#")]
    dbk = float(head[1].split("dbk=")[1].split()[0])
    dt = float(head[3].split("dt=")[1].split()[0])
    table = np.array([[float(x) for x in line.split()] for line in lines if line.strip() and not line.startswith("#")])
    bk, density = table[:, 0], table[:, 1:]
    nt = density.shape[1]
    edges = np.concatenate([bk - dbk / 2, [bk[-1] + dbk / 2]])
    q = (np.arange(2000) + 0.5) / 2000
    traj = np.zeros((nt, len(q)))
    for j in range(nt):
        p = np.maximum(density[:, j], 0.0)
        if p.sum() > 0:
            traj[j] = np.interp(q, np.concatenate([[0.0], np.cumsum(p)]) / p.sum(), edges)
    beta = ((1 + traj) ** 2 - 1) / ((1 + traj) ** 2 + 1)
    return lt.SampledKernel(dt * np.arange(nt), beta)


def _ahkin_cascade() -> la.Feeding:
    lines = (AHKIN / "in_casc.dat").read_text().splitlines()
    n = int(lines[0])
    taus = [float(x) for x in lines[1].split()]
    b = [[float(x) for x in lines[2 + i].split()] for i in range(n)]
    feeding = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            feeding[j, i] = b[i][j]  # row i of AhKin's matrix: the levels that feed level i
    pops = [float(x) for x in lines[2 + n].split()]
    return la.Feeding(taus, feeding, pops, level=n - 1)


def test_fit_lineshape_recovers_ahkin_on_its_example():
    ref = __import__("json").loads((AHKIN / "result.json").read_text())
    inp = [line.split("#")[0].split() for line in (AHKIN / "in_AhKin.dat").read_text().splitlines()]
    inp = [x for x in inp if x]
    e_g = float(inp[0][0])
    cal0, cal1 = float(inp[1][0]), float(inp[1][1])
    r0, r1, r2 = (float(x) for x in inp[3])
    bkg = float(inp[5][0])
    spec = np.array([[float(x) for x in line.split()] for line in (AHKIN / "EL1192GFAA53.txt").read_text().splitlines()
                     if line.strip()])
    ch, counts = spec[:, 0], spec[:, 1]
    # AhKin's channel c holds energies from cal0 + c·cal1 to cal0 + (c + 1)·cal1 (its peak channel is
    # floor((E_γ − cal0)/cal1)); with these edges our line at its τ matches its best line shape to < 1%.
    edges = cal0 + cal1 * np.concatenate([ch, [ch[-1] + 1.0]])
    fwhm = r0 + r1 * e_g + r2 * e_g ** 2  # AhKin evaluates its resolution at E_γ
    kern = _ahkin_kernel()
    # As AhKin: the background fixed, χ² with the counts as variances, the cascade's feeders fixed.
    r = la.fit_lineshape(counts, edges, e_g, kern, 0.0, _ahkin_cascade(), bkg, fwhm_kev=fwhm, statistic="neyman",
                         stopping_uncertainty=0.0)
    assert ref["tau_ps"] - ref["minus_ps"] <= r.tau_ps <= ref["tau_ps"] + ref["plus_ps"], r.tau_ps
    # AhKin's interval is where χ²/n rises by 1 (Δχ² = n, here 56 channels); ours is Δχ² = 1, so about √56 smaller.
    assert r.stat_ps < 0.5 * ref["plus_ps"]
    # At AhKin's τ, our line (lifetimes.lineshape on AhKin's kernel and cascade, folded, scaled to the data) is
    # AhKin's best line shape (ech0.131308.dat, background included) to 1.5% of the peak.
    ech_lines = (AHKIN / ref["best_lineshape"]).read_text().splitlines()
    ech = np.array([[float(x) for x in line.split()] for line in ech_lines if line.strip()])
    feed = _ahkin_cascade()
    taus = list(feed.taus_ps)
    taus[feed.level] = ref["tau_ps"]
    fine, fold = la._response(edges, e_g, fwhm, 0.95 * e_g, 1.05 * e_g)
    shape = fold @ lt.lineshape(e_g, taus, feed.feeding, feed.populations, kern, 0.0, fine)[feed.level]
    area = ((counts - bkg) @ shape) / (shape @ shape)  # AhKin's normalisation: least squares, unweighted
    np.testing.assert_allclose(area * shape + bkg, ech[:, 1], atol=0.015 * ech[:, 1].max())


# ---------------------------------------------------------------------------------------------------------------
# B from τ


def _scheme():
    levels = [Level(Value(0.0), 0.0, 1, "0+"), Level(Value(652.4), 2.0, 1, "2+"), Level(Value(1397.6), 4.0, 1, "4+"),
              Level(Value(1414.0), 2.0, 1, "2+")]
    transitions = [Transition(1, 0, Value(652.4), Value(100.0), "E2"),
                   Transition(2, 1, Value(745.2), Value(100.0), "E2"),
                   Transition(3, 1, Value(761.6), Value(80.0, 4.0), "M1+E2", Value(-2.0)),
                   Transition(3, 0, Value(1414.0), Value(20.0, 2.0), "E2")]
    return LevelScheme("98Ru", levels, transitions, reference="typed for the test")


def test_b_from_result_inverts_tau_from_b():
    scheme = _scheme()
    b_true = 1300.0
    tau = lt.tau_from_b(652.4, "E2", b_true)
    res = la.TauResult(tau, 0.02 * tau, la._budget(0.02 * tau), "decay curve")
    b, unc = la.b_from_result(res, scheme, "2+₁")
    assert b == pytest.approx(b_true, rel=1e-12)
    assert unc == pytest.approx(0.02 * b_true, rel=1e-9)
    assert la.b_from_result(res, scheme, 1)[0] == pytest.approx(b_true)
    # A level with two branches: the branching comes from the intensities, the mixing from the scheme.
    tau3 = 1.7
    res3 = la.TauResult(tau3, 0.2, la._budget(0.2), "decay curve")
    b_m = la.b_from_result(res3, scheme, "2+₂")  # the strongest branch, 761.6 keV M1+E2: B(E2) with δ = −2
    assert b_m[0] == pytest.approx(lt.b_from_tau(761.6, "E2", tau3, 0.8, 0.0, -2.0), rel=1e-12)
    b_g = la.b_from_result(res3, scheme, "2+₂", final="0+")
    assert b_g[0] == pytest.approx(lt.b_from_tau(1414.0, "E2", tau3, 0.2), rel=1e-12)
    assert b_g[1] / b_g[0] > 0.2 / tau3  # the branching's uncertainty adds to τ's
    with pytest.raises(ValueError):
        la.b_from_result(res3, scheme, "6+")


# ---------------------------------------------------------------------------------------------------------------
# Deorientation


def _ru98_table_one() -> dict:
    setups = tomllib.loads((REF / "published" / "setups.toml").read_text(encoding="utf-8"))
    text = setups["ru98_rdds_coulex"]["analysis"]
    out = {}
    for name in ("alpha_2", "alpha_4", "Gamma_2", "Gamma_4"):
        m = re.search(name + r" ([\d.]+)\((\d+)\)", text)
        value, digits = m.group(1), m.group(2)
        decimals = len(value.split(".")[1]) if "." in value else 0
        out[name] = (float(value), int(digits) * 10.0 ** -decimals)
    return out


def test_ru98_table_one_is_read():
    t = _ru98_table_one()
    assert t["alpha_2"] == (0.18, pytest.approx(0.11))
    assert t["Gamma_2"] == (0.142, pytest.approx(0.031))
    assert t["Gamma_4"] == (0.314, pytest.approx(0.110))


def test_deorientation_correction():
    a = {2: 0.5, 4: -0.2}
    # No attenuation, or the same for both components: no correction.
    assert la.deorientation_correction(40.0, a, {2: 1.0, 4: 1.0}, {2: 1.0, 4: 1.0}) == pytest.approx(1.0)
    g = {2: 0.6, 4: 0.4}
    assert la.deorientation_correction(40.0, a, g, g) == pytest.approx(1.0)
    # By hand at 0°: P_k(1) = 1.
    f = la.deorientation_correction(0.0, a, {2: 0.6, 4: 0.4}, {2: 0.9, 4: 0.8})
    assert f == pytest.approx((1 + 0.5 * 0.6 - 0.2 * 0.4) / (1 + 0.5 * 0.9 - 0.2 * 0.8))
    # Arrays per distance broadcast, with the attenuations of lifetimes.py.
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    al, ga = {2: 0.18, 4: 0.24}, {2: 0.142, 4: 0.314}
    f = la.deorientation_correction(36.5, a, lt.deorientation(t, al, ga), lt.deorientation_shifted(t, TAU_RU, al, ga))
    assert f.shape == t.shape
    assert f[0] == pytest.approx(1.0, abs=0.01) and np.all(f < 1.0 + 1e-12)


def _angular_data(rng, a_k, alpha, gamma, rings, counts=4e5):
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    g_u = lt.deorientation(t, alpha, gamma)
    g_s = lt.deorientation_shifted(t, TAU_RU, alpha, gamma)
    u = np.exp(-t / TAU_RU)
    x = np.cos(np.radians(rings))
    p = {2: 0.5 * (3 * x ** 2 - 1), 4: (35 * x ** 4 - 30 * x ** 2 + 3) / 8}
    out = []
    for share, g in ((u, g_u), (1 - u, g_s)):
        w = 1 + sum(a_k[k] * np.asarray(g[k])[:, None] * p[k][None, :] for k in a_k)
        mu = counts / len(rings) * share[:, None] * w
        y = rng.poisson(mu).astype(float)
        out.append((y, sqrt_counts(y)))
    return out


def test_fit_deorientation_recovers_ru98_table_one():
    table = _ru98_table_one()
    alpha = {2: table["alpha_2"][0], 4: table["alpha_4"][0]}
    gamma = {2: table["Gamma_2"][0], 4: table["Gamma_4"][0]}
    a_k = {2: 0.6, 4: -0.4}  # a strongly aligned Coulex 2⁺ (illustrative)
    rings = np.array([36.5, 46.5, 85.0, 95.0, 133.5, 143.5])  # SPEEDY's rings
    rng = np.random.default_rng(71)
    unsh, shif = _angular_data(rng, a_k, alpha, gamma, rings)
    fa, fg, cov = la.fit_deorientation(D_RU, unsh, shif, a_k, TAU_RU, BETA_RU, theta_deg=rings)
    err = np.sqrt(np.diag(cov))  # (α₂, α₄, Γ₂, Γ₄)
    assert abs(fa[2] - alpha[2]) <= 2 * err[0], (fa, err)
    assert abs(fg[2] - gamma[2]) <= 2 * err[2], (fg, err)
    # With this many counts the fit is more precise than Table I.
    assert err[0] < table["alpha_2"][1] and err[2] < table["Gamma_2"][1]
    # The coefficient form: a_k G_k per distance (as an angular-distribution fit per distance gives them).
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    g_u = lt.deorientation(t, alpha, gamma)
    g_s = lt.deorientation_shifted(t, TAU_RU, alpha, gamma)
    coef_u = {k: (a_k[k] * g_u[k] + rng.normal(0, 0.01, len(t)), np.full(len(t), 0.01)) for k in (2, 4)}
    coef_s = {k: (a_k[k] * g_s[k] + rng.normal(0, 0.01, len(t)), np.full(len(t), 0.01)) for k in (2, 4)}
    fa2, fg2, cov2 = la.fit_deorientation(D_RU, coef_u, coef_s, a_k, TAU_RU, BETA_RU)
    err2 = np.sqrt(np.diag(cov2))
    assert abs(fa2[2] - alpha[2]) <= 2.5 * err2[0] and abs(fg2[2] - gamma[2]) <= 2.5 * err2[2], (fa2, fg2, err2)


def test_deorientation_corrected_ddcm():
    """⁹⁸Ru's 2⁺₁ with deorientation: the measured unshifted share at a forward ring is biased; corrected with
    W_U/W_S from the fitted α, Γ, the decay curve gives τ back."""
    alpha, gamma = {2: 0.18, 4: 0.24}, {2: 0.142, 4: 0.314}
    a_k = {2: 0.6, 4: -0.4}
    t = lt.flight_times(D_RU, BETA_RU, proper=True)
    u = np.exp(-t / TAU_RU)
    f = la.deorientation_correction(36.5, a_k, lt.deorientation(t, alpha, gamma),
                                    lt.deorientation_shifted(t, TAU_RU, alpha, gamma))
    measured = u * f / (u * f + (1 - u))
    biased = la.decay_curve(D_RU, measured, np.full_like(u, 0.002), BETA_RU)
    assert abs(biased.tau_ps - TAU_RU) > 0.1
    corrected = measured / (measured + f * (1 - measured))
    good = la.decay_curve(D_RU, corrected, np.full_like(u, 0.002), BETA_RU)
    assert good.tau_ps == pytest.approx(TAU_RU, rel=1e-6)
