"""Lifetime physics (planning/lifetimes.md, wave 1B): B ↔ τ, Bateman decay curves, RDDS shares, DSAM attenuation and
line shapes, deorientation in vacuum, against wave 0's references in tests/reference/nuclear/lifetimes/."""

import csv
import math
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import lifetimes as lt
from physim.nuclear.gamma import doppler_energy, level_names
from physim.nuclear.levels import Level, LevelScheme, Transition, Value, rate_per_b, weisskopf_unit
from physim.nuclear.orientation import simple_scheme

REF = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "lifetimes"


def rows(name: str) -> list:
    with open(REF / name, newline="", encoding="utf-8") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines))


def simpson(y, x):
    """Composite Simpson's rule on an odd number of equally spaced points."""
    h = x[1] - x[0]
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def close(actual, expected, rel=1e-9):
    """Relative agreement; values that underflow to 0 in the reference must be (nearly) 0 here too."""
    return actual == pytest.approx(expected, rel=rel, abs=1e-300)


# ---------------------------------------------------------------------------------------------------------------
# RDDS: one level and Bateman


def test_one_level_unshifted_share_matches_the_closed_form():
    table = rows("analytic/rdds_one_level.csv")
    assert len(table) > 100
    for r in table:
        tau, beta, d = float(r["tau_ps"]), float(r["beta"]), float(r["distance_um"])
        t_lab = lt.flight_times(d, beta)
        assert close(t_lab, float(r["flight_time_ps"]))
        # Proper flight time (time dilation included) and without it.
        assert close(lt.unshifted([tau], None, [1.0], lt.flight_times(d, beta, proper=True))[0], float(r["unshifted"]))
        assert close(lt.unshifted([tau], None, [1.0], t_lab)[0], float(r["unshifted_no_dilation"]))


def test_one_level_arrays_and_shapes():
    d = np.array([1.0, 10.0, 100.0])
    t = lt.flight_times(d, 0.05)
    assert lt.unshifted([10.0], None, [1.0], t).shape == (1, 3)
    assert lt.bateman([10.0], None, [2.0], t).shape == (1, 3)
    assert lt.bateman([10.0], None, [2.0], 5.0) == pytest.approx([2.0 * math.exp(-0.5)], rel=1e-13)
    assert lt.flight_times(d, 0.05, proper=True) == pytest.approx(t / lt.lorentz_gamma(0.05), rel=1e-14)
    # A recoil at 30° to the axis flies further to cross the same gap.
    assert lt.flight_times(100.0, 0.05, math.cos(math.radians(30))) == pytest.approx(
        lt.flight_times(100.0, 0.05) / math.cos(math.radians(30)))


def test_two_level_cascade_matches_bateman():
    feeding = [[0.0, 1.0], [0.0, 0.0]]
    table = rows("analytic/rdds_two_levels.csv")
    assert {(r["tau1_ps"], r["tau2_ps"]) for r in table} >= {("10", "10.5")}  # the near-degenerate pair is there
    for r in table:
        taus = [float(r["tau1_ps"]), float(r["tau2_ps"])]
        pops = [float(r["p1"]), float(r["p2"])]
        t = lt.flight_times(float(r["distance_um"]), float(r["beta"]))  # the file's time column is rounded
        assert close(t, float(r["flight_time_ps"]))
        upper, lower = lt.unshifted(taus, feeding, pops, t)
        assert close(upper, float(r["unshifted_upper"]))
        assert close(lower, float(r["unshifted_lower"]))
        n1, _ = lt.bateman(taus, feeding, pops, t)
        assert close(n1 / pops[0], float(r["unshifted_upper"]))


@pytest.mark.parametrize("delta", [1e-3, 1e-6, 1e-9, 1e-13, 0.0])
def test_bateman_is_stable_when_the_lifetimes_are_nearly_equal(delta):
    """N₂ = p₁ τ₂/(τ₁ − τ₂)(e^{−t/τ₁} − e^{−t/τ₂}) loses every digit as τ₂ → τ₁; written with expm1 it does not."""
    tau1, tau2 = 10.0, 10.0 * (1 + delta)
    t = np.array([0.0, 0.1, 1.0, 10.0, 50.0, 300.0, 2000.0])
    d = 1 / tau2 - 1 / tau1
    with np.errstate(invalid="ignore", divide="ignore"):
        stable = np.where(d != 0, -np.expm1(-t * d) / np.where(d != 0, d * tau1, 1.0), t / tau1) * np.exp(-t / tau1)
    n = lt.bateman([tau1, tau2], [[0, 1], [0, 0]], [1.0, 0.0], t)
    assert n[0] == pytest.approx(np.exp(-t / tau1), rel=1e-12)
    assert n[1] == pytest.approx(stable, rel=1e-11, abs=1e-300)


def test_cascade_with_branching_conserves_nuclei():
    taus = [3.0, 0.5, 20.0, np.inf]  # the last level does not decay (the ground state)
    feeding = [[0, 0.7, 0.2, 0.1],
               [0, 0, 0.6, 0.4],
               [0, 0, 0, 1.0],
               [0, 0, 0, 0]]
    t = np.geomspace(1e-3, 500, 40)
    n = lt.bateman(taus, feeding, [1.0, 0.3, 0.0, 0.0], t)
    assert n.sum(axis=0) == pytest.approx(1.3, rel=1e-12)
    assert np.all(n >= 0)
    u = lt.unshifted(taus, feeding, [1.0, 0.3, 0.0, 0.0], t)
    assert np.all(np.diff(u[:3], axis=1) <= 1e-15)  # fewer decays after a longer flight
    assert u[3] == pytest.approx(1.0)  # the ground state never decays: all of it "after"
    # The level-3 share: decays of level 3 after t over all its decays, by direct integration of λ₃ N₃.
    def decays(start):
        x = np.linspace(start, 1000.0, 40_001)
        return simpson(lt.bateman(taus, feeding, [1.0, 0.3, 0.0, 0.0], x)[2] / 20.0, x)

    total = decays(0.0)
    for k in (5, 20, 35):
        assert u[2, k] == pytest.approx(decays(t[k]) / total, rel=1e-6)


def test_feeding_matrix_is_checked():
    with pytest.raises(ValueError):
        lt.bateman([1.0, 2.0], [[0, 1.2], [0, 0]], [1, 0], 1.0)
    with pytest.raises(ValueError):
        lt.unshifted([1.0, 2.0], [[0, 1], [1, 0]], [1, 0], 1.0)  # a closed loop: no decay ever leaves
    with pytest.raises(ValueError):
        lt.bateman([1.0, 0.0], None, [1, 0], 1.0)
    with pytest.raises(ValueError):
        lt.unshifted([1.0], None, [1.0], -1.0)


# ---------------------------------------------------------------------------------------------------------------
# B ↔ τ


def test_rates_per_unit_b_are_the_textbook_ones():
    # T(E2) = 1.223e9 E⁵ B(E2) [e² fm⁴], T(M1) = 1.758e13 E³ B(M1) [µ_N²], E in MeV (Bohr & Mottelson).
    assert lt.rate_per_unit_b("E2", 1000.0) == pytest.approx(1.223e9, rel=3e-3)
    assert lt.rate_per_unit_b("M1", 1000.0) == pytest.approx(1.758e13, rel=3e-3)
    assert lt.rate_per_unit_b("E1", 500.0) == rate_per_b(1, 500.0)


@pytest.mark.parametrize("lam, b, kwargs", [
    ("E2", 1300.0, {}),
    ("E2", 50.0, {"branching": 0.3, "conversion": 0.02}),
    ("E2", 20.0, {"branching": 0.6, "conversion": 0.004, "mixing": -2.5}),
    ("M1", 0.05, {"branching": 0.6, "conversion": 0.004, "mixing": -2.5}),
    ("E1", 1e-4, {}),
    ("E3", 2e4, {"conversion": 0.1}),
])
def test_tau_and_b_are_inverse(lam, b, kwargs):
    tau = lt.tau_from_b(800.0, lam, b, **kwargs)
    assert tau > 0
    assert lt.b_from_tau(800.0, lam, tau, **kwargs) == pytest.approx(b, rel=1e-13)


def test_tau_from_b_for_a_known_state():
    # ⁹⁸Ru 2⁺₁ (652 keV): τ = 8.36 ps (Radeck et al.) is B(E2↓) ≈ 830 e² fm⁴, about 31 W.u.
    b = lt.b_from_tau(652.0, "E2", 8.36, conversion=0.0016)
    assert b == pytest.approx(1.0 / (8.36e-12 * 1.0016) / (1.223e9 * 0.652**5), rel=3e-3)
    assert 28 < b / weisskopf_unit(2, 98) < 34
    # Mixing: the E2 part of an M1+E2 γ ray carries δ²/(1 + δ²) of it.
    pure = lt.tau_from_b(800.0, "E2", 30.0)
    assert lt.tau_from_b(800.0, "E2", 30.0, mixing=1.0) == pytest.approx(pure / 2, rel=1e-13)
    with pytest.raises(ValueError):
        lt.tau_from_b(800.0, "E1", 1e-3, mixing=0.2)
    with pytest.raises(ValueError):
        lt.tau_from_b(800.0, "E2", 30.0, mixing=0.0)


def ensdf_like() -> LevelScheme:
    """A ⁹⁸Ru-like scheme typed as ENSDF gives it: half-lives, γ-ray intensities, conversion and an E2/M1 mixing."""
    def v(x, unc=None):
        return Value(x, unc, "ensdf")

    s = LevelScheme("98Ru", reference="typed for the test")
    s.levels = [Level(v(0.0), 0, 1, "0+"),
                Level(v(652.4), 2, 1, "2+", v(8.36e-12 * math.log(2), 0.2e-12)),
                Level(v(1397.8), 4, 1, "4+", v(2.31e-12 * math.log(2))),
                Level(v(1414.0), 2, 1, "2+", v(1.7e-12 * math.log(2)))]
    s.transitions = [Transition(1, 0, v(652.4), v(100.0), "E2", conversion=v(0.0016)),
                     Transition(2, 1, v(745.4), v(100.0), "E2", conversion=v(0.0011)),
                     Transition(3, 1, v(761.6), v(60.0), "M1+E2", mixing=v(-2.1), conversion=v(0.0013)),
                     Transition(3, 0, v(1414.0), v(40.0), "E2", conversion=v(0.0003))]
    s.derive()
    return s


def test_b_from_tau_inverts_level_scheme_derive():
    s = ensdf_like()
    taus = lt.level_taus(s)
    assert set(taus) == {"2+₁", "4+", "2+₂"}
    assert taus["2+₁"] == (pytest.approx(8.36, rel=1e-12), "half-life (ensdf)")
    branch_total = 60.0 * 1.0013 + 40.0 * 1.0003
    cases = [(1, 0, "2+₁", 1.0, 0.0016, None), (2, 1, "4+", 1.0, 0.0011, None),
             (3, 1, "2+₂", 60.0 * 1.0013 / branch_total, 0.0013, -2.1),
             (3, 0, "2+₂", 40.0 * 1.0003 / branch_total, 0.0003, None)]
    for i, f, name, branching, alpha, mixing in cases:
        e = s.levels[i].energy.value - s.levels[f].energy.value
        b = lt.b_from_tau(e, "E2", taus[name][0], branching, alpha, mixing)
        assert b == pytest.approx(s.b(i, f, "E2"), rel=1e-12)


def test_level_taus_from_b_when_there_is_no_half_life():
    s = ensdf_like()
    expected = {name: tau for name, (tau, _) in lt.level_taus(s).items()}
    for level in s.levels:
        level.half_life = None  # the derived matrix elements stay
    taus = lt.level_taus(s)
    for name, tau in expected.items():
        assert taus[name][0] == pytest.approx(tau, rel=1e-12)
        assert taus[name][1] == "B(E2) (derived)"


def test_level_taus_from_a_typed_b_e2():
    s = simple_scheme("58Ni", "1.4542 MeV", "E2", "700 e2fm4")
    name = level_names(s.levels)[1]  # "2": a typed level without a parity written
    tau, source = lt.level_taus(s)[name]
    assert source == "B(E2) (user)"
    assert tau == pytest.approx(lt.tau_from_b(1454.2, "E2", s.b(1, 0, "E2")), rel=1e-12)
    # With the half-life typed too, it wins.
    s.levels[1].half_life = Value(1e-12, source="user")
    assert lt.level_taus(s)[name] == (pytest.approx(1.0 / math.log(2), rel=1e-12), "half-life (user)")


def test_level_taus_with_unknown_branching_is_an_upper_limit():
    s = ensdf_like()
    for level in s.levels:
        level.half_life = None
    s.transitions[2].intensity = None  # the 2⁺₂ → 2⁺₁ intensity is unknown
    tau, source = lt.level_taus(s)["2+₂"]
    assert source.endswith("other branches unknown")
    # Both branches have a B (derived before the intensity was lost), so their rates add up to the true total.
    assert tau == pytest.approx(1.7, rel=1e-12)


# ---------------------------------------------------------------------------------------------------------------
# Planning: sensitive range, distances, peak separation


@pytest.mark.parametrize("case, tau, beta, used", [
    ("98Ru 2+1 (Radeck et al. 2012)", 8.36, 13.49 / lt.C_UM_PER_PS, (3.0, 450.0)),
    ("50Cr 2+1 (Beckers et al. 2022)", 12.2, 0.039, (17.0, 750.0)),
])
def test_sensitive_range_and_distances_agree_with_published_choices(case, tau, beta, used):
    d_min, d_max = lt.sensitive_range(tau, beta)
    v_tau = beta * lt.C_UM_PER_PS * tau
    assert d_min < v_tau < d_max
    # The published series reach the short end and end within a factor 1.5 of the long end of the range.
    assert used[0] <= d_min
    assert used[1] / 1.5 <= d_max <= used[1] * 1.5
    distances = lt.recommend_distances(tau, beta)
    assert len(distances) == 8
    assert distances == sorted(distances)
    assert used[0] <= distances[0] and distances[-1] <= 1.5 * used[1]
    ratios = np.diff(np.log(distances))
    assert ratios == pytest.approx(np.full(7, ratios.mean()), rel=0.15)  # log-spaced (rounded to 2 figures)


def test_sensitive_range_is_where_the_precision_is_within_the_factor():
    d_min, d_max = lt.sensitive_range(10.0, 0.05)
    length = lt.lorentz_gamma(0.05) * 0.05 * lt.C_UM_PER_PS * 10.0
    best = lt._precision(lt._best_x())
    assert best <= min(lt._precision(x) for x in np.linspace(0.5, 3, 2001))
    assert lt._precision(d_min / length) == pytest.approx(2 * best, rel=1e-9)
    assert lt._precision(d_max / length) == pytest.approx(2 * best, rel=1e-9)
    wide = lt.sensitive_range(10.0, 0.05, factor=3.0)
    assert wide[0] < d_min and wide[1] > d_max
    assert lt.recommend_distances(10.0, 0.05, n=1) == [pytest.approx(lt._round_sig(1.5936 * length), rel=1e-3)]


def test_peak_separation_is_the_relativistic_doppler_shift():
    e0, beta = 652.0, 0.045
    for theta in (0.0, 36.5, 90.0, 143.5, 180.0):
        expected = doppler_energy(e0, beta, math.cos(math.radians(theta))) - e0
        assert lt.peak_separation(e0, beta, theta) == pytest.approx(expected, rel=1e-14)
    assert lt.peak_separation(e0, beta, 36.5) > 20  # keV: separates easily from a 2-keV germanium peak
    assert lt.peak_separation(e0, beta, 143.5) < -20
    assert lt.peak_separation(e0, beta, 90.0) == pytest.approx(-e0 * beta**2 / 2, rel=1e-2)  # second order only


# ---------------------------------------------------------------------------------------------------------------
# DSAM


def test_dsam_attenuation_for_stopping_proportional_to_speed():
    table = rows("analytic/dsam_exponential.csv")
    assert len(table) >= 20
    kernels = {}
    for r in table:
        alpha, tau = float(r["alpha_ps"]), float(r["tau_ps"])
        k = kernels.setdefault(alpha, lt.exponential_kernel(0.04, alpha))
        assert isinstance(k, lt.Kernel)
        assert lt.dsam_attenuation(tau, k) == pytest.approx(float(r["attenuation_F"]), rel=1e-6)


def test_dsam_line_shape_gives_the_distribution_of_speeds_at_emission():
    """At 0° each v/v₀ bin maps to one energy bin; the shares must match the reference's binned distribution."""
    beta0, e0 = 0.04, 1000.0
    x = np.linspace(0.0, 1.0, 21)
    edges = doppler_energy(e0, beta0 * x, 1.0)
    for r in rows("analytic/dsam_exponential.csv"):
        alpha, tau = float(r["alpha_ps"]), float(r["tau_ps"])
        expected = np.array([float(v) for v in r["share_v_over_v0_in_20_bins"].split(";")])
        shape = lt.lineshape(e0, [tau], None, [1.0], lt.exponential_kernel(beta0, alpha), 0.0, edges)
        assert shape.shape == (1, 20)
        assert shape[0].sum() == pytest.approx(1.0, abs=1e-9)
        assert shape[0] == pytest.approx(expected, abs=2e-5)


def test_dsam_line_shape_of_a_fed_level_follows_its_decay_curve():
    """A level fed from a longer-lived one: the share of its γ rays with v/v₀ in a bin is the share of its decays
    in the times the recoil spends at those speeds (from :func:`unshifted`)."""
    alpha, beta0, e0 = 0.5, 0.03, 800.0
    taus, feeding, pops = [0.8, 0.3], [[0, 1.0], [0, 0]], [0.7, 0.3]
    x = np.linspace(0.05, 1.0, 20)
    edges = doppler_energy(e0, beta0 * x, 1.0)
    shape = lt.lineshape([1200.0, e0], taus, feeding, pops, lt.exponential_kernel(beta0, alpha), 0.0, edges)
    after = lt.unshifted(taus, feeding, pops, alpha * np.log(1.0 / x))[1]  # x = e^{−t/α}
    assert shape[1] == pytest.approx(np.diff(after), abs=2e-5)
    # The mean speed at emission is the product of the two attenuations (a convolution of exponentials).
    fine = np.linspace(0.0, 1.0, 2001)
    rows_ = lt.lineshape(e0, taus, feeding, [1.0, 0.0], lt.exponential_kernel(beta0, alpha), 0.0,
                         doppler_energy(e0, beta0 * fine, 1.0))
    mean_x = rows_[1] @ (0.5 * (fine[1:] + fine[:-1]))
    assert mean_x == pytest.approx(alpha / (alpha + 0.8) * alpha / (alpha + 0.3), rel=1e-4)


def test_dsam_line_shape_from_a_sampled_kernel_with_spread():
    """Recoils at several angles: each ring sees the projection; the stopped peak holds the late decays."""
    rng = np.random.default_rng(3)
    cos = rng.uniform(0.6, 1.0, 400)
    t = np.linspace(0.0, 4.0, 401)
    kernel = lt.SampledKernel(t, 0.035 * np.clip(1 - t / 1.2, 0, None)[:, None] * cos[None, :])
    edges = np.linspace(990.0, 1040.0, 251)
    forward = lt.lineshape(1000.0, [0.6], None, [1.0], kernel, 30.0, edges)[0]
    backward = lt.lineshape(1000.0, [0.6], None, [1.0], kernel, 150.0, np.linspace(960.0, 1010.0, 251))[0]
    stopped_share = lt.unshifted([0.6], None, [1.0], 1.2)[0]  # decays after all recoils have stopped
    assert forward.sum() == pytest.approx(1.0, abs=1e-9)
    assert backward.sum() == pytest.approx(1.0, abs=1e-9)
    stopped_bin = np.searchsorted(edges, 1000.0, side="right") - 1
    # The stopped bin (0.2 keV wide) also holds the decays of the last moments before stopping.
    assert stopped_share <= forward[stopped_bin] <= stopped_share + 0.005
    centres = 0.5 * (edges[1:] + edges[:-1])
    assert forward @ centres > 1000.0
    assert backward @ (0.5 * (np.linspace(960.0, 1010.0, 251)[1:] + np.linspace(960.0, 1010.0, 251)[:-1])) < 1000.0


# ---------------------------------------------------------------------------------------------------------------
# Deorientation


RADECK = {"alpha": {2: 0.18, 4: 0.24}, "gamma_per_ps": {2: 0.142, 4: 0.314}}  # PRC 85, 014301, Table I


def test_deorientation_closed_form():
    t = np.array([0.0, 1.0, 10.0, 1e4])
    g = lt.deorientation(t, **RADECK)
    assert g[2] == pytest.approx(0.18 + 0.82 * np.exp(-0.142 * t), rel=1e-15)
    assert g[4][0] == 1.0 and g[4][-1] == pytest.approx(0.24)
    assert lt.deorientation(5.0, {2: 0.3}, {2: 0.0})[2] == 1.0


@pytest.mark.parametrize("tau", [0.5, 8.36, 100.0, math.inf])
def test_deorientation_shifted_matches_numeric_integration(tau):
    flights = np.array([1e-9, 0.05, 1.0, 8.0, 33.0, 200.0])
    shifted = lt.deorientation_shifted(flights, tau, **RADECK)
    for k in (2, 4):
        for n, flight in enumerate(flights):
            x = np.linspace(0.0, flight, 20001)
            weight = np.exp(-x / tau) if math.isfinite(tau) else np.ones_like(x)
            gk = lt.deorientation(x, **RADECK)[k]
            numeric = simpson(gk * weight, x) / simpson(weight, x)
            assert shifted[k][n] == pytest.approx(numeric, rel=1e-8)
    # No flight: no deorientation; the shifted share is less attenuated than the unshifted one at the stopper.
    assert lt.deorientation_shifted(0.0, tau, **RADECK)[2] == 1.0
    g_unshifted = lt.deorientation(flights, **RADECK)
    assert np.all(shifted[2][1:] > g_unshifted[2][1:])


def test_dsam_line_shape_binning_is_kept_between_calls(monkeypatch):
    """The kernel binned once (for fits over τ) gives the same shape as binning it with each τ's weights."""
    rng = np.random.default_rng(5)
    t = np.linspace(0.0, 2.0, 81)
    kernel = lt.SampledKernel(t, 0.03 * np.clip(1 - t / 1.0, 0, None)[:, None] * rng.uniform(0.5, 1.0, 300)[None, :],
                              weights=rng.uniform(0.5, 1.5, 300))
    edges = np.linspace(995.0, 1035.0, 161)
    args = ([1000.0, 1000.0], [0.4, 0.2], [[0, 0.8], [0, 0]], [1.0, 0.2], kernel, 20.0, edges)
    kept = lt.lineshape(*args)
    again = lt.lineshape(*args)
    monkeypatch.setattr(lt, "_MATRIX_LIMIT", 0)
    direct = lt.lineshape(*args)
    assert again == pytest.approx(kept, rel=1e-14, abs=1e-18)
    assert kept == pytest.approx(direct, rel=1e-10, abs=1e-15)
