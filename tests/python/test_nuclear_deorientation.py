"""Deorientation in vacuum (the two-state model with GOSIA's parameters): its limits, the time-dependent G_k against the
average over the decay, density matrices, arrays against numbers, and GOSIA's own printed coefficients."""

from __future__ import annotations

import math

import numpy as np
import pytest

from physim.nuclear import deorientation as deo
from physim.nuclear import lifetimes, validation
from physim.nuclear.angular import wigner_6j
from physim.nuclear.deorientation import DEFAULTS, Parameters, TwoState, atomic_spin, charge_states, hard_core
from physim.nuclear.gamma_events import _attenuate

#: Other parameters than the defaults (all of them).
OTHER = Parameters(atomic_spin=2.5, width_per_ps=0.05, rate_per_ps=0.02, correlation_ps=1.5, field_tesla=900.0,
                   field_power=0.45, g_factor=0.8)
LEVELS = (("110Cd", 2), ("110Cd", 4), ("166Er", 6), ("26Mg", 2), ("57Fe", 1.5))


def values(g: dict) -> np.ndarray:
    return np.array([g[k] for k in sorted(g)])


# -- the pieces -----------------------------------------------------------------------------------------------------


def test_the_defaults_are_gosias():
    p = DEFAULTS
    assert (p.atomic_spin, p.width_per_ps, p.rate_per_ps, p.correlation_ps, p.field_tesla, p.field_power,
            p.g_factor) == (3.0, 0.02, 0.0345, 3.5, 600.0, 0.6, None)
    assert TwoState("166Er", 2).g_factor == pytest.approx(68 / 166)
    assert TwoState("166Er", 2, Parameters(g_factor=0.3)).g_factor == 0.3
    for wrong in ({"atomic_spin": 0}, {"atomic_spin": 2.2}, {"width_per_ps": -1}, {"rate_per_ps": -1},
                  {"correlation_ps": 0.5}, {"g_factor": math.nan}):
        with pytest.raises(ValueError):
            Parameters(**wrong)


def test_the_atomic_spins_are_gosias_function_ats():
    """Its branches, line by line (odd numbers of electrons by m = N/2 + 1, even ones by N/2)."""
    odd = {0.5: (1, 2, 3, 6, 7, 10, 15, 16, 19, 24, 25, 28, 31, 35, 37, 40, 41, 44),
           1.5: (4, 5, 8, 9, 11, 17, 18, 20, 26, 27, 36, 42, 43, 45), 2.5: (12, 14, 21, 23, 32, 39), 4.5: (13, 22, 38),
           3.5: (29, 30, 48), 7.5: (33,), 6.5: (34,), 5.5: (46, 47)}
    even = {2.0: (4, 8, 17, 26, 28, 30, 32, 42, 45, 48), 1.0: (10, 36), 3.0: (12, 21, 37),
            4.0: (13, 22, 29, 31, 34, 38, 47), 8.0: (33,), 6.0: (46,)}
    for n in range(-1, 100):
        want = 0.0
        if 0 < n <= 96:
            m, table = (n // 2 + 1, odd) if n % 2 else (n // 2, even)
            want = next((j for j, ms in table.items() if m in ms), 0.0)
        assert atomic_spin(n) == want, n


def test_a_static_field_leaves_the_hard_core():
    assert hard_core(2.0, 0.0, 2) == pytest.approx(1.0) and hard_core(2.0, 3.0, 0) == pytest.approx(1.0)
    assert hard_core(2.0, 3.0, 6) == 0.0 == hard_core(0.0, 3.0, 2)          # no such tensors
    for spin, j, k in ((2.0, 0.5, 2), (2.0, 3.0, 4), (4.0, 2.5, 6), (1.5, 1.0, 2)):
        fs = np.arange(abs(spin - j), spin + j + 0.5)
        # With the interference of different F the sum is the start, 1; the hard core is the part without.
        total = sum((2 * f + 1) * (2 * g + 1) / (2 * j + 1) * wigner_6j(f, g, k, spin, spin, j) ** 2
                    for f in fs for g in fs)
        assert total == pytest.approx(1.0) and 0 < hard_core(spin, j, k) < 1


def test_the_charge_states_are_gosias():
    for z, beta in ((48, 0.04), (68, 0.02), (12, 0.1), (92, 0.06)):
        h = 1 / (1 + (z**0.45 * 0.012008 / beta) ** 1.666667)
        centre = z * h**0.6
        width = math.sqrt(centre * (1 - h)) / 2
        want = {q: math.exp(-(((centre - q) / width) ** 2) / 2)
                for q in range(max(int(centre - 3 * width - 0.5), 1), min(int(centre + 3 * width + 0.5), z) + 1)}
        q, share = charge_states(z, beta)
        got = {int(a): s for a, s in zip(q[0], share[0]) if s > 0}
        assert set(got) == set(want) and share.sum() == pytest.approx(1.0)
        for a, s in want.items():
            assert got[a] == pytest.approx(s / sum(want.values()))
    # Slower than any charge by that rule, and nearly bare: one state is left, not none.
    for z, beta, only in ((28, 1e-4, 1), (2, 0.9, 2)):
        q, share = charge_states(z, beta)
        assert share.sum() == pytest.approx(1.0) and q[0][np.argmax(share[0])] == only


# -- limits -----------------------------------------------------------------------------------------------------------


def test_nothing_happens_without_time_moment_field_tensor_or_switch():
    for g in (TwoState("110Cd", 4).g(0.0, 0.05), TwoState("110Cd", 4).average(0.0, 0.05),
              TwoState("110Cd", 4).average(30.0, 0.0),
              TwoState("110Cd", 4, Parameters(g_factor=0.0)).average(30.0, 0.05),
              TwoState("110Cd", 4, Parameters(g_factor=0.0)).g(np.inf, 0.05),
              TwoState("110Cd", 4, Parameters(field_tesla=0.0)).g(30.0, 0.05),
              TwoState("110Cd", 4, Parameters(rate_per_ps=0.0)).average(30.0, 0.05)):      # Λ* = 0: GOSIA's switch
        assert g == {2: 1.0, 4: 1.0, 6: 1.0}
    assert TwoState("110Cd", 0).average(30.0, 0.05) == {} == TwoState("110Cd", 0.5).g(1.0, 0.05)
    assert [TwoState("110Cd", s).ranks for s in (0, 1, 2, 3, 4)] == [(), (2,), (2, 4), (2, 4, 6), (2, 4, 6)]
    assert TwoState("57Fe", 1.5, ranks="all").ranks == (1, 2, 3)
    # The sign of g does not matter, and it starts flat.
    a, b = (TwoState("110Cd", 2, Parameters(g_factor=g)).g(7.0, 0.05) for g in (0.4, -0.4))
    assert a == b
    assert values(TwoState("110Cd", 2).g(1e-4, 0.05)) == pytest.approx(1.0, abs=1e-6)
    m = TwoState("110Cd", 2)
    for call in (m.g, m.average):
        for args in ((-1.0, 0.05), (math.nan, 0.05), (1.0, -0.1), (1.0, 1.0)):
            with pytest.raises(ValueError):
                call(*args)


@pytest.mark.parametrize("parameters", [None, OTHER])
def test_at_long_times_the_hard_core_is_left(parameters):
    """What reaches the static state, Λ* ∫ G_k^fl(t) e^{−Λ*t} dt, times the hard core of the static field."""
    for nuclide, spin in LEVELS:
        m = TwoState(nuclide, spin, parameters)
        t = np.linspace(0.0, 2000.0, 400_001)
        lam, p = (m.rates(0.06)[x][0] for x in ("lam", "p"))
        tc, big = m.parameters.correlation_ps, m.parameters.rate_per_ps
        c = (1 - lam * tc, lam / p + 2 * lam * tc, -(lam / p + lam * tc))
        fluctuating = sum(np.exp(-np.outer(t, lam + n * p)) * c[n] for n in range(3))
        late = m.alpha_static(0.06)[0] * np.trapezoid(big * fluctuating * np.exp(-big * t)[:, None], t, axis=0)
        assert values(m.g(np.inf, 0.06)) == pytest.approx(late, abs=1e-6)
        assert values(m.g(1e6, 0.06)) == pytest.approx(late, abs=1e-6)
        assert values(m.average(np.inf, 0.06)) == pytest.approx(values(m.g(np.inf, 0.06)), abs=1e-12)
        assert values(m.average(1e8, 0.06)) == pytest.approx(late, abs=1e-5)


# -- in time --------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("parameters", [None, OTHER, Parameters(width_per_ps=0.0), Parameters(correlation_ps=1.0),
                                        Parameters(atomic_spin=0.5, correlation_ps=40.0)])
def test_the_average_over_the_decay_is_gosias_closed_form(parameters):
    x, weight = np.polynomial.legendre.leggauss(12)
    for nuclide, spin in LEVELS:
        m = TwoState(nuclide, spin, parameters, ranks="all")
        for tau in (0.05, 1.0, 7.4, 60.0, 4000.0):
            for beta in (0.015, 0.05, 0.1):
                # ∫₀¹ G(−τ ln s) ds, in pieces: the fast start and the long tail both resolved.
                edges = np.r_[0.0, np.geomspace(1e-6, 1.0, 40), np.linspace(1.25, 40.0, 60)]
                lo, hi = edges[:-1, None], edges[1:, None]
                t = (lo + hi) / 2 + (hi - lo) / 2 * x
                w = (hi - lo) / 2 * weight * np.exp(-t)
                g = m.g(tau * t, beta)
                closed = m.average(tau, beta)
                for k in m.ranks:
                    assert np.sum(w * g[k]) == pytest.approx(closed[k], abs=1e-6), (nuclide, spin, tau, beta, k)


def test_it_only_fades_and_higher_ranks_faster():
    t = np.geomspace(0.01, 1e4, 400)
    for nuclide, spin in LEVELS:
        for parameters in (None, OTHER):
            m = TwoState(nuclide, spin, parameters)
            for beta in (0.02, 0.06):
                g, mean = values(m.g(t, beta)), values(m.average(t, beta))
                assert np.all(np.diff(g, axis=1) <= 1e-12) and np.all((g > 0) & (g <= 1))
                assert np.all(np.diff(mean, axis=1) <= 1e-12) and np.all(mean >= values(m.g(np.inf, beta))[:, None])
                assert np.all(np.diff(mean, axis=0) < 0)


def test_arrays_give_what_numbers_give():
    rng = np.random.default_rng(3)
    t, beta = rng.exponential(20.0, 300), rng.uniform(0.0, 0.12, 300)
    t[:3], beta[3:6] = 0.0, 0.0
    for nuclide, spin in LEVELS:
        m = TwoState(nuclide, spin, ranks="all")
        for call in (m.g, m.average):
            many = call(t, beta)
            for i in range(0, 300, 7):
                one = call(float(t[i]), float(beta[i]))
                assert all(isinstance(v, float) and v == pytest.approx(many[k][i], abs=1e-12) for k, v in one.items())
    m = TwoState("110Cd", 2)
    assert m.g(t.reshape(30, 10), 0.05)[2].shape == (30, 10) and m.g(3.0, beta)[4].shape == (300,)
    assert m.g(np.zeros(0), np.zeros(0))[2].shape == (0,)
    # The same shape of answer as the plunger's static model.
    assert set(m.g(t, 0.05)) == set(lifetimes.deorientation(t, {2: 0.2, 4: 0.1}, {2: 0.1, 4: 0.2}))
    # More nuclei than are worked on at once.
    deo_chunk = deo._CHUNK
    try:
        deo._CHUNK = 70
        assert np.array_equal(m.g(t, beta)[4], TwoState("110Cd", 2).g(t, beta)[4])
        parts = m.g(t, beta)[4]
    finally:
        deo._CHUNK = deo_chunk
    assert parts == pytest.approx(m.g(t, beta)[4], abs=1e-14)
    assert deo.model("110Cd", 2) is deo.model("110Cd", 2.0)


# -- density matrices -------------------------------------------------------------------------------------------------


def test_density_matrices_stay_positive():
    """With every rank attenuated: up to spin 2 (to 1e-6); above, no state is negative by more than 0.005, which is the
    model's own doing (each rank has its own rate). The trace stays 1."""
    rng = np.random.default_rng(1)
    t = np.r_[0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, np.inf]
    for spin, least in ((1, -1e-6), (1.5, -1e-6), (2, -1e-6), (4, -0.005), (6, -0.005)):
        d = int(2 * spin + 1)
        x = rng.normal(size=(60, d)) + 1j * rng.normal(size=(60, d))
        x[:d] = np.eye(d)
        x /= np.linalg.norm(x, axis=1)[:, None]
        rho = x[:, :, None] * x[:, None, :].conj()
        for nuclide in ("26Mg", "110Cd", "166Er", "238U"):
            m = TwoState(nuclide, spin, ranks="all")
            for call in (m.g, m.average):
                for beta in (0.02, 0.06, 0.1):
                    g = call(t, beta)
                    for i in range(len(t)):
                        out = _attenuate(rho, {k: np.full(len(rho), v[i]) for k, v in g.items()})
                        assert np.linalg.eigvalsh(out).min() > least, (nuclide, spin, beta, t[i])
                        assert np.trace(out, axis1=1, axis2=2) == pytest.approx(1.0)


# -- GOSIA ------------------------------------------------------------------------------------------------------------


def test_it_is_gosias_in_blind_runs():
    """Blind GOSIA runs (tests/reference/nuclear/coulex/gosia_gk.csv, the app's numbers stored before them): G_2, G_4
    and G_6 of 13 levels of ⁶²Ni, ¹⁹⁴Pt, ¹⁵²Sm, ⁴⁸Ti and ⁹⁸Ru at v/c of 0.02, 0.05 and 0.08, with GOSIA's parameters
    and with others, within 0.002 and half the last printed digit (largest difference 0.0006)."""
    ref = validation.reference_dir()
    path = ref / "coulex" / "gosia_gk.csv" if ref else None
    if path is None or not path.exists():
        pytest.skip("no GOSIA reference")
    judged, beyond = 0, 0
    for r in validation.read_reference(path)[1]:
        par = Parameters(width_per_ps=float(r["gamma_per_ps"]), rate_per_ps=float(r["lambda_per_ps"]),
                         correlation_ps=float(r["tau_c_ps"]), g_factor=float(r["g"]) if r["g"] else None,
                         field_tesla=float(r["k_1e8t"]) * 1e8, field_power=float(r["x"]))
        assert (par == DEFAULTS) == (r["vac"] == "default")
        spin = float(r["spin"])
        tau, beta = float(r["tau_ps"]), float(r["beta"])
        assert tau == pytest.approx(float(r["tau_ps_gosia"]), rel=3e-4) and beta == pytest.approx(
            float(r["beta_gosia"]), abs=1e-6)
        own = TwoState(r["nuclide"], spin, par).average(tau, beta)
        for k, g in TwoState(r["nuclide"], spin, par, ranks=(2, 4, 6)).average(tau, beta).items():
            assert abs(g - float(r[f"g{k}"])) <= 0.002 + 0.00005, (r["nuclide"], r["level"], r["beta"], r["vac"], k)
            # What GOSIA prints for a rank above 2I acts on nothing: the app has no such coefficient.
            assert own.get(k) == (g if k <= 2 * spin else None)
            judged, beyond = judged + (k <= 2 * spin), beyond + (k > 2 * spin)
    assert (judged, beyond) == (174, 60)
