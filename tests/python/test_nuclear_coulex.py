"""Coulomb excitation and γ-ray Doppler shifts (backlog item 43)."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment, SetupError
from physim.nuclear.coulex import HBARC_MEV_FM, Coulex, orbit_integrals, parse_b, ylm_equator
from physim.nuclear.gamma import doppler_energy, doppler_table
from physim.nuclear.rutherford import E2_MEV_FM, Rutherford

trapezoid = getattr(np, "trapezoid", None) or np.trapz  # NumPy 1.x has only trapz


def ni58(**kw):
    args = dict(excite="target", energy="1.454 MeV", multipolarity="E2", b_up="0.0695 e2b2")
    args.update(kw)
    return Coulex("16O", "58Ni", kw.pop("beam_energy", 30.0), **{k: v for k, v in args.items()
                                                                if k != "beam_energy"})


@pytest.mark.parametrize("lam", [1, 2, 3])
def test_spherical_harmonics_on_the_orbit_plane(lam):
    assert sum(ylm_equator(lam, mu) ** 2 for mu in range(-lam, lam + 1)) == pytest.approx(
        (2 * lam + 1) / (4 * math.pi), rel=1e-12)
    for mu in range(-lam, lam + 1):
        if (lam + mu) % 2:
            assert ylm_equator(lam, mu) == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize("lam, exact", [(1, 2.0), (2, 2 / 3), (3, 4 / 15)])
def test_backscattering_in_the_sudden_limit_is_a_closed_form(lam, exact):
    """At θ = 180° (ε = 1) and ξ = 0 the orbit is a straight line in and out: every I_μ is ∫ dw/(1 + cosh w)^λ,
    which is 2, 2/3 and 4/15 for λ = 1, 2, 3."""
    np.testing.assert_allclose(orbit_integrals(lam, 1.0, 0.0), exact, rtol=1e-8)


@pytest.mark.parametrize("lam, xi", [(1, 0.381), (2, 0.381), (2, 0.112), (3, 1.203)])
def test_orbit_integrals_agree_with_a_brute_force_integration(lam, xi):
    """The panels, the half-orbit symmetry and the endpoint expansion of the tail against a plain fine
    trapezoidal integration of the whole orbit, to the stated 1e-9 (absolute; the integrals are of order 1)."""
    from physim.nuclear.coulex import _orbit_pieces, orbit_integrals

    for eps in (1.0, 1.5, 4.0, 20.0):
        w = np.linspace(0.0, 16.0, 2_000_001)
        g, phase, _ = _orbit_pieces(lam, eps, xi, w)
        e = np.exp(1j * phase)
        brute = np.array([2 * np.trapezoid(e * (g[mu] if mu >= 0 else np.conj(g[-mu])), w).real
                          for mu in range(-lam, lam + 1)])
        got = orbit_integrals(lam, eps, xi)
        assert np.max(np.abs(got - brute)) < 5e-9, (eps, got, brute)
        assert np.all(got.imag == 0)


def test_a_probability_table_is_fast():
    import time

    from physim.nuclear.coulex import orbit_integrals

    orbit_integrals.cache_clear()
    c = Coulex("16O", "58Ni", 50.0, energy=1.454, multipolarity="E1", b_up="0.01 e2b")  # the slowest case before
    t = time.perf_counter()
    c.probability(np.linspace(1.0, 180.0, 180))
    assert time.perf_counter() - t < 2.0


def test_excitation_probability_closed_form_at_180_degrees():
    c = ni58(energy="1e-9 MeV")  # ξ → 0
    assert c.xi < 1e-8
    s = 5 / (4 * math.pi)
    pref = (4 * math.pi * 8 * E2_MEV_FM / (HBARC_MEV_FM * 5)) ** 2 * c.b_up / 5
    closed = pref * (2 / 3) ** 2 * s / (c.beta * c.a**2) ** 2
    assert c.probability(180.0, exact=True) == pytest.approx(closed, rel=1e-6)


@pytest.mark.parametrize("theta", [60.0, 120.0, 170.0])
def test_agrees_with_integration_along_engine_orbits(theta):
    """The same first-order amplitude, integrated in time along a Coulomb orbit that physim's engine integrates
    (adaptive Dormand–Prince), with no use of the hyperbolic parametrisation."""
    r = Rutherford("16O", "58Ni", 40.0)
    m1, m2 = r.kinematics.m1, r.kinematics.m2
    mu = m1 * m2 / (m1 + m2)
    v = math.sqrt(2 * r.e_cm / mu)
    a = r.k / (mu * v * v)
    orbit = r.trajectories([r.impact_parameter(theta)], distance=4000 * r.d0, samples=200_001, rtol=1e-12)[0]
    x, y, t = orbit.pos[:, 0], orbit.pos[:, 1], orbit.t
    rr, phi = np.hypot(x, y), np.arctan2(y, x)
    for xi in (0.0, 0.4):
        omega = xi * v / a
        engine = sum(abs(ylm_equator(2, m) * trapezoid(np.exp(1j * (omega * t + m * phi)) / rr**3, t)) ** 2
                     for m in range(-2, 3))
        integrals = orbit_integrals(2, 1 / math.sin(math.radians(theta) / 2), xi)
        mine = sum(abs(ylm_equator(2, m) * integrals[m + 2]) ** 2 for m in range(-2, 3)) / (v * a**2) ** 2
        assert engine == pytest.approx(mine, rel=1e-6), (theta, xi)


def test_adiabatic_suppression():
    """Excitation dies away once the collision is slower than the nuclear period (ξ ≳ 1)."""
    probs, xis = [], []
    for e_star in (0.01, 0.5, 1.0, 2.0, 4.0, 6.0):
        c = ni58(energy=e_star)
        xis.append(c.xi)
        # Same B: compare P / (prefactor at that velocity), i.e. the orbit factor alone.
        probs.append(c.probability(180.0, exact=True) * (c.beta * c.a**2) ** 2)
    assert np.all(np.diff(xis) > 0)
    assert np.all(np.diff(probs) < 0)
    assert probs[-1] / probs[0] < 0.01


def test_scaling_with_b_and_with_the_exciting_charge():
    base = ni58()
    assert ni58(b_up="0.139 e2b2").probability(120.0) == pytest.approx(2 * base.probability(120.0), rel=1e-12)
    assert parse_b("695 e2fm4", 2) == pytest.approx(695.0, rel=1e-12)
    assert parse_b("0.0695 e2b2", 2) == pytest.approx(695.0, rel=1e-12)
    # The same orbit and state, but the projectile excited (by Z = 28) instead of the target (by Z = 8).
    proj = ni58(excite="projectile")
    assert proj.xi == pytest.approx(base.xi) and proj.a == pytest.approx(base.a)
    assert proj.probability(120.0) == pytest.approx(base.probability(120.0) * (28 / 8) ** 2, rel=1e-12)
    with pytest.raises(ValueError, match="unit"):
        parse_b("0.05 e2b", 2)


def test_interpolated_probability_is_accurate():
    c = ni58()
    th = np.array([33.7, 77.25, 120.4, 155.55, 179.6])
    np.testing.assert_allclose(c.probability(th), c.probability(th, exact=True), rtol=1e-5)
    assert c.probability(0.0) == 0.0


def test_cross_sections():
    c = ni58()
    assert c.cross_section_cm(120.0) == pytest.approx(c.probability(120.0) * c.rutherford_cm(120.0))
    lab, second = c.cross_section_lab(100.0)
    assert second is None and lab > 0
    total = c.total()
    assert 0 < total < 100  # mb


def test_safe_distance():
    safe = ni58()  # 30 MeV: safe at every angle
    assert safe.max_safe_angle() == 180.0 and not any("safe distance" in w for w in safe.warnings())
    unsafe = ni58(beam_energy=48.0)
    th = unsafe.max_safe_angle()
    assert 0 < th < 180
    assert unsafe.closest_approach(th) == pytest.approx(unsafe.safe_distance, rel=1e-9)
    assert unsafe.safe_distance == pytest.approx(1.25 * (16 ** (1 / 3) + 58 ** (1 / 3)) + 5)
    assert any("safe distance" in w for w in unsafe.warnings())
    assert not any("safe distance" in w for w in unsafe.warnings(theta_cm_max=th - 1))


def test_setup_file():
    exp = Experiment.example("coulex_ni58")
    assert exp.reaction == "coulex" and exp.excitation.energy_mev == pytest.approx(1.454)
    assert exp.excitation.b_up_e2fm == pytest.approx(695.0)
    assert len(exp.gamma_detectors) == 4
    assert Experiment.from_toml(exp.to_toml()).to_dict() == exp.to_dict()
    d = exp.to_dict()
    del d["reaction"]["b_up"]
    with pytest.raises(SetupError, match="b_up"):
        Experiment.from_dict(d)
    d = exp.to_dict()
    d["reaction"]["b_up"] = "0.05 e2b"
    with pytest.raises(SetupError, match="not a unit of B"):
        Experiment.from_dict(d)
    d = exp.to_dict()
    d["reaction"] = {"type": "elastic", "energy": "1 MeV"}
    with pytest.raises(SetupError, match="only apply"):
        Experiment.from_dict(d)


def test_doppler_formula_is_a_lorentz_boost():
    """E_γ from boosting the photon four-vector (E₀, E₀ n̂') of the rest frame into the lab."""
    rng = np.random.default_rng(1)
    for _ in range(20):
        beta = rng.uniform(0, 0.5)
        cos_rest = rng.uniform(-1, 1)
        g = 1 / math.sqrt(1 - beta**2)
        e_lab = g * (1.0 + beta * cos_rest)  # E₀ = 1
        cos_lab = (cos_rest + beta) / (1 + beta * cos_rest)  # aberration
        assert doppler_energy(1.0, beta, cos_lab) == pytest.approx(e_lab, rel=1e-12)


def test_doppler_table():
    rows = {(r["particle_detector"], r["gamma_detector"]): r for r in doppler_table(Experiment.example("coulex_ni58"))}
    cd45, cd135, cd90 = rows[("CD", "Ge45")], rows[("CD", "Ge135")], rows[("CD", "Ge90")]
    beta = cd45["beta_mean"]
    # Backscattered 16O in the CD (134°–167°): the 58Ni recoils go forward, 6°–23° from the beam, all around it, so
    # the shift is E₀ β ⟨cos α⟩ with α spread about 45°: between E₀ β cos 68° and E₀ β cos 22°, and near cos 45°.
    assert 1454 * beta * math.cos(math.radians(68)) < cd45["shift_kev"] < 1454 * beta * math.cos(math.radians(22))
    assert cd45["shift_kev"] == pytest.approx(1454 * beta * math.cos(math.radians(45)), rel=0.08)
    assert cd135["shift_kev"] == pytest.approx(-cd45["shift_kev"], rel=0.05)
    assert abs(cd90["shift_kev"]) < 0.05 * cd45["shift_kev"]
    for r in rows.values():
        assert r["fwhm_kev"] >= r["doppler_fwhm_kev"] > 0
