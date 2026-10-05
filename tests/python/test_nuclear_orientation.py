"""Orientation of excited states, particle–γ correlation and decay (backlog 54)."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment, angular
from physim.nuclear.coulex import Coulex, orbit_integrals
from physim.nuclear.gamma import correlation_table, excitation_of
from physim.nuclear.levels import Level, LevelScheme, Transition, Value
from physim.nuclear.orientation import Excitation, Populated, fed_density, simple_scheme, tensors_of
from physim.nuclear.planner import Planner

NI = ("16O", "58Ni", "30 MeV")


def ni(**kw) -> Excitation:
    return Excitation(*NI, simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2"), **kw)


def sphere(n: int = 24):
    """Directions and weights that integrate over all directions exactly for the low orders used here."""
    x, w = np.polynomial.legendre.leggauss(n)
    phi = np.linspace(0, 2 * np.pi, 2 * n + 1)[:-1]
    theta, ph = np.meshgrid(np.arccos(x), phi, indexing="ij")
    return theta, ph, np.repeat(w[:, None], len(phi), axis=1) * 2 * np.pi / len(phi)


def cascade_scheme(alpha: float = 0.0, branch_to_ground: float = 0.0, mixing=None) -> LevelScheme:
    """0⁺, 2⁺ (500 keV), 4⁺ (1200 keV, not excited directly) and a second 2⁺ (1500 keV) that feeds the first."""
    s = LevelScheme("58Ni")
    v = lambda x: Value(x, source="user")  # noqa: E731
    s.levels = [Level(v(0.0), 0.0, 1), Level(v(500.0), 2.0, 1), Level(v(1200.0), 4.0, 1), Level(v(1500.0), 2.0, 1)]
    s.transitions = [Transition(1, 0, v(500.0), v(100.0), "E2", conversion=v(alpha) if alpha else None),
                     Transition(3, 1, v(1000.0), v(100.0), "M1+E2" if mixing is not None else "E2",
                                mixing=v(mixing) if mixing is not None else None)]
    if branch_to_ground:
        s.transitions.append(Transition(3, 0, v(1500.0), v(branch_to_ground), "E2"))
    s.set_matrix_element(0, 1, "E2", 30.0)
    s.set_matrix_element(0, 3, "E2", 12.0)
    return s


# -- the algebra ----------------------------------------------------------------------------------------------------


def test_coupling_coefficients_have_their_textbook_values():
    assert angular.clebsch_gordan(2, 0, 2, 0, 2, 0) == pytest.approx(-math.sqrt(2 / 7))
    assert angular.clebsch_gordan(0.5, 0.5, 0.5, -0.5, 0, 0) == pytest.approx(1 / math.sqrt(2))
    assert angular.clebsch_gordan(1, 1, 1, -1, 2, 0) == pytest.approx(1 / math.sqrt(6))
    assert angular.wigner_3j(2, 2, 0, 1, -1, 0) == pytest.approx(-1 / math.sqrt(5))
    assert angular.wigner_6j(1, 1, 1, 1, 1, 1) == pytest.approx(1 / 6)
    assert angular.wigner_6j(2, 2, 2, 2, 2, 2) == pytest.approx(-3 / 70)
    assert angular.wigner_3j(2, 2, 2, 1, 1, 1) == 0 and angular.wigner_6j(1, 1, 3, 1, 1, 1) == 0
    # Orthogonality: Σ_m ⟨j1 m1 j2 m2 | J M⟩² over all (m1, m2) is 1.
    total = sum(angular.clebsch_gordan(1.5, m1, 1, m2, 2.5, 0.5) ** 2
                for m1 in (-1.5, -0.5, 0.5, 1.5) for m2 in (-1, 0, 1))
    assert total == pytest.approx(1.0)


def test_f_and_u_coefficients_match_the_tables():
    # Krane and Steffen's tables.
    assert angular.f_coefficient(2, 2, 2, 0, 2) == pytest.approx(-0.5976, abs=1e-4)
    assert angular.f_coefficient(4, 2, 2, 0, 2) == pytest.approx(-1.0690, abs=1e-4)
    assert angular.f_coefficient(2, 1, 1, 0, 1) == pytest.approx(0.7071, abs=1e-4)
    assert angular.f_coefficient(2, 2, 2, 2, 4) == pytest.approx(-0.4477, abs=1e-4)
    assert angular.f_coefficient(0, 2, 2, 0, 2) == pytest.approx(1.0)
    assert angular.u_coefficient(0, 4, 2, 2) == pytest.approx(1.0)
    assert angular.u_coefficient(2, 4, 2, 2) == pytest.approx(0.7491, abs=1e-4)
    assert angular.u_coefficient(4, 4, 2, 2) == pytest.approx(0.2847, abs=1e-4)
    # A mixed transition: A_k = [F(LL) + 2δF(LL′) + δ²F(L′L′)] / (1 + δ²).
    d = 0.5
    expected = (angular.f_coefficient(2, 1, 1, 2, 2) + 2 * d * angular.f_coefficient(2, 1, 2, 2, 2)
                + d * d * angular.f_coefficient(2, 2, 2, 2, 2)) / (1 + d * d)
    assert angular.distribution_coefficient(2, 2, 2, 1, 2, d) == pytest.approx(expected)


def test_spherical_harmonics_are_orthonormal():
    theta, phi, w = sphere(12)
    for (l1, m1), (l2, m2) in (((2, 1), (2, 1)), ((4, -3), (4, -3)), ((2, 1), (4, 1)), ((3, 2), (3, -2))):
        integral = np.sum(angular.spherical_harmonic(l1, m1, theta, phi)
                          * np.conj(angular.spherical_harmonic(l2, m2, theta, phi)) * w)
        assert integral == pytest.approx(1.0 if (l1, m1) == (l2, m2) else 0.0, abs=1e-12)
    th = np.linspace(0, np.pi, 9)
    assert np.allclose(np.abs(angular.spherical_harmonic(2, 1, th, 0.3)) ** 2,
                       15 / (8 * np.pi) * np.sin(th) ** 2 * np.cos(th) ** 2)


# -- excitation -----------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("multipolarity, b_up", [("E2", "0.0695 e2b2"), ("E3", "0.02 e2b3"), ("E1", "0.01 e2b")])
def test_summed_over_substates_the_probability_is_the_first_order_one(multipolarity, b_up):
    scheme = simple_scheme("58Ni", "1.454 MeV", multipolarity, b_up)
    ex = Excitation(*NI, scheme)
    cx = Coulex(*NI, energy="1.454 MeV", multipolarity=multipolarity, b_up=b_up)
    for theta in (20.0, 60.0, 120.0, 180.0):
        assert ex.probability(1, theta) == pytest.approx(cx.probability(theta, exact=True), rel=1e-8)
    assert ex.cross_sections()["direct"][1] == pytest.approx(cx.total(), rel=5e-3)


def test_any_ground_state_spin_gives_the_same_total_for_the_same_b():
    """B(E2↑) fixes the excitation probability whatever the spins: the average over the ground state's
    substates and the sum over the final ones see to that."""
    even = ni().probability(1, 130.0)
    for spin in (0.5, 1.5, 2.0):
        scheme = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2", ground_spin=spin)
        assert scheme.levels[1].spin == spin + 2
        assert Excitation(*NI, scheme).probability(1, 130.0) == pytest.approx(even, rel=1e-10)


def test_projectile_excitation_uses_the_target_charge():
    scheme = simple_scheme("16O", "6.917 MeV", "E2", "0.004 e2b2")
    ex = Excitation("16O", "58Ni", "60 MeV", scheme, excite="projectile")
    cx = Coulex("16O", "58Ni", "60 MeV", excite="projectile", energy="6.917 MeV", multipolarity="E2",
                b_up="0.004 e2b2")
    assert ex.probability(1, 150.0) == pytest.approx(cx.probability(150.0, exact=True), rel=1e-8)


def test_levels_that_cannot_be_excited_are_left_out_with_a_note():
    s = cascade_scheme()
    s.set_matrix_element(0, 2, "E2", 5.0)     # 0⁺ → 4⁺ by E2: the spins do not allow it
    s.levels.append(Level(Value(40000.0, source="user"), 2.0, 1))
    s.set_matrix_element(0, 4, "E2", 5.0)     # above the energy available
    ex = Excitation(*NI, s)
    assert sorted(ex.paths) == [1, 3]
    assert any("E2 cannot connect" in n for n in ex.notes) and any("above the energy" in n for n in ex.notes)


def test_orbit_integrals_without_energy_loss_have_their_closed_form():
    """For ξ = 0 the integrals follow from the hyperbola in polar form: with cos φ₀ = 1/ε,
    I₀ ∝ 2(ε sin φ₀ − φ₀) and I±₂ ∝ ε(sin φ₀ + sin 3φ₀ / 3) − sin 2φ₀."""
    for theta in (40.0, 90.0, 150.0):
        eps = 1 / math.sin(math.radians(theta) / 2)
        p0 = math.acos(1 / eps)
        i = orbit_integrals(2, eps, 0.0)
        c0 = 2 * (eps * math.sin(p0) - p0)
        c2 = eps * (math.sin(p0) + math.sin(3 * p0) / 3) - math.sin(2 * p0)
        assert i[2].real / i[4].real == pytest.approx(c0 / c2, rel=1e-6)
        assert i[0] == pytest.approx(i[4]) and np.abs(i.imag).max() < 1e-12


# -- the angular correlation ----------------------------------------------------------------------------------------


@pytest.mark.parametrize("theta_cm", [30.0, 75.0, 120.0, 165.0])
def test_the_correlation_of_0_2_0_is_the_radiation_pattern_of_the_amplitudes(theta_cm):
    """For 0⁺ → 2⁺ → 0⁺ the γ ray carries the 2⁺ state's amplitudes as quadrupole radiation: its distribution is
    |Σ a_M X_2M|² with the vector spherical harmonics (Jackson). The tensor formula must give the same."""
    ex = ni()
    state = ex.at(theta_cm)
    a = ex.amplitudes(theta_cm)[1][:, 0]
    theta, phi, _ = sphere(10)
    assert np.allclose(state.w(1, 0, theta, phi), angular.multipole_pattern(a, theta, phi), rtol=1e-10, atol=1e-18)
    # Only M = 0 and ±2 are reached, in the frame whose z is perpendicular to the orbit, and with real amplitudes.
    assert a[1] == 0 and a[3] == 0 and abs(a[0]) > 0 and np.abs(a.imag).max() < 1e-12 * np.abs(a).max()


def test_the_closed_form_in_the_sudden_limit():
    """Without energy loss (ξ = 0) the amplitudes are a_M ∝ Y_2M(π/2, 0) I_−M with the closed-form integrals,
    so the whole correlation is known in closed form. A state at 1 keV is as good as ξ = 0."""
    ex = Excitation(*NI, simple_scheme("58Ni", "0.001 MeV", "E2", "0.0695 e2b2"))
    theta, phi, _ = sphere(10)
    for theta_cm in (50.0, 110.0, 160.0):
        eps = 1 / math.sin(math.radians(theta_cm) / 2)
        p0 = math.acos(1 / eps)
        c0 = 2 * (eps * math.sin(p0) - p0)
        c2 = eps * (math.sin(p0) + math.sin(3 * p0) / 3) - math.sin(2 * p0)
        y0, y2 = -math.sqrt(5 / (16 * math.pi)), math.sqrt(15 / (32 * math.pi))
        a = np.array([y2 * c2, 0, y0 * c0, 0, y2 * c2])
        closed = angular.multipole_pattern(a, theta, phi) / np.sum(a**2)
        state = ex.at(theta_cm)
        assert np.allclose(state.w(1, 0, theta, phi) / state.gamma_yield(1, 0), closed, rtol=2e-4, atol=5e-6)  # ξ is 5e-4, not 0


def test_backscattering_gives_the_m0_pattern_about_the_beam():
    """At 180° the orbit is along the beam, so only m = 0 along the beam is populated and the γ rays of
    2⁺ → 0⁺ follow (15/8π) sin²α cos²α = 1 + (5/7)P₂ − (12/7)P₄ over 4π, α from the beam."""
    state = ni().at(180.0)
    alpha = np.linspace(0.05, 3.1, 9)
    for phi in (0.0, 1.3, 4.0):
        dirs = np.c_[np.sin(alpha) * np.cos(phi), np.sin(alpha) * np.sin(phi), np.cos(alpha)]
        w = state.w_lab(1, 0, dirs) / state.gamma_yield(1, 0)
        assert np.allclose(w, 15 / (8 * np.pi) * np.sin(alpha) ** 2 * np.cos(alpha) ** 2, rtol=1e-9, atol=1e-12)
        c = np.cos(alpha)
        p2, p4 = (3 * c**2 - 1) / 2, (35 * c**4 - 30 * c**2 + 3) / 8
        assert np.allclose(4 * np.pi * w, 1 + 5 / 7 * p2 - 12 / 7 * p4, atol=1e-9)


def test_the_tensor_formula_holds_for_any_state():
    """Complex amplitudes test the phase conventions that real ones cannot."""
    rng = np.random.default_rng(1)
    a = rng.normal(size=5) + 1j * rng.normal(size=5)
    state = Populated(ni(), 90.0, {1: 1.0}, {1: tensors_of(np.outer(a, a.conj()), 2.0)})
    theta, phi, _ = sphere(10)
    assert np.allclose(state.w(1, 0, theta, phi), angular.multipole_pattern(a, theta, phi), rtol=1e-10)


@pytest.mark.parametrize("theta_cm", [45.0, 140.0])
def test_integrated_over_all_directions_the_yield_is_unchanged(theta_cm):
    theta, phi, w = sphere()
    for ex in (ni(), Excitation(*NI, cascade_scheme(branch_to_ground=30.0, mixing=0.4))):
        state = ex.at(theta_cm)
        for i, f in ex.transitions():
            assert np.sum(state.w(i, f, theta, phi) * w) == pytest.approx(state.gamma_yield(i, f), rel=1e-10)
    # In the laboratory too: the forward throw moves γ rays, it does not make or lose any.
    state = ni().at(theta_cm)
    dirs = np.stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)], axis=-1).reshape(-1, 3)
    for velocity in (None, (0.03, -0.02, 0.05)):
        lab = state.w_lab(1, 0, dirs, phi_particle_deg=25.0, velocity=velocity)
        assert np.sum(lab * w.ravel()) == pytest.approx(state.gamma_yield(1, 0), rel=1e-6)


def test_isotropic_emission_is_a_switch():
    state = ni(isotropic=True).at(120.0)
    theta, phi, _ = sphere(6)
    assert np.allclose(state.w(1, 0, theta, phi), state.gamma_yield(1, 0) / (4 * np.pi))
    assert state.gamma_yield(1, 0) == pytest.approx(ni().at(120.0).gamma_yield(1, 0))


def test_the_orbit_frame_in_the_laboratory():
    state = ni().at(60.0)
    for phi in (0.0, 70.0, 200.0):
        ax = state.axes(phi)
        assert np.allclose(ax.T @ ax, np.eye(3), atol=1e-12) and np.linalg.det(ax) == pytest.approx(1.0)
        beam, out = ax.T @ [0, 0, 1], ax.T @ [math.sin(math.radians(60)) * math.cos(math.radians(phi)),
                                               math.sin(math.radians(60)) * math.sin(math.radians(phi)), 0.5]
        assert np.allclose(beam, [-0.5, math.cos(math.radians(30)), 0], atol=1e-12)
        assert np.allclose(out, [0.5, math.cos(math.radians(30)), 0], atol=1e-12)
    # Turning the particle about the beam turns the γ-ray pattern with it.
    d = np.array([[0.3, 0.5, 0.81]])
    turned = np.array([[0.3 * math.cos(1.0) - 0.5 * math.sin(1.0), 0.3 * math.sin(1.0) + 0.5 * math.cos(1.0), 0.81]])
    assert state.w_lab(1, 0, d, 0.0) == pytest.approx(state.w_lab(1, 0, turned, math.degrees(1.0)))
    # For projectile excitation the exciting nucleus is on the other side.
    ex = Excitation("16O", "58Ni", "60 MeV", simple_scheme("16O", "6.917 MeV", "E2", "0.004 e2b2"),
                    excite="projectile")
    assert np.allclose(ex.at(60.0).axes(0.0)[:, :2], -state.axes(0.0)[:, :2])


def test_the_forward_throw_of_a_moving_emitter():
    state = ni(isotropic=True).at(100.0)
    beta = 0.05
    along, against, across = (state.w_lab(1, 0, [d], velocity=(0, 0, beta))[0] for d in ([0, 0, 1], [0, 0, -1],
                                                                                         [1, 0, 0]))
    flat = state.gamma_yield(1, 0) / (4 * np.pi)
    assert along / flat == pytest.approx((1 + beta) / (1 - beta))
    assert against / flat == pytest.approx((1 - beta) / (1 + beta))
    assert across / flat == pytest.approx(1 - beta**2)


# -- decay and feeding ----------------------------------------------------------------------------------------------


def test_conversion_takes_its_share_of_the_gamma_rays():
    alpha = 0.25
    plain, converted = Excitation(*NI, cascade_scheme()), Excitation(*NI, cascade_scheme(alpha=alpha))
    for theta in (60.0, 150.0):
        a, b = plain.at(theta), converted.at(theta)
        assert b.population(1) == pytest.approx(a.population(1))
        assert b.gamma_yield(1, 0) == pytest.approx(a.gamma_yield(1, 0) / (1 + alpha))
    cs = converted.cross_sections()
    assert cs["gamma"][(1, 0)] == pytest.approx(cs["populated"][1] / (1 + alpha))


def test_feeding_adds_the_population_of_the_levels_above():
    ex = Excitation(*NI, cascade_scheme(branch_to_ground=25.0))
    state = ex.at(140.0)
    to_first = 100.0 / 125.0
    assert state.population(1) == pytest.approx(state.direct[1] + to_first * state.direct[3])
    assert state.gamma_yield(3, 0) == pytest.approx(0.2 * state.direct[3])
    assert state.gamma_yield(3, 1) == pytest.approx(to_first * state.direct[3])
    cs = ex.cross_sections()
    assert cs["populated"][1] == pytest.approx(cs["direct"][1] + to_first * cs["direct"][3])
    assert 2 not in cs["direct"]                                  # the 4⁺ state is not reached in one step


def test_feeding_carries_the_orientation_down_with_u_k():
    rng = np.random.default_rng(3)
    m = rng.normal(size=(9, 9)) + 1j * rng.normal(size=(9, 9))
    rho = m @ m.conj().T
    parent = tensors_of(rho, 4.0)
    child = tensors_of(fed_density(rho, 4.0, 2.0, 2), 2.0)
    for k in (0, 2, 4):
        for q in range(-k, k + 1):
            assert child[(k, q)] == pytest.approx(angular.u_coefficient(k, 4, 2, 2) * parent[(k, q)])
    # In a scheme: the first 2⁺ state fed only from the second has the second's tensors times U_k(2 2 2).
    s = cascade_scheme()
    s.set_matrix_element(0, 1, "E2", 1e-9)
    state = Excitation(*NI, s).at(130.0)
    for k in (2, 4):
        ratio = state.tensors[1][(k, 0)] / state.tensors[3][(k, 0)]
        assert ratio == pytest.approx(angular.u_coefficient(k, 2, 2, 2), rel=1e-6)


def test_a_mixed_transition_uses_its_mixing_ratio():
    theta, phi, _ = sphere(8)
    pure = Excitation(*NI, cascade_scheme()).at(130.0)
    mixed = Excitation(*NI, cascade_scheme(mixing=0.6)).at(130.0)
    assert not np.allclose(pure.w(3, 1, theta, phi), mixed.w(3, 1, theta, phi), rtol=1e-3)
    a2 = mixed.coefficients(3, 1)[(2, 0)] * mixed.tensors[3][(0, 0)] / mixed.tensors[3][(2, 0)]
    assert a2 == pytest.approx(angular.distribution_coefficient(2, 2, 2, 1, 2, 0.6))


# -- in an experiment -----------------------------------------------------------------------------------------------


def test_the_correlation_table_of_the_example():
    exp = Experiment.example("coulex_ni58")
    rows = {(r["particle_detector"], r["gamma_detector"]): r["factor"] for r in correlation_table(exp)}
    assert len(rows) == 12
    # Behind the target the particle comes straight back: γ rays at 90° to the beam are few, those near 45° and
    # 135° many (the sin² cos² pattern, washed out by the size of the detectors).
    assert rows[("CD", "Ge90")] < 0.6 and rows[("CD", "Ge45")] > 1.3 and rows[("CD", "Ge135")] > 1.3
    # Left and right detectors, mirror images of each other, see the same.
    for g in ("Ge90", "Ge45", "Ge135", "Ge0"):
        assert rows[("DSSD-L", g)] == pytest.approx(rows[("DSSD-R", g)], rel=0.02)
    d = exp.to_dict()
    d["reaction"]["emission"] = "isotropic"
    flat = correlation_table(Experiment.from_dict(d))
    assert all(abs(r["factor"] - 1) < 0.06 for r in flat)        # only the forward throw is left
    assert excitation_of(exp).isotropic is False and excitation_of(Experiment.from_dict(d)).isotropic is True


def test_particle_rates_do_not_change_and_coincidences_use_the_correlation():
    p = Planner.example("coulex_ni58")
    r = p.rates()
    q = Planner.example("coulex_ni58")
    assert q.set("reaction", "emission", "isotropic")
    flat = q.rates()
    for row, other in zip(r["rows"], flat["rows"]):
        assert row["rate_per_s"] == other["rate_per_s"] and row["excitation_per_s"] == other["excitation_per_s"]
        assert row["coincidence_per_s"] == pytest.approx(row["excitation_per_s"] * row["gamma_efficiency"])
        assert other["gamma_efficiency"] == pytest.approx(flat["gamma_efficiency"], rel=0.03)
    cd = next(row for row in r["rows"] if row["detector"] == "CD")
    assert cd["gamma_efficiency"] != pytest.approx(r["gamma_efficiency"], rel=0.05)
    g = p.gamma()
    assert g["emission"] == "correlated" and len(g["correlation"]) == 12
    assert not q.set("reaction", "emission", "sideways")


def test_populations_of_a_level_scheme_in_the_planner():
    p = Planner.example("coulex_ni58")
    assert p.populations() == {"available": False, "reason": "Look up a level scheme first."}
    d = p.experiment.to_dict()
    d["levels"] = {"target": cascade_scheme(alpha=0.1, branch_to_ground=25.0).to_dict()}
    p = Planner(Experiment.from_dict(d))
    pop = p.populations()
    assert pop["available"] and pop["role"] == "target" and pop["nuclide"] == "58Ni"
    by = {x["index"]: x for x in pop["levels"]}
    assert set(by) == {1, 3} and by[1]["populated_mb"] > by[1]["direct_mb"] > by[3]["direct_mb"] > 0
    first = next(x for x in pop["gammas"] if (x["initial"], x["final"]) == (1, 0))
    assert first["gamma_mb"] == pytest.approx(by[1]["populated_mb"] / 1.1)
    assert first["energy_kev"] == pytest.approx(500.0)
