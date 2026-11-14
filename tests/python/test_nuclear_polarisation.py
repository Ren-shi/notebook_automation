"""The E1 polarisation correction (GOSIA's form, manual eq. 3.27–3.28, after Alder and Winther, appendix J): every
E2 collision function is multiplied by (1 − z a/r) along the orbit, z = d E_cm A / Z² of the excited nucleus. It
vanishes with d = 0; its size matches a case worked by hand; it weakens the reorientation of the ⁶²Ni 2⁺ state
less than the excitation; the setup's switch (on by default, with GOSIA's d = 0.0038) reaches the events, the plan
and the analysis, and the analysis gives back the B(E2) a run was made with when both have the correction. The
numbers below are for the d of GOSIA's manual, 0.005, unless they say otherwise."""

import copy

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.analysis import Analysis, Settings
from physim.nuclear.coulex import (GOSIA_DIPOLE_POLARISATION, GOSIA_MANUAL_DIPOLE_POLARISATION, Coulex,
                                   polarisation_strength)
from physim.nuclear.coupled import CoupledChannels
from physim.nuclear.gamma import (cascade_excitation, e1_polarisation, excitation_of, other_excitation,
                                  state_transition)
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.multistep import Multistep
from physim.nuclear.orientation import Excitation, simple_scheme
from physim.nuclear.planner import Planner
from physim.nuclear import data
from physim.nuclear.rates import Rates, channels, coulex_for, stack

D = GOSIA_MANUAL_DIPOLE_POLARISATION  # 0.005, GOSIA's manual
OFF = 0.0
NI58 = ("16O", "58Ni", 30.0)
NI62 = ("62Ni", "194Pt", 237.0)


def setup(name: str, d=None, order=None) -> Experiment:
    raw = Experiment.example(name).to_dict()
    if d is not None:
        raw["reaction"]["e1_polarisation"] = d
    if order is not None:
        raw["reaction"]["order"] = order
    return Experiment.from_dict(raw)


# -- the model ------------------------------------------------------------------------------------------------------


def test_z_is_gosias_for_target_and_projectile_excitation():
    o16, ni58, ni62, pt194 = (data.nuclide(n) for n in ("16O", "58Ni", "62Ni", "194Pt"))
    # Target excitation, GOSIA's manual eq. 3.28: z = d E_p A₂ / (Z₂² (1 + A₁/A₂)).
    assert polarisation_strength(D, o16, ni58, 30.0, "target") == pytest.approx(D * 30 * 58 / (28**2 * (1 + 16 / 58)))
    assert Coulex(*NI58, energy=1.454, b_up=695.0, e1_polarisation=D).polarisation == pytest.approx(0.0086976, rel=1e-4)
    # Projectile excitation, GOSIA's subroutine LOAD (IZ1 < 0): d E_p A_p / (Z_p² (1 + A_p/A_t)), that is
    # d E_cm A / Z² of the projectile, as for the target.
    z = polarisation_strength(D, ni62, pt194, 237.0, "projectile")
    assert z == pytest.approx(D * 237 * 62 / (28**2 * (1 + 62 / 194))) and z == pytest.approx(0.071016, rel=1e-4)
    assert polarisation_strength(D, ni62, pt194, 237.0, "target") == pytest.approx(D * 237 * 194
                                                                                  / (78**2 * (1 + 62 / 194)))
    assert polarisation_strength(0.0, ni62, pt194, 237.0, "projectile") == 0.0


def test_the_correction_vanishes_with_the_parameter_at_zero_and_leaves_e1_and_e3_alone():
    scheme = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2")
    for d in (0.0, None):
        assert Coulex(*NI58, energy=1.454, b_up=695.0, e1_polarisation=d).probability(150.0, exact=True) == \
            Coulex(*NI58, energy=1.454, b_up=695.0).probability(150.0, exact=True)
        assert Excitation(*NI58, scheme, e1_polarisation=d).probability(1, 150.0) == \
            Excitation(*NI58, scheme).probability(1, 150.0)
        cc = CoupledChannels(*NI58, scheme, e1_polarisation=d)
        assert cc.polarisation == 0.0 and cc.probability(1, 150.0) == CoupledChannels(*NI58, scheme).probability(1, 150.0)
    # Only E2 has it: GOSIA's QE applies the factor to λ = 2 alone.
    for mult, b in (("E1", 1e-3), ("E3", 1e4)):
        with_d = Coulex(*NI58, energy=3.0, multipolarity=mult, b_up=b, e1_polarisation=0.05)
        assert with_d.polarisation > 0
        assert with_d.probability(150.0, exact=True) == Coulex(*NI58, energy=3.0, multipolarity=mult,
                                                               b_up=b).probability(150.0, exact=True)
    assert Coulex(*NI58, energy=1.454, b_up=695.0, e1_polarisation=D).probability(150.0, exact=True) < \
        Coulex(*NI58, energy=1.454, b_up=695.0).probability(150.0, exact=True)


def test_its_size_in_sudden_backscattering_is_the_hand_worked_one():
    """At 180° (ε = 1) only μ = 0 couples, and with ξ → 0 the first-order amplitude is ∫ dw / (cosh w + 1)²
    = 2/3; the correction takes z ∫ dw / (cosh w + 1)³ = 4/15 z from it, so the probability falls by
    (1 − 0.4 z)². A level at 1 keV has ξ = 5 × 10⁻⁴."""
    scheme = simple_scheme("58Ni", "0.001 MeV", "E2", "1 e2fm4")
    for d in (D, 0.05):
        c0, c1 = (Coulex(*NI58, energy=0.001, b_up=1.0, e1_polarisation=x) for x in (0.0, d))
        expected = (1 - 0.4 * c1.polarisation) ** 2
        assert c0.xi < 1e-3
        assert c1.probability(180.0, exact=True) / c0.probability(180.0, exact=True) == pytest.approx(expected, abs=1e-6)
        e0, e1 = (Excitation(*NI58, scheme, e1_polarisation=x) for x in (0.0, d))
        assert e1.probability(1, 180.0) / e0.probability(1, 180.0) == pytest.approx(expected, abs=1e-6)
        k0, k1 = (CoupledChannels(*NI58, scheme, tolerance=1e-6, e1_polarisation=x) for x in (0.0, d))
        assert k1.probability(1, 180.0) / k0.probability(1, 180.0) == pytest.approx(expected, abs=1e-5)
        # Away from 180° the coupled equations still agree with first order for this weak state.
        assert k1.probability(1, 90.0) / k0.probability(1, 90.0) == pytest.approx(
            c1.probability(90.0, exact=True) / c0.probability(90.0, exact=True), abs=1e-5)


def test_the_reorientation_of_the_62ni_2plus_state_is_weakened_less_than_its_excitation():
    """The template's 2⁺ state of ⁶²Ni has ⟨2‖E2‖2⟩ = +6.6 e fm² (Q = +5 e fm²): the reorientation effect, the
    population with it over the population without it, is +8.1% at 170°. The correction weakens the diagonal
    coupling with the transitional one, so the effect keeps its sign and shrinks by only a few percent of itself,
    while the population falls by 6%."""
    scheme = setup("coulex_ni62_pt194").levels["beam"]
    assert scheme.matrix_element(1, 1, "E2").value == pytest.approx(6.597, abs=1e-3)
    flat = copy.deepcopy(scheme)
    flat.set_matrix_element(1, 1, "E2", 0.0, source="user")
    solve = {(d, s is flat): CoupledChannels(*NI62, s, excite="projectile", tolerance=1e-6, e1_polarisation=d)
             for d in (0.0, D) for s in (scheme, flat)}
    for theta in (150.0, 170.0):
        p = {k: cc.probability(1, theta) for k, cc in solve.items()}
        effect = {d: p[(d, False)] / p[(d, True)] - 1 for d in (0.0, D)}
        assert 0.07 < effect[0.0] < 0.09
        assert 0.97 < effect[D] / effect[0.0] < 0.995, (theta, effect)
        assert 0.92 < p[(D, False)] / p[(0.0, False)] < 0.95
    # The correction falls off one power of r faster than the excitation, so it weighs the backward angles more:
    # the ratio of the yields at 160° and 60° falls by about 2%, about what a diagonal element a third smaller
    # does without the correction.
    ratio = {d: solve[(d, False)].probability(1, 160.0) / solve[(d, False)].probability(1, 60.0) for d in (0.0, D)}
    assert -0.025 < ratio[D] / ratio[0.0] - 1 < -0.015


def test_with_a_diagonal_element_the_solver_and_second_order_theory_agree_on_the_correction():
    """The Runge–Kutta integration applies the factor as the nested quadrature of second-order perturbation
    theory does (an independent integration of the same couplings), reorientation included."""
    s = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0064 e2b2")
    s.set_matrix_element(1, 1, "E2", 26.4, source="user")
    p = {}
    for d in (0.0, 0.05):
        cc = CoupledChannels(*NI58, s, tolerance=1e-6, e1_polarisation=d)
        a = cc.solve(170.0)[cc._level_slices[1], 0]
        first, second = cc.second_order_check(1, 170.0)
        p[d] = (np.sum(np.abs(a) ** 2), np.sum(np.abs(first + second) ** 2))
    change_solver, change_theory = p[0.05][0] / p[0.0][0], p[0.05][1] / p[0.0][1]
    assert change_solver < 0.95
    assert change_solver == pytest.approx(change_theory, abs=2e-3)


# -- the setup's switch ---------------------------------------------------------------------------------------------


def test_the_setup_switch_reaches_the_solutions_the_plan_and_the_gosia_file():
    p = Planner.example("coulex_ni62_pt194")
    # On by default, with GOSIA's own default: a setup that does not set it has d = 0.0038.
    assert p.experiment.excitation.e1_polarisation is None
    assert e1_polarisation(p.experiment) == GOSIA_DIPOLE_POLARISATION == 0.0038
    assert cascade_excitation(p.experiment).e1_polarisation == 0.0038
    assert "DIP,3.800" in Multistep(p.experiment, p.experiment.levels["beam"], role="beam").gosia_input().splitlines()
    assert not p.set("reaction", "e1_polarisation", 5) and "DIP," in p.problems[0]
    assert not p.set("reaction", "e1_polarisation", -0.001)
    assert p.set("reaction", "e1_polarisation", D)
    on, off = p.experiment, setup("coulex_ni62_pt194", OFF)
    assert Experiment.from_toml(on.to_toml()).excitation.e1_polarisation == D
    assert Experiment.from_toml(off.to_toml()).excitation.e1_polarisation == 0.0 and e1_polarisation(off) == 0.0
    # The solutions: the switch is in the cache key, and each nucleus gets its own z.
    a, b = cascade_excitation(off), cascade_excitation(on)
    assert a is not b and a.polarisation == 0.0 and isinstance(b, CoupledChannels)
    assert b.polarisation == pytest.approx(D * b.beam_energy * 194 / 256 * 62 / 28**2, rel=2e-3)
    pt = other_excitation(on)
    assert pt.polarisation == pytest.approx(D * pt.beam_energy * 194 / 256 * 194 / 78**2, rel=2e-3)
    first = cascade_excitation(setup("coulex_ni62_pt194", D, order="first"))
    assert type(first) is Excitation and first.e1_polarisation == D
    assert b.near(150.0).direct[1] / a.near(150.0).direct[1] == pytest.approx(0.937, abs=0.005)
    # The plan's rates follow the populations: the nickel is excited about 6% less often in the silicon.
    ra, rb = Rates(off), Rates(on)
    for det in (g.name for g in ra.array):
        if ra.rate(det, what="excitations") > 0:
            assert 0.92 < rb.rate(det, what="excitations") / ra.rate(det, what="excitations") < 0.96, det
    # One state, no scheme: the first-order Coulex of the rates and the angular correlation have it too.
    one, none = setup("coulex_ni58", D), setup("coulex_ni58", OFF)
    layers = stack(one)
    ch = next(c for c in channels(layers, one) if c.excitation is not None)
    c = coulex_for(one, ch, layers)  # at the beam energy in the middle of the target, 28.9 MeV
    assert c.polarisation == pytest.approx(D * c.kinematics.beam_energy_mev * 58 / 74 * 58 / 28**2)
    assert coulex_for(none, ch, layers).polarisation == 0.0
    assert excitation_of(one).e1_polarisation == D and excitation_of(none).e1_polarisation == 0.0
    # The GOSIA input sets DIP to the setup's (1000 d); off, a millionth, since GOSIA takes 0 as "use 0.0038".
    assert "DIP,5.000" in Multistep(on, on.levels["beam"], role="beam").gosia_input().splitlines()
    assert "DIP,0.001" in Multistep(off, off.levels["beam"], role="beam").gosia_input().splitlines()


def test_the_analysis_gives_back_the_b_e2_of_a_run_made_with_the_correction():
    on, off = setup("coulex_ni62_pt194", D), setup("coulex_ni62_pt194", OFF)
    g = simulate_gammas(on, 150_000, seed=3)
    with_d, without = Analysis(on, g), Analysis(off, g)
    # The all-orders fit solves the same equations as the events, correction included: at the B(E2) the run was
    # made with it expects the γ rays of the events' solution (to the rounding of B(E2) in the setup's text).
    level, ground = state_transition(on)
    truth = on.excitation.b_up_e2fm
    grid, y = with_d._scheme_yield(truth, level, 120.0, 170.0)
    events = cascade_excitation(on)
    expected = np.array([events._populated(float(t), events._densities(events.solve(float(t)))).gamma_yield(level, ground)
                         for t in grid])
    assert np.allclose(y[(level, ground)], expected, rtol=1e-5)
    _, y_off = without._scheme_yield(truth, level, 120.0, 170.0)
    assert 0.92 < float(np.mean(y[(level, ground)] / y_off[(level, ground)])) < 0.96
    # The first-order reading has it too.
    s = Settings(particle_gate="all")
    keep = np.ones(len(g), dtype=bool)
    assert 0.92 < with_d._mean_probability(s, keep) / without._mean_probability(s, keep) < 0.95
    # B(E2) back: with the correction on both sides, the value put in, within the run's statistics (as without
    # it, test_nuclear_both_nuclei); leaving it out of the analysis reads about 10% less from the same run.
    r, r_off = with_d.run(s), without.run(s)
    assert r.b_e2fm4 == pytest.approx(truth, rel=0.08)
    assert 0.87 < r_off.b_e2fm4 / r.b_e2fm4 < 0.93
