"""Multi-step Coulomb excitation, reorientation and the shape (backlog 58)."""

import math
import time

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.coupled import CoupledChannels
from physim.nuclear.coulex import Coulex
from physim.nuclear.levels import Level, LevelScheme, Transition, Value, quadrupole_factor
from physim.nuclear.multistep import Multistep
from physim.nuclear.orientation import Excitation, simple_scheme
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates

NI = ("16O", "58Ni", "30 MeV")


def v(x: float) -> Value:
    return Value(x, source="user")


def two_levels(me01: float = 26.4, q_efm2: float = 0.0) -> LevelScheme:
    s = LevelScheme("58Ni")
    s.levels = [Level(v(0.0), 0.0, 1), Level(v(1454.2), 2.0, 1)]
    s.transitions = [Transition(1, 0, v(1454.2), v(100.0), "E2")]
    s.set_matrix_element(0, 1, "E2", me01)
    if q_efm2:
        s.set_matrix_element(1, 1, "E2", q_efm2 / quadrupole_factor(2.0))
    return s


def band(n_levels: int = 25) -> LevelScheme:
    """A rotational-like band with couplings between neighbours, for the size and speed checks."""
    s = LevelScheme("194Pt")
    s.levels = [Level(v(0.0), 0.0, 1)] + [Level(v(300.0 + 120 * i), float(2 * ((i % 4) + 1)) if i % 4 else 2.0, 1)
                                        for i in range(n_levels - 1)]
    s.transitions = [Transition(i + 1, 0, v(300.0 + 120 * i), v(100.0), "E2") for i in range(n_levels - 1)]
    s.transitions += [Transition(i + 1, i, v(120.0), v(100.0), "E2") for i in range(1, n_levels - 1)]
    for i in range(n_levels - 1):
        if s.levels[i + 1].spin == 2.0:
            s.set_matrix_element(0, i + 1, "E2", 100.0 / (1 + i))
        if i > 0:
            s.set_matrix_element(i, i + 1, "E2", 80.0)
    return s


# -- the coupled equations ------------------------------------------------------------------------------------------


def test_with_weak_coupling_the_result_is_first_order():
    weak = simple_scheme("58Ni", "1.454 MeV", "E2", "0.00000695 e2b2")   # a ten-thousandth of the real strength
    ex, cc = Excitation(*NI, weak), CoupledChannels(*NI, weak, tolerance=1e-6)
    for theta in (30.0, 90.0, 170.0):
        assert cc.probability(1, theta) == pytest.approx(ex.probability(1, theta), rel=2e-6), theta
    # The deviation from first order grows with the strength: the real ⁵⁸Ni state is 1.5% below first order at
    # 170°, from the ground state's depletion and the paths back and forth.
    real = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2")
    ratio = CoupledChannels(*NI, real).probability(1, 170.0) / Excitation(*NI, real).probability(1, 170.0)
    assert 0.97 < ratio < 0.995


def test_the_probabilities_of_all_levels_sum_to_one():
    cc = CoupledChannels("16O", "194Pt", "60 MeV", band(8))
    for theta in (40.0, 100.0, 160.0):
        p = cc.probabilities(theta)
        assert sum(p.values()) == pytest.approx(1.0, abs=1e-8), theta
        assert all(x >= 0 for x in p.values()) and p[0] < 1.0


def test_a_diagonal_matrix_element_reproduces_second_order_perturbation_theory():
    """With a quadrupole moment the second-order path 0 → 2 → 2 interferes with the first-order amplitude: the
    reorientation effect. The solver's amplitude must equal first plus second order, by nested quadrature on
    the same orbit, to the size of the third-order term."""
    # A weaker transition keeps the third order (the ground state's depletion) out of the way; the size of the
    # reorientation effect depends on the moment alone.
    s = two_levels(8.0, q_efm2=20.0)
    cc = CoupledChannels(*NI, s, tolerance=1e-6)
    for theta in (120.0, 170.0):
        a = cc.solve(theta)[cc._level_slices[1], 0]
        first, second = cc.second_order_check(1, theta)
        assert np.max(np.abs(second)) > 0.02 * np.max(np.abs(first))     # the effect is there
        p_solver = np.sum(np.abs(a) ** 2)
        p_12 = np.sum(np.abs(first + second) ** 2)
        p_1 = np.sum(np.abs(first) ** 2)
        assert abs(p_solver - p_12) < 0.05 * abs(p_12 - p_1), (theta, p_solver, p_12, p_1)
    # The sign of the moment sets the sign of the effect: prolate (negative Q_s) lowers the yield here, oblate
    # raises it, against Q = 0.
    p = {q: CoupledChannels(*NI, two_levels(26.4, q)).probability(1, 150.0) for q in (-20.0, 0.0, 20.0)}
    assert p[-20.0] < p[0.0] < p[20.0]
    assert (p[20.0] - p[-20.0]) / p[0.0] > 0.05


def test_two_step_excitation_reaches_a_level_first_order_cannot():
    s = band(3)                       # 0⁺, 2⁺, 4⁺: the 4⁺ only through the 2⁺
    s.matrix_elements = [me for me in s.matrix_elements if not (me.a == 0 and me.b == 2)]
    cc, ex = CoupledChannels("16O", "194Pt", "60 MeV", s), Excitation("16O", "194Pt", "60 MeV", s)
    assert ex.probability(2, 160.0) == 0.0
    assert cc.probability(2, 160.0) > 1e-4
    state = cc.at(160.0)
    assert state.population(2) > 0 and state.gamma_yield(2, 1) == pytest.approx(state.population(2) / 2)
    assert state.population(1) > state.direct[1] * 0.99     # feeding from the 4⁺ adds to the 2⁺


def test_a_25_level_case_is_fast_enough():
    cc = CoupledChannels("16O", "194Pt", "60 MeV", band(25), tolerance=1e-4)
    assert len(cc.basis) == 265
    t = time.perf_counter()
    p = cc.probabilities(150.0)
    assert time.perf_counter() - t < 5.0
    assert sum(p.values()) == pytest.approx(1.0, abs=1e-8)


# -- the experiment ---------------------------------------------------------------------------------------------


def test_first_order_yields_of_the_experiment_match_the_analytic_rates():
    """The same first-order excitation, integrated over the detectors the Multistep way, against Rates: equal
    but for the Rutherford factor (symmetrised orbit here, elastic there, 6% at this energy)."""
    exp = Experiment.example("coulex_ni58")
    y = Multistep(exp, two_levels(math.sqrt(695.0)), "target", angle_step=6, energies=2, first_order=True).yields()
    r = Rates(exp)
    for d in y.detectors:
        ratio = y.detectors[d][(1, 0)] / r.rate(d, what="excitations", counted=False)
        assert 1.0 < ratio < 1.12, (d, ratio)
    assert set(y.rings["CD"]) == set(range(exp.detectors[0].rings))
    assert sum(y.rings["CD"][k][(1, 0)] for k in y.rings["CD"]) == pytest.approx(y.detectors["CD"][(1, 0)])


def test_all_orders_differ_from_first_order_and_the_shapes_are_told_apart():
    exp = Experiment.example("coulex_ni58")
    ms = Multistep(exp, two_levels(math.sqrt(695.0)), "target", angle_step=8, energies=1)
    full, first = ms.yields(), Multistep(exp, ms.scheme, "target", 8, 1, first_order=True).yields()
    for d in full.detectors:
        assert 0.95 < full.detectors[d][(1, 0)] / first.detectors[d][(1, 0)] < 1.0
    sh = ms.shapes(1)
    assert sh["transition"] == (1, 0) and sh["q_efm2"]["prolate"] < 0 < sh["q_efm2"]["oblate"]
    assert sh["totals"]["prolate"] < sh["totals"]["spherical"] < sh["totals"]["oblate"]
    assert sh["separable"] and sh["total_difference_sigma"] > 3
    assert all(r["difference_sigma"] >= 0 for r in sh["rings"]) and len(sh["rings"]) > 10


def test_the_fit_recovers_the_matrix_elements_put_in():
    exp = Experiment.example("coulex_ni58")
    truth = two_levels(26.4, q_efm2=15.0)
    ms = Multistep(exp, truth, "target", angle_step=10, energies=1)
    t = 86400.0
    counts = ms.yields()
    measured = {d: {(1, 0): (counts.detectors[d][(1, 0)] * t, 0.01 * counts.detectors[d][(1, 0)] * t)}
                for d in counts.detectors}
    start = two_levels(22.0, q_efm2=3.0)
    r = Multistep(exp, start, "target", angle_step=10, energies=1).fit(measured, [(0, 1, "E2"), (1, 1, "E2")])
    assert r["values"][(0, 1, "E2")] == pytest.approx(26.4, rel=0.02)
    assert r["values"][(1, 1, "E2")] == pytest.approx(15.0 / quadrupole_factor(2.0), rel=0.2)
    assert r["uncertainties"][(0, 1, "E2")] < 0.03 * 26.4 and r["chi2"] < 3 * max(r["degrees_of_freedom"], 1)
    with pytest.raises(ValueError, match="one to three"):
        ms.fit(measured, [])


def test_the_gosia_input_is_written_from_the_scheme():
    exp = Experiment.example("coulex_ni58")
    ms = Multistep(exp, two_levels(26.4, q_efm2=15.0), "target")
    text = ms.gosia_input()
    assert text.startswith("OP,TITL") and "OP,GOSI" in text and "LEVE" in text and "ME" in text
    assert "1 +1 0 0.000000" in text and "2 +1 2 1.454200" in text
    assert "1 2 0.26400" in text                               # e fm² → e b
    assert "EXPT" in text and "-8 16 30.000" in text          # target excited: the projectile's Z negative
    assert text.rstrip().endswith("OP,EXIT")


def test_the_planner_solves_all_orders_with_a_scheme():
    p = Planner.example("coulex_ni58")
    assert not p.multistep()["available"]
    d = p.experiment.to_dict()
    d["levels"] = {"target": two_levels(math.sqrt(695.0), q_efm2=-10.0).to_dict()}
    p = Planner(Experiment.from_dict(d))
    m = p.multistep(angle_step=10, energies=1)
    assert m["available"] and m["detectors"] == ["CD", "DSSD-L", "DSSD-R"]
    x = m["transitions"][0]
    assert x["transition"] == (1, 0) and all(0.9 < x["all_orders"][dd] / x["first_order"][dd] < 1.0
                                            for dd in m["detectors"])
    assert m["shapes"]["level"] == 1 and set(m["shapes"]["q_efm2"]) == {"prolate", "spherical", "oblate"}
    assert "OP,GOSI" in m["gosia_input"] and "OP,INTI" in m["gosia_input"]
    assert p.multistep(angle_step=10, energies=1) is m
