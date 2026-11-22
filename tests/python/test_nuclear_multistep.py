"""Multi-step Coulomb excitation, reorientation and the shape (backlog 58)."""

import math
import time

import numpy as np
import pytest

from physim.external import gosia, programs
from physim.nuclear import Experiment
from physim.nuclear.coupled import CoupledChannels
from physim.nuclear.coulex import Coulex
from physim.nuclear.levels import Level, LevelScheme, Transition, Value, quadrupole_factor
from physim.nuclear.multistep import GOSIA_4PI_FORT9, Multistep
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


def test_all_orders_agree_with_gosia():
    """GOSIA 20110524.14 on the same case to all orders (tests/reference/nuclear/coulex_gosia_all_orders.csv): within
    1% at five angles, and its step from first order to all orders is physim's within 0.1 percentage point."""
    from physim.nuclear import validation

    ref = validation.reference_dir()
    path = ref / "coulex_gosia_all_orders.csv" if ref else None
    if path is None or not path.exists():
        pytest.skip("no GOSIA reference")
    rows = validation.read_reference(path)[1]
    first = {float(r["theta_cm_deg"]): float(r["probability"])
             for r in validation.read_reference(ref / "coulex_gosia.csv")[1]}
    real = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2")
    cc, ex = CoupledChannels(*NI, real, tolerance=1e-7), Excitation(*NI, real)
    for r in rows:
        theta, gosia = float(r["theta_cm_deg"]), float(r["probability"])
        ours = cc.probability(1, theta)
        assert ours == pytest.approx(gosia, rel=0.01), theta
        step_ours = ours / ex.probability(1, theta) - 1
        step_gosia = gosia / first[theta] - 1
        assert step_ours == pytest.approx(step_gosia, abs=1e-3), theta


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


#: ⁶²Ni at 111.6 MeV on ~1 mg/cm² carbon, carbon recoils counted in two full rings at 14–28° and 28–44°
#: (Allmond et al., PRC 90, 034309 (2014)): projectile excitation in inverse kinematics, the beam never reaching
#: the rings.
NI62_ON_C = """
schema = "physim.experiment/1"
title = "62Ni on 12C, recoils detected"
[reaction]
type = "coulex"
excite = "projectile"
energy = "1172.91 keV"
multipolarity = "E2"
b_up = "884.111 e2fm4"
[beam]
nuclide = "62Ni"
energy = "111.6 MeV"
current = "5 ppA"
charge_state = 10
[target]
material = "12C"
thickness = "1.03741 mg/cm2"
[run]
beam_time = "24 h"
[[detectors]]
name = "Ring2"
shape = "annular"
theta = "0 deg"
distance = "50 mm"
inner_radius = "12.4664 mm"
outer_radius = "26.5855 mm"
rings = 7
sectors = 10
thickness = "1000 um"
[[detectors]]
name = "Ring3"
shape = "annular"
theta = "0 deg"
distance = "50 mm"
inner_radius = "26.5855 mm"
outer_radius = "48.2844 mm"
rings = 7
sectors = 12
thickness = "1000 um"
"""


def ni62_scheme() -> LevelScheme:
    s = LevelScheme("62Ni")
    s.levels = [Level(v(0.0), 0.0, 1), Level(v(1172.91), 2.0, 1)]
    s.transitions = [Transition(1, 0, v(1172.91), v(100.0), "E2")]
    s.set_matrix_element(0, 1, "E2", 29.734)       # ⟨0⁺‖E2‖2⁺⟩ = 0.29734 e b, Q(2⁺) = 0
    return s


def expt_lines(text: str) -> list:
    lines = text.splitlines()
    i = lines.index("EXPT")
    return lines[i + 2:i + 2 + int(lines[i + 1].split()[0])]


def test_the_gosia_input_is_written_from_the_scheme():
    exp = Experiment.example("coulex_ni58")
    ms = Multistep(exp, two_levels(26.4, q_efm2=15.0), "target")
    text = ms.gosia_input(gamma="4pi")
    assert text.startswith("OP,TITL") and "OP,GOSI" in text and "LEVE" in text and "ME" in text
    assert "1 +1 0 0.000000" in text and "2 +1 2 1.454200" in text
    assert "1 2 0.264000" in text                              # e fm² → e b
    assert "OP,YIEL" in text and "OP,INTI" in text and text.rstrip().endswith("OP,EXIT")


def test_gosia_sign_rule_for_the_target_scheme():
    """GOSIA (20110524, CMLAB and INVKIN): a positive partner Z means the TARGET is the investigated nucleus. The
    app wrote it negative before, which GOSIA reads as ⁵⁸Ni the projectile on ¹⁶O and stops ("maximum scattering
    angle is 14.06 degrees")."""
    exp = Experiment.example("coulex_ni58")
    text = Multistep(exp, two_levels(26.4), "target").gosia_input(gamma="4pi")
    lines = text.splitlines()
    assert lines[lines.index("EXPT") + 1].split()[1:] == ["28", "58"]          # the investigated nucleus: ⁵⁸Ni
    expt = expt_lines(text)
    assert expt and all(line.split()[:2] == ["8", "16"] for line in expt)    # the beam, Z positive
    assert all(float(line.split()[2]) == pytest.approx(30.0) for line in expt)


def test_gosia_sign_rule_for_the_beam_scheme():
    """The beam's scheme: the partner (the target) is written with a negative Z, the projectile being the
    investigated nucleus."""
    exp = Experiment.from_toml(NI62_ON_C)
    text = Multistep(exp, ni62_scheme(), "beam").gosia_input(gamma="4pi")
    lines = text.splitlines()
    assert lines[lines.index("EXPT") + 1].split()[1:] == ["28", "62"]          # the investigated nucleus: ⁶²Ni
    assert all(line.split()[:2] == ["-6", "12"] for line in expt_lines(text))


def test_gosia_cont_options_and_the_blank_line_after_end():
    exp = Experiment.from_toml(NI62_ON_C)
    text = Multistep(exp, ni62_scheme(), "beam").gosia_input(gamma="4pi")
    lines = text.splitlines()
    cont = lines[lines.index("CONT"):lines.index("OP,YIEL")]
    assert cont == ["CONT", "ACC,7.", "DIP,3.8", "SPL,1.", "END,", ""]      # splines; GOSIA needs the blank line


def test_gosia_recoil_detection_in_inverse_kinematics():
    """⁶²Ni on carbon with the carbon recoils counted: the beam (at most ~10.8°) cannot reach the rings, so each
    ring is one experiment with the recoil detected: a negative angle on the EXPT line that is the BEAM's lab
    angle, and OP,INTI's mesh in recoil angles (OP,INTG would want the beam's angles, two-valued in ring 3)."""
    exp = Experiment.from_toml(NI62_ON_C)
    ms = Multistep(exp, ni62_scheme(), "beam")
    xs = ms.gosia_experiments()
    assert [(x["detector"], x["particle"]) for x in xs] == [("Ring2", "recoil"), ("Ring3", "recoil")]
    assert all(x["phi"] is None for x in xs)                                  # full rings: axial symmetry
    assert xs[0]["theta"] == pytest.approx((14.0, 28.0), abs=1e-3)
    assert xs[1]["theta"] == pytest.approx((28.0, 44.0), abs=1e-3)
    text = ms.gosia_input(gamma="4pi")
    for line, x in zip(expt_lines(text), xs):
        f = line.split()
        theta = float(f[3])
        assert theta < 0 and -theta < 11.2                                    # the beam's angle, below its maximum
        assert f[6:9] == ["0", "0.00", "360.00"] and f[9] == str(x["ikin"])
    inti = text.splitlines()[text.splitlines().index("OP,INTI") + 1:]
    head = inti[0].split()
    assert head[:2] == ["9", "9"] and float(head[4]) == pytest.approx(14.0, abs=1e-3)
    assert float(head[2]) < 75.0 and float(head[3]) == pytest.approx(111.6)    # exit to entry energy
    assert [float(t) for t in inti[2].split()][-1] == pytest.approx(28.0, abs=1e-3)


def test_gosia_side_detectors_cover_their_azimuth():
    """Square DSSDs at the side: IAX = 1, a negative angle count in OP,INTI and the azimuths covered at each angle;
    forward of 90° the ⁵⁸Ni recoil reaches them too, so each has a second experiment with the recoil detected."""
    exp = Experiment.example("coulex_ni58")
    ms = Multistep(exp, two_levels(26.4), "target")
    xs = ms.gosia_experiments()
    assert [(x["detector"], x["particle"]) for x in xs] == [
        ("CD", "ejectile"), ("DSSD-L", "ejectile"), ("DSSD-L", "recoil"), ("DSSD-R", "ejectile"),
        ("DSSD-R", "recoil")]
    assert xs[0]["phi"] is None and all(x["phi"] is not None for x in xs[1:])
    expt = expt_lines(ms.gosia_input(gamma="4pi"))
    assert expt[0].split()[6] == "0" and all(line.split()[6] == "1" for line in expt[1:])
    assert all(float(line.split()[3]) < 0 for line in (expt[2], expt[4]))


def test_gosia_input_with_the_gamma_array_runs_gdet_first():
    exp = Experiment.example("coulex_ni58")
    if not exp.gamma_detectors:
        with pytest.raises(ValueError, match="no γ detectors"):
            Multistep(exp, two_levels(26.4), "target").gosia_input()
        return
    text = Multistep(exp, two_levels(26.4), "target").gosia_input()
    lines = text.splitlines()
    assert lines.index("OP,GDET") < lines.index("OP,GOSI")
    assert lines[lines.index("OP,GDET") + 1] == str(len(exp.gamma_detectors))


def test_gosia_integrated_output_is_read_back():
    out = """
     INTEGRATED RUTHERFORD CROSS SECTION=.5746E+03  FOR EXP. 1

                                                  INTEGRATED YIELDS

     EXPERIMENT  1  DETECTOR  1
     ENERGY RANGE   69.209--- 111.600 MEV   SCATTERING ANGLE RANGE  14.000--- 28.000 DEG

     NI     NF     II     IF     YIELD     NORMALIZED YIELD

      2      1    2.0    0.0   0.41474E-01   0.10000E+01
"""
    res = gosia.integrated(out)
    assert res[1]["rutherford"] == pytest.approx(574.6)
    assert res[1]["yields"][1][(2, 1)] == pytest.approx(0.041474)


def _gosia_exe():
    folder = programs.find("gosia")
    exe = folder / "gosia.exe" if folder is not None else None
    return exe if exe is not None and exe.is_file() else None


@pytest.mark.skipif(_gosia_exe() is None, reason="GOSIA is not on this machine")
def test_gosia_run_reproduces_the_allmond_ratios(tmp_path):
    """The app's input run through gosia.exe gives the γ yield over GOSIA's (inelastic) integrated Rutherford of an
    independent hand-written GOSIA input for the same setup (SRIM stopping and thickness; 9.076e-4 in ring 2 and
    6.647e-4 in ring 3) within 0.5 %. The app's thickness is tuned with its own stopping to the same energy loss."""
    exp = Experiment.from_toml(NI62_ON_C)
    ms = Multistep(exp, ni62_scheme(), "beam")
    (tmp_path / "fort.9").write_text(GOSIA_4PI_FORT9)
    res = gosia.integrated(gosia.run(_gosia_exe(), ms.gosia_input(gamma="4pi"), tmp_path).output)
    ratios = [4 * math.pi * res[k]["yields"][1][(2, 1)] / res[k]["rutherford"] for k in (1, 2)]
    assert ratios == pytest.approx([9.076e-4, 6.647e-4], rel=5e-3)


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
