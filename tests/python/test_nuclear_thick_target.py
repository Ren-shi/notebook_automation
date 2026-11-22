"""Coulomb excitation in a thick target: the excitation probability integrated over the depth, and Q_s for a single
0⁺ → 2⁺ transition.

The reference case is Allmond et al., PRC 90, 034309 (2014): 62Ni at 111.6 MeV on 1.04 mg/cm² of carbon, the carbon
recoils counted in two rings (14–28° and 28–44° lab). The beam loses 42 MeV (38 %) in the target and P(θ) changes
150-fold through it. GOSIA (all orders, integrated over the target, elastic Rutherford) gives σ_C/σ_R = 1.0080e-3
(ring 2) and 7.699e-4 (ring 3) for <0⁺‖E2‖2⁺> = 0.29734 eb and Q = 0.

The Coulomb-excitation cross section carries Alder and Winther's inelastic Rutherford factor, the elastic one times
v_i/v_f (GOSIA's convention). Against those GOSIA numbers the app is then about 1.3 % (ring 2) and 1.2 % (ring 3)
high, and the whole of it is two simplifications of GOSIA's, measured by rerunning GOSIA through the app's own
writer (:data:`GOSIA_MATCHED`): its integer mass numbers (A = 62 and 12; the app uses the nuclear masses, and P,
steep in the energy here, moves ~0.6 %) and its E1 polarisation (DIP 3.8, −0.5 %), which the app does not model.
With the nuclear masses and DIP off GOSIA and the app agree within 0.2 %.
"""

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import rates as rates_module
from physim.nuclear.coulex import Coulex
from physim.nuclear.events import generator_config, simulate
from physim.nuclear.gamma import cascade_excitation, has_static_moment
from physim.nuclear.levels import quadrupole_factor
from physim.nuclear.rates import Rates, beam_energy_at, channels, depth_probability, stack
from physim.nuclear.rutherford import Rutherford

RINGS = ("Ring2", "Ring3")
#: GOSIA, all orders, integrated over the target, divided by the elastic Rutherford cross section.
GOSIA = {"Ring2": 1.0080e-3, "Ring3": 7.699e-4}
#: All orders against first order at Q = 0, from the coupled equations at mid-target (Round 2's builder). Through
#: the depth the app finds about −0.7 % and −0.4 %: the front of the target, where P and the higher orders are
#: largest, weighs most.
ALL_VS_FIRST = {"Ring2": -0.0049, "Ring3": -0.0028}
#: Against :data:`GOSIA` the app is +1.3 % (ring 2) and +1.2 % (ring 3): GOSIA's mass numbers and E1 polarisation.
GOSIA_TOLERANCE = 0.02
#: GOSIA (20110524.14) run on the app's own input (``Multistep.gosia_input(gamma="4pi")``, the app's stopping and
#: thickness) with the nuclear masses of ⁶²Ni and ¹²C (61.91298 u and 11.99671 u, for 62 and 12) and DIP,0.001 (E1
#: polarisation off): γ yield ×4π over the app's elastic Rutherford integral (R_el) and over GOSIA's own, inelastic,
#: integrated Rutherford (R_G).
GOSIA_MATCHED = {"Ring2": 1.01929e-3, "Ring3": 7.77906e-4}
GOSIA_MATCHED_RG = {"Ring2": 9.19657e-4, "Ring3": 6.73100e-4}
#: GOSIA's integrated Rutherford (mb·mg/cm², full azimuth) on the app's unmodified input (mass numbers, DIP 3.8):
#: the inelastic kinematics with its Rutherford factor (CMLAB), against which the app's own factor is checked.
GOSIA_RUTHERFORD = {"Ring2": 574.6, "Ring3": 1766.0}
#: The ring edges (recoil lab angles, degrees).
RING_DEG = {"Ring2": (14.0, 28.0), "Ring3": (28.0, 44.0)}
#: GOSIA's change of the ratio for Q_s = +0.1 eb (−0.1 eb is the same, of the other sign), in percent.
GOSIA_Q_EFFECT = {"Ring2": 2.7, "Ring3": 1.9}

SETUP = """
schema = "physim.experiment/1"
title = "62Ni on natC (Allmond 2014)"

[reaction]
type = "coulex"
excite = "projectile"
energy = "1172.91 keV"
multipolarity = "E2"
b_up = "884.111 e2fm4"
{order}

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
counts_wanted = 1000

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
dead_layer = "0.5 um"
resolution = "500 keV"
threshold = "1 MeV"

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
dead_layer = "0.5 um"
resolution = "500 keV"
threshold = "1 MeV"
{extra}
"""

LEVELS = """
[levels.beam]
nuclide = "62Ni"

[[levels.beam.level]]
energy = "0 keV"
jpi = "0+"

[[levels.beam.level]]
energy = "1172.91 keV"
jpi = "2+"
half_life = "1.442e-12 s"

[[levels.beam.transition]]
from = 1
to = 0
energy = "1172.89 keV"
multipolarity = "E2"

[[levels.beam.matrix_element]]
from = 0
to = 1
multipolarity = "E2"
value = "29.734 efm2"
{diagonal}
"""

GAMMA = """
[[gamma_detectors]]
name = "G90"
theta = "90 deg"
distance = "100 mm"
crystal_diameter = "60 mm"
crystal_length = "70 mm"
resolution = "2.1 keV"
"""


def allmond(q_eb=None, order=None, gamma=False) -> Experiment:
    """The setup; with ``q_eb`` a {0⁺, 2⁺} level scheme for the beam whose 2⁺ has Q_s = ``q_eb``."""
    extra = ""
    if q_eb is not None:
        me = q_eb * 100.0 / quadrupole_factor(2.0)
        extra = LEVELS.format(diagonal=f"""
[[levels.beam.matrix_element]]
from = 1
to = 1
multipolarity = "E2"
value = "{me:.12g} efm2"
""")
    if gamma:
        extra += GAMMA
    return Experiment.from_toml(SETUP.format(order=f'order = "{order}"' if order else "", extra=extra))


def ratios(exp, **kw) -> dict:
    rt = Rates(exp, **kw)
    return {r: rt.rate(r, counted=False, what="excitations") / rt.rate(r, counted=False) for r in RINGS}


def excited_channel(exp):
    layers = stack(exp)
    return next(ch for ch in channels(layers, exp) if ch.excitation is not None), layers


@pytest.fixture(scope="module")
def first_order():
    return ratios(allmond())


@pytest.fixture(scope="module")
def all_orders():
    return {q: ratios(allmond(q)) for q in (0.0, 0.1, -0.1)}


# -- the depth-integrated excitation probability ---------------------------------------------------------------


def test_depth_probability_follows_the_energy_through_the_layer():
    """The interpolated P at any energy in the layer is the first-order P solved at that energy, and the depth
    average is the brute-force average over the depth with the Rutherford weight (∝ 1/E_cm²)."""
    exp = allmond()
    ch, layers = excited_channel(exp)
    dp = depth_probability(exp, ch, layers)
    assert not dp.thin and dp.e_front == pytest.approx(111.6) and 69.0 < dp.e_back < 69.5
    theta = np.array([60.0, 90.0, 120.0, 150.0, 178.0])

    def exact(e):
        return Coulex("62Ni", "12C", e, excite="projectile", energy=1.17291, multipolarity="E2",
                      b_up=884.111).probability(theta)

    for e in (111.0, 100.0, 90.0, 80.0, 70.0):
        np.testing.assert_allclose(dp.at(e).probability(theta), exact(e), rtol=1e-3)
    zs = (np.arange(48) + 0.5) / 48 * dp.reached
    es = beam_energy_at(exp, 0, zs, layers)
    w = np.array([Rutherford("62Ni", "12C", float(e)).cross_section_cm(90.0) for e in es])
    brute = (w[:, None] * np.array([exact(float(e)) for e in es])).sum(axis=0) / w.sum()
    np.testing.assert_allclose(dp.probability(theta), brute, rtol=2e-3)
    # The mid-target P (what the app used before) is about 20 % short at the backward angles.
    mid = rates_module.layer_energy(exp, 0, layers)
    assert dp.probability(150.0) > 1.15 * float(exact(mid)[3])


def test_slices_average_to_the_whole_layer():
    exp = allmond()
    ch, layers = excited_channel(exp)
    dp = depth_probability(exp, ch, layers)
    n = dp.slice_count()
    assert 4 <= n <= rates_module.MAX_SLICES
    theta = np.array([90.0, 150.0])
    parts = dp.slices(n, theta)
    assert sum(t for _, t, _ in parts) == pytest.approx(dp.thickness)
    # Each slice weighs as its Rutherford integral; together they are the whole layer's average.
    wts = []
    for lo, t, _ in parts:
        x, w = np.polynomial.legendre.leggauss(8)
        es = beam_energy_at(exp, 0, lo + (x + 1) / 2 * t, layers)
        wts.append(float(np.sum(w * [Rutherford("62Ni", "12C", float(e)).cross_section_cm(90.0) for e in es])) * t)
    combined = sum(wt * p for wt, (_, _, p) in zip(wts, parts)) / sum(wts)
    np.testing.assert_allclose(combined, dp.probability(theta), rtol=1e-3)
    # The front slices, at the higher energies, excite more.
    assert parts[0][2][1] > 3 * parts[-1][2][1]


def test_rates_integrate_p_over_the_depth(first_order, monkeypatch):
    """Rates with P through the depth: converged with the default four depth points; with the elastic Rutherford
    factor in place of Alder and Winther's, the 40-slice workaround of Round 2 (9.903e-4 and 7.529e-4) is
    reproduced, and the inelastic factor (v_i/v_f, 1.04–1.10 through this target) raises both rings by 3.8 %."""
    fine = ratios(allmond(), depth_points=16)
    for r in RINGS:
        assert first_order[r] == pytest.approx(fine[r], rel=2e-4)
    assert first_order["Ring2"] == pytest.approx(1.0285e-3, rel=2e-3)
    assert first_order["Ring3"] == pytest.approx(7.820e-4, rel=2e-3)
    monkeypatch.setattr(rates_module, "rutherford_factor", lambda e_cm, e_star: 1.0)
    rates_module._DEPTH_KEPT.clear()
    elastic = ratios(allmond())
    rates_module._DEPTH_KEPT.clear()
    assert elastic["Ring2"] == pytest.approx(9.903e-4, rel=2e-3)
    assert elastic["Ring3"] == pytest.approx(7.529e-4, rel=2e-3)
    for r in RINGS:
        assert first_order[r] / elastic[r] == pytest.approx(1.038, abs=0.003)


@pytest.mark.parametrize("ring", RINGS)
def test_first_order_against_gosia(first_order, ring):
    """First order, corrected to all orders, against GOSIA: +1.1 to +1.5 %, GOSIA's mass numbers and E1
    polarisation (module docstring)."""
    ours = first_order[ring] * (1 + ALL_VS_FIRST[ring])
    assert ours == pytest.approx(GOSIA[ring], rel=GOSIA_TOLERANCE)
    assert 0.008 < ours / GOSIA[ring] - 1 < 0.02


def ring_integrals(exp, depth_points: int = 24, angle_points: int = 200) -> dict:
    """Per ring, over the target's depth and the ring's recoil angles (the first branch of the inelastic
    kinematics, the one GOSIA integrates): (∫σ_Coulex, ∫ Alder–Winther's Rutherford factor alone, ∫σ_R elastic),
    mb·mg/cm²."""
    layers = stack(exp)
    ch, _ = excited_channel(exp)
    dp = depth_probability(exp, ch, layers)
    x, w = np.polynomial.legendre.leggauss(depth_points)
    zs, wz = (x + 1) / 2 * dp.reached, w * dp.reached / 2
    es = beam_energy_at(exp, 0, zs, layers)
    xt, wt = np.polynomial.legendre.leggauss(angle_points)
    out = {}
    for r in RINGS:
        lo, hi = RING_DEG[r]
        th = lo + (xt + 1) / 2 * (hi - lo)
        dom = 2 * np.pi * np.sin(np.radians(th)) * wt * np.radians(hi - lo) / 2
        sums = np.zeros(3)
        for e, we in zip(es, wz):
            sums[2] += we * np.sum(np.nan_to_num(Rutherford("62Ni", "12C", float(e)).cross_section_lab(th, "recoil")[0])
                                   * dom)
            cx = Coulex("62Ni", "12C", float(e), excite="projectile", energy=1.17291, b_up=884.111)
            pt = cx.kinematics.at_lab(th, "recoil")[0]
            tej = 180.0 - np.asarray(pt.theta_cm)
            ok = ~np.isnan(tej)
            aw = np.where(ok, np.asarray(cx.rutherford_cm(np.where(ok, tej, 90.0)))
                          * np.nan_to_num(np.asarray(pt.jacobian)), 0.0) * dom
            sums[0] += we * np.sum(aw * np.asarray(dp.at(float(e)).probability(np.where(ok, tej, 90.0))))
            sums[1] += we * np.sum(aw)
        out[r] = sums
    return out


@pytest.fixture(scope="module")
def integrals_q0():
    return ring_integrals(allmond(0.0))


@pytest.mark.parametrize("ring", RINGS)
def test_rutherford_factor_is_gosias(integrals_q0, ring):
    """The app's Coulomb-excitation Rutherford factor (Alder and Winther: elastic × v_i/v_f, on the inelastic
    kinematics), integrated over the ring and the target, is GOSIA's integrated Rutherford on the same input; and
    it is 11 % (ring 2) and 16 % (ring 3) above the elastic integral, which is why a ratio must say which it
    divides by."""
    _, aw, el = integrals_q0[ring]
    assert aw == pytest.approx(GOSIA_RUTHERFORD[ring], rel=2e-3)
    assert aw / el == pytest.approx({"Ring2": 1.111, "Ring3": 1.158}[ring], abs=0.003)


@pytest.mark.parametrize("ring", RINGS)
def test_all_orders_against_gosia_with_matched_masses(all_orders, integrals_q0, ring):
    """Against GOSIA run on the app's input with the nuclear masses and without E1 polarisation, the app agrees
    within 0.5 % (0.2 % found), both divided by the elastic Rutherford integral (R_el, from the rates) and by
    GOSIA's own inelastic one (R_G)."""
    assert all_orders[0.0][ring] == pytest.approx(GOSIA_MATCHED[ring], rel=5e-3)
    c, aw, _ = integrals_q0[ring]
    assert c / aw == pytest.approx(GOSIA_MATCHED_RG[ring], rel=5e-3)


def test_thin_layers_keep_the_mid_layer_probability():
    """A layer the beam barely slows in keeps the one solution at its middle: the old results, exactly."""
    exp = Experiment.from_toml(SETUP.format(order="", extra="").replace('"1.03741 mg/cm2"', '"0.002 mg/cm2"'))
    ch, layers = excited_channel(exp)
    dp = depth_probability(exp, ch, layers)
    assert dp.thin
    mid = rates_module.excitation_probability(exp, ch, layers)
    theta = np.linspace(10, 179, 50)
    np.testing.assert_array_equal(dp.probability(theta), mid.probability(theta))
    np.testing.assert_array_equal(dp.at(111.0).probability(theta), mid.probability(theta))
    cfg, _ = generator_config(exp, pilot=False)
    assert len(cfg["layers"]) == 1


def test_coulex_ni58_moves_by_the_depth_average_only():
    """The 58Ni example (a 7 % energy loss) changes by the physics of the depth average, about +0.6 %: the same
    as averaging the exact first-order P over the depth by brute force."""
    exp = Experiment.example("coulex_ni58")
    ch, layers = excited_channel(exp)
    dp = depth_probability(exp, ch, layers)
    theta = np.array([100.0, 140.0, 170.0])
    e = ch.excitation
    zs = (np.arange(48) + 0.5) / 48 * dp.reached
    es = beam_energy_at(exp, 0, zs, layers)
    beam = rates_module.beam_ion(exp)
    w = np.array([Rutherford(beam, ch.nuclide.name, float(x)).cross_section_cm(90.0) for x in es])
    p = np.array([Coulex(beam, ch.nuclide.name, float(x), excite=e.excite, energy=e.energy_mev,
                         multipolarity=e.multipolarity, b_up=e.b_up_e2fm).probability(theta) for x in es])
    np.testing.assert_allclose(dp.probability(theta), (w[:, None] * p).sum(axis=0) / w.sum(), rtol=2e-4)
    mid = rates_module.excitation_probability(exp, ch, layers).probability(theta)
    assert np.all(np.abs(dp.probability(theta) / mid - 1) < 0.01)


def test_events_agree_with_the_rates(first_order):
    """The generator cuts the target into slices with P averaged over each: its excitation rates match the
    analytic ones, and the excitations sit at the front of the target, where the beam is fastest."""
    exp = allmond()
    cfg, used = generator_config(exp)
    assert len(cfg["layers"]) > 1 and sum(lay["thickness"] for lay in cfg["layers"]) == pytest.approx(1.03741, 1e-3)
    ev = simulate(exp, events=400_000, seed=11)
    rt = Rates(exp)
    exc = next(label for label in ev.channels if "excited" in label)
    for r in RINGS:
        got, err = ev.rate(r, counted=False, channel=exc)
        want = rt.rate(r, counted=False, what="excitations")
        assert abs(got - want) < 4 * err + 2e-3 * want
    c = ev.columns
    excited = c["channel"] == ev.channels.index(exc)
    elastic = ~excited
    mean_depth = np.average(c["depth"][excited], weights=c["weight"][excited])
    assert mean_depth < 0.9 * np.average(c["depth"][elastic], weights=c["weight"][elastic])
    assert np.all(c["depth"] <= 1.03741 * 1.001)


# -- Q_s for a single 0⁺ → 2⁺ transition -----------------------------------------------------------------------


def test_a_static_moment_keeps_the_coupled_equations_for_one_transition():
    assert cascade_excitation(allmond()) is None  # no scheme: the [reaction] state alone
    solved = cascade_excitation(allmond(0.1))
    assert solved is not None and has_static_moment(solved.scheme)
    assert type(solved).__name__ == "CoupledChannels"
    # First order has no reorientation: the one-γ-ray chain stays.
    assert cascade_excitation(allmond(0.1, order="first")) is None
    # Without the diagonal element a single transition keeps the one-γ-ray chain, as before.
    plain = Experiment.from_toml(SETUP.format(order="", extra=LEVELS.format(diagonal="")))
    assert cascade_excitation(plain) is None


@pytest.mark.parametrize("ring", RINGS)
def test_all_orders_against_gosia(all_orders, first_order, ring):
    """The {0⁺, 2⁺} scheme with Q = 0, to all orders through the depth, against GOSIA directly (+1.3 % and +1.2 %:
    GOSIA's mass numbers and E1 polarisation); and against first order: a few tenths of a percent lower, more so
    for ring 2 (larger CM angles, larger P)."""
    ours = all_orders[0.0][ring]
    assert ours == pytest.approx(GOSIA[ring], rel=GOSIA_TOLERANCE)
    assert 0.008 < ours / GOSIA[ring] - 1 < 0.018
    assert -0.010 < ours / first_order[ring] - 1 < -0.002
    if ring == "Ring2":
        assert all_orders[0.0]["Ring2"] / first_order["Ring2"] < all_orders[0.0]["Ring3"] / first_order["Ring3"]


@pytest.mark.parametrize("ring", RINGS)
def test_quadrupole_moment_changes_the_yield_as_in_gosia(all_orders, ring):
    """Q_s = ±0.1 eb changes the ratio by about ±2.7 % (ring 2) and ±1.9 % (ring 3), of GOSIA's sign: a positive
    Q_s raises the yield at these backward CM angles."""
    zero = all_orders[0.0][ring]
    up = 100 * (all_orders[0.1][ring] / zero - 1)
    down = 100 * (all_orders[-0.1][ring] / zero - 1)
    assert up == pytest.approx(GOSIA_Q_EFFECT[ring], abs=0.4)
    assert down == pytest.approx(-GOSIA_Q_EFFECT[ring], abs=0.4)


def test_order_first_ignores_the_static_moment(first_order):
    got = ratios(allmond(0.1, order="first"))
    for r in RINGS:
        assert got[r] == pytest.approx(first_order[r], rel=1e-9)


def test_gamma_events_follow_the_single_transition_scheme():
    """With Q_s set the [reaction] nucleus decays through its scheme: the events hold its one 2⁺ → 0⁺ line."""
    from physim.nuclear.gamma import transition_lines
    from physim.nuclear.gamma_events import simulate_gammas

    exp = allmond(0.1, gamma=True)
    lines = transition_lines(exp)
    assert [(ln["initial"], ln["final"]) for ln in lines] == [(1, 0)]
    g = simulate_gammas(exp, events=40_000, seed=2)
    assert len(g.columns["energy0"]) > 0
    np.testing.assert_allclose(np.unique(np.round(g.columns["energy0"], 4)), [1.17289], atol=2e-4)
