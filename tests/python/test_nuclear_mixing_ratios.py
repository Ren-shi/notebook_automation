"""ENSDF's mixing ratios and multipolarities as physim reads them. δ is the amplitude of the higher multipole over the
lower (E2/M1, L + 1 over L) with Krane and Steffen's sign, whatever order the multipolarity is written in, and E0 (no γ
ray, only conversion) is left out of the γ-ray multipoles: "E2+M1+E0" is an E2/M1 γ ray. The bands were fixed before
the change:

1. δ gives the same result whatever order the multipolarity is written in, for every caller;
2. A₂ and A₄ of a mixed 2 → 2 and 4 → 4 line against Krane and Steffen's F coefficients as exact surds, to 1e-6, for
   both signs of δ;
3. the E2 share of every mixed line of the templates and of ENSDF's ⁶²Ni, ⁹⁸Ru, ¹⁹⁴Pt and ¹⁵²Sm is δ²/(1+δ²);
4. B(E2) from the lifetime of a mixed line agrees with the lifetime analysis's to 1e-9."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from physim.nuclear import Experiment, angular, ensdf, lifetimes, nucleus, results
from physim.nuclear.correlation import Emission
from physim.nuclear.deorientation import Parameters
from physim.nuclear.experiment import example_names
from physim.nuclear.gamma import Faded
from physim.nuclear.levels import (LN2, Level, LevelScheme, MatrixElement, Transition, Value, components, fractions,
                                   leading, rate_per_b)
from physim.nuclear.lifetime_analysis import b_from_result
from physim.nuclear.orientation import Excitation, tensors_of
from physim.nuclear.residues import _multipoles

needs_ensdf = pytest.mark.skipif(not ensdf.available(), reason="no local copy of ENSDF (scripts/fetch_ensdf.py)")

#: One M1+E2 line written every way ENSDF writes it.
SPELLINGS = ["M1+E2", "E2+M1", "E2(+M1)", "(E2(+M1))", "(M1+E2)", "[M1+E2]", "M1+E2+E0", "E2+M1+E0", "E0+M1+E2"]
DELTAS = [3.1, -3.1, 0.3, -0.3, 15.0, -9.3]


def v(x: float) -> Value:
    return Value(x, source="user")


def scheme(spelling: str, delta, j_upper: float = 2.0) -> LevelScheme:
    """A J → J line mixed as written, from the level :func:`upper` to :func:`lower`. J = 2: 0⁺, 2⁺₁ at 1000 keV and
    2⁺₂ at 2000 keV (with a half-life), which decays to the 2⁺₁ (mixed) and to the ground state (E2). J = 4: 0⁺, 2⁺ at
    1000 keV, 4⁺₁ at 1500 keV and 4⁺₂ at 2000 keV (with a half-life), whose one line goes to the 4⁺₁."""
    s = LevelScheme("62Ni")
    mixing = None if delta is None else v(delta)
    if j_upper == 2.0:
        s.levels = [Level(v(0.0), 0.0, 1, "0+"), Level(v(1000.0), 2.0, 1, "2+", v(1.4e-12)),
                    Level(v(2000.0), 2.0, 1, "2+", v(0.5e-12))]
        s.transitions = [Transition(1, 0, v(1000.0), v(100.0), "E2", conversion=v(0.001)),
                         Transition(2, 1, v(1000.0), v(60.0), spelling, mixing, v(0.002)),
                         Transition(2, 0, v(2000.0), v(40.0), "E2", conversion=v(0.0005))]
        s.matrix_elements = [MatrixElement(0, 1, "E2", v(40.0)), MatrixElement(1, 2, "E2", v(30.0)),
                             MatrixElement(0, 2, "E2", v(5.0))]
    else:
        s.levels = [Level(v(0.0), 0.0, 1, "0+"), Level(v(1000.0), 2.0, 1, "2+", v(1.4e-12)),
                    Level(v(1500.0), 4.0, 1, "4+", v(3e-12)), Level(v(2000.0), 4.0, 1, "4+", v(0.5e-12))]
        s.transitions = [Transition(1, 0, v(1000.0), v(100.0), "E2", conversion=v(0.001)),
                         Transition(2, 1, v(500.0), v(100.0), "E2", conversion=v(0.01)),
                         Transition(3, 2, v(500.0), v(60.0), spelling, mixing, v(0.002))]
        s.matrix_elements = [MatrixElement(0, 1, "E2", v(40.0)), MatrixElement(1, 2, "E2", v(60.0)),
                             MatrixElement(2, 3, "E2", v(30.0))]
    return s


def upper(j: float = 2.0) -> int:
    return 2 if j == 2.0 else 3


def lower(j: float = 2.0) -> int:
    return upper(j) - 1



def excitation(s: LevelScheme) -> Excitation:
    return Excitation("12C", s.nuclide, 30.0, s, excite="target")


# -- band 1: the order written does not matter, for every caller ----------------------------------------------------


def test_components_are_the_gamma_ray_multipoles_lower_first():
    for text in SPELLINGS:
        assert components(text) == ["M1", "E2"], text
    assert components("M2+E1") == components("E1+M2") == ["E1", "M2"]
    assert components("E2+E0") == components("E0+E2") == ["E2"]
    assert components("E0") == ["E0"] and components("(E2)") == ["E2"]
    for text in ("D+Q", "M1,E2", "", "E1+M2+E3", "M1+M1"):
        assert components(text) is None, text


@pytest.mark.parametrize("delta", DELTAS)
def test_fractions_give_the_higher_multipole_delta_squared_over_one_plus_delta_squared(delta):
    want = {"M1": 1 / (1 + delta**2), "E2": delta**2 / (1 + delta**2)}
    for text in SPELLINGS:
        assert fractions(text, delta) == pytest.approx(want, rel=1e-12, abs=0), text
    assert fractions("E1+M2", delta) == pytest.approx({"E1": 1 / (1 + delta**2), "M2": delta**2 / (1 + delta**2)})
    assert fractions("E2+M1", None) == {}


@pytest.mark.parametrize("delta", [3.1, -0.3])
def test_b_from_a_half_life_reads_delta_the_same_whatever_the_order(delta):
    """LevelScheme's B from the level's half-life and its branching (``derive``)."""
    values = set()
    for text in SPELLINGS:
        s = scheme(text, delta)
        s.derive()
        b_e2 = s._b_down(s.transitions[1], "E2")[0]
        b_m1 = s._b_down(s.transitions[1], "M1")[0]
        total = 60.0 * 1.002 + 40.0 * 1.0005
        rate = LN2 / 0.5e-12 * 60.0 / total
        assert b_e2 == pytest.approx(rate * delta**2 / (1 + delta**2) / rate_per_b(2, 1000.0), rel=1e-12), text
        assert b_m1 == pytest.approx(rate / (1 + delta**2) / rate_per_b(1, 1000.0), rel=1e-12), text
        values.add((round(b_e2, 9), round(b_m1, 12), round(s.b(2, 1, "E2"), 9)))
    assert len(values) == 1


@pytest.mark.parametrize("delta", [3.1, -0.3])
def test_tau_from_a_matrix_element_reads_delta_the_same_whatever_the_order(delta):
    """lifetimes: τ of the 4⁺₂ level from the B(E2) of its one, mixed, line to the 4⁺₁."""
    taus = set()
    for text in SPELLINGS:
        s = scheme(text, delta, j_upper=4.0)
        tau, how = lifetimes._taus_from_matrix_elements(s, 3, s.branches(3))
        b = s.b(3, 2, "E2")
        want = 1e12 / (lifetimes.rate_per_unit_b("E2", 500.0) * b / (delta**2 / (1 + delta**2)) * 1.002)
        assert tau == pytest.approx(want, rel=1e-12), text
        taus.add(round(tau, 12))
    assert len(taus) == 1


@pytest.mark.parametrize("delta", DELTAS)
def test_the_decay_reads_delta_the_same_whatever_the_order(delta):
    """orientation.Excitation: L < L′, δ as given, its A_k and its γ-ray share; residues._multipoles the same."""
    seen = set()
    for text in SPELLINGS:
        s = scheme(text, delta)
        ex = excitation(s)
        d = ex.decay(2, 1)
        assert (d.l1, d.l2, d.mixing) == (1, 2, delta), text
        a = tuple(round(ex.distribution_coefficient(2, 1, k), 12) for k in (2, 4))
        seen.add((d.branch, d.gamma, a))
        assert _multipoles(s.transitions[1], 2.0, 2.0) == (1, 2, delta), text
    assert len(seen) == 1
    # Without δ the shares are not known: the line is taken as the multipole written first (ENSDF writes the dominant
    # one first, "(+M1)" a possible admixture), E0 left out.
    for text in SPELLINGS:
        first = 2 if text.lstrip("([").startswith("E2") or text.startswith("E0+E2") else 1
        s = scheme(text, None)
        d = excitation(s).decay(2, 1)
        assert (d.l1, d.l2, d.mixing) == (first, 0, 0.0), text
        assert _multipoles(s.transitions[1], 2.0, 2.0) == (first, 0, 0.0), text
        assert fractions(text, None) == {}
    assert leading("E2+M1+E0") == leading("(E2(+M1))") == "E2" and leading("E0+M1+E2") == leading("M1+E2") == "M1"


def test_the_deorientation_hands_on_the_same_orientation_whatever_the_order():
    """gamma.Faded: the unobserved feeding's U-coefficient shares come from the decay's multipoles and δ."""
    feeds = []
    for text in SPELLINGS:
        ex = excitation(scheme(text, 3.1))
        f = Faded(ex, 1, [0.0, 2.0, 0.7], Parameters())
        feeds.append({k: m.copy() for k, m in f._feeds.items()})
    for other in feeds[1:]:
        assert all(np.array_equal(other[k], feeds[0][k]) for k in feeds[0])
    # And it is the shares δ²/(1+δ²) of E2 and 1/(1+δ²) of M1.
    u = {lam: angular.u_coefficient(2, 2.0, 2.0, lam) for lam in (1, 2)}
    share = 60.0 * 1.002 / (60.0 * 1.002 + 40.0 * 1.0005)
    want = share * (u[1] / (1 + 3.1**2) + u[2] * 3.1**2 / (1 + 3.1**2))
    ex = excitation(scheme("E2+M1+E0", 3.1))
    f = Faded(ex, 1, [0.0, 2.0, 0.7], Parameters())
    a, b = f.above.index(2), f.above.index(1)
    assert f._feeds[2][b, a] == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("delta", [3.1, -0.3])
def test_b_from_a_measured_lifetime_reads_delta_the_same_whatever_the_order(delta):
    """lifetime_analysis.b_from_result."""
    bs = set()
    for text in SPELLINGS:
        b, err = b_from_result(SimpleNamespace(tau_ps=0.7, total_ps=0.07), scheme(text, delta), 2, final=1)
        bs.add(round(b, 9))
    assert len(bs) == 1


def test_the_conversion_between_b_and_tau_reads_delta_the_same_whatever_the_order(monkeypatch):
    """nucleus.decay: the inputs of a B ↔ τ conversion through a level's one γ ray (here a mixed 4⁺₂ → 4⁺₁ line)."""
    seen = []
    for text in SPELLINGS + ["M1+E2+E0"]:
        s = scheme(text, 3.1, j_upper=4.0)
        keys = tuple(results.level_key("62Ni", lv.energy.value, lv.spin, lv.parity) for lv in s.levels)
        data = SimpleNamespace(levels=[SimpleNamespace(dropped=0) for _ in s.levels])
        monkeypatch.setattr(nucleus, "adopted_scheme", lambda nuclide, found=(data, s, keys): found)
        key, inputs, why = nucleus.decay("62Ni", keys[3])
        assert why == "" and key.multipolarity == "E2" and inputs["mixing"] == 3.1, text
        seen.append((key.text, tuple(sorted(inputs.items()))))
    assert len(set(seen)) == 1


# -- band 2: A₂ and A₄ against Krane and Steffen ---------------------------------------------------------------------

#: F_k(L L′ J J) for J → J from Krane and Steffen's formula as exact surds (exact Racah sums; their tables give
#: F₂(1 1 2 2) = −0.4183 and F₂(1 2 2 2) = −0.6124). F₄ of a dipole is zero.
F = {(2, 2): {(2, 1, 1): -math.sqrt(7 / 40), (2, 1, 2): -math.sqrt(3 / 8), (2, 2, 2): math.sqrt(45 / 2744),
              (4, 2, 2): -math.sqrt(32 / 343)},
     (4, 4): {(2, 1, 1): -math.sqrt(77 / 400), (2, 1, 2): -math.sqrt(9 / 80), (2, 2, 2): math.sqrt(4225 / 60368),
              (4, 2, 2): -math.sqrt(936 / 3773)}}
BOTH_SIGNS = [3.0, -3.0, 0.3, -0.3, 13.0, -13.0]


def closed_form(k: int, j: float, delta: float) -> float:
    f = F[(int(j), int(j))]
    return (f.get((k, 1, 1), 0.0) + 2 * delta * f.get((k, 1, 2), 0.0) + delta**2 * f[(k, 2, 2)]) / (1 + delta**2)


def test_the_f_coefficients_are_krane_and_steffens():
    for (j, _), table in F.items():
        for (k, l1, l2), want in table.items():
            assert angular.f_coefficient(k, l1, l2, j, j) == pytest.approx(want, abs=1e-12)
            assert angular.f_coefficient(k, l2, l1, j, j) == pytest.approx(want, abs=1e-12)
    assert angular.f_coefficient(2, 2, 2, 0, 2) == pytest.approx(-0.5976, abs=5e-5)
    assert angular.f_coefficient(4, 2, 2, 0, 2) == pytest.approx(-1.0690, abs=5e-5)


@pytest.mark.parametrize("j", [2.0, 4.0])
@pytest.mark.parametrize("delta", BOTH_SIGNS)
def test_a2_a4_of_a_mixed_line_are_krane_and_steffens(j, delta):
    """Through the app's own path: the decay of a scheme written E2+M1 or M1+E2 with ENSDF's δ."""
    for text in ("E2+M1", "M1+E2", "E2+M1+E0"):
        ex = excitation(scheme(text, delta, j_upper=j))
        for k in (2, 4):
            a = ex.distribution_coefficient(upper(j), lower(j), k)
            assert a == pytest.approx(closed_form(k, j, delta), abs=1e-6), (text, k)
    a2 = closed_form(2, j, delta)
    # The sign of δ matters: the interference term 2δ F₂(1 2) flips with it.
    assert closed_form(2, j, -delta) != pytest.approx(a2, abs=1e-3)


@pytest.mark.parametrize("j", [2.0, 4.0])
@pytest.mark.parametrize("delta", BOTH_SIGNS)
def test_the_emitted_gamma_ray_has_krane_and_steffens_distribution(j, delta):
    """correlation.Emission (the event generator's): W from a level aligned along z against 1 + Σ B_k A_k P_k with the
    closed-form A_k, for either sign of δ."""
    d = int(round(2 * j + 1))
    p = np.exp(-0.5 * (np.arange(d) - j) ** 2 / 1.1**2)
    rho = np.diag(p / p.sum()).astype(complex)
    t = tensors_of(rho, j)
    theta = np.radians(np.linspace(0.0, 180.0, 37))
    dirs = np.c_[np.sin(theta), np.zeros_like(theta), np.cos(theta)]
    w = Emission(j, j, 1, 2, delta).weights(np.repeat(rho[None], len(dirs), axis=0), dirs)
    want = np.ones_like(theta)
    for k in (2, 4):
        legendre = np.polynomial.legendre.legval(np.cos(theta), [0] * k + [1])
        want += closed_form(k, j, delta) * t[(k, 0)].real / t[(0, 0)].real * legendre
    assert np.allclose(w, want, atol=1e-6)


# -- band 3: the E2 share of every mixed line ------------------------------------------------------------------------


def e2_share(ex, t) -> float:
    d = ex.decay(t.initial, t.final)
    shares = {d.l1: 1.0} if not d.l2 else {d.l1: 1 / (1 + d.mixing**2), d.l2: d.mixing**2 / (1 + d.mixing**2)}
    return shares.get(2, 0.0)


def mixed_lines(s: LevelScheme) -> list:
    return [t for t in s.transitions if t.mixing is not None and components(t.multipolarity or "") is not None
            and len(components(t.multipolarity)) == 2]


def check_shares(s: LevelScheme) -> int:
    ex = excitation(s)
    lines = mixed_lines(s)
    for t in lines:
        d = t.mixing.value
        lower, higher = components(t.multipolarity)
        want = d * d / (1 + d * d)
        assert fractions(t.multipolarity, d)[higher] == pytest.approx(want, rel=1e-12), t
        if higher == "E2":
            assert e2_share(ex, t) == pytest.approx(want, rel=1e-12), t
    return len(lines)


def test_the_e2_share_of_every_mixed_line_of_the_templates():
    counted = {}
    for name in example_names():
        for role, s in Experiment.example(name).levels.items():
            counted[name, role] = check_shares(s)
    # The templates' mixed lines: 62Ni's and 194Pt's 2+2 -> 2+1 ("M1+E2+E0", "E2+M1+E0") and 194Pt's 4+2 -> 4+1
    # ("(E2(+M1))") in both Pt templates; 143Sm's seven M1+E2 lines.
    assert counted["coulex_ni62_pt194", "beam"] == 1 and counted["coulex_ni62_pt194", "target"] == 2
    assert counted["coulex_ni62_pt194_clovers", "beam"] == 1 and counted["coulex_ni62_pt194_clovers", "target"] == 2
    assert counted["fusion_sm143_dsam", "143Sm"] == 7


@needs_ensdf
@pytest.mark.parametrize("nuclide, at_least", [("62Ni", 1), ("98Ru", 1), ("194Pt", 2), ("152Sm", 8)])
def test_the_e2_share_of_every_mixed_line_of_ensdf(nuclide, at_least):
    assert check_shares(LevelScheme.from_ensdf(nuclide, max_energy_kev=3000)) >= at_least


@needs_ensdf
def test_ensdfs_two_plus_two_lines_are_mostly_e2():
    """What the old reading got wrong: ⁹⁸Ru's 2⁺₂ → 2⁺₁ "E2+M1" (δ = 13) was 99.4 % M1, ¹⁹⁴Pt's 4⁺₂ → 4⁺₁ "(E2(+M1))"
    (δ = 3) 90 % M1, and ⁶²Ni's "M1+E2+E0" pure M1."""
    for nuclide, energy, delta in (("98Ru", 1414.4, 13.0), ("194Pt", 1229.5, 3.0), ("62Ni", 2301.8, 3.1)):
        s = LevelScheme.from_ensdf(nuclide, max_energy_kev=3000)
        t = next(t for t in mixed_lines(s) if abs(s.levels[t.initial].energy.value - energy) < 0.5)
        assert t.mixing.value == delta
        assert e2_share(excitation(s), t) == pytest.approx(delta**2 / (1 + delta**2), rel=1e-12)


# -- band 4: B(E2) from a lifetime -----------------------------------------------------------------------------------


@pytest.mark.parametrize("text", SPELLINGS)
@pytest.mark.parametrize("delta", [3.1, -0.3, 15.0])
def test_b_e2_from_a_lifetime_agrees_with_the_lifetime_analysis(text, delta):
    """The level scheme's B(E2) from the half-life (as for ENSDF) and the lifetime analysis's from τ = T½/ln 2 of the
    same mixed line, with its branching and conversion."""
    s = scheme(text, delta)
    mine = s._b_down(s.transitions[1], "E2")[0]
    tau_ps = 0.5e-12 / LN2 * 1e12
    theirs, _ = b_from_result(SimpleNamespace(tau_ps=tau_ps, total_ps=0.1 * tau_ps), s, 2, final=1)
    assert mine == pytest.approx(theirs, rel=1e-9)
