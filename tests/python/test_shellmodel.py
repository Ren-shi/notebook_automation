"""Shell-model results: KSHELL's logs read into numbers, and chosen values put into a level scheme.

The logs are the wave 0 fixtures in ``tests/reference/nuclear/shellmodel/kshell`` (²⁸Si with USDB, ⁴⁸Cr with KB3),
read as the bytes KSHELL wrote. The cross-checks are KSHELL's own ``collect_logs.py`` summary of the same logs and
the published numbers in the fixtures' README (Frycz et al., arXiv:2404.14506).
"""

import json
import math
import re
from pathlib import Path

import pytest

from physim.nuclear import Experiment, shellmodel as sm
from physim.nuclear.levels import Level, LevelScheme, MatrixElement, Transition, Value, quadrupole_factor

FIXTURES = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "shellmodel" / "kshell"
SI28 = FIXTURES / "si28_usdb" / "logs"
SI28_EFF = FIXTURES / "si28_usdb" / "logs_eff_1.31_0.46"
CR48 = FIXTURES / "cr48_kb3" / "log_cr48_m0.txt"


def raw(path: Path) -> str:
    """A log exactly as written: line endings untouched."""
    return path.read_bytes().decode("utf-8")


@pytest.fixture(scope="module")
def si28():
    return sm.read_folder(SI28)


@pytest.fixture(scope="module")
def si28_eff():
    """The 1.31/0.46 transit runs, with the states from the kshell logs they were run on."""
    return sm.read_logs([raw(p) for p in sorted(SI28.glob("log_Si28_usdb_j*p.txt"))],
                        [raw(p) for p in sorted(SI28_EFF.glob("log_tr_*.txt"))])


# -- kshell.exe logs -------------------------------------------------------------------------------------------------


def test_kshell_log_states_and_header():
    log = sm.parse_kshell_log(raw(SI28 / "log_Si28_usdb_j4p.txt"))
    assert [s.label for s in log] == ["2+1", "2+2", "2+3"]
    assert [s.energy_mev for s in log] == [-133.92904, -128.33707, -127.95966]
    assert log[1].ex_mev == pytest.approx(-128.33707 + 133.92904)
    assert all(s.t2 == 0 for s in log)
    first = log[0]
    assert (first.q_p, first.q_n, first.q_s) == (10.454, 10.454, 20.907)
    assert (first.l_p, first.l_n, first.s_p, first.s_n, first.mu) == (0.959, 0.959, 0.041, 0.041, 1.031)
    assert (log.dimension, log.mtot, log.parity) == (81122, 4, 1)
    assert (log.eff_charge, log.gl, log.gs) == ((1.5, 0.5), (1.0, 0.0), (5.585, -3.826))
    assert (log.fn_int, log.wave) == ("usdb.snt", "Si28_usdb_j4p.wav")
    assert (log.z_val, log.n_val, log.core_z, log.core_n) == (6, 6, 8, 8)
    assert log.orbits == ["0d3/2", "0d5/2", "1s1/2"]
    assert log.elapsed_s == pytest.approx(12.894)
    assert log.notes == []


def test_moments_in_the_kshell_log_are_charges_and_g_factors_times_the_parts():
    """<eQ> = e_p Q_p + e_n Q_n and μ = <gm> = g_l,p L_p + g_l,n L_n + g_s,p s_p + g_s,n s_n."""
    for name in ("log_Si28_usdb_j4p.txt", "log_Si28_usdb_j8p.txt"):
        log = sm.parse_kshell_log(raw(SI28 / name))
        (ep, en), (glp, gln), (gsp, gsn) = log.eff_charge, log.gl, log.gs
        for s in log:
            assert s.q_s == pytest.approx(ep * s.q_p + en * s.q_n, abs=2e-3)
            assert s.mu == pytest.approx(glp * s.l_p + gln * s.l_n + gsp * s.s_p + gsn * s.s_n, abs=2e-3)


def test_magnetic_moment_agrees_with_the_diagonal_m1_line_of_transit(si28):
    kshell = {s.label: s for name in ("log_Si28_usdb_j4p.txt", "log_Si28_usdb_j8p.txt")
              for s in sm.parse_kshell_log(raw(SI28 / name))}
    diagonal = [t for t in si28.transitions if t.kind == "M1" and t.diagonal]
    assert len(diagonal) == 5
    for t in diagonal:
        label = sm.label(*t.initial)
        assert t.moment == pytest.approx(kshell[label].mu, abs=6e-4)
        # μ = √(4π/3) ⟨J J 1 0|J J⟩ / √(2J+1) × Mred(M1), with ⟨J J 1 0|J J⟩ = √(J/(J+1))
        j = t.initial[0] / 2
        assert t.moment == pytest.approx(math.sqrt(4 * math.pi / 3) * math.sqrt(j / (j + 1)) / math.sqrt(2 * j + 1)
                                         * t.mred, abs=2e-4)


def test_a_log_that_is_not_j_projected_mixes_j_values():
    log = sm.parse_kshell_log(raw(CR48))
    assert len(log) == 32 and (log.dimension, log.mtot, log.parity) == (1963461, 0, 1)
    assert log[0].energy_mev == -32.95367 and log[0].label == "0+1"
    assert [s.label for s in log[:9]] == ["0+1", "2+1", "4+1", "6+1", "2+2", "2+3", "4+2", "3+1", "0+2"]
    assert log[1].q_s == -29.140 and log[1].q_p == -14.570
    assert all(s.mu is None for s in log)  # an M = 0 run prints no M1 parts
    assert log[0].q_s is None
    assert log.gs == (5.586, -3.826) and log.orbits == ["0f7/2", "1p3/2", "0f5/2", "1p1/2"]
    labels = [s.label for s in log]
    assert len(set(labels)) == 32
    for (j2, parity) in {(s.j2, s.parity) for s in log}:
        same = [s for s in log if (s.j2, s.parity) == (j2, parity)]
        assert [s.index for s in sorted(same, key=lambda s: s.energy_mev)] == list(range(1, len(same) + 1))
    # <JJ> = 6.00003 for the last state is close enough to J(J+1) = 6
    assert log.notes == []


def test_crlf_lf_and_mixed_line_endings_read_the_same():
    text = raw(SI28 / "log_Si28_usdb_tr_j4p_j8p.txt")
    assert "\r\n" in text
    lf = text.replace("\r\n", "\n")
    lines = lf.split("\n")
    mixed = "".join(line + ("\r\n" if i % 3 else "\n") for i, line in enumerate(lines))
    reference = sm.parse_transit_log(text)
    assert len(reference) > 0
    assert sm.parse_transit_log(lf) == reference == sm.parse_transit_log(mixed)
    k = raw(SI28 / "log_Si28_usdb_j8p.txt")
    k_mixed = "".join(line + ("\n" if i % 2 else "\r\n") for i, line in enumerate(k.replace("\r\n", "\n").split("\n")))
    assert sm.parse_kshell_log(k) == sm.parse_kshell_log(k_mixed)


# -- transit.exe logs ------------------------------------------------------------------------------------------------


def test_transit_lines_run_from_the_higher_state_to_the_lower():
    log = sm.parse_transit_log(raw(SI28 / "log_Si28_usdb_tr_j0p_j4p.txt"))
    assert (log.wave_l, log.wave_r, log.parity_l, log.parity_r) == ("Si28_usdb_j0p.wav", "Si28_usdb_j4p.wav", 1, 1)
    assert log.eff_charge == (1.5, 0.5) and log.kinds == ["E2"] and len(log) == 9
    # " 0(   1) -135.861 4(   1) -133.929   1.932  -22.3433  499.2211   99.8442": the 2+ (right) is higher
    t = log[0]
    assert (t.kind, t.initial, t.final, t.waves) == ("E2", (4, 1, 1), (0, 1, 1),
                                                     ("Si28_usdb_j4p.wav", "Si28_usdb_j0p.wav"))
    assert (t.mred, t.b_down, t.b_up) == (-22.3433, 99.8442, 499.2211)
    assert t.b_down == pytest.approx(t.mred**2 / 5, abs=1e-3)
    # " 0(   2) -131.024 4(   1) -133.929  -2.905    9.0754   82.3630   16.4726": the 0+ (left) is higher
    t = log[3]
    assert (t.initial, t.final, t.mred, t.b_down, t.b_up) == ((0, 1, 2), (4, 1, 1), 9.0754, 82.3630, 16.4726)
    assert t.moment is None


def test_transit_from_one_file_to_itself_keeps_each_pair_once_with_the_moment():
    log = sm.parse_transit_log(raw(SI28_EFF / "log_tr_j4p_j4p.txt"))
    assert log.eff_charge == (1.31, 0.46) and log.gl == (1.0, 0.0) and log.gs == (5.585, -3.826)
    e2 = [t for t in log if t.kind == "E2"]
    assert len(e2) == 6  # 3 diagonal + 3 pairs, each printed in both orders
    diag = {t.initial: t for t in e2 if t.diagonal}
    assert diag[(4, 1, 1)].moment == 18.5028 and diag[(4, 1, 1)].mred == 24.4122
    assert diag[(4, 1, 1)].moment == pytest.approx(quadrupole_factor(2) * 24.4122, abs=2e-3)
    assert log.notes == []


def test_turning_a_matrix_element_round_applies_the_phase_for_odd_spin_changes():
    """⟨a‖M‖b⟩ = (−1)^(J_a − J_b) ⟨b‖M‖a⟩: a 3+ → 2+ line printed with the higher state on the left."""
    text = "\n".join([
        "&INPUT", ' FN_LOAD_WAVE_L="x_j6p.wav    ",', ' FN_LOAD_WAVE_R="x_j4p.wav    ",', " /",
        " E2 transition  e^2 fm^4  eff_charge=  1.5000  0.5000 parity  1  1",
        " 2xJf      Ef      2xJi     Ei       Ex       Mred.    B(EM )->   B(EM)<-   Mom.",
        " 6(   1)  -10.000 4(   1)  -11.000  -1.000    3.0000    1.2857    1.8000    0.0000",
        " 6(   1)  -10.000 4(   2)   -9.000   1.000    2.0000    0.5714    0.8000    0.0000",
        "", " M1 transition  mu_N^2  gl,gs=  1.0000  0.0000  5.5850 -3.8260 parity  1  1",
        " 2xJf      Ef      2xJi     Ei       Ex       Mred.    B(EM )->   B(EM)<-   Mom.",
        " 6(   1)  -10.000 4(   1)  -11.000  -1.000    0.5000    0.0357    0.0500    0.0000", ""])
    log = sm.parse_transit_log(text)
    e2_down, e2_up, m1 = log
    assert (e2_down.initial, e2_down.final, e2_down.mred) == ((6, 1, 1), (4, 1, 1), -3.0)
    assert (e2_down.b_down, e2_down.b_up) == (1.2857, 1.8)  # B(3+ → 2+) = 9/7
    assert (e2_up.initial, e2_up.final, e2_up.mred) == ((4, 1, 2), (6, 1, 1), 2.0)  # already ⟨f‖M‖i⟩
    assert (m1.kind, m1.mred) == ("M1", -0.5)
    assert log.gl == (1.0, 0.0)


def row(j2l, nl, el, j2r, nr, er, *columns) -> str:
    """One table line in transit.exe's format (2(i2,"(",i4,")",f9.3), f8.3, 5f10.4); a column given as None is
    printed as Fortran prints an overflow."""
    out = f"{j2l:2d}({nl:4d}){el:9.3f}{j2r:2d}({nr:4d}){er:9.3f}{er - el:8.3f}"
    return out + "".join("*" * 10 if x is None else f"{x:10.4f}" for x in columns)


def test_overflowed_columns_are_worked_out_or_left_out():
    assert row(0, 1, -135.861, 4, 1, -133.929, -22.3433, 499.2211, 99.8442, 0.0) == \
        " 0(   1) -135.861 4(   1) -133.929   1.932  -22.3433  499.2211   99.8442    0.0000"
    text = "\n".join([
        " E2 transition  e^2 fm^4  eff_charge=  1.5000  0.5000 parity  1  1",
        " 2xJf      Ef      2xJi     Ei       Ex       Mred.    B(EM )->   B(EM)<-   Mom.",
        row(0, 1, -135.861, 4, 1, -133.929, -400.0, None, None, 0.0),
        row(0, 1, -135.861, 4, 2, -133.0, None, None, None, 0.0), ""])
    log = sm.parse_transit_log(text)
    assert len(log) == 1 and log[0].b_down == pytest.approx(160000 / 5) and log[0].b_up == pytest.approx(160000)
    assert any("overflowed" in n for n in log.notes)


# -- one result from a folder ---------------------------------------------------------------------------------------


def summary():
    """KSHELL's collect_logs.py summary: the level table and the B(E2) and B(M1) tables."""
    levels, tables, kind = [], {"E2": [], "M1": []}, None
    for line in raw(SI28 / "summary_Si28_usdb.txt").splitlines():
        if m := re.match(r"\s*\d+\s+(\d+)\s+([+-])\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+log", line):
            levels.append((2 * int(m.group(1)), 1 if m.group(2) == "+" else -1, int(m.group(3)), int(m.group(4)),
                           float(m.group(5)), float(m.group(6))))
        elif m := re.match(r"B\((E2|M1)\)", line):
            kind = m.group(1)
        elif kind and (m := re.match(r"\s*(\d+)([+-])\(\s*(\d+)\)\s+(\S+)\s+(\d+)([+-])\(\s*(\d+)\)\s+(\S+)\s+(\S+)"
                                     r"\s+(\S+)\(\s*\S+\)\s+(\S+)\(", line)):
            g = m.groups()
            tables[kind].append(((2 * int(g[0]), 1 if g[1] == "+" else -1, int(g[2])),
                                 (2 * int(g[4]), 1 if g[5] == "+" else -1, int(g[6])),
                                 float(g[3]), float(g[7]), float(g[9]), float(g[10])))
    return levels, tables


def test_the_folder_reproduces_kshells_own_summary(si28):
    levels, tables = summary()
    assert len(levels) == len(si28.states) == 8
    for j2, parity, index, t2, energy, ex in levels:
        s = si28.state(j2, parity, index)
        assert s.energy_mev == pytest.approx(energy, abs=6e-4) and s.ex_mev == pytest.approx(ex, abs=6e-4)
        assert s.t2 == 2 * t2
    assert [s.label for s in si28.states] == ["0+1", "2+1", "4+1", "0+2", "4+2", "0+3", "2+2", "2+3"]
    assert len(tables["E2"]) == 19 and len(tables["M1"]) == 4
    for kind, rows in tables.items():
        for initial, final, ex_i, ex_f, b_down, b_up in rows:
            t = si28.transition(kind, initial, final)
            assert t is not None and (t.initial, t.final) == (initial, final), (kind, initial, final)
            assert si28.state(*initial).ex_mev == pytest.approx(ex_i, abs=1.5e-3)
            assert si28.state(*final).ex_mev == pytest.approx(ex_f, abs=1.5e-3)
            digits = 0.05 if kind == "E2" else 5e-4
            assert t.b_down == pytest.approx(b_down, abs=digits) and t.b_up == pytest.approx(b_up, abs=digits)


def test_the_folder_result_and_what_produced_it(si28):
    assert (si28.nuclide, si28.interaction, si28.kshell_version) == ("28Si", "usdb", "unknown")
    assert si28.space == "0d3/2 0d5/2 1s1/2 outside 16O"
    assert (si28.eff_charge, si28.gl, si28.gs) == ((1.5, 0.5), (1.0, 0.0), (5.585, -3.826))
    assert si28.dimensions == {(0, 1): 93710, (4, 1): 81122, (8, 1): 52175}
    assert si28.elapsed_s > 20 and si28.notes == []
    ground = si28.state(0, 1, 1)
    assert ground.energy_mev == pytest.approx(-135.861, abs=5e-4)  # the README's benchmark, −135.9 published
    # every transition: B↓ = Mred² / (2J_i + 1), B↑ = Mred² / (2J_f + 1), to KSHELL's printed digits
    for t in si28.transitions:
        assert t.b_down == pytest.approx(t.mred**2 / (t.initial[0] + 1), abs=2e-3 * abs(t.mred) + 1e-4)
        assert t.b_up == pytest.approx(t.mred**2 / (t.final[0] + 1), abs=2e-3 * abs(t.mred) + 1e-4)
        assert t.waves == ()
        assert t.diagonal or si28.state(*t.initial).energy_mev > si28.state(*t.final).energy_mev
    # moments from transit's diagonal lines (4 decimals)
    assert si28.state(4, 1, 1).q_s == 20.9071 and si28.state(4, 1, 1).mu == 1.0312


def test_the_readme_benchmarks_with_effective_charges_131_046(si28_eff):
    assert si28_eff.eff_charge == (1.31, 0.46)
    assert round(si28_eff.transition("E2", (4, 1, 1), (0, 1, 1)).b_down, 2) == 78.20
    assert round(si28_eff.transition("E2", (8, 1, 1), (4, 1, 1)).b_down, 2) == 110.02
    assert round(si28_eff.state(4, 1, 1).q_s, 2) == 18.50
    # the 4+ moments have no diagonal line in these runs: recomputed from the kshell parts with 1.31/0.46
    s = si28_eff.state(8, 1, 1)
    assert s.q_s == pytest.approx(1.31 * s.q_p + 0.46 * s.q_n)
    assert any("Q_s" in n and "worked out again" in n for n in si28_eff.notes)


def test_a_folder_of_transit_logs_alone_still_gives_the_states(si28_eff):
    alone = sm.read_folder(SI28_EFF)
    assert [s.label for s in alone.states] == [s.label for s in si28_eff.states]
    assert any("known only from the transit logs" in n for n in alone.notes)
    assert alone.nuclide == "28Si"
    for t in si28_eff.transitions:
        u = alone.transition(t.kind, t.initial, t.final)
        assert u is not None and u.mred == t.mred and u.b_down == t.b_down


def test_transit_logs_with_different_charges_are_refused():
    with pytest.raises(ValueError, match="different eff charges"):
        sm.read_logs([], [raw(SI28 / "log_Si28_usdb_tr_j0p_j4p.txt"), raw(SI28_EFF / "log_tr_j0p_j4p.txt")])
    with pytest.raises(ValueError, match="no kshell.exe or transit.exe logs"):
        sm.read_folder(FIXTURES.parent / "snt")


def test_result_survives_json(si28):
    request = sm.ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3), (8, 1, 2)))
    si28.request = request
    try:
        d = json.loads(json.dumps(si28.to_dict()))
        again = sm.ShellModelResult.from_dict(d)
    finally:
        si28.request = None
    assert again.request == request and again.request.key() == request.key()
    assert again.states == si28.states and again.transitions == si28.transitions
    assert again.dimensions == si28.dimensions and again.eff_charge == si28.eff_charge
    assert again.to_dict() == d


# -- the request ------------------------------------------------------------------------------------------------------


def test_request_key_is_stable_and_covers_what_changes_the_result():
    a = sm.ShellModelRequest("Si28", "USDB", [(4, 1, 3), (0, 1, 3)], ["E2", "M1"], (1.5, 0.5))
    b = sm.ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3)), ("M1", "E2"), (1.5, 0.50), threads=4,
                             allow_over_limit=True)
    assert a.nuclide == "28Si" and a.interaction == "usdb" and a.states == ((0, 1, 3), (4, 1, 3))
    assert a.key() == b.key() and len(a.key()) == 16
    assert a.key() == sm.ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3)), eff_charge=(1.5, 0.5)).key()
    assert a.key() != sm.ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 3)), eff_charge=(1.31, 0.46)).key()
    assert a.key() != sm.ShellModelRequest("28Si", "usdb", ((0, 1, 3), (4, 1, 2)), eff_charge=(1.5, 0.5)).key()
    assert a.key() != sm.ShellModelRequest("28Si", "usda", ((0, 1, 3), (4, 1, 3)), eff_charge=(1.5, 0.5)).key()
    assert a.key() != a.key(build="kshell 2019, gfortran 16.2")
    # pinned: a cache folder named by this key must not move between versions of physim
    assert sm.ShellModelRequest("28Si", "usdb", ((0, 1, 3),)).key() == "a87ebbe0a870961a"
    with pytest.raises(AttributeError):
        a.nuclide = "30Si"


@pytest.mark.parametrize("kwargs, message", [
    ({"nuclide": "28Xx"}, "not an element"),
    ({"interaction": " "}, "name the interaction"),
    ({"states": ((3, 1, 1),)}, "2J must be even"),
    ({"states": ((4, 0, 1),)}, "parity"),
    ({"states": ((4, 1, 0),)}, "at least one"),
    ({"states": ((4, 1, 1), (4, 1, 2))}, "asked for twice"),
    ({"transitions": ("E3",)}, "transitions are"),
    ({"eff_charge": 1.5}, "two numbers"),
])
def test_request_problems_are_named(kwargs, message):
    args = {"nuclide": "28Si", "interaction": "usdb", "states": ((0, 1, 1),)} | kwargs
    with pytest.raises(ValueError, match=message):
        sm.ShellModelRequest(**args)


def test_labels():
    assert sm.label(4, 1, 1) == "2+1" and sm.label(3, -1, 2) == "3/2-2" and sm.label(0, 1) == "0+"
    assert sm.parse_label("2+1") == (4, 1, 1) and sm.parse_label("3/2-") == (3, -1, None)
    with pytest.raises(ValueError):
        sm.parse_label("2/2+1")


# -- into a level scheme ----------------------------------------------------------------------------------------------


def si28_scheme() -> LevelScheme:
    """A small ²⁸Si scheme typed for the test, in the shape ENSDF gives (values near the evaluated ones): 0+, 2+,
    4+, 0+, 3+, with an ENSDF γ ray and a derived 2+ → 0+ matrix element."""
    def e(x):
        return Value(x, 0.1, "ensdf")

    s = LevelScheme("28Si", [Level(e(0.0), 0, 1, "0+"), Level(e(1779.0), 2, 1, "2+"), Level(e(4617.9), 4, 1, "4+"),
                             Level(e(4979.9), 0, 1, "0+"), Level(e(6276.2), 3, 1, "3+")],
                    [Transition(1, 0, e(1779.0), e(100.0), "E2"), Transition(2, 1, e(2838.9), e(100.0), "E2")],
                    reference="typed for the test")
    s.matrix_elements.append(MatrixElement(0, 1, "E2", Value(18.0, 0.5, "derived", "from B(E2) in ENSDF")))
    return s


def test_states_are_matched_by_spin_parity_then_order(si28):
    pairs = sm.match(si28, si28_scheme())
    got = {p.state.label: p.level for p in pairs}
    assert got == {"0+1": 0, "2+1": 1, "4+1": 2, "0+2": 3, "4+2": None, "0+3": None, "2+2": None, "2+3": None}
    two = next(p for p in pairs if p.state.label == "2+1")
    assert two.known_kev == 1779.0 and two.delta_kev == pytest.approx(1931.69 - 1779.0, abs=0.01)
    assert [p.state.label for p in pairs] == [s.label for s in sorted(si28.states, key=lambda s: s.ex_mev)]


def test_apply_fills_gaps_and_keeps_ensdf_and_derived_values(si28_eff):
    s = si28_scheme()
    changes = sm.apply(si28_eff, s)
    # ENSDF energies and the derived 2+ → 0+ matrix element are untouched
    assert [lv.energy.source for lv in s.levels] == ["ensdf"] * 5
    assert s.matrix_element(0, 1).value == 18.0 and s.matrix_element(0, 1).source == "derived"
    # the 4+ → 2+ matrix element is new, with KSHELL's sign; ⟨4+‖E2‖2+⟩ from ⟨2+‖E2‖4+⟩ (even ΔJ: same sign)
    v = s.matrix_element(1, 2)
    t = si28_eff.transition("E2", (8, 1, 1), (4, 1, 1))
    assert (v.value, v.source) == (t.mred, "shell model") and v.value < 0
    assert "usdb" in v.note and "e_p = 1.31 e, e_n = 0.46 e" in v.note
    assert s.b(2, 1) == pytest.approx(110.02, abs=0.01)
    # quadrupole moments become diagonal matrix elements
    assert s.quadrupole_moment(1) == pytest.approx(18.5028, abs=1e-4)
    assert s.quadrupole_moment(2) == pytest.approx(si28_eff.state(8, 1, 1).q_s)
    assert s.matrix_element(1, 1).source == "shell model"
    # 0+2 → 2+1 (no E2 between the two 0+)
    assert s.matrix_element(1, 3).value == pytest.approx(si28_eff.transition("E2", (0, 1, 2), (4, 1, 1)).mred)
    assert s.matrix_element(0, 3) is None
    # the 3+ has no shell-model partner in this result
    assert all(m.a != 4 and m.b != 4 for m in s.matrix_elements)
    assert any(c.startswith("kept 3 ensdf energies, 1 derived E2 matrix element (") for c in changes)
    assert sm.apply(si28_eff, s)[:-1] == []  # applying again changes nothing


def test_apply_replaces_protected_values_only_when_asked(si28_eff):
    s = si28_scheme()
    s.set_matrix_element(1, 2, "E2", "-0.3 eb")  # a user value
    sm.apply(si28_eff, s, ["E2"])
    assert s.matrix_element(1, 2).source == "user" and s.matrix_element(0, 1).source == "derived"
    changes = sm.apply(si28_eff, s, ["E2"], replace=True)
    t = si28_eff.transition("E2", (4, 1, 1), (0, 1, 1))
    assert s.matrix_element(0, 1).value == t.mred and s.matrix_element(0, 1).source == "shell model"
    assert s.matrix_element(1, 2).source == "shell model"
    assert any("was +18" in c and "derived" in c for c in changes)
    sm.apply(si28_eff, s, ["energy"], replace=True)
    assert s.levels[1].energy.value == pytest.approx(1931.69, abs=0.01) and s.levels[1].energy.source == "shell model"
    assert s.levels[0].energy.value == 0.0 and s.levels[0].energy.source == "ensdf"


def test_apply_to_chosen_states_and_quantities(si28_eff):
    s = si28_scheme()
    changes = sm.apply(si28_eff, s, ["Q", "2+1"])
    assert [(m.a, m.b) for m in s.matrix_elements if m.value.source == "shell model"] == [(1, 1)]
    s = si28_scheme()
    sm.apply(si28_eff, s, ["E2", "4+1"])
    assert {(m.a, m.b) for m in s.matrix_elements if m.value.source == "shell model"} == {(1, 2)}  # 0+ to 4+ is not E2
    assert any("B(M1)" in c for c in sm.apply(si28_eff, si28_scheme(), ["M1"]))
    assert "B(M1)" not in " ".join(changes)
    with pytest.raises(ValueError, match="cannot apply 'banana'"):
        sm.apply(si28_eff, s, ["banana"])
    with pytest.raises(ValueError, match="is for 28Si, the level scheme for 30Si"):
        sm.apply(si28_eff, LevelScheme("30Si", [Level(Value(0.0), 0, 1)]))


def test_unmatched_states_become_new_levels_in_order(si28_eff):
    s = si28_scheme()
    changes = sm.apply(si28_eff, s, ["new"])
    assert [lv.energy.value for lv in s.levels] == pytest.approx(
        [0.0, 1779.0, 4617.9, 4979.9, 6276.2, 7004.95, 7326.75, 7523.66, 7901.07], abs=0.01)
    assert [lv.label for lv in s.levels] == ["0+", "2+", "4+", "0+", "3+", "4+", "0+", "2+", "2+"]
    assert [lv.energy.source for lv in s.levels[5:]] == ["shell model"] * 4
    assert sum(c.startswith("new level") for c in changes) == 4
    # the existing transitions and matrix elements still point at the same levels
    assert [(t.initial, t.final) for t in s.transitions] == [(1, 0), (2, 1)]
    assert s.matrix_element(0, 1).source == "derived"
    # the new levels are joined by shell-model matrix elements, so Coulomb excitation reaches them
    assert s.matrix_element(1, 7).value == pytest.approx(si28_eff.transition("E2", (4, 1, 2), (4, 1, 1)).mred)
    assert {5, 6, 7, 8} <= s.connected()
    # naming one unmatched state adds that one only
    one = si28_scheme()
    sm.apply(si28_eff, one, ["2+2"])
    assert len(one.levels) == 6 and one.levels[5].label == "2+"


def test_signs_follow_the_level_order_for_odd_spin_changes():
    """A level scheme keeps ⟨higher level‖M‖lower level⟩ and the result ⟨final‖M‖initial⟩ = ⟨lower‖M‖higher⟩, so an
    odd change of spin turns the sign, ⟨a‖M‖b⟩ = (−1)^(J_a − J_b) ⟨b‖M‖a⟩, unless the scheme has the two levels the
    other way round from the shell model."""
    states = [sm.ShellState(0, 1, 1, -10.0, 0.0), sm.ShellState(6, 1, 1, -9.0, 1.0), sm.ShellState(4, 1, 1, -8.8, 1.2)]
    transitions = [sm.ShellTransition("E2", (4, 1, 1), (6, 1, 1), 3.0, 9 / 5, 9 / 7),
                   sm.ShellTransition("E2", (4, 1, 1), (0, 1, 1), 4.0, 16 / 5, 16.0)]
    result = sm.ShellModelResult(None, "", "test", "unknown", (1.5, 0.5), None, None, states, transitions,
                                 nuclide="28Si")
    # as in the shell model: 3+ (level 1) below 2+ (level 2); stored ⟨2+‖E2‖3+⟩ = (−1)^(2−3) ⟨3+‖E2‖2+⟩
    same = LevelScheme("28Si", [Level(Value(0.0), 0, 1, "0+"), Level(Value(900.0), 3, 1, "3+"),
                                Level(Value(1200.0), 2, 1, "2+")])
    sm.apply(result, same, ["E2"])
    assert same.matrix_element(1, 2).value == -3.0
    assert same.matrix_element(0, 2).value == 4.0  # 2+ → 0+: even change, same sign
    assert same.b(2, 1) == pytest.approx(9 / 5)  # B(E2; 2+ → 3+)
    # the other way round: 2+ (level 1) below 3+ (level 2); stored ⟨3+‖E2‖2+⟩ is the result's value as it is
    swapped = LevelScheme("28Si", [Level(Value(0.0), 0, 1, "0+"), Level(Value(900.0), 2, 1, "2+"),
                                   Level(Value(1200.0), 3, 1, "3+")])
    sm.apply(result, swapped, ["E2"])
    assert swapped.matrix_element(1, 2).value == 3.0
    assert swapped.matrix_element(0, 1).value == 4.0


def test_shell_model_values_survive_the_setup_file(si28_eff):
    exp = Experiment.example("coulex_ni58")
    exp.target.material = "28Si"
    s = si28_scheme()
    sm.apply(si28_eff, s, ["all", "new"])
    exp.levels["target"] = s
    again = Experiment.from_toml(exp.to_toml())
    t = again.levels["target"]
    assert len(t.levels) == 9 and t.levels[6].energy.source == "shell model"
    assert t.levels[6].energy.note == s.levels[6].energy.note
    assert "KSHELL shell model, usdb" in t.levels[6].energy.note
    v, w = t.matrix_element(1, 2), s.matrix_element(1, 2)
    assert (v.value, v.source, v.note) == (pytest.approx(w.value), "shell model", w.note)
    assert t.matrix_element(0, 1).source == "derived"
    assert again.to_toml() == exp.to_toml()


def test_derive_keeps_shell_model_values(si28_eff):
    s = si28_scheme()
    s.transitions[0].weisskopf["E2"] = Value(13.0, 1.0, "ensdf")
    sm.apply(si28_eff, s, ["E2"], replace=True)
    s.derive()
    assert s.matrix_element(0, 1).source == "shell model"
    assert s.matrix_element(1, 2).source == "shell model"


def test_a_scheme_from_the_shell_model_alone(si28):
    s = sm.to_level_scheme(si28, max_energy_kev=7400)
    assert [lv.label for lv in s.levels] == ["0+", "2+", "4+", "0+", "4+", "0+"]
    assert s.levels[0].energy.value == 0.0 and {lv.energy.source for lv in s.levels} == {"shell model"}
    assert s.b(1, 0) == pytest.approx(99.8442, abs=1e-3)
    assert s.quadrupole_moment(1) == pytest.approx(20.9071, abs=1e-3)
    problems = []
    again = LevelScheme.from_dict(s.to_dict(), "levels.target", problems)
    assert problems == [] and again.matrix_element(1, 2).value == s.matrix_element(1, 2).value
