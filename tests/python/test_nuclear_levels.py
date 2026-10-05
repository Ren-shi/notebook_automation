"""Level schemes and the ENSDF reader (backlog item 50).

The ENSDF text here is typed by hand in ENSDF's card layout, with values close to the evaluated ones for ¹⁹⁴Pt. It
is not a copy of the evaluation. Tests marked ``needs_ensdf`` read the real local copy and are skipped without it.
"""

import math

import pytest

from physim.nuclear import Experiment, SetupError, ensdf
from physim.nuclear.levels import (LevelScheme, components, fractions, parse_matrix_element, quadrupole_factor,
                                   rate_per_b, weisskopf_unit)

needs_ensdf = pytest.mark.skipif(not ensdf.available(), reason="no local copy of ENSDF (scripts/fetch_ensdf.py)")


def card(nucid, kind, cols=(), cont=" ", comment=" "):
    """One 80-column record: ``cols`` are (first column counting from 1, text)."""
    line = list(f"{nucid:<5}{cont}{comment}{kind} ".ljust(80))
    for start, text in cols:
        line[start - 1:start - 1 + len(text)] = text
    return "".join(line)


def level(e, de="", j="", t="", dt=""):
    return card("194PT", "L", [(10, e), (20, de), (22, j), (40, t), (50, dt)])


def gamma(e, de="", ri="", m="", mr="", cc=""):
    return card("194PT", "G", [(10, e), (20, de), (22, ri), (32, m), (42, mr), (56, cc)])


PT194 = [
    card("194OS", " ", [(10, "ADOPTED LEVELS, GAMMAS")]),
    card("194OS", "L", [(10, "0.0"), (22, "0+")]),
    card("194OS", "L", [(10, "218.5"), (22, "2+")]),
    "",
    card("194PT", " ", [(10, "ADOPTED LEVELS, GAMMAS"), (75, "202109")]),
    card("194PT", "H", [(10, "TYP=FUL$AUT=A. Nother$CIT=NDS 177, 1 (2021)$")]),
    card("194PT", "L", [(10, "A comment that is not data")], comment="c"),
    level("0.0", j="0+", t="STABLE"),
    level("328.464", "12", "2+", "41.7 PS", "4"),
    card("194PT", "L", [(10, "MOME2=+0.48 14")], cont="2"),
    gamma("328.464", "12", "100", "E2", cc="0.0635"),
    level("622.004", "20", "2+", "35 PS", "4"),
    gamma("293.541", "14", "100", "M1+E2", "+14", "0.105"),
    card("194PT", "G", [(10, "BE2W=89 11$BM1W=7.7E-5 14")], cont="2"),
    gamma("622.00", "3", "9.3", "E2", cc="0.0126"),
    card("194PT", "G", [(10, "BE2W=0.29 4")], cont="2"),
    level("811.25", "3", "4+", "3.7 PS", "3"),
    card("194PT", "L", [(10, "MOMM1=1.1 2$MOME2=0.5 2 (2016St14)")], cont="2"),
    gamma("482.78", "3", "100", "E2", cc="0.0228"),
    card("194PT", "G", [(10, "FL=328.46+X")], cont="2"),
    level("1229.5+X", j="(5-)"),
    gamma("300.0", ri="100", m="E1"),
    level("1432.5", "1", "(3-)"),
    gamma("1104.0", "1", "100", "E1"),
    gamma("555.5", "1", "40", "D+Q"),
    "",
    card("194PT", " ", [(10, "194IR B- DECAY")]),
    card("194PT", "L", [(10, "999.0"), (22, "7+")]),
]


def pt194():
    return LevelScheme.from_adopted(ensdf.parse(iter(PT194), 78, 194))


# -- reading ENSDF fields -------------------------------------------------------------------------------------------


@pytest.mark.parametrize("text, unc, expected", [
    ("328.464", "12", (328.464, 0.012)), ("100", "5", (100.0, 5.0)), ("1.5E3", "2", (1500.0, 200.0)),
    ("0.0635", "", (0.0635, None)), ("2.5", "+5-3", (2.5, 0.5)), ("7", "LT", (7.0, None)),
])
def test_numbers_with_uncertainty_in_the_last_digits(text, unc, expected):
    value, err = ensdf.number(text, unc)
    assert value == pytest.approx(expected[0])
    assert err == (None if expected[1] is None else pytest.approx(expected[1]))


def test_fields_that_are_not_numbers():
    assert ensdf.number("1229.5+X") is None
    assert ensdf.number("") is None


def test_half_lives_and_widths():
    assert ensdf.half_life("41.7 PS", "4") == pytest.approx((41.7e-12, 0.4e-12))
    assert ensdf.half_life("STABLE") is None
    # A width of 1 eV is a half-life of ħ ln2 / Γ = 4.562e-16 s.
    assert ensdf.half_life("1.0 EV")[0] == pytest.approx(4.5624e-16, rel=1e-4)


@pytest.mark.parametrize("text, expected", [
    ("2+", (2.0, 1, False)), ("(2+)", (2.0, 1, True)), ("3/2-", (1.5, -1, False)), ("(3/2)-", (1.5, -1, True)),
    ("4", (4.0, None, False)), ("2+,3+", (None, None, True)), ("", (None, None, False)),
])
def test_spin_and_parity(text, expected):
    assert ensdf.spin_parity(text) == expected


def test_only_the_adopted_dataset_of_the_nuclide_is_read():
    data = ensdf.parse(iter(PT194), 78, 194)
    assert data.nuclide == "194Pt" and data.reference == "NDS 177, 1 (2021)" and data.date == "202109"
    # The level with energy "1229.5+X" cannot be placed, so it and its γ ray are left out.
    assert [lv.energy[0] for lv in data.levels] == [0.0, 328.464, 622.004, 811.25, 1432.5]
    assert data.levels[1].energy[1] == pytest.approx(0.012)
    assert data.levels[0].stable and not data.levels[1].stable
    assert data.levels[1].moments["E2"] == pytest.approx((0.48, 0.14))
    assert ensdf.parse(iter(PT194), 76, 194).levels[1].energy[0] == 218.5
    assert ensdf.parse(iter(PT194), 79, 194) is None


def test_gamma_rays_find_their_final_levels():
    data = ensdf.parse(iter(PT194), 78, 194)
    finals = [[g.final for g in lv.gammas] for lv in data.levels]
    # 1432.5 − 555.5 = 877 keV matches no level, so that γ ray is left out and counted.
    assert finals == [[], [0], [1, 0], [1], [1]]
    assert data.unplaced == 1
    g = data.levels[2].gammas[0]
    assert g.multipolarity == "M1+E2" and g.mixing[0] == 14 and g.conversion[0] == pytest.approx(0.105)
    assert g.weisskopf["E2"] == pytest.approx((89.0, 11.0))
    assert g.weisskopf["M1"][0] == pytest.approx(7.7e-5)


# -- formulas -------------------------------------------------------------------------------------------------------


def test_weisskopf_units_and_rates_match_the_textbook_constants():
    # B_W(E2) = 0.0594 A^{4/3} e²fm⁴ and B_W(E1) = 0.0645 A^{2/3} e²fm².
    assert weisskopf_unit(2, 194) == pytest.approx(0.0594 * 194 ** (4 / 3), rel=1e-3)
    assert weisskopf_unit(1, 194) == pytest.approx(0.0645 * 194 ** (2 / 3), rel=2e-3)
    # T(E1) = 1.587e15 E³ B, T(E2) = 1.223e9 E⁵ B, T(E3) = 5.698e2 E⁷ B (E in MeV, B in e² fm^{2λ}). The textbook
    # constants were computed with older values of e² and ħc, and differ by 0.2% from CODATA 2018.
    for lam, constant in ((1, 1.587e15), (2, 1.223e9), (3, 5.698e2)):
        assert rate_per_b(lam, 1000.0) == pytest.approx(constant, rel=3e-3)
    assert rate_per_b(2, 500.0) == pytest.approx(rate_per_b(2, 1000.0) / 32)


def test_quadrupole_factor():
    assert quadrupole_factor(2) == pytest.approx(0.7579, abs=1e-4)
    assert quadrupole_factor(0) == 0 and quadrupole_factor(0.5) == 0


def test_multipolarity_shares():
    assert components("(E2)") == ["E2"] and components("M1+E2") == ["M1", "E2"]
    assert components("D+Q") is None and components("M1,E2") is None and components("") is None
    assert fractions("E2", None) == {"E2": 1.0}
    assert fractions("M1+E2", 2.0) == pytest.approx({"M1": 0.2, "E2": 0.8})
    assert fractions("M1+E2", None) == {}


def test_matrix_element_units():
    assert parse_matrix_element("1.28 eb", 2) == pytest.approx(128.0)
    assert parse_matrix_element("128 efm2", 2) == 128.0
    assert parse_matrix_element("0.5 eb1.5", 3) == pytest.approx(500.0)
    with pytest.raises(ValueError, match="with its unit"):
        parse_matrix_element("1.28", 2)


# -- deriving matrix elements ---------------------------------------------------------------------------------------


def test_matrix_element_from_the_half_life():
    s = pt194()
    s.derive()
    # By hand: rate = ln2 / 41.7 ps / (1 + 0.0635); B(E2↓) = rate / (1.2252e9 × 0.328464⁵) e²fm⁴.
    b_down = math.log(2) / 41.7e-12 / 1.0635 / (1.2252e9 * 0.328464**5)
    assert s.b(1, 0, "E2") == pytest.approx(b_down, rel=1e-3)
    assert s.b(0, 1, "E2") == pytest.approx(5 * b_down, rel=1e-3)  # B↑ = (2J_f + 1)/(2J_i + 1) B↓
    assert s.b_weisskopf(1, 0) == pytest.approx(50.0, rel=0.02)
    v = s.matrix_element(0, 1, "E2")
    assert v.source == "derived" and "half-life" in v.note and "sign assumed" in v.note
    assert v.unc / v.value == pytest.approx(0.5 * 0.4 / 41.7, rel=1e-6)
    assert s.matrix_element(1, 0, "E2") is v


def test_weisskopf_value_is_preferred_and_mixing_is_respected():
    s = pt194()
    s.derive()
    # 622 → 328 (level 2 → 1) has B(E2) = 89 W.u. in the file; 622 → 0 has 0.29 W.u.
    assert s.b_weisskopf(2, 1) == pytest.approx(89.0)
    assert s.b_weisskopf(2, 0) == pytest.approx(0.29)
    assert "W.u." in s.matrix_element(1, 2).note
    # Without the W.u. value, the half-life route gives the E2 share δ²/(1+δ²) of the 293.5 keV branch.
    s2 = pt194()
    for t in s2.transitions:
        t.weisskopf.clear()
    s2.derive()
    total = 100 * 1.105 + 9.3 * 1.0126
    rate = math.log(2) / 35e-12 * 100 / total * (196 / 197)
    assert s2.b(2, 1, "E2") == pytest.approx(rate / rate_per_b(2, 293.541), rel=1e-9)


def test_quadrupole_moment_becomes_a_diagonal_matrix_element():
    s = pt194()
    assert s.quadrupole_moment(1) == pytest.approx(48.0)  # 0.48 b in e fm²
    v = s.matrix_element(1, 1, "E2")
    assert v.value == pytest.approx(48.0 / 0.7579, rel=1e-4) and v.source == "derived"
    assert s.quadrupole_moment(0) is None
    # A moment written without a sign: only its size is known, and the note says so.
    assert s.quadrupole_moment(3) == pytest.approx(50.0)
    assert "no sign" in s.matrix_element(3, 3).note and "no sign" not in v.note


def test_a_nuclide_without_placed_levels_is_reported(tmp_path, monkeypatch):
    monkeypatch.setenv(ensdf.ENV, str(tmp_path))
    lines = [card("194AU", " ", [(10, "ADOPTED LEVELS")]), card("194AU", "L", [(10, "0+X"), (22, "1-")])]
    (tmp_path / "ensdf.194").write_text("\n".join(lines))
    ensdf._adopted.cache_clear()
    with pytest.raises(ensdf.EnsdfMissing, match="no level of 194Au with a known energy"):
        LevelScheme.from_ensdf("194Au")
    ensdf._adopted.cache_clear()


def test_the_users_value_takes_precedence():
    s = pt194()
    s.set_matrix_element(0, 1, "E2", "1.20 eb")
    s.derive()
    v = s.matrix_element(0, 1, "E2")
    assert (v.value, v.source) == (pytest.approx(120.0), "user")
    assert s.b(0, 1) == pytest.approx(120.0**2)
    with pytest.raises(ValueError, match="no level 9"):
        s.set_matrix_element(0, 9, "E2", "1 eb")
    with pytest.raises(ValueError, match="E1, E2, E3"):
        s.set_matrix_element(0, 1, "M1", 1.0)


def test_levels_without_what_is_needed_get_no_matrix_element():
    s = pt194()
    s.derive()
    # The (3-) level has an E1 γ ray but no half-life and no W.u. value.
    assert s.matrix_element(1, 4, "E1") is None
    assert any("left out" in n for n in s.notes)


# -- truncation -----------------------------------------------------------------------------------------------------


def test_select_keeps_reachable_levels_and_renumbers():
    s = pt194()
    s.derive()
    assert s.connected() == {0, 1, 2, 3}
    cut = s.select(700.0)
    assert [lv.energy.value for lv in cut.levels] == [0.0, 328.464, 622.004]
    assert {(t.initial, t.final) for t in cut.transitions} == {(1, 0), (2, 1), (2, 0)}
    assert cut.b(0, 1) == pytest.approx(s.b(0, 1))
    assert any("3 of 5 levels kept" in n for n in cut.notes)
    everything = s.select(None, connected_only=False)
    assert len(everything.levels) == 5


def test_first_excitation():
    s = pt194()
    s.derive()
    n, energy, b_up = s.first_excitation()
    assert (n, energy) == (1, 328.464) and b_up == pytest.approx(s.b(0, 1))
    assert LevelScheme("4He", s.levels[:1]).first_excitation() is None


# -- the setup file -------------------------------------------------------------------------------------------------


def test_setup_file_round_trip_keeps_values_and_sources():
    exp = Experiment.example("coulex_ni58")
    s = pt194()
    s.derive()
    s.set_matrix_element(0, 1, "E2", "1.28 eb", unc=2.0)
    exp.target.material = "194Pt"
    exp.levels["target"] = s
    again = Experiment.from_toml(exp.to_toml())
    t = again.levels["target"]
    assert t.nuclide == "194Pt" and "NDS 177" in t.reference
    assert [lv.label for lv in t.levels] == ["0+", "2+", "2+", "4+", "3-"]
    assert t.levels[4].jpi == "(3-)"
    assert t.levels[1].energy == s.levels[1].energy
    assert t.levels[1].half_life.value == pytest.approx(41.7e-12) and t.levels[1].half_life.source == "ensdf"
    v = t.matrix_element(0, 1)
    assert (v.value, v.unc, v.source) == (pytest.approx(128.0), 2.0, "user")
    assert t.matrix_element(1, 2).source == "derived" and "W.u." in t.matrix_element(1, 2).note
    assert t.transitions[1].mixing.value == 14 and t.transitions[1].weisskopf["E2"].value == 89
    assert again.to_toml() == exp.to_toml()


def test_setup_files_without_levels_are_unchanged():
    exp = Experiment.example("coulex_ni58")
    assert exp.levels == {} and "levels" not in exp.to_dict()


def test_problems_in_a_level_scheme_are_named():
    d = Experiment.example("coulex_ni58").to_dict()
    d["levels"] = {
        "beam": {"nuclide": "18O", "level": [{"energy": "0 keV", "jpi": "0+"}]},
        "target": {"nuclide": "58Ni",
                   "level": [{"energy": "5 keV", "jpi": "0+"}, {"energy": "1454", "jpi": "2+"},
                             {"energy": "3 keV", "source": "guess"}],
                   "transition": [{"from": 0, "to": 1, "energy": "1454 keV"}, {"from": 7, "to": 0, "energy": "1 keV"}],
                   "matrix_element": [{"from": 0, "to": 1, "value": "0.26"},
                                      {"from": 0, "to": 1, "multipolarity": "M1", "value": "1 eb"}]},
        "detector": {},
    }
    with pytest.raises(SetupError) as e:
        Experiment.from_dict(d)
    text = "\n".join(e.value.problems)
    for expected in ("levels.beam: nuclide 18O is not the beam (16O)",
                     "levels.target level 1: energy", "has no unit",
                     "source must be one of ensdf, derived, assumed, user, got 'guess'",
                     "the first level must be the ground state",
                     "levels.target transition 0: a transition goes from a higher level to a lower one",
                     "levels.target transition 1: 'from' and 'to' must be level numbers",
                     "levels.target matrix element 0: value: write an E2 matrix element with its unit",
                     "levels.target matrix element 1: multipolarity must be one of E1, E2, E3",
                     "levels: unknown section 'detector'"):
        assert expected in text, expected


def test_target_scheme_must_be_a_nuclide_of_the_target():
    d = Experiment.example("coulex_ni58").to_dict()
    d["levels"] = {"target": {"nuclide": "60Ni", "level": [{"energy": "0 keV", "jpi": "0+"}]}}
    with pytest.raises(SetupError, match="60Ni is not in the target"):
        Experiment.from_dict(d)


# -- the local copy -------------------------------------------------------------------------------------------------


def test_a_missing_copy_says_what_to_do(tmp_path, monkeypatch):
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "nothing"))
    assert not ensdf.available()
    with pytest.raises(ensdf.EnsdfMissing, match="fetch_ensdf.py"):
        LevelScheme.from_ensdf("194Pt")


def test_local_copy_is_read_from_files_and_zips(tmp_path, monkeypatch):
    import zipfile

    monkeypatch.setenv(ensdf.ENV, str(tmp_path))
    (tmp_path / "README.md").write_text("about this copy")
    assert not ensdf.available()
    with zipfile.ZipFile(tmp_path / "ensdf_test.zip", "w") as zf:
        zf.writestr("ensdf.193", "\n".join(PT194).replace("194PT", "193PT"))
        zf.writestr("ensdf.194", "\n".join(PT194))
    ensdf._adopted.cache_clear()
    assert ensdf.available()
    s = LevelScheme.from_ensdf("194Pt", max_energy_kev=1000)
    assert [lv.label for lv in s.levels] == ["0+", "2+", "2+", "4+"]
    assert s.b_weisskopf(1, 0) == pytest.approx(50.0, rel=0.02)
    with pytest.raises(ensdf.EnsdfMissing, match="no adopted levels for 196Pt"):
        LevelScheme.from_ensdf("196Pt")
    ensdf._adopted.cache_clear()


@needs_ensdf
@pytest.mark.parametrize("nuclide, energy, jpi, wu", [
    ("194Pt", 328.5, "2+", (40, 60)), ("58Ni", 1454.2, "2+", (5, 15)), ("208Pb", 2614.5, "3-", (25, 45)),
    ("20Ne", 1633.7, "2+", (15, 25)), ("152Sm", 121.8, "2+", (120, 170)),
])
def test_first_excited_states_from_the_real_copy(nuclide, energy, jpi, wu):
    s = LevelScheme.from_ensdf(nuclide, max_energy_kev=3000)
    mult = "E3" if jpi == "3-" else "E2"
    n, e, _ = s.first_excitation(mult)
    assert e == pytest.approx(energy, abs=0.5) and s.levels[n].label == jpi
    assert wu[0] < s.b_weisskopf(n, 0, mult) < wu[1]


@needs_ensdf
def test_an_odd_mass_nucleus_from_the_real_copy():
    s = LevelScheme.from_ensdf("195Pt", max_energy_kev=500)
    assert s.levels[0].label == "1/2-" and len(s.levels) > 3
    assert all(m.value.value > 0 for m in s.matrix_elements if m.a != m.b)


# -- the planner ----------------------------------------------------------------------------------------------------


@pytest.fixture
def copy_with_ni58(tmp_path, monkeypatch):
    """A local 'ENSDF' folder holding a hand-typed scheme for 58Ni (values close to the evaluated ones)."""
    def c(kind, cols, cont=" "):
        return card(" 58NI", kind, cols, cont)

    lines = [c(" ", [(10, "ADOPTED LEVELS, GAMMAS")]),
             c("L", [(10, "0.0"), (22, "0+"), (40, "STABLE")]),
             c("L", [(10, "1454.21"), (20, "9"), (22, "2+"), (40, "0.652 PS"), (50, "23")]),
             c("G", [(10, "1454.20"), (20, "9"), (22, "100"), (32, "E2")]),
             c("G", [(10, "BE2W=10.0 4")], cont="2"),
             c("L", [(10, "2459.21"), (20, "14"), (22, "4+"), (40, "1.5 PS")]),
             c("G", [(10, "1005.0"), (22, "100"), (32, "E2")])]
    monkeypatch.setenv(ensdf.ENV, str(tmp_path))
    (tmp_path / "ensdf.058").write_text("\n".join(lines))
    ensdf._adopted.cache_clear()
    yield
    ensdf._adopted.cache_clear()


def test_planner_looks_up_a_scheme_and_plans_a_state(copy_with_ni58):
    from physim.nuclear.planner import Planner

    p = Planner.example("coulex_ni58")
    info = p.levels()
    assert info["ensdf"] and info["target"] is None
    assert info["nuclides"] == {"beam": "16O", "target": "58Ni", "target_choices": ["58Ni"]}
    assert p.lookup_levels("target", "3 MeV")
    t = p.levels()["target"]
    assert [lv["jpi"] for lv in t["levels"]] == ["0+", "2+", "4+"]
    row = t["matrix_elements"][0]
    assert row["source"] == "derived" and row["b_up_wu"] == pytest.approx(50.0)  # 5 × 10 W.u. for 0+ → 2+
    assert t["transitions"][0]["branching"] == 1.0
    # Planning the 2+ state fills in the reaction from the scheme.
    assert p.use_state("target", 1)
    exc = p.experiment.excitation
    assert exc.excite == "target" and exc.energy_mev == pytest.approx(1.45421)
    assert exc.b_up_e2fm == pytest.approx(50.0 * weisskopf_unit(2, 58), rel=1e-5)
    assert p.gamma()["available"]
    # The user's matrix element survives a second look-up, and is marked as theirs.
    assert p.set_matrix_element("target", 0, 1, "E2", "0.25 eb")
    assert p.lookup_levels("target", "2 MeV")
    t = p.levels()["target"]
    assert len(t["levels"]) == 2
    assert (t["matrix_elements"][0]["value_efm"], t["matrix_elements"][0]["source"]) == (pytest.approx(25.0), "user")
    assert "[levels.target]" in p.to_toml() and Planner(Experiment.from_toml(p.to_toml())).levels()["target"]
    assert p.remove_levels("target") and "levels" not in p.to_toml()


def test_planner_reports_a_missing_copy_or_nuclide(copy_with_ni58):
    from physim.nuclear.planner import Planner

    p = Planner.example("coulex_ni58")
    with pytest.raises(ensdf.EnsdfMissing, match="no adopted levels for 16O"):
        p.lookup_levels("beam")
    with pytest.raises(ValueError, match="no E2 matrix element"):
        p.lookup_levels("target")
        p.use_state("target", 2)


def test_level_diagram_renders():
    pytest.importorskip("plotly")
    from physim.nuclear import app

    s = pt194()
    s.derive()
    fig = app.figure_levels(s.select(1000.0).table())
    assert len(fig.data) == 4  # one line per level
    assert len([a for a in fig.layout.annotations if a.showarrow]) == 4  # one arrow per transition
    for theme in app.THEMES:
        app.themed(fig, theme).to_json()
