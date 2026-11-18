"""Shell-model results run elsewhere: KSHELL folders and logs, KSHELL's summary, and physim's table format.

The KSHELL logs and summary are the wave 0 fixtures in ``tests/reference/nuclear/shellmodel/kshell``; the tables
are in ``tests/reference/nuclear/shellmodel/import`` (its README says what each is). Cross-checks: the summary
against the full logs it summarises, the hand-typed table against the logs it was typed from, the template
against the table (the same numbers in other units), and round trips result → table → result.
"""

import math
import shutil
from pathlib import Path

import pytest

from physim.nuclear import shellmodel as sm
from physim.nuclear import shellmodel_import as si
from physim.nuclear.levels import Level, LevelScheme, MatrixElement, Value, quadrupole_factor, weisskopf_unit

REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "shellmodel"
SI28 = REFERENCE / "kshell" / "si28_usdb" / "logs"
SI28_EFF = REFERENCE / "kshell" / "si28_usdb" / "logs_eff_1.31_0.46"
SUMMARY = SI28 / "summary_Si28_usdb.txt"
CR48 = REFERENCE / "kshell" / "cr48_kb3" / "log_cr48_m0.txt"
TABLE = REFERENCE / "import" / "si28_usdb_table.csv"
TEMPLATE = REFERENCE / "import" / "table_template.csv"


def raw(path: Path) -> str:
    return path.read_bytes().decode("utf-8")


@pytest.fixture(scope="module")
def full():
    return si.import_path(SI28)


@pytest.fixture(scope="module")
def summary():
    return si.import_path(SUMMARY)


@pytest.fixture(scope="module")
def eff():
    """The 1.31/0.46 transit runs with the states of the kshell logs: what the table was typed from."""
    return sm.read_logs([raw(p) for p in sorted(SI28.glob("log_Si28_usdb_j*p.txt"))],
                        [raw(p) for p in sorted(SI28_EFF.glob("log_tr_*.txt"))])


def off_diagonal(result):
    return {(t.kind, t.initial, t.final): t for t in result.transitions if not t.diagonal}


# ---------------------------------------------------------------------------------------------------------------
# KSHELL logs


def test_a_kshell_folder_reads_as_read_folder_does_with_its_source_noted(full):
    plain = sm.read_folder(SI28)
    assert full.to_dict()["states"] == plain.to_dict()["states"]
    assert full.to_dict()["transitions"] == plain.to_dict()["transitions"]
    assert (full.nuclide, full.interaction, full.eff_charge) == ("28Si", "usdb", (1.5, 0.5))
    p = si.provenance(full)
    assert p["format"] == "KSHELL logs" and Path(p["path"]) == SI28.resolve() and p["reference"] is None
    assert "read 3 kshell.exe and 5 transit.exe logs" in full.notes
    # the summary beside the logs is not used, but checked against them
    assert any("summary_Si28_usdb.txt" in n and "agrees" in n for n in full.notes)
    assert full.notes[1:] == [n for n in full.notes if not n.startswith("imported from")]


def test_a_list_of_logs_and_a_single_log(full):
    logs = sorted(SI28.glob("log_*.txt"))
    listed = si.import_path(logs)
    assert listed.to_dict()["states"] == full.to_dict()["states"]
    assert len(listed.transitions) == len(full.transitions)
    one = si.import_path(SI28 / "log_Si28_usdb_j4p.txt")
    assert [s.label for s in one.states] == ["2+1", "2+2", "2+3"]
    assert si.provenance(one)["path"] == str((SI28 / "log_Si28_usdb_j4p.txt").resolve())


def test_transit_logs_alone_still_give_the_states():
    r = si.import_path(SI28_EFF)
    assert r.nuclide == "28Si" and r.eff_charge == (1.31, 0.46)
    assert r.transition("E2", (4, 1, 1), (0, 1, 1)).b_down == pytest.approx(78.2005)
    assert "read 0 kshell.exe and 3 transit.exe logs" in r.notes


def test_cr48_log_that_is_not_j_projected():
    r = si.import_path(CR48)
    assert (r.nuclide, r.interaction) == ("48Cr", "kb3")  # from the core and valence nucleons, and FN_INT
    assert len(r.states) == 32 and r.dimensions == {(0, 1): 1963461}
    assert r.states[0].label == "0+1" and r.states[0].energy_mev == pytest.approx(-32.95367)
    assert {s.j2 for s in r.states} >= {0, 4, 8, 12}  # an M = 0 run mixes J values
    assert si.provenance(r)["format"] == "KSHELL logs"


def test_nuclide_and_interaction_given_are_checked_or_take_over():
    with pytest.raises(ValueError, match="are for 48Cr, not 50Cr"):
        si.import_path(CR48, nuclide="50Cr")
    r = si.import_path(CR48, nuclide="Cr-48", interaction="KB3")
    assert r.nuclide == "48Cr" and r.interaction == "KB3"
    assert not any("called" in n for n in r.notes)
    r = si.import_path(CR48, interaction="kb3g")
    assert r.interaction == "kb3g"
    assert "the files name the interaction 'kb3'; called 'kb3g' as asked" in r.notes


def test_a_summary_that_disagrees_with_the_logs_beside_it_is_named(tmp_path):
    for p in SI28.glob("log_*.txt"):
        shutil.copy(p, tmp_path)
    (tmp_path / "VERSION").write_text("KSHELL test build\n", encoding="utf-8")
    text = raw(SUMMARY).replace("-133.929    1.932", "-133.900    1.961")
    (tmp_path / "summary_Si28_usdb.txt").write_bytes(text.encode("utf-8"))
    (tmp_path / "Si28_usdb_j0p.wav").write_bytes(b"\0\1\2")  # other files of a run folder are skipped
    r = si.import_path(tmp_path)
    assert r.kshell_version == "KSHELL test build"
    assert any("disagrees with the logs for 2+1" in n for n in r.notes)


# ---------------------------------------------------------------------------------------------------------------
# KSHELL's summary


def test_summary_alone_matches_the_full_logs_to_its_rounding(summary, full):
    assert (summary.nuclide, summary.interaction, summary.space) == ("28Si", "usdb", full.space)
    assert summary.eff_charge is None and summary.gl is None and summary.request is None
    assert [s.key for s in summary.states] == [s.key for s in sorted(full.states, key=lambda s: s.energy_mev)]
    for s in summary.states:
        f = full.state(*s.key)
        assert s.energy_mev == pytest.approx(f.energy_mev, abs=0.0005 + 1e-9)
        assert s.ex_mev == pytest.approx(f.ex_mev, abs=0.001 + 1e-9)
        assert s.t2 == f.t2 == 0
        assert s.q_s is None and s.mu is None
    got, want = off_diagonal(summary), off_diagonal(full)
    assert set(got) == set(want)
    for k, t in got.items():
        w = want[k]
        h = 0.05 if t.kind == "E2" else 0.0005  # the summary prints B(E2) to 0.1, B(M1) to 0.001
        assert t.b_down == pytest.approx(w.b_down, abs=h + 1e-9)
        assert t.b_up == pytest.approx(w.b_up, abs=h + 1e-9)
        assert t.mred >= 0 and t.mred ** 2 == pytest.approx(w.mred ** 2, abs=h * (min(k[1][0], k[2][0]) + 1) + 1e-6)
    assert si.unsigned(summary) == set(got)
    assert si.provenance(summary) == {"path": str(SUMMARY.resolve()), "format": "KSHELL summary", "reference": None}


def test_summary_numbers_are_kshell_s_printed_ones(summary):
    t = summary.transition("E2", (4, 1, 1), (0, 1, 1))
    assert (t.initial, t.final, t.b_down, t.b_up) == ((4, 1, 1), (0, 1, 1), 99.8, 499.2)
    # the 0+2 at 4.837 MeV lies above the 2+1: B(E2; 0+2 → 2+1) is the "->" column
    t = summary.transition("E2", (0, 1, 2), (4, 1, 1))
    assert (t.initial, t.b_down, t.b_up) == ((0, 1, 2), 82.4, 16.5)
    assert t.mred == pytest.approx(math.sqrt(82.4))  # from B↓ × 1, more precise than B↑ × 5
    assert summary.state(0, 1, 2).ex_mev == 4.836  # the level table's value, not the B table's 4.837


def test_summary_in_a_folder_of_its_own_and_without_names(tmp_path):
    shutil.copy(SUMMARY, tmp_path)
    r = si.import_path(tmp_path)
    assert si.provenance(r)["format"] == "KSHELL summary" and r.nuclide == "28Si"
    bare = "\n".join(line.split("log_")[0] for line in raw(SUMMARY).splitlines())
    (tmp_path / "out.txt").write_text(bare, encoding="utf-8")
    with pytest.raises(ValueError, match="pass nuclide="):
        si.import_path(tmp_path / "out.txt")
    with pytest.raises(ValueError, match="pass interaction="):
        si.import_path(tmp_path / "out.txt", nuclide="28Si")
    r = si.import_path(tmp_path / "out.txt", nuclide="28Si", interaction="USDB")
    assert (r.nuclide, r.interaction, len(r.states)) == ("28Si", "usdb", 8)
    with pytest.raises(ValueError, match="mass = 28"):
        si.parse_summary(bare, nuclide="30Si", interaction="usdb")  # the B tables give the mass
    with pytest.raises(ValueError, match="several KSHELL summaries"):
        shutil.copy(SUMMARY, tmp_path / "summary_Si28_usdb_again.txt")
        si.import_path(tmp_path)


def test_parse_summary_mass_check_against_names():
    text = raw(SUMMARY).replace("mass = 28", "mass = 30")
    with pytest.raises(ValueError, match="name 28Si but the B tables say mass = 30"):
        si.parse_summary(text)


# ---------------------------------------------------------------------------------------------------------------
# physim's table format


def test_the_typed_table_matches_the_logs_it_was_typed_from(eff):
    r = si.import_path(TABLE)
    assert (r.nuclide, r.interaction, r.eff_charge, r.gl, r.gs) == ("28Si", "usdb", (1.31, 0.46), (1.0, 0.0),
                                                                     (5.585, -3.826))
    assert r.kshell_version.startswith("KSHELL 2019")
    p = si.provenance(r)
    assert p["format"] == "physim shell-model table" and p["path"] == str(TABLE.resolve())
    assert p["reference"].startswith("physim's KSHELL fixtures") and "logs_eff_1.31_0.46/" in p["reference"]
    assert any("Frycz et al." in n and "78.3 e2fm4" in n for n in r.notes)
    assert [s.label for s in r.states] == ["0+1", "2+1", "4+1"]
    for s in r.states:
        e = eff.state(*s.key)
        assert s.energy_mev == pytest.approx(e.energy_mev, abs=0.0006)
        assert s.ex_mev == pytest.approx(e.ex_mev, abs=0.0006)
    assert r.state(4, 1, 1).q_s == 18.50 and eff.state(4, 1, 1).q_s == pytest.approx(18.50, abs=0.005)
    assert r.state(4, 1, 1).mu == eff.state(4, 1, 1).mu == 1.0312
    assert r.state(0, 1, 1).t2 == 0 and r.state(4, 1, 1).t2 is None
    for k in (((4, 1, 1), (0, 1, 1)), ((8, 1, 1), (4, 1, 1))):
        t, e = r.transition("E2", *k), eff.transition("E2", *k)
        assert (t.initial, t.final) == (e.initial, e.final)
        assert t.mred == e.mred
        assert t.b_down == pytest.approx(e.b_down, abs=0.005)
        assert t.b_up == pytest.approx(e.b_up, abs=0.05)
    assert r.transition("E2", (4, 1, 1), (0, 1, 1)).b_down == 78.20  # kept as typed
    assert si.unsigned(r) == set()


def test_the_template_gives_the_same_numbers_in_other_units():
    t, x = si.import_path(TEMPLATE), si.import_path(TABLE)
    assert (t.nuclide, t.interaction, t.eff_charge) == (x.nuclide, x.interaction, x.eff_charge)
    assert si.provenance(t)["reference"].startswith("TEMPLATE")
    for a, b in zip(t.states, x.states):
        assert a.key == b.key
        assert a.energy_mev == pytest.approx(b.energy_mev, abs=1e-9)  # keV → MeV
        assert a.ex_mev == pytest.approx(b.ex_mev, abs=1e-9)
        assert a.q_s == pytest.approx(b.q_s) if b.q_s is not None else a.q_s is None  # e b → e fm²
    up = t.transition("E2", (4, 1, 1), (0, 1, 1))  # typed as B(E2; 0+1 → 2+1) in e²b²
    assert up.b_up == pytest.approx(391.0) and up.b_down == pytest.approx(78.2005, abs=1e-3)
    assert up.mred == pytest.approx(-19.7738)  # e b → e fm²
    wu = t.transition("E2", (8, 1, 1), (4, 1, 1))  # B(E2) in W.u., M(E2) written from the lower state up
    assert wu.b_down == pytest.approx(21.78 * weisskopf_unit(2, 28))
    assert wu.mred == pytest.approx(-31.4671)  # (−1)^(J_a − J_b) = +1 for ΔJ = 2
    assert wu.b_up == pytest.approx(31.4671 ** 2 / 5)  # B↓ is kept as given, B↑ follows from M


def test_round_trip_kshell_result_to_table_and_back(full, tmp_path):
    text = si.to_table(full, tmp_path / "si28.csv")
    assert (tmp_path / "si28.csv").read_text(encoding="utf-8") == text
    back = si.import_path(tmp_path / "si28.csv")
    assert (back.nuclide, back.interaction, back.space, back.eff_charge, back.gl, back.gs, back.dimensions) == \
        (full.nuclide, full.interaction, full.space, full.eff_charge, full.gl, full.gs, full.dimensions)
    assert back.kshell_version == full.kshell_version == "unknown"
    assert back.notes[1:] == full.notes  # the old provenance stays, after the new one
    assert si.provenance(back)["path"] == str((tmp_path / "si28.csv").resolve())
    assert len(back.states) == len(full.states)
    for s in back.states:
        f = full.state(*s.key)
        for name in ("energy_mev", "ex_mev", "q_s", "mu"):
            got, want = getattr(s, name), getattr(f, name)
            assert (got is None and want is None) or got == pytest.approx(want, rel=1e-9, abs=1e-9), name
        assert s.t2 == f.t2
    got, want = off_diagonal(back), off_diagonal(full)
    assert set(got) == set(want)
    for k, t in got.items():
        w = want[k]
        assert t.mred == w.mred and t.b_down == w.b_down
        assert t.b_up == pytest.approx(w.b_up, rel=2e-5, abs=2e-4)  # from M, printed to 4 decimals
    # and once more: the table of the table is the same table, apart from the new provenance note
    again = si.to_table(back)
    assert again.replace("# notes: " + back.notes[0] + "\n", "") == text


def test_round_trip_keeps_unknown_signs_unknown(summary):
    back = si.parse_table(si.to_table(summary))
    assert si.unsigned(back) == si.unsigned(summary) != set()
    assert "M(E2)" not in si.to_table(summary)
    assert si.provenance(back)["path"] == "text given directly"
    for k, t in off_diagonal(back).items():
        w = off_diagonal(summary)[k]
        assert (t.b_down, t.b_up, t.mred) == (w.b_down, w.b_up, w.mred)


def test_round_trip_of_the_typed_table(tmp_path):
    first = si.import_path(TABLE)
    back = si.import_path(Path(si.to_table(first, tmp_path / "t.csv") and tmp_path / "t.csv"))
    assert si.provenance(back)["reference"] == si.provenance(first)["reference"]
    assert [(s.key, s.energy_mev, s.ex_mev, s.q_s, s.mu, s.t2) for s in back.states] == \
        [(s.key, s.energy_mev, s.ex_mev, s.q_s, s.mu, s.t2) for s in first.states]
    assert [(t.kind, t.initial, t.final, t.mred, t.b_down) for t in back.transitions] == \
        [(t.kind, t.initial, t.final, t.mred, t.b_down) for t in first.transitions]


def test_imported_table_goes_into_a_level_scheme_through_apply():
    r = si.import_path(TABLE)
    s = LevelScheme("28Si", [Level(Value(0.0, 0.1, "ensdf"), 0, 1, "0+"), Level(Value(1779.0, 0.1, "ensdf"), 2, 1, "2+")])
    s.matrix_elements.append(MatrixElement(0, 1, "E2", Value(18.0, 0.5, "derived", "from B(E2) in ENSDF")))
    changes = sm.apply(r, s, ["E2", "Q", "new"])
    assert s.levels[1].energy.source == "ensdf"  # never replaced
    assert s.matrix_element(0, 1, "E2").source == "derived"
    q = s.matrix_element(1, 1, "E2")
    assert q.source == "shell model" and q.value * quadrupole_factor(2) == pytest.approx(18.50)
    four = next(i for i, lv in enumerate(s.levels) if lv.spin == 4)
    assert s.levels[four].energy.value == pytest.approx(4607.0)
    assert s.matrix_element(1, four, "E2").value == pytest.approx(-31.4671)
    assert any("kept" in c for c in changes)


def test_odd_a_labels_turned_round_m1_and_unsigned_b():
    text = """# format: physim shell-model table 1
# nuclide: 25Mg
# interaction: a model of the test
## synthetic numbers, only to test the reading
kind,state,final,value,unit,comment
Ex,5/2+,,0,MeV,
Ex,1/2+1,,0.6,MeV,
Ex,3/2+1,,0.9,MeV,
Ex,7/2+1,,1.6,MeV,
M(M1),5/2+,7/2+1,0.5,mu_N,written upwards
B(M1),7/2+1,5/2+,0.01745,W.u.,downwards: 0.5²/8 μ_N² in W.u.
B(E2),3/2+1,1/2+1,40,e2fm4,no M row: sign unknown
M(E2),5/2+1,5/2+1,10,e fm2,diagonal: Q_s
"""
    r = si.parse_table(text)
    assert r.interaction == "a model of the test" and r.nuclide == "25Mg"
    assert [s.label for s in r.states] == ["5/2+1", "1/2+1", "3/2+1", "7/2+1"]
    assert "the table gives excitation energies only; energy_mev holds the excitation energy" in r.notes
    m1 = r.transition("M1", (7, 1, 1), (5, 1, 1))
    assert (m1.initial, m1.final) == ((7, 1, 1), (5, 1, 1))
    assert m1.mred == pytest.approx(-0.5)  # ⟨5/2‖M1‖7/2⟩ = (−1)^(5/2 − 7/2) ⟨7/2‖M1‖5/2⟩
    assert m1.b_down == pytest.approx(0.01745 * si.WEISSKOPF_M1) and m1.b_up == pytest.approx(0.25 / 6)
    assert si.WEISSKOPF_M1 == pytest.approx(1.790, abs=5e-4)
    assert r.state(5, 1, 1).q_s == pytest.approx(10 * quadrupole_factor(2.5))
    assert si.unsigned(r) == {("E2", (3, 1, 1), (1, 1, 1))}
    assert r.transition("E2", (3, 1, 1), (1, 1, 1)).mred == pytest.approx(math.sqrt(160))


# ---------------------------------------------------------------------------------------------------------------
# Errors


GOOD_HEAD = "# format: physim shell-model table 1\n# nuclide: 28Si\n# interaction: usdb\n"
COLUMNS = "kind,state,final,value,unit\n"
LEVELS = "Ex,0+1,,0,MeV\nEx,2+1,,1.8,MeV\n"


@pytest.mark.parametrize("text, message", [
    ("# nuclide: 28Si\n# interaction: usdb\n" + COLUMNS + LEVELS, "format: physim shell-model table 1"),
    (GOOD_HEAD.replace("table 1", "table 2") + COLUMNS + LEVELS, "version 2"),
    (GOOD_HEAD + "# efective charges: 1.5, 0.5\n" + COLUMNS + LEVELS, "'efective charges' is not a header entry"),
    (GOOD_HEAD + "# effective charges: 1.5\n" + COLUMNS + LEVELS, "two numbers"),
    (GOOD_HEAD + "# nuclide: 30Si\n" + COLUMNS + LEVELS, "'nuclide' is given twice"),
    (GOOD_HEAD.replace("# nuclide: 28Si\n", "") + COLUMNS + LEVELS, "no nuclide"),
    (GOOD_HEAD.replace("# interaction: usdb\n", "") + COLUMNS + LEVELS, "no interaction"),
    (GOOD_HEAD, "no values"),
    (GOOD_HEAD + "kind,state,value,unit\n" + LEVELS, "missing final"),
    (GOOD_HEAD + "kind,state,final,value,unit,sigma\n" + LEVELS, "'sigma' not known"),
    (GOOD_HEAD + COLUMNS + "Energy,0+1,,0,MeV\n", "'Energy' is not a kind"),
    (GOOD_HEAD + COLUMNS + "Ex,0+1,,0,eV\n", "'eV' is not a unit for Ex"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,0+1,78,e fm2\n", "not a unit for B\\(E2\\)"),
    (GOOD_HEAD + COLUMNS + "Ex,two plus,,0,MeV\n", "cannot read 'twoplus' as a state"),
    (GOOD_HEAD + COLUMNS + "Ex,5/2+1,,0,MeV\n", "5/2\\+ cannot be a state of 28Si"),
    (GOOD_HEAD + COLUMNS + "Ex,2+0,,0,MeV\n", "indices start at 1"),
    (GOOD_HEAD + COLUMNS + "Ex,0+1,,zero,MeV\n", "'zero' is not a number"),
    (GOOD_HEAD + COLUMNS + "Ex,0+1,,,MeV\n", "no value"),
    (GOOD_HEAD + COLUMNS + ",,,0,MeV\n", "'' is not a kind"),
    (GOOD_HEAD + COLUMNS + "Ex,,,0,MeV\n", "no state"),
    (GOOD_HEAD + COLUMNS + LEVELS + "Ex,2+1,,1.9,MeV\n", "Ex of 2\\+1 is given twice \\(line 6\\)"),
    (GOOD_HEAD + COLUMNS + LEVELS + "Q,2+1,0+1,10,e fm2\n", "leave 'final' empty"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,,78,e2fm4\n", "needs a final state"),
    (GOOD_HEAD + COLUMNS + LEVELS + "T,0+1,,0,MeV\n", "isospin has no unit"),
    (GOOD_HEAD + COLUMNS + LEVELS + "T,0+1,,0.3,\n", "isospin is 0, 1/2"),
    (GOOD_HEAD + COLUMNS + "Ex,0+1,,0,MeV\nEx,2+,,1.8,MeV\nEx,2+2,,4,MeV\n", "without an index means the only 2\\+ state"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),4+1,2+1,110,e2fm4\n", "no energy for 4\\+1"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,0+1,78,e2fm4\nM(E2),2+1,0+1,-30,e fm2\n", "disagree"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,0+1,78,e2fm4\nB(E2),2+1,0+1,79,e2fm4\n", "given twice \\(line 7\\)"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,0+1,78,e2fm4\nB(E2),0+1,2+1,300,e2fm4\n", "disagree"),
    (GOOD_HEAD + COLUMNS + LEVELS + "M(E2),2+1,0+1,-19.75,e fm2\nM(E2),0+1,2+1,19.75,e fm2\n", "given twice"),
    (GOOD_HEAD + COLUMNS + "Ex,0+1,,0,MeV\nEx,2+,,1.8,MeV\nEx,2+1,,1.9,MeV\n", "Ex of 2\\+1 is given twice"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(M1),2+1,0+1,1,mu_N2\n", "breaks the selection rules"),
    (GOOD_HEAD + COLUMNS + LEVELS + "Ex,3-1,,5,MeV\nB(E2),3-1,2+1,1,e2fm4\n", "breaks the selection rules"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,0+1,-78,e2fm4\n", "cannot be negative"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,2+1,7,e2fm4\n", "diagonal M\\(E2\\)"),
    (GOOD_HEAD + COLUMNS + LEVELS + "M(E2),0+1,0+1,7,e fm2\n", "no quadrupole moment"),
    (GOOD_HEAD + COLUMNS + LEVELS + "Q,2+1,,10,e fm2\nM(E2),2+1,2+1,20,e fm2\n", "gives Q_s"),
    (GOOD_HEAD + COLUMNS + "E,0+1,,-100.000,MeV\nEx,0+1,,0.000,MeV\nE,2+1,,-98.000,MeV\nEx,2+1,,1.500,MeV\n", "E - Ex"),
    (GOOD_HEAD + COLUMNS + "E,0+1,,-100,MeV\nEx,2+1,,1.5,MeV\n", "same kind"),
    (GOOD_HEAD + COLUMNS + LEVELS + "B(E2),2+1,0+1,78,e2fm4,extra\n", "more cells than columns"),
])
def test_malformed_tables_say_what_is_wrong(text, message):
    with pytest.raises(ValueError, match=message):
        si.parse_table(text)


def test_nuclide_argument_against_the_table_header():
    text = GOOD_HEAD + COLUMNS + LEVELS
    with pytest.raises(ValueError, match="is for 28Si, not 30Si"):
        si.parse_table(text, nuclide="30Si")
    r = si.parse_table(text.replace("# nuclide: 28Si\n", "").replace("# interaction: usdb\n", ""),
                       nuclide="28Si", interaction="USDA")
    assert (r.nuclide, r.interaction) == ("28Si", "usda")
    r = si.parse_table(text, interaction="usda")
    assert r.interaction == "usda" and any("called 'usda' as asked" in n for n in r.notes)


def test_files_that_are_not_shell_model_results(tmp_path):
    with pytest.raises(FileNotFoundError):
        si.import_path(tmp_path / "missing")
    (tmp_path / "notes.txt").write_text("just some text\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cannot tell what"):
        si.import_path(tmp_path / "notes.txt")
    with pytest.raises(ValueError, match="holds no KSHELL logs"):
        si.import_path(tmp_path)
    (tmp_path / "garbage.txt").write_text(" Energy levels\n\n    N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file\n",
                                          encoding="utf-8")
    with pytest.raises(ValueError, match="no energy-level table"):
        si.import_path(tmp_path / "garbage.txt")
    with pytest.raises(ValueError, match="no files given"):
        si.import_path([])
    mixed = tmp_path / "mixed"
    mixed.mkdir()
    shutil.copy(CR48, mixed)
    shutil.copy(TABLE, mixed)
    with pytest.raises(ValueError, match="both KSHELL logs and a physim table"):
        si.import_path(mixed)
    with pytest.raises(ValueError, match="are for 48Cr, not 47Cr"):
        si.import_path(CR48, nuclide="47Cr")


def test_nushellx_output_is_deferred_with_what_is_needed(tmp_path):
    (tmp_path / "si28usdb.lpt").write_text("not a real NuShellX file\n", encoding="utf-8")
    with pytest.raises(NotImplementedError, match="real NuShellX output"):
        si.import_path(tmp_path / "si28usdb.lpt")
    with pytest.raises(NotImplementedError, match="physim shell-model table"):
        si.import_path(tmp_path)
