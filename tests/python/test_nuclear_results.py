"""The results ledger, the registry of observables and the comparison engine (wave 0 of the results pages): rows
that round-trip and are only ever appended, one name per level and transition, and the comparison of the ⁶²Ni
B(E2; 0⁺ → 2⁺₁) values by a rule fixed beforehand. Data model and plumbing: no run is simulated here."""

import json

import pytest

from physim.nuclear import Experiment, compare, ensdf, lifetimes, observables, results
from physim.nuclear.gamma import level_names
from physim.nuclear.levels import Level, LevelScheme, Value
from physim.nuclear.runs import ExperimentFolder
from physim.nuclear.workbench import KINDS

needs_ensdf = pytest.mark.skipif(not ensdf.available(), reason="no local copy of ENSDF (scripts/fetch_ensdf.py)")

KEY = "62Ni:2+1>0+1:E2"
NDS = results.reference("Nichols, Singh and Tuli", 2012, where="Nucl. Data Sheets 113, 973, adopted levels")
#: The first adopted levels of ⁶²Ni as (energy keV, 2J, parity, ordinal), for anchoring without a local copy.
NI62 = ((0.0, 0, 1, 1), (1172.91, 4, 1, 1), (2048.55, 0, 1, 2), (2301.85, 4, 1, 2), (2336.46, 8, 1, 1))


def b_e2(value, method, evidence, origin=NDS, **more):
    return results.Row("62Ni", "B(E2)", KEY, value, "e2fm4", method, evidence, direction="up", origin=origin, **more)


@pytest.fixture
def no_ensdf(tmp_path, monkeypatch):
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))


@pytest.fixture
def ni62_levels(monkeypatch):
    monkeypatch.setattr(results, "ensdf_levels", lambda nuclide: NI62 if nuclide == "62Ni" else ())


def scheme():
    """⁶²Ni cut so that its first 2⁺ and its second 0⁺ are missing: the scheme's own numbering is not ENSDF's."""
    def level(energy, spin):
        return Level(Value(energy, source="user"), float(spin), +1, f"{spin}+")

    return LevelScheme("62Ni", [level(0.0, 0), level(2301.8, 2), level(2336.5, 4), level(3500.0, 2)])


# -- rows -----------------------------------------------------------------------------------------------------------


def test_a_row_round_trips_and_its_id_is_its_content():
    row = b_e2(806, "nrf", "measured", plus=69, minus=59, detail="thin target", depends_on=["abc"])
    again = results.Row.from_dict(json.loads(json.dumps(row.to_dict())))
    assert again == row and again.id == row.id and len(row.id) == 16
    assert row.code == results.code() and set(row.code) == {"version", "commit"}
    assert row.created.endswith("+00:00") and row.depends_on == ("abc",)
    assert b_e2(807, "nrf", "measured", plus=69, minus=59, created=row.created).id != row.id
    assert (row.unc("plus"), row.unc("minus"), row.unc()) == (69, 59, 64)
    stat_sys = b_e2(918.6, "coulex-yield", "projected", stat=12.9, sys=63.4, truth=884.1,
                    origin=results.from_run("62Ni-on-194Pt", "03-beam-24h", {"gate": "inelastic"}))
    assert stat_sys.unc() == pytest.approx(64.70, abs=0.01) and stat_sys.relative == pytest.approx(0.0704, abs=1e-4)
    assert not stat_sys.independent and row.independent
    assert set(results.INDEPENDENT) == set(results.EVIDENCE) - {results.PROJECTED}


@pytest.mark.parametrize("change, message", [
    ({"method": "guess"}, "not a method"),
    ({"evidence": "rumour"}, "not a kind of evidence"),
    ({"observable": "B(E7)"}, "not an observable"),
    ({"evidence": "projected"}, "needs 'truth'"),
    ({"truth": 884.1}, "projected rows only"),
    ({"direction": None}, "direction stated"),
    ({"unit": "W.u."}, "always the ↓ value"),
    ({"unit": "ps"}, "not a unit of B"),
    ({"total": 5.0, "stat": 3.0}, "one form"),
    ({"plus": 5.0}, "both plus and minus"),
    ({"total": -1.0}, "not negative"),
    ({"origin": {"note": "from memory"}}, "origin needs a 'kind'"),
    ({"key": "62Ni:2+1"}, "as a transition"),
    ({"nuclide": "60Ni"}, "is not of 60Ni"),
    ({"observable": "tau", "unit": "ps", "key": "62Ni:2+1"}, "has no direction"),
])
def test_a_row_is_checked_when_it_is_made(change, message):
    good = {"nuclide": "62Ni", "observable": "B(E2)", "key": KEY, "value": 906.0, "unit": "e2fm4",
            "method": "coulex-fit", "evidence": "measured", "direction": "up", "origin": NDS}
    results.Row(**good)
    with pytest.raises(ValueError, match=message):
        results.Row(**{**good, **change})


# -- the files ------------------------------------------------------------------------------------------------------


def test_the_ledger_is_only_ever_appended_to(tmp_path):
    folder = tmp_path / "experiments" / "62Ni-on-194Pt"
    first = b_e2(918.6, "coulex-yield", "projected", stat=12.9, sys=63.4, truth=884.1,
                 origin=results.from_run(folder.name, "03-beam-24h"), created="2026-11-09T22:14:36+00:00")
    assert results.append(folder, first) and not results.append(folder, first)
    before = (folder / results.RESULTS_FILE).read_text(encoding="utf-8")
    second = first.revised(value=912.0, created="2026-11-10T08:00:00+00:00")
    assert second.supersedes == first.id and second.id != first.id and second.truth == first.truth
    assert results.append(folder, second)
    after = (folder / results.RESULTS_FILE).read_text(encoding="utf-8")
    assert after.startswith(before) and len(after.splitlines()) == 2
    rows = results.of_experiment(folder)
    assert rows == [first, second]
    assert results.current(rows) == [second] and results.history(rows, second) == [second, first]
    # A line that cannot be read is passed over and said, not fatal.
    with (folder / results.RESULTS_FILE).open("a", encoding="utf-8") as f:
        f.write('{"id": "x", "observable": "B(E9)"}\n')
    problems = []
    assert results.read(folder / results.RESULTS_FILE, problems) == rows and "line 3" in problems[0]


def test_a_nucleus_gathers_its_experiments_and_its_own_points(tmp_path, monkeypatch, no_ensdf):
    root, home = tmp_path / "experiments", tmp_path / "nuclei"
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(root))
    monkeypatch.delenv(results.ENV, raising=False)
    assert results.nuclei_home() == home
    lab = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=root)
    other = ExperimentFolder.create(Experiment.example("coulex_ni58"), root=root)
    ExperimentFolder.create(Experiment.example("alpha_on_gold"), root=root)  # studies no nucleus
    projected = b_e2(918.6, "coulex-yield", "projected", stat=12.9, sys=63.4, truth=884.1,
                     origin=results.from_run(lab.path.name, "03-beam-24h"), created="2026-11-09T22:14:36+00:00")
    evaluated = b_e2(884.1, "evaluation", "evaluated", total=13.1, created="2026-11-01T00:00:00+00:00")
    results.append(lab.path, projected)
    assert results.add_point(evaluated) and (home / "62Ni" / results.POINTS_FILE).is_file()
    assert results.of_nuclide("Ni-62") == [evaluated, projected]
    assert results.select(results.of_nuclide("62Ni"), "B(E2)", KEY, "evaluated") == [evaluated]
    assert results.of_nuclide("58Ni") == []
    found = {e["nuclide"]: e for e in results.nuclei()}
    assert list(found) == ["58Ni", "62Ni"]
    assert found["62Ni"] == {"nuclide": "62Ni", "rows": 2, "experiments": [lab.path.name]}
    assert found["58Ni"] == {"nuclide": "58Ni", "rows": 0, "experiments": [other.path.name]}
    monkeypatch.setenv(results.ENV, str(tmp_path / "elsewhere"))
    assert results.nuclei_home() == tmp_path / "elsewhere" and results.of_nuclide("62Ni") == [projected]


def test_backfill_makes_projected_rows_once(tmp_path, no_ensdf):
    lab = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=tmp_path)
    run = lab.path / "runs" / "03-beam-24h"
    run.mkdir()
    lab.setup().save(run / "setup.toml")
    (run / "summary.json").write_text(json.dumps({"number": 3, "kind": "beam"}), encoding="utf-8")
    analysis = {"settings": {"particle_gate": "inelastic"}, "normalisation": "rutherford", "gate": "inelastic",
                "b_e2fm4": 918.6, "statistical": 0.014, "systematic": 0.069, "truth_e2fm4": 884.105,
                "budget": {"efficiency": 0.05, "angular correlation": 0.03}, "saved": "2026-11-09T22:14:36+00:00"}
    (run / "analysis.json").write_text(json.dumps(analysis), encoding="utf-8")
    empty = lab.path / "runs" / "04-beam-10min"  # a run that was never analysed
    empty.mkdir()
    (empty / "summary.json").write_text("{}", encoding="utf-8")

    (row,) = results.backfill(lab.path)
    assert (row.nuclide, row.observable, row.unit, row.direction) == ("62Ni", "B(E2)", "e2fm4", "up")
    assert (row.method, row.evidence, row.truth) == ("coulex-yield", "projected", 884.105)
    assert row.key == "62Ni:2+@1173>0+@0:E2"  # unanchored: no ENSDF copy here
    assert row.stat == pytest.approx(12.86, abs=0.01) and row.sys == pytest.approx(63.38, abs=0.01)
    assert row.budget["efficiency"] == pytest.approx(45.93, abs=0.01)
    assert row.created == analysis["saved"] and row.code == {"version": None, "commit": None}
    assert row.origin["run"] == "03-beam-24h" and row.origin["settings"] == analysis["settings"]
    assert row.detail == "rutherford normalisation, inelastic gate"
    assert results.backfill(lab.path) == [] and results.of_experiment(lab.path) == [row]
    assert json.loads((run / "analysis.json").read_text(encoding="utf-8")) == analysis  # left as it was

    analysis.update(b_e2fm4=900.0, saved="2026-11-12T09:00:00+00:00")  # analysed again
    (run / "analysis.json").write_text(json.dumps(analysis), encoding="utf-8")
    (newer,) = results.backfill(lab.path)
    assert newer.supersedes == row.id and results.current(results.of_experiment(lab.path)) == [newer]
    assert results.backfill(lab.path) == []


# -- names ----------------------------------------------------------------------------------------------------------


def test_keys_read_back_from_their_text():
    for text in ("62Ni:2+1", "62Ni:2+@1173", "57Ni:3/2-(2)", "62Ni:?@2336", "62Ni:10+1"):
        assert results.LevelKey.parse(text).text == text
    key = results.LevelKey.parse("57Ni:3/2-2")
    assert (key.nuclide, key.spin, key.parity, key.ordinal, key.anchored) == ("57Ni", 1.5, -1, 2, True)
    assert not results.LevelKey.parse("62Ni:2+@1173").anchored and not results.LevelKey.parse("62Ni:2+(1)").anchored
    t = results.TransitionKey.parse(KEY)
    assert (t.text, t.initial.text, t.final.text, t.multipolarity, t.anchored) == (KEY, "62Ni:2+1", "62Ni:0+1", "E2",
                                                                                 True)
    for bad in ("62Ni", "62Ni:2+", "2+1", KEY):
        with pytest.raises(ValueError, match="as a level"):
            results.LevelKey.parse(bad)
    with pytest.raises(ValueError, match="as a transition"):
        results.TransitionKey.parse("62Ni:2+1>0+1")


def test_the_four_conventions_map_onto_one_key(ni62_levels):
    s = scheme()
    # The scheme's index and its level name: ENSDF's ordinal, however the scheme was cut.
    assert level_names(s.levels) == ["0+", "2+₁", "4+", "2+₂"]
    assert results.key_of_level(s, 1).text == "62Ni:2+2" == results.key_of_name(s, "2+₁").text
    assert results.key_of_level(s, 1).energy_kev == 2301.85 and results.key_of_level(s, 1).anchored
    assert results.key_of_level(s, 2).text == "62Ni:4+1"
    # A level ENSDF does not have keeps its energy and Jπ.
    lost = results.key_of_level(s, 3)
    assert (lost.text, lost.anchored, lost.ordinal) == ("62Ni:2+@3500", False, None)
    assert results.level_key("62Ni", 1175.5, 2, +1).text == "62Ni:2+@1176"  # more than 2 keV away
    assert results.level_key("62Ni", 1172.0, 3, -1).text == "62Ni:3-@1172"  # another Jπ at that energy
    assert results.level_key("Ni-62", 2048.0).text == "62Ni:0+2"  # by its energy alone
    # A shell-model state (2J, π, n), and one beyond the levels ENSDF has.
    assert results.key_of_shell_state("62Ni", (4, 1, 2)).text == "62Ni:2+2"
    assert results.key_of_shell_state("62Ni", (8, 1, 3)).text == "62Ni:4+(3)"
    # The [reaction] energy.
    assert results.key_of_reaction(Experiment.example("coulex_ni62_pt194")).text == KEY
    with pytest.raises(ValueError, match="not Coulomb excitation"):
        results.key_of_reaction(Experiment.example("alpha_on_gold"))
    # A transition: the upper level first, whichever way it is given.
    assert results.key_of_transition(s, 0, 1).text == "62Ni:2+2>0+1:E2" == results.key_of_transition(s, 1, 0).text
    with pytest.raises(ValueError, match="one nuclide"):
        results.transition_key(results.level_key("62Ni", 0.0), results.level_key("58Ni", 0.0))


def test_without_ensdf_every_level_is_unanchored(no_ensdf):
    s = scheme()
    assert [results.key_of_level(s, i).text for i in range(4)] == ["62Ni:0+@0", "62Ni:2+@2302", "62Ni:4+@2336",
                                                                  "62Ni:2+@3500"]
    assert results.key_of_shell_state("62Ni", (4, 1, 1)).text == "62Ni:2+(1)"
    assert results.key_of_reaction(Experiment.example("coulex_ni62_pt194")).text == "62Ni:2+@1173>0+@0:E2"


@needs_ensdf
def test_levels_are_anchored_to_the_local_ensdf_copy():
    s = scheme()
    assert [results.key_of_level(s, i).text for i in range(3)] == ["62Ni:0+1", "62Ni:2+2", "62Ni:4+1"]
    assert results.key_of_reaction(Experiment.example("coulex_ni62_pt194")).text == KEY


# -- the registry ---------------------------------------------------------------------------------------------------


def test_the_registry_names_conventions_routes_and_conversions():
    assert set(observables.OBSERVABLES) == {"B(E2)", "tau", "Qs", "B(M1)", "mu", "E_level"}
    b = observables.get("B(E2)")
    assert (b.belongs_to, b.unit, b.direction, b.symbol) == ("transition", "e2fm4", "up", "B(E2)↑")
    assert [observables.get(k).belongs_to for k in ("tau", "Qs", "E_level")] == ["level"] * 3
    kinds = {k.key: k.live for k in KINDS}
    for o in observables.OBSERVABLES.values():
        assert o.routes and o.used_for and all(u.text and u.reference for u in o.used_for)
        for route in o.routes:
            assert route.method in results.METHODS
            assert route.live is (None if route.kind is None else kinds[route.kind])
        for conversion in o.conversions:
            assert callable(getattr(observables, conversion.function))
    assert {r.kind for r in observables.routes("B(E2)", live=True)} == {"coulex", "lifetimes"}
    assert {r.kind for r in observables.routes("B(E2)", live=False)} == {"inelastic"}
    # The conversions are the existing functions: B↑ = 5 B↓ for 0⁺ → 2⁺, W.u. always of B↓, B ↔ τ.
    assert observables.b_up_over_down(0, 2) == pytest.approx(5.0)
    wu = observables.weisskopf(884.1 / 5, "62Ni")
    assert wu == pytest.approx(12.13, abs=0.01)
    assert observables.to_canonical("B(E2)", wu, "W.u.", "down", nuclide="62Ni", j_lower=0,
                                    j_upper=2) == pytest.approx(884.1)
    assert observables.to_canonical("B(E2)", 0.08841, "e2b2", "up") == pytest.approx(884.1)
    assert observables.to_canonical("tau", 2080.0, "fs") == pytest.approx(2.08)
    assert observables.matrix_element(884.1, 0) ** 2 == pytest.approx(884.1)
    assert observables.quadrupole_moment(observables.diagonal_matrix_element(-27.5, 2), 2) == pytest.approx(-27.5)
    with pytest.raises(ValueError, match="always the ↓ value"):
        observables.canonical_factor("B(E2)", "W.u.", "up", nuclide="62Ni")
    with pytest.raises(ValueError, match="spins of both levels"):
        observables.canonical_factor("B(E2)", "e2fm4", "down")


# -- the comparison -------------------------------------------------------------------------------------------------


def ni62_rows():
    evaluated = b_e2(884.1, "evaluation", "evaluated", total=13.1)
    return [
        evaluated,
        b_e2(906, "coulex-fit", "measured", total=37, origin=results.reference("Allmond et al.", 2014)),
        b_e2(806, "nrf", "measured", plus=69, minus=59),
        b_e2(408, "shell-model", "theory", origin=results.shell_model_job("jun45", "f5pg9", {"p": 1.5, "n": 0.5})),
        b_e2(918.6, "coulex-yield", "projected", stat=12.9, sys=63.4, truth=884.1,
             origin=results.from_run("62Ni-on-194Pt", "03-beam-24h")),
    ]


def test_the_62ni_values_beside_the_evaluation():
    rows = ni62_rows()
    c = compare.compare(rows)
    assert (c.observable, c.label, c.unit, c.key) == ("B(E2)", "B(E2)↑", "e2fm4", KEY)
    reference, coulex, nrf, theory = c.entries
    assert c.reference is reference and reference.is_reference and reference.pull is None
    assert reference.uncertainty_percent == pytest.approx(1.48, abs=0.01)
    assert (coulex.pull, coulex.verdict) == (pytest.approx(0.56, abs=0.005), "agree")
    assert coulex.difference_percent == pytest.approx(2.48, abs=0.01) and coulex.independent_of_reference
    # The asymmetric uncertainty counts with the side that faces the reference: +69, the value lying below it.
    assert (nrf.pull, nrf.verdict) == (pytest.approx(-1.11, abs=0.005), "agree")
    assert (nrf.plus, nrf.minus) == (69, 59)
    # No uncertainty: a ratio and no verdict.
    assert theory.ratio == pytest.approx(0.46, abs=0.005) and (theory.pull, theory.verdict) == (None, None)
    # The projected row is closed against its truth and set beside what is known; it is not among the entries.
    (projected,) = c.projected
    assert projected.closure == pytest.approx(0.53, abs=0.005) and projected.closure_verdict == "agree"
    assert projected.reach_percent == pytest.approx(7.0, abs=0.05)
    assert projected.known_percent == pytest.approx(1.48, abs=0.01)
    assert projected.best_measured_percent == pytest.approx(4.08, abs=0.01) and c.best_measured is coulex
    assert not hasattr(projected, "pull")
    plain = c.to_dict()
    assert json.loads(json.dumps(plain)) == plain and plain["entries"][1]["row"]["id"] == rows[1].id


def test_the_rule_the_reference_and_what_is_not_independent():
    assert [compare.verdict(p) for p in (0.0, -1.99, 2.0, -3.0, 3.01, None)] == ["agree", "agree", "tension",
                                                                                "tension", "differ", None]
    rows = ni62_rows()
    evaluated, coulex = rows[0], rows[1]
    adopting = b_e2(884.1, "evaluation", "evaluated", total=13.1, depends_on=[coulex.id],
                    created="2026-11-20T00:00:00+00:00")
    low = b_e2(745, "rdds", "measured", total=30)
    c = compare.compare(rows + [evaluated.revised(value=884.1, depends_on=[coulex.id]), low])
    by_id = {e.id: e for e in c.entries}
    assert evaluated.id not in by_id  # superseded: only what holds now is compared
    assert by_id[coulex.id].independent_of_reference is False and by_id[rows[2].id].independent_of_reference
    assert by_id[low.id].verdict == "differ" and by_id[low.id].pull == pytest.approx(-4.25, abs=0.01)
    # Any row can be the reference, but not a projected one; without one the rows are only listed.
    c = compare.compare(rows, reference=coulex)
    assert c.reference.id == coulex.id and c.entries[0].pull == pytest.approx(-0.56, abs=0.005)
    assert c.projected[0].known_percent == pytest.approx(4.08, abs=0.01)
    with pytest.raises(ValueError, match="cannot be the reference"):
        compare.compare(rows, reference=rows[4])
    listed = compare.compare(rows[1:])
    assert listed.reference is None and all(e.pull is None and e.ratio is None for e in listed.entries)
    assert listed.projected[0].known_percent is None and listed.projected[0].closure is not None
    assert adopting.depends_on == (coulex.id,)
    with pytest.raises(ValueError, match="different observables or of different levels"):
        compare.compare(rows + [results.Row("62Ni", "B(E2)", "62Ni:2+2>0+1:E2", 20.0, "e2fm4", "evaluation",
                                            "evaluated", direction="up", origin=NDS)])


def test_rows_in_other_conventions_and_a_lifetime_join_the_comparison():
    evaluated = b_e2(884.1, "evaluation", "evaluated", total=13.1)
    # B↓ in Weisskopf units: 12.13 W.u. is 884.1 e²fm⁴ ↑, and the comparison says how it got there.
    wu = results.Row("62Ni", "B(E2)", KEY, 12.13, "W.u.", "evaluation", "measured", direction="down", total=0.5,
                     origin=NDS)
    entry = compare.compare([evaluated, wu]).entries[1]
    assert entry.value == pytest.approx(884.1, rel=1e-3) and entry.uncertainty_percent == pytest.approx(4.12, abs=0.01)
    assert "W.u." in entry.converted and "5 ×" in entry.converted
    # τ of the 2⁺₁ level becomes B(E2)↑ through the existing lifetimes.b_from_tau.
    tau_ps = lifetimes.tau_from_b(1172.91, "E2", 884.1 / 5)
    tau = results.Row("62Ni", "tau", "62Ni:2+1", 1e3 * tau_ps, "fs", "rdds", "measured", plus=100.0, minus=50.0,
                      origin=results.own_measurement("2026-11-24", "a plunger run"))
    b = compare.convert(tau, "B(E2)", KEY, energy_kev=1172.91)
    assert (b.observable, b.key, b.unit, b.direction, b.method, b.evidence) == ("B(E2)", KEY, "e2fm4", "up",
                                                                               "converted", "measured")
    assert b.value == pytest.approx(884.1) and b.depends_on == (tau.id,) and b.origin["from"] == tau.id
    # B falls when τ rises: the sides swap, the relative sizes stay.
    assert b.plus / b.value == pytest.approx(50.0 / tau.value) and b.minus / b.value == pytest.approx(100.0 / tau.value)
    back = compare.convert(b, "tau", energy_kev=1172.91)
    assert (back.key, back.unit) == ("62Ni:2+1", "ps") and back.value == pytest.approx(tau_ps)
    assert back.depends_on == (b.id, tau.id)
    joined = compare.compare([evaluated, b]).entries[1]
    assert joined.pull == pytest.approx(0.0, abs=1e-6) and joined.verdict == "agree"
    with pytest.raises(ValueError, match="does not start at the level"):
        compare.convert(tau, "B(E2)", "62Ni:2+2>0+1:E2", energy_kev=2301.85)
    with pytest.raises(ValueError, match="no conversion"):
        compare.convert(tau, "Qs", energy_kev=1.0)
