"""The nucleus above the experiments (wave 1 of the results pages): the evaluation as rows, data points that are
corrected and withdrawn but never deleted, and the nucleus views on a temporary folder with the ⁶²Ni numbers of
the mock-ups. Plumbing on synthetic rows: no run is simulated here."""

import json
from pathlib import Path

import pytest

from physim.nuclear import Experiment, compare, ensdf, lifetimes, logbook, nucleus, results
from physim.nuclear.levels import LevelScheme
from physim.nuclear.runs import ExperimentFolder

needs_ensdf = pytest.mark.skipif(not ensdf.available(), reason="no local copy of ENSDF (scripts/fetch_ensdf.py)")

KEY = "62Ni:2+1>0+1:E2"
#: The first adopted levels of ⁶²Ni as (energy keV, 2J, parity, ordinal), for anchoring without a local copy.
NI62 = ((0.0, 0, 1, 1), (1172.91, 4, 1, 1), (2048.55, 0, 1, 2), (2301.85, 4, 1, 2), (2336.46, 8, 1, 1))
ENSDF_ORIGIN = results.reference("ENSDF, adopted levels and γ rays", 2025, where="NDS 204, 1 (2025)",
                                 evaluation=nucleus.ENSDF, dataset_date="202507")
ALLMOND = {"authors": "Allmond et al.", "year": 2014, "doi": "10.1103/PhysRevC.90.034309", "where": "Table II"}


def card(kind, cols=(), cont=" "):
    """One 80-column ENSDF record of ⁶²Ni: ``cols`` are (first column counting from 1, text)."""
    line = list(f" 62NI{cont} {kind} ".ljust(80))
    for start, text in cols:
        line[start - 1:start - 1 + len(text)] = text
    return "".join(line)


def level_card(e, de, j, t="", dt=""):
    return card("L", [(10, e), (20, de), (22, j), (40, t), (50, dt)])


def gamma_card(e, ri="100"):
    return card("G", [(10, e), (22, ri), (32, "E2")])


#: A small adopted data set in ENSDF's format: a stated B(E2), one to derive, a half-life that is only a limit, and
#: a level with a γ ray that fits no level.
TINY_ENSDF = [
    card(" ", [(10, "ADOPTED LEVELS, GAMMAS"), (75, "202507")]),
    card("H", [(10, "TYP=FUL$CIT=NDS 204, 1 (2025)$")]),
    level_card("0.0", "", "0+", "STABLE"),
    level_card("1172.91", "3", "2+", "1.442 PS", "21"),
    card("L", [(10, "MOME2=+0.05 12")], cont="2"),
    gamma_card("1172.89"),
    card("G", [(10, "BE2W=12.13 18")], cont="2"),
    level_card("2048.55", "8", "0+", "1.40 PS", "51"),
    gamma_card("875.6"),
    level_card("2890.63", "15", "0+", "3.1 PS", "GT"),
    gamma_card("1717.8"),
    level_card("3500.0", "5", "2+", "0.50 PS", "10"),
    gamma_card("2327.1"),
    gamma_card("999.9", ri="20"),
]


@pytest.fixture
def lab(tmp_path, monkeypatch):
    """Empty experiments and nuclei folders, no local ENSDF copy, and the ⁶²Ni levels to anchor to: the same on
    every computer. The 2⁺₁ level decays by its one γ ray, as ENSDF has it."""
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path / "experiments"))
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setattr(results, "ensdf_levels", lambda nuclide: NI62 if nuclide == "62Ni" else ())

    def decay(nuclide, level_key):
        if level_key.text != "62Ni:2+1":
            return None, None, "no γ ray known"
        return (results.TransitionKey.parse(KEY),
                {"energy_kev": 1172.9, "branching": 1.0, "conversion": 0.0, "mixing": None}, "")

    monkeypatch.setattr(nucleus, "decay", decay)
    return tmp_path


def evaluated(b_wu=12.13):
    """A small evaluation of ⁶²Ni, as :func:`nucleus.evaluation` makes it: the B(E2) depends on the τ row."""
    def row(observable, key, value, unit, **more):
        return results.Row("62Ni", observable, key, value, unit, "evaluation", "evaluated", origin=ENSDF_ORIGIN,
                           **more)

    tau = row("tau", "62Ni:2+1", 2.080, "ps", total=0.030)
    return [row("E_level", "62Ni:2+1", 1172.91, "keV", total=0.03), tau,
            row("Qs", "62Ni:2+1", 0.05, "eb", total=0.12),
            row("B(E2)", KEY, b_wu, "W.u.", direction="down", total=0.18, depends_on=(tau.id,))]


# -- the evaluation as rows ---------------------------------------------------------------------------------------


def test_limits_and_unplaced_gamma_rays_are_not_made_into_values(tmp_path, monkeypatch):
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "ensdf"))
    (tmp_path / "ensdf").mkdir()
    (tmp_path / "ensdf" / "ensdf.062").write_text("\n".join(TINY_ENSDF), encoding="latin-1")
    assert [lv.dropped for lv in ensdf.adopted("62Ni").levels] == [0, 0, 0, 0, 1]
    by = {(r.observable, r.key): r for r in nucleus.evaluation("62Ni")}
    assert sorted(k for o, k in by if o == "E_level") == ["62Ni:0+2", "62Ni:0+3", "62Ni:2+1", "62Ni:2+2"]
    # "3.1 PS GT" is a limit: no τ of 0+3 and no B(E2) from it. Every τ row has its uncertainty.
    assert sorted(k for o, k in by if o == "tau") == ["62Ni:0+2", "62Ni:2+1", "62Ni:2+2"]
    assert all(r.total > 0 for r in by.values() if r.observable in ("tau", "B(E2)"))
    assert by["tau", "62Ni:2+1"].value == pytest.approx(2.0804, abs=1e-4)
    assert by["tau", "62Ni:2+1"].total == pytest.approx(0.0303, abs=1e-4)
    moment = by["Qs", "62Ni:2+1"]
    assert (moment.value, moment.unit, moment.total) == (0.05, "eb", 0.12)
    # B(E2): ENSDF's own where it states one, the level scheme's where the half-life and one γ ray give it; none
    # for the level whose second γ ray fits no level (its branching is not complete).
    assert sorted(k for o, k in by if o == "B(E2)") == ["62Ni:0+2>2+1:E2", KEY]
    stated, derived = by["B(E2)", KEY], by["B(E2)", "62Ni:0+2>2+1:E2"]
    assert (stated.value, stated.total, stated.unit, stated.direction) == (12.13, 0.18, "W.u.", "down")
    scheme = LevelScheme.from_adopted(ensdf.adopted("62Ni"))
    scheme.derive()
    assert derived.value == pytest.approx(scheme.b(2, 1, "E2")) and derived.unit == "e2fm4"
    assert derived.depends_on == (by["tau", "62Ni:0+2"].id,) and "half-life" in derived.detail
    assert all(r.origin["dataset_date"] == "202507" and "NDS 204" in nucleus.source(r) for r in by.values())

    # B(E2) ↔ τ goes through a level's one γ ray; a level with a γ ray that could not be placed is refused.
    key, inputs, why = nucleus.decay("62Ni", nucleus.level("62Ni", "2+1"))
    assert (key.text, inputs["energy_kev"], inputs["branching"], why) == (KEY, 1172.89, 1.0, "")
    assert nucleus.decay("62Ni", nucleus.level("62Ni", "2+2"))[1:] == (
        None, "ENSDF lists 1 more γ ray from the level 2+₂ that could not be placed: its branching is not known")
    stuck = results.Row("62Ni", "tau", "62Ni:2+2", 0.7, "ps", "dsam", "measured", total=0.1,
                        origin=results.reference("Someone", 2001))
    assert "could not be placed" in nucleus.gives("62Ni", stuck)[0]["missing"]
    assert nucleus.with_converted("62Ni", [stuck]) == ([stuck], {})
    # A converted row is the same row each time it is made: its source's time and code, so the same id.
    kenn = results.Row("62Ni", "tau", "62Ni:2+1", 2.01, "ps", "dsam", "measured", total=0.12,
                       origin=results.reference("Kenn et al.", 2001), created="2001-05-01T00:00:00+00:00")
    (_, first), (_, again) = (nucleus.with_converted("62Ni", [kenn])[0] for _ in range(2))
    assert (first.key, first.created, first.code, first.depends_on) == (KEY, kenn.created, kenn.code, (kenn.id,))
    assert first.id == again.id == compare.convert(kenn, "B(E2)", key, **inputs).id
    # A lifetime analysis writes its own B(E2) beside its τ (the B depends on the τ): neither is converted again.
    run = results.from_run("62Ni-plunger", "01-beam-24h")
    tau = results.Row("62Ni", "tau", "62Ni:2+1", 2.0, "ps", "rdds", "projected", stat=0.1, truth=2.08, origin=run)
    b = results.Row("62Ni", "B(E2)", KEY, 180.0, "e2fm4", "rdds", "projected", direction="down", total=9.0,
                    truth=176.8, depends_on=(tau.id,), origin=run)
    assert nucleus.with_converted("62Ni", [tau, b]) == ([tau, b], {})


@needs_ensdf
def test_the_adopted_values_come_as_evaluated_rows(tmp_path):
    rows = nucleus.evaluation("Ni-62", max_energy_kev=2400)
    assert rows and all((r.evidence, r.method, r.nuclide) == ("evaluated", "evaluation", "62Ni") for r in rows)
    assert all(nucleus.is_evaluation(r) and r.origin["dataset_date"] for r in rows)
    assert all("ENSDF" in nucleus.source(r) for r in rows)
    by = {(r.observable, r.key): r for r in rows}
    assert len(by) == len(rows) and {r.observable for r in rows} == {"E_level", "tau", "Qs", "B(E2)", "B(M1)", "mu"}
    energy, tau, b = by["E_level", "62Ni:2+1"], by["tau", "62Ni:2+1"], by["B(E2)", KEY]
    assert energy.value == pytest.approx(1172.9, abs=0.5) and 0 < energy.total < 1
    # τ is what the lifetime code makes of ENSDF's half-life, and B(E2) is ENSDF's own number in Weisskopf units.
    scheme = LevelScheme.from_adopted(ensdf.adopted("62Ni"))
    assert tau.value == lifetimes.level_taus(scheme)["2+₁"][0] and tau.unit == "ps" and tau.total > 0
    line = next(t for t in scheme.transitions if (t.initial, t.final) == (1, 0))
    given = line.weisskopf["E2"]
    assert (b.value, b.total, b.unit, b.direction) == (given.value, given.unc, "W.u.", "down")
    assert ("E_level", "62Ni:0+1") not in by and max(r.value for r in rows if r.observable == "E_level") <= 2400
    assert all(r.total is not None for r in rows if r.observable == "tau")  # no limit is made into a value
    home = tmp_path / "nuclei"
    first = nucleus.snapshot_evaluation("62Ni", home, max_energy_kev=2400)
    again = nucleus.snapshot_evaluation("62Ni", home, max_energy_kev=2400)
    assert (len(first["added"]), first["unchanged"]) == (len(rows), 0)
    assert (again["added"], again["superseded"], again["unchanged"]) == ([], [], len(rows))


def test_the_snapshot_is_taken_once_and_a_changed_evaluation_supersedes(lab, monkeypatch):
    monkeypatch.setattr(nucleus, "evaluation", lambda nuclide, max_energy_kev=None: evaluated())
    first = nucleus.snapshot_evaluation("62Ni")
    file = lab / "nuclei" / "62Ni" / results.POINTS_FILE
    before = file.read_text(encoding="utf-8")
    assert len(first["added"]) == 4 and first["reference"] == "ENSDF adopted, NDS 204, 1 (2025)"
    again = nucleus.snapshot_evaluation("62Ni")  # rows made anew (new ids, new times): nothing is written
    assert (again["added"], again["superseded"], again["unchanged"]) == ([], [], 4)
    assert file.read_text(encoding="utf-8") == before
    held = results.of_nuclide("62Ni")
    b = results.select(held, "B(E2)", KEY)[0]
    assert b.depends_on == (results.select(held, "tau")[0].id,)
    # A newer ENSDF copy with another B(E2): one new row, which supersedes the old one; the old one stays.
    monkeypatch.setattr(nucleus, "evaluation", lambda nuclide, max_energy_kev=None: evaluated(b_wu=12.4))
    newer = nucleus.snapshot_evaluation("62Ni")
    assert (len(newer["added"]), len(newer["superseded"]), newer["unchanged"]) == (0, 1, 3)
    assert newer["superseded"][0].supersedes == b.id and file.read_text(encoding="utf-8").startswith(before)
    held = results.of_nuclide("62Ni")
    assert [r.value for r in results.select(results.current(held), "B(E2)", KEY)] == [12.4]
    assert [h["status"] for h in nucleus.history(held, "B(E2)", KEY)] == ["holds", "superseded"]


# -- data points --------------------------------------------------------------------------------------------------


def test_a_level_and_a_transition_as_one_says_them(lab):
    for said in ("2+1", "2⁺₁", "62Ni:2+1", "2+ at 1172.9 keV", "1.1729 MeV 2+", (1172.9, "2+"), 1172.9,
                 {"energy": "1173 keV", "jpi": "2+"}):
        assert nucleus.level("62Ni", said).text == "62Ni:2+1", said
    assert nucleus.level("62Ni", "gs").text == "62Ni:0+1" and nucleus.level("62Ni", "0+ gs").anchored
    assert nucleus.level("62Ni", "3- at 3757 keV").text == "62Ni:3-@3757"  # not among the levels: unanchored
    for said in ("2+1 > 0+1", "0+1 -> 2+1", "gs → 2+ at 1173 keV", "2+1 to gs", KEY, {"from": "2+1", "to": "gs"}):
        assert nucleus.transition("62Ni", said).text == KEY, said
    assert nucleus.label("B(E2)", KEY) == "B(E2; 0+₁ → 2+₁)↑" and nucleus.label("tau", "62Ni:2+1") == "τ(2+₁)"
    for said, message in (("2+9", "cannot be placed"), ("the first one", "cannot read"),
                          ("60Ni:2+1", "not of 62Ni"), ("2+ at -5 keV", "cannot be negative"),
                          ("x+ at 100 keV", "spin and parity")):
        with pytest.raises(ValueError, match=message):
            nucleus.level("62Ni", said)
    with pytest.raises(ValueError, match="two levels"):
        nucleus.transition("62Ni", "2+1")
    assert nucleus.level("62Ni", "2+9", check=False).text == "62Ni:2+9"  # for looking rows up


def test_a_data_point_is_added_corrected_and_withdrawn_but_never_deleted(lab):
    own = {"date": "2027-06-14", "note": "iThemba LABS, runs 12 to 31"}
    point = nucleus.add_data_point("62Ni", "B(E2)", "0+ gs -> 2+ at 1173 keV", 0.0906, "e2b2", direction="up",
                                   stat=0.0020, sys=0.0031, method="coulex-yield", own=own)
    assert (point.key, point.evidence, point.origin["kind"]) == (KEY, "measured", "own measurement")
    assert nucleus.cell(point)["text"] == "906 ± 37 e2fm4↑" and nucleus.is_data_point(point)
    file = lab / "nuclei" / "62Ni" / results.POINTS_FILE
    with pytest.raises(nucleus.PointProblem, match="in the ledger of 62Ni already"):
        nucleus.add_data_point("62Ni", "B(E2)", KEY, 0.0906, "e2b2", direction="up", stat=0.0020, sys=0.0031,
                               method="coulex-yield", own=own)
    # A correction is a new row that supersedes the old one: by naming it, or by saying what changes.
    fixed = nucleus.add_data_point("62Ni", "B(E2)", KEY, 0.0916, "e2b2", direction="up", stat=0.0020, sys=0.0031,
                                   method="coulex-yield", own=own, corrects=point.id[:8])
    again = nucleus.correct_data_point("62Ni", fixed.id, uncertainty=0.0040, detail="thick target")
    assert fixed.supersedes == point.id and again.supersedes == fixed.id
    assert (again.value, again.total, again.stat, again.origin) == (0.0916, 0.0040, None, fixed.origin)
    assert results.current(results.of_nuclide("62Ni")) == [again]
    with pytest.raises(ValueError, match="replaced by " + fixed.id):
        nucleus.correct_data_point("62Ni", point.id, value=1.0)
    # A retraction too: the point leaves what holds now, and every row is still in the file.
    gone = nucleus.retract_data_point("62Ni", again.id, "entered for the wrong nucleus")
    assert gone.supersedes == again.id and gone.retracted and results.current(results.of_nuclide("62Ni")) == []
    assert len(file.read_text(encoding="utf-8").splitlines()) == 4
    assert [h["status"] for h in nucleus.history(results.of_nuclide("62Ni"))] == \
        ["retraction", "retracted", "superseded", "superseded"]
    with pytest.raises(nucleus.PointProblem, match="needs its reason"):
        nucleus.retract_data_point("62Ni", gone.id, " ")
    with pytest.raises(ValueError, match="withdrawn by"):
        nucleus.retract_data_point("62Ni", again.id, "twice")
    # Only data points are corrected here: the evaluation is replaced by a new snapshot.
    energy = evaluated()[0]
    results.add_point(energy)
    with pytest.raises(ValueError, match="not a data point"):
        nucleus.retract_data_point("62Ni", energy.id, "no")


def test_a_bad_data_point_is_refused_with_everything_that_is_wrong(lab):
    with pytest.raises(nucleus.PointProblem) as refused:
        nucleus.add_data_point("62Ni", "B(E2)", "2+1", "many", "barn", uncertainty=-1, plus=2, method="guess",
                               evidence="guessed")
    said = "; ".join(refused.value.problems)
    for part in ("a transition is two levels", "the value must be a number", "'barn' is not a unit of B(E2)",
                 "needs its direction", "one form", "'uncertainty' must be a number that is not negative",
                 "a data point is measured, outside program, theory or evaluated, not 'guessed'", "the method is one of",
                 "where the point is from"):
        assert part in said, part
    assert len(refused.value.problems) == 9 and not (lab / "nuclei").exists()
    good = {"direction": "up", "uncertainty": 37.0, "method": "coulex-fit", "reference": ALLMOND}
    for change, message in (
            ({"uncertainty": None}, "a measured point needs its uncertainty"),
            ({"method": "shell-model"}, "not a method of a point that is measured"),
            ({"method": "evaluation", "evidence": "evaluated", "reference": {
                "authors": "B. Singh et al.", "year": 2025, "where": "Nucl. Data Sheets 204, 1 (2025)"}},
             "with the snapshot of the evaluation"),
            ({"method": "evaluation", "evidence": "evaluated", "reference": None, "own": {"date": "2027"}},
             "an evaluation is a published one"),
            ({"evidence": "projected"}, "simulated run"),
            ({"unit": "W.u."}, "always the ↓ value"),
            ({"reference": {"authors": "", "year": 1066}}, "needs its authors.*year it was published"),
            ({"reference": None, "own": {"date": "last June"}}, "written as 2027-06-14"),
            ({"own": {"date": "2027-06"}}, "one of the two"),
            ({"corrects": "abc"}, "no row 'abc'")):
        with pytest.raises(nucleus.PointProblem, match=message):
            nucleus.add_data_point("62Ni", "B(E2)", KEY, 906.0, **{**good, **change})
    for value, change, message in (
            (-5.0, {}, r"B\(E2\) cannot be negative, got -5"), (True, {}, "must be a number, got True"),
            (906.0, {"uncertainty": 0}, "an uncertainty above 0"),
            (906.0, {"uncertainty": True}, "'uncertainty' must be"),
            (906.0, {"uncertainty": None, "stat": 0.0, "sys": 0}, "an uncertainty above 0")):
        with pytest.raises(nucleus.PointProblem, match=message):
            nucleus.add_data_point("62Ni", "B(E2)", KEY, value, **{**good, **change})
    with pytest.raises(nucleus.PointProblem, match="τ must be above 0, got 0"):
        nucleus.add_data_point("62Ni", "tau", "2+1", 0.0, uncertainty=0.1, method="dsam", reference=ALLMOND)
    assert not (lab / "nuclei").exists()
    theory = nucleus.add_data_point("62Ni", "Qs", "2+1", -6.7, method="shell-model", evidence="theory",
                                    reference={"authors": "Honma et al.", "year": 2009})
    assert (theory.unit, theory.direction, theory.unc()) == ("efm2", None, None)  # theory may have no uncertainty


def test_a_row_that_was_not_written_is_not_returned(lab, monkeypatch):
    monkeypatch.setattr(results, "_now", lambda: "2026-11-24T12:00:00+00:00")  # everything in one second
    point = {"direction": "up", "uncertainty": 37.0, "method": "coulex-fit", "reference": ALLMOND}
    first = nucleus.add_data_point("62Ni", "B(E2)", KEY, 906.0, **point)
    nucleus.retract_data_point("62Ni", first.id, "entered twice")
    with pytest.raises(nucleus.PointProblem, match="written this very second"):  # it would have the same id
        nucleus.add_data_point("62Ni", "B(E2)", KEY, 906.0, **point)
    assert results.current(results.of_nuclide("62Ni")) == []


def test_the_backfill_keeps_to_its_method_and_says_nothing_twice(tmp_path, monkeypatch):
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    folder = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=tmp_path)
    run = folder.path / "runs" / "03-beam-24h"
    run.mkdir()
    folder.setup().save(run / "setup.toml")
    (run / "summary.json").write_text("{}", encoding="utf-8")
    analysis = {"settings": {"particle_gate": "inelastic"}, "b_e2fm4": 918.6, "statistical": 0.014,
                "systematic": 0.069, "truth_e2fm4": 884.105, "budget": {"efficiency": 0.05},
                "saved": "2026-11-09T22:14:36+00:00"}

    def analysed(**changes):
        analysis.update(changes)
        (run / "analysis.json").write_text(json.dumps(analysis), encoding="utf-8")
        return results.backfill(folder.path)

    (first,) = analysed()
    # The all-orders fit of the same run has a row of its own, by another method.
    fit = results.Row("62Ni", "B(E2)", first.key, 905.0, "e2fm4", "coulex-fit", "projected", direction="up",
                      stat=10.0, sys=40.0, truth=884.105, origin=results.from_run(folder.path.name, run.name),
                      created="2026-11-10T00:00:00+00:00")
    results.append(folder.path, fit)
    # Analysed again with the same settings and result: only the stamp is new, and nothing is added.
    assert analysed(saved="2026-11-11T08:00:00+00:00") == []
    # Analysed again with another result: the yield's row is superseded, the fit's is left alone.
    (newer,) = analysed(b_e2fm4=900.0, saved="2026-11-12T09:00:00+00:00")
    assert newer.supersedes == first.id and newer.method == "coulex-yield"
    assert results.current(results.of_experiment(folder.path)) == [fit, newer]
    assert analysed(settings={"particle_gate": "all"}, saved="2026-11-13T09:00:00+00:00")[0].supersedes == newer.id


# -- the views ----------------------------------------------------------------------------------------------------


def test_the_nucleus_an_experiment_studies():
    assert results.studied_nuclide(Experiment.example("coulex_ni62_pt194")) == "62Ni"
    assert nucleus.partner_nuclide(Experiment.example("coulex_ni62_pt194")) == "194Pt"
    assert results.studied_nuclide(Experiment.example("fusion_cr50_plunger")) == "50Cr"  # [run] channel = "4p2n"
    assert results.studied_nuclide(Experiment.example("fusion_dy158")) is None  # no studied channel marked
    assert results.studied_nuclide(Experiment.example("alpha_on_gold")) is None


def test_the_nuclei_and_one_nucleus_as_plain_data(lab):
    root = lab / "experiments"
    coulex = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=root)
    fusion = ExperimentFolder.create(Experiment.example("fusion_dy158"), root=root)
    ExperimentFolder.create(Experiment.example("alpha_on_gold"), root=root)
    for run in ("01-beam-1h", "03-beam-24h"):
        (coulex.path / "runs" / run).mkdir()
        (coulex.path / "runs" / run / "summary.json").write_text("{}", encoding="utf-8")
    # The ⁶²Ni numbers of the mock-up: the evaluation, a Coulex and a DSAM point, the shell model, a simulated run.
    projected = results.Row("62Ni", "B(E2)", KEY, 918.6, "e2fm4", "coulex-yield", "projected", direction="up",
                            stat=12.9, sys=63.4, truth=884.1,
                            origin=results.from_run(coulex.path.name, "03-beam-24h"))
    results.append(coulex.path, projected)
    adopted = results.Row("62Ni", "B(E2)", KEY, 884.1, "e2fm4", "evaluation", "evaluated", direction="up",
                          total=13.1, origin=ENSDF_ORIGIN)
    jun45 = results.Row("62Ni", "B(E2)", KEY, 408.0, "e2fm4", "shell-model", "theory", direction="up",
                        origin=results.shell_model_job("JUN45", "f5pg9", {"p": 1.5, "n": 1.1}))
    results.add_points([adopted, jun45])
    allmond = nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 906.0, direction="up", uncertainty=37.0,
                                     method="coulex-fit", reference=ALLMOND, detail="GOSIA")
    kenn = nucleus.add_data_point("62Ni", "tau", "2+1", 2.01, "ps", uncertainty=0.12, method="dsam",
                                  reference={"authors": "Kenn et al.", "year": 2001})

    found = nucleus.nuclei()
    json.dumps(found)
    assert [e["nuclide"] for e in found["nuclei"]] == ["62Ni", "194Pt"]
    assert {u["folder"]: u["reason"].split()[0:2] for u in found["unassigned"]} == {
        fusion.path.name: ["the", "setup"], "4He-on-197Au": ["elastic", "scattering"]}
    assert "[run] channel" in next(u["reason"] for u in found["unassigned"] if u["folder"] == fusion.path.name)
    ni, pt = found["nuclei"]
    [mine] = ni["experiments"]
    assert [mine[k] for k in ("folder", "role", "kind", "runs", "analysed_runs", "rows")] == [
        coulex.path.name, "studied", "coulex", 2, 1, 1]
    assert mine["latest"]["text"] == "919 ± 65 e2fm4↑"
    assert (ni["runs"], ni["analysed_runs"], ni["rows"], ni["points"], ni["evaluation"]["rows"]) == (2, 1, 5, 2, 1)
    assert [(x["role"], x["rows"], x["latest"]) for x in pt["experiments"]] == [("partner", 0, None)]
    assert (pt["rows"], pt["table"], pt["open"]) == (0, [], [])

    # One row per observable and level or transition; the two groups of evidence in separate entries.
    b, tau = ni["table"]
    assert (b["observable"], b["key"], b["label"], b["unit"]) == ("B(E2)", KEY, "B(E2; 0+₁ → 2+₁)↑", "e2fm4")
    assert b["evaluated"]["id"] == adopted.id and b["evaluated"]["text"] == "884 ± 13 e2fm4↑"
    assert (b["measured"]["count"], b["measured"]["best"]["id"]) == (2, allmond.id)  # with the B(E2) the τ gives
    assert b["measured"]["max"] == pytest.approx(914.6, abs=0.1)
    assert (b["theory"]["count"], b["theory"]["newest"]["value"], b["outside_program"]) == (1, 408.0, None)
    assert b["projected"]["best"]["uncertainty_percent"] == pytest.approx(7.04, abs=0.01)
    assert b["open"] == ["coulex-fit", "rdds", "dsam"]  # a published value does not close a route of one's own
    assert tau["label"] == "τ(2+₁)" and tau["evaluated"] is None and tau["measured"]["count"] == 2
    assert [o["what"] for o in ni["open"]] == ["route", "route", "runs"] and ni["open"][-1]["count"] == 1
    assert nucleus.view("62Ni")["table"] == ni["table"] and nucleus.view("58Ni")["table"] == []

    # One observable of one transition: the comparison, the routes, what the reference gives, the history.
    v = nucleus.observable_view("62Ni", "B(E2)", "0+1 -> 2+1")
    c = v["comparison"]
    assert v["key"] == KEY and c.reference.id == adopted.id
    assert v["reference"]["source"].startswith("ENSDF adopted")
    entries = {v["cells"][e.id]["source"]: e for e in c.entries}
    assert entries["Allmond et al. 2014"].pull == pytest.approx(0.558, abs=0.001)
    assert entries["Allmond et al. 2014"].verdict == "agree"
    theory = entries["shell model, JUN45"]  # no uncertainty: its ratio, and no verdict
    assert (theory.pull, theory.verdict, theory.ratio) == (None, None, pytest.approx(0.4615, abs=1e-4))
    made = entries["Kenn et al. 2001"]  # the DSAM lifetime, as a B(E2): made here, not written to the ledger
    assert v["cells"][made.id]["via"]["id"] == kenn.id and made.method == "converted" and made.verdict == "agree"
    assert made.id not in {r.id for r in results.of_nuclide("62Ni")}
    [p] = c.projected  # the simulated run: closed against what the simulation was given, never pulled
    assert (p.id, p.truth) == (projected.id, 884.1) and p.closure == pytest.approx(0.533, abs=0.001)
    assert (p.reach_percent, p.known_percent, p.best_measured_percent) == (
        pytest.approx(7.04, abs=0.01), pytest.approx(1.48, abs=0.01), pytest.approx(4.08, abs=0.01))
    assert projected.id not in {e.id for e in c.entries}
    routes = {(r["method"], r["kind"]): r for r in v["routes"]}
    assert [(r["status"], r["own"], r["open"]) for r in (
        routes["coulex-yield", "coulex"], routes["coulex-fit", "coulex"], routes["rdds", "lifetimes"],
        routes["dsam", "lifetimes"])] == [("projected", 1, False), ("done", 0, True), ("in physim", 0, True),
                                          ("done", 0, True)]
    # The yield route has a simulated run and no measurement: "projected", not "done". The groups stay apart.
    assert [(routes[r, "coulex"]["independent"], routes[r, "coulex"]["projected"])
            for r in ("coulex-yield", "coulex-fit")] == [(0, 1), (1, 0)]
    assert routes["inelastic-scattering", "inelastic"]["status"] == "later"
    assert routes["shell-model", None]["status"] == "done"
    assert (routes["nrf", None]["status"], routes["nrf", None]["means"]) == ("no experiment", "literature")
    assert all(r["reason"] for r in v["routes"])  # each is tested, or says what is not checked
    gives = {g["to"]: g for g in v["gives"]}
    assert gives["B(E2)↓ in W.u."]["text"] == "12.13 ± 0.18 W.u.↓"
    assert gives["B(E2)↓"]["text"] == "176.8 ± 2.6 e2fm4↓"
    assert gives["tau"]["value"] == pytest.approx(2.08, abs=0.01) and gives["tau"]["key"] == "62Ni:2+1"
    element = gives["matrix element"]  # a square root: its uncertainty has no function in the code
    assert (element["value"], element["plus"], element["unit"]) == (pytest.approx(29.73, abs=0.01), None, "efm2")
    assert [h["status"] for h in v["history"]] == ["holds"] * 4 and v["history"][0]["id"] == allmond.id
    # Another reference; a projected row cannot be one.
    against = nucleus.observable_view("62Ni", "B(E2)", KEY, reference=allmond.id[:6])["comparison"]
    assert against.reference.id == allmond.id
    assert nucleus.gives("62Ni", kenn)[0]["value"] == pytest.approx(914.6, abs=0.1)
    with pytest.raises(ValueError, match="projected row cannot be the reference"):
        nucleus.observable_view("62Ni", "B(E2)", KEY, reference=projected.id)
    assert nucleus.observable_view("62Ni", "Qs", "2+1")["comparison"] is None
    assert [c["evidence"] for c in nucleus.holding("62Ni", "B(E2)")] == [
        "measured", "measured", "evaluated", "theory", "projected"]


# -- a value counts once, other evaluations, and what the screens draw ---------------------------------------------


LITERATURE = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "literature_ni62_cut.json"
KENN = {"authors": "Kenn et al.", "year": 2001, "where": "DSAM"}


def point(value, unit="e2fm4", **more):
    return nucleus.add_data_point("62Ni", "B(E2)", KEY, value, unit, **{
        "direction": "up", "uncertainty": 0.04 * value, "method": "coulex-fit", "reference": ALLMOND, **more})


def add_ledger(*rows):
    results.add_points(list(rows))


def test_a_value_counts_once_whatever_its_unit_or_detail(lab):
    typed = point(906.0, uncertainty=37.0)
    kenn = nucleus.add_data_point("62Ni", "tau", "2+1", 2.01, "ps", uncertainty=0.12, method="dsam", reference=KENN)
    # The same two values in a file: in e²b² and with a long detail (the DOI says so), and with stat and sys (the
    # first author and the year say so, as one of the two has no DOI).
    found = nucleus.import_points(LITERATURE, "62Ni")
    assert {x["where"]: x["id"] for x in found["already"]} == {"quantities.BE2_up[0]": typed.id,
                                                              "quantities.tau_2plus[0]": kenn.id}
    with pytest.raises(nucleus.PointThere, match="is in the ledger of 62Ni already: 906 ± 37 e2fm4↑.*counts once") as e:
        point(12.43, "W.u.", direction="down", uncertainty=0.5, detail="GOSIA")  # 906 e²fm⁴ to three digits
    assert e.value.row == typed
    with pytest.raises(nucleus.PointThere, match="counts once"):
        point(0.0906, "e2b2", uncertainty=0.0037,
              reference={"authors": "Allmond, J. M.", "year": 2014, "doi": "https://doi.org/10.1103/PHYSREVC.90.034309"})
    # Another method, another paper, another year or another value is another point; so is a correction of the same.
    for more in ({"method": "coulex-yield"}, {"reference": {**ALLMOND, "doi": "10.1000/other"}},
                 {"reference": {"authors": "Allmond et al.", "year": 2015}}):
        assert point(906.0, **more).value == 906.0
    assert point(918.0).value == 918.0
    assert nucleus.correct_data_point("62Ni", typed.id, uncertainty=40.0).total == 40.0
    a, b = (results.Row("62Ni", "B(E2)", KEY, v, "e2fm4", "coulex-yield", "measured", direction="up", total=30.0,
                        origin=results.reference(who, 1971, where=where))
            for v, who, where in ((880.0, "NSR 1971ChZF", ""), (880.0, "Roe", "Ann. Phys. 1, 1 (1971) [1971ChZF]")))
    assert nucleus.repeats(a, b)  # by the NSR key number


def test_the_duplicates_of_an_older_ledger_are_withdrawn_once(lab):
    def row(value, unit, created, authors="Allmond et al.", **more):
        return results.Row("62Ni", "B(E2)", KEY, value, unit, "coulex-fit", "measured", direction="up",
                           origin=results.reference(authors, 2014, "10.1103/PhysRevC.90.034309", "Table II"),
                           created=created, **more)

    first = row(906.0, "e2fm4", "2026-11-24T10:00:00+00:00", total=37.0)
    twice = row(0.0906, "e2b2", "2026-11-24T11:00:00+00:00", "J. M. Allmond et al.", total=0.0037)
    fixed = row(906.0, "e2fm4", "2026-11-24T12:00:00+00:00", total=36.0, supersedes=first.id)  # came in first, still
    other = row(860.0, "e2fm4", "2026-11-24T13:00:00+00:00", total=37.0)
    add_ledger(*evaluated(), first, twice, fixed, other)
    file = lab / "nuclei" / "62Ni" / results.POINTS_FILE
    before = file.read_text(encoding="utf-8")
    [gone] = nucleus.retract_duplicates()
    assert gone.supersedes == twice.id and gone.retracted.startswith(f"repeats row {fixed.id} (906 ± 36 e2fm4↑")
    assert file.read_text(encoding="utf-8").startswith(before)  # nothing is deleted
    after = file.read_text(encoding="utf-8")
    assert nucleus.retract_duplicates("62Ni") == [] and nucleus.retract_duplicates() == []
    assert file.read_text(encoding="utf-8") == after
    view = nucleus.observable_view("62Ni", "B(E2)", KEY)
    assert view["counts"]["measured"]["count"] == 2 and "1 evaluated, 2 measured" in logbook.nucleus_html("62Ni")


def test_another_evaluation_is_compared_without_a_verdict_and_is_never_the_reference(lab):
    adopted = evaluated()
    add_ledger(*adopted)
    ensdf_row = adopted[-1]
    allmond = point(906.0, uncertainty=37.0)
    raman = next(x["row"] for x in nucleus.import_points(LITERATURE, "62Ni")["added"]
                 if x["where"] == "quantities.BE2_up[2]")  # written after ENSDF's
    assert nucleus.is_data_point(raman) and not nucleus.is_evaluation(raman)
    view = nucleus.observable_view("62Ni", "B(E2)", KEY)
    entries = {e.id: e for e in view["comparison"].entries}
    assert view["comparison"].reference.id == ensdf_row.id and entries[allmond.id].verdict == "agree"
    e = entries[raman.id]
    assert e.difference_percent == pytest.approx(0.67, abs=0.02) and (e.pull, e.verdict) == (None, None)
    assert e.independent_of_reference is False and entries[allmond.id].independent_of_reference is True
    [group] = [g for g in view["groups"] if g["id"] == "evaluation"]
    assert [i["id"] for i in group["entries"]] == [raman.id]
    assert group["entries"][0]["verdict_text"] == "not judged: not independent"
    counts = view["counts"]
    assert counts["evaluations"] == 1 and raman.id not in {
        i["id"] for g in view["groups"] if g["label"].startswith("Measured") for i in g["entries"]}
    row = nucleus.view("62Ni")["table"][0]
    assert row["evaluated"]["id"] == ensdf_row.id and row["other_evaluations"]["count"] == 1
    # Picked by hand it is the reference, the measurements are judged against it (not independent), ENSDF is not.
    against = nucleus.observable_view("62Ni", "B(E2)", KEY, reference=raman.id)["comparison"]
    entries = {e.id: e for e in against.entries}
    assert entries[allmond.id].verdict == "agree" and entries[allmond.id].independent_of_reference is False
    assert entries[ensdf_row.id].verdict is None and entries[ensdf_row.id].difference is not None
    assert nucleus.correct_data_point("62Ni", raman.id, uncertainty=0.0026).origin == raman.origin


def test_an_observable_ready_to_draw(lab):
    add_ledger(*evaluated(), results.Row(
        "62Ni", "B(E2)", KEY, 918.6, "e2fm4", "coulex-yield", "projected", direction="up", stat=12.9, sys=63.4,
        truth=884.105, origin=results.from_run("62Ni-on-194Pt", "03-beam-24h")))
    allmond = point(906.0, uncertainty=37.0)
    li = point(618.0, uncertainty=42.0, method="electron-scattering", detail="(e,e') at 150 MeV",
               reference={"authors": "NSR 1972Li28", "year": 1972})
    kenn = nucleus.add_data_point("62Ni", "tau", "2+1", 2.01, "ps", uncertainty=0.12, method="dsam", reference=KENN)
    model = point(15.5, "W.u.", direction="down", uncertainty=None, method="shell-model", evidence="theory",
                  reference={"authors": "A. Chakraborty et al.", "year": 2011})
    v = nucleus.observable_view("62Ni", "B(E2)", KEY)
    drawn = {k: v[k] for k in ("groups", "counts", "stands_out", "axis", "units", "orders")}
    json.dumps(drawn)
    groups = {g["id"]: g for g in v["groups"]}
    assert list(groups) == ["reference", "coulex", "electron-scattering", "lifetime", "theory", "projected"]
    [made] = groups["lifetime"]["entries"]
    assert (made["family"], made["year"], made["method"], made["made_from"]) == (
        "lifetime", 2001, "DSAM (line shape)", "from τ = 2.01 ± 0.12 ps") and made["id"] != kenn.id
    assert groups["lifetime"]["note"] == "1 value, converted from the τ it reports"
    assert groups["lifetime"]["label"] == "Measured · from the lifetime τ(2⁺₁)"
    assert "given to its simulation (884.1)" in groups["projected"]["note"]
    [run] = groups["projected"]["entries"]
    assert (run["pull_text"], run["verdict_text"], run["known_text"], run["stat_sys_text"]) == (
        "+0.53σ", "closes", "reaches ±7.0 %", "± 13 (stat) ± 63 (sys)")
    counts = v["counts"]
    assert counts["measured"]["text"] == "2 agree · 1 differs" and counts["measured"]["count"] == 3
    assert counts["theory"]["text"] == "× 1.28 · above" and counts["best_measured"]["id"] == allmond.id
    assert counts["projected"]["closing_text"] == "closes" and "of which 6.9 % is systematic" in counts["projected"]["text"]
    differs, theory = v["stands_out"]
    assert (differs["kind"], differs["ids"], differs["note"]) == ("differ", [li.id], "(e,e') at 150 MeV")
    assert "618 ± 42 e²fm⁴↑, a pull of −6.05: it differs" in differs["text"]
    assert (theory["kind"], theory["ids"]) == ("theory", [model.id]) and "× 1.28" in theory["text"]
    # Everything is on one axis, in each unit: the places are fractions of it.
    axis = v["axis"]
    items = [i for g in v["groups"] for i in g["entries"]]
    assert all(0 <= i["x_low"] <= i["x"] <= i["x_high"] <= 1 for i in items)
    assert axis["low"] < 576 and axis["high"] > 1129 and axis["band"]["x_low"] < axis["band"]["x"]
    place = {i["id"]: (i["value"] - axis["low"]) / (axis["high"] - axis["low"]) for i in items}
    assert all(i["x"] == pytest.approx(place[i["id"]]) for i in items)
    assert [(u["unit"], u["label"]) for u in v["units"]] == [("e2fm4", "e²fm⁴↑"), ("e2b2", "e²b²↑"), ("W.u.", "W.u.↓")]
    # A value without an uncertainty keeps the three digits it was given with, in every unit.
    bare = next(i for i in items if i["id"] == model.id)
    assert (bare["text"], bare["texts"]) == ("1130", {"e2fm4": "1130", "e2b2": "0.113", "W.u.": "15.5"})
    assert v["orders"]["pull"][1] == li.id and v["orders"]["year"][1] == allmond.id
    assert all(sorted(o) == sorted(i["id"] for i in items) for o in v["orders"].values())
    # The card of the landing page says the same; the evaluation's own rows have none.
    table = {(t["observable"], t["key"]): t for t in nucleus.view("62Ni")["table"]}
    card_ = table["B(E2)", KEY]
    assert card_["counts"] == counts and card_["strip"]["axis"] == axis and card_["strip"]["units"] == v["units"]
    assert [d["id"] for d in card_["strip"]["dots"]] == v["orders"]["method"]
    assert next(d for d in card_["strip"]["dots"] if d["id"] == model.id)["texts"] == bare["texts"]
    assert table["E_level", "62Ni:2+1"]["counts"] is None and table["E_level", "62Ni:2+1"]["strip"] is None
    # The τ of the same nucleus: what Coulex gives is converted, what a lifetime gives is not.
    tau = {g["id"]: g for g in nucleus.observable_view("62Ni", "tau", "2+1")["groups"]}
    assert tau["theory"]["entries"][0]["texts"] == {"ps": "1.63", "fs": "1630"}  # by way of the B(E2) it was given as
    assert tau["coulex"]["note"] == "1 value, converted from the B(E2) it reports" and tau["lifetime"]["note"] == "1 value"
    empty = nucleus.observable_view("62Ni", "mu", "2+1")
    assert (empty["groups"], empty["axis"], empty["counts"]["measured"]["count"]) == ([], None, 0)


def test_a_value_made_from_one_that_was_replaced_stands_out(lab):
    run = results.from_run("62Ni-plunger", "01-beam-24h")
    tau = results.Row("62Ni", "tau", "62Ni:2+1", 2.0, "ps", "rdds", "measured", total=0.1, origin=run)
    b = results.Row("62Ni", "B(E2)", KEY, 180.0, "e2fm4", "rdds", "measured", direction="down", total=9.0,
                    depends_on=(tau.id,), origin=run)
    add_ledger(*evaluated(), tau, b, tau.revised(value=2.2))
    found = nucleus.observable_view("62Ni", "B(E2)", KEY)["stands_out"]
    [gone] = [s for s in found if s["kind"] == "superseded"]
    assert gone["ids"] == [b.id] and "replaced by 2.20 ± 0.10 ps" in gone["text"]


def test_the_lowest_levels_come_from_the_evaluation(lab):
    nothing = nucleus.lowest_levels("62Ni")
    assert nothing["levels"] == [] and "no local copy of ENSDF" in nothing["reason"]

    def energy(key, value):
        return results.Row("62Ni", "E_level", key, value, "keV", "evaluation", "evaluated", total=0.03,
                           origin=ENSDF_ORIGIN)

    add_ledger(*evaluated(), energy("62Ni:4+1", 2336.46), energy("62Ni:0+2", 2048.55), energy("62Ni:2+2", 2301.85))
    point(906.0, uncertainty=37.0)
    found = nucleus.lowest_levels("62Ni", count=4)
    json.dumps(found)
    assert [(x["key"], x["jpi"], x["short"]) for x in found["levels"]] == [
        ("62Ni:0+1", "0⁺₁", "0"), ("62Ni:2+1", "2⁺₁", "1172.9 keV"), ("62Ni:0+2", "0⁺₂", "2048.6 keV"),
        ("62Ni:2+2", "2⁺₂", "2301.8 keV")]
    ground, two = found["levels"][:2]
    assert two["text"] == "1172.910 ± 0.030 keV" and ground["ground"] and found["reason"] == ""
    assert [(o["observable"], o["role"], o["headline"]) for o in two["observables"]] == [
        ("B(E2)", "from", True), ("tau", "level", True), ("Qs", "level", False)]
    [b] = ground["observables"]
    assert (b["key"], b["role"], b["other"]) == (KEY, "to", {"key": "62Ni:2+1", "jpi": "2⁺₁"})
    assert b["value"]["text"] == "884 ± 13 e2fm4↑"
    assert found["ratio"]["value"] == pytest.approx(1.992, abs=0.001) and len(nucleus.lowest_levels("62Ni")["levels"]) == 5
