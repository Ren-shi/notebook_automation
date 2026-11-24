"""The nucleus above the experiments (wave 1 of the results pages): the evaluation as rows, data points that are
corrected and withdrawn but never deleted, and the nucleus views on a temporary folder with the ⁶²Ni numbers of
the mock-ups. Plumbing on synthetic rows: no run is simulated here."""

import json

import pytest

from physim.nuclear import Experiment, compare, ensdf, lifetimes, nucleus, results
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
    assert len(by) == len(rows) and {r.observable for r in rows} == {"E_level", "tau", "Qs", "B(E2)"}
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
                               evidence="evaluated")
    said = "; ".join(refused.value.problems)
    for part in ("a transition is two levels", "the value must be a number", "'barn' is not a unit of B(E2)",
                 "needs its direction", "one form", "'uncertainty' must be a number that is not negative",
                 "not 'evaluated': the evaluation is taken from ENSDF", "the method is one of",
                 "where the point is from"):
        assert part in said, part
    assert len(refused.value.problems) == 9 and not (lab / "nuclei").exists()
    good = {"direction": "up", "uncertainty": 37.0, "method": "coulex-fit", "reference": ALLMOND}
    for change, message in (
            ({"uncertainty": None}, "a measured point needs its uncertainty"),
            ({"method": "shell-model"}, "not a method of a measured point"),
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
    assert routes["other", "inelastic"]["status"] == "later" and routes["shell-model", None]["status"] == "done"
    assert (routes["nrf", None]["status"], routes["nrf", None]["means"]) == ("no experiment", "literature")
    assert all((r["applies"], r["reason"]) == (None, "") for r in v["routes"])  # the place a later wave fills
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
