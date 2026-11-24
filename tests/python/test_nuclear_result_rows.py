"""The planner's ledger rows (wave 1 of the results pages): what each analysis returns becomes rows that are
appended, a second analysis supersedes the first, the shell-model copy adds a calculation once, and a ledger that
cannot be written never fails the analysis. Built on made-up results: no run is simulated here."""

import json
import math
from types import SimpleNamespace

import pytest

from physim.nuclear import Experiment, ensdf, result_rows, results
from physim.nuclear.analysis import PeakFit, Result, save_result
from physim.nuclear.levels import Level, LevelScheme, Value, quadrupole_factor
from physim.nuclear.planner import Planner
from physim.nuclear.runs import ExperimentFolder
from physim.nuclear.shellmodel import ShellModelResult, ShellState, ShellTransition

E2 = "62Ni:2+@1173>0+@0:E2"  # unanchored: these tests have no ENSDF copy


@pytest.fixture(autouse=True)
def apart(tmp_path, monkeypatch):
    """No ENSDF copy (the keys do not depend on the computer), and the nuclei home in the test's own folder."""
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))


def analysis_result(b=918.6, **more) -> Result:
    """A Coulex analysis as the Analysis step returns it, with made-up numbers."""
    fit = PeakFit(area=5200.0, area_unc=75.0, centroid=1.1729, sigma=0.0012, background=40.0, chi2=1.1, bins=24)
    d = dict(settings={"particle_gate": "inelastic", "particle_energy": (20.0, 30.0)}, steps=[], fit=fit,
             area=5200.0, efficiency=0.021, correlation=1.02, branch=1.0, normalisation="rutherford",
             norm_counts=1.0e7, b_e2fm4=b, statistical=0.014, systematic=0.069, monte_carlo=0.004, b_e2b2=b * 1e-4, b_wu=12.4,
             budget={"efficiency": 0.05, "angular correlation": 0.03}, truth_e2fm4=884.105, pull=0.5, shape={},
             counts_per_shift=900.0, hours_for_precision=30.0, gate="inelastic")
    return Result(**{**d, **more})


@pytest.fixture
def lab(tmp_path):
    """An experiment folder with one beam run, and a stand-in for the planner that is open on it."""
    folder = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=tmp_path / "experiments")
    run = folder.path / "runs" / "03-beam-24h"
    run.mkdir(parents=True)
    folder.setup().save(run / "setup.toml")
    (run / "summary.json").write_text(json.dumps({"number": 3, "kind": "beam"}), encoding="utf-8")
    current = SimpleNamespace(folder=run, label=run.name)
    planner = SimpleNamespace(experiment_folder=folder.path, _beam_run=lambda: current,
                              _data_experiment=folder.setup)
    return SimpleNamespace(folder=folder.path, run=run, planner=planner)


def lifetime_analysis(tau=2.10, method="ddcm", **more) -> dict:
    """What Planner.analyse_lifetimes returns for the first 2⁺ of ⁶²Ni, with made-up numbers."""
    budget = {"statistics": 0.08, "distance": 0.03, "velocity": 0.04, "feeding": 0.0, "stopping": 0.0}
    out = {"available": True, "level": "2+", "method": method, "tau_ps": tau, "stat_ps": 0.08, "budget": budget,
           "total_ps": math.sqrt(sum(v * v for v in budget.values())), "b_down": [176.0, 8.5], "b_unit": "e2fm4",
           "b_multipole": "E2", "runs": ["04-beam-10um", "05-beam-30um", "06-beam-100um"],
           "compare": [{"source": "level scheme: derived", "tau_ps": 2.09, "b_down": 176.8},
                       {"source": "ENSDF (local copy)", "tau_ps": 2.01, "b_down": 183.9}],
           "transition": {"nuclide": "62Ni", "gamma_kev": 1172.9,
                          "initial": {"index": 1, "energy_kev": 1172.9, "spin": 2.0, "parity": 1},
                          "final": {"index": 0, "energy_kev": 0.0, "spin": 0.0, "parity": 1}}}
    return {**out, **more}


def shell_result(b_up=100.0) -> ShellModelResult:
    states = [ShellState(0, 1, 1, -203.1), ShellState(4, 1, 1, -201.9, ex_mev=1.2, q_s=-12.0, mu=0.4),
              ShellState(8, 1, 1, -200.7, ex_mev=2.4)]  # (2J, parity, index), KSHELL's energy, MeV
    transitions = [ShellTransition("E2", (4, 1, 1), (0, 1, 1), 10.0, b_down=b_up / 5, b_up=b_up),
                   ShellTransition("E2", (4, 1, 1), (4, 1, 1), -15.8, moment=-12.0),
                   ShellTransition("M1", (4, 1, 1), (4, 1, 1), 0.9, moment=0.4)]
    return ShellModelResult(None, "f5pg9", "jun45", "R2", (1.5, 0.5), None, None, states, transitions,
                            nuclide="62Ni")


# -- Coulex ---------------------------------------------------------------------------------------------------------


def test_a_coulex_analysis_gives_the_row_the_backfill_would(lab):
    r = analysis_result()
    save_result(lab.run, r)
    (row,) = result_rows.write_coulex(lab.planner, r)
    assert (row.nuclide, row.observable, row.key) == ("62Ni", "B(E2)", E2) and (row.unit, row.direction) == (
        "e2fm4", "up")
    assert (row.method, row.evidence, row.truth) == ("coulex-yield", "projected", 884.105)  # never measured
    assert row.stat == pytest.approx(918.6 * 0.014) and row.sys == pytest.approx(918.6 * 0.069)
    assert row.budget == {"efficiency": pytest.approx(45.93), "angular correlation": pytest.approx(27.558)}
    assert row.origin["experiment"] == lab.folder.name and row.origin["run"] == "03-beam-24h"
    assert row.origin["settings"]["particle_energy"] == [20.0, 30.0] and row.code == results.code()
    assert results.of_experiment(lab.folder) == [row]
    assert results.backfill(lab.folder) == []  # the kept file's own row is there already

    (lab.folder / results.RESULTS_FILE).unlink()
    (old,) = results.backfill(lab.folder)  # the same analysis, read from the kept file instead
    assert result_rows._says(old) == result_rows._says(row)


def test_a_second_analysis_supersedes_the_first(lab):
    (first,) = result_rows.write_coulex(lab.planner, analysis_result())
    assert result_rows.write_coulex(lab.planner, analysis_result()) == [first]  # nothing new said: no new row
    changed = analysis_result(b=900.0, normalisation="target")
    (second,) = result_rows.write_coulex(lab.planner, changed)
    assert second.supersedes == first.id and second.detail == "target normalisation, inelastic gate"
    rows = results.of_experiment(lab.folder)
    assert rows == [first, second] and results.current(rows) == [second]
    other = SimpleNamespace(folder=lab.run.with_name("07-beam-1h"), label="07-beam-1h")
    lab.planner._beam_run = lambda: other
    (third,) = result_rows.write_coulex(lab.planner, analysis_result())  # another run: beside it, not instead
    assert third.supersedes is None
    assert results.current(results.of_experiment(lab.folder)) == [second, third]


def test_an_old_experiment_is_backfilled_when_it_is_opened(lab):
    save_result(lab.run, analysis_result())
    assert results.of_experiment(lab.folder) == []
    Planner.open_experiment(lab.folder)
    (row,) = results.of_experiment(lab.folder)
    assert row.origin["backfilled"] is True and row.value == 918.6
    assert result_rows.backfill_if_new(lab.folder) == []  # once: the ledger is there now


# -- lifetimes ------------------------------------------------------------------------------------------------------


def test_a_lifetime_analysis_gives_tau_and_the_b_that_depends_on_it(lab):
    a = lifetime_analysis()
    tau, b = result_rows.write_lifetimes(lab.planner, a)
    assert (tau.observable, tau.key, tau.method) == ("tau", "62Ni:2+@1173", "rdds")
    assert (tau.value, tau.unit) == (2.10, "ps")
    assert (tau.evidence, tau.truth, tau.stat) == ("projected", 2.09, 0.08)
    assert tau.unc() == pytest.approx(a["total_ps"])
    assert set(tau.budget) == {"distance", "velocity", "feeding", "stopping"}
    assert tau.origin["runs"] == a["runs"] and tau.origin["settings"] == {"method": "ddcm", "level": "2+"}
    assert (b.observable, b.key, b.value, b.direction, b.total) == ("B(E2)", E2, 176.0, "down", 8.5)
    assert (b.method, b.evidence, b.truth, b.depends_on) == ("rdds", "projected", 176.8, (tau.id,))

    assert result_rows.write_lifetimes(lab.planner, a) == [tau, b]  # analysed again, the same: nothing new
    assert len(results.of_experiment(lab.folder)) == 2
    tau2, b2 = result_rows.write_lifetimes(lab.planner, lifetime_analysis(tau=2.30, b_down=[160.7, 8.0]))
    assert (tau2.supersedes, b2.supersedes, b2.depends_on) == (tau.id, b.id, (tau2.id,))
    assert results.current(results.of_experiment(lab.folder)) == [tau2, b2]

    dsam, _ = result_rows.write_lifetimes(lab.planner, lifetime_analysis(method="lineshape"))
    assert dsam.method == "dsam" and dsam.supersedes is None  # another method: beside the plunger's rows
    # What has no row: an M1 line's B (not in the registry yet), and a level whose lifetime the setup lacks.
    assert result_rows.lifetime_rows(lifetime_analysis(b_multipole="M1"), "x", tau) == []
    assert result_rows.write_lifetimes(lab.planner, lifetime_analysis(compare=[])) == []
    assert result_rows.write_lifetimes(lab.planner, {"available": False, "reason": "No plunger runs"}) == []


def test_the_last_lifetime_analysis_is_kept_with_the_experiment(lab):
    p = Planner(ExperimentFolder(lab.folder).setup())
    p.folder = ExperimentFolder(lab.folder)
    assert p.kept_lifetimes() is None
    p._keep("lifetime.json", lifetime_analysis(), p.experiment_folder)
    assert Planner.open_experiment(lab.folder).kept_lifetimes()["tau_ps"] == 2.10


# -- the fit to all orders ------------------------------------------------------------------------------------------


def test_a_fit_gives_b_and_qs_and_leaves_out_what_the_registry_lacks():
    scheme = LevelScheme("62Ni", [Level(Value(0.0, source="user"), 0.0, +1, "0+"),
                                  Level(Value(1172.9, source="user"), 2.0, +1, "2+")])
    scheme.set_matrix_element(0, 1, "E2", 29.7)
    scheme.set_matrix_element(1, 1, "E2", 12.0)
    fit = {"values": {(0, 1, "E2"): 30.0, (1, 1, "E2"): 10.0, (0, 1, "E3"): 55.0},
           "uncertainties": {(0, 1, "E2"): 0.5, (1, 1, "E2"): 4.0, (0, 1, "E3"): 9.0}, "chi2": 1.2,
           "degrees_of_freedom": 3}
    b, q = result_rows.fit_rows(fit, scheme, "ni-on-pt", "03-beam-24h", {"role": "target"})
    assert (b.observable, b.key, b.direction) == ("B(E2)", E2, "up")
    assert (b.method, b.evidence) == ("coulex-fit", "projected")
    assert (b.value, b.truth) == (pytest.approx(900.0), pytest.approx(29.7 ** 2))
    assert b.stat == pytest.approx(900.0 * 2 * 0.5 / 30.0) and "= 30 ± 0.5 e fm²" in b.detail
    assert b.origin["matrix_element"] == {"from": 0, "to": 1, "value": 30.0, "uncertainty": 0.5}
    assert (q.observable, q.key, q.unit) == ("Qs", "62Ni:2+@1173", "efm2")
    factor = quadrupole_factor(2.0)
    assert (q.value, q.stat, q.truth) == tuple(pytest.approx(x * factor) for x in (10.0, 4.0, 12.0))


# -- the shell model ------------------------------------------------------------------------------------------------


def test_a_shell_model_calculation_is_copied_once(tmp_path):
    added = result_rows.copy_shell_model(shell_result())
    assert sorted((r.observable, r.key, r.value) for r in added) == [
        ("B(E2)", "62Ni:2+(1)>0+(1):E2", 100.0), ("E_level", "62Ni:2+(1)", 1200.0),
        ("E_level", "62Ni:4+(1)", 2400.0), ("Qs", "62Ni:2+(1)", -12.0)]
    assert {(r.method, r.evidence, r.truth, r.unc()) for r in added} == {("shell-model", "theory", None, None)}
    assert added[0].origin == results.shell_model_job("jun45", "f5pg9", {"proton": 1.5, "neutron": 0.5},
                                                      program="KSHELL R2")
    assert results.of_nuclide("62Ni", root=tmp_path / "none") == sorted(added, key=lambda r: r.created)

    assert result_rows.copy_shell_model(shell_result()) == []  # the same calculation, whichever way it comes
    assert result_rows.copy_shell_model(shell_result(), job="0123456789abcdef") == []
    assert result_rows.read_shell_model(shell_result(), "0123456789abcdef") == []
    assert len(results.read(tmp_path / "nuclei" / "62Ni" / results.POINTS_FILE)) == 4

    (again,) = result_rows.copy_shell_model(shell_result(b_up=104.0))  # run again, and one value moved
    old = next(r for r in added if r.observable == "B(E2)")
    assert (again.value, again.supersedes) == (104.0, old.id)


# -- a ledger that cannot be written --------------------------------------------------------------------------------


def test_a_failed_write_never_fails_the_analysis(lab, tmp_path, monkeypatch, caplog):
    blocked = tmp_path / "blocked"
    blocked.write_text("a file where the folder should be", encoding="utf-8")
    lab.planner.experiment_folder = blocked  # results.jsonl cannot be made under a file
    assert result_rows.write_coulex(lab.planner, analysis_result()) == []
    assert result_rows.write_lifetimes(lab.planner, lifetime_analysis()) == []
    monkeypatch.setenv(results.ENV, str(blocked))
    assert result_rows.copy_shell_model(shell_result()) == []
    assert result_rows.backfill_if_new(blocked / "nowhere") == []
    # A level that cannot be named: no row, and the analysis's own result is as it was.
    lab.planner.experiment_folder = lab.folder
    unnamed = lifetime_analysis()
    unnamed["transition"] = {**unnamed["transition"], "nuclide": "not a nuclide"}
    before = json.dumps(unnamed, sort_keys=True)
    assert result_rows.write_lifetimes(lab.planner, unnamed) == []
    assert json.dumps(unnamed, sort_keys=True) == before and results.of_experiment(lab.folder) == []
    assert "no ledger row for the lifetime analysis" in caplog.text
