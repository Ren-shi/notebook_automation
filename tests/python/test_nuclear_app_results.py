"""The results and their comparison on an experiment's own screens (results wave 2b): the nucleus in the rail, the
short comparison under a result, the lifetime analysis kept from an earlier session, and where the steps of an
experiment without an analysis lead. Synthetic rows in a temporary folder and no simulated run: seconds."""

import asyncio
import copy
from types import SimpleNamespace

import pytest

from physim.nuclear import ensdf, record, results
from physim.nuclear import workbench as wb
from physim.nuclear.app import tab_analysis as ta
from physim.nuclear.app import tab_lifetimes as tl
from physim.nuclear.experiment import Experiment
from physim.nuclear.planner import Planner
from physim.nuclear.record import _g, provenance, record_html

KEY = "62Ni:2+1>0+1:E2"


@pytest.fixture(autouse=True)
def apart(tmp_path, monkeypatch):
    """The experiments and the nuclei of this test, and no local ENSDF copy (as in CI)."""
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path / "experiments"))
    monkeypatch.setenv("PHYSIM_NUCLEI", str(tmp_path / "nuclei"))
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))


def _row(value, method, evidence, origin, **more):
    return results.Row("62Ni", "B(E2)", KEY, value, "e2fm4", method, evidence, origin=origin, direction="up", **more)


def _experiment(tmp_path):
    """An experiment's folder with two projected rows of one run: the analysis, and a later fit to all orders."""
    folder = tmp_path / "experiments" / "62Ni-on-194Pt"
    folder.mkdir(parents=True)
    (folder / "experiment.toml").write_text(Experiment.example("coulex_ni62_pt194").to_toml(), encoding="utf-8")
    run = results.from_run(folder.name, "03-beam-24h")
    results.append(folder, _row(918.6, "coulex-yield", "projected", run, stat=13.0, sys=63.7, truth=884.1))
    results.append(folder, _row(901.0, "coulex-fit", "projected", run, stat=20.0, truth=884.1))
    return SimpleNamespace(experiment_folder=folder, experiment=Experiment.example("coulex_ni62_pt194"))


def _lines(p, value=918.6, missing=""):
    view = ta.comparison_of(p, "B(E2)", KEY)
    mine = ta.own_projection(view, p.experiment_folder.name, value, run="03-beam-24h")
    short = ta.short_comparison(view, mine, missing)
    return short, dict(short["lines"])


def test_the_short_comparison_with_a_reference(tmp_path):
    p = _experiment(tmp_path)
    results.add_point(_row(884.0, "evaluation", "evaluated", results.reference("ENSDF", 2012), total=13.0))
    results.add_point(_row(890.0, "coulex-fit", "measured", results.reference("A. Author", 2014), total=30.0))
    short, said = _lines(p)
    assert short["title"] == "B(E2; 0⁺₁ → 2⁺₁)↑ of ⁶²Ni"
    assert "simulated and analysed back" in said["Projected"] and "not a measurement of ⁶²Ni" in said["Projected"]
    assert said["Reference"].startswith("884 ± 13 e²fm⁴ ↑, evaluated (ENSDF")
    # Closed against what the simulation was given, and set beside how well the value is known: no pull.
    assert "gives back the value the simulation was given (884.1 e²fm⁴) within 0.5σ" in said["Closure"]
    assert "This setup reaches ±7.1%; the value is known to ±1.5%" in said["Reach"].replace(" %", "%")
    assert "the best single measurement reaches ±3.4%" in said["Reach"].replace(" %", "%")
    assert "pull" not in " ".join(said.values())
    assert said["Other routes"].startswith("Of the other routes, ") and "hold a value" in said["Other routes"]
    assert "still open" in said["Other routes"] and "Evaluation" not in said["Other routes"].split("still open")[1]
    assert short["link"] == {"text": "Compare on the ⁶²Ni page",
                             "href": "/?nucleus=62Ni&observable=B%28E2%29&key=62Ni%3A2%2B1%3E0%2B1%3AE2"}
    # The value shown is not the ledger's: the row is of another analysis, and is not shown as this one's.
    short, said = _lines(p, 950.0, "Not in the results ledger yet.")
    assert said["This result"] == "Not in the results ledger yet." and "Closure" not in said and "Reach" not in said


def test_the_short_comparison_without_a_reference(tmp_path):
    short, said = _lines(_experiment(tmp_path))
    assert said["Reference"] == "None: no evaluated value of B(E2; 0⁺₁ → 2⁺₁)↑ is entered for ⁶²Ni."
    assert "within 0.5σ" in said["Closure"]
    assert "no evaluated value is entered" in said["Reach"] and "no measured point" in said["Reach"]
    assert "none holds a value yet" in said["Other routes"]
    # The same run's fit to all orders is another route: simulated only, and not this result.
    assert "1 has simulated runs only" in said["Other routes"]


def test_the_nucleus_in_the_rail():
    assert wb.studied(Experiment.example("coulex_ni62_pt194")) == {
        "nuclide": "62Ni", "label": "⁶²Ni", "href": "/?nucleus=62Ni"}
    fusion = Experiment.example("fusion_dy158").to_dict()
    fusion["run"].pop("channel", None)
    assert wb.studied(Experiment.from_dict(fusion)) is None  # no channel is chosen: no nucleus
    fusion["run"]["channel"] = "158Dy"
    assert wb.studied(Experiment.from_dict(fusion))["href"] == "/?nucleus=158Dy"
    assert wb.studied(Experiment.example("alpha_on_gold")) is None


@pytest.mark.parametrize("example, kind, needs", [
    ("fusion_dy158", "fusion", "cross sections through the efficiency"),
    ("alpha_on_gold", "beam", "Rutherford cross section")])
def test_the_steps_without_an_analysis_lead_to_the_report(example, kind, needs):
    p = Planner.example(example)
    assert {s["key"]: s for s in wb.steps(p)}["analysis"]["state"] == "waiting"
    p.run = SimpleNamespace(kind=kind, number=1, experiment=p.experiment)
    data, analysis = wb.next_step(p, "data"), wb.next_step(p, "analysis")
    assert (data["label"], data["to"]) == ("Next: the report", "report") and needs in data["hint"]
    assert "does not analyse a run" in data["hint"] and "the plan's rates" in analysis["hint"]
    assert analysis["to"] == "report" and wb.next_step(p, "report")["to"] is None
    # A run of a Coulomb excitation is still analysed after the setup became another kind.
    p.run = SimpleNamespace(kind="beam", number=1, experiment=Experiment.example("coulex_ni62_pt194"))
    assert wb.next_step(p, "data")["to"] == "analysis"


ANALYSIS = {
    "available": True, "reason": None, "level": "2+", "method": "ddcm", "tau_ps": 7.6, "stat_ps": 0.31,
    "budget": {"statistics": 0.31, "stopping": 0.12}, "points": [],
    "curve": {"x": [1.0, 10.0, 100.0], "y": [0.97, 0.75, 0.32], "kind": "decay_curve"},
    "b_down": [610.0, 30.0], "b_unit": "e2fm4", "b_multipole": "E2",
    "compare": [{"source": "level scheme (ENSDF)", "tau_ps": 8.4, "b_down": 600.0},
                {"source": "Coulex B(E2) of the setup", "tau_ps": 8.3, "b_down": 605.0}],
    "transition": {"nuclide": "98Ru", "gamma_kev": 652.4,
                   "initial": {"index": 1, "energy_kev": 652.4, "spin": 2.0, "parity": 1},
                   "final": {"index": 0, "energy_kev": 0.0, "spin": 0.0, "parity": 1}},
    "runs": ["01-beam-1h"], "notes": ["Fixture note."]}


def test_the_kept_lifetime_analysis_is_shown_after_reopening(monkeypatch):
    pytest.importorskip("plotly")
    pytest.importorskip("nicegui")
    from nicegui import ui
    from nicegui.testing import user_simulation

    from physim.nuclear import result_rows
    from physim.nuclear.app.page import build_page

    p = Planner.example("coulex_ru98_plunger")
    folder = p.create_experiment("ru98")
    p._keep("lifetime.json", copy.deepcopy(ANALYSIS), p.experiment_folder)
    result_rows.write_lifetimes(p, ANALYSIS)
    p = Planner.open_experiment(folder)
    status = {s["key"]: s for s in wb.steps(p)}["analysis"]
    assert (status["status"], status["state"]) == ("τ(2+) 7.6 ± 0.31 ps", "done")
    # A run taken since, at a distance, is not in it; a run without a distance is not a plunger's.
    runs = [{"label": "01-beam-1h", "kind": "beam", "distance_um": 3.0}, {"label": "02", "kind": "beam"},
            {"label": "03", "kind": "beam", "distance_um": 9.0}]
    assert wb.runs_not_read(p, ANALYSIS, runs) == ["03"]
    # The choice of level and method is the experiment's: another one begins afresh.
    ctx = SimpleNamespace(state={"planner": p})
    tl.analysis_choice(ctx)["session"] = True
    assert tl.analysis_choice(ctx)["session"]
    ctx.state["planner"] = Planner.example("fusion_cr50_plunger")
    assert not tl.analysis_choice(ctx)["session"]

    def run_again(self, level=None, series=None, method=None):
        raise AssertionError("the kept analysis is shown, not run again")

    monkeypatch.setattr(Planner, "analyse_lifetimes", run_again)

    async def main():
        async with user_simulation(lambda: build_page(prebuild=True, experiment=str(folder))) as user:
            await user.open("/")
            await user.should_see("Kept from an earlier session")
            texts = " | ".join(str(e.text) for e in user.find(kind=ui.label).elements)
            links = [(str(e.text), e.props["href"]) for e in user.find(kind=ui.link).elements]
            for part in ("7.6 ± 0.31 ps", "τ(2⁺ (652 keV)) of ⁹⁸Ru", "B(E2; 0⁺ → 2⁺ (652 keV))↑ of ⁹⁸Ru",
                         "PROJECTED", "the value the simulation was given (8.4 ps)"):
                assert part in texts, part
            # The step's own table keeps only what the comparison above does not say.
            table = [r["s"] for t in user.find(kind=ui.table).elements for r in t.rows if set(r) == {"s", "t", "b"}]
            assert table == ["Coulex B(E2) of the setup"]
            assert "could not be shown" not in texts
            assert ("⁹⁸Ru", "/?nucleus=98Ru") in links
            assert any(v.startswith("/?nucleus=98Ru&observable=tau&key=98Ru") for _, v in links)
            await asyncio.sleep(0.2)

    asyncio.run(main())


def test_the_record_says_what_the_setup_is(monkeypatch):
    assert (_g(130.0, 2), _g(437100.0), _g(None)) == ("130", "437,100", "—")  # it was "1.3e+02"
    ni = Planner.example("coulex_ni62_pt194")
    d = ni.experiment.to_dict()
    d["reaction"]["b_up"] = "0.1 e2b2"
    own = Experiment.from_dict(d)
    assert "Reaction: Fusion-evaporation" in record_html(Planner.example("fusion_dy158"))
    assert "Reaction: Elastic scattering" in record_html(Planner.example("alpha_on_gold"))
    # The B(E2) that the level scheme's matrix element gives is sourced as that element is; one's own is "user".
    source = [next(r for r in provenance(x) if r[0].startswith("State at"))[2] for x in (ni.experiment, own)]
    assert source[0].startswith("derived (from B(E2) = 12.13 W.u. in ENSDF") and source[1].startswith("user (")
    # The ratio of B↑ and B↓ is the spins': 5 for 0⁺ → 2⁺, and nothing without the spins.
    ratio, upper, lower = record._decay(ni.experiment)
    assert (round(ratio, 9), upper, lower) == (5.0, "2⁺ (1173 keV)", "0⁺")  # (no ENSDF: no ordinal)
    unknown = results.TransitionKey.parse("62Ni:?@1173>0+1:E2")
    monkeypatch.setattr(results, "key_of_reaction", lambda experiment: unknown)
    assert record._decay(ni.experiment) is None
