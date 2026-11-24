"""The results tools of the MCP server (physim mcp): the nuclei, the rows of one nucleus, a data point added,
corrected and withdrawn, and the comparison, on a temporary folder with synthetic rows. No experiment is opened."""

import asyncio
import json

import pytest

pytest.importorskip("mcp.server.mcpserver", reason="the MCP extra is not installed")

from mcp.client import Client  # noqa: E402

from physim.mcp.server import INSTRUCTIONS, build_server  # noqa: E402
from physim.mcp.session import Session  # noqa: E402
from physim.mcp.tools_results import ResultTools  # noqa: E402
from physim.nuclear import Experiment, ensdf, nucleus, results  # noqa: E402
from physim.nuclear.runs import ExperimentFolder  # noqa: E402

KEY = "62Ni:2+1>0+1:E2"
NI62 = ((0.0, 0, 1, 1), (1172.91, 4, 1, 1), (2048.55, 0, 1, 2))
ALLMOND = {"authors": "Allmond et al.", "year": 2014, "doi": "10.1103/PhysRevC.90.034309", "where": "Table II"}


@pytest.fixture
def tools(tmp_path, monkeypatch):
    """The tools on an experiments folder with one Coulex experiment and its simulated run's row; the nuclei's own
    rows beside it; no local ENSDF copy, and the ⁶²Ni levels to anchor to."""
    monkeypatch.setenv(results.ENV, str(tmp_path / "elsewhere"))  # the session's own folder wins over this
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setattr(results, "ensdf_levels", lambda nuclide: NI62 if nuclide == "62Ni" else ())
    monkeypatch.setattr(nucleus, "decay", lambda nuclide, level_key: (None, None, "no γ ray known"))
    root = tmp_path / "experiments"
    coulex = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=root)
    (coulex.path / "runs" / "03-beam-24h").mkdir()
    (coulex.path / "runs" / "03-beam-24h" / "summary.json").write_text("{}", encoding="utf-8")
    results.append(coulex.path, results.Row(
        "62Ni", "B(E2)", KEY, 918.6, "e2fm4", "coulex-yield", "projected", direction="up", stat=12.9, sys=63.4,
        truth=884.1, origin=results.from_run(coulex.path.name, "03-beam-24h")))
    made = ResultTools(Session(root))
    assert made.home == tmp_path / "nuclei" and ResultTools(Session()).home == tmp_path / "elsewhere"
    results.add_point(results.Row(
        "62Ni", "B(E2)", KEY, 12.13, "W.u.", "evaluation", "evaluated", direction="down", total=0.18,
        origin=results.reference("ENSDF, adopted levels and γ rays", 2025, where="NDS 204, 1 (2025)",
                                 evaluation=nucleus.ENSDF, dataset_date="202507")), made.home)
    return made


def _data(result):
    assert not result.is_error, result.content[0].text
    if result.structured_content is not None:
        return result.structured_content
    return json.loads(result.content[0].text)


def test_a_data_point_through_the_tools(tools):
    added = tools.add_data_point("Ni-62", "B(E2)", "2+1 > 0+1", 906.0, "coulex-fit", direction="up",
                                 uncertainty=37.0, reference=ALLMOND, detail="GOSIA")
    assert added["added"] == "B(E2; 0+₁ → 2+₁)↑ = 906 ± 37 e2fm4↑" and added["key"] == KEY
    assert added["source"] == "Allmond et al. 2014" and "supersedes" not in added
    with pytest.raises(nucleus.PointProblem, match="needs its direction.*a measured point needs its uncertainty"):
        tools.add_data_point("62Ni", "B(E2)", KEY, 906.0, "coulex-fit", reference=ALLMOND)
    fixed = tools.add_data_point("62Ni", "B(E2)", KEY, 896.0, "coulex-fit", direction="up", uncertainty=37.0,
                                 reference=ALLMOND, corrects=added["id"])
    assert fixed["supersedes"] == added["id"]

    rows = tools.list_results("62Ni", "B(E2)", "2+1 > 0+1")["rows"]
    json.dumps(rows)
    assert [(r["evidence"], r["value"]) for r in rows] == [
        ("measured", "896 ± 37 e2fm4↑"), ("evaluated", "884 ± 13 e2fm4↑"), ("projected", "919 ± 65 e2fm4↑")]
    assert rows[1]["as_given"].startswith("12.13 W.u.↓") and rows[2]["simulation_was_given"] == "884.1 e2fm4"
    with pytest.raises(ValueError, match="Name the observable"):
        tools.list_results("62Ni", key="2+1")

    gone = tools.retract_data_point("62Ni", fixed["id"][:8], "the paper was withdrawn")
    assert gone["retracted"] == fixed["id"] and gone["reason"] == "the paper was withdrawn"
    assert [r["evidence"] for r in tools.list_results("62Ni", "B(E2)")["rows"]] == ["evaluated", "projected"]
    past = tools.list_results("62Ni", "B(E2)", KEY, history=True)["history"]
    assert [h["status"] for h in past] == ["retraction", "retracted", "superseded", "holds", "holds"]
    assert past[1]["reason"] == "the paper was withdrawn" and past[2]["superseded_by"] == fixed["id"]


def test_the_nuclei_and_the_comparison_through_the_tools(tools):
    tools.add_data_point("62Ni", "B(E2)", KEY, 906.0, "coulex-fit", direction="up", uncertainty=37.0,
                         reference=ALLMOND)
    tools.add_data_point("62Ni", "B(E2)", KEY, 408.0, "shell-model", evidence="theory", direction="up",
                         own_measurement={"date": "2026-11-22", "note": "JUN45, KSHELL"})
    listed = tools.list_nuclei()
    json.dumps(listed)
    ni, pt = listed["nuclei"]
    assert (ni["nucleus"], ni["rows_that_hold"], ni["data_points"]) == ("62Ni", 4, 2)
    assert ni["experiments"][0].endswith("studied: 1 run, 1 analysed; latest B(E2; 0+₁ → 2+₁)↑ = 919 ± 65 e2fm4↑ "
                                         "(projected)")
    assert ni["observables"] == [
        f"B(E2; 0+₁ → 2+₁)↑ [{KEY}]: evaluated 884 ± 13 e2fm4↑; 1 measured: 906 ± 37 e2fm4↑ (Allmond et al. "
        "2014); 1 from theory: 408 e2fm4↑ (own measurement, 2026-11-22); projected by 1 simulated run: "
        "reaches ±7 %"]
    assert ni["open"] == [
        "B(E2; 0+₁ → 2+₁)↑: Coulex fitted to all orders, RDDS (decay curve, DDCM), DSAM (line shape)"]
    assert pt["nucleus"] == "194Pt" and "partner" in pt["experiments"][0] and listed["not_in_a_nucleus_yet"] == []
    assert tools.list_results("62Ni")["observables"] == ni["observables"]

    c = tools.compare_results("62Ni", "B(E2)", "0+1 -> 2+1")
    json.dumps(c)
    assert c["reference"].startswith("ENSDF adopted, NDS 204, 1 (2025): 884 ± 13 e2fm4↑")
    against = {r["source"]: r.get("against_reference") for r in c["independent"]}
    assert against["Allmond et al. 2014"] == "+2.5 %, pull +0.56: agree"
    assert against["own measurement, 2026-11-22"] == "× 0.461 (no uncertainty: not judged)"
    [p] = c["projected"]  # apart from the independent rows, and closed against what the simulation was given
    assert p["closure"].startswith("+0.53 standard deviations") and p["closure"].endswith(": agree")
    assert (p["setup_reaches"], p["known_to"], p["best_single_measurement"]) == ("±7 %", "±1.5 %", "±4.1 %")
    assert "against_reference" not in p and "never pulled against" in c["rule"]
    routes = {r["route"]: r for r in c["routes"]}
    assert routes["Coulex from the yield (first order)"] == {  # a simulated run, nothing measured
        "route": "Coulex from the yield (first order)", "experiment": "Coulomb excitation", "status": "projected",
        "independent": 0, "projected": 1, "own": 1, "note": "the γ-ray yield per scattered particle, first order",
        "applies": "yes", "reason": "62Ni (3.63 % of natural Ni) can be the target or the beam"}
    fit = routes["Coulex fitted to all orders"]  # a published value: done, and still open to plan oneself
    assert [fit[k] for k in ("status", "independent", "projected", "own", "open")] == ["done", 1, 0, 0, True]
    assert routes["RDDS (decay curve, DDCM)"] == {
        "route": "RDDS (decay curve, DDCM)", "experiment": "Lifetimes (RDDS, DSAM)", "status": "in physim",
        "open": True, "note": "through τ, the γ-ray energy, the branching and the conversion",
        "applies": "not tested", "reason": "no τ of the level 2+₁ is entered (take the evaluation or add a data "
                                           "point): whether it is in the range of a plunger or of DSAM is not tested"}
    gives = {g["gives"]: g for g in c["reference_gives"]}
    assert gives["B(E2)↓ in W.u."]["value"] == "12.13 ± 0.18 W.u.↓" and gives["tau"] == {
        "gives": "tau", "missing": "no γ ray known"}
    rotor, theory = c["checks"]
    assert rotor["rows"] and "No Q_s(2⁺₁) is entered" in rotor["says"] and rotor["how_to_read_it"]
    assert theory["says"].startswith("0.46 = 408 e²fm⁴ ↑") and len(theory["rows"]) == 2
    assert all(u["text"] and u["reference"] for u in c["used_for"])
    empty = tools.compare_results("62Ni", "tau", "2+1")
    assert "note" in empty and "independent" not in empty and empty["routes"]
    assert not (tools.home.parent / "elsewhere").exists()


def test_a_reference_of_zero_is_compared_in_the_unit(tools):
    results.add_point(results.Row("62Ni", "Qs", "62Ni:2+1", 0.0, "efm2", "evaluation", "evaluated", total=1.0,
                                  origin=results.reference("ENSDF", 2025, evaluation=nucleus.ENSDF)), tools.home)
    tools.add_data_point("62Ni", "Qs", "2+1", 5.0, "coulex-fit", uncertainty=3.0, reference=ALLMOND)
    tools.add_data_point("62Ni", "Qs", "2+1", -6.7, "shell-model", evidence="theory", reference=ALLMOND)
    c = tools.compare_results("62Ni", "Qs", "2+1")  # no difference in percent and no ratio against 0
    assert {r["evidence"]: r["against_reference"] for r in c["independent"]} == {
        "evaluated": "the reference", "measured": "+5 efm2, pull +1.58: agree",
        "theory": "-6.7 efm2 (no uncertainty: not judged)"}
    json.dumps([tools.list_nuclei(), tools.list_results("62Ni"), tools.list_results("62Ni", "Qs")])


def test_the_server_has_the_results_tools(tools):
    server = build_server(session=tools.session)
    assert "list_nuclei" in INSTRUCTIONS and "never mixed" in INSTRUCTIONS

    async def talk():
        async with Client(server) as c:
            listed = {t.name: t for t in (await c.list_tools()).tools}
            assert {name for name, _ in ResultTools.TOOLS} <= set(listed)
            assert listed["list_nuclei"].annotations.read_only_hint
            assert listed["compare_results"].annotations.read_only_hint
            assert not listed["add_data_point"].annotations.read_only_hint
            assert "own_measurement" in listed["add_data_point"].input_schema["properties"]
            # No experiment is open, and none is needed.
            out = _data(await c.call_tool("list_nuclei", {}))
            assert [n["nucleus"] for n in out["nuclei"]] == ["62Ni", "194Pt"]
            out = _data(await c.call_tool("add_data_point", {
                "nucleus": "62Ni", "observable": "tau", "key": {"energy": "1172.9 keV", "jpi": "2+"},
                "value": 2.01, "unit": "ps", "uncertainty": 0.12, "method": "dsam",
                "reference": {"authors": "Kenn et al.", "year": 2001}}))
            assert out["key"] == "62Ni:2+1" and out["added"] == "τ(2+₁) = 2.01 ± 0.12 ps"
            out = _data(await c.call_tool("compare_results",
                                          {"nucleus": "62Ni", "observable": "B(E2)", "key": KEY}))
            assert len(out["independent"]) == 1 and len(out["projected"]) == 1
            refused = await c.call_tool("add_data_point", {
                "nucleus": "62Ni", "observable": "tau", "key": "2+1", "value": 2.0, "method": "dsam",
                "evidence": "projected"})
            assert refused.is_error and "simulated run" in refused.content[0].text
            assert "where the point is from" in refused.content[0].text
            refused = await c.call_tool("snapshot_evaluation", {"nucleus": "62Ni"})
            assert refused.is_error and "no local copy of ENSDF" in refused.content[0].text
            assert not listed["import_data_points"].annotations.read_only_hint

    asyncio.run(talk())
    assert server.session.planner is None


def test_points_imported_from_a_file(tools, tmp_path):
    points = tmp_path / "points.json"
    points.write_text(json.dumps([
        {"observable": "tau", "key": "2+1", "value": 1.6, "unit": "ps", "plus": 0.4, "minus": 0.6,
         "method": "dsam", "reference": {"authors": "Ohnuma et al.", "year": 1978}},
        {"observable": "tau", "key": "2+1", "value": 2.0, "method": "dsam", "evidence": "evaluated",
         "reference": {"authors": "Someone", "year": 2020}}]), encoding="utf-8")
    out = tools.import_data_points(str(points), "62Ni")
    json.dumps(out)
    assert out["summary"] == "1 added, 0 already there, 1 refused" and out["shape"] == "points"
    assert out["added"][0]["as_stored"] == "τ(2+₁) = 1.60 +0.40 −0.60 ps" and out["added"][0]["where"] == "points[0]"
    assert "an evaluation, not a measurement" in out["refused"][0]["reasons"][0]
    again = tools.import_data_points(str(points), "62Ni")
    assert again["summary"] == "0 added, 1 already there, 1 refused"
    assert again["already_there"][0]["id"] == out["added"][0]["id"]


def test_the_tools_make_the_rows_of_an_experiment_the_app_never_opened(tools):
    folder = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=tools.session.root)
    run = folder.path / "runs" / "01-beam-12h"
    run.mkdir(parents=True)
    folder.setup().save(run / "setup.toml")
    (run / "summary.json").write_text("{}", encoding="utf-8")
    (run / "analysis.json").write_text(json.dumps({
        "settings": {}, "b_e2fm4": 930.0, "statistical": 0.02, "systematic": 0.07, "truth_e2fm4": 884.1,
        "saved": "2026-11-09T22:14:36+00:00"}), encoding="utf-8")
    # A tool that reads writes nothing: the rows are made once, when the tools are (the server starts).
    tools.compare_results("62Ni", "B(E2)", KEY)
    assert not (folder.path / results.RESULTS_FILE).exists()
    started_again = type(tools)(tools.session)
    assert (folder.path / results.RESULTS_FILE).exists()
    rows = started_again.compare_results("62Ni", "B(E2)", KEY)["projected"]
    assert any(r["value"].startswith("930") for r in rows)
