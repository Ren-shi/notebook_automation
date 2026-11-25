"""Data points imported from a file, the report of a nucleus and the route that holds the reference (wave 3 of the
results pages), on temporary folders with a handful of synthetic ⁶²Ni rows and a cut-down literature file. No run
is simulated and nothing is computed but the views."""

import json
from pathlib import Path

import pytest

from physim.nuclear import Experiment, ensdf, logbook, nucleus, results
from physim.nuclear.runs import ExperimentFolder

KEY = "62Ni:2+1>0+1:E2"
NI62 = ((0.0, 0, 1, 1), (1172.91, 4, 1, 1), (2048.55, 0, 1, 2), (2336.46, 8, 1, 1))
LITERATURE = Path(__file__).resolve().parents[1] / "reference" / "nuclear" / "literature_ni62_cut.json"
ENSDF_ORIGIN = results.reference("ENSDF, adopted levels and γ rays", 2025, where="NDS 204, 1 (2025)",
                                 evaluation=nucleus.ENSDF)


@pytest.fixture
def lab(tmp_path, monkeypatch):
    """Empty experiments and nuclei folders, no local ENSDF copy, the ⁶²Ni levels to anchor to, and no γ ray to
    convert B(E2) and τ through: the same on every computer."""
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path / "experiments"))
    monkeypatch.setenv(results.ENV, str(tmp_path / "nuclei"))
    monkeypatch.setenv(ensdf.ENV, str(tmp_path / "no-ensdf"))
    monkeypatch.setattr(results, "ensdf_levels", lambda nuclide: NI62 if nuclide == "62Ni" else ())
    monkeypatch.setattr(nucleus, "decay", lambda nuclide, level_key: (None, None, "no γ ray known"))
    return tmp_path


def _points(lab) -> list:
    return results.read(lab / "nuclei" / "62Ni" / results.POINTS_FILE)


# -- data points from a file ----------------------------------------------------------------------------------------


def test_a_literature_file_is_imported_once_with_its_methods_and_refusals(lab):
    found = nucleus.import_points(LITERATURE, "62Ni")
    assert found["shape"] == "literature"
    added = {x["where"]: x["row"] for x in found["added"]}
    assert len(added) == 9 and not found["already"]
    allmond = added["quantities.BE2_up[0]"]
    assert (allmond.key, allmond.value, allmond.unit, allmond.direction, allmond.method, allmond.total) == (
        KEY, 0.0906, "e2b2", "up", "coulex-fit", 0.0037)
    assert allmond.detail.startswith("Coulomb excitation, inverse kinematics") and "GOSIA" in allmond.detail
    assert allmond.detail.endswith("a note of the file, kept in the row's detail")  # "superseded by ..." stays
    assert nucleus.source(allmond) == "J. M. Allmond et al. 2014" and allmond.origin["doi"].endswith("034309")
    # A lifetime stays a τ, as the paper gives it; stat and sys, and an asymmetric uncertainty, are kept.
    kenn, ohnuma = added["quantities.tau_2plus[0]"], added["quantities.tau_2plus[1]"]
    assert (kenn.observable, kenn.unit, kenn.stat, kenn.sys, kenn.method) == ("tau", "ps", 0.07, 0.1, "dsam")
    assert (ohnuma.plus, ohnuma.minus, ohnuma.total) == (0.4, 0.6, None) and nucleus.source(kenn) == "O. Kenn et al. 2001"
    assert added["quantities.tau_2plus[2]"].method == "inelastic-neutron"  # DSAM after (n,n'γ)
    assert added["quantities.tau_2plus[3]"].method == "nrf"
    assert (added["quantities.Qs_2plus[0]"].unit, added["quantities.Qs_2plus[0]"].method) == ("eb", "coulex-fit")
    model = added["quantities.shell_model[0]"]
    assert (model.evidence, model.method, model.unit, model.direction, model.key) == (
        "theory", "shell-model", "W.u.", "down", KEY)
    assert (added["quantities.other[0]"].observable, added["quantities.other[0]"].method) == ("mu", "moment")
    # Another evaluation is taken, marked as an average of measurements.
    raman = added["quantities.BE2_up[2]"]
    assert (raman.method, raman.evidence, raman.origin["not_independent_of"]) == (
        "evaluation", "evaluated", "measurements")
    # ENSDF's values and a quantity without a value are refused, each with its reason; the others go on.
    refused = {x["where"]: " ".join(x["reasons"]) for x in found["refused"]}
    assert set(refused) == {"quantities.BE2_up[1]", "quantities.BE2_down[0]", "quantities.shell_model[1]"}
    assert "with the snapshot of the evaluation" in refused["quantities.BE2_down[0]"]
    assert refused["quantities.shell_model[1]"].startswith("the file gives no value for it")
    assert len(_points(lab)) == 9
    # The same file again adds nothing: every point is there already, with its row.
    again = nucleus.import_points(LITERATURE, "62Ni")
    assert not again["added"] and len(again["already"]) == 9 and len(again["refused"]) == 3
    assert {x["id"] for x in again["already"]} == {r.id for r in added.values()} and len(_points(lab)) == 9


def test_points_from_json_and_csv_files(lab):
    good = {"observable": "B(E2)", "key": "2+1 > 0+1", "value": 900, "unit": "e2fm4", "direction": "up",
            "plus": 60, "minus": 50, "method": "coulex-fit", "reference": {"authors": "Roe", "year": 2020}}
    (lab / "points.json").write_text(json.dumps({"nucleus": "62Ni", "points": [
        good, {**good, "value": "many", "method": None}, {**good, "nucleus": "58Ni"}, {**good, "colour": "red"}]}),
        encoding="utf-8")
    found = nucleus.import_points(lab / "points.json", "62Ni")
    assert found["shape"] == "points" and [x["where"] for x in found["added"]] == ["points[0]"]
    row = found["added"][0]["row"]
    assert (row.plus, row.minus, row.total) == (60.0, 50.0, None)
    reasons = {x["where"]: " ".join(x["reasons"]) for x in found["refused"]}
    assert "the value must be a number" in reasons["points[1]"] and "pick the method: " in reasons["points[1]"]
    assert "Coulex fitted to all orders: coulex-fit" in reasons["points[1]"] and "None" not in reasons["points[1]"]
    assert reasons["points[2]"] == "the point is of 58Ni, not of 62Ni"
    assert "'colour' is not a field of a data point" in reasons["points[3]"]
    (lab / "points.csv").write_text(
        "observable,key,value,unit,uncertainty,method,authors,year,date,note\n"
        "tau,2+1,2.01,ps,0.12,dsam,Kenn et al.,2001,,\n"
        "tau,2+1,2.2,ps,0.2,rdds,,,2027-06-14,runs 12 to 31\n"
        "tau,2+1,2.01,ps,0.12,dsam,Kenn et al.,2001,,\n", encoding="utf-8")
    found = nucleus.import_points(lab / "points.csv", "62Ni")
    assert found["shape"] == "csv" and len(found["added"]) == 2 and not found["refused"]
    assert [x["where"] for x in found["already"]] == ["line 4"]  # the same point twice in one file
    own = found["added"][1]["row"]
    assert (own.method, own.origin["kind"], own.origin["date"]) == ("rdds", "own measurement", "2027-06-14")
    assert nucleus.method_of("Inelastic electron scattering (e,e'), 45-65 MeV") == ("electron-scattering", "measured")
    assert nucleus.method_of("ENSDF 2025 adopted; 4+1 T1/2 = 0.88(21) ps")[1] == "evaluated"
    with pytest.raises(ValueError, match="holds no data points"):
        (lab / "other.json").write_text("{\"a\": 1}", encoding="utf-8")
        nucleus.import_points(lab / "other.json")


# -- the route that holds the reference -----------------------------------------------------------------------------


def test_the_route_of_the_reference_is_flagged(lab):
    results.add_points([results.Row("62Ni", "B(E2)", KEY, 884.1, "e2fm4", "evaluation", "evaluated",
                                    direction="up", total=13.1, origin=ENSDF_ORIGIN)])
    view = nucleus.observable_view("62Ni", "B(E2)", KEY)
    routes = {r["method"]: r for r in view["routes"]}
    assert routes["evaluation"]["is_reference"] and routes["evaluation"]["others"] == 0
    assert not routes["coulex-fit"]["is_reference"] and routes["coulex-fit"]["others"] == 0
    # The same rows given to both views: nothing is read twice, and they say the same.
    rows = results.of_nuclide("62Ni")
    assert nucleus.view("62Ni", rows=rows)["table"] == nucleus.view("62Ni")["table"]
    assert nucleus.observable_view("62Ni", "B(E2)", KEY, rows=rows)["routes"] == view["routes"]


# -- the report of a nucleus ----------------------------------------------------------------------------------------


def _files(folder: Path) -> dict:
    return {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in folder.rglob("*") if p.is_file()}


def test_the_report_of_a_nucleus(lab):
    coulex = ExperimentFolder.create(Experiment.example("coulex_ni62_pt194"), root=lab / "experiments")
    results.append(coulex.path, results.Row(
        "62Ni", "B(E2)", KEY, 918.6, "e2fm4", "coulex-yield", "projected", direction="up", stat=12.9, sys=63.4,
        truth=884.1, origin=results.from_run(coulex.path.name, "03-beam-24h")))
    results.add_points([results.Row("62Ni", "B(E2)", KEY, 884.1, "e2fm4", "evaluation", "evaluated",
                                    direction="up", total=13.1, origin=ENSDF_ORIGIN),
                        results.Row("62Ni", "E_level", "62Ni:2+1", 1172.91, "keV", "evaluation", "evaluated",
                                    total=0.03, origin=ENSDF_ORIGIN)])
    nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 906.0, direction="up", uncertainty=37.0,
                           method="coulex-fit", reference={"authors": "Allmond et al.", "year": 2014})
    before = _files(lab)
    page = logbook.nucleus_html("62Ni")
    assert _files(lab) == before  # building it writes nothing
    assert page.startswith("<!doctype html>") and "<style>" in page and "⁶²Ni: the report of the nucleus" in page
    # Results first, then the experiments and the history; a quantity with the evaluation alone has no section.
    order = [page.index(x) for x in ('<h2 id="q1">1. B(E2; 0⁺₁ → 2⁺₁)↑', "The values, side by side",
                                     "The routes to it", "Still to run", '<h2 id="experiments">A. The experiments',
                                     '<h2 id="history">B. The history of the rows')]
    assert order == sorted(order) and page.count('<h2 id="q') == 1
    assert "B(E2)↓ in W.u. = 12.13 ± 0.18 W.u.↓" in page and "e²fm⁴" in page and "e2fm4" not in page
    # The measured point is pulled against the reference; the simulated run is closed, never pulled.
    independent = page[page.index("<h3>Independent evidence"):page.index("<h3>Projected")]
    assert "pull +" in independent and "agree" in independent
    projected = page[page.index("<th>Simulated run</th>"):]
    projected = projected[:projected.index("</table>")]
    assert "closes" in projected and "pull" not in projected
    assert "62Ni on 194Pt" in page and "Allmond et al. 2014" in page.split('id="history"')[1]


def test_the_decay_curve_of_a_ddcm_is_drawn_against_the_shifted_share():
    pytest.importorskip("plotly")
    from physim.nuclear.app import figures as f

    points = [{"distance_um": d, "unshifted": u, "shifted": 1 - u, "error": 0.02}
              for d, u in ((10.0, 0.2), (100.0, 0.5), (1000.0, 0.9))]
    ddcm = {"points": points, "curve": {"x": [10.0, 100.0, 1000.0], "y": [0.8, 0.5, 0.1],
                                        "kind": "speed-weighted shifted intensity (DDCM fit)"}}
    fig = f.figure_decay_curve(ddcm)
    assert list(fig.data[0].y) == pytest.approx([0.8, 0.5, 0.1]) and "shifted share" in fig.data[1].name
    assert fig.layout.yaxis.title.text.startswith("shifted share")
    curve = dict(ddcm["curve"], kind="unshifted share (decay-curve fit)")
    fig = f.figure_decay_curve({"points": points, "curve": curve})
    assert list(fig.data[0].y) == pytest.approx([0.2, 0.5, 0.9]) and fig.data[1].name == "fit"
    other = dict(ddcm["curve"], kind="decays before each foil (DDCM fit), against x₁ + x₂ β/β′")
    fig = f.figure_decay_curve({"points": points, "curve": other})
    assert fig.data[1].yaxis == "y2" and fig.layout.yaxis2.title.text == other["kind"]
