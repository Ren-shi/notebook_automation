"""The Analysis reads the run and the gates (backlog item 68): B(E2) from a run with a gate made on the Data tab and
no simulation of its own, predicted against measured, and the alignment and multistep blocks on run data."""

import math

import pytest

from physim.nuclear import dataviews as dv
from physim.nuclear import workbench as wb
from physim.nuclear.orientation import simple_scheme
from physim.nuclear.planner import Planner


def test_a_ten_minute_run_gives_b_e2_back_with_a_data_tab_gate(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="10 min")
    assert p.run.scale == 1.0
    p.put_gate(dv.Gate("every particle", None, None, "any"))
    out = p.analyse(gate="every particle")
    assert out["available"], out.get("reason")
    assert out["settings"]["particle_gate"] == "all"
    b, truth = out["b_e2fm4"], out["truth_e2fm4"]
    assert math.isfinite(b) and abs(out["pull"]) < 2.0, (b, truth, out["pull"])
    assert "run's real events" in out["statistics_from"] and out["monte_carlo"] == 0.0
    assert out["whole_run_statistical"] == pytest.approx(out["statistical"])  # all 10 min are real
    # No simulation of its own: only the run was read.
    assert not any(isinstance(k, tuple) and k[0] in ("events", "gammas") for k in p._cache)


@pytest.mark.parametrize("key", [t[0] for t in wb.TEMPLATES])
def test_predicted_and_measured_rates_agree_for_every_template(key, tmp_path):
    p = Planner(wb.template(key))
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="1 min", budget="20 s")
    rows = p.compare()["rows"]
    rates = [r for r in rows if r["quantity"] == "Rate"]
    assert len(rates) == len(p.experiment.detectors)
    for r in rates:
        assert r["agree"], r
    for r in rows:
        if r["quantity"] == "Coincidences in the peak":
            assert r["agree"] is not False, r


def test_the_efficiency_is_measured_by_a_source_run(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    p.start_run("source", duration="2 h", source="152Eu", activity="100 kBq")
    p.start_run("beam", duration="10 s")
    eff = [r for r in p.compare()["rows"] if r["quantity"] == "Efficiency at the transition"]
    assert eff and all(r["measured"] is not None for r in eff)
    for r in eff:
        assert r["agree"], r


def _close_and_with_a_scheme() -> Planner:
    p = Planner.example("coulex_ni58")
    d = p.draft
    d["detectors"] = [x for x in d["detectors"] if x["name"] == "CD"]
    d["gamma_detectors"] = [gd for gd in d["gamma_detectors"] if gd["name"] != "Ge135"]
    for gd in d["gamma_detectors"]:
        gd["distance"] = "60 mm"
    d["beam"]["current"] = "2 pnA"
    exc = d["reaction"]
    d["levels"] = {"target": simple_scheme("58Ni", exc["energy"], "E2", exc["b_up"]).to_dict()}
    assert p._apply(d), p.problems
    return p


def test_alignment_and_multistep_read_the_run(tmp_path):
    p = _close_and_with_a_scheme()
    p.create_experiment(root=tmp_path)
    run = p.start_run("alignment", duration="5 min", offset_mm=1.0)
    out = p.alignment(fit=False)
    assert out["available"] and out["offset_mm"] == 1.0
    assert p.gamma_events() is run.gammas()  # the run's own γ rays
    view = dv.gamma_vs_crystal(run, "recoil", offset_mm=1.0)
    assert view["counts"].sum() > 0
    fit = p.fit_matrix_elements()
    (key, value), = fit["values"].items()
    unc = fit["uncertainties"][key]
    truth = p.experiment.levels["target"].matrix_element(0, 1, "E2").value
    assert math.isfinite(value) and unc > 0
    assert abs(abs(value) - abs(truth)) < 4 * unc + 0.1 * abs(truth), (value, unc, truth)
