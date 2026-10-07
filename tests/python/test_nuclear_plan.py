"""The Plan tab's model (backlog item 65): the answers without any button pressed, and the predictions kept with the
next run."""

import math

import pytest

from physim.nuclear import workbench as wb
from physim.nuclear.planner import Planner


@pytest.mark.parametrize("key", [t[0] for t in wb.TEMPLATES])
def test_plan_answers_for_every_template_without_simulating(key):
    p = Planner(wb.template(key))
    pl = p.plan()
    assert pl["beam_time_s"] > 0 and pl["detectors"]
    for d in pl["detectors"]:
        assert d["rate_per_s"] >= 0 and 0 <= d["dead_time_fraction"] < 1
        assert len(d["safe_rings"]) + len(d["unsafe_rings"]) == d["rings"] >= 1
    if pl["counts_wanted"]:
        assert pl["beam_time_needed_s"] > 0 and pl["limiting_detector"] in {d["detector"] for d in pl["detectors"]}
        assert pl["enough"] == (pl["beam_time_needed_s"] <= pl["beam_time_s"])
    exp = p.experiment
    if exp.excitation is not None and exp.gamma_detectors:
        assert len(pl["gamma_detectors"]) == len(exp.gamma_detectors)
        assert all(0 < g["efficiency"] < 1 for g in pl["gamma_detectors"])
        # The matrix adds up to each particle detector's coincidence rate.
        for d in pl["detectors"]:
            assert sum(pl["coincidences"][d["detector"]].values()) == pytest.approx(d["coincidence_per_s"], rel=1e-9)
        shares = [d["excitation_share"] for d in pl["detectors"]]
        assert sum(shares) == pytest.approx(1.0)
    # Nothing was simulated to get here.
    assert not any(isinstance(k, tuple) and k[0] in ("events", "gammas") for k in p._cache)


def test_dead_time_and_unsafe_rings():
    p = Planner.example("coulex_ni58")
    p.set("run", "dead_time", "10 us")
    pl = p.plan()
    cd = next(d for d in pl["detectors"] if d["detector"] == "CD")
    assert cd["dead_time_fraction"] == pytest.approx(cd["rate_per_s"] * 1e-5 / (1 + cd["rate_per_s"] * 1e-5))
    total = sum(d["rate_per_s"] for d in pl["detectors"])
    assert pl["live_fraction"] == pytest.approx(1 / (1 + 1e-5 * total))
    p.set("beam", "energy", "64 MeV")  # above the safe energy at backward angles
    assert any(d["unsafe_rings"] for d in p.plan()["detectors"])


def test_the_next_run_keeps_the_plan(tmp_path):
    p = Planner.example("alpha_on_gold")
    p.create_experiment(root=tmp_path)
    before = p.plan()
    run = p.start_run("beam", duration="5 s")
    kept = run.summary["plan"]
    assert [d["detector"] for d in kept["detectors"]] == [d["detector"] for d in before["detectors"]]
    for a, b in zip(kept["detectors"], before["detectors"]):
        assert a["rate_per_s"] == pytest.approx(b["rate_per_s"])
    assert kept["beam_time_needed_s"] == pytest.approx(before["beam_time_needed_s"])
    p.set("beam", "energy", "6 MeV")
    assert p.plan()["detectors"][0]["rate_per_s"] != pytest.approx(kept["detectors"][0]["rate_per_s"])
    assert p.folder.load_run(1).summary["plan"]["detectors"][0]["rate_per_s"] == pytest.approx(
        before["detectors"][0]["rate_per_s"])
    assert not math.isnan(kept["live_fraction"])
