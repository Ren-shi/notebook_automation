"""The planner app's model (backlog item 41): every tab renders for every example, edits are checked, sweeps."""

import math

import numpy as np
import pytest

from physim.nuclear import SetupError
from physim.nuclear.planner import EXPLAIN, TABS, Planner


@pytest.mark.parametrize("name", Planner.examples())
def test_every_tab_renders_for_every_example(name):
    p = Planner.example(name)
    for tab in TABS:
        out = p.tab(tab)
        assert isinstance(out, dict) and out, tab
        e = p.explain(tab)
        assert {"title", "formula", "assumptions", "limits", "register"} <= set(e)
    geo = p.geometry()
    assert len(geo["detectors"]) == len(p.experiment.detectors)
    rows = p.rates()["rows"]
    assert all(r["rate_per_s"] >= 0 for r in rows)
    spectra = p.spectra(events=20_000)["spectra"]
    assert set(spectra) == {d.name for d in p.experiment.detectors}
    orbits = p.trajectories()["orbits"]
    assert all(o["xy"].shape[1] == 2 for o in orbits)
    with pytest.raises(KeyError):
        p.tab("nope")


def test_rates_tab_matches_the_beam_time_arithmetic():
    p = Planner.example("alpha_on_gold")
    t = p.rates()
    for row in t["rows"]:
        assert row["counts_in_run"] == pytest.approx(row["rate_per_s"] * t["beam_time_s"])
        assert row["beam_time_s"] == pytest.approx(t["counts_wanted"] / row["rate_per_s"])
        assert row["relative_error"] == pytest.approx(1 / math.sqrt(row["counts_in_run"]))
    assert sum(t["strips"]["A45"].values()) == pytest.approx(p._rates().rate("A45"))


def test_editing_keeps_the_last_valid_setup():
    p = Planner.example("alpha_on_gold")
    assert not p.set("beam", "energy", "-3 MeV")
    assert p.problems and "energy" in p.problems[0]
    assert p.experiment.beam.energy_mev == 5.5  # results still show the last valid setup
    assert p.warnings()[0].level == "error"
    with pytest.raises(SetupError):
        p.save("never-written.toml")
    assert p.set("beam", "energy", "6 MeV") and not p.problems
    assert p.experiment.beam.energy_mev == 6.0
    assert p.set("backing", "material", "C") is False  # a backing needs a thickness too
    assert p.set("backing", "thickness", "10 ug/cm2")
    assert p.experiment.target.backing.material == "C"
    assert p.set("backing", "material", None) is False and p.set("backing", "thickness", None)
    assert p.experiment.target.backing is None
    assert p.set("title", "", "My run")
    assert p.experiment.title == "My run"


def test_detector_table_operations(tmp_path):
    p = Planner.example("alpha_on_gold")
    n = len(p.experiment.detectors)
    assert p.duplicate_detector("A45", theta="50 deg")
    assert [d.name for d in p.experiment.detectors][-1] == "A45-2"
    assert p.experiment.detectors[-1].theta.to("deg") == 50.0
    assert p.duplicate_detector("A45")
    assert p.experiment.detectors[-1].name == "A45-3"
    assert p.set("detector 1", "distance", "120 mm")
    assert p.experiment.detectors[0].distance.to("mm") == 120.0
    assert p.remove_detector("A45-2") and p.remove_detector(len(p.experiment.detectors) - 1)
    assert len(p.experiment.detectors) == n
    assert p.add_detector(name="A75", shape="circle", theta="75 deg", distance="100 mm", radius="2.5 mm",
                          thickness="300 um")
    assert "A75" in p.rates()["strips"]
    with pytest.raises(KeyError):
        p.remove_detector("nope")
    path = tmp_path / "setup.toml"
    p.save(path)
    again = Planner.load(path)
    assert [d.name for d in again.experiment.detectors] == [d.name for d in p.experiment.detectors]


def test_sweeps():
    p = Planner.example("alpha_on_gold")
    s = p.sweep("beam energy", ["4 MeV", "6 MeV"], "rate", detector="A45")
    # Rutherford: rate ∝ 1/E² (slightly more, as the beam slows in the target).
    assert s["y"][0] / s["y"][1] == pytest.approx((6 / 4) ** 2, rel=0.02)
    t = p.sweep("target thickness", ["0.25 mg/cm2", "0.5 mg/cm2"], detector="A45")
    # Twice the nuclei, and the beam slows a little more, which raises the cross section by ~1%.
    assert t["y"][1] / t["y"][0] == pytest.approx(2.0, rel=0.02)
    a = p.sweep("detector angle", ["40 deg", "50 deg"], "peak energy", detector="A45")
    assert a["y"][0] > a["y"][1] and a["x"] == [40.0, 50.0]
    b = p.sweep("beam energy", ["5 MeV"], "beam time", detector="A45")
    assert b["y"][0] > 0
    w = p.sweep("target thickness", ["0.1 mg/cm2", "1 mg/cm2"], "peak width", detector="A60")
    assert w["y"][1] > w["y"][0]
    with pytest.raises(ValueError):
        p.sweep("magnet", ["1 T"])
    with pytest.raises(ValueError):
        p.sweep("beam energy", ["5 MeV"], "colour")


def test_warnings_banner_and_explain_texts():
    p = Planner.example("oxygen_on_lead_array")
    w = p.warnings()
    assert any(x.level == "warning" and "Coulomb barrier" in x.text for x in w)
    assert any("pile-up" in x.text for x in w)
    levels = [x.level for x in w]
    assert levels == sorted(levels, key=["error", "warning", "note"].index)
    assert set(EXPLAIN) == set(TABS)


def test_energy_loss_tab():
    p = Planner.example("oxygen_on_lead_array")
    el = p.energy_loss()
    target, backing = el["layers"]
    assert target["energy_out_mev"] == pytest.approx(backing["energy_in_mev"])
    assert target["loss_mev"] > 0 and backing["loss_mev"] > 0
    assert np.all(np.diff(el["beam_energy_mev"]) <= 1e-12)
    assert el["depth_mg_cm2"][-1] == pytest.approx(0.22)
    for d in el["detectors"]:
        assert d["after_dead_layer_mev"] < d["ejectile_energy_mev"] < d["punch_through_mev"]
