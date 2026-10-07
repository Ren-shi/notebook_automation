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


def test_gamma_tab():
    assert Planner.example("alpha_on_gold").gamma() == {
        "available": False, "reason": 'The setup has no excited state ([reaction] type = "coulex").'}
    g = Planner.example("coulex_ni58").gamma()
    assert g["available"] and g["max_safe_angle"] == 180.0
    assert set(g["rates"]) == {"CD", "DSSD-L", "DSSD-R"} and all(v > 0 for v in g["rates"].values())
    assert len(g["doppler"]) == 12
    labels = [c["label"] for c in Planner.example("coulex_ni58").kinematics()["curves"]]
    assert any("58Ni* 1454 keV" in lab for lab in labels)


def test_coulomb_excitation_is_editable(tmp_path):
    """Everything a Coulomb-excitation setup needs can be set from the planner (and so from the app)."""
    p = Planner.example("alpha_on_gold")
    # Switching on Coulomb excitation fills in E2 excitation of the target, and its first state from the local
    # ENSDF copy when there is one; without a copy the state must be given.
    if not p.set("reaction", "type", "coulex"):
        assert any("'energy' is missing" in x for x in p.problems) and any("'b_up' is missing" in x for x in p.problems)
    assert p.set("reaction", "energy", "0.547 MeV") or p.problems
    assert p.set("reaction", "b_up", "0.3 e2b2")
    assert p.experiment.excitation.multipolarity == "E2" and p.experiment.excitation.excite == "target"
    assert p.gamma()["available"] and p.gamma()["doppler"] == []
    # γ-ray detectors: add, edit, duplicate, remove.
    assert p.add_gamma_detector(name="Ge1", theta="90 deg", distance="100 mm", radius="30 mm")
    assert p.set("gamma detector 1", "theta", "120 deg")
    assert p.duplicate_gamma_detector("Ge1", phi="180 deg")
    assert [g.name for g in p.experiment.gamma_detectors] == ["Ge1", "Ge1-2"]
    # Every particle detector that records excitation events is paired with both γ detectors.
    pairs = {(r["particle_detector"], r["gamma_detector"]) for r in p.gamma()["doppler"]}
    seen = {d for d, _ in pairs}
    assert seen and pairs == {(d, g) for d in seen for g in ("Ge1", "Ge1-2")}
    assert not p.set("gamma detector 2", "radius", "-1 mm") and p.problems
    assert p.set("gamma detector 2", "radius", "30 mm")
    # Back to elastic and again: the state's settings come back.
    assert p.set("reaction", "type", "elastic") and p.experiment.excitation is None
    assert p.set("reaction", "type", "coulex") and p.experiment.excitation.energy_mev == pytest.approx(0.547)
    assert p.remove_gamma_detector(0) and p.remove_gamma_detector("Ge1-2")
    assert "gamma_detectors" not in p.draft
    p.save(tmp_path / "s.toml")
    assert Planner.load(tmp_path / "s.toml").experiment.excitation.b_up == "0.3 e2b2"


def test_gamma_detectors_in_the_geometry():
    g = Planner.example("coulex_ni58").geometry()
    assert [x["name"] for x in g["gamma_detectors"]] == ["Ge90", "Ge45", "Ge135", "Ge0"]
    ge = g["gamma_detectors"][0]
    # The outline is a circle of the crystal's radius around the detector centre, facing the target.
    r = np.linalg.norm(ge["outline"] - ge["centre"], axis=1)
    assert np.allclose(r, 35.0)
    assert np.allclose((ge["outline"] - ge["centre"]) @ (ge["centre"] / 120.0), 0.0, atol=1e-9)
    assert g["extent"] >= 1.25 * 150.0 - 1e-9  # 187.5 up to rounding (187.49999999999997 on some numpy builds)
    assert Planner.example("alpha_on_gold").geometry()["gamma_detectors"] == []


def test_particle_energies_for_coulomb_excitation():
    """The energies each particle detector sees, against the elastic kinematic factor by hand."""
    rows = Planner.example("coulex_ni58").gamma()["particles"]
    by = {(r["detector"], r["particle"], r["where"]): r for r in rows}
    assert {(d, w) for d, _, w in by} == {(d, w) for d in ("CD", "DSSD-L", "DSSD-R")
                                         for w in ("min", "centre", "max")}
    # The beam cannot reach backward angles as a recoil partner: no 58Ni rows in the CD.
    assert not any(k[0] == "CD" and k[1] == "recoil" for k in by)
    m1, m2, e = 16.0, 58.0, 30.0  # mass numbers are enough at the 1% level
    for key, r in by.items():
        assert r["difference_mev"] > 0  # exciting the state always costs energy at a fixed angle
        assert 0 < r["beta_excited"] < 0.05
        th = math.radians(r["theta_lab"])
        if key[1] == "ejectile":
            k = ((m1 * math.cos(th) + math.sqrt(m2**2 - (m1 * math.sin(th)) ** 2)) / (m1 + m2)) ** 2
        else:
            k = 4 * m1 * m2 / (m1 + m2) ** 2 * math.cos(th) ** 2
        assert r["elastic_mev"] == pytest.approx(k * e, rel=0.01)
    # Detecting the excited 58Ni itself: β follows from its own energy.
    r = by[("DSSD-L", "recoil", "centre")]
    m = 57.935 * 931.494
    assert r["beta_excited"] == pytest.approx(math.sqrt(r["excited_mev"] * (r["excited_mev"] + 2 * m))
                                              / (r["excited_mev"] + m), rel=1e-3)  # nuclear vs atomic mass


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
