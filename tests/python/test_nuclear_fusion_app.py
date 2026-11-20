"""The app redesign for fusion-evaporation: the kinds of experiment, the ancillary detectors in the setup file and
the scene, the channel studied, and what the fusion screens read from the planner (planning/app-redesign)."""

import pytest

from physim.nuclear import scene
from physim.nuclear import workbench as wb
from physim.nuclear.experiment import Experiment, SetupError
from physim.nuclear.planner import Planner


def _fusion(**changes) -> dict:
    d = Experiment.example("fusion_dy158").to_dict()
    for key, value in changes.items():
        if value is None:
            d.pop(key, None)
        else:
            d[key] = value
    return d


def test_fusion_needs_no_particle_detector_and_round_trips():
    d = _fusion(detectors=None)
    exp = Experiment.from_dict(d)
    assert exp.detectors == []
    again = Experiment.from_toml(exp.to_toml())
    assert again.detectors == [] and "detectors" not in exp.to_dict()
    # Elastic scattering still needs one.
    d = Experiment.example("alpha_on_gold").to_dict()
    d.pop("detectors")
    with pytest.raises(SetupError, match="no detectors"):
        Experiment.from_dict(d)


def test_ancillary_sections():
    d = _fusion(charged_particle_ball={"elements": 40, "radius": "60 mm", "material": "CsI"},
                recoil_separator={"theta_acceptance": "3 deg", "kind": "gas-filled"},
                gamma_calorimeter={"elements": 60, "radius": "120 mm", "material": "BGO", "fold_min": 6})
    exp = Experiment.from_dict(d)
    assert exp.charged_particle_ball.elements == 40 and exp.recoil_separator.kind == "gas-filled"
    assert exp.gamma_calorimeter.fold_min == 6
    again = Experiment.from_toml(exp.to_toml())
    assert again.recoil_separator.theta_acceptance == exp.recoil_separator.theta_acceptance
    with pytest.raises(SetupError, match="material must be 'Si' or 'CsI'"):
        Experiment.from_dict(_fusion(charged_particle_ball={"elements": 4, "radius": "60 mm", "material": "Ge"}))
    with pytest.raises(SetupError, match="theta_min must be below theta_max"):
        Experiment.from_dict(_fusion(charged_particle_ball={"elements": 4, "radius": "60 mm", "theta_min": "90 deg",
                                                            "theta_max": "30 deg"}))
    elastic = Experiment.example("alpha_on_gold").to_dict()
    elastic["recoil_separator"] = {"theta_acceptance": "3 deg"}
    with pytest.raises(SetupError, match="recoil_separator: .*fusion-evaporation"):
        Experiment.from_dict(elastic)


def test_studied_channel():
    d = _fusion()
    d["run"]["channel"] = "158Dy"
    assert Experiment.from_dict(d).run.channel == "158Dy"
    d["run"]["channel"] = "2n a"  # the particles in any order
    Experiment.from_dict(d)
    d["run"]["channel"] = "6n"
    with pytest.raises(SetupError, match="channel '6n' is not in the channel table"):
        Experiment.from_dict(d)
    elastic = Experiment.example("alpha_on_gold").to_dict()
    elastic["run"]["channel"] = "4n"
    with pytest.raises(SetupError, match="run: channel"):
        Experiment.from_dict(elastic)


def test_plan_beam_time_follows_the_studied_channel():
    p = Planner.example("fusion_dy158")
    pl = p.plan()
    rows = pl["residues"]["rows"]
    largest = max(range(len(rows)), key=lambda i: rows[i]["cross_section_mb"])
    assert pl["study"]["index"] == largest and pl["beam_time_needed_s"] == rows[largest]["beam_time_s"]
    assert pl["limiting_detector"] is None and not pl["study"]["chosen"]
    assert p.set("run", "channel", "5n")
    pl = p.plan()
    assert pl["study"]["channel"] == "5n" and pl["study"]["chosen"]
    assert pl["beam_time_needed_s"] == rows[2]["beam_time_s"]
    assert wb.steps(p)[1]["status"].startswith("5n:")


def test_ancillaries_survive_a_switch_of_reaction():
    p = Planner.example("fusion_dy158")
    assert p.add_ancillary("recoil_separator") and p.add_ancillary("gamma_calorimeter")
    assert p.set("reaction", "type", "elastic") and p.experiment.recoil_separator is None
    assert p.set("reaction", "type", "fusion-evaporation")
    assert p.experiment.recoil_separator is not None and p.experiment.gamma_calorimeter is not None
    assert p.remove_ancillary("gamma_calorimeter") and p.experiment.gamma_calorimeter is None
    # A plunger is a lifetime section, added like an ancillary (backlog: lifetimes, wave 1A).
    assert p.add_ancillary("plunger") and p.experiment.lifetime.method == "rdds"
    with pytest.raises(KeyError):
        p.add_ancillary("nope")


def test_kinds_and_templates():
    live = [k.key for k in wb.KINDS if k.live]
    assert live == ["coulex", "elastic", "fusion-evaporation", "lifetimes"]
    assert all(not k.templates for k in wb.KINDS if not k.live)
    for key, *_ , ancillaries in wb.FUSION_TEMPLATES:
        exp = wb.template(key)
        assert exp.fusion is not None and wb.kind_of(exp)[0] == "fusion-evaporation", key
        assert all(getattr(exp, name) is not None for name in ancillaries), key
    assert wb.template("fusion-blank").detectors == []
    assert wb.kind_of(Experiment.example("coulex_ni62_pt194"))[1] == "Coulex · beam"
    assert wb.kind_of(Experiment.example("alpha_on_gold"))[0] == "elastic"
    exp = wb.new_setup("16O", "80 MeV", "120Sn", "1 mg/cm2", "fusion-evaporation", detectors_from="fusion-wall")
    assert exp.neutron_wall is not None and exp.fusion.channels


def test_ancillary_cards():
    p = Planner(wb.template("fusion-all"))
    keys = [d.key for d in wb.decisions(p)]
    assert [k for k in keys if k.startswith("anc:")] == [f"anc:{k}" for k in wb.ANCILLARY_CARDS]
    assert keys[-1] == "run"
    share = wb.ancillary_numbers(p)["recoil_separator"]
    assert 0 < share["all"] < 1 and set(share["channels"]) == {"4n"}


def test_band_and_evaporation():
    p = Planner.example("fusion_dy158")
    b = p.band()
    assert b["available"] and b["residue"] == "158Dy"
    first = next(r for r in b["rows"] if r["initial"] == "2+")
    assert first["b_e2_wu"] == pytest.approx(186, rel=0.02)  # from T½ = 1.66 ns, α = 2.82
    assert all(r["rate_per_s"] > 0 for r in b["rows"])
    ev = p.evaporation()
    assert ev["channel"] == "4n" and ev["multiplicity"]["n"]["mean"] == pytest.approx(4.0, abs=0.2)
    assert 0 < ev["cone"]["half_angle_deg"] < 30 and "tof" not in ev
    p.add_ancillary("neutron_wall")
    assert "tof" in p.evaporation()


def test_scene_draws_the_ancillaries():
    exp = wb.template("fusion-all")
    shapes = scene.ancillaries(exp, cone_deg=10.0)
    groups = [s["group"] for s in shapes]
    assert groups.count("neutron_wall") == exp.neutron_wall.cells
    assert {"charged_particle_ball", "gamma_calorimeter", "recoil_separator", "residues"} <= set(groups)
    assert scene.ancillary_reach_mm(shapes) >= 600
    assert scene.ancillaries(Experiment.example("alpha_on_gold")) == []
    # The advice beside the scene speaks of fusion, not of Coulomb excitation.
    text = " ".join(scene.advice(exp))
    assert "excitation is strongest" not in text and "monitor" in text


def test_fusion_warnings_use_one_barrier():
    p = Planner.example("fusion_dy158")
    texts = [w.text for w in p.warnings()]
    assert not any("below the Coulomb barrier" in t for t in texts)
    assert all(w.level == "note" for w in p.warnings() if w.text.startswith("For the beam monitor"))


def test_start_page_helpers():
    from physim.nuclear.app.start import experiment_kind, fusion_preview

    text = fusion_preview("12C", "60 MeV", "150Nd")
    assert text.startswith("162Dy*") and "above the fusion barrier" in text
    assert fusion_preview("nonsense", "?", "?") == ""
    assert experiment_kind("/no/such/folder") == ""


def test_fusion_figures():
    pytest.importorskip("plotly")
    from physim.nuclear.app import figures

    p = Planner.example("fusion_dy158")
    figures.figure_fusion_gammas(p, events=20_000).to_json()
    figures.figure_fusion_gammas(p, "Ge90", corrected=False, only="4n", events=20_000).to_json()
    for f in figures.figure_evaporation(p.evaporation()):
        f.to_json()
    figures.figure_doppler(p, 320.5).to_json()
    figures.figure_band(p.band()).to_json()
