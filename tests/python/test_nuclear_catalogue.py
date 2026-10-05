"""The detector catalogue, clover crystals, dead material and the chamber (backlog item 51)."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment, SetupError, catalogue
from physim.nuclear.detectors import Array
from physim.nuclear.events import simulate
from physim.nuclear.gamma import doppler_table
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates

S3 = {"name": "S3", "model": "S3", "theta": "180 deg", "distance": "30 mm", "resolution": "50 keV",
      "threshold": "500 keV"}
CLOVER = {"name": "Clo", "model": "clover", "theta": "135 deg", "phi": "90 deg", "distance": "200 mm"}
LABR = {"name": "La", "model": "LaBr3_2x2", "theta": "90 deg", "phi": "0 deg", "distance": "150 mm"}


def setup(detectors=(S3,), gammas=(CLOVER, LABR), **sections):
    d = Experiment.example("coulex_ni58").to_dict()
    d["detectors"] = [dict(x) for x in detectors]
    d["gamma_detectors"] = [dict(x) for x in gammas]
    if not gammas:
        del d["gamma_detectors"]
    d.update(sections)
    return d


# -- the catalogue --------------------------------------------------------------------------------------------------


def test_models_and_their_sources():
    assert catalogue.names("particle") == ["S3"]
    assert catalogue.names("gamma") == ["clover", "LaBr3_2x2"]
    for name in catalogue.names():
        m = catalogue.model(name)
        assert m.title and m.sources and m.notes, name
        # Every field marked as typical is a field of the model or one of its blocking parts.
        known = set(m.fields) | {b["name"] for b in m.blocking}
        assert set(m.typical) <= known, name
    with pytest.raises(ValueError, match="'S4' is not a particle detector model; the models are S3"):
        catalogue.model("S4", "particle")
    with pytest.raises(ValueError, match="not a γ-ray detector model"):
        catalogue.model("S3", "gamma")


def test_s3_matches_its_data_sheet():
    m = catalogue.model("S3")
    f = m.fields
    # Micron's data sheet: 24 rings, 32 sectors, active area 22 to 70 mm in diameter, chip 20 to 76 mm.
    assert (f["rings"], f["sectors"], f["inner_radius"], f["outer_radius"]) == (24, 32, "11 mm", "35 mm")
    assert m.facts["ring_pitch"] == "1 mm" and m.facts["ring_strip_width"] == "0.886 mm"
    assert (35 - 11) / 24 == 1.0
    assert (m.facts["chip_inner_diameter"], m.facts["chip_outer_diameter"]) == ("20 mm", "76 mm")
    assert {"thickness", "dead_layer", "board"} == set(m.typical)


def test_clover_and_labr3_match_their_documents():
    c = catalogue.model("clover").fields
    assert (c["crystals"], c["crystal_diameter"], c["crystal_length"]) == (4, "50 mm", "70 mm")
    assert catalogue.model("clover").facts["crystal_gap_max"] == "0.7 mm"
    la = catalogue.model("LaBr3_2x2").fields
    assert (la["crystals"], la["crystal_diameter"], la["crystal_length"]) == (1, "50.8 mm", "50.8 mm")
    # 2.1% of 1332 keV.
    assert float(la["resolution"].split()[0]) == pytest.approx(0.021 * 1332, abs=0.5)


# -- a model in a setup file ----------------------------------------------------------------------------------------


def test_a_model_fills_in_the_fields_and_can_be_overridden():
    exp = Experiment.from_dict(setup())
    s3 = exp.detectors[0]
    assert (s3.model, s3.shape, s3.rings, s3.sectors) == ("S3", "annular", 24, 32)
    assert str(s3.thickness) == "1000 um"
    thin = Experiment.from_dict(setup([dict(S3, thickness="300 um", rings=12)]))
    assert str(thin.detectors[0].thickness) == "300 um" and thin.detectors[0].rings == 12
    # Saved with every field written out, and read back the same.
    text = exp.to_toml()
    assert 'model = "S3"' in text and "rings = 24" in text and 'model = "clover"' in text
    assert Experiment.from_toml(text).to_toml() == text


def test_unknown_models_and_missing_sizes_are_named():
    with pytest.raises(SetupError) as e:
        Experiment.from_dict(setup([dict(S3, model="S9")],
                                   [dict(CLOVER, model="S3"), {"theta": "90 deg", "distance": "100 mm"},
                                    {"theta": "90 deg", "distance": "100 mm", "crystals": 4,
                                     "crystal_diameter": "50 mm"},
                                    {"theta": "90 deg", "distance": "100 mm", "crystals": 3,
                                     "crystal_diameter": "50 mm"}]))
    text = "\n".join(e.value.problems)
    for expected in ("detector 1 (S3): model 'S9' is not a particle detector model",
                     "gamma detector 1: model 'S3' is not a γ-ray detector model",
                     "gamma detector 2: give 'radius', or a 'model' (one of clover, LaBr3_2x2)",
                     "gamma detector 3: four crystals need 'crystal_pitch'",
                     "gamma detector 4: crystals must be 1 or 4"):
        assert expected in text, expected


def test_setups_without_models_are_unchanged():
    exp = Experiment.example("coulex_ni58")
    assert all(d.model is None and d.blocking() == [] for d in exp.detectors)
    assert Array.from_experiment(exp).blockers == []
    g = exp.gamma_detectors[0]
    assert g.elements() == [("", pytest.approx(tuple(120 * c for c in g.direction())), 35.0)]
    assert g.geometric_efficiency() == pytest.approx((1 - math.cos(math.atan2(35, 120))) / 2)


# -- S3 geometry ----------------------------------------------------------------------------------------------------


def test_each_s3_ring_covers_the_solid_angle_of_a_hand_calculation():
    g = Array.from_experiment(Experiment.from_dict(setup()))["S3"]
    z = 30.0
    per_ring = np.zeros(24)
    for (ring, _), omega in g.segment_solid_angles().items():
        per_ring[ring] += omega
    for k in range(24):
        r1, r2 = 11.0 + k, 12.0 + k
        hand = 2 * math.pi * (z / math.hypot(z, r1) - z / math.hypot(z, r2)) * 1e3  # msr
        assert per_ring[k] == pytest.approx(hand, rel=1e-3), k
    assert g.solid_angle() == pytest.approx(per_ring.sum(), rel=1e-6)


def test_the_s3s_own_board_does_not_hide_it():
    exp = Experiment.from_dict(setup())
    array = Array.from_experiment(exp)
    assert [b.name for b in array.blockers][:2] == ["S3 (chip edge)", "S3 (board)"]
    assert array.shadowing() == {} and array.warnings() == []
    # The beam passes through the hole in the board.
    bare = Experiment.from_dict(setup([{k: v for k, v in S3.items() if k != "model"}
                                       | {"shape": "annular", "inner_radius": "11 mm", "outer_radius": "35 mm",
                                          "rings": 24, "sectors": 32, "thickness": "1000 um",
                                          "dead_layer": "0.5 um"}]))
    assert Rates(exp).rate("S3") == pytest.approx(Rates(bare).rate("S3"), rel=1e-12)


def test_a_board_hides_a_detector_behind_it():
    # A pad 10 mm behind the S3, at 44 mm from the axis: behind the board (10 to 50 mm), outside the silicon.
    pad = {"name": "P", "shape": "circle", "position": ["44 mm", "0 mm", "-40 mm"], "facing": [0.0, 0.0, 1.0],
           "radius": "3 mm", "thickness": "300 um", "threshold": "500 keV"}
    exp = Experiment.from_dict(setup([S3, pad]))
    array = Array.from_experiment(exp)
    assert array.shadowing()[("S3 (board)", "P")] == pytest.approx(1.0)
    assert any("S3 (board) hides 100% of P" in w for w in array.warnings())
    assert Rates(exp).rate("P") == 0.0
    events = simulate(exp, 200_000, seed=3)
    assert events.detectors == ["S3", "P"]
    assert events.select("P", counted=False).sum() == 0 and events.select("S3").sum() > 1000
    # Without the model, nothing is in the way.
    no_board = Experiment.from_dict(setup([{k: v for k, v in S3.items() if k != "model"}
                                           | catalogue.model("S3").fields, pad]))
    assert Rates(no_board).rate("P") > 0


# -- γ-ray detectors ------------------------------------------------------------------------------------------------


def test_clover_crystals_sit_around_the_axis():
    exp = Experiment.from_dict(setup())
    clover, labr = exp.gamma_detectors
    crystals = clover.elements()
    assert [c[0] for c in crystals] == ["A", "B", "C", "D"] and all(c[2] == 25.0 for c in crystals)
    centres = np.array([c[1] for c in crystals])
    axis = np.array(clover.direction())
    # Their mean is on the axis at the stated distance, and each is half a diagonal of the pitch from it.
    assert centres.mean(axis=0) == pytest.approx(200 * axis)
    assert np.linalg.norm(centres - 200 * axis, axis=1) == pytest.approx(45 / math.sqrt(2))
    assert (centres - 200 * axis) @ axis == pytest.approx(0, abs=1e-9)
    assert np.linalg.norm(centres[0] - centres[1]) == pytest.approx(45.0)
    assert labr.elements()[0][2] == 25.4
    # Four crystals cover four times one crystal's share, to first order.
    one = (1 - math.cos(math.atan2(25, math.hypot(200, 45 / math.sqrt(2))))) / 2
    assert clover.geometric_efficiency() == pytest.approx(4 * one, rel=1e-9)


def test_doppler_table_has_a_row_per_crystal():
    rows = doppler_table(Experiment.from_dict(setup()), order=8, depth_points=2)
    assert [r["gamma_detector"] for r in rows] == ["Clo A", "Clo B", "Clo C", "Clo D", "La"]
    shifts = {r["gamma_detector"]: r["shift_kev"] for r in rows}
    # The target recoils forward, so the clover at 135° sees a red shift; crystals nearer the beam axis see more.
    assert all(shifts[k] < 0 for k in ("Clo A", "Clo B", "Clo C", "Clo D"))
    assert shifts["Clo A"] == pytest.approx(shifts["Clo B"], rel=1e-6)
    assert shifts["Clo A"] != pytest.approx(shifts["Clo C"], rel=1e-3)
    # One crystal has a smaller opening angle than the whole clover, so its Doppler width is smaller.
    disc = dict(CLOVER, name="Clo")
    disc.pop("model")
    envelope = math.tan(math.radians(Experiment.from_dict(setup()).gamma_detectors[0].half_angle_deg())) * 200
    whole = doppler_table(Experiment.from_dict(setup(gammas=[dict(disc, radius=f"{envelope} mm")])), order=8,
                          depth_points=2)[0]
    assert all(r["doppler_fwhm_kev"] < whole["doppler_fwhm_kev"] for r in rows[:4])


def test_a_clover_housing_hides_a_detector_in_the_rates_and_the_events():
    pad = {"name": "P", "shape": "circle", "theta": "45 deg", "phi": "0 deg", "distance": "100 mm",
           "radius": "5 mm", "thickness": "300 um", "threshold": "500 keV"}
    in_front = dict(CLOVER, theta="45 deg", phi="0 deg", distance="60 mm")
    hidden = Experiment.from_dict(setup([pad], [in_front]))
    clear = Experiment.from_dict(setup([pad], [CLOVER]))
    assert Rates(clear).rate("P") > 0 and Rates(hidden).rate("P") == 0.0
    assert any("Clo (housing) hides 100% of P" in w for w in Array.from_experiment(hidden).warnings())
    assert simulate(clear, 50_000, seed=1).select("P", counted=False).sum() > 0
    assert len(simulate(hidden, 50_000, seed=1)) == 0  # every particle aimed at P stops in the housing


def test_dead_material_in_the_beam_is_reported():
    downstream = dict(CLOVER, theta="0 deg", distance="150 mm")
    warnings = Array.from_experiment(Experiment.from_dict(setup(gammas=[downstream]))).warnings()
    assert any("Clo (housing) sits in the beam path" in w for w in warnings)


def test_monte_carlo_agrees_with_the_rates_for_an_s3():
    exp = Experiment.from_dict(setup(gammas=[]))
    events = simulate(exp, 400_000, seed=11)
    rate, error = events.rate("S3", counted=False)
    assert abs(rate - Rates(exp).rate("S3", counted=False)) < 4 * error


# -- the chamber ----------------------------------------------------------------------------------------------------


def test_chamber_keeps_particle_detectors_in_and_gamma_detectors_out():
    chamber = {"radius": "120 mm", "wall_thickness": "3 mm", "wall_material": "Al", "beam_pipe_radius": "20 mm"}
    exp = Experiment.from_dict(setup(chamber=chamber))
    assert str(exp.chamber.radius) == "120 mm" and exp.chamber.wall_material == "Al"
    assert Experiment.from_toml(exp.to_toml()).chamber == exp.chamber
    with pytest.raises(SetupError) as e:
        Experiment.from_dict(setup([dict(S3, distance="110 mm")], [dict(CLOVER, distance="125 mm")],
                                   chamber=chamber))
    text = "\n".join(e.value.problems)
    assert "detector 1 (S3): reaches 145 mm from the target, outside the chamber (radius 120 mm)" in text
    assert "gamma detector 1 (Clo): its front is 120 mm from the target, inside the chamber wall" in text
    with pytest.raises(SetupError, match="chamber: 'radius' is missing"):
        Experiment.from_dict(setup(chamber={"wall_material": "Al"}))


# -- the planner ----------------------------------------------------------------------------------------------------


def test_planner_geometry_draws_crystals_and_boards():
    p = Planner(Experiment.from_dict(setup()))
    g = p.geometry()
    assert [b["name"] for b in g["blocking"]] == ["S3 (chip edge)", "S3 (board)"]
    clover, labr = g["gamma_detectors"]
    assert (clover["model"], clover["crystals"], labr["crystals"]) == ("clover", 4, 1)
    # The clover's outline reaches the back of its 70 mm crystals.
    depth = np.nanmax(clover["outline"] @ np.array(p.experiment.gamma_detectors[0].direction()))
    assert depth == pytest.approx(270.0)
    pytest.importorskip("plotly")
    from physim.nuclear import app

    app.figure_geometry(p).to_json()
