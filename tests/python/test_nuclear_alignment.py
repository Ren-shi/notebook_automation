"""A movable target, misalignment and sources placed anywhere (backlog 60)."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment, response
from physim.nuclear.alignment import diagnostic, fit_offset, flatness, overlay, with_offset
from physim.nuclear.analysis import Analysis, Settings
from physim.nuclear.detectors import Array
from physim.nuclear.experiment import SetupError
from physim.nuclear.gamma_events import recorrect, simulate_gammas
from physim.nuclear.planner import Planner


def moved(mm: float) -> Experiment:
    d = Experiment.example("coulex_ni58").to_dict()
    d["target"]["position"] = f"{mm} mm"
    return Experiment.from_dict(d)


def test_moving_the_target_changes_each_ring_angle_as_a_hand_calculation_gives():
    exp0, exp5 = Experiment.example("coulex_ni58"), moved(5.0)
    cd0, cd5 = Array.from_experiment(exp0).geometries[0], Array.from_experiment(exp5).geometries[0]
    det = exp0.detectors[0]
    r_in, r_out = 9.0, 41.0
    n_rings = det.rings
    for k in range(n_rings):
        lo, hi = r_in + k * (r_out - r_in) / n_rings, r_in + (k + 1) * (r_out - r_in) / n_rings
        # The CD stands 30 mm behind the chamber's centre; the target moved 5 mm forward makes it 35 mm.
        for g, dist in ((cd0, 30.0), (cd5, 35.0)):
            a, b = g.theta_range((k, 0))
            assert b == pytest.approx(180 - math.degrees(math.atan(lo / dist)), abs=0.01)
            assert a == pytest.approx(180 - math.degrees(math.atan(hi / dist)), abs=0.01)
    assert exp5.detectors[0].position_mm() == pytest.approx((0.0, 0.0, -35.0), abs=1e-9)
    assert exp5.detectors[0].chamber_position_mm() == pytest.approx((0.0, 0.0, -30.0), abs=1e-9)
    # γ-ray detectors follow too: the one at 90° and 120 mm is now at 120.1 mm and a little backward.
    g90 = exp5.gamma_detectors[0]
    assert g90.distance_mm() == pytest.approx(math.hypot(120.0, 5.0))
    assert g90.direction()[2] == pytest.approx(-5.0 / math.hypot(120.0, 5.0))
    # The setup file keeps the chamber's coordinates and round-trips.
    again = Experiment.from_toml(exp5.to_toml())
    assert again.target.position_mm == 5.0 and again.detectors[0].position_mm() == exp5.detectors[0].position_mm()
    assert "_origin" not in exp5.to_toml()


def test_a_setup_with_the_target_at_the_origin_is_unchanged():
    p0 = Planner.example("coulex_ni58")
    p1 = Planner(moved(0.0))
    assert p0.rates()["rows"] == p1.rates()["rows"]
    d = Experiment.example("coulex_ni58").to_dict()
    assert "position" not in d["target"]


def test_a_target_ladder_picks_the_target_in_the_beam():
    d = Experiment.example("coulex_ni58").to_dict()
    d["target"] = {"ladder": [["58Ni", "0.5 mg/cm2"], ["208Pb", "1 mg/cm2"]], "selected": 2}
    exp = Experiment.from_dict(d)
    assert exp.target.material == "208Pb" and str(exp.target.thickness) == "1 mg/cm2"
    again = Experiment.from_toml(exp.to_toml())
    assert again.target.material == "208Pb" and again.target.ladder[0][0] == "58Ni"
    d["target"]["selected"] = 3
    with pytest.raises(SetupError, match="the ladder has 2 targets"):
        Experiment.from_dict(d)
    d["target"] = {"thickness": "1 mg/cm2"}
    with pytest.raises(SetupError, match="material"):
        Experiment.from_dict(d)


def test_a_source_moved_off_centre_changes_each_crystal_rate_with_its_solid_angle():
    exp = Experiment.example("coulex_ni58")
    gd = exp.gamma_detectors[0]                       # Ge90: 35 mm radius at 120 mm along +y
    on, off = response.Response(exp, gd), response.Response(exp, gd, source=(0.0, 60.0, 0.0))
    assert on.coverage[0] == pytest.approx((1 - math.cos(math.atan(35 / 120))) / 2)
    assert off.coverage[0] == pytest.approx((1 - math.cos(math.atan(35 / 60))) / 2, rel=1e-6)
    away = response.Response(exp, gd, source=(0.0, -60.0, 0.0))
    assert away.coverage[0] == pytest.approx((1 - math.cos(math.atan(35 / 180))) / 2, rel=1e-6)
    sideways = response.Response(exp, gd, source=(60.0, 0.0, 0.0))
    assert sideways.coverage[0] < on.coverage[0]
    run0 = response.source_run(exp, "60Co", "10 kBq", "10 min")
    run1 = response.source_run(exp, "60Co", "10 kBq", "10 min", position=(0.0, 60.0, 0.0))
    for name in ("Ge90", "Ge45"):
        g = next(x for x in exp.gamma_detectors if x.name == name)
        r0, r1 = response.Response(exp, g), response.Response(exp, g, source=(0.0, 60.0, 0.0))
        assert run1.expected[name].sum() / run0.expected[name].sum() == pytest.approx(
            r1.coverage[0] / r0.coverage[0], rel=1e-9)


@pytest.fixture(scope="module")
def ni():
    exp = Experiment.example("coulex_ni58")
    return exp, simulate_gammas(exp, 600_000, seed=1)


def test_the_correction_with_the_true_geometry_is_what_the_events_carry(ni):
    exp, g = ni
    again = recorrect(g, exp)
    assert np.allclose(again["corrected_recoil"], g["corrected_recoil"])
    wrong = recorrect(g, with_offset(exp, 3.0))
    assert not np.allclose(wrong["corrected_recoil"], g["corrected_recoil"])
    # The target moved 3 mm forward: a detector at 30 mm behind is now 33 mm away.
    assert with_offset(exp, 3.0).detectors[0].position_mm()[2] == pytest.approx(-33.0)
    assert with_offset(with_offset(exp, 3.0), -3.0).target.position is None


def test_the_diagnostic_is_flat_when_aligned_and_the_fit_recovers_an_offset(ni):
    exp, g = ni
    rows = diagnostic(g, exp)
    assert len(rows) > 100 and all(r["error_kev"] > 0 and r["n"] >= 3 for r in rows)
    # Against the pattern a simulation with the right geometry gives, the plot is flat within statistics.
    reference = diagnostic(simulate_gammas(exp, 600_000, seed=7), exp)
    chi2 = flatness(rows, reference=reference)
    assert chi2 < 1.5 * len(rows)
    aligned = fit_offset(g, exp)
    assert abs(aligned["offset_mm"]) < 2.5 * aligned["uncertainty_mm"] and not aligned["at_edge"]
    # The analysis assumes the target 2 mm from where it is: the fit gives the 2 mm back, within its uncertainty.
    found = fit_offset(g, with_offset(exp, 2.0))
    assert abs(found["offset_mm"] - (-2.0)) < 1.5 * found["uncertainty_mm"], found
    assert found["against"] == "simulation" and len(found["scan"]) == 17
    # Judged against flat lines instead, the finite crystals leave a pattern of their own.
    assert flatness(rows) > chi2


def test_the_overlay_shows_the_shift_and_broadening(ni):
    exp, g = ni
    o = overlay(g, exp, with_offset(exp, 6.0))
    assert len(o["true"]) == len(o["assumed"]) == len(o["edges"]) - 1
    assert o["true"].sum() == pytest.approx(o["assumed"].sum(), rel=0.05)
    assert abs(o["shift_kev"]) < 5 and o["fwhm_assumed_kev"] > 0 and "broadening_kev" in o


def test_the_analysis_can_assume_a_misplaced_target_and_budgets_it(ni):
    exp, g = ni
    a = Analysis(exp, g)
    right = a.run(Settings(particle_gate="all"))
    off = a.run(Settings(particle_gate="all", target_offset_mm=6.0))
    assert off.settings["target_offset_mm"] == 6.0
    d_mm = np.mean([gd.distance_mm() for gd in exp.gamma_detectors])
    assert right.budget["detector positions"] >= 2 * 1.0 / d_mm - 1e-9
    assert abs(off.b_e2fm4 / right.b_e2fm4 - 1) < 0.1


def test_the_planner_and_the_figures():
    pytest.importorskip("plotly")
    from physim.nuclear import app

    p = Planner.example("coulex_ni58")
    r = p.alignment(2.0, 150_000, 1, fit=False)
    assert r["available"] and r["offset_mm"] == 2.0 and "fit" not in r and r["diagnostic"]
    app.figure_alignment(p, r).to_json()
    app.figure_overlay(r).to_json()
    assert not Planner.example("alpha_on_gold").alignment()["available"]
