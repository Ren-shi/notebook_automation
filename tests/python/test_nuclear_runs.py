"""Experiments and runs (backlog item 63): unweighted events for the real part of a run, the rest scaled; run
folders that reload; the setup locked during a run; runs marked stale when the physics changes."""

import math
import time

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.planner import Planner
from physim.nuclear.rates import Rates
from physim.nuclear.runs import ExperimentFolder, RunData, RunInProgress, physics_changes, slug


def _within(measured, expected, sigmas=4.0):
    return abs(measured - expected) <= sigmas * math.sqrt(max(expected, 1.0))


def test_ten_minute_run_matches_the_rates(tmp_path):
    exp = Experiment.example("coulex_ni58")
    folder = ExperimentFolder.create(exp, root=tmp_path)
    t = time.time()
    run = folder.take_run("beam", duration="10 min")
    took = time.time() - t
    r = Rates(exp)
    for d in run.summary["detectors"]:
        assert _within(run.counts(d), r.rate(d) * 600.0), d
    assert run.real_s == pytest.approx(600.0)
    assert took < 60.0, f"took {took:.0f} s"
    assert (run.folder / "events.npz").stat().st_size < 50e6
    ev = run.events()
    assert np.all(ev["weight"] * ev.beam_time_s == pytest.approx(1.0))


def test_run_folder_reloads_the_same(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    run = p.start_run("beam", duration="20 s")
    assert run.label == "01-beam-20s"
    spectra = p.spectra()["spectra"]
    gammas = p.gamma_spectra()

    q = Planner.open_experiment(p.experiment_folder)
    assert q.run.label == run.label
    assert q.run.summary["counts"] == run.summary["counts"]
    assert q.experiment.to_dict() == p.experiment.to_dict()
    assert q.run.setup == run.setup
    for name, s in q.spectra()["spectra"].items():
        assert np.array_equal(s["counts"], spectra[name]["counts"])
    g2 = q.gamma_spectra()
    for name in gammas["spectra"]:
        assert np.array_equal(g2["spectra"][name]["measured"], gammas["spectra"][name]["measured"])
    assert q.runs()[0]["describe"].startswith("Beam run, 20s")


def test_physics_change_marks_the_run_stale_and_nothing_resimulates(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    run = p.start_run("beam", duration="10 s")
    n_gammas = len(p.gamma_events())
    assert not p.run_status()["stale"]
    p.set("title", None, "Another title")
    p.set("detector 1", "name", "Renamed")
    assert not p.run_status()["stale"]
    p.set("target", "thickness", "1.0 mg/cm2")
    status = p.run_status()
    assert status["stale"] and any("target.thickness" in c for c in status["changes"])
    assert p.runs()[0]["stale"]
    # The views still read the run, taken with the old target; nothing simulates on its own.
    assert len(p.gamma_events()) == n_gammas
    assert p.gamma_spectra()["run"]["stale"]
    assert p._data_experiment().to_dict() == run.experiment.to_dict()


def test_long_run_is_scaled_and_labelled(tmp_path):
    exp = Experiment.example("alpha_on_gold")
    folder = ExperimentFolder.create(exp, root=tmp_path)
    run = folder.take_run("beam", duration="24 h", budget="60 s")
    assert run.real_s == pytest.approx(60.0) and run.duration_s == pytest.approx(86400.0)
    assert run.summary["scaled"] and "scaled" in run.describe()
    r = Rates(exp)
    for d in run.summary["detectors"]:
        expected = r.rate(d) * 86400.0
        real = run.counts(d)
        # The full run's count is the real part's, scaled: its uncertainty is that of the real part.
        assert abs(run.counts_in_run(d) - expected) <= 4 * run.scale * math.sqrt(max(real, 1.0)), d
    h = np.array([0.0, 10.0, 100.0])
    scaled = run.scaled(h)
    assert scaled[0] == 0 and abs(scaled[2] - 100 * run.scale) < 6 * math.sqrt(100 * run.scale)


def test_setup_is_locked_while_a_run_is_taken(tmp_path):
    p = Planner.example("alpha_on_gold")
    p.create_experiment(root=tmp_path)
    bg = p.start_run("beam", duration="10 min", budget="10 min", background=True)
    assert p.running
    with pytest.raises(RunInProgress):
        p.set("beam", "energy", "6 MeV")
    run = p.stop_run()
    assert not p.running
    assert run.summary["status"] == "stopped" and run.duration_s == pytest.approx(run.real_s)
    assert run.duration_s < 600.0
    assert p.set("beam", "energy", "6 MeV")
    assert bg.result is run or bg.result.label == run.label


def test_extend_continues_the_same_stream(tmp_path):
    p = Planner.example("alpha_on_gold")
    p.create_experiment(root=tmp_path)
    first = p.start_run("beam", duration="10 s")
    n0 = first.summary["particles"]
    ev0 = first.events()["event"].copy()
    longer = p.extend_run("10 s")
    assert longer.duration_s == pytest.approx(20.0) and longer.real_s == pytest.approx(20.0)
    ev = longer.events()["event"]
    assert len(ev) > n0
    assert np.array_equal(ev[:n0], ev0)
    assert np.all(np.diff(ev) >= 0) and ev[n0] > ev0[-1]
    assert longer.summary["extended"] == 1


def test_source_and_alignment_runs_are_listed(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    src = p.start_run("source", duration="10 min", source="152Eu")
    assert src.kind == "source" and src.label == "01-152Eu-source"
    assert src.source().spectra and sum(src.summary["gamma_counts"].values()) > 0
    al = p.start_run("alignment", duration="10 s", offset_mm=1.0)
    assert al.summary["assumed_offset_mm"] == 1.0
    kinds = [r["kind"] for r in p.runs()]
    assert kinds == ["source", "alignment"]
    out = p.alignment(fit=False)
    assert out["available"] and out["offset_mm"] == 1.0


def test_analysis_reads_unit_weight_runs(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="60 s")
    g = p.gamma_events()
    assert len(g) > 0
    assert np.allclose(g["weight"] * g.events.beam_time_s, 1.0)
    out = p.analyse()
    assert "available" in out


def test_labels_are_not_physics():
    a = Experiment.example("coulex_ni58").to_dict()
    b = Experiment.example("coulex_ni58").to_dict()
    b["title"] = "x"
    b["detectors"][0]["name"] = "y"
    assert physics_changes(a, b) == []
    b["beam"]["energy"] = "70 MeV"
    assert physics_changes(a, b) == [f"beam.energy: {a['beam']['energy']} → 70 MeV"]


def test_experiment_names_and_folders(tmp_path):
    assert slug("⁵⁸Ni on ²⁰⁸Pb") == "58Ni-on-208Pb"
    exp = Experiment.example("alpha_on_gold")
    one = ExperimentFolder.create(exp, root=tmp_path)
    two = ExperimentFolder.create(exp, root=tmp_path)
    assert one.path != two.path and one.name == two.name
    listed = Planner.experiments(tmp_path)
    assert {x["path"] for x in listed} == {str(one.path), str(two.path)}
    assert isinstance(RunData, type)


def test_segment_directions_match_each_rows_own_segment():
    """The table over the hit segments gives every row the mean angles of its own particle in its own segment: the
    recoils' where a recoil hit, the scattered beam's (at the setup's masses) where that did."""
    from physim.nuclear.detectors import Array, target_mass_ratio
    from physim.nuclear.runs import _segment_directions

    exp = Planner.example("coulex_ni58").experiment
    array = Array.from_experiment(exp)
    ratio = target_mass_ratio(exp)
    assert ratio == pytest.approx(16 / 58, rel=0.01)
    rng = np.random.default_rng(3)
    det = rng.integers(0, len(array.geometries), 400)
    si = np.array([rng.integers(0, len({s[0] for s in array.geometries[d].segments})) for d in det])
    sj = np.array([rng.integers(0, len({s[1] for s in array.geometries[d].segments})) for d in det])
    recoil = rng.integers(0, 2, 400).astype(bool)
    theta, phi = _segment_directions(exp, {"detector": det, "segment_i": si, "segment_j": sj, "recoil": recoil})
    for k in range(0, 400, 7):
        th, ph, _ = array.geometries[int(det[k])].segment_angles(
            (int(si[k]), int(sj[k])), "recoil" if recoil[k] else "ejectile", ratio)
        assert theta[k] == th and phi[k] == ph
    # One segment, both particles: two different angles.
    i = [g.name for g in array.geometries].index("DSSD-L")
    one = {"detector": np.full(2, i), "segment_i": np.zeros(2), "segment_j": np.full(2, 8),
           "recoil": np.array([False, True])}
    th, _ = _segment_directions(exp, one)
    assert abs(th[1] - th[0]) > 0.01
    empty = _segment_directions(exp, {"detector": det[:0], "segment_i": si[:0], "segment_j": sj[:0],
                                      "recoil": recoil[:0]})
    assert len(empty[0]) == 0 and len(empty[1]) == 0


@pytest.mark.parametrize("name", ["coulex_ni62_c_allmond", "coulex_ni58", "coulex_ru98_plunger"])
def test_stored_angles_are_the_mean_angles_of_the_rings(name):
    """What a stored run takes for its particles' angles is, ring by ring, the mean angle of the particles that hit:
    within 0.02° (and the statistics of the simulation) for scattered beam and recoils, in normal and inverse
    kinematics, with sectors, strips and a whole ring."""
    from physim.nuclear.events import simulate
    from physim.nuclear.runs import _segment_directions

    exp = Experiment.example(name)
    c = simulate(exp, 1_000_000, 5).columns
    theta, _ = _segment_directions(exp, c)
    seen = 0
    for d in range(len(exp.detectors)):
        for kind in (False, True):
            m = c["counted"] & (c["detector"] == d) & (c["recoil"] == kind)
            for ring in np.unique(c["segment_i"][m]):
                r = m & (c["segment_i"] == ring)
                if r.sum() < 1000:
                    continue
                w, off = c["weight"][r], theta[r] - c["theta"][r]
                mean = np.average(off, weights=w)
                err = math.sqrt(np.sum((w * (off - mean)) ** 2)) / w.sum()
                assert abs(mean) < 0.02 + 4 * err, (d, kind, ring, mean, err)
                seen += 1
    assert seen >= 2


def test_a_ring_without_sectors_corrects_for_the_motion_along_the_beam_only():
    """A ring all round the beam says nothing about the azimuth: γ rays of one energy in crystals at one polar angle
    get one corrected energy, whatever their azimuth, and the particle's θ (not 0) sets the emitter's speed."""
    from physim.nuclear.detectors import Array
    from physim.nuclear.gamma_events import doppler_correct

    exp = Experiment.example("coulex_ru98_plunger")
    array = Array.from_experiment(exp)
    th, _, known = array.geometries[0].segment_angles((0, 0), "recoil")
    assert known == 0.0 and 15.0 < th < 25.0
    elements = [(i, k) for i, gd in enumerate(exp.gamma_detectors) for k, _ in enumerate(gd.elements())]
    names = [gd.name for gd in exp.gamma_detectors]
    for group in (("F1", "F2", "F3", "F4"), ("B1", "B2", "B3", "B4")):
        crystals = np.array([c for c, (i, k) in enumerate(elements) if names[i] in group and k == 0])
        n = len(crystals)
        cols = {"detector": np.zeros(n), "segment_i": np.zeros(n), "segment_j": np.zeros(n),
                "recoil": np.ones(n, dtype=bool)}
        out = doppler_correct(exp, cols, np.arange(n), crystals, np.full(n, 0.65246))["ejectile"]
        assert out == pytest.approx(np.full(n, out[0]), rel=1e-9)
        assert abs(out[0] / 0.65246 - 1) > 0.01


def test_kinematic_lines_are_at_the_angles_of_their_own_particle():
    from physim.nuclear.dataviews import kinematic_lines
    from physim.nuclear.detectors import Array, target_mass_ratio

    exp = Experiment.example("coulex_ni58")
    g = Array.from_experiment(exp)["CD"]
    ratio = target_mass_ratio(exp)
    for line in kinematic_lines(exp, "CD"):
        assert line["particle"] == "ejectile"
        assert line["theta"] == pytest.approx([g.segment_angles((r, 0), "ejectile", ratio)[0] for r in line["rings"]])
    exp = Experiment.example("coulex_ni62_c_allmond")
    g = Array.from_experiment(exp)["Ring2"]
    lines = [x for x in kinematic_lines(exp, "Ring2") if x["particle"] == "recoil"]
    assert {x["group"] for x in lines} == {"elastic", "excited"}
    for line in lines:
        assert line["theta"] == pytest.approx([g.segment_angles((r, 0), "recoil")[0] for r in line["rings"]])
