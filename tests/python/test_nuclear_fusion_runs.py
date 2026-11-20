"""Simulated fusion-evaporation runs (lifetimes wave 2E): residue events and their γ rays, triggered on γ-ray fold
and the ancillaries' gates, stored like particle runs, started, stopped, extended and loaded like a Coulex run,
with γγ coincidences that show the band."""

import json
import math

import numpy as np
import pytest

from physim.nuclear import Experiment, fusion_runs, runs
from physim.nuclear.residues import simulate_residue_gammas, simulate_residues
from physim.nuclear.runs import ExperimentFolder, RunData

EXAMPLE = "fusion_dy158"


def _within(measured, expected, sigmas=5.0):
    return abs(measured - expected) <= sigmas * math.sqrt(max(expected, 1.0))


def _take(taker):
    while not taker.done:
        taker.step()
    return taker.write()


def _ground_band(exp, residue="158Dy"):
    """{(J_i, J_f): energy keV} of the yrast E2 transitions of the residue's ground band."""
    scheme = exp.levels[residue]
    out = {}
    for t in scheme.transitions:
        ji, jf = scheme.levels[t.initial].spin, scheme.levels[t.final].spin
        if ji is not None and jf is not None and ji - jf == 2 and t.energy is not None:
            out.setdefault((ji, jf), float(t.energy.value))
    return out


def _counts_at(counts, edges, energy_kev, half=2.5):
    centres = 0.5 * (edges[1:] + edges[:-1])
    return float(counts[np.abs(centres - energy_kev) <= half].sum())


@pytest.fixture(scope="module")
def taken(tmp_path_factory):
    """One fusion_dy158 run of an hour, the first 8 s simulated event by event (γ-ray fold ≥ 2)."""
    folder = ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path_factory.mktemp("fusion"))
    run = _take(fusion_runs.start(folder, duration="1 h", budget="8 s", seed=1))
    return folder, run


def test_run_is_stored_like_a_particle_run(taken):
    folder, run = taken
    assert isinstance(run, fusion_runs.FusionRunData)
    assert run.label == "01-fusion-1h"
    for name in ("setup.toml", "events.npz", "summary.json"):
        assert (run.folder / name).is_file()
    s = json.loads((run.folder / "summary.json").read_text(encoding="utf-8"))
    assert s["kind"] == "fusion" and s["status"] == "complete"
    assert s["distance_um"] is None
    assert s["trigger"]["fold"] == fusion_runs.DEFAULT_FOLD
    assert run.real_s == pytest.approx(8.0)
    assert run.duration_s == pytest.approx(3600.0)
    assert run.scale == pytest.approx(450.0)
    assert run.describe().startswith("Fusion run, 1h, γ-ray fold ≥ 2")
    assert run.experiment.to_dict() == folder.setup().to_dict()
    # The residues are Poisson about σ_total × luminosity × real time.
    rate = simulate_residues(run.experiment, 1, 1).expected()
    assert _within(s["residues"], rate * run.real_s)
    # The trigger kept exactly the residues of fold ≥ 2, and each kept event has its γ rays.
    assert s["triggered"] == sum(s["folds"][2:]) == len(run.events()["event"]) > 0
    ev = run.events()
    assert np.all(ev["fold"] >= 2)
    assert np.all(np.diff(ev["event"]) > 0)
    assert np.all(ev["weight"] * run.real_s == pytest.approx(1.0))
    g = run.gammas()
    per_event = np.bincount(g["row"], minlength=len(ev["event"]))
    assert np.all(per_event >= ev["fold"])
    assert g["channel"].max() < len(run.channels)
    assert run.counts() == s["triggered"]


def test_loads_like_a_coulex_run(taken):
    folder, run = taken
    generic = folder.load_run(run.number)
    loaded = fusion_runs.load(generic)
    assert isinstance(loaded, fusion_runs.FusionRunData)
    assert loaded.summary == run.summary
    for k, v in run.events().items():
        assert np.array_equal(loaded.events()[k], v), k
    assert np.array_equal(loaded.gammas()["corrected"], run.gammas()["corrected"])
    assert fusion_runs.load(run.folder).label == run.label
    assert folder.runs()[0]["kind"] == "fusion"


def test_planner_particle_views_handle_a_fusion_run(tmp_path):
    """Planner.spectra, run_yields and compare with a fusion run current (no particle events): the particle
    spectra are the setup's simulated ones with a note, the multi-step yields say why there are none, and compare
    lists what a fusion run measures instead of comparing the monitor's rate with 0 (lifetime wave 4B)."""
    from physim.nuclear.planner import Planner

    p = Planner(Experiment.example(EXAMPLE))
    p.create_experiment(root=tmp_path)
    run = p.start_run(kind="fusion", duration="1 h", budget="1 s", seed=1)
    assert run.kind == "fusion" and p.run is run
    s = p.spectra(events=20_000)
    assert set(s["spectra"]) == {d.name for d in p.experiment.detectors}
    assert s["run"] is None and "fusion run" in s["note"]
    with pytest.raises(ValueError, match="fusion run"):
        p.run_yields()
    c = p.compare()
    quantities = {r["quantity"] for r in c["rows"]}
    assert "Rate" not in quantities  # the plan's particle rates do not apply
    residues = next(r for r in c["rows"] if r["quantity"] == "Residues made")
    assert residues["agree"] and residues["measured"] > 0
    line = next(r for r in c["rows"] if r["quantity"] == "Full-energy counts in the line")
    assert line["agree"] is None and line["measured"] > 0 and line["predicted"] > 0
    assert any("no particle detectors" in n for n in c["notes"])


def test_gated_matrix_shows_the_band(taken):
    _, run = taken
    band = _ground_band(run.experiment)
    gate = band[(4.0, 2.0)]
    others = [band[(2.0, 0.0)], band[(6.0, 4.0)], band[(8.0, 6.0)]]
    counts, edges = run.gated(gate, 2.5, bins=2000, range_kev=(0.0, 1000.0))
    singles, s_edges = run.gammas().spectrum(bins=2000, range_kev=(0.0, 1000.0))
    at_gate = _counts_at(counts, edges, gate)
    peaks = [_counts_at(counts, edges, e) for e in others]
    # Each coincident band member stands out of the gated spectrum's background...
    background = np.median(counts[(edges[:-1] > 500) & (edges[:-1] < 1000)]) * 5  # ±2.5 keV in 0.5 keV bins
    assert min(peaks) >= 8 and min(peaks) > 4 * (background + 1)
    # ...and the gate line itself is suppressed: a cascade emits it once (what is left there comes from other
    # residues' lines and Compton events near the gate; over seeds 1–3 it is 5–25 % of the band peaks).
    assert at_gate < 0.4 * np.mean(peaks)
    assert _counts_at(singles, s_edges, gate) > 5 * at_gate + 50
    # The matrix is symmetric and holds every pair both ways.
    m, m_edges = run.matrix(bins=500, range_kev=(0.0, 1000.0))
    assert np.array_equal(m, m.T)
    a, b = run.gammas().pairs()
    assert m.sum() == np.count_nonzero((a < 1000.0) & (b < 1000.0)) > 0


def test_stop_ends_the_run_early(tmp_path):
    folder = ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path)
    taker = fusion_runs.start(folder, duration="1 h", budget="60 s", seed=3)
    taker.step()
    taker.step()
    p = taker.progress()
    assert p["kind"] == "fusion" and 0 < p["fraction"] < 1 and p["residues"] > 0
    taker.stop = True
    assert taker.done
    run = taker.write()
    assert run.summary["status"] == "stopped"
    assert run.duration_s == pytest.approx(run.real_s) and run.real_s < 60.0
    assert run.describe().endswith("(stopped)")


def test_extend_continues_the_stream(tmp_path):
    folder = ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path)
    first = _take(fusion_runs.start(folder, duration="1.5 s", budget="10 s", seed=5))
    before = {k: v.copy() for k, v in first.events().items()}
    gen = first.summary["generated"]
    taker = fusion_runs.extend(first, "1.5 s")
    assert taker.folder == first.folder and taker.number == first.number
    run = _take(taker)
    assert run.summary["extended"] == 1
    assert run.real_s == pytest.approx(3.0) and run.duration_s == pytest.approx(3.0)
    assert run.summary["generated"] > gen
    ev = run.events()
    n = len(before["event"])
    assert len(ev["event"]) >= n
    for k, v in before.items():
        if k != "weight":  # 1 / real seconds: the real part grew
            assert np.array_equal(ev[k][:n], v), k
    assert np.all(np.diff(ev["event"]) > 0)
    assert len(folder.runs()) == 1
    rate = simulate_residues(run.experiment, 1, 1).expected()
    assert _within(run.summary["residues"], rate * 3.0)
    with pytest.raises(ValueError, match="not a fusion run"):
        fusion_runs.extend(RunData(run.folder, dict(run.summary, kind="beam"), run.setup), "1 s")


def test_trigger_fold_changes_the_rate_as_expected(tmp_path):
    exp = Experiment.example(EXAMPLE)
    folder = ExperimentFolder.create(exp, root=tmp_path)
    by_fold = {f: _take(fusion_runs.start(folder, duration="2 s", budget="2 s", seed=11, fold=f)) for f in (1, 2, 3)}
    one, two, three = by_fold[1], by_fold[2], by_fold[3]
    # The same seed makes the same residues: a higher fold keeps exactly the events of the lower one's that pass.
    assert one.summary["residues"] == two.summary["residues"] == three.summary["residues"]
    e1, e2 = one.events(), two.events()
    assert np.array_equal(e1["event"][e1["fold"] >= 2], e2["event"])
    assert three.summary["triggered"] == int(np.count_nonzero(e1["fold"] >= 3))
    assert one.summary["triggered"] > 10 * two.summary["triggered"] > 0
    assert two.summary["triggered"] > three.summary["triggered"]
    # Fold ≥ 1 at the rate the residues' γ rays reach the crystals: P(fold ≥ 1) from an independent sample.
    n = 200_000
    g = simulate_residue_gammas(exp, seed=99, residue_events=simulate_residues(exp, n, 99), gammas_per_event=1)
    c = np.asarray(g["counted"], dtype=bool)
    p1 = len(np.unique(np.asarray(g["event"])[c])) / n
    expected = one.summary["residues"] * p1
    assert _within(one.summary["triggered"], expected)
    assert one.summary["rates_per_s"]["trigger"] == pytest.approx(one.summary["triggered"] / one.real_s)


def _with_ancillaries(**sections):
    d = Experiment.example(EXAMPLE).to_dict()
    d.update(sections)
    return Experiment.from_dict(d)


def test_ancillaries_gate_and_are_recorded(tmp_path):
    exp = _with_ancillaries(
        recoil_separator={"theta_acceptance": "1.5 deg", "kind": "gas-filled"},
        gamma_calorimeter={"elements": 60, "radius": "120 mm", "material": "BGO", "fold_min": 8},
        neutron_wall={"cells": 50, "distance": "510 mm", "theta_min": "4 deg"},
        charged_particle_ball={"elements": 40, "radius": "60 mm", "material": "CsI", "threshold": "2 MeV"})
    folder = ExperimentFolder.create(exp, root=tmp_path)
    gated = _take(fusion_runs.start(folder, duration="1 s", budget="1 s", seed=21, fold=1))
    free = _take(fusion_runs.start(folder, duration="1 s", budget="1 s", seed=21, fold=1, gates=False))
    tr = gated.summary["trigger"]
    assert tr["separator"] and tr["calorimeter_fold"] == 8 and tr["gates"]
    ev, fr = gated.events(), free.events()
    assert np.all(ev["separator"]) and np.all(ev["theta"] <= 1.5 + 1e-4)
    assert np.all(ev["calorimeter"] >= 8)
    assert 0 < len(ev["event"]) < len(fr["event"])
    assert not np.all(fr["separator"]) and fr["calorimeter"].min() < 8
    # Neutrons and charged particles are recorded; they gate only when asked for.
    assert fr["neutrons"].max() >= 1 and fr["charged"].max() >= 1
    assert gated.summary["ancillaries"]["separator_accepted"] < gated.summary["residues"]
    n_gate = _take(fusion_runs.start(folder, duration="1 s", budget="1 s", seed=21, fold=1, gates=False,
                                     neutrons_min=1))
    assert np.all(n_gate.events()["neutrons"] >= 1)
    assert np.array_equal(n_gate.events()["event"], fr["event"][fr["neutrons"] >= 1])
    with pytest.raises(ValueError):
        fusion_runs.start(ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path / "b"),
                          duration="1 s", neutrons_min=1)


def test_runs_within_half_the_free_ram(tmp_path, monkeypatch):
    folder = ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path)
    monkeypatch.setattr(fusion_runs, "_memory_limit", lambda: 20_000 * fusion_runs.BYTES_PER_EVENT)
    taker = fusion_runs.start(folder, duration="1 s", budget="1 s", seed=2)
    assert taker.runner.block == 20_000 and taker.memory_limit == 20_000 * fusion_runs.BYTES_PER_EVENT
    monkeypatch.setattr(fusion_runs, "_memory_limit",
                        lambda: fusion_runs.MIN_BLOCK_EVENTS * fusion_runs.BYTES_PER_EVENT - 1)
    taker = fusion_runs.start(folder, duration="1 s", budget="1 s", seed=2)
    taker.step()
    assert taker.done and taker.stopped_for == "memory"
    run = taker.write()
    assert run.summary["status"] == "stopped" and run.summary["stopped_for"] == "memory"
    assert any("RAM" in n for n in run.summary["notes"])


def test_distance_is_stored_and_passed_when_set(tmp_path, monkeypatch):
    import functools

    from physim.nuclear import residues

    # The 50Cr plunger with the bottom of its yrast band typed in from ENSDF (NDS 157, 1 (2019)), so the test needs
    # no local ENSDF copy: 2+ 783.31 keV (T½ 9.08 ps), 4+ 1881.42 keV (T½ 2.2 ps).
    d = Experiment.example("fusion_cr50_plunger").to_dict()
    d.setdefault("levels", {})["50Cr"] = {
        "nuclide": "50Cr",
        "level": [{"energy": "0 keV", "jpi": "0+", "source": "ensdf"},
                  {"energy": "783.31 keV", "jpi": "2+", "half_life": "9.08 ps", "source": "ensdf"},
                  {"energy": "1881.42 keV", "jpi": "4+", "half_life": "2.2 ps", "source": "ensdf"}],
        "transition": [{"from": 1, "to": 0, "energy": "783.3 keV", "intensity": 100, "multipolarity": "E2"},
                       {"from": 2, "to": 1, "energy": "1098.1 keV", "intensity": 100, "multipolarity": "E2"}]}
    d["reaction"]["channels"] = [ch for ch in d["reaction"]["channels"] if ch["residue"] == "50Cr"]
    plunger = Experiment.from_dict(d)
    distance = plunger.lifetime.distances[1].to("um")              # 42 um, one of the setup's own
    made = []                                                      # the γ rays the run's generator made
    original = residues.simulate_residue_gammas

    @functools.wraps(original)                                     # keeps the signature fusion_runs inspects
    def spy(*args, **kwargs):
        g = original(*args, **kwargs)
        made.append(g)
        return g

    monkeypatch.setattr(residues, "simulate_residue_gammas", spy)
    folder = ExperimentFolder.create(plunger, root=tmp_path)
    run = _take(fusion_runs.start(folder, duration="0.02 s", budget="0.02 s", distance_um=distance))
    assert run.summary["distance_um"] == distance and run.distance_um == distance
    assert json.loads((run.folder / "summary.json").read_text(encoding="utf-8"))["distance_um"] == distance
    assert run.label.endswith(f"-{distance:g}um")
    assert made and all(np.all(g["distance_um"] == distance) for g in made)
    stored = fusion_runs.load(run.folder).gammas()
    assert len(stored["stopped"]) > 0 and stored["stopped"].any()   # some γ rays from residues in the stopper

    # A distance on a setup without a plunger is refused.
    monkeypatch.setattr(residues, "simulate_residue_gammas", original)
    folder = ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path)
    with pytest.raises(ValueError, match="plunger"):
        _take(fusion_runs.start(folder, duration="0.02 s", budget="0.02 s", distance_um=50.0))


def test_needs_a_fusion_setup_with_gamma_detectors(tmp_path):
    coulex = ExperimentFolder.create(Experiment.example("coulex_ni58"), root=tmp_path)
    with pytest.raises(ValueError, match="fusion"):
        fusion_runs.start(coulex, duration="1 s")
    folder = ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path)
    with pytest.raises(ValueError):
        fusion_runs.start(folder, duration="1 s", fold=0)
    est = fusion_runs.estimate(folder.setup(), duration="1 h", budget="10 s")
    assert est["real_s"] == 10.0 and est["scaled_s"] == 3590.0 and est["residues"] > 0


@pytest.mark.skipif("fusion" not in runs.KINDS, reason="runs.start has no fusion kind yet (wave 2A)")
def test_runs_start_takes_a_fusion_run(tmp_path):
    folder = ExperimentFolder.create(Experiment.example(EXAMPLE), root=tmp_path)
    taker = runs.start(folder, kind="fusion", duration="1 h", budget="0.5 s", seed=1)
    run = _take(taker)
    assert run.kind == "fusion" and run.summary["trigger"]["fold"] == fusion_runs.DEFAULT_FOLD
    assert run.real_s == pytest.approx(0.5)
    loaded = fusion_runs.load(folder.load_run(run.number))
    assert loaded.summary["triggered"] == run.summary["triggered"]
