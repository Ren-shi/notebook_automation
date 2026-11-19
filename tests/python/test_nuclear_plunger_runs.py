"""Runs of a plunger setup (lifetimes wave 2A): one run per target-to-stopper distance makes a distance series,
stored and reloaded through the experiment folder; and the "fusion" kind, handed to ``fusion_runs``."""

import json
import sys
import types

import numpy as np
import pytest

import physim.nuclear as nuclear
from physim.nuclear import Experiment
from physim.nuclear import runs as runs_module
from physim.nuclear.runs import KINDS, DistanceSeries, ExperimentFolder, RunData


@pytest.fixture(scope="module")
def ru98():
    return Experiment.example("coulex_ru98_plunger")


def test_a_distance_series_round_trips_through_the_folder(tmp_path, ru98):
    folder = ExperimentFolder.create(ru98, root=tmp_path)
    taken = [folder.take_run("beam", duration="1 h", budget="0.4 s", distance_um=d) for d in (10.0, 100.0)]
    # A run without a distance takes the setup's first.
    first = folder.take_run("beam", duration="1 h", budget="0.4 s")
    assert [r.distance_um for r in taken + [first]] == [10.0, 100.0, 3.0]
    assert [r.label for r in taken + [first]] == ["01-beam-1h-10um", "02-beam-1h-100um", "03-beam-1h-3um"]
    assert json.loads((taken[1].folder / "summary.json").read_text(encoding="utf-8"))["distance_um"] == 100.0
    assert taken[0].describe().startswith("Beam run, 1h at 10 µm")

    # Reloaded from disk: the distance, and where each γ ray was emitted.
    shares = []
    for r in taken:
        again = folder.load_run(r.number)
        assert again.distance_um == r.distance_um
        g = again.gammas()
        assert g.distance_um == r.distance_um and np.all(g["distance_um"] == r.distance_um)
        assert g.segments == ["target", "gap", "stopper"]
        assert set(np.unique(g["segment"])) <= {0, 1, 2} and g["stopped"].dtype == bool
        assert np.all(g["segment"][g["stopped"]] == 2)
        shares.append(float(np.mean(g["segment"] == 2)))
        with np.load(again.folder / "events.npz") as z:
            assert z["g_stopped"].dtype == np.int8 and len(z["g_segment"]) == len(z["g_particle"])
    assert shares[0] > shares[1] + 0.2  # further from the stopper, fewer decays in it

    series = folder.distance_series()
    assert len(series) == 1 and isinstance(series[0], DistanceSeries)
    s = series[0]
    assert s.distances == [3.0, 10.0, 100.0] and len(s) == 3
    assert [r.number for r in s.runs[100.0]] == [2]
    assert 3.0 not in s.missing() and 450.0 in s.missing()
    assert folder.distance_series(ru98)[0].distances == s.distances

    # A second run at one distance joins it; a run with another stopper starts a series of its own.
    folder.take_run("beam", duration="1 h", budget="0.2 s", distance_um=10.0)
    thicker = Experiment.from_dict(ru98.to_dict())
    thicker.lifetime.stopper.thickness = "20 mg/cm2"
    folder.take_run("beam", duration="1 h", budget="0.2 s", experiment=thicker, distance_um=50.0)
    series = folder.distance_series()
    assert [x.distances for x in series] == [[3.0, 10.0, 100.0], [50.0]]
    assert [r.number for r in series[0].runs[10.0]] == [1, 4]
    assert [x.distances for x in folder.distance_series(thicker)] == [[50.0]]


def test_extending_and_retaking_keep_the_distance(tmp_path, ru98):
    folder = ExperimentFolder.create(ru98, root=tmp_path)
    run = folder.take_run("beam", duration="1 h", budget="0.3 s", distance_um=25.0)
    n = len(run.gammas())
    t = runs_module.extend(run, "1 h")
    t.budget_s = 0.6
    t.real_target = 0.6
    while not t.done:
        t.step()
    longer = t.write()
    assert longer.distance_um == 25.0 and longer.number == run.number and len(longer.gammas()) > n
    assert np.all(longer.gammas()["distance_um"] == 25.0)
    again = runs_module.retake(folder, longer)
    assert again.distance_um == 25.0
    with pytest.raises(ValueError, match="negative"):
        folder.take_run("beam", duration="1 h", budget="0.1 s", distance_um=-5.0)
    assert folder.take_run("beam", duration="1 h", budget="0.1 s", distance_um="60 um").distance_um == 60.0


def test_runs_without_a_plunger_are_as_before(tmp_path):
    exp = Experiment.example("coulex_ni58")
    folder = ExperimentFolder.create(exp, root=tmp_path)
    run = folder.take_run("beam", duration="10 s")
    assert run.distance_um is None and run.summary["distance_um"] is None and run.label == "01-beam-10s"
    g = run.gammas()
    assert g.distance_um is None and np.all(np.isnan(g["distance_um"])) and not g["stopped"].any()
    assert folder.distance_series() == []
    with pytest.raises(ValueError, match="plunger"):
        folder.take_run("beam", duration="10 s", distance_um=10.0)

    # A run stored before the events said where each γ ray was emitted: no such columns. It loads, and extends.
    with np.load(run.folder / "events.npz") as z:
        arrays = {k: z[k] for k in z.files if not k.endswith(("_stopped", "_segment"))}
    np.savez_compressed(run.folder / "events.npz", **arrays)
    summary = json.loads((run.folder / "summary.json").read_text(encoding="utf-8"))
    summary.pop("distance_um")
    (run.folder / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    old = RunData.load(run.folder)
    assert old.distance_um is None and not old.gammas()["stopped"].any()
    t = runs_module.extend(old, "10 s")
    while not t.done:
        t.step()
    longer = t.write()
    with np.load(longer.folder / "events.npz") as z:
        assert len(z["g_segment"]) == len(z["g_particle"]) and len(z["g_stopped"]) == len(z["g_particle"])
    assert longer.duration_s == pytest.approx(20.0) and longer.gammas() is not None


def test_a_fusion_run_is_handed_to_fusion_runs(tmp_path, monkeypatch):
    """``runs.start(folder, kind="fusion", ...)`` imports ``fusion_runs`` (written by another session) only then and
    returns what its ``start`` returns; here a stand-in takes its place."""
    assert "fusion" in KINDS
    seen = {}

    class Taker:  # RunTaker-like: stepped until done, then written
        done = False

        def __init__(self):
            self.extra = {}

        def step(self):
            self.done = True

        def progress(self):
            return {}

        def write(self):
            return "written"

    marker = Taker()

    def start(folder, **kwargs):
        seen["folder"], seen["kwargs"] = folder, kwargs
        marker.done = False
        return marker

    stand_in = types.ModuleType("physim.nuclear.fusion_runs")
    stand_in.start = start
    monkeypatch.setitem(sys.modules, "physim.nuclear.fusion_runs", stand_in)
    monkeypatch.setattr(nuclear, "fusion_runs", stand_in, raising=False)
    folder = ExperimentFolder.create(Experiment.example("coulex_ni58"), root=tmp_path)
    out = runs_module.start(folder, kind="fusion", duration="1 h", budget="30 s", seed=4, distance_um=50.0,
                            gate="fold >= 2")
    assert out is marker and seen["folder"] is folder
    assert seen["kwargs"] == {"duration": "1 h", "budget": "30 s", "seed": 4, "experiment": None, "plan": None,
                              "distance_um": 50.0, "gate": "fold >= 2"}
    # The folder's own take_run goes the same way.
    assert folder.take_run("fusion", duration="2 h") == "written"
    assert "distance_um" not in seen["kwargs"] and seen["kwargs"]["duration"] == "2 h"
    # Options a particle run does not know are a mistake there.
    with pytest.raises(TypeError, match="gate"):
        runs_module.start(folder, kind="beam", duration="10 s", gate="fold >= 2")

    # A stored fusion run loads as fusion_runs' own RunData (which may itself call RunData.load), and extending or
    # retaking it goes to fusion_runs too.
    class FusionRunData(RunData):
        @classmethod
        def load(cls, path):
            base = RunData.load(path)
            return cls(base.folder, base.summary, base.setup)

    stand_in.FusionRunData = FusionRunData
    stand_in.extend = lambda run, more: ("extended", run.number, more)
    path = folder.path / "runs" / "01-fusion-1h-50um"
    path.mkdir(parents=True)
    folder.setup().save(path / "setup.toml")
    (path / "summary.json").write_text(json.dumps({"number": 1, "kind": "fusion", "seed": 9, "duration_s": 3600.0,
                                                   "real_s": 60.0, "distance_um": 50.0,
                                                   "trigger": {"fold": 3, "neutrons_min": 1, "charged_min": None,
                                                               "gates": True}}), encoding="utf-8")
    run = folder.load_run(1)
    assert type(run) is FusionRunData and run.distance_um == 50.0
    assert run.describe().startswith("Fusion run, 1h at 50 µm")
    assert runs_module.extend(run, "1 h") == ("extended", 1, "1 h")
    again = runs_module.retake(folder, run)
    assert again is marker and again.extra == {"retake_of": 1}
    assert seen["kwargs"]["seed"] == 9 and seen["kwargs"]["distance_um"] == 50.0
    # The retake keeps the run's trigger.
    assert [seen["kwargs"][k] for k in ("fold", "neutrons_min", "charged_min", "gates")] == [3, 1, None, True]
    assert [s.distances for s in folder.distance_series()] == [[50.0]]
