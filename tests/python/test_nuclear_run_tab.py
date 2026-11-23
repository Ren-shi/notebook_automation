"""The Run tab's model (backlog item 66): progress in experiment time, live counters that end at the summary's
numbers, Stop and Extend, the run list, and what is not simulated."""

import threading
import time

import numpy as np
import pytest

from physim.nuclear import gamma_events
from physim.nuclear import runs as _runs
from physim.nuclear import workbench as wb
from physim.nuclear.planner import Planner


def test_progress_is_shown_throughout_and_ends_at_the_summary(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    seen = []
    run = p.start_run("beam", duration="24 h", budget="2 min", progress=seen.append)
    assert len(seen) >= 2, "progress should come more than once"
    times = [g["real_s"] for g in seen]
    assert times == sorted(times) and times[-1] == pytest.approx(120.0)
    last = seen[-1]
    assert last["counts"] == run.summary["counts"]
    assert last["coincidences"] == run.summary["coincidences"]
    assert last["scaled_s"] == pytest.approx(86400.0 - 120.0)
    assert last["real_done"]
    # Coincidences in the full-energy peak are a part of those at any energy; the projection scales them up.
    for k, v in last["peak_coincidences"].items():
        assert v <= last["coincidences"][k]
    cd = sum(v for k, v in last["peak_coincidences"].items() if k.startswith("CD|"))
    assert last["projected"]["CD"] == pytest.approx(cd * 86400.0 / 120.0)


def test_estimate_before_the_run():
    p = Planner.example("alpha_on_gold")
    e = p.estimate_run("2 h", "10 min")
    assert e["real_s"] == 600 and e["scaled_s"] == pytest.approx(6600)
    total = sum(r["rate_per_s"] for r in p.rates()["rows"])
    assert e["particles"] == pytest.approx(600 * total)
    assert e["cpu_s"] > 0 and e["size_mb"] > 0
    assert p.estimate_run("5 min", "1 h")["scaled_s"] == 0


def test_stop_keeps_what_is_accumulated_and_the_list_shows_it(tmp_path):
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    bg = p.start_run("beam", duration="1 h", budget="1 h", background=True)
    deadline = time.time() + 60
    while time.time() < deadline and (p.run_progress() or {}).get("real_s", 0) <= 0:
        time.sleep(0.2)
    assert p.run_sample(5), "the last particles are there while the run is taken"
    run = p.stop_run()
    assert bg.result is not None and run.summary["status"] == "stopped"
    longer = p.extend_run("1 min")
    assert longer.real_s == pytest.approx(run.real_s + 60.0)
    rows = p.runs()
    assert len(rows) == 1 and rows[0]["extended"] == 1


@pytest.mark.parametrize("where", ["particles", "gammas"])
def test_progress_shows_a_chunk_only_when_its_particles_are_stored(tmp_path, monkeypatch, where):
    """A step is held where it turns the generator's chunk into stored columns (its beam time is known to the
    generator by then), or in the γ chain: the progress, the counters and the sample still show the chunks before it,
    all of the same beam time, and all of the new chunk once the step is over."""
    p = Planner.example("coulex_ni58")
    folder = _runs.ExperimentFolder.create(p.experiment, root=tmp_path)
    taker = _runs.start(folder, "beam", duration="1 h", budget="1 h")
    held, go = threading.Event(), threading.Event()
    module, name = (_runs, "_store_particles") if where == "particles" else (gamma_events, "simulate_gammas")
    real = getattr(module, name)

    def wait(*args, **kw):
        out = real(*args, **kw)
        held.set()
        assert go.wait(60)
        return out

    monkeypatch.setattr(module, name, wait)
    for _ in range(2):
        shown, last = taker.progress(), taker.sample(5)
        held.clear()
        go.clear()
        t = threading.Thread(target=taker.step)
        t.start()
        assert held.wait(60)
        during, sample, ahead = taker.progress(), taker.sample(5), taker.runner.real_s
        go.set()
        t.join(60)
        assert ahead > during["real_s"], "the generator is ahead of what is stored"
        for k in ("real_s", "events", "counts", "gamma", "coincidences", "peak_coincidences"):
            assert during[k] == shown[k], k
        assert sample == last
        after = taker.progress()
        assert after["real_s"] == ahead
        assert after["events"] == sum(after["counts"].values()) > during["events"]
        assert len(taker.sample(5)) == 5
    # What is stored goes with the beam time that is stored.
    taker.stop = True
    run = taker.write()
    assert run.summary["generator"]["real_s"] == run.real_s == after["real_s"]
    assert run.summary["particles"] == after["events"]


def test_a_watched_run_has_the_events_of_its_seed(tmp_path):
    """The same seed gives the same events whether the run is taken at once or watched from another thread."""
    p = Planner.example("coulex_ni58")
    p.create_experiment(root=tmp_path)
    plain = p.start_run("beam", duration="1 min", seed=5)
    bg = p.start_run("beam", duration="1 min", seed=5, background=True)
    while bg.running:
        p.run_progress()
        p.run_sample(5)
    watched = bg.join()
    a, b = plain._arrays(), watched._arrays()
    assert set(a) == set(b) and all(np.array_equal(a[k], b[k]) for k in a)
    assert plain.summary["generator"] == watched.summary["generator"]
    assert plain.summary["counts"] == watched.summary["counts"]


def test_the_not_simulated_list_names_items_71_to_74():
    items = {n for n, _, _ in wb.NOT_SIMULATED}
    assert items == {71, 72, 73, 74}
    names = " ".join(what for _, what, _ in wb.NOT_SIMULATED).lower()
    for word in ("cascades", "summing", "pile-up", "lifetimes", "contaminant", "halo"):
        assert word in names
