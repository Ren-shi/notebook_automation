"""Particle and γ-ray tracks in the scene (backlog 59)."""

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.detectors import Array
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.planner import Planner
from physim.nuclear.tracks import describe, sample_tracks


@pytest.fixture(scope="module")
def ni():
    exp = Experiment.example("coulex_ni58")
    return exp, simulate_gammas(exp, 150_000, seed=2)


def test_each_track_ends_in_the_segment_the_event_names(ni):
    exp, g = ni
    array = Array.from_experiment(exp)
    tracks = sample_tracks(exp, g, n=40, select="all", seed=3)
    assert len(tracks) == 40
    checked = 0
    for t in tracks:
        for p in t.paths:
            h = p["hit"]
            if h and "segment" in h:
                end = np.array(p["points"][1])
                hit = array.geometries[h["index"]].hit(end / np.linalg.norm(end))
                assert hit is not None and list(hit.segment) == h["segment"] and hit.detector == h["detector"]
                assert hit.distance == pytest.approx(np.linalg.norm(end), rel=1e-9)
                checked += 1
            if h and "crystal" in h:
                end = np.array(p["points"][1])
                name, centre, radius = next(c for c in _crystals(exp) if c[0] == h["crystal"])
                assert np.linalg.norm(end - centre) <= radius + 1e-6
                checked += 1
        assert t.paths[0]["what"] == "beam" and t.paths[0]["points"][1] == [0.0, 0.0, 0.0]
        assert {p["what"] for p in t.paths} >= {"beam", "ejectile", "recoil"}
    assert checked > 40


def _crystals(exp):
    out = []
    for i, gd in enumerate(exp.gamma_detectors):
        for label, centre, radius in gd.elements():
            out.append(((gd.name or f"G{i + 1}") + (f" {label}" if label else ""), np.array(centre), radius))
    return out


def test_with_the_coincidence_filter_every_track_has_a_particle_and_a_gamma_hit(ni):
    exp, g = ni
    tracks = sample_tracks(exp, g, n=25, select="coincidences", seed=1)
    assert len(tracks) == 25
    for t in tracks:
        assert t.coincidence and t.hits and t.crystal_hits and t.gammas
        assert any(p["what"] == "gamma" and p["hit"] for p in t.paths)
        # The γ ray belongs to this event and the crystal named is the one hit.
        assert all(int(row["event"]) == t.event for row in t.gammas)
    plain = sample_tracks(exp, g, n=25, select="all", seed=1)
    assert not all(t.coincidence for t in plain)


def test_a_channel_filter_and_the_sampling_choice(ni):
    exp, g = ni
    excited = next(ch for ch in g.events.channels if "excited" in ch)
    tracks = sample_tracks(exp, g, n=15, select=excited, seed=4)
    assert all(t.channel == excited for t in tracks)
    with pytest.raises(ValueError, match="no channel"):
        sample_tracks(exp, g, n=5, select="nothing")
    # Weighted by rate the sample is dominated by the frequent (elastic, forward) events; as generated it is not.
    weighted = sample_tracks(exp, g, n=60, weighted=True, seed=5)
    plain = sample_tracks(exp, g, n=60, weighted=False, seed=5)
    assert np.mean([t.weight for t in weighted]) > 1.5 * np.mean([t.weight for t in plain])
    text = describe(weighted, "all", True, g.events.n_events)
    assert "60 of 150000" in text and "looks like a run" in text
    assert "over-represented" in describe(plain, "all", False, g.events.n_events)
    assert describe([], "coincidences", True, 1) == "No event matches the filter."
    # The same seed gives the same sample.
    assert [t.event for t in sample_tracks(exp, g, n=10, seed=9)] == [t.event for t in sample_tracks(exp, g, n=10, seed=9)]


def test_tracks_from_particle_events_alone():
    exp = Experiment.example("alpha_on_gold")
    p = Planner(exp)
    r = p.tracks(n=12, select="all", weighted=False, seed=1, events=50_000)
    assert len(r["tracks"]) == 12 and not r["gammas"] and r["channels"] == list(r["channels"])
    for t in r["tracks"]:
        assert not t.coincidence and t.hits
        assert all(p["what"] != "gamma" for p in t.paths)
    with pytest.raises(ValueError, match="give the"):
        sample_tracks(exp)


def test_the_scene_view_names_tracks_for_selection():
    pytest.importorskip("nicegui")
    from physim.nuclear import scene_view

    assert scene_view.parse_name("track:3") == (None, None)
    assert "__ITEMS__" in scene_view._ANIMATION and "requestAnimationFrame" in scene_view._ANIMATION
    for theme in scene_view.PALETTES.values():
        assert {"ejectile", "recoil", "gamma", "track_beam", "hit"} <= set(theme)
