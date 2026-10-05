"""The planner's guidance (backlog item 46): help for every input, the guided steps, and the readings of every result."""

import pytest

from physim.nuclear import guide
from physim.nuclear.planner import TABS, Planner


def _app():
    pytest.importorskip("plotly")
    from physim.nuclear import app

    return app


def test_every_input_has_help():
    app = _app()
    fields = [("title", "")]
    fields += [("beam", f) for f, _, _ in app.BEAM_FIELDS]
    fields += [("target", f) for f, _, _ in app.TARGET_FIELDS]
    fields += [("backing", f) for f, _, _ in app.BACKING_FIELDS]
    fields += [("run", f) for f, _, _ in app.RUN_FIELDS]
    fields += [("reaction", f) for f, _, _ in app.REACTION_FIELDS] + [("reaction", f) for f in
                                                                      ("type", "excite", "multipolarity")]
    fields += [("detector 1", f) for f, _, _ in app.DETECTOR_FIELDS] + [("detector 3", "shape")]
    fields += [("gamma detector 2", f) for f, _, _ in app.GAMMA_FIELDS]
    missing = [f for f in fields if guide.help_for(*f) is None]
    assert not missing, missing
    for h in guide.HELP.values():
        assert h.what and h.typical and h.effect
        assert h.text().count("\n") == (1 if h.effect == "—" else 2)


@pytest.mark.parametrize("name", Planner.examples())
def test_every_tab_has_a_reading(name):
    p = Planner.example(name)
    for tab in TABS:
        lines = guide.reading(p, tab)
        if tab == "gamma" and p.experiment.excitation is None:
            assert lines == []
            continue
        assert lines and all(isinstance(s, str) and s.endswith((".", ")")) for s in lines), (tab, lines)
        assert not any("e+0" in s for s in lines), lines  # numbers are written for people


def test_readings_use_this_setups_numbers():
    p = Planner.example("alpha_on_gold")
    rows = p.rates()["rows"]
    fastest = max(rows, key=lambda r: r["rate_per_s"])["detector"]
    text = " ".join(guide.reading(p, "rates"))
    assert text.startswith(f"{fastest} counts fastest") and "pile-up" in text  # A20 is above 5000/s
    loss = p.energy_loss()["layers"][0]["loss_mev"]
    assert f"{loss * 1e3:.3g} keV" in guide.reading(p, "energy_loss")[0]
    # The Coulomb-excitation peaks in the CD overlap over the whole detector, but not ring by ring.
    c = Planner.example("coulex_ni58")
    spectra = " ".join(guide.reading(c, "spectra"))
    assert "In CD" in spectra and "overlap" in spectra and "ring by ring" in spectra
    assert "keV less than an elastic one" in " ".join(guide.reading(c, "gamma"))
    # Halving the beam time below what is needed changes the verdict.
    c.set("run", "beam_time", "1 s")
    assert "is not enough" in " ".join(guide.reading(c, "rates"))


def test_guided_steps():
    p = Planner.example("alpha_on_gold")
    keys = [s["step"].key for s in guide.steps_for(p)]
    assert keys == ["goal", "beam", "target", "detectors", "rates", "spectra", "report"]
    assert [s["step"].key for s in guide.steps_for(Planner.example("coulex_ni58"))] == [
        "goal", "beam", "target", "detectors", "gamma", "rates", "spectra", "report"]
    assert all(t in TABS for s in guide.STEPS for t in s.tabs)
    assert set(guide.GOALS) == {"elastic", "coulex"}
    assert all(ex in Planner.examples() for _, _, ex in guide.GOALS.values())

    def problems(planner):
        return {s["step"].key: s["problems"] for s in guide.steps_for(planner) if s["problems"]}

    assert problems(p) == {}
    p.set("beam", "energy", "abc")
    assert list(problems(p)) == ["beam"]
    p.set("beam", "energy", "5.5 MeV")
    p.set("detector 2", "theta", "200 deg")
    assert list(problems(p)) == ["detectors"]
    p.set("detector 2", "theta", "30 deg")
    p.set("backing", "material", "C")  # no thickness yet
    assert list(problems(p)) == ["target"]
    p.set("backing", "material", None)
    p.set("run", "beam_time", "12")
    assert list(problems(p)) == ["rates"]
    p.set("run", "beam_time", "12 h")
    p.set("reaction", "type", "coulex")  # the state is not given yet
    assert list(problems(p)) == ["goal"]
    c = Planner.example("coulex_ni58")
    c.set("gamma detector 1", "radius", "x")
    assert list(problems(c)) == ["gamma"]


def test_placements():
    p = Planner.example("alpha_on_gold")
    names = [d.name for d in p.experiment.detectors]
    for pl in guide.PLACEMENTS:
        fields = guide.placement(pl.key, names)
        assert fields["name"] not in names
        assert p.add_detector(**fields), (pl.key, p.problems)
        names.append(fields["name"])
    ranges = {g["name"]: g["theta_range"] for g in p.geometry()["detectors"]}
    assert ranges["F1"][1] < 90 and ranges["B1"][0] > 90
    lo, hi = ranges["CD1"]
    assert 120 < lo < hi < 170  # the ring sits around the beam, backward
    assert guide.placement("ring", ["CD1", "CD2"])["name"] == "CD3"


def test_coulomb_excitation_beam_time_is_set_by_coincidences():
    """Coulomb excitation is counted in particle–γ coincidences, so its beam time is hours, not seconds."""
    p = Planner.example("coulex_ni58")
    r = p.rates()
    assert r["measured"] == "coincidences" and r["gamma_efficiency_typical"]
    rates = p._rates()
    e0 = p.experiment.excitation.energy_mev
    eff = sum(g.peak_efficiency(e0, p.experiment) for g in p.experiment.gamma_detectors)
    assert r["gamma_efficiency"] == pytest.approx(eff)
    for row in r["rows"]:
        exc = rates.rate(row["detector"], what="excitations")
        assert row["excitation_per_s"] == pytest.approx(exc) and exc < 0.01 * row["rate_per_s"]
        # In coincidence with each particle detector the γ rays are not emitted evenly: its own efficiency.
        own = row["gamma_efficiency"]
        assert 0.5 * eff < own < 1.5 * eff
        assert row["coincidence_per_s"] == pytest.approx(exc * own)
        assert row["beam_time_s"] == pytest.approx(r["counts_wanted"] / (exc * own))
        assert row["beam_time_s"] > 3600
    # A measured photopeak efficiency replaces the typical response.
    for i in range(len(p.experiment.gamma_detectors)):
        assert p.set(f"gamma detector {i + 1}", "efficiency", "1 %")
    r2 = p.rates()
    assert r2["gamma_efficiency"] == pytest.approx(0.04) and not r2["gamma_efficiency_typical"]
    assert not p.set("gamma detector 1", "efficiency", "150 %") and "at most 100 %" in p.problems[0]
    text = " ".join(guide.reading(Planner.example("coulex_ni58"), "rates"))
    assert "particle–γ coincidences" in text and "typical response" in text and " h." in text
    # Without γ detectors, the excitation events themselves.
    q = Planner.example("coulex_ni58")
    for _ in range(4):
        q.remove_gamma_detector(0)
    assert q.rates()["measured"] == "excitations"
    # Elastic setups are unchanged.
    a = Planner.example("alpha_on_gold")
    assert a.rates()["measured"] == "all" and "excitation_per_s" not in a.rates()["rows"][0]


def test_geometry_reading_points_backward():
    p = Planner.example("alpha_on_gold")
    for i in range(len(p.experiment.detectors) - 1, -1, -1):
        if p.geometry()["detectors"][i]["theta_range"][1] > 90:
            p.remove_detector(i)
    assert "Backward pad" in " ".join(guide.reading(p, "geometry"))
    c = " ".join(guide.reading(Planner.example("coulex_ni58"), "geometry"))
    assert "CD has the largest share of excitation events" in c
