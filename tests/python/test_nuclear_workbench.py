"""The workbench's structure (backlog item 64): one mode, the setup as decisions with their consequences, the status
strip, new experiments and templates, the seven stages."""

import pytest

from physim.nuclear import ensdf, guide
from physim.nuclear import workbench as wb
from physim.nuclear.planner import Planner


def _card(p, key):
    return next(d for d in wb.decisions(p) if d.key == key)


def test_new_experiment_without_a_form_of_every_field(tmp_path):
    exp = wb.new_setup("16O", "64 MeV", "208Pb", measure="elastic")
    p = Planner(exp)
    p.create_experiment(root=tmp_path)
    # The detectors come after, from the Add menus.
    names = [d.get("name") for d in p.draft["detectors"]]
    assert p.add_detector(**guide.placement("forward", names))
    assert p.add_gamma_detector(**guide.gamma_placement("clover"))
    strip = {x["label"]: x["text"] for x in wb.status_strip(p)}
    assert strip["Experiment"] == "16O on 208Pb"
    assert strip["Beam"] == "16O at 64 MeV" and strip["Target"].startswith("208Pb")
    assert strip["Detectors"] == "2 particle + 1 γ"
    assert strip["Last run"] == "no run yet"
    # A light beam on a heavy target gets its first detector behind the target.
    assert exp.detectors[0].theta.to("deg") == 180
    assert wb.new_setup("208Pb", "1 GeV", "12C").detectors[0].theta.to("deg") == 45


def test_a_template_brings_only_its_detectors():
    exp = wb.new_setup("62Ni", "237 MeV", "194Pt", measure="elastic", detectors_from="coulex-cd-clovers")
    d = exp.to_dict()
    assert d["beam"]["nuclide"] == "62Ni" and d["beam"]["energy"] == "237 MeV"
    assert d["target"]["material"] == "194Pt"
    template = wb.template("coulex-cd-clovers").to_dict()
    assert [x["name"] for x in d["detectors"]] == [x["name"] for x in template["detectors"]]
    assert [x["name"] for x in d["gamma_detectors"]] == [x["name"] for x in template["gamma_detectors"]]
    assert len(wb.new_setup("62Ni", "237 MeV", "194Pt", detectors_from="blank").detectors) == 1


def test_an_elastic_run_with_gamma_detectors_has_no_gamma_rays(tmp_path):
    # A template's clovers on a setup that excites nothing: the run's data hold particles only.
    p = Planner(wb.new_setup("62Ni", "237 MeV", "208Pb", measure="elastic", detectors_from="coulex-cd-clovers"))
    p.create_experiment(root=tmp_path)
    run = p.start_run("beam", duration="1 h", budget="5 s")
    assert run.gammas() is None and run.gammas(plain=True) is None
    assert len(run.events().columns["event"])


@pytest.mark.skipif(not ensdf.available(), reason="no local copy of ENSDF (scripts/fetch_ensdf.py)")
def test_choosing_coulex_on_the_card_fills_in_the_first_2plus_state():
    p = Planner(wb.new_setup("62Ni", "237 MeV", "184Pt", measure="elastic"))
    assert p.set("reaction", "type", "coulex")
    exc = p.experiment.excitation
    assert exc.excite == "target" and exc.energy_mev == pytest.approx(0.163, abs=1e-3)
    assert p.experiment.levels["target"].nuclide == "184Pt"
    # The other nucleus: its own state, not the target's.
    assert p.set("reaction", "excite", "projectile")
    assert p.experiment.excitation.energy_mev == pytest.approx(1.1729, abs=1e-3)
    assert p.experiment.levels["beam"].nuclide == "62Ni"
    # Values the user typed are kept when switching to elastic and back.
    assert p.set("reaction", "energy", "1.2 MeV") and p.set("reaction", "type", "elastic")
    assert p.set("reaction", "type", "coulex")
    assert p.experiment.excitation.energy_mev == pytest.approx(1.2)


def test_the_pickers_offer_elements_mass_numbers_and_compounds():
    assert wb.elements()["Ni"].startswith("Ni · nickel")
    masses = wb.mass_numbers("Ni")
    assert masses[62] == "62 · 3.63 %" and 56 in masses  # stable with its abundance, and radioactive ones
    assert wb.most_abundant("Pt") == 195 and wb.most_abundant("Tc") is None
    assert wb.split_nuclide("62Ni") == ("Ni", 62) and wb.split_nuclide("Pt") == ("Pt", None)
    assert wb.split_nuclide("CD2") is None and wb.split_nuclide("") is None


@pytest.mark.skipif(not ensdf.available(), reason="no local copy of ENSDF (scripts/fetch_ensdf.py)")
def test_the_state_picker_lists_states_with_a_known_b_up():
    p = Planner(wb.new_setup("62Ni", "237 MeV", "184Pt", measure="coulex-target"))
    choices = wb.states(p.experiment.levels["target"])
    first = next(iter(choices))
    assert choices[first].startswith("2+ · 163.0 keV · B(E2↑) 39,477")
    # A new target: its own state, and the old target's scheme is gone.
    assert p.set("target", "material", "194Pt")
    assert p.experiment.levels["target"].nuclide == "194Pt"
    assert p.experiment.excitation.energy_mev == pytest.approx(0.3285, abs=1e-3)
    # The beam is not the excited nucleus: changing it leaves the target's state alone.
    assert p.set("beam", "nuclide", "58Ni")
    assert p.experiment.excitation.energy_mev == pytest.approx(0.3285, abs=1e-3)


def test_every_card_has_a_consequence_that_follows_its_fields():
    p = Planner.example("coulex_ni58")
    cards = wb.decisions(p)
    assert [c.group for c in cards][:2] == ["Beam", "Target and reaction"]
    groups = list(dict.fromkeys(c.group for c in cards))
    assert groups == list(wb.GROUPS)
    for c in cards:
        assert c.state and c.consequence and c.consequence != "—", c.key
    before = _card(p, "target").consequence
    p.set("target", "thickness", "1.0 mg/cm2")
    after = _card(p, "target").consequence
    assert before != after and "loses" in after
    cd = _card(p, "detector:0").consequence
    p.set("detector 1", "distance", "60 mm")
    assert _card(p, "detector:0").consequence != cd
    assert _card(p, "detector:0").state.startswith("CD: annular at 60 mm")


def test_every_field_of_every_card_has_help():
    for name in Planner.examples():
        p = Planner.example(name)
        for c in wb.decisions(p):
            sec = "detector" if c.section.startswith("detector") else (
                "gamma" if c.section.startswith("gamma") else c.section)
            missing = [f for f, _, _ in c.fields + c.more if guide.help_for(sec, f) is None]
            assert not missing, (name, c.key, missing)
    for f, _, _ in wb.BACKING_FIELDS:
        assert guide.help_for("backing", f) is not None


def test_checks_say_when_the_setup_changed_since_the_run(tmp_path):
    p = Planner.example("alpha_on_gold")
    p.create_experiment(root=tmp_path)
    p.start_run("beam", duration="5 s")
    assert not any("changed since run" in c["text"] for c in wb.checks(p))
    p.set("beam", "energy", "6 MeV")
    assert "changed since run 1" in wb.checks(p)[0]["text"]
    strip = {x["label"]: x["text"] for x in wb.status_strip(p)}
    assert strip["Last run"].startswith("Beam run, 5s") and "changed since" in strip["Last run"]


def test_stages_and_templates():
    assert [label for _, label in wb.STAGES] == ["Setup", "Plan", "Run", "Data", "Analysis", "Report", "Physics"]
    for key, label, example in wb.TEMPLATES:
        exp = wb.template(key)
        exp.validate()
        if example:
            assert exp.to_dict() == Planner.example(example).experiment.to_dict()
    with pytest.raises(KeyError):
        wb.template("nope")


def test_the_app_opens_on_the_beam_excitation_at_backward_angles(tmp_path, monkeypatch):
    """With no experiment yet the app opens the first template: 62Ni (6+, 237 MeV) on 194Pt, the projectile
    excited, every particle detector behind 90°."""
    from physim.nuclear.app.page import _starting_planner
    from physim.nuclear.detectors import Array

    assert wb.TEMPLATES[0][2] == "coulex_ni62_pt194"
    exp = wb.template(wb.TEMPLATES[0][0])
    assert (exp.beam.nuclide, exp.beam.charge_state, exp.beam.energy_mev) == ("62Ni", 6, 237.0)
    assert exp.target.material == "194Pt" and exp.excitation.excite == "projectile"
    assert exp.excitation.energy_mev == pytest.approx(1.17291)
    # The silicon behind the target stands 300 mm from it, and no particle detector reaches forward of 90°:
    # there the scattered beam and the slow 194Pt recoils of distant collisions arrive.
    cd = exp.detectors[0]
    assert cd.position_mm() == pytest.approx((0.0, 0.0, -300.0), abs=1e-9)
    ranges = {g.name: g.theta_range() for g in Array.from_experiment(exp)}
    assert all(lo > 100 for lo, _ in ranges.values()), ranges
    assert ranges["CD"][0] > 170
    # Both nuclei have their level scheme: 62Ni to its 4+, 194Pt with its first 2+.
    beam, target = exp.levels["beam"], exp.levels["target"]
    assert [(lv.label, round(lv.energy.value)) for lv in beam.levels[1:5:3]] == [("2+", 1173), ("4+", 2336)]
    assert (target.nuclide, target.levels[1].label, round(target.levels[1].energy.value)) == ("194Pt", "2+", 328)
    p = Planner(exp)
    assert not [w for w in p.warnings() if w.level != "note"]
    assert set(p.levels()) >= {"beam", "target"}
    # No experiment and nothing asked for: that template (an existing experiment would be opened instead).
    monkeypatch.setenv("PHYSIM_EXPERIMENTS", str(tmp_path))
    first, offer = _starting_planner(None, None, None)
    assert first.experiment.to_dict() == exp.to_dict() and offer


def test_coulex_needs_a_state_or_says_so():
    exp = wb.new_setup("16O", "64 MeV", "Xe", measure="coulex-beam")
    if exp.excitation is None:
        notes = wb.new_setup_notes(exp)
        assert notes and "B(E2)" in notes[0]
    else:
        assert exp.excitation.excite == "projectile"


def test_changes_between_the_applied_setup_and_a_draft():
    """What the app's "not applied" bar lists, and which cards it marks."""
    from physim.nuclear.planner import Planner

    p = Planner.example("coulex_ni58")
    draft = Planner(p.experiment)
    draft.adopt(p)
    assert wb.changes(p._draft, draft._draft) == []
    assert draft.set("beam", "energy", "31 MeV") and draft.set("run", "counts_wanted", 888)
    assert draft.set("detector 2", "distance", "90 mm")
    found = wb.changes(p._draft, draft._draft)
    assert found == ["Beam: energy 30 MeV → 31 MeV", "DSSD-L: distance 100 mm → 90 mm",
                     "Run conditions: counts wanted 500 → 888"]
    # Detectors are matched by name: with the first one removed, the others are not changes.
    assert draft.remove_detector(0) and draft.add_detector(**{**p._draft["detectors"][2], "name": "New"})
    found = wb.changes(p._draft, draft._draft)
    assert "CD: removed" in found and "New: added" in found and not any(x.startswith("DSSD-R") for x in found)
    assert wb.card_setup(draft._draft, "detector:0")["name"] == "DSSD-L"
    assert wb.card_setup(p._draft, "detector:0", name="DSSD-L") == p._draft["detectors"][1]
    assert wb.card_setup(draft._draft, "detector:0") != wb.card_setup(p._draft, "detector:0", name="DSSD-L")
    assert wb.card_setup(draft._draft, "detector:1") == wb.card_setup(p._draft, "detector:1", name="DSSD-R")
    assert wb.card_setup(p._draft, "detector:9") is None and wb.card_setup(p._draft, "gamma:0", name="nope") is None
    assert wb.card_setup(draft._draft, "target") == wb.card_setup(p._draft, "target")
    assert draft.set("reaction", "type", "elastic")
    assert wb.card_setup(draft._draft, "target") != wb.card_setup(p._draft, "target")
    assert any(x.startswith("Reaction: type") for x in wb.changes(p._draft, draft._draft))


def test_the_cards_of_a_draft_compute_nothing():
    from physim.nuclear.planner import Planner

    p = Planner.example("coulex_ni58")
    cards = wb.draft_decisions(p)
    assert [c.key for c in cards] == [c.key for c in wb.decisions(p)]
    assert "rates" in p._cache  # from decisions(); the draft's cards did not need them:
    q = Planner.example("coulex_ni58")
    cards = wb.draft_decisions(q)
    assert not q._cache and all(c.consequence == "—" for c in cards)
    assert cards[0].state == wb.decisions(q)[0].state and cards[0].values == q._draft["beam"]
