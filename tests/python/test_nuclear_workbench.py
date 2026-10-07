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


def test_coulex_needs_a_state_or_says_so():
    exp = wb.new_setup("16O", "64 MeV", "Xe", measure="coulex-beam")
    if exp.excitation is None:
        notes = wb.new_setup_notes(exp)
        assert notes and "B(E2)" in notes[0]
    else:
        assert exp.excitation.excite == "projectile"
