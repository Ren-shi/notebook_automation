"""The lifetime section of a setup (plunger RDDS and DSAM): reading, writing, problems and the planner's draft
methods (wave 1A of the lifetime work)."""

import copy

import pytest

from physim.nuclear import Experiment, SetupError
from physim.nuclear.experiment import Degrader2, Deorientation, Lifetime, Stopper
from physim.nuclear.planner import Planner

LIFETIME_EXAMPLES = ("coulex_ru98_plunger", "fusion_cr50_plunger", "fusion_sm143_dsam")


def setup(name):
    return Experiment.example(name).to_dict()


def problems_of(data):
    with pytest.raises(SetupError) as info:
        Experiment.from_dict(data)
    return info.value.problems


def one_problem(data, *words):
    """The setup has exactly one problem, and it says all of ``words``."""
    problems = problems_of(data)
    assert len(problems) == 1, problems
    assert all(w in problems[0] for w in words), problems[0]
    return problems[0]


# -- the examples ------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("name", LIFETIME_EXAMPLES)
def test_examples_load_validate_and_round_trip(name, tmp_path):
    exp = Experiment.example(name)
    exp.validate()
    assert exp.lifetime is not None
    assert Experiment.from_toml(exp.to_toml()) == exp
    assert Experiment.from_dict(exp.to_dict()) == exp
    exp.save(tmp_path / "again.toml")
    assert (tmp_path / "again.toml").read_text(encoding="utf-8") == exp.to_toml()


def test_ru98_plunger_is_the_published_setup():
    exp = Experiment.example("coulex_ru98_plunger")
    lt = exp.lifetime
    assert exp.reaction == "coulex" and exp.excitation.excite == "projectile"
    assert exp.beam.energy_mev == 300 and exp.target.material == "24Mg"
    assert lt.method == "rdds" and isinstance(lt.stopper, Stopper) and lt.degrader2 is None
    assert (lt.stopper.kind, lt.stopper.material, str(lt.stopper.thickness)) == ("stopper", "Cu", "15.7 mg/cm2")
    assert lt.distances_um[0] == 3 and lt.distances_um[-1] == 450 and lt.offset_um == 8
    assert lt.flight_distances_um[0] == 11
    assert lt.deorientation.alpha == {2: 0.18, 4: 0.24}
    assert lt.deorientation.gamma_per_ps == pytest.approx({2: 0.142, 4: 0.314})
    assert [str(g.theta) for g in exp.gamma_detectors].count("41.5 deg") == 4
    # The recoil at the paper's speed (v/c = 4.50(6) % at the shifted peak's centroid), within the model's 5 %.
    nuclide, made, left = exp.lifetime_recoil()
    beta = (2 * left / (98 * 931.494)) ** 0.5
    assert nuclide == "98Ru" and made > left > 0 and beta == pytest.approx(0.045, rel=0.05)


def test_cr50_is_a_three_foil_plunger():
    lt = Experiment.example("fusion_cr50_plunger").lifetime
    assert lt.stopper.kind == "differential" and isinstance(lt.degrader2, Degrader2)
    assert (lt.degrader2.material, str(lt.degrader2.thickness), lt.degrader2.gap_um) == ("Ta", "6.5 mg/cm2", 100)
    assert lt.offset is None and lt.offset_um == 0


def test_sm143_is_dsam_on_gold():
    exp = Experiment.example("fusion_sm143_dsam")
    assert exp.lifetime == Lifetime(method="dsam")
    assert exp.target.backing.material == "Au"
    assert exp.lifetime_recoil()[0] == "143Sm" and exp.lifetime_recoil()[2] == 0  # stops in the backing


# -- problems ----------------------------------------------------------------------------------------------------


def test_lifetime_needs_coulex_or_fusion():
    d = setup("alpha_on_gold")
    d["lifetime"] = {"method": "dsam"}
    d["target"]["backing"] = {"material": "C", "thickness": "20 mg/cm2"}
    one_problem(d, "lifetime", "\"coulex\" or \"fusion-evaporation\"")


def test_method_must_be_known():
    d = setup("fusion_sm143_dsam")
    d["lifetime"]["method"] = "doppler"
    one_problem(d, "lifetime: method must be 'rdds' or 'dsam'")
    del d["lifetime"]["method"]
    one_problem(d, "lifetime: 'method' is missing")


@pytest.mark.parametrize("distances", [None, []])
def test_a_plunger_needs_distances(distances):
    d = setup("coulex_ru98_plunger")
    if distances is None:
        del d["lifetime"]["distances"]
    else:
        d["lifetime"]["distances"] = distances
    one_problem(d, "needs its target-to-stopper distances")


@pytest.mark.parametrize("bad, words", [
    ("0 um", ("must all be positive", "number 2 is '0 um'")),
    ("-5 um", ("must all be positive", "'-5 um'")),
    ("5 MeV", ("not a length",)),
    (5, ("with a unit",)),
])
def test_distances_must_be_positive_lengths(bad, words):
    d = setup("coulex_ru98_plunger")
    d["lifetime"]["distances"] = ["3 um", bad, "10 um"]
    one_problem(d, "lifetime: distances", *words)


def test_offset_must_not_be_negative():
    d = setup("coulex_ru98_plunger")
    d["lifetime"]["offset"] = "-1 um"
    one_problem(d, "lifetime: offset must not be negative")


def test_a_plunger_needs_a_stopper():
    d = setup("coulex_ru98_plunger")
    del d["lifetime"]["stopper"]
    one_problem(d, "needs the foil after the gap", "[lifetime.stopper]")
    d = setup("coulex_ru98_plunger")
    d["lifetime"]["stopper"]["kind"] = "wall"
    one_problem(d, "lifetime.stopper: kind must be")
    d = setup("coulex_ru98_plunger")
    del d["lifetime"]["stopper"]["thickness"]
    one_problem(d, "lifetime.stopper: 'thickness' is missing")


def test_the_second_foil_goes_with_a_differential_plunger():
    d = setup("fusion_cr50_plunger")
    del d["lifetime"]["degrader2"]
    one_problem(d, "differential plunger needs its second foil", "[lifetime.degrader2]")
    d = setup("fusion_cr50_plunger")
    d["lifetime"]["stopper"]["kind"] = "degrader"
    one_problem(d, "lifetime.degrader2: a second foil only applies", "differential")
    d = setup("fusion_cr50_plunger")
    del d["lifetime"]["degrader2"]["gap"]
    one_problem(d, "lifetime.degrader2: 'gap' is missing")


def test_dsam_has_no_plunger_fields():
    d = setup("fusion_sm143_dsam")
    d["lifetime"]["distances"] = ["10 um"]
    d["lifetime"]["stopper"] = {"kind": "stopper", "material": "Au", "thickness": "10 mg/cm2"}
    one_problem(d, "distances, stopper only apply to method = \"rdds\"")


def test_dsam_needs_something_to_stop_the_recoil():
    d = setup("fusion_sm143_dsam")
    del d["target"]["backing"]
    one_problem(d, "DSAM measurement needs the recoil to slow down and stop", "[target.backing]")
    d["target"]["thickness"] = "20 mg/cm2"  # a thick target stops it by itself
    Experiment.from_dict(d)
    d = setup("fusion_sm143_dsam")
    d["target"]["backing"]["thickness"] = "0.5 mg/cm2"
    one_problem(d, "the backing (0.5 mg/cm2 of Au) does not stop the recoil", "143Sm")


def test_a_stopper_too_thin_to_stop_the_recoil():
    d = setup("coulex_ru98_plunger")
    d["lifetime"]["stopper"]["thickness"] = "2 mg/cm2"
    one_problem(d, "lifetime.stopper: 2 mg/cm2 of Cu does not stop the recoil", "98Ru", "make it thicker")


def test_a_degrader_must_let_the_recoil_through():
    d = setup("coulex_ru98_plunger")
    d["lifetime"]["stopper"]["kind"] = "degrader"
    one_problem(d, "the degrader foil (15.7 mg/cm2 of Cu) stops the recoil", "kind = \"stopper\"")
    d["lifetime"]["stopper"]["thickness"] = "1 mg/cm2"
    Experiment.from_dict(d)


def test_both_foils_of_a_differential_plunger():
    d = setup("fusion_cr50_plunger")
    d["lifetime"]["stopper"]["thickness"] = "20 mg/cm2"
    one_problem(d, "the differential foil (20 mg/cm2 of Ta) stops the recoil")
    d = setup("fusion_cr50_plunger")
    d["lifetime"]["degrader2"]["thickness"] = "1 mg/cm2"
    one_problem(d, "lifetime.degrader2: 1 mg/cm2 of Ta does not stop the recoil", "50Cr", "after the first foil")


def test_a_recoil_that_never_leaves_the_target():
    d = setup("fusion_sm143_dsam")
    d["lifetime"] = {"method": "rdds", "distances": ["10 um"],
                     "stopper": {"kind": "stopper", "material": "Au", "thickness": "10 mg/cm2"}}
    one_problem(d, "stops in the target and its backing and never reaches the plunger's gap")


def test_deorientation_is_one_form_or_the_other():
    d = setup("coulex_ru98_plunger")
    deor = d["lifetime"]["deorientation"]
    d["lifetime"]["deorientation"] = {"g": 0.47, "c2": "3.3 ps", "c4": "1.5 ps", "alpha2": 0.18, "alpha4": 0.24}
    exp = Experiment.from_dict(d)
    assert exp.lifetime.deorientation.gamma_per_ps == pytest.approx({2: 0.47 / 3.3, 4: 0.47 / 1.5})
    assert Experiment.from_toml(exp.to_toml()) == exp
    d["lifetime"]["deorientation"] = dict(deor, g=0.47)
    one_problem(d, "lifetime.deorientation: gamma2, gamma4, g mix the two forms")
    d["lifetime"]["deorientation"] = {"alpha2": 0.18, "gamma2": "0.142 /ps"}
    one_problem(d, "lifetime.deorientation: alpha4, gamma4 missing")
    d["lifetime"]["deorientation"] = {"alpha2": 0.18, "alpha4": 0.2, "g": 0.4, "c2": "3 ps"}
    one_problem(d, "lifetime.deorientation: c4 missing")
    d["lifetime"]["deorientation"] = dict(deor, gamma2="0.142 ps")
    one_problem(d, "lifetime.deorientation: gamma2", "wrong unit")
    d["lifetime"]["deorientation"] = dict(deor, alpha2=1.5)
    one_problem(d, "lifetime.deorientation: alpha2 must be between 0 and 1")
    d["lifetime"]["deorientation"] = dict(deor, beta2=0.1)
    one_problem(d, "lifetime.deorientation: unknown field 'beta2'")


def test_deorientation_only_with_a_plunger():
    d = setup("fusion_sm143_dsam")
    d["lifetime"]["deorientation"] = {"alpha2": 0.18, "gamma2": "0.142 /ps", "alpha4": 0.24, "gamma4": "0.314 /ps"}
    one_problem(d, "deorientation only apply to method = \"rdds\"")


def test_edits_in_python_are_validated():
    exp = Experiment.example("coulex_ru98_plunger")
    exp.lifetime.distances = ["10 um", "20 um"]
    exp.lifetime.deorientation = Deorientation(alpha2=0.1, alpha4=0.2, gamma2="0.1 /ps", gamma4="0.2 /ps")
    exp.validate()
    exp.lifetime.stopper.thickness = "1 mg/cm2"
    with pytest.raises(SetupError, match="does not stop the recoil"):
        exp.validate()


# -- the planner's draft methods ---------------------------------------------------------------------------------


def test_add_a_plunger_to_a_fusion_setup():
    p = Planner.example("fusion_dy158")
    assert p.add_ancillary("plunger"), p.problems
    lt = p.experiment.lifetime
    assert lt.method == "rdds" and lt.distances_um == [10, 20, 50, 100, 200, 500]
    assert (lt.stopper.kind, lt.stopper.material) == ("stopper", "Au")
    assert p.remove_lifetime() and p.experiment.lifetime is None
    with pytest.raises(KeyError, match="plunger"):
        p.add_ancillary("nope")


def test_the_default_stopper_stops_the_recoil():
    p = Planner.example("coulex_ru98_plunger")
    assert p.remove_lifetime() and p.add_ancillary("plunger"), p.problems
    nuclide, _, left = p.experiment.lifetime_recoil()
    from physim.nuclear.stopping import Stopping

    assert Stopping(nuclide, "Au").range(left) < p.experiment.lifetime.stopper.material_data().areal_density_mg_cm2(
        p.experiment.lifetime.stopper.thickness)


def test_set_lifetime_makes_and_edits_the_section():
    p = Planner.example("coulex_ni62_pt194")
    assert p.experiment.lifetime is None
    assert p.set_lifetime(distances=[5, 15.5, "1 mm"], offset=2), p.problems
    lt = p.experiment.lifetime
    assert [str(x) for x in lt.distances] == ["5 um", "15.5 um", "1 mm"] and str(lt.offset) == "2 um"
    deor = {"alpha2": 0.2, "gamma2": "0.1 /ps", "alpha4": 0.3, "gamma4": "0.2 /ps"}
    assert p.set_lifetime(deorientation=deor) and p.experiment.lifetime.deorientation.alpha[4] == 0.3
    assert p.set_lifetime(deorientation={}) and p.experiment.lifetime.deorientation is None
    assert not p.set_lifetime(distances=["0 um"]) and "must all be positive" in p.problems[0]
    assert p.set_lifetime(distances=["10 um"]) and not p.problems


def test_dsam_and_back_keeps_the_plunger():
    p = Planner.example("coulex_ru98_plunger")
    plunger = copy.deepcopy(p.experiment.lifetime)
    assert not p.set_lifetime("dsam")  # no backing, and the thin target does not stop the recoil
    assert "needs the recoil to slow down and stop" in p.problems[0]
    assert set(p.draft["lifetime"]) == {"method"}
    assert p.set("backing", "material", "Au") is False  # a backing needs its thickness too
    assert p.set("backing", "thickness", "20 mg/cm2"), p.problems
    assert p.experiment.lifetime == Lifetime(method="dsam")
    # Back to the plunger: the backing now stops the recoil before the gap, until it goes.
    assert not p.set_lifetime("rdds") and "never reaches the plunger's gap" in p.problems[0]
    p.set("backing", "thickness", None)
    assert p.set("backing", "material", None), p.problems
    assert p.experiment.lifetime == plunger


def test_stoppers_and_the_second_foil():
    p = Planner.example("fusion_cr50_plunger")
    foil2 = copy.deepcopy(p.experiment.lifetime.degrader2)
    assert p.set_stopper("stopper", "Ta", "10 mg/cm2"), p.problems
    assert p.experiment.lifetime.degrader2 is None and p.experiment.lifetime.stopper.kind == "stopper"
    assert p.set_stopper("differential", "Ta", 4.4), p.problems  # the second foil comes back
    assert p.experiment.lifetime.degrader2 == foil2
    assert str(p.experiment.lifetime.stopper.thickness) == "4.4 mg/cm2"
    assert p.set_degrader2("Au", "8 mg/cm2", 150), p.problems
    assert p.experiment.lifetime.degrader2.gap_um == 150 and p.experiment.lifetime.degrader2.material == "Au"
    assert p.set_stopper("degrader", "Ta", "2 mg/cm2") and p.experiment.lifetime.degrader2 is None
    assert p.set_degrader2("Ta", "6.5 mg/cm2", "100 um") and p.experiment.lifetime.stopper.kind == "differential"


def test_set_stopper_with_a_density():
    p = Planner.example("fusion_dy158")
    assert p.set_stopper("stopper", "Au", "5 um", density="19.3 g/cm3"), p.problems
    lt = p.experiment.lifetime
    assert lt.method == "rdds" and str(lt.stopper.density) == "19.3 g/cm3" and lt.distances


def test_the_lifetime_waits_while_the_reaction_is_elastic():
    p = Planner.example("coulex_ru98_plunger")
    plunger = copy.deepcopy(p.experiment.lifetime)
    assert p.set("reaction", "type", "elastic"), p.problems
    assert p.experiment.lifetime is None
    assert p.set("reaction", "type", "coulex"), p.problems
    assert p.experiment.lifetime == plunger
    q = Planner.example("fusion_dy158")  # no studied channel in [run], which only a fusion setup may have
    assert q.add_ancillary("plunger")
    plunger = copy.deepcopy(q.experiment.lifetime)
    assert q.set("reaction", "type", "elastic") and q.experiment.lifetime is None
    assert q.set("reaction", "type", "fusion-evaporation"), q.problems
    assert q.experiment.lifetime == plunger


def test_a_draft_is_applied_or_discarded():
    p = Planner.example("fusion_dy158")
    draft = Planner(p.experiment)
    assert draft.add_ancillary("plunger") and p.experiment.lifetime is None
    draft.adopt(p)  # discard: the draft takes the planner's setup back
    assert draft.experiment.lifetime is None
    assert draft.set_lifetime(distances=["30 um"]) and p.adopt(draft)
    assert p.experiment.lifetime.distances_um == [30]


def test_a_sweep_leaves_the_lifetime_section_out():
    # At 4 MeV the residue cannot leave the target, but the monitor's rate is still worth sweeping.
    p = Planner.example("fusion_cr50_plunger")
    s = p.sweep("beam energy", ["4 MeV", "100 MeV"])
    assert len(s["y"]) == 2 and p.experiment.lifetime is not None
