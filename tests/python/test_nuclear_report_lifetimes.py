"""The beam-time report's lifetime section (lifetime wave 3A): the plan per level, the caveats, the limits this
setup reaches, and the lifetime checks picked by capability (so checks added later appear without a change in the
report). The setups are the tests' own; the validation results are synthetic, so nothing here waits on the suite."""

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from physim.nuclear import Experiment, report  # noqa: E402
from physim.nuclear import validation as v  # noqa: E402

RDDS = "Lifetimes: recoil distance (RDDS)"
DSAM = "Lifetimes: Doppler-shift attenuation (DSAM)"


def _coulex(**lifetime):
    d = Experiment.example("coulex_ru98_plunger").to_dict()
    d["lifetime"] = {"method": "rdds", "distances": ["5 um", "20 um", "60 um", "200 um"],
                     "stopper": {"kind": "stopper", "material": "Cu", "thickness": "20 mg/cm2"}, **lifetime}
    return Experiment.from_dict(d)


@pytest.fixture(scope="module")
def rep():
    return report.build(_coulex(), events=5_000, validation=False)


def test_the_report_has_a_lifetime_section_for_a_lifetime_setup(rep):
    assert rep.lifetime is not None
    plan = rep.lifetime["plan"]
    assert plan["available"] and plan["levels"][0]["method"] == "rdds"
    parts = rep.parts(dpi=40)
    html = parts["lifetime"]
    assert "Lifetime measurement" in html and "Plan per level" in html
    assert "5, 20, 60, 200 µm" in html and "20 mg/cm2 of Cu" in html
    assert parts["lifetime_plain"].startswith("<h3>")
    page = rep.to_html()
    assert page.index("Lifetime measurement") < page.index("<h2>Energy loss</h2>")
    numbers = rep.numbers()["lifetime"]
    assert numbers["levels"][0]["tau_ps"] == pytest.approx(plan["levels"][0]["tau_ps"])
    assert len(numbers["flags"]) == len(report.LIFETIME_LIMITS)
    assert rep.lifetime["analysis"] is None  # no runs: nothing to analyse


def test_caveats_say_what_this_setup_reaches(rep):
    flags = {lim["key"]: lim["flag"] for lim in rep.lifetime["limits"]}
    # A Coulex plunger of a state longer than 1 ps without deorientation, analytic spectra in the plan.
    assert flags["dynamic_deorientation"] and "No deorientation" in flags["dynamic_deorientation"]
    assert flags["coulex_spectra"] and flags["typical_recoil"]
    assert flags["sideways"] is None  # a DSAM limit: not this plunger's
    html = rep.parts(dpi=40)["lifetime"]
    assert "Caveats" in html and "Static deorientation only" in html
    with_dor = report.lifetime_limits(_coulex(deorientation={"alpha2": 0.2, "gamma2": "0.1 /ps", "alpha4": 0.2,
                                                             "gamma4": "0.3 /ps"}), rep.lifetime["plan"])
    # The analysis corrects for deorientation where it has the alignment: the caveat no longer says it does not.
    assert "does not correct" not in with_dor[0]["flag"]
    assert "corrects the shifted and unshifted intensities" in with_dor[0]["flag"]


def test_a_setup_without_lifetime_has_no_section():
    r = report.build(Experiment.example("alpha_on_gold"), events=5_000, validation=False)
    assert r.lifetime is None and r.parts(dpi=40)["lifetime"] == "" and "lifetime" not in r.numbers()


def test_lifetime_capabilities_are_picked_by_name_and_method():
    caps = ["Rutherford cross section", RDDS, DSAM, "Lifetimes: fast timing (LaBr3)"]
    assert report.lifetime_capabilities(caps) == [RDDS, DSAM, "Lifetimes: fast timing (LaBr3)"]
    assert report.lifetime_capabilities(caps, "rdds") == [RDDS, "Lifetimes: fast timing (LaBr3)"]
    assert report.lifetime_capabilities(caps, "dsam") == [DSAM, "Lifetimes: fast timing (LaBr3)"]
    assert set(report.lifetime_capabilities(v.CAPABILITIES)) >= {RDDS, DSAM}


def test_checks_appear_by_capability_with_known_failures_marked():
    results = [v.Result(RDDS, "literature", "rdds_new_check_from_3c", "a paper", "5%", "pass", 0.01, detail="τ 8.4 ps"),
               v.Result(DSAM, "tool", "dsam_lineshape_vs_apcad", "APCAD 1.0", "10%", "fail", 0.17,
                        detail="first-order shift physim 17% larger"),
               v.Result(DSAM, "tool", "dsam_other", "AhKin", "10%", "fail", 0.3, detail="too wide"),
               v.Result("Rutherford cross section", "literature", "geiger", "x", "1%", "pass", 0.0)]
    rows = report.lifetime_checks(results, known={"dsam_lineshape_vs_apcad"})
    by = {r["name"]: r for r in rows}
    assert set(by) == {"rdds_new_check_from_3c", "dsam_lineshape_vs_apcad", "dsam_other"}
    assert by["rdds_new_check_from_3c"]["status"] == "pass"
    assert by["dsam_lineshape_vs_apcad"]["status"] == "known failure"
    assert "17% larger" in by["dsam_lineshape_vs_apcad"]["detail"]  # the suite's own numbers, not the report's
    assert by["dsam_other"]["status"] == "fail"
    assert [r["capability"] for r in rows] == [RDDS, DSAM, DSAM]
    assert {r["name"] for r in report.lifetime_checks(results, known=set(), method="rdds")} == {
        "rdds_new_check_from_3c"}


def test_known_failures_reach_the_caveats_and_the_models_table(rep):
    rep.lifetime["checks"] = report.lifetime_checks(
        [v.Result(RDDS, "tool", "rdds_known_one", "NaPyTau 1.0", "5%", "fail", 0.2, detail="physim 20% off")],
        known={"rdds_known_one"})
    parts = rep.parts(dpi=40)
    assert "rdds_known_one" in parts["lifetime"] and "physim 20% off" in parts["lifetime"]
    assert "Lifetime checks" in parts["models"] and "known failure" in parts["models"]
    rep.lifetime["checks"] = []


def test_a_dsam_report_names_its_backing():
    d = Experiment.example("fusion_sm143_dsam").to_dict()
    d["reaction"]["side_feeding_levels"] = 0
    d.setdefault("levels", {})["143Sm"] = {
        "nuclide": "143Sm",
        "level": [{"energy": "0 keV", "jpi": "13/2-"}, {"energy": "530 keV", "jpi": "15/2-", "half_life": "0.3 ps"}],
        "transition": [{"from": 1, "to": 0, "energy": "530 keV", "intensity": 100, "multipolarity": "M1"}]}
    r = report.build(Experiment.from_dict(d), events=5_000, validation=False)
    html = r.parts(dpi=40)["lifetime"]
    assert "Doppler-shift attenuation" in html and "13 mg/cm2 of Au" in html
    flags = {lim["key"]: lim["flag"] for lim in r.lifetime["limits"]}
    assert flags["stopping"] and flags["sideways"]  # DSAM, with detectors at 90°
    assert r.lifetime["plan"]["levels"][0]["attenuation"] is not None
