"""The reproduced Coulomb-excitation experiment: Allmond et al., PRC 90, 034309 (2014), ⁶²Ni on carbon
(physim.nuclear.coulex_reproduction). The check itself (``reproduction_allmond2014``) runs in
test_nuclear_validation.py; these test its parts, that every number but two is inside its band, and that those two
(GOSIA as run, 1 %) miss it by the known, explained amount.

Check 4 rebuilds the paper's per-ring ratios from its matrix elements through GOSIA's ratio(M): the paper gives no
yields, so it tests whether physim extracts the same B(E2) from the same data as GOSIA did."""

import math
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear import coulex_reproduction as cr

REFERENCE = Path(__file__).resolve().parents[2] / "tests" / "reference" / "nuclear"
needs_reference = pytest.mark.skipif(not (REFERENCE / cr.REFERENCE_FILE).is_file(),
                                     reason="reference files are not shipped with the wheel")


@pytest.fixture(scope="module")
def refd():
    if not (REFERENCE / cr.REFERENCE_FILE).is_file():
        pytest.skip("no reference file")
    return cr.reference(REFERENCE)


def test_the_example_is_the_papers_setup():
    exp = Experiment.example(cr.EXAMPLE)
    assert exp.beam.energy_mev == pytest.approx(111.6)
    assert exp.target.material == "12C"
    rings = {d.name: d for d in exp.detectors}
    assert set(rings) == {"Ring2", "Ring3"}
    for name, (lo, hi) in {"Ring2": (14, 28), "Ring3": (28, 44)}.items():
        d = rings[name]
        dist = float(str(d.distance).split()[0])
        inner = float(str(d.inner_radius).split()[0])
        outer = float(str(d.outer_radius).split()[0])
        assert math.degrees(math.atan(inner / dist)) == pytest.approx(lo, abs=0.01)
        assert math.degrees(math.atan(outer / dist)) == pytest.approx(hi, abs=0.01)
    thetas = sorted(float(str(g.theta).split()[0]) for g in exp.gamma_detectors)
    assert thetas == [90.0] * 5 + [132.0] * 3 + [154.0]
    assert "Allmond" in exp.description and "assumed" in exp.description.lower()


def test_the_ratio_form_inverts():
    a, b = 1.1479762535634853e-2, -8.85933123935079e-4
    for m in (0.26, 0.29734, 0.3038, 0.34):
        assert cr.me_of(a, b, cr.ratio_of(a, b, m)) == pytest.approx(m, rel=1e-12)
    pts = {m: cr.ratio_of(a, b, m) for m in (0.28, 0.29734, 0.32)}
    assert cr.ratio_fit(pts) == pytest.approx((a, b), rel=1e-8)


@needs_reference
def test_the_reference_holds_the_papers_numbers_and_the_bands(refd):
    p = refd["paper"]
    assert (p["ring2"]["me_eb"], p["ring2"]["me_stat_eb"]) == (0.3038, 0.0037)
    assert (p["ring3"]["me_eb"], p["ring3"]["me_stat_eb"]) == (0.3002, 0.0024)
    assert (p["final"]["b_e2b2"], p["final"]["b_unc_e2b2"]) == (0.0906, 0.0037)
    assert refd["bands"]["paper"]["rel"] == pytest.approx(0.0037 / 0.0906, abs=1e-3)
    # The plan's 1 % for GOSIA as run (it was 2 % for a time, written after the answer was known); no band is raised.
    bands = refd["bands"]
    assert (bands["gosia"]["matched_rel"], bands["gosia"]["as_run_rel"]) == (0.005, 0.01)
    assert (bands["closure"]["share_of_paper_stat"], bands["closure"]["mc_sigmas"]) == (0.25, 3)
    assert (bands["stored"]["probability_rel"], bands["sensitivity"]["pp"]) == (0.001, 0.5)
    assert len(set(bands["closure"]["seeds"])) >= 4 and 11 not in bands["closure"]["seeds"]
    # A quarter of the paper's statistical error of B(E2) per ring, and of the weighted mean of both (0.0012).
    assert [100 * cr._closure_band(refd, r) for r in cr.ROWS] == pytest.approx([0.609, 0.400, 0.335], abs=1e-3)
    # GOSIA's fit goes through GOSIA's own run at 0.29734 e b.
    for ring in ("ring2", "ring3"):
        f = refd["gosia"]["fit"][ring]
        assert cr.ratio_of(f["a"], f["b"], 0.29734) == pytest.approx(refd["gosia"]["as_run"][ring]["ratio_elastic"],
                                                                     rel=1e-4)


@needs_reference
def test_the_thickness_reproduces_the_measured_losses(refd):
    c = cr.check_thickness(refd)
    assert c["passed"], c
    assert c["thickness_mg_cm2"] == pytest.approx(1.03741, rel=1e-4)
    loss = {x["ion"]: x["loss_mev"] for x in c["rows"]}
    assert loss["62Ni"] == pytest.approx(42.4, abs=0.1)


@needs_reference
def test_the_comparisons_catch_what_they_should(refd):
    """Each comparison fails when physim is off: 1 % on the ratio against GOSIA matched, 6 % on the B(E2) against
    the paper, a wrong sign of the reorientation effect."""
    good = {r: refd["gosia"]["matched"][cr.REF_RING[r]]["ratio_elastic"] for r in cr.RINGS}

    def rows(ratios):
        return {(x["setting"], x["ring"]): x for x in cr.check_gosia(refd, ratios)["rows"]}

    assert rows(good)["matched", "Ring2"]["passed"]
    off = rows({r: 1.01 * v for r, v in good.items()})
    assert not any(x["passed"] for x in off.values())
    # A miss of the recorded size is the known one; half a percent more, or an agreement, is not; nor is a miss of the
    # other band, or one on the other side.
    as_run = {r: refd["gosia"]["as_run"][cr.REF_RING[r]]["ratio_elastic"] for r in cr.RINGS}
    size = refd["bands"]["gosia"]["known_miss"]
    for extra, known, passed in ((0.0, True, False), (0.0004, True, False), (0.005, False, False),
                                 (-0.0006, False, False), (-0.005, False, True), (-0.026, False, False)):
        got = rows({r: (1 + size[cr.REF_RING[r]] + extra) * v for r, v in as_run.items()})
        for r in cr.RINGS:
            x = got["as_run", r]
            assert (x["known_miss"], x["passed"]) == (known, passed), (extra, r)
            assert not got["matched", r]["known_miss"]
    # physim's curve equal to GOSIA's: B(E2) is the paper's exactly; 6 % lower everywhere fails the 4.1 %.
    fits = refd["gosia"]["fit"]
    scan = {m: {r: cr.ratio_of(fits[cr.REF_RING[r]]["a"], fits[cr.REF_RING[r]]["b"], m) for r in cr.RINGS}
            for m in (0.28, 0.29734, 0.32)}
    same = cr.check_paper(refd, scan)
    assert same["passed"] and all(abs(x["deviation"]) < 1e-9 for x in same["rows"])
    high = {m: {r: v * 1.06 for r, v in row.items()} for m, row in scan.items()}
    assert not cr.check_paper(refd, high)["passed"]
    s = refd["gosia"]["sensitivity"]
    base = {r: 1.0 for r in cr.RINGS}
    right = cr.check_sensitivities(refd, base, {r: 1 + s[cr.REF_RING[r]]["q_plus"] for r in cr.RINGS},
                                   {r: 1 + s[cr.REF_RING[r]]["q_minus"] for r in cr.RINGS},
                                   {r: s[cr.REF_RING[r]]["feeding"] for r in cr.RINGS})
    assert right["passed"]
    flipped = cr.check_sensitivities(refd, base, {r: 1 + s[cr.REF_RING[r]]["q_minus"] for r in cr.RINGS},
                                     {r: 1 + s[cr.REF_RING[r]]["q_plus"] for r in cr.RINGS},
                                     {r: s[cr.REF_RING[r]]["feeding"] for r in cr.RINGS})
    assert not flipped["passed"]


def checks_with(rows) -> dict:
    """Checks as :func:`cr.outcome` takes them, the GOSIA rows given as (deviation, band, known)."""
    gosia = [{"setting": "as_run", "ring": r, "deviation": dev, "band": band, "passed": abs(dev) <= band,
              "known_miss": known} for r, (dev, band, known) in zip(cr.RINGS, rows)]
    other = {"rows": [{"passed": True}], "passed": True, "worst": 0.6}
    worst = max(abs(x["deviation"]) / x["band"] for x in gosia)
    return {"gosia": {"rows": gosia, "passed": all(x["passed"] for x in gosia), "worst": worst}, "closure": other}


def test_a_known_miss_is_no_pass():
    known = cr.outcome(checks_with([(0.0132, 0.01, True), (0.0116, 0.01, True)]))
    assert not known["passed"] and known["only_known_misses"] and len(known["known_misses"]) == 2
    assert known["worst"] == pytest.approx(1.32) and known["worst_other"] == pytest.approx(0.6)
    text = cr.verdict(known)
    assert text.startswith("NOT PASSED") and "+1.32 % (Ring2)" in text and "worst 1.32" in text and "0.60" in text
    new = cr.outcome(checks_with([(0.0132, 0.01, True), (0.02, 0.01, False)]))
    assert not new["passed"] and not new["only_known_misses"]
    assert cr.verdict(dict(new, checks=checks_with([(0.0132, 0.01, True), (0.02, 0.01, False)]))).startswith("FAILED")
    good = cr.outcome(checks_with([(0.002, 0.01, False), (0.001, 0.01, False)]))
    assert good["passed"] and not good["only_known_misses"] and "every number inside" in cr.verdict(good)


@pytest.fixture(scope="module")
def reproduced():
    """The full reproduction (about 110 s on one core, of which the closure's four samples take 70)."""
    r = cr.run(REFERENCE)
    if r is None:
        pytest.skip("no reference file")
    return r


def test_the_reproduction_misses_one_band_by_the_known_amount(reproduced):
    """The band of 1 % against GOSIA as run is missed, +1.32 % and +1.16 %, and only that: a change of either, or a
    number that leaves its band or an agreement that is new, fails here. The check is reported as not passed."""
    r = reproduced
    assert not r["passed"] and r["only_known_misses"], r["detail"]
    assert set(r["checks"]) == {"thickness", "gosia", "closure", "paper", "sensitivities", "stored"}
    outside = [(name, x.get("setting"), x.get("ring")) for name, c in r["checks"].items() for x in c["rows"]
               if not x["passed"]]
    assert outside == [("gosia", "as_run", "Ring2"), ("gosia", "as_run", "Ring3")], r["detail"]
    miss = {x["ring"]: x["deviation"] for x in r["known_misses"]}
    assert miss["Ring2"] == pytest.approx(0.0132, abs=5e-4) and miss["Ring3"] == pytest.approx(0.0116, abs=5e-4)
    assert r["worst"] == pytest.approx(1.32, abs=0.05) and r["worst_other"] <= 1.0
    # What GOSIA's masses and E1 polarisation explain of it, and what is left: the matched rows.
    g = {(x["setting"], x["ring"]): x for x in r["checks"]["gosia"]["rows"]}
    for ring, left in (("Ring2", 0.0020), ("Ring3", 0.0012)):
        x = g["as_run", ring]
        assert (1 + x["deviation"]) / (1 + x["explained"]) - 1 == pytest.approx(g["matched", ring]["deviation"])
        assert g["matched", ring]["deviation"] == pytest.approx(left, abs=3e-4)
    table = cr.table(r)
    assert table.count("MISS (known, explained)") == 2 and "FAIL" not in table and "The known miss: GOSIA" in table
    assert "NOT PASSED" in r["detail"] and "worst 1.32 of its band" in r["detail"]


def test_all_orders_are_judged_and_first_order_is_the_second_line(reproduced):
    c = reproduced["checks"]
    assert {x["order"] for x in c["gosia"]["rows"]} == {"all"}
    assert {x["order"] for x in c["gosia"]["second_line"]} == {"first"}
    # First order is above all orders in the ratio (it has no saturation), so it reads B(E2) lower.
    for a, b, more in zip(c["gosia"]["rows"], c["gosia"]["second_line"], (0.0070, 0.0040) * 2):
        assert b["physim"] / a["physim"] - 1 == pytest.approx(more, abs=3e-4)
    for a, b in zip(c["paper"]["rows"], c["paper"]["second_line"]["rows"]):
        assert b["b_physim_e2b2"] < a["b_physim_e2b2"]
        assert abs(b["physim_b"]) < 0.02 * abs(a["physim_b"])  # in first order the ratio goes as M²
    for x, more in zip(c["closure"]["rows"], (-0.0070, -0.0040, -0.0048)):
        assert x["first_order"] - x["deviation"] == pytest.approx(more, abs=6e-4)
    paper = {x["ring"]: x for x in c["paper"]["rows"]}
    # GOSIA as run is ~1.3 % above physim (integer masses, E1 polarisation), so physim reads B(E2) ~1.2 % lower.
    for ring in cr.RINGS:
        assert -0.025 < paper[ring]["deviation"] < 0.0
    table = cr.table(reproduced)
    assert table.count("(1st order)") == 13 and "reported" in table


def test_the_other_parts_of_the_reproduction(reproduced):
    c = reproduced["checks"]
    # A positive Q_s raises the ratio at these backward CM angles, in physim as in GOSIA.
    sens = {(x["what"], x["ring"]): x for x in c["sensitivities"]["rows"]}
    assert all(sens["q_plus", ring]["physim_percent"] > 0 > sens["q_minus", ring]["physim_percent"] for ring in cr.RINGS)
    closure = c["closure"]
    assert len(closure["seeds"]) == 4 and all(len(x["each"]) == 4 and x["scatter"] > 0 for x in closure["rows"])
    assert all(x["monte_carlo"] < 0.008 for x in closure["rows"])
    stored = c["stored"]["rows"]
    assert [x["rings"] for x in stored if x["what"] == "B(E2)"] == ["Ring2", "Ring3", "both"]
    assert all(x["error"] < 0.02 for x in stored) and all(x["counts"] > 8000 for x in stored if x["what"] == "B(E2)")
    assert np.isfinite(reproduced["worst"])


# -- a run taken in the app's way ---------------------------------------------------------------------------------------


def test_b_e2_comes_back_from_a_run_taken_in_the_app(tmp_path):
    """A run of the app's normal kind (24 h, the first 10 min event by event), stored, loaded and analysed, ring by
    ring: segments for angles, the Doppler correction and the segment averages of the stored columns. B(E2) is 50 times
    the real one (in first order, as the events are, it cancels) so that ten minutes give 500 and 1200 counts. Band and
    seed were fixed before it was run: a quarter of the paper's statistical error and three statistical errors.

    Eight runs of an hour of each kind, taken off line (seeds 8001–8008, 8101–8108, 8201–8208), gave for ring 2, ring 3
    and both +3.1(47), +5.2(30), +4.8(26) % with the real B(E2); −0.3(7), +0.3(4), 0.0(4) % with 50 times it; and
    −0.4(7), −0.4(4), −0.5(4) % with that and one ring per detector, as the paper's has."""
    from physim.nuclear.analysis import Analysis, Settings
    from physim.nuclear.runs import ExperimentFolder

    d = Experiment.example(cr.EXAMPLE).to_dict()
    truth = 50 * Experiment.example(cr.EXAMPLE).excitation.b_up_e2fm
    d["reaction"]["b_up"] = f"{truth!r} e2fm4"
    lab = ExperimentFolder.create(Experiment.from_dict(d), root=tmp_path)
    lab.take_run("beam", duration="24 h", budget=600, seed=9002)
    run = ExperimentFolder(lab.path).load_run(1)
    assert run.real_s == pytest.approx(600) and run.scale == pytest.approx(144)
    a = Analysis(run.experiment, run.gammas(), run=run)
    for row, (dets, bias) in {"Ring2": (["Ring2"], 0.0061), "Ring3": (["Ring3"], 0.0040), "both": (None, 0.0033)}.items():
        r = a.run(Settings(detectors=dets, particle_gate="all"))
        assert r.fit.failed is None and r.statistical < 0.055, row
        assert abs(r.b_e2fm4 / truth - 1) < bias + 3 * r.statistical, (row, r.b_e2fm4 / truth - 1, r.statistical)
