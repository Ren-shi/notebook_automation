"""The reproduced Coulomb-excitation experiment: Allmond et al., PRC 90, 034309 (2014), ⁶²Ni on carbon
(physim.nuclear.coulex_reproduction). The check itself (``reproduction_allmond2014``) runs in
test_nuclear_validation.py; these test its parts and that the five checks pass.

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
    assert (refd["bands"]["gosia"]["matched_rel"], refd["bands"]["gosia"]["as_run_rel"]) == (0.005, 0.02)
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
    assert cr.check_gosia(refd, good)["rows"][0]["passed"]
    assert not cr.check_gosia(refd, {r: 1.01 * v for r, v in good.items()})["passed"]
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


@needs_reference
def test_the_reproduction_passes_its_five_checks():
    """The full reproduction (about a minute: the closure run of 500 000 events dominates)."""
    r = cr.run(REFERENCE)
    assert r["passed"], r["detail"]
    assert set(r["checks"]) == {"thickness", "gosia", "closure", "paper", "sensitivities"}
    paper = {x["ring"]: x for x in r["checks"]["paper"]["rows"]}
    # GOSIA as run is ~1.3 % above physim (integer masses, E1 polarisation), so physim reads B(E2) ~1.2 % lower.
    for ring in cr.RINGS:
        assert -0.025 < paper[ring]["deviation"] < 0.0
    # A positive Q_s raises the ratio at these backward CM angles, in physim as in GOSIA.
    sens = {(x["what"], x["ring"]): x for x in r["checks"]["sensitivities"]["rows"]}
    assert all(sens["q_plus", ring]["physim_percent"] > 0 > sens["q_minus", ring]["physim_percent"] for ring in cr.RINGS)
    assert "paper" in cr.table(r) or "B(E2)" in cr.table(r)
    assert np.isfinite(r["worst"]) and r["worst"] <= 1.0
