"""The published lifetime reproductions and the simulated runs of the lifetime checks (lifetime wave 3C):
physim.nuclear.lifetime_reproduction, the three examples made after the papers, and the simulated series NaPyTau
analysed (tests/reference/nuclear/lifetimes/). The reproductions themselves run in test_nuclear_validation.py."""

import math
from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import gamma
from physim.nuclear import lifetime_analysis as la
from physim.nuclear import lifetime_reproduction as lr
from physim.nuclear import lifetimes as lt
from physim.nuclear.experiment import Experiment
from physim.nuclear.gamma import level_names

REFERENCE = Path(__file__).resolve().parents[1] / "reference" / "nuclear"
needs_reference = pytest.mark.skipif(not (REFERENCE / "lifetimes").is_dir(), reason="no reference files")


# -- the examples, as the papers give them ------------------------------------------------------------------------


def test_ru98_example_takes_the_papers_lifetime_not_ensdf():
    exp = Experiment.example("coulex_ru98_plunger")
    taus = lt.level_taus(gamma.excitation_of(exp).scheme)
    assert list(taus.values())[0][0] == pytest.approx(8.36, abs=0.005)


def test_cr50_example_level_scheme():
    scheme = Experiment.example("fusion_cr50_plunger").levels["50Cr"]
    taus = lt.level_taus(scheme)
    assert taus["2+"][0] == pytest.approx(12.2, rel=1e-4) and taus["4+"][0] == pytest.approx(2.4, rel=1e-4)
    assert taus["6+"][0] == pytest.approx(0.8 / math.log(2), rel=1e-4)  # ENSDF's 0.80 ps half-life
    assert [t.energy.value for t in scheme.transitions] == [783.3, 1098.1, 1282.5]


def test_sm143_example_is_the_papers_level_scheme():
    """Rajbanshi et al. 2018: DB I on the 43/2⁻ level at 8614 keV, DB II forked at 49/2⁻; energies as sums of the
    transitions, Table I's lifetimes and branchings."""
    exp = Experiment.example("fusion_sm143_dsam")
    scheme = exp.levels["143Sm"]
    labels = lr._sm143_labels(scheme)
    energy = {lab: scheme.levels[i].energy.value for i, lab in enumerate(labels) if lab}
    assert energy["45/2- (DB I)"] == pytest.approx(8614 + 239.5)
    assert energy["51/2- (DB I)"] == pytest.approx(8614 + 239.5 + 339.6 + 444.8 + 578.2)
    assert energy["51/2- (DB II)"] == pytest.approx(8614 + 239.5 + 339.6 + 444.8 + 527.5)
    assert energy["57/2(-) (DB II)"] == pytest.approx(energy["51/2- (DB II)"] + 584.9 + 646.2 + 698.0)
    taus = lr.lt_level_taus_array(scheme)
    index = {lab: i for i, lab in enumerate(labels) if lab}
    paper = {r["state"]: r["tau_ps"] for r in lr.published(REFERENCE)["lifetimes"]
             if r["case"] == "sm143_dsam_fusion"} if (REFERENCE / "lifetimes").is_dir() else {}
    for lab, tau in paper.items():
        assert taus[index[lab]] == pytest.approx(tau, rel=1e-4), lab
    # Every crossover E2 matches the sum of the two dipoles it bypasses (Table I) within 0.5 keV.
    e = {t.initial: {} for t in scheme.transitions}
    for t in scheme.transitions:
        e[t.initial][t.final] = t.energy.value
    for t in scheme.transitions:
        if t.multipolarity == "E2":
            lev_e = scheme.levels[t.initial].energy.value - scheme.levels[t.final].energy.value
            assert t.energy.value == pytest.approx(lev_e, abs=0.5)
    assert exp.fusion.channels[1].entry_spin == 28.5


# -- the pieces ----------------------------------------------------------------------------------------------------


def test_recoil_alignment_matches_the_correlation_table_for_both_particles():
    exp = Experiment.example("coulex_ru98_plunger")
    both = lr.coulex_alignment(exp, "both")
    iso = Experiment.example("coulex_ru98_plunger")
    iso.excitation.emission = "isotropic"
    corr = {(r["gamma_detector"], r["crystal"]): r["factor"] for r in gamma.correlation_table(exp)}
    flat = {(r["gamma_detector"], r["crystal"]): r["factor"] for r in gamma.correlation_table(iso)}
    f1 = np.mean([corr[k] / flat[k] for k in corr if k[0].startswith("F")])
    assert both["w"][41.5] == pytest.approx(f1, rel=1e-9)
    rec = lr.coulex_alignment(exp, "recoil")
    # Gated on the ²⁴Mg recoils (head-on collisions) the 2⁺ is much more aligned: W(90°) about 0.5, not 0.8.
    assert 0.4 < rec["w"][90.0] < 0.6 and 0.75 < both["w"][90.0] < 0.85
    assert rec[2] > both[2] > 0 and rec[4] < both[4] < 0


def test_ring_cos_is_below_the_nominal_angle():
    exp = Experiment.example("coulex_ru98_plunger")
    c = lr.ring_cos(exp, (41.5, 138.5))
    assert 0.70 < c[41.5] < math.cos(math.radians(41.5)) and c[138.5] == pytest.approx(-c[41.5], rel=1e-6)


def test_window_areas_and_efficiency_ratio():
    """A plunger spectrum built by hand: a stopped peak at E₀, a flat-topped shifted line, a flat background and a
    shifted line seen 7% more efficiently. The windows give the components, the efficiency ratio comes back."""
    e0, beta, ring = 652.0, 0.045, 41.5
    edges = np.arange(e0 - 45.0, e0 + 45.0 + 1e-9, 0.5)
    c = 0.5 * (edges[1:] + edges[:-1])
    shares = np.array([0.9, 0.6, 0.3, 0.1])
    k = 1.07
    spectra = []
    for u in shares:
        stopped = 4000 * u * np.exp(-0.5 * ((c - e0) / 0.9) ** 2) / (0.9 * math.sqrt(2 * math.pi)) * 0.5
        flat = (np.abs(c - (e0 + 20.0)) <= 8.0) * (4000 * k * (1 - u) / 32.0)
        spectra.append(stopped + flat + 3.0)
    run = {"edges_kev": edges, "e0_kev": e0, "spectra": {ring: np.array(spectra)}}
    a = lr.plunger_areas(run, ring, beta)
    assert a["unshifted"] == pytest.approx(4000 * shares, rel=1e-3)
    assert a["shifted"] == pytest.approx(4000 * k * (1 - shares), rel=1e-3)
    assert lr.efficiency_ratio(a["unshifted"], a["shifted"], a["err_unshifted"], a["err_shifted"]) == \
        pytest.approx(k, rel=1e-3)
    with pytest.raises(ValueError):
        lr.plunger_areas(run, ring, beta, method="other")


def test_feeders_and_the_sm143_cascade():
    scheme = Experiment.example("fusion_sm143_dsam").levels["143Sm"]
    side = {lab: 0.15 for lab in lr.SM143_INTENSITY}
    cas = lr.sm143_cascade(scheme, side, steps=1)
    n = len(scheme.levels)
    assert len(cas["taus_ps"]) == n + 9 and cas["populations"].sum() == pytest.approx(1.0)
    assert np.all(cas["populations"][:n] == 0)  # every band level is fed through its side-feeding level
    assert cas["e_gamma"]["49/2- (DB I)"] == 444.8 and cas["e_gamma"]["51/2- (DB II)"] == 527.5
    top = cas["index"]["57/2(-) (DB II)"]
    assert list(lr._feeders(cas["feeding"], top)) == [top, n + list(cas["index"]).index("57/2(-) (DB II)")]
    low = lr._feeders(cas["feeding"], cas["index"]["45/2- (DB I)"])
    assert len(low) == 2 * 9  # every band level and every side-feeding level reaches the 45/2⁻
    five = lr.sm143_cascade(scheme, side, steps=5)
    assert len(five["taus_ps"]) == n + 45 and five["taus_ps"][n] == pytest.approx(0.03)


def test_pair_stack_sets_both_gaps():
    exp = Experiment.example("fusion_cr50_plunger")
    segs = lr._pair_stack(exp, 42.0, 40.0)
    names = [s.name for s in segs]
    assert names == ["target", "gap", "degrader", "gap2", "degrader2"]
    assert segs[names.index("gap2")].gap_um == 40.0
    assert segs[names.index("gap")].gap_um == pytest.approx(42.0 + exp.lifetime.offset_um)


# -- the reference files -------------------------------------------------------------------------------------------


@needs_reference
def test_published_values_are_read():
    pub = lr.published(REFERENCE)
    assert lr._paper(pub, "ru98_rdds_coulex", "2+1")["tau_ps"] == 8.36
    assert lr._paper(pub, "cr50_rdds_fusion", "4+1")["plus_ps"] == 0.3
    assert pub["setups"]["sm143_dsam_fusion"]["side_feeding_ps"]["51/2- (DB II)"] == 0.15
    assert lr.published(REFERENCE / "nowhere") is None


@needs_reference
def test_napytau_on_the_simulated_series():
    """NaPyTau's τ(x) on physim's simulated areas, and physim's DDCM in NaPyTau's procedure on the same areas."""
    sim = lr.read_napytau_simulated(REFERENCE)
    assert len(sim["distances_um"]) == 11 and sim["offset_um"] == 8.0
    assert sim["tau_ps"] == pytest.approx(8.36, abs=1e-3) and 0.04 < sim["beta"] < 0.055
    res = lr.ddcm_vs_napytau(REFERENCE)
    assert res["worst"] < 1e-4
    # A deliberately different procedure (weighted fit) is visibly not NaPyTau's.
    r = la.ddcm(sim["distances_um"], sim["unshifted"], sim["shifted"], (sim["err_unshifted"], sim["err_shifted"]),
                sim["beta"], polynomial=3, fit="shifted", weighted=True, proper=False, offset_um=sim["offset_um"])
    ours = np.array([p["tau_ps"] for p in r.per_distance])
    assert np.nanmax(np.abs(ours / sim["napytau"][3] - 1)) > 0.01


@needs_reference
def test_one_level_closure():
    c = lr.closure_one_level(REFERENCE)
    for tau, stat in (c["decay_curve"], c["ddcm"]):
        assert abs(tau - c["tau_ps"]) <= 2 * stat
    assert 1.0 < c["efficiency_ratio"] < 1.15  # the forward ring sees the shifted line boosted
