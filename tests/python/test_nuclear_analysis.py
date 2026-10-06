"""Automatic analysis: from peak areas to B(E2) and the shape (backlog 56)."""

import math

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.analysis import Analysis, Settings, fit_peak
from physim.nuclear.gamma_events import simulate_gammas
from physim.nuclear.levels import weisskopf_unit
from physim.nuclear.planner import Planner

EVENTS = 400_000


@pytest.fixture(scope="module")
def ni():
    exp = Experiment.example("coulex_ni58")
    return exp, Analysis(exp, simulate_gammas(exp, EVENTS, seed=1))


def projectile_case() -> Experiment:
    d = Experiment.example("coulex_ni58").to_dict()
    d["beam"]["energy"] = "70 MeV"
    d["reaction"] = {"type": "coulex", "excite": "projectile", "energy": "6.917 MeV", "multipolarity": "E2",
                     "b_up": "0.004 e2b2"}
    return Experiment.from_dict(d)


# -- the peak fit ---------------------------------------------------------------------------------------------------


def test_the_peak_fit_recovers_a_known_peak():
    rng = np.random.default_rng(2)
    edges = np.arange(1.40, 1.51, 0.0005)
    centres = (edges[:-1] + edges[1:]) / 2
    sigma = 0.0025
    expected = 3000 * 0.0005 / (sigma * math.sqrt(2 * math.pi)) * np.exp(-0.5 * ((centres - 1.454) / sigma) ** 2)
    expected += 40 + 200 * (centres - 1.454)
    counts = rng.poisson(expected).astype(float)
    fit = fit_peak(counts, edges, 1.454, sigma * 0.8, half_width=4.0)
    assert fit.area == pytest.approx(3000, abs=3 * fit.area_unc + 30)
    assert fit.centroid == pytest.approx(1.454, abs=3e-4) and fit.sigma == pytest.approx(sigma, rel=0.15)
    assert fit.gauss_area == pytest.approx(3000, rel=0.08)
    assert 0.5 * fit.bins < fit.chi2 < 2 * fit.bins
    with pytest.raises(ValueError, match="too few bins"):
        fit_peak(counts[:3], edges[:4], 1.454, sigma)


# -- B(E2) ----------------------------------------------------------------------------------------------------------


def test_b_e2_agrees_with_the_input_with_the_rutherford_normalisation(ni):
    exp, a = ni
    truth = exp.excitation.b_up_e2fm
    for settings in (Settings(), Settings(particle_gate="all")):
        r = a.run(settings)
        assert abs(r.pull) < 2.5, (settings.particle_gate, r.b_e2fm4, truth, r.statistical, r.systematic)
        assert r.b_e2fm4 == pytest.approx(truth, rel=0.12)
        assert r.b_e2b2 == pytest.approx(r.b_e2fm4 * 1e-4) and r.b_wu == pytest.approx(r.b_e2fm4 / weisskopf_unit(2, 58))
        assert r.truth_e2fm4 == truth and r.normalisation == "rutherford"
        assert [s["step"] for s in r.steps] == ["particle gate", "Doppler correction", "peak fit",
                                                "random coincidences", "yield", "normalisation"]
        assert r.fit.randoms > 0 and r.area < r.fit.area
    assert a.history[-1].statistical < a.history[-2].statistical    # the gate costs counts


def test_b_e2_agrees_with_the_input_with_the_target_normalisation():
    """The ¹⁶O projectile is excited and the ⁵⁸Ni 2⁺ state (known B(E2)) is the reference."""
    exp = projectile_case()
    a = Analysis(exp, simulate_gammas(exp, EVENTS, seed=2))
    ref = {"energy": "1454.2 keV", "b_up": "0.0695 e2b2", "unc": 0.03}
    both = [a.run(Settings(normalisation="target", reference=ref)), a.run(Settings())]
    for r in both:
        assert abs(r.pull) < 2.5, (r.normalisation, r.b_e2fm4, r.statistical, r.systematic)
    assert both[0].budget["reference B(E2)"] == pytest.approx(0.03)
    assert any("reference peak" in n for n in both[0].notes)
    # The two normalisations agree with each other more closely than with the truth: the same peak area.
    assert both[0].b_e2fm4 == pytest.approx(both[1].b_e2fm4, rel=0.05)
    with pytest.raises(ValueError, match="needs Settings.reference"):
        a.run(Settings(normalisation="target"))


def test_the_scatter_over_seeds_matches_the_quoted_statistical_uncertainty():
    exp = Experiment.example("coulex_ni58")
    values, quoted, mc = [], [], []
    for seed in range(1, 7):
        r = Analysis(exp, simulate_gammas(exp, 150_000, seed=seed)).run(Settings(particle_gate="all"))
        values.append(r.b_e2fm4)
        quoted.append(r.statistical * r.b_e2fm4)
        mc.append(r.monte_carlo * r.b_e2fm4)
    scatter = float(np.std(values, ddof=1))
    # The scatter between seeds is the Monte Carlo's own uncertainty, which the result quotes beside the
    # experiment's statistical one (six seeds give the scatter to about 30%).
    assert 0.5 * np.mean(mc) < scatter < 2.0 * np.mean(mc), (values, mc)
    # The experiment's statistical uncertainty is that of the counts in the run, smaller than a 150 000-event
    # sample's; the two would meet with a sample as large as the run.
    assert np.mean(quoted) < np.mean(mc)
    r = Analysis(exp, simulate_gammas(exp, 150_000, seed=1)).run(Settings(particle_gate="all"))
    assert r.monte_carlo == pytest.approx(np.mean(mc) / r.b_e2fm4, rel=0.3)


def test_each_systematic_is_found_by_changing_its_input(ni):
    exp, a = ni
    r = a.run(Settings(particle_gate="all"))
    assert set(r.budget) == {"efficiency", "angular correlation", "beam energy", "detector positions",
                             "matrix elements assumed"}
    assert r.budget["efficiency"] == 0.05 and r.budget["angular correlation"] == 0.03
    assert r.systematic == pytest.approx(math.sqrt(sum(v**2 for v in r.budget.values())))
    # Beam energy: ⟨P⟩/B at E(1 + δ) against E, as the budget computes it.
    p0 = a._mean_probability(Settings(particle_gate="all"), None)
    p1 = a._mean_probability(Settings(particle_gate="all"), None, beam_energy=a._e_mid * 1.005)
    assert r.budget["beam energy"] == pytest.approx(abs(p1 / p0 - 1))
    assert 0.01 < r.budget["beam energy"] < 0.1
    # Doubling the assumed efficiency uncertainty doubles its entry and grows the total.
    r2 = a.run(Settings(particle_gate="all", efficiency_unc=0.10))
    assert r2.budget["efficiency"] == 0.10 and r2.systematic > r.systematic
    assert r2.total_unc == pytest.approx(math.hypot(r2.statistical, r2.systematic))


def test_a_changed_gate_or_window_changes_the_result_and_is_recorded(ni):
    exp, a = ni
    base = a.run(Settings())
    narrow = a.run(Settings(fit_half_width=2.0))
    assert narrow.settings["fit_half_width"] == 2.0
    assert any("Changed from the previous run" in n and "fit_half_width" in n for n in narrow.notes)
    rings = a.run(Settings(rings=list(range(8, 16))))
    assert rings.steps[0]["text"].startswith("inelastic group") and "rings [8" in rings.steps[0]["text"]
    assert rings.area < base.area and rings.norm_counts < base.norm_counts
    assert any("'rings'" in n for n in rings.notes)
    only_cd = a.run(Settings(detectors=["CD"], particle_gate="all"))
    assert only_cd.norm_counts < base.norm_counts and abs(only_cd.pull) < 3
    assert len(a.history) >= 4 and a.history[-1] is only_cd


def test_the_doppler_correction_for_the_wrong_nucleus_is_a_choice(ni):
    exp, a = ni
    right = a.run(Settings(correction="recoil"))
    wrong = a.run(Settings(correction="projectile"))
    assert right.steps[1]["text"] == "for the recoil" and wrong.steps[1]["text"] == "for the projectile"
    assert wrong.fit.fwhm_kev > right.fit.fwhm_kev


def test_shape_readings_and_beam_time(ni):
    exp, a = ni
    r = a.run(Settings(particle_gate="all"))
    b = r.b_e2fm4
    z, mass = 28, 58
    assert r.shape["beta2"] == pytest.approx(4 * math.pi / (3 * z * (1.2 * mass ** (1 / 3)) ** 2) * math.sqrt(b))
    assert r.shape["q0_efm2"] == pytest.approx(math.sqrt(16 * math.pi * b / 5))
    assert r.shape["qs_2plus_efm2"] == pytest.approx(-2 / 7 * r.shape["q0_efm2"])
    assert "rotor" in r.shape["note"]
    t = a.gammas.events.beam_time_s * a.gammas.live_fraction
    assert r.counts_per_shift == pytest.approx(r.area / t * 8 * 3600)
    assert r.hours_for_precision == pytest.approx((1 / 0.05**2) / (r.area / t) / 3600)
    precise = a.run(Settings(particle_gate="all", wanted_precision=0.01))
    assert precise.hours_for_precision == pytest.approx(25 * r.hours_for_precision)


def test_the_planner_runs_the_analysis():
    p = Planner.example("coulex_ni58")
    r = p.analyse({"particle_gate": "all"}, 150_000, 1)
    assert r["available"] and abs(r["pull"]) < 3 and r["runs"] == 1
    r2 = p.analyse(Settings(fit_half_width=3.0), 150_000, 1)
    assert r2["runs"] == 2 and any("Changed" in n for n in r2["notes"])
    assert not p.analyse({"normalisation": "target"}, 150_000, 1)["available"]
    assert not Planner.example("alpha_on_gold").analyse()["available"]
