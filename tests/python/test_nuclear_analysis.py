"""Automatic analysis: from peak areas to B(E2) and the shape (backlog 56)."""

import json
import math
from dataclasses import asdict, replace
from types import SimpleNamespace

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.analysis import Analysis, PeakFit, Result, Settings, fit_peak, fit_tails
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


@pytest.fixture(scope="module")
def projectile():
    """A sample that the fit failed on when it had a straight line for its background (seed 2, as seeds 1, 4 and 6)."""
    exp = projectile_case()
    return exp, Analysis(exp, simulate_gammas(exp, EVENTS, seed=2))


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


def peak_with_tails(n: float, share: float, edges: np.ndarray) -> np.ndarray:
    """The counts expected from a peak of ``n`` at 1172.9 keV, ``share`` of it 2.3 times as wide as the rest (σ =
    2.8 keV), on a background that its own scattered γ rays leave below it."""
    x, width = 1e3 * (edges[:-1] + edges[1:]) / 2 - 1172.9, 1e3 * (edges[1] - edges[0])
    shape = sum(p * np.exp(-0.5 * (x / s) ** 2) / (s * math.sqrt(2 * math.pi)) for p, s in
                ((1 - share, 2.8), (share, 6.4))) * width
    return n * shape + 1.7e-3 * n * width * (1 - np.cumsum(shape)) + 1e-4 * n * width


def fitted(counts, edges, variance=None) -> PeakFit:
    first = fit_peak(counts, edges, 1.1729, 0.003, 12.0, max_sigma=0.035, variance=variance)
    return fit_peak(counts, edges, 1.1729, first.sigma, 4.0, max_sigma=0.035, variance=variance)


def test_a_peak_with_tails_keeps_its_area():
    """The Gaussian's ±3σ leave out more than 0.27 % of such a peak, and its line takes the tails for background."""
    edges = np.arange(1.10, 1.25, 0.0002)
    expected = peak_with_tails(5000, 0.12, edges)
    rng = np.random.default_rng(4)
    plain, whole = [], []
    for _ in range(20):
        counts = rng.poisson(expected).astype(float)
        fit = fitted(counts, edges)
        full = fit_tails(counts, edges, fit)
        assert full.tails is not None and 0.965 < full.tails["inside"] < 0.99
        assert full.centroid == fit.centroid and full.sigma == fit.sigma   # the window stays that of the Gaussian
        plain.append(fit.area)
        whole.append(full.area)
    assert np.mean(plain) < 0.985 * 5000
    assert np.mean(whole) == pytest.approx(5000, abs=3 * np.std(whole) / math.sqrt(20) + 15)
    assert np.std(whole) == pytest.approx(full.area_unc, rel=0.5)
    # The curve is the fit's: it holds the peak and its background; and it survives being kept as JSON.
    e = (edges[:-1] + edges[1:]) / 2
    assert float(np.sum(full.curve(e) - full.background_curve(e))) == pytest.approx(full.tails["area"], rel=2e-3)
    inside = np.abs(e - full.centroid) < 3 * full.sigma
    assert float(full.background_curve(e)[inside].sum()) == pytest.approx(full.background, rel=0.02)
    assert PeakFit.of(dict(json.loads(json.dumps(asdict(full))), later="unknown")) == full


def test_a_weighted_sample_with_tails_keeps_its_area():
    """Entries that stand for many counts each are fitted with their own variance."""
    edges = np.arange(1.10, 1.25, 0.0002)
    rng = np.random.default_rng(5)
    areas = []
    for _ in range(10):
        entries = rng.poisson(peak_with_tails(20_000, 0.12, edges)).astype(float)
        full = fit_tails(50 * entries, edges, fitted(50 * entries, edges, 2500 * entries), variance=2500 * entries)
        assert full.tails is not None and full.tails["gain"] < 2000   # that of 20 000 entries, not of 10⁶ counts
        areas.append(full.area)
    assert np.mean(areas) == pytest.approx(1e6, rel=0.006)


@pytest.mark.parametrize("n", [300, 3000, 30_000])
def test_a_gaussian_has_no_tails(n):
    """A peak without tails, whatever its background, keeps the Gaussian's fit; so does one of too few counts."""
    edges = np.arange(1.10, 1.25, 0.0002)
    rng = np.random.default_rng(6)
    for _ in range(10):
        counts = rng.poisson(peak_with_tails(n, 0.0, edges)).astype(float)
        fit = fitted(counts, edges)
        assert fit_tails(counts, edges, fit) is fit


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
        assert [s["step"] for s in r.steps] == ["particle gate", "Doppler correction", "peak search", "peak fit",
                                                "random coincidences", "yield", "normalisation"]
        assert r.fit.randoms > 0 and r.area < r.fit.area
    assert a.history[-1].statistical < a.history[-2].statistical    # the gate costs counts


def test_b_e2_agrees_with_the_input_with_the_target_normalisation(projectile):
    """The ¹⁶O projectile is excited and the ⁵⁸Ni 2⁺ state (known B(E2)) is the reference. The line's Compton edge is in
    the fit window. Over 24 other seeds, of samples this large, area over true counts scatters by 5 % (half of it the true
    counts of the peak's 400 effective entries, which the area of the peak and its continuum does not follow), so a
    single one is judged by 3 % and three times that; with a straight line it was 26 to 33 % too large, or failed."""
    exp, a = projectile
    g = a.gammas
    ref = {"energy": "1454.2 keV", "b_up": "0.0695 e2b2", "unc": 0.03}
    both = [a.run(Settings(normalisation="target", reference=ref)), a.run(Settings())]
    full = a._gate(Settings())[0] & (np.abs(g["deposited"] - g["energy_lab"]) < 1e-6)
    w = g["weight"][full] * g.events.beam_time_s * g.live_fraction
    x = g["corrected_projectile"][full]
    centre = float(np.average(x, weights=w))
    neff = float(w.sum() ** 2 / (w**2).sum())
    for r in both:
        assert r.fit.failed is None and r.fit.below is not None and r.fit.shape.endswith("Compton continuum")
        assert abs(r.fit.area / w.sum() - 1) < 0.03 + 3 / math.sqrt(neff), (r.fit.area, w.sum(), neff)
        assert abs(r.fit.centroid - centre) < 0.005 + 3 * r.fit.sigma / math.sqrt(neff)
        assert abs(r.b_e2fm4 / r.truth_e2fm4 - 1) < 0.03 + 3 * max(r.monte_carlo, 1 / math.sqrt(neff)), r.b_e2fm4
        assert abs(r.pull) < 2.5, (r.normalisation, r.b_e2fm4, r.statistical, r.systematic)
        # What the response gives below the peak for a count in it has its line in the budget.
        assert 0.02 < r.budget["detector response below the peak"] < 0.06
        assert any("Compton edge" in n for n in r.notes)
    assert both[0].budget["reference B(E2)"] == pytest.approx(0.03)
    assert any("reference peak" in n for n in both[0].notes)
    # The two normalisations agree with each other more closely than with the truth: the same peak area.
    assert both[0].b_e2fm4 == pytest.approx(both[1].b_e2fm4, rel=0.05)
    with pytest.raises(ValueError, match="or Settings.reference"):
        a.run(Settings(normalisation="target"))


def test_a_fit_with_the_edge_is_kept_and_drawn_as_it_was_fitted(projectile):
    _, a = projectile
    r = a.run(Settings())
    again = Result.from_dict(json.loads(json.dumps(r.to_dict())))
    assert again.fit == r.fit and again.fit.below["area"] == r.fit.area
    edges, counts = np.array(r.spectrum["edges"]), np.array(r.spectrum["counts"])
    e = (edges[:-1] + edges[1:]) / 2
    half = 4 * 2.3548 * r.fit.sigma
    window = np.abs(e - a.e0) < half
    assert r.fit.curve(e[window]).sum() == pytest.approx(counts[window].sum(), rel=2e-3)
    under = r.fit.background_curve(e) - r.fit.level - r.fit.slope * (e - r.fit.origin)
    # The continuum below the peak, and none above it; more below the edge than between it and the peak.
    flat, below = (float(under[np.abs(e - x) < 0.01].mean()) for x in (6.76, 6.60))
    assert 0 < flat < below and abs(under[e > r.fit.centroid + 4.5 * r.fit.sigma]).max() < 1e-3 * flat
    line = r.lines[0]
    found = [pk for pk in r.search["projectile"]["peaks"] if pk["label"]]
    assert len(found) == 1 and found[0]["area"] == pytest.approx(r.fit.area, rel=0.03)  # the search has it too
    # A result kept before there was such a field has a straight line.
    old = r.to_dict()
    del old["fit"]["below"]
    assert Result.from_dict(old).fit.below is None
    pytest.importorskip("plotly")
    from physim.nuclear.app.figures import figure_peak

    traces = {t.name: t for t in figure_peak(line).data}
    drawn = traces["the background under the peak"]
    assert np.allclose(np.asarray(drawn.y, dtype=float), 4 * r.fit.background_curve(np.asarray(drawn.x) * 1e-3))
    assert np.allclose(np.asarray(traces[f"fit: {r.fit.shape}"].y, dtype=float),
                       4 * r.fit.curve(np.asarray(drawn.x) * 1e-3))


def test_a_run_of_few_counts_is_fitted_with_the_edge(projectile):
    """24 h of beam, 150 counts in the peak, each once: the Poisson fit, and the simulated peak of the same setup."""
    exp, a = projectile
    g = a.gammas
    real = g.events.beam_time_s
    rows = np.repeat(np.arange(len(g)), np.random.default_rng(5).poisson(g["weight"] * real))
    cols = {k: np.asarray(v)[rows] for k, v in g.columns.items()}
    cols["weight"] = np.full(len(rows), 1.0 / real)
    b = Analysis(exp, replace(g, columns=cols), run=SimpleNamespace(summary={"one_particle": True}, scale=1.0,
                                                                    real_s=real, duration_s=real))
    s = Settings(subtract_randoms=False)
    keep = b._gate(s)[0]
    true = int((keep & (np.abs(cols["deposited"] - cols["energy_lab"]) < 1e-6)).sum())
    counts, edges, _, sigma, _ = b._spectrum(s, keep)
    assert b._variance is None and 100 < true < 200
    fit = b._peak(s, 0, b.e0, counts, edges, sigma, None)
    assert fit.failed is None and fit.below is not None
    # Scatter over 48 runs: 10 % of the true counts, 5 % of the area itself, as it prints.
    assert fit.area == pytest.approx(true, rel=0.03 + 3 * 0.10) and 0.03 < fit.area_unc / fit.area < 0.08
    # A Gaussian of any starting width: the detectors' resolution (a single entry), or that of the edge.
    for start in (0.0023, 0.15):
        assert b._peak(s, 0, b.e0, counts, edges, start, None).area == pytest.approx(fit.area, rel=1e-3)


def test_where_no_edge_is_in_the_window_the_fit_is_what_it_was(ni, monkeypatch):
    """The fit of a 1454 keV line is the two passes and the tails, with a line and a step, to the last digit: the edge is
    217 keV below it."""
    _, a = ni
    for gate in ("inelastic", "all"):
        s = Settings(particle_gate=gate)
        counts, edges, _, sigma, _ = a._spectrum(s, a._gate(s)[0])
        var = a._variance
        first = fit_peak(counts, edges, a.e0, sigma, 12.0, max_sigma=0.03 * a.e0, variance=var)
        second = fit_peak(counts, edges, a.e0, first.sigma, 4.0, max_sigma=0.03 * a.e0, variance=var)
        expected = fit_tails(counts, edges, second, variance=var)
        fit = a._peak(s, 0, a.e0, counts, edges, sigma, var)
        assert fit == expected and fit.below is None and not fit.failed
        assert a._edge_fit(s, 0, a.e0, counts, edges, var, fit) is None
        assert a._below_change(s, 0, fit) == 0.0
    # A spectrum that does not reach down to the edge: the detectors' response is not even asked.
    monkeypatch.setattr(Analysis, "_below", lambda *args: pytest.fail("not needed"))
    assert a._edge_fit(s, 0, a.e0, counts, edges, var, fit) is None


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
    assert set(r.budget) == {"efficiency", "angular correlation", "beam energy", "detector positions"}
    assert r.budget["efficiency"] == 0.05 and r.budget["angular correlation"] == 0.03
    assert r.systematic == pytest.approx(math.sqrt(sum(v**2 for v in r.budget.values())))
    # Beam energy: ⟨P⟩/B with the energy the beam enters with raised by δ, against nominal, as the budget computes it.
    p0 = a._mean_probability(Settings(particle_gate="all"), exact=True)
    p1 = a._mean_probability(Settings(particle_gate="all"), energy_scale=1.005, exact=True)
    assert r.budget["beam energy"] == pytest.approx(abs(p1 / p0 - 1))
    assert 0.01 < r.budget["beam energy"] < 0.1
    # Doubling the assumed efficiency uncertainty doubles its entry and grows the total.
    r2 = a.run(Settings(particle_gate="all", efficiency_unc=0.10))
    assert r2.budget["efficiency"] == 0.10 and r2.systematic > r.systematic
    assert r2.total_unc == pytest.approx(math.hypot(r2.statistical, r2.systematic))


def test_an_input_without_an_uncertainty_has_no_term_and_the_budget_says_so(ni):
    """Round 3, fix 12: the target's thickness and the matrix elements that are not measured."""
    exp, a = ni
    r = a.run(Settings(particle_gate="all"))
    assert r.settings["thickness_unc"] is None and "target thickness" not in r.budget
    assert "matrix elements assumed" not in r.budget       # it was there as 0, with nothing known about it
    assert r.budget_notes[0].startswith("target thickness: no uncertainty given")
    assert r.budget_notes[1].startswith("matrix elements assumed: the setup has no level scheme of 58Ni")
    thick = a.run(Settings(particle_gate="all", thickness_unc=0.1))
    assert 0 < thick.budget["target thickness"] < 0.1 and thick.systematic > r.systematic
    assert thick.systematic == pytest.approx(math.sqrt(sum(v**2 for v in thick.budget.values())))
    assert thick.steps[-1]["step"] == "normalisation" and "thicker" in thick.steps[-1]["text"]
    assert "thicker" not in r.steps[-1]["text"] and any("thicker" in n for n in thick.budget_notes)
    # A result kept by an earlier version has neither the setting nor the notes, and still loads.
    d = json.loads(json.dumps(r.to_dict()))
    del d["budget_notes"], d["settings"]["thickness_unc"]
    d["budget"]["matrix elements assumed"] = 0.0
    old = Result.from_dict(d)
    assert old.budget_notes is None and old.total_unc == r.total_unc
    assert Settings(**old.settings).thickness_unc is None
    assert Result.from_dict(json.loads(json.dumps(thick.to_dict()))).budget_notes == thick.budget_notes


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
