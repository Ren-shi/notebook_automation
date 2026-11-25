"""The peak fit of the analysis at few counts and at many (Allmond 2014 verification, round 3, fixes 3 and 4): Poisson
likelihood, bounds and failed fits, the tail shape of a simulated peak, and the area at high statistics.

Before: least squares with the variances 1/max(counts, 1) and a Gaussian counted within ±3σ lost 1 to 4 % of a peak
of 60 to 1400 counts with the tails of a thick target (one fit in a thousand ran away at 60 counts), the share of two
Gaussians within ±3σ was 0.2 % too large for one of a million, and a peak on a background gained several percent."""

import json
import math
from dataclasses import asdict

import numpy as np
import pytest

from physim.nuclear import Experiment
from physim.nuclear.analysis import SHAPE_REACH, TAILS_WORTH, PeakFit, _biweight, fit_peak, fit_tails, tail_shape
from physim.nuclear.response import ELECTRON_MEV, Response, compton_edge

E0, EDGES = 1.1729, np.arange(1.10, 1.25, 0.0002)
X = 1e3 * (EDGES[:-1] + EDGES[1:]) / 2 - 1e3 * E0  # keV from the transition


def gaussian(sigma: float, centre: float = 0.0) -> np.ndarray:
    """The share of a Gaussian (keV) in each bin."""
    return np.exp(-0.5 * ((X - centre) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi)) * 0.2


def box(half: float, sigma: float) -> np.ndarray:
    """The share of a flat top ±``half`` keV wide, smeared with σ, in each bin."""
    from math import erf

    f = np.vectorize(lambda x: (erf((x + half) / sigma / math.sqrt(2)) - erf((x - half) / sigma / math.sqrt(2))) / 2)
    return f(X) / (2 * half) * 0.2


def thick_target_peak(n: float) -> np.ndarray:
    """The counts expected of a peak as a thick target and ring segments leave it: a flat-topped core, and tails that
    end (12 % in a wider flat top, 0.9 keV below), on the step its own scattered γ rays make and little else."""
    shape = 0.88 * box(2.6, 2.3) + 0.12 * np.interp(X + 0.9, X, box(8.5, 2.6))
    shape = shape / shape.sum()
    return n * shape + 1.7e-3 * 0.2 * n * (1 - np.cumsum(shape) + shape / 2)


def gauss_fit(counts, variance=None) -> PeakFit:
    """As Analysis does: two passes, the window of the second from the width of the first."""
    first = fit_peak(counts, EDGES, E0, 0.003, 12.0, max_sigma=0.035, variance=variance)
    return fit_peak(counts, EDGES, E0, 0.003 if first.failed else first.sigma, 4.0, max_sigma=0.035, variance=variance)


def peak_fit(counts, shape=None) -> PeakFit:
    """As Analysis does: with its own tails, else the simulated ones."""
    fit = gauss_fit(counts)
    own = fit_tails(counts, EDGES, fit)
    return own if own.tails is not None or shape is None else fit_tails(counts, EDGES, fit, shape=shape)


@pytest.fixture(scope="module")
def shape() -> dict:
    """The tails of a "simulated" peak of 15 000 counts (its own draw), as a weighted sample has them."""
    entries = np.random.default_rng(101).poisson(thick_target_peak(15_000)).astype(float)
    fit = fit_tails(0.1 * entries, EDGES, gauss_fit(0.1 * entries, 0.01 * entries), variance=0.01 * entries)
    assert fit.tails is not None and not fit.tails["simulated"] and fit.failed is None
    return tail_shape(fit)


# -- at few counts --------------------------------------------------------------------------------------------------


@pytest.mark.parametrize("n", [60, 140, 1400, 5000])
def test_the_area_of_few_counts_is_right_with_the_simulated_tails(shape, n):
    """Within 0.5 %, or the toys' own error; no fit runs away and none fails; the quoted error is the scatter."""
    rng = np.random.default_rng(n)
    mu = thick_target_peak(n)
    toys = [peak_fit(rng.poisson(mu).astype(float), shape) for _ in range(250)]
    area = np.array([f.area for f in toys]) / n
    err = area.std() / math.sqrt(len(area))
    assert all(f.failed is None and f.tails is not None for f in toys)
    assert abs(area.mean() - 1) < max(0.005, 2.5 * err), (area.mean(), err)
    assert np.all(np.abs(area - 1) < 6 / math.sqrt(n))
    assert area.std() == pytest.approx(np.mean([f.area_unc for f in toys]) / n, rel=0.12)
    if n < 1000:
        assert np.mean([f.tails["simulated"] for f in toys]) > 0.98


def test_the_area_has_no_step_where_the_tails_come_from_the_counts(shape):
    """From all simulated to all their own: the mean area stays."""
    own, mean, rng = [], [], np.random.default_rng(8)
    for n in (700, 1500, 3000, 6000):
        mu = thick_target_peak(n)
        toys = [peak_fit(rng.poisson(mu).astype(float), shape) for _ in range(120)]
        own.append(np.mean([not f.tails["simulated"] for f in toys]))
        mean.append(np.mean([f.area for f in toys]) / n)
        assert abs(mean[-1] - 1) < 0.005 + 2.5 / math.sqrt(n * len(toys)), (n, mean, own)
    assert own[0] < 0.5 < own[-1]


@pytest.mark.parametrize("n, level", [(140, 0.0), (140, 0.5), (1400, 0.5), (1400, 20.0)])
def test_a_gaussian_peak_on_a_background_keeps_its_area(n, level):
    """The line beside the peak is right at any level: with the variances of the counts it was low where the bins hold
    few (a peak of 1400 on 20 per bin read 6 % high)."""
    rng = np.random.default_rng(n)
    toys = [gauss_fit(rng.poisson(n * gaussian(2.9) + level).astype(float)) for _ in range(300)]
    area = np.array([f.area for f in toys]) / n
    err = area.std() / math.sqrt(len(area))
    assert abs(area.mean() - 1) < 0.005 + 2.5 * err, (area.mean(), err)
    assert all(f.tails is None and f.failed is None for f in toys)


def test_a_gaussian_peak_fitted_with_a_gaussian_is_the_same_as_before():
    """Many counts: likelihood and least squares agree far within the statistical error."""
    counts = np.random.default_rng(3).poisson(50_000 * gaussian(2.9) + 100).astype(float)
    ml, ls = gauss_fit(counts), gauss_fit(counts, np.maximum(counts, 1.0))
    assert ml.area == pytest.approx(ls.area, abs=0.2 * ml.area_unc)
    assert ml.sigma == pytest.approx(ls.sigma, rel=2e-3) and ml.area_unc == pytest.approx(ls.area_unc, rel=0.02)


# -- at many counts -------------------------------------------------------------------------------------------------


def test_the_area_of_many_counts_is_that_of_the_whole_peak():
    """A flat top and tails that end are no sum of two Gaussians; the area is the counts over the background of the
    whole window, which does not need them to be. (On the simulated spectrum itself: in the closure tests.)"""
    n = 1_000_000
    fit = peak_fit(thick_target_peak(n))  # the expected counts themselves: no scatter
    assert fit.tails is not None and not fit.tails["simulated"]
    assert fit.area == pytest.approx(n, rel=1e-3)
    # The line of the background is nowhere below zero.
    t, e = fit.tails, 1e-3 * X + E0
    window = np.abs(e - fit.centroid) < 4 * 2.3548 * fit.sigma
    assert np.all((t["level"] + t["slope"] * (e - t["origin"]))[window] > -1e-6 * n)
    # The curve is the fit's, and it survives being kept.
    assert float(np.sum(fit.curve(e) - fit.background_curve(e))) == pytest.approx(fit.area, rel=2e-3)
    assert PeakFit.of(json.loads(json.dumps(asdict(fit)))) == fit


def test_the_wider_component_ends():
    """The tails of a Doppler-corrected peak end; a bell of standard deviation 1 that ends at ±√7."""
    z = np.linspace(-4, 4, 80_001)
    density, above = _biweight(z)
    step = z[1] - z[0]
    assert density.sum() * step == pytest.approx(1, abs=1e-6)
    assert (density * z**2).sum() * step == pytest.approx(1, abs=1e-5)
    assert np.all(density[np.abs(z) > np.sqrt(7)] == 0) and np.all(density >= 0)
    assert above[0] == pytest.approx(1) and above[-1] == pytest.approx(0, abs=1e-12)
    assert 1 - np.cumsum(density) * step == pytest.approx(above, abs=1e-4)
    assert peak_fit(thick_target_peak(1_000_000)).tails["wide_shape"] == "biweight"


def test_a_fit_kept_before_the_wider_component_ended_has_a_gaussian_there():
    fit = peak_fit(thick_target_peak(100_000))
    kept = asdict(fit)
    del kept["tails"]["wide_shape"]
    old = PeakFit.of(json.loads(json.dumps(kept)))
    t, e = old.tails, 1e-3 * X + E0
    gauss = lambda c, s: np.exp(-0.5 * ((e - c) / s) ** 2) / (s * np.sqrt(2 * np.pi))  # noqa: E731
    peak = (1 - t["share"]) * gauss(t["centre"], t["sigma"]) + t["share"] * gauss(t["wide_centre"], t["wide_sigma"])
    assert old.curve(e) - old.background_curve(e) == pytest.approx(t["area"] * t["width"] * peak)
    assert not np.allclose(old.curve(e), fit.curve(e), rtol=1e-6)


# -- bounds and failures --------------------------------------------------------------------------------------------


def test_a_fit_at_its_bound_is_failed():
    rng = np.random.default_rng(9)
    broad = rng.poisson(300 * gaussian(40.0) + 0.05).astype(float)
    fit = fit_peak(broad, EDGES, E0, 0.003, 12.0, max_sigma=0.010)
    assert fit.failed and "width" in fit.failed and fit.sigma == pytest.approx(0.010)
    beside = rng.poisson(2000 * gaussian(2.9, 10.0) + 0.05).astype(float)
    fit = fit_peak(beside, EDGES, E0, 0.003, 4.0, max_sigma=0.006)
    assert fit.failed
    nothing = rng.poisson(np.where(X < -20, 3.0, 0.5)).astype(float)
    nothing[np.abs(X) < 8] = 0.0
    fit = fit_peak(nothing, EDGES, E0, 0.003, 4.0)
    assert fit.failed and "no counts above" in fit.failed
    assert gauss_fit(rng.poisson(300 * gaussian(2.9)).astype(float)).failed is None


def test_the_simulated_shape_bounds_the_width_and_the_centre(shape):
    rng = np.random.default_rng(10)
    broad = rng.poisson(400 * gaussian(14.0)).astype(float)
    fit = fit_tails(broad, EDGES, gauss_fit(broad), shape=shape)
    assert fit.failed and "width" in fit.failed
    assert fit.tails["sigma"] <= SHAPE_REACH * shape["core_sigma"] * (1 + 1e-9)
    beside = rng.poisson(400 * gaussian(2.9, 14.0)).astype(float)
    fit = fit_tails(beside, EDGES, gauss_fit(beside), shape=shape)
    assert fit.failed
    assert fit_tails(0 * broad, EDGES, gauss_fit(broad), shape=shape).tails is None  # nothing to fit


def test_the_simulated_tails_are_those_of_the_shape(shape):
    counts = np.random.default_rng(11).poisson(thick_target_peak(300)).astype(float)
    fit = fit_tails(counts, EDGES, gauss_fit(counts), shape=shape)
    t = fit.tails
    assert t["simulated"] and t["share"] == shape["share"]
    assert t["wide_sigma"] / t["sigma"] == pytest.approx(shape["ratio"])
    assert (t["wide_centre"] - t["centre"]) / t["sigma"] == pytest.approx(shape["offset"])
    # The width that is quoted is the Gaussian's of the simulated peak, stretched as the core is.
    assert fit.sigma / t["simulated_sigma"] == pytest.approx(t["sigma"] / shape["core_sigma"])
    assert fit.shape == "Gaussian with simulated tails on a line with a step"
    assert PeakFit.of(dict(json.loads(json.dumps(asdict(fit))), later="unknown")) == fit


# -- decays inside the target ---------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def one_sided() -> tuple:
    """A spread (MeV) with a tail of 15 % on one side, −20 to 0 keV: what is left of the shift when some nuclei decay
    inside the target."""
    at = np.linspace(-20, 4, 241)
    share = 0.85 * np.exp(-0.5 * (at / 1.4) ** 2) / np.exp(-0.5 * (at / 1.4) ** 2).sum()
    tail = np.where(at < 0, np.exp(at / 6) * (at + 20), 0.0)
    return 1e-3 * at, share + 0.15 * tail / tail.sum()


def spread_peak(n: float, spread: tuple, sigma: float = 0.9) -> np.ndarray:
    shape = sum(p * gaussian(sigma, 1e3 * d) for d, p in zip(*spread))
    return n * shape + 2e-3 * 0.2 * n * np.maximum(1 - np.cumsum(shape) + shape / 2, 0.0)


def test_the_spread_of_the_setup_gives_the_whole_peak(one_sided):
    n = 1_000_000
    counts = spread_peak(n, one_sided)
    plain = fit_tails(counts, EDGES, gauss_fit(counts))
    fit = fit_tails(counts, EDGES, gauss_fit(counts), spread=one_sided)
    t, e = fit.tails, 1e-3 * X + E0
    assert fit.area == pytest.approx(n, rel=3e-4) and abs(plain.area / n - 1) > 2e-3
    assert t["share"] == 0 and 1e3 * t["sigma"] == pytest.approx(0.9, abs=0.02) and not t["simulated"]
    assert t["centre"] == pytest.approx(E0, abs=2e-5)
    assert sum(t["spread"]["share"]) == pytest.approx(1) and fit.shape.startswith("Gaussian spread")
    # The curve is the fit's, and it survives being kept.
    assert float(np.sum(fit.curve(e) - fit.background_curve(e))) == pytest.approx(fit.area, rel=2e-3)
    assert fit.curve(e)[np.abs(X) < 20] == pytest.approx(counts[np.abs(X) < 20], abs=1e-3 * counts.max())
    assert PeakFit.of(json.loads(json.dumps(asdict(fit)))) == fit
    # Without a spread nothing changes, and a Gaussian peak stays one whatever the spread.
    assert fit_tails(counts, EDGES, gauss_fit(counts), spread=None) == plain
    few = np.random.default_rng(5).poisson(60 * gaussian(2.9) + 0.02).astype(float)
    assert fit_tails(few, EDGES, gauss_fit(few), spread=one_sided).tails is None


def test_few_counts_take_the_spread_unless_the_gaussian_is_clearly_better(one_sided):
    """With a negative worth, as a weighted sample has it: too few counts to show the spread take it all the same and
    keep the area (their Gaussian loses the tail), and many counts that are no peak of that shape do not."""
    rng = np.random.default_rng(8)
    few = [rng.poisson(spread_peak(140, one_sided)).astype(float) for _ in range(60)]
    toys = [fit_tails(c, EDGES, gauss_fit(c), spread=one_sided, worth=-TAILS_WORTH) for c in few]
    assert all("spread" in f.tails for f in toys)
    area = np.array([f.area for f in toys]) / 140
    assert abs(area.mean() - 1) < 3 * area.std(ddof=1) / math.sqrt(len(area))
    assert sum(fit_tails(c, EDGES, gauss_fit(c), spread=one_sided).tails is None for c in few) > 30
    wide = rng.poisson(100_000 * gaussian(2.9) + 2.0).astype(float)
    assert fit_tails(wide, EDGES, gauss_fit(wide), spread=one_sided, worth=-TAILS_WORTH).tails is None


def test_few_counts_take_the_spread_with_the_shape(one_sided):
    """The shape holds the spread, and the counts stretch it with the core."""
    entries = np.random.default_rng(6).poisson(spread_peak(15_000, one_sided)).astype(float)
    shape = tail_shape(fit_tails(entries, EDGES, gauss_fit(entries), spread=one_sided))
    rng = np.random.default_rng(7)
    toys = [fit_tails(c, EDGES, gauss_fit(c), shape=shape)
            for c in (rng.poisson(spread_peak(140, one_sided)).astype(float) for _ in range(60))]
    assert all(f.failed is None and f.tails["simulated"] for f in toys)
    area = np.array([f.area for f in toys]) / 140
    assert abs(area.mean() - 1) < 3 * area.std(ddof=1) / math.sqrt(len(area))
    t = toys[0].tails
    wide = np.ptp(t["spread"]["shift"]) / np.ptp(shape["spread"][0])
    assert wide == pytest.approx(t["sigma"] / shape["core_sigma"], rel=0.03)


# -- another line beside the peak ------------------------------------------------------------------------------------


def blended_peak(n: float, ratio: float = 0.025, at: float = -9.36) -> np.ndarray:
    """The peak of ``n`` counts as the clovers of ``coulex_ni62_pt194_clovers`` see it (σ 4.5 keV, with tails that end),
    and a line of the same shape with ``ratio`` of its counts 9.36 keV below (4⁺ → 2⁺ of ⁶²Ni)."""
    def one(centre):
        z = (X - centre) / 3.0
        return 0.8 * gaussian(4.0, centre) + 0.2 * _biweight(z / 3.2)[0] / 3.0 / 3.2 * 0.2

    both = one(0.0) + ratio * one(at)
    return n * both + 4e-4 * n * 0.2 * np.ones(len(X))


BESIDE = ({"shift": -9.36e-3, "ratio": 0.025, "transition": (4, 1), "counted": True},)


def test_the_peak_and_the_line_beside_it_keep_their_area():
    """With the line in the fit, at its ratio to the peak, the area is that of both. (The Gaussian on a line, which is
    what a peak with another line in the window got, counts them low.)"""
    mu = blended_peak(1e6)
    fit = gauss_fit(mu, mu)
    with_line = fit_tails(mu, EDGES, fit, variance=mu, lines=BESIDE)
    assert with_line.area == pytest.approx(1.025e6, rel=0.004)
    assert fit.area / 1.025e6 - 1 < -0.01
    ln = with_line.tails["lines"]
    assert [x["counted"] for x in ln] == [True, True] and ln[1]["transition"] == [4, 1]
    assert sum(x["share"] for x in ln) == pytest.approx(1.0)
    assert with_line.tails["inside"] < 1 and json.loads(json.dumps(asdict(with_line)))["tails"]["lines"] == ln
    # The curve on its background is the fitted counts, whichever way the fit is drawn or kept.
    kept = PeakFit.of(json.loads(json.dumps(asdict(with_line))))
    x = (EDGES[:-1] + EDGES[1:]) / 2
    m = np.abs(x - E0) < 0.04
    assert kept.curve(x)[m].sum() == pytest.approx(mu[m].sum(), rel=0.003)
    assert np.sum((kept.curve(x) - kept.background_curve(x))[m]) == pytest.approx(with_line.area, rel=0.01)
    assert "another line" in with_line.shape


def test_a_line_outside_the_window_is_in_the_fit_but_not_in_the_area():
    """A line 40 keV (9 σ) away is in the fit, at the edge of its window, and not in the area."""
    mu = blended_peak(1e6, ratio=0.05, at=-40.0)
    apart = ({"shift": -40e-3, "ratio": 0.05, "transition": (3, 1), "counted": False},)
    fit = gauss_fit(mu, mu)
    with_line = fit_tails(mu, EDGES, fit, variance=mu, lines=apart)
    assert with_line.area == pytest.approx(1e6, rel=0.006)
    assert [x["counted"] for x in with_line.tails["lines"]] == [True, False]


def test_a_gaussian_peak_with_a_line_beside_it_has_no_tails():
    mu = 1e4 * gaussian(4.0) + 300 * gaussian(4.0, -9.36)
    y = np.random.default_rng(5).poisson(mu).astype(float)
    fit = fit_tails(y, EDGES, gauss_fit(y), lines=({"shift": -9.36e-3, "ratio": 0.03, "transition": (4, 1),
                                                    "counted": True},))
    assert fit.tails["share"] == 0.0 and fit.tails["lines"][1]["transition"] == [4, 1]
    assert fit.area == pytest.approx(10_300, abs=4 * fit.area_unc) and "Gaussian (and another" in fit.shape


# -- the line's own Compton edge in the window ----------------------------------------------------------------------------

HIGH, HIGH_EDGES = 6.917, np.arange(6.571, 7.2635, 0.0006)
HIGH_SIGMA = 0.057
#: A germanium crystal's response: 16 incomplete energies for each full one, a seventh of them above the Compton edge
#: 246 keV below the line and so under the peak's foot.
GE = Response(Experiment.example("coulex_ni58"), Experiment.example("coulex_ni58").gamma_detectors[0])


def below_the_line() -> np.ndarray:
    n = int(0.9 / 0.0006)
    return GE.incomplete(HIGH, HIGH + 0.0006 * np.arange(-n, 1))


def line_and_continuum(n: float, sigma: float = HIGH_SIGMA) -> np.ndarray:
    """The counts expected of a line of ``n`` full-energy counts: the peak and what it leaves below it, both spread by
    ``sigma`` (MeV), by summing over the energies that are left."""
    below = below_the_line()
    x = (HIGH_EDGES[:-1] + HIGH_EDGES[1:]) / 2
    at = np.r_[HIGH, HIGH + 0.0006 * (np.arange(len(below)) - len(below) + 0.5)]
    each = np.r_[1.0, below]
    out = np.zeros(len(x))
    for e, c in zip(at, each):
        out += c * np.exp(-0.5 * ((x - e) / sigma) ** 2)
    return n * out * 0.0006 / (sigma * math.sqrt(2 * math.pi))


def test_the_incomplete_energies_are_those_of_the_response():
    edges = np.arange(0.0, 7.4, 0.001)
    r = GE
    counts = r.incomplete(HIGH, edges)
    pt = float(r.peak_to_total(HIGH))
    assert counts.sum() == pytest.approx((1 - pt) / pt, rel=1e-6)
    whole = r.shape(HIGH, edges)  # with the resolution, and the peak
    low = edges[1:] < 6.60
    assert counts[low].sum() * pt == pytest.approx(whole[low].sum(), rel=1e-3)
    edge = compton_edge(HIGH)
    at = lambda e: float(counts[np.abs((edges[:-1] + edges[1:]) / 2 - e) < 0.02].mean())  # noqa: E731
    assert at(edge - 0.03) > 1.1 * at(edge + 0.03) > 0 and at(HIGH + 0.05) == 0
    assert counts[np.argmin(np.abs(edges - (HIGH - ELECTRON_MEV))) - 1:][:2].max() > 10 * at(6.3)  # the single escape
    # Add-back puts the Compton deposits into the peak, a shield takes every incomplete energy away.
    exp = Experiment.example("coulex_ni62_pt194_clovers")
    clover = Response(exp, exp.gamma_detectors[0])
    pt = float(clover.peak_to_total(1.1729))
    assert clover.incomplete(1.1729, edges).sum() == pytest.approx((1 - pt) / pt, rel=1e-3)


def test_a_line_whose_compton_edge_is_in_the_window_keeps_its_area():
    """A straight line for the background reads these expected counts 13 % too large and 0.7 σ low, and in samples the fit
    often ran to its bound on the edge."""
    mu = line_and_continuum(150.0)
    below = below_the_line()
    plain = fit_peak(mu, HIGH_EDGES, HIGH, HIGH_SIGMA, 4.0, max_sigma=0.03 * HIGH)
    assert plain.failed or (plain.area > 1.1 * 150 and plain.centroid < HIGH - 0.03)
    for start in (0.003, 0.02, HIGH_SIGMA, 0.15):  # and the first pass finds it from any width
        fit = fit_peak(mu, HIGH_EDGES, HIGH, start, 0.35 / (2.3548 * start) + 1, max_sigma=0.03 * HIGH, below=below)
        assert fit.failed is None and fit.area == pytest.approx(150, rel=2e-3)
        assert abs(fit.centroid - HIGH) < 2e-4 and fit.sigma == pytest.approx(HIGH_SIGMA, rel=5e-3)
    x = (HIGH_EDGES[:-1] + HIGH_EDGES[1:]) / 2
    assert np.allclose(fit.curve(x), mu, atol=2e-3 * mu.max())  # what is drawn
    assert fit.background_curve(x).sum() == pytest.approx((mu - line_and_continuum_peak(150.0)).sum(), rel=3e-3)
    assert float(fit.curve(HIGH)) == pytest.approx(float(np.interp(HIGH, x, mu)), rel=5e-3)
    kept = json.loads(json.dumps(asdict(fit)))
    assert PeakFit.of(kept) == fit
    kept.pop("below")
    assert PeakFit.of(kept).below is None and PeakFit.of(kept).shape == "Gaussian on a line"


def line_and_continuum_peak(n: float) -> np.ndarray:
    x = (HIGH_EDGES[:-1] + HIGH_EDGES[1:]) / 2
    return n * np.exp(-0.5 * ((x - HIGH) / HIGH_SIGMA) ** 2) * 0.0006 / (HIGH_SIGMA * math.sqrt(2 * math.pi))


@pytest.mark.parametrize("weight", [1.0, 0.25])
def test_few_counts_on_the_edge_are_not_biased_and_scatter_as_printed(weight):
    """150 counts in the peak, each once (a run) or from entries that stand for a quarter of one (a sample)."""
    rng = np.random.default_rng(31)
    mu, below = line_and_continuum(150.0), below_the_line()
    areas, uncs, centres = [], [], []
    for _ in range(60):
        counts = weight * rng.poisson(mu / weight)
        fit = fit_peak(counts, HIGH_EDGES, HIGH, HIGH_SIGMA, 4.0, max_sigma=0.03 * HIGH, below=below,
                       variance=None if weight == 1.0 else weight * counts)
        assert fit.failed is None
        areas.append(fit.area)
        uncs.append(fit.area_unc)  # of as many counts
        centres.append(fit.centroid)
    scatter = float(np.std(areas, ddof=1))
    assert abs(np.mean(areas) / 150 - 1) < 3 * scatter / 150 / math.sqrt(60)
    assert 0.7 < scatter / (math.sqrt(weight) * np.mean(uncs)) < 1.4
    assert abs(np.mean(centres) - HIGH) < 3 * np.std(centres) / math.sqrt(60)


def test_a_peak_with_tails_on_its_own_edge_keeps_its_area():
    """The tails of the peak (fit_tails) and the continuum spread as they are."""
    wide = 0.85 * line_and_continuum(20_000.0, 0.045) + 0.15 * line_and_continuum(20_000.0, 0.09)
    fit = fit_peak(wide, HIGH_EDGES, HIGH, 0.05, 4.0, max_sigma=0.03 * HIGH, below=below_the_line())
    full = fit_tails(wide, HIGH_EDGES, fit)
    assert full.tails is not None and full.below is not None and full.tails["step"] == 0
    # The wider component is a bell that ends, and this one is a Gaussian (the area of the Gaussian alone is 2.6 % off).
    assert full.area == pytest.approx(20_000, rel=0.01) and abs(full.centroid - HIGH) < 1e-3
    assert abs(full.area / 20_000 - 1) < 0.5 * abs(fit.area / 20_000 - 1)
    x = (HIGH_EDGES[:-1] + HIGH_EDGES[1:]) / 2
    window = np.abs(x - fit.centroid) < 4 * 2.3548 * fit.sigma
    assert full.curve(x[window]).sum() == pytest.approx(wide[window].sum(), rel=1e-3)
    assert full.below["area"] == full.area and full.shape.endswith("the line's own Compton continuum")
    # With the shape of that peak, a few counts.
    few = np.random.default_rng(8).poisson(wide / 100).astype(float)
    got = fit_tails(few, HIGH_EDGES, fit_peak(few, HIGH_EDGES, HIGH, 0.05, 4.0, max_sigma=0.03 * HIGH,
                                              below=below_the_line()), shape=tail_shape(full))
    assert got.failed is None and got.tails["simulated"] and got.area == pytest.approx(200, rel=0.2)


def test_another_background_under_a_line_on_its_edge_is_the_straight_lines():
    """Twice as many counts in the window that are not the line's as that are: the area is still not biased, and scatters
    as printed. And no counts at all are a failed fit."""
    rng = np.random.default_rng(17)
    below = below_the_line()
    level = 2 * line_and_continuum(150.0).sum() / (len(HIGH_EDGES) - 1)
    areas, uncs = [], []
    for _ in range(60):
        counts = rng.poisson(level + line_and_continuum(150.0)).astype(float)
        fit = fit_peak(counts, HIGH_EDGES, HIGH, HIGH_SIGMA, 4.0, max_sigma=0.03 * HIGH, below=below)
        assert fit.failed is None
        areas.append(fit.area)
        uncs.append(fit.area_unc)
    scatter = float(np.std(areas, ddof=1))
    assert abs(np.mean(areas) / 150 - 1) < 3 * scatter / 150 / math.sqrt(60)
    assert 0.7 < scatter / np.mean(uncs) < 1.4
    nothing = np.zeros(len(HIGH_EDGES) - 1)
    fit = fit_peak(nothing, HIGH_EDGES, HIGH, HIGH_SIGMA, 4.0, below=below)
    assert fit.failed.startswith("there are no counts above")
