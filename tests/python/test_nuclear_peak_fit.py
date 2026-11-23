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

from physim.nuclear.analysis import SHAPE_REACH, PeakFit, fit_peak, fit_tails, tail_shape

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
