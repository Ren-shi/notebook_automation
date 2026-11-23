"""Automatic analysis of the simulated experiment: from the γ-ray peak areas to B(E2), with uncertainties, the
shape, and what the planned beam time determines.

The steps are those of a Coulomb-excitation analysis, each shown and adjustable::

    from physim.nuclear import Experiment
    from physim.nuclear.analysis import Analysis, Settings
    from physim.nuclear.gamma_events import simulate_gammas

    exp = Experiment.example("coulex_ni58")
    a = Analysis(exp, simulate_gammas(exp, 400_000, seed=1))
    result = a.run()                                   # the default settings
    result.b_e2fm4, result.statistical, result.systematic
    result.budget                                      # each systematic contribution
    a.run(Settings(rings=range(4, 16), fit_half_width=3.0))   # changed settings; the change is recorded

1. **Particle gate:** the particle must sit in the inelastic group of its detector: its measured energy nearer
   the kinematic energy of a particle that excited the state than that of an elastic one, and within
   ``gate_width`` of the first. Rings (or strips) can be chosen. When the two groups overlap through the depth of
   the target, the result warns and says to use the gate "all" (:attr:`Result.gate_warning`).
2. **Doppler correction** for the emitting nucleus (or the other one, to see what happens).
3. **Peak search:** each corrected spectrum (for the projectile, and for the recoil) is searched for peaks over
   its whole range (:func:`find_peaks`), each peak found is fitted, and those that sit on a transition of the
   nucleus that correction is right for are named (:func:`identify`). The other nucleus's lines are smeared by
   the wrong correction: each nucleus shows its peaks in its own spectrum.
4. **Peak fit:** a Gaussian on a straight line, by maximum likelihood within ``fit_half_width`` resolutions
   of the transition energy, with the random coincidences expected in the window subtracted. Where the peak has
   tails that a Gaussian leaves out (a thick target: the emitter's speed, and so the Doppler shift, depends on how
   deep the reaction was), the area is that of a peak with them (:func:`fit_tails`). A run of too few counts to fix
   the tails takes their share, width and offset from a simulation of the same setup
   (:meth:`Analysis.simulated_peak`), and the budget says what that is worth. A fit that ends at a bound is reported
   as failed (:attr:`PeakFit.failed`), with no area.
5. **Yield:** the area over the full-energy-peak efficiency (each crystal's at the Doppler-shifted energy it sees
   the line at), the γ-ray branch (internal conversion) and the angular-correlation factor.
6. **Normalisation**, chosen by the user: ``"rutherford"``, to the particles counted in the same rings,
   so the beam current and the target thickness cancel; or ``"target"``, to a transition of the other nucleus of
   the collision whose strength is known. When the events hold that nucleus's γ rays (the setup has its level
   scheme), its peak is found, fitted and integrated in its own corrected spectrum like the first one, and its
   yield is compared with what its level scheme gives to all orders (:attr:`Result.lines`); otherwise the
   reference is ``Settings.reference``, its peak taken from its analytic rate with Poisson noise (noted in the
   result).
7. **B(E2)** from the first-order proportionality of the excitation probability to B(E2), the excitations
   expected per particle taken from the forward model's cross sections through the depth of the target and with
   the excited particle's kinematics (:meth:`Analysis._mean_probability`), in e²fm⁴, e²b² and
   Weisskopf units, with the statistical uncertainty of the fit and a budget of systematic contributions, each
   found by changing its input by one standard deviation and running again. When the events follow a level
   scheme to all orders, first order is only the start: the scheme's coupled equations are solved with the
   matrix element to the state as the unknown, until they give the measured γ rays per particle, every
   transition inside the peak's window counted (the feeding of the state from above, and a line that the peak
   blends with). ``Settings.order = "first"`` keeps the first-order value.
8. **Shape:** β₂, the strength in Weisskopf units, E(4⁺)/E(2⁺) if a level scheme is there, and the quadrupole
   moment a rigid rotor would have, with the warning that these readings depend on the rotor model.

The result also holds the comparison with the value that was put in, and the beam time for a chosen precision.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, fields, replace
from typing import Optional

import numpy as np

from . import data
from .coulex import Coulex
from .gamma import correlation_table, emitters, transition_lines
from .kinematics import TwoBody
from .levels import weisskopf_unit
from .quantity import Quantity
from .rates import (Rates, SchemeProbability, _q, beam_energy_at, beam_ion, channels, coulex_for, depth_probability,
                    exit_energy, reaction_at, stack, stopping, tilt_deg)
from .response import FWHM_PER_SIGMA, Response

#: Hours in a shift, for the counts per shift.
SHIFT_H = 8.0
#: The expected excitations per particle (:meth:`Analysis._mean_probability`) are tabulated against the particle's
#: lab angle in steps of this many degrees, from the forward cross sections integrated over the depth of each layer
#: with this many Gauss–Legendre points.
FORWARD_STEP_DEG = 0.05
FORWARD_DEPTH_POINTS = 24
#: The elastic and inelastic groups of a particle detector overlap at a lab direction when the energies of the
#: elastic group, from every depth of the target, cover more than this share of the inelastic group's
#: (:meth:`Analysis._groups_overlap`); the inelastic particle gate is called unreliable when they overlap for more
#: than ``GATE_OVERLAP_WARN`` of the excitations in the chosen rings.
GATE_GROUP_COVER = 0.5
GATE_OVERLAP_WARN = 0.5
#: A peak is taken to show its tails (:func:`fit_tails`) when they improve the fit's deviance by more than this: the
#: three parameters they add do that by chance about once in 60 000 fits. A run whose peak has too few counts to
#: reach it takes the tails of a simulated one (:meth:`Analysis.simulated_peak`).
TAILS_WORTH = 25.0
#: With a simulated tail shape, how many times wider or narrower than the simulated peak's core the fitted core may
#: be, and how many of its σ its centre may lie away from the simulated one.
SHAPE_REACH = 2.0
#: The simulated peak: its seed (one more if the run has that seed), the reactions generated, the γ rays each of its
#: particles emits, the full-energy counts it is made for and the most γ rays emitted for them.
SHAPE_SEED = 1_618_033
SHAPE_EVENTS = 600_000
SHAPE_COPIES = 10
SHAPE_COUNTS = 15_000
SHAPE_EMISSIONS = 1_000_000
#: A simulated tail share is never taken to be known better than this part of what it changes the area by (nor
#: worse than all of it): how well the simulated width matches the measured one says where in between.
SHAPE_KNOWN = 0.1


def _with_beam_energy(experiment, energy_mev: float):
    """The setup with another beam energy (nothing else of it changes)."""
    return replace(experiment, beam=replace(experiment.beam, energy=Quantity(energy_mev, "MeV")))


class _Zero:
    """P = 0 at every angle: below the excitation threshold."""

    def probability(self, theta_cm):
        return 0.0 if np.ndim(theta_cm) == 0 else np.zeros(np.shape(theta_cm))


@dataclass
class Settings:
    """The adjustable choices of the analysis."""

    #: Which particle detectors to use (all if None) and which rings or strips of them (all if None).
    detectors: Optional[list] = None
    rings: Optional[list] = None
    #: "inelastic" keeps particles in the inelastic group, "all" takes every detected particle.
    particle_gate: str = "inelastic"
    #: A window on the particle's measured energy (MeV), as drawn on the Data tab's energy-against-ring view.
    particle_energy: Optional[tuple] = None
    #: Half-width of the inelastic gate, in resolutions (FWHM) of the particle detector.
    gate_width: float = 3.0
    #: Doppler correction: "emitter" (the nucleus that was excited), "projectile" or "recoil".
    correction: str = "emitter"
    #: Half-width of the fit window, in resolutions (FWHM) of the γ-ray peak.
    fit_half_width: float = 4.0
    subtract_randoms: bool = True
    #: "rutherford" or "target".
    normalisation: str = "rutherford"
    #: For "target": the reference transition, {"energy": "328 keV", "b_up": "1.65 e2b2", "unc": 0.03}, its
    #: B(E2↑) known to the relative uncertainty ``unc``. The reference is in the other nucleus of the collision.
    reference: Optional[dict] = None
    #: The target's offset along the beam assumed in the Doppler correction, mm from its true place (0: the
    #: analysis knows the geometry); and the uncertainty of the target's place, mm, for the budget.
    target_offset_mm: float = 0.0
    #: Relative systematic uncertainties of the inputs.
    efficiency_unc: float = 0.05
    beam_energy_unc: float = 0.005
    position_unc_mm: float = 1.0
    correlation_unc: float = 0.03
    #: The relative statistical uncertainty the beam time is asked for.
    wanted_precision: float = 0.05
    #: How B(E2) follows from the yield: "auto" (to all orders when the events follow a level scheme to all
    #: orders, else first order) or "first".
    order: str = "auto"

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}


@dataclass
class PeakFit:
    """A Gaussian on a straight line fitted to a spectrum."""

    area: float
    area_unc: float
    centroid: float
    sigma: float
    background: float
    chi2: float
    bins: int
    randoms: float = 0.0
    #: The area of the fitted Gaussian alone, for comparison with the summed area.
    gauss_area: float = 0.0
    #: The fitted curve, for drawing it: height × exp(−(E − centroid)²/2σ²) + level + slope × (E − origin),
    #: in counts per bin with E in MeV.
    height: float = 0.0
    level: float = 0.0
    slope: float = 0.0
    origin: float = 0.0
    #: The peak with tails that the area comes from (:func:`fit_tails`), or None for the Gaussian alone: the
    #: "centre" and "sigma" (MeV) of its core, the "share" of its area in the wider Gaussian, the "wide_centre" and
    #: "wide_sigma" of that one, the total "area" (counts) and the bin "width" (MeV); the "level" and "slope" (about
    #: "origin") of the line above it and the "step" the background takes below it (counts per bin); "inside", the
    #: share of the peak within ±3σ of the Gaussian; "gain", what the tails took off the fit's deviance; and whether
    #: they are "simulated" (the share, the offset and the width of the wider Gaussian, relative to the core, are
    #: those of a simulated peak: :meth:`Analysis.simulated_peak`). With simulated tails also the σ of the Gaussian
    #: fitted to that peak ("simulated_sigma", MeV), the share of it that the area of that Gaussian loses
    #: ("simulated_loss") and the area of the Gaussian fitted to these counts ("gaussian_area", or None).
    tails: Optional[dict] = None
    #: Why the fit is no measurement of a peak (its width or its centre ended at a bound, or it found no counts above
    #: the background), or None.
    failed: Optional[str] = None

    @classmethod
    def of(cls, d: dict) -> PeakFit:
        """The fit kept as a dict (``dataclasses.asdict``), whatever version of physim kept it."""
        return cls(**{k: v for k, v in d.items() if k in {f.name for f in fields(cls)}})

    def _components(self) -> tuple:
        """(share of the area, centre, σ) of the two Gaussians of a peak with tails."""
        t = self.tails
        return (1 - t["share"], t["centre"], t["sigma"]), (t["share"], t["wide_centre"], t["wide_sigma"])

    def background_curve(self, energy) -> np.ndarray:
        """The fitted background under the peak at ``energy`` (MeV), counts per bin."""
        e = np.asarray(energy, dtype=float)
        t = self.tails
        if t is None:
            return self.level + self.slope * (e - self.origin)
        above = sum(share * 0.5 * np.vectorize(math.erfc, otypes=[float])((e - centre) / sigma / math.sqrt(2))
                    for share, centre, sigma in self._components())
        return t["level"] + t["slope"] * (e - t["origin"]) + t["step"] * above

    def curve(self, energy) -> np.ndarray:
        """The fitted peak on its background at ``energy`` (MeV), counts per bin."""
        e = np.asarray(energy, dtype=float)
        t = self.tails
        if t is None:
            return self.height * np.exp(-0.5 * ((e - self.centroid) / self.sigma) ** 2) + self.background_curve(e)
        peak = sum(share * np.exp(-0.5 * ((e - centre) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
                   for share, centre, sigma in self._components())
        return t["area"] * t["width"] * peak + self.background_curve(e)

    @property
    def fwhm_kev(self) -> float:
        return 1e3 * FWHM_PER_SIGMA * self.sigma

    @property
    def shape(self) -> str:
        """What was fitted, in words."""
        if self.tails is None:
            return "Gaussian on a line"
        return f"Gaussian with {'simulated ' if self.tails.get('simulated') else ''}tails on a line with a step"


def fit_peak(counts: np.ndarray, edges: np.ndarray, energy: float, sigma0: float, half_width: float = 4.0,
             iterations: int = 40, max_sigma: Optional[float] = None, variance=None) -> PeakFit:
    """Fit a Gaussian on a straight line to ``counts`` within ``half_width`` × FWHM of ``energy``, by maximum
    likelihood (the counts of each bin are Poisson; the line is nowhere negative), then count the peak within ±3σ of
    its centre over the line. Returns the area in counts. A fit whose width or centre ends at the bound it may not
    pass (``max_sigma`` and 0.3 of ``sigma0``, ±2 ``sigma0``), or which finds no counts above the line, is marked
    :attr:`PeakFit.failed`: it describes no peak, and its area means nothing.

    ``variance`` is the variance of each bin when it is not its counts: a spectrum of weighted events scaled to
    a beam time holds far more counts than it has entries, and a fit that took the counts for their own variance
    would trust every bin too much and follow the scatter of the few entries. The fit is then least squares with
    these variances (Gauss–Newton). Only the fit's weights use it; the uncertainty of the area stays that of the
    counts, which is what the beam time gives."""
    centres = (edges[:-1] + edges[1:]) / 2
    width = float(edges[1] - edges[0])
    half = half_width * FWHM_PER_SIGMA * sigma0
    m = np.abs(centres - energy) <= half
    x, y = centres[m] - energy, counts[m].astype(float)
    if m.sum() < 6:
        raise ValueError("the fit window holds too few bins; widen it or use finer bins")
    widest = max_sigma if max_sigma else 4.0 * sigma0
    if variance is None:
        p, chi2 = _gaussian(x, y, sigma0, width, widest)
    else:
        var = np.asarray(variance, dtype=float)[m]
        # A bin without entries is as uncertain as one entry.
        filled = var[var > 0]
        w = 1 / np.maximum(var, float(np.median(filled)) if filled.size else 1.0)
        p, chi2 = _gaussian_weighted(x, y, w, sigma0, widest, iterations)
    failed = None
    gauss_area = p[0] * p[2] * math.sqrt(2 * math.pi) / width
    if not (p[0] > 0 if variance is not None else gauss_area > 1):
        failed = "there are no counts above the background where the peak is looked for"
    else:
        fwhm = 1e3 * FWHM_PER_SIGMA * p[2]
        if p[2] >= widest * (1 - 1e-9):
            failed = f"the peak's width ran to the largest the fit allows (FWHM {fwhm:.1f} keV)"
        elif p[2] <= 0.3 * sigma0 * (1 + 1e-9):
            failed = f"the peak's width ran to the smallest the fit allows (FWHM {fwhm:.1f} keV)"
        elif abs(p[1]) >= 2 * sigma0 * (1 - 1e-9):
            failed = f"the peak's centre ran to the edge of where it is looked for ({1e3 * (energy + p[1]):.1f} keV)"
    # The area is the counts within ±3σ of the fitted centre less the fitted line under them: a peak that is
    # not quite Gaussian (the Doppler correction leaves a flat-topped one) is counted all the same. The 0.9973
    # restores the Gaussian tails outside ±3σ.
    inside = np.abs(x - p[1]) < 3 * p[2]
    under = max(float(np.sum(p[3] + p[4] * x[inside])), 0.0)
    area = (float(y[inside].sum()) - under) / 0.9973
    # Counting statistics of the window, plus the background's own from the bands the line was fitted to.
    outside = ~inside
    bg_var = under * under / max(float(y[outside].sum()), 1.0) if outside.any() else under
    area_unc = math.sqrt(max(float(y[inside].sum()), 0.0) + bg_var)
    return PeakFit(float(area), float(area_unc), float(energy + p[1]), float(p[2]), under, chi2, int(m.sum()),
                   gauss_area=float(gauss_area), height=float(p[0]), level=float(p[3]), slope=float(p[4]),
                   origin=float(energy), failed=failed)


def _gaussian(x: np.ndarray, y: np.ndarray, sigma0: float, width: float, widest: float) -> tuple:
    """(height, centre, σ, level and slope of the line; deviance) of the likeliest Gaussian on a line that is nowhere
    negative to the Poisson counts ``y`` of bins at ``x`` (MeV, from where the peak is expected)."""
    z, b = x / sigma0, width / sigma0
    u = (z - z[0]) / (z[-1] - z[0])

    def model(p):
        return p[0] * b / (math.sqrt(2 * math.pi) * p[2]) * np.exp(-0.5 * ((z - p[1]) / p[2]) ** 2) \
            + p[3] * (1 - u) + p[4] * u

    # Area (counts), centre, σ (in σ₀), and the line at the window's two ends.
    y = np.maximum(y, 0.0)
    side = max(len(y) // 8, 2)
    ends = [float(y[:side].mean()), float(y[-side:].mean())]
    start = np.array([max(float(y.sum()) - 0.5 * len(y) * sum(ends), 1.0), 0.0, 1.0, *ends])
    low = np.array([0.0, -2.0, 0.3, 0.0, 0.0])
    high = np.array([np.inf, 2.0, widest / sigma0, np.inf, np.inf])
    p, deviance, _ = _likeliest(model, y, 1.0, start, low, high, np.ones(5, dtype=bool))
    slope = (p[4] - p[3]) / (x[-1] - x[0])
    return np.array([p[0] * b / (math.sqrt(2 * math.pi) * p[2]), p[1] * sigma0, p[2] * sigma0,
                     p[3] - slope * x[0], slope]), deviance


def _gaussian_weighted(x: np.ndarray, y: np.ndarray, w: np.ndarray, sigma0: float, widest: float,
                       iterations: int) -> tuple:
    """The same parameters by least squares with the weights ``w`` (Gauss–Newton), and χ²."""
    # Parameters: height, mean (relative to energy), sigma, background level and slope.
    p = np.array([max(float(y.max() - np.median(y)), 1.0), 0.0, sigma0, float(np.median(y)), 0.0])

    def model(q):
        return q[0] * np.exp(-0.5 * ((x - q[1]) / q[2]) ** 2) + q[3] + q[4] * x

    for _ in range(iterations):
        r = y - model(p)
        jac = np.empty((len(x), 5))
        g = np.exp(-0.5 * ((x - p[1]) / p[2]) ** 2)
        jac[:, 0] = g
        jac[:, 1] = p[0] * g * (x - p[1]) / p[2] ** 2
        jac[:, 2] = p[0] * g * (x - p[1]) ** 2 / p[2] ** 3
        jac[:, 3] = 1.0
        jac[:, 4] = x
        a = jac.T @ (jac * w[:, None])
        b = jac.T @ (w * r)
        try:
            step = np.linalg.solve(a + 1e-12 * np.eye(5), b)
        except np.linalg.LinAlgError:
            break
        p = p + step
        # Keep the peak a peak: its width between a half and three times the resolution, its centre in the
        # window, its height positive.
        p[2] = min(max(abs(p[2]), 0.3 * sigma0), widest)
        p[1] = min(max(p[1], -2 * sigma0), 2 * sigma0)
        p[0] = max(p[0], 0.0)
        if np.max(np.abs(step[:3]) / (np.abs(p[:3]) + 1e-12)) < 1e-7:
            break
    r = y - model(p)
    return p, float(np.sum(w * r * r))


def _likeliest(model, y: np.ndarray, scale: float, start: np.ndarray, low: np.ndarray, high: np.ndarray,
               free: np.ndarray) -> tuple:
    """(parameters, deviance, covariance) of the bounded fit of ``model(parameters)`` to the counts ``y`` of bins whose
    variance is ``scale`` × their expected counts, by maximum Poisson likelihood. (Least squares with the variances of
    the counts would count bins of few entries low.) Only the ``free`` parameters move; Fisher scoring with
    Levenberg–Marquardt damping. The covariance is that of as many counts (not of as many entries), zero for a
    parameter that does not move or ends at its bound."""
    p = np.clip(np.asarray(start, dtype=float), low, high)
    idx = np.flatnonzero(free)
    floor = 1e-4 * scale  # no bin expects less than that share of an entry
    some = y > 0
    const = float(np.sum(y[some] * np.log(y[some]) - y[some])) if some.any() else 0.0

    def expect(q):
        return np.maximum(model(q), floor)

    def deviance(e):
        return 2 * (float(e.sum() - np.sum(y[some] * np.log(e[some]))) + const) / scale

    def jacobian(p):
        jac = np.empty((len(y), len(idx)))
        base = model(p)
        for j, k in enumerate(idx):
            h = 1e-5 * max(abs(p[k]), 1e-2)
            h = -h if p[k] + h > high[k] else h
            q = p.copy()
            q[k] += h
            jac[:, j] = (model(q) - base) / h
        return jac

    e = expect(p)
    cost, lam = deviance(e), 1e-3
    for _ in range(200):
        jac = jacobian(p)
        a, b = jac.T @ (jac / e[:, None]), jac.T @ ((y - e) / e)
        # A parameter at its bound that would leave it stays there.
        go = ~(((p[idx] <= low[idx]) & (b < 0)) | ((p[idx] >= high[idx]) & (b > 0)))
        if not go.any():
            break
        better = False
        for _ in range(12):
            m = a[np.ix_(go, go)] + lam * np.diag(np.diag(a)[go] + 1e-30)
            q = p.copy()
            q[idx[go]] += np.linalg.lstsq(m, b[go], rcond=None)[0]
            q = np.clip(q, low, high)
            eq = expect(q)
            dq = deviance(eq)
            if dq <= cost:
                better = True
                break
            lam *= 5
        if not better:
            break
        gained, moved = cost - dq, float(np.max(np.abs(q - p) / (np.abs(p) + 1e-3)))
        p, e, cost, lam = q, eq, dq, max(lam / 5, 1e-9)
        if gained < 1e-7 and moved < 1e-4:
            break
    jac = jacobian(p)
    a = jac.T @ (jac / e[:, None])
    off = (p[idx] > low[idx]) & (p[idx] < high[idx]) & (np.diag(a) > 1e-20)
    cov = np.zeros((len(p), len(p)))
    if off.any():
        a = a[np.ix_(off, off)]
        d = np.sqrt(np.diag(a))
        k = idx[off]
        cov[np.ix_(k, k)] = np.linalg.pinv(a / np.outer(d, d)) / np.outer(d, d)
    return p, cost, cov


def fit_tails(counts: np.ndarray, edges: np.ndarray, fit: PeakFit, half_width: float = 4.0, variance=None,
              shape: Optional[dict] = None, worth: float = TAILS_WORTH) -> PeakFit:
    """The peak of ``fit`` (:func:`fit_peak`) with the tails it has beyond a Gaussian, or ``fit`` itself if it has none
    worth fitting.

    A Doppler-corrected peak is the sum of what every crystal, particle and depth of the target leave of the shift:
    a core with wider tails. From a thick target 2 % of it lies beyond ±3σ of the Gaussian that fits its core, and
    a straight line fitted beside it takes those tails for background. Here the peak is two Gaussians (the wider one
    1.25 to 3.5 times as wide, its centre within a σ of the other's) on a straight line, nowhere negative, with a
    step: the peak's own γ rays that scatter on their way are below it, so the background steps down as the peak is
    crossed. The area is that of the fitted peak. At the maximum of the likelihood the model holds as many counts as
    the window does, so it is the counts of the whole window less the background under them. (Counting within ±3σ
    and dividing by the share of the model there takes that share from the model: a peak with a flat top, and tails
    that end, is no sum of two Gaussians, and 0.2–0.3 % was lost on one.)

    The fit is by maximum likelihood (:func:`_likeliest`), and the tails are kept when they fit better than a Gaussian
    on the same background by more than ``worth`` (:data:`TAILS_WORTH`); without them the result is ``fit``.
    ``variance`` is each bin's, if not its counts (its sum over theirs is the counts that one entry stands for).

    With a ``shape`` (:func:`tail_shape` of a peak that was fitted with its tails: one from a simulation of the setup)
    the share, the centre and the width of the wider Gaussian, relative to the core's, are those of that peak. The
    counts fix only the area, the centre and width of the core and the background, within the window of that peak.
    The core may be :data:`SHAPE_REACH` times as wide or as narrow as the shape's, and its centre as many of its σ
    away; a fit that ends there is marked :attr:`PeakFit.failed`."""
    centres = (edges[:-1] + edges[1:]) / 2
    width = float(edges[1] - edges[0])
    origin, s0 = (fit.centroid, fit.sigma) if shape is None else (shape["centroid"], shape["sigma"])
    m = np.abs(centres - origin) <= half_width * FWHM_PER_SIGMA * s0
    y = np.maximum(np.asarray(counts, dtype=float)[m], 0.0)
    if len(y) < 20 or not y.sum() > 0 or (shape is None and not fit.area > 0):
        return fit
    scale = 1.0 if variance is None else float(np.asarray(variance, dtype=float)[m].sum()) / float(y.sum())
    if not scale > 0:
        return fit
    x, b = (centres[m] - origin) / s0, width / s0  # in σ of the Gaussian
    u = (x - x[0]) / (x[-1] - x[0])

    def peak(p):  # the share of the peak in each bin
        core = np.exp(-0.5 * ((x - p[1]) / p[2]) ** 2) / p[2]
        wide = np.exp(-0.5 * ((x - p[1] - p[4] * p[2]) / (p[2] * p[5])) ** 2) / (p[2] * p[5])
        return ((1 - p[3]) * core + p[3] * wide) * b / math.sqrt(2 * math.pi)

    def background(p, share):
        return p[6] * (1 - u) + p[7] * u + p[8] * (1 - np.cumsum(share) + share / 2)

    def model(p):
        share = peak(p)
        return p[0] * share + background(p, share)

    # Area; centre and σ of the core; share, centre (from the core's, in its σ) and width (over the core's) of the
    # wider Gaussian; the line at the window's two ends; the step.
    side = max(len(y) // 8, 2)
    below, above = float(y[:side].mean()), float(y[-side:].mean())
    start = np.array([fit.area, 0.0, 0.95, 0.1, 0.0, 2.0, above, above, max(below - above, 0.0)])
    low = np.array([0.0, -1.0, 0.3, 0.0, -1.0, 1.25, 0.0, 0.0, 0.0])
    high = np.array([np.inf, 1.0, 1.5, 0.4, 1.0, 3.5, np.inf, np.inf, np.inf])
    free = np.ones(9, dtype=bool)
    alone = free.copy()
    alone[3:6] = False
    failed, gain = None, 0.0
    try:
        if shape is None:
            gaussian = start.copy()
            gaussian[3] = 0.0
            _, plain, _ = _likeliest(model, y, scale, gaussian, low, high, alone)
            p, deviance, cov = _likeliest(model, y, scale, start, low, high, free)
            gain = plain - deviance
            if not gain > worth:
                return fit
        else:
            core = shape["core_sigma"] / s0
            start[:6] = [max(float(y.sum()) - 0.5 * len(y) * (below + above), 1.0), (shape["centre"] - origin) / s0,
                         core, shape["share"], shape["offset"], shape["ratio"]]
            low[:6], high[:6] = -np.inf, np.inf
            low[0] = 0.0
            low[1], high[1] = start[1] - SHAPE_REACH * core, start[1] + SHAPE_REACH * core
            low[2], high[2] = core / SHAPE_REACH, core * SHAPE_REACH
            p, deviance, cov = _likeliest(model, y, scale, start, low, high, alone)
            if not p[0] > 1:
                failed = "there are no counts above the background where the setup puts the peak"
            elif not low[2] < p[2] < high[2]:
                failed = (f"the peak's width ran to {SHAPE_REACH:g} times "
                          f"{'more' if p[2] > core else 'less'} than the setup's resolution and Doppler spread give")
            elif not low[1] < p[1] < high[1]:
                failed = "the peak's centre ran to the edge of where the setup puts it"
    except np.linalg.LinAlgError:
        return fit
    share = peak(p)
    slope = float(p[7] - p[6]) / float(x[-1] - x[0])
    tails = {"area": float(p[0]), "width": width, "centre": origin + s0 * float(p[1]), "sigma": s0 * float(p[2]),
             "share": float(p[3]), "wide_centre": origin + s0 * float(p[1] + p[4] * p[2]),
             "wide_sigma": s0 * float(p[2] * p[5]), "level": float(p[6]) - slope * float(x[0]),
             "slope": slope / s0, "origin": origin, "step": float(p[8]), "gain": float(gain),
             "simulated": shape is not None}
    if shape is not None:  # the Gaussian of a peak of that shape: the window, and the width that is quoted
        stretch = float(p[2]) * s0 / shape["core_sigma"]
        fit = replace(fit, sigma=stretch * shape["sigma"], failed=failed,
                      centroid=tails["centre"] + stretch * (shape["centroid"] - shape["centre"]))
        tails["simulated_sigma"] = shape["sigma"]
    inside = np.abs(centres[m] - fit.centroid) < 3 * fit.sigma
    tails["inside"] = float(share[inside].sum())
    unc = math.sqrt(max(float(cov[0, 0]), float(p[0])))
    return replace(fit, area=float(p[0]), area_unc=unc, background=float(background(p, share)[inside].sum()),
                   chi2=deviance, bins=int(m.sum()), tails=tails)


def tail_shape(fit: PeakFit) -> Optional[dict]:
    """What :func:`fit_tails` needs to give another peak the tails of ``fit``: the share of the wider Gaussian, its
    centre from the core's ("offset", in the core's σ) and its σ over the core's ("ratio"), the core's "centre" and
    "core_sigma", and the "centroid" and "sigma" of the Gaussian that was fitted to it. None for a peak without
    tails."""
    t = fit.tails
    if t is None:
        return None
    return {"share": t["share"], "offset": (t["wide_centre"] - t["centre"]) / t["sigma"],
            "ratio": t["wide_sigma"] / t["sigma"], "centre": t["centre"], "core_sigma": t["sigma"],
            "centroid": fit.centroid, "sigma": fit.sigma}


def find_peaks(counts, edges, fwhm, threshold: float = 5.0, contrast: float = 0.1, variance=None,
               widths: tuple = (1.0, 2.0, 4.0, 8.0), widest: float = 0.025, most: int = 30) -> list:
    """Search a spectrum for peaks, without knowing where they are.

    At every bin, and for peaks ``widths`` times as wide as the detector's resolution (``fwhm``, a function of
    the energy in MeV that gives the FWHM in MeV; a Doppler-corrected peak is wider than the resolution, but
    never wider than ``widest`` × its energy), the counts in a window of ±2σ are compared with the level of two
    bands beside it, one on each side, a σ away. A peak is where that excess is largest, if it is at least
    ``threshold`` standard deviations and at least ``contrast`` × the level under it (a Compton edge or a bump
    left by the wrong Doppler correction is a small bend on a high level). ``variance`` is the variance of each
    bin if not its counts (a spectrum with something subtracted).

    Returns the peaks, most significant first, as dicts: "energy" (the centre of gravity of the excess, MeV),
    "sigma" (its r.m.s. width, MeV), "significance", "counts" (the excess) and "level" (under the window)."""
    y = np.asarray(counts, dtype=float)
    n = len(y)
    if n < 12:
        return []
    var = y if variance is None else np.asarray(variance, dtype=float)
    centres = (np.asarray(edges[:-1]) + np.asarray(edges[1:])) / 2
    step = float(edges[1] - edges[0])
    res = np.maximum(np.asarray(fwhm(centres), dtype=float), 2 * step)
    cy, cv = np.r_[0.0, np.cumsum(y)], np.r_[0.0, np.cumsum(np.maximum(var, 0.0))]
    idx = np.arange(n)
    best = np.zeros(n)
    best_h = np.ones(n, dtype=np.int64)
    best_g = np.ones(n, dtype=np.int64)
    for k in widths:
        f = np.minimum(k * res, np.maximum(widest * centres, res))
        h = np.maximum(np.rint(0.85 * f / step).astype(np.int64), 1)   # the window: ±2σ of such a peak
        g = np.maximum(np.rint(0.45 * f / step).astype(np.int64), 1)   # a gap of about a σ
        w = 2 * h + 1
        lo, hi = idx - h, idx + h + 1
        l0, l1, r0, r1 = lo - g - w, lo - g, hi + g, hi + g + w
        ok = (l0 >= 0) & (r1 <= n)

        def total(c, a, b):
            return c[np.clip(b, 0, n)] - c[np.clip(a, 0, n)]

        inside = total(cy, lo, hi)
        level = (total(cy, l0, l1) + total(cy, r0, r1)) / 2
        excess = inside - level
        noise = np.sqrt(np.maximum(total(cv, lo, hi) + (total(cv, l0, l1) + total(cv, r0, r1)) / 4, 1.0))
        sig = np.where(ok & (excess >= contrast * np.maximum(level, 0.0)), excess / noise, 0.0)
        better = sig > best
        best[better], best_h[better], best_g[better] = sig[better], h[better], g[better]
    out = []
    taken = np.zeros(n, dtype=bool)
    for i in np.argsort(-best, kind="stable"):
        if best[i] < threshold or len(out) >= most:
            break
        h, g = int(best_h[i]), int(best_g[i])
        lo, hi = i - h, i + h + 1
        # A window that reaches into a peak already found sees that peak again, off its centre: not another one.
        if taken[lo:hi].any():
            continue
        taken[max(lo - g, 0):hi + g] = True
        w = 2 * h + 1
        level = (y[lo - g - w:lo - g].sum() + y[hi + g:hi + g + w].sum()) / (2 * w)
        net = np.maximum(y[lo:hi] - level, 0.0)
        if net.sum() <= 0:
            continue
        mean = float(np.average(centres[lo:hi], weights=net))
        rms = math.sqrt(float(np.average((centres[lo:hi] - mean) ** 2, weights=net)))
        out.append({"energy": mean, "sigma": max(rms, step / 2), "significance": float(best[i]),
                    "counts": float((y[lo:hi] - level).sum()), "level": float(level * w)})
    return out


def identify(peaks: list, lines: list, correction: str) -> list:
    """Name the peaks of a search (:func:`find_peaks`, each with its fitted "energy" and "sigma" in MeV) in a
    spectrum corrected for ``correction`` ("projectile" or "recoil"): a peak within half its width (and 1.5 keV)
    of a transition of the nucleus that correction is right for (``lines`` from
    :func:`physim.nuclear.gamma.transition_lines`) is that transition, the nearest if several fit. A peak that
    fits none but lies within 6% of a line of the other nucleus is given as "near" it: that is where the wrong
    correction leaves its counts. Adds "label", "nuclide", "line_kev" and "near" to each peak, in place."""
    for pk in peaks:
        e_kev = 1e3 * pk["energy"]
        reach = max(0.5 * FWHM_PER_SIGMA * 1e3 * pk["sigma"], 1.5)
        own = [ln for ln in lines if ln["correction"] == correction and abs(ln["energy_kev"] - e_kev) <= reach]
        pk.update(label=None, nuclide=None, line_kev=None, near=None)
        if own:
            ln = min(own, key=lambda ln: abs(ln["energy_kev"] - e_kev))
            pk.update(label=ln["label"], nuclide=ln["nuclide"], line_kev=ln["energy_kev"])
            continue
        far = [ln for ln in lines if ln["correction"] != correction
               and abs(ln["energy_kev"] - e_kev) <= 0.06 * ln["energy_kev"]]
        if far:
            pk["near"] = min(far, key=lambda ln: abs(ln["energy_kev"] - e_kev))["label"]
    return peaks


def _half_maximum_sigma(counts: np.ndarray, edges: np.ndarray, energy: float) -> float:
    """A first idea of a peak's width from the counts: the full width at half maximum of the highest point near
    ``energy``, over a baseline from the window's ends, as σ. Zero if there is no peak to see."""
    centres = (edges[:-1] + edges[1:]) / 2
    near = np.abs(centres - energy) <= 0.02 * energy
    if not near.any():
        return 0.0
    # Smoothed at several scales; the scale at which the peak stands out most against the counting noise wins,
    # so a wide, low peak is not mistaken for a spike of noise.
    best = None
    for k in (5, 15, 45, 135):
        if k >= len(counts) // 3:
            break
        smooth = np.convolve(counts, np.ones(k) / k, mode="same")
        edge_bins = max(len(smooth) // 10, 2)
        base = float(np.mean(np.r_[smooth[:edge_bins], smooth[-edge_bins:]]))
        top = int(np.flatnonzero(near)[np.argmax(smooth[near])])
        height = smooth[top] - base
        significance = height * math.sqrt(k) / math.sqrt(max(base, 1.0))
        if height > 0 and (best is None or significance > best[0]):
            best = (significance, smooth, base, top, height)
    if best is None:
        return 0.0
    _, smooth, base, top, height = best
    half = base + height / 2
    left, right = top, top
    while left > 0 and smooth[left] > half:
        left -= 1
    while right < len(smooth) - 1 and smooth[right] > half:
        right += 1
    return float((centres[right] - centres[left]) / FWHM_PER_SIGMA)


@dataclass
class Result:
    """What the analysis found."""

    settings: dict
    #: The steps in words and numbers, in order.
    steps: list
    fit: PeakFit
    #: The γ-ray yield per collision-equivalent and the quantities that went into it.
    area: float
    efficiency: float
    correlation: float
    branch: float
    normalisation: str
    norm_counts: float
    #: B(E2↑) in e²fm⁴ with its statistical and total systematic relative uncertainties, and in other units.
    b_e2fm4: float
    statistical: float
    systematic: float
    #: The relative uncertainty from the size of the Monte Carlo sample itself (not of the experiment): the
    #: scatter between seeds. Simulate more events to shrink it.
    monte_carlo: float
    b_e2b2: float
    b_wu: float
    budget: dict
    #: The value put into the setup and the pull (difference over the total uncertainty).
    truth_e2fm4: float
    pull: float
    shape: dict
    #: Counts in the peak per shift and the beam time for the wanted precision, hours.
    counts_per_shift: float
    hours_for_precision: float
    notes: list = field(default_factory=list)
    #: For the analysis of a run (backlog 63): where the statistical uncertainty comes from ("the run's real
    #: events"), and what it would be for the whole run, scaled beyond the real part.
    statistics_from: str = "the Monte Carlo sample"
    whole_run_statistical: Optional[float] = None
    #: The gated, Doppler-corrected spectrum the peak was fitted to, for drawing it: {"edges" (MeV), "counts",
    #: "randoms" (expected per bin), "correction" ("projectile" or "recoil"), "energy" (the transition, MeV)}.
    spectrum: Optional[dict] = None
    #: The peak search (step 3), by Doppler correction ("projectile", "recoil"): {"edges" (MeV), "counts" (random
    #: coincidences subtracted), "gate" (the particles in words), "nuclide" (the nucleus this correction is right
    #: for), "peaks": [{"energy_kev", "fwhm_kev", "area", "area_unc", "significance", "label", "nuclide",
    #: "line_kev", "near"}]}.
    search: Optional[dict] = None
    #: The peaks the analysis integrated, the ``[reaction]`` transition first and, when the events hold the other
    #: nucleus, its first transition: dicts with "nuclide", "label", "energy_kev", "correction", "gate", "spectrum"
    #: (as :attr:`spectrum`), "fit" (a :class:`PeakFit` as a dict), "area" and "area_unc" (random coincidences
    #: subtracted), "randoms", "efficiency", "correlation", "gate_efficiency", "yield" (γ rays emitted),
    #: "per_particle" (the yield per particle in the same rings), "expected_per_particle" (what the level
    #: scheme gives to all orders; the other nucleus only), "ratio" (measured over expected), "found" (the peak
    #: of the search it is, or None) and "blends" (other transitions inside the ±3σ window).
    lines: Optional[list] = None
    #: What first-order theory makes of the same yield, e²fm⁴ (equal to :attr:`b_e2fm4` when first order is all
    #: there is: one state, or first-order events).
    b_first_order: Optional[float] = None
    #: The particle gate the result used ("inelastic" or "all"), and for the inelastic gate the share of the
    #: excitations in the rings at whose lab direction the elastic and inelastic groups overlap through the depth
    #: of the target; above :data:`GATE_OVERLAP_WARN` :attr:`gate_warning` says to use "all" (the first note too).
    gate: Optional[str] = None
    gate_overlap: Optional[float] = None
    gate_warning: Optional[str] = None

    @property
    def total_unc(self) -> float:
        return math.hypot(self.statistical, self.systematic)

    def to_dict(self) -> dict:
        """The result as plain data, to keep with a run (:func:`save_result`)."""
        return json.loads(json.dumps(asdict(self), default=_plain))

    @classmethod
    def from_dict(cls, d: dict) -> Result:
        d = {k: v for k, v in d.items() if k in {f.name for f in fields(cls)}}
        d["fit"] = PeakFit.of(d.get("fit") or {})
        return cls(**d)


#: The file in a run's folder that keeps its last analysis.
RESULT_FILE = "analysis.json"


def _plain(x):
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (set, tuple)):
        return list(x)
    return str(x)


def save_result(folder, result: Result) -> None:
    """Keep ``result`` in a run's folder, so that the reports find it after the setup changes or the app closes."""
    from datetime import datetime, timezone
    from pathlib import Path

    d = result.to_dict()
    d["saved"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    (Path(folder) / RESULT_FILE).write_text(json.dumps(d, indent=1), encoding="utf-8")


def load_result(folder) -> Optional[Result]:
    """The analysis kept in a run's folder (:func:`save_result`), or None if there is none or it cannot be read."""
    from pathlib import Path

    p = Path(folder) / RESULT_FILE
    if not p.is_file():
        return None
    try:
        return Result.from_dict(json.loads(p.read_text(encoding="utf-8")))
    except (ValueError, TypeError, KeyError):
        return None


class Analysis:
    """The analysis of one simulated experiment (``gammas`` from
    :func:`physim.nuclear.gamma_events.simulate_gammas`). Every :meth:`run` is kept in :attr:`history` with
    its settings, so a change and its effect are on record."""

    def __init__(self, experiment, gammas, run=None):
        self.experiment = experiment
        self.gammas = gammas
        self.events = gammas.events
        #: The run the γ rays were taken in (:class:`physim.nuclear.runs.RunData`), or None for a weighted sample.
        self.taken = run
        self.history: list = []
        exc = experiment.excitation
        if exc is None:
            raise ValueError("the analysis needs a Coulomb-excitation setup")
        self.e0 = exc.energy_mev
        self.rates = Rates(experiment)
        self._correlation = {(r["particle_detector"], r["gamma_detector"]): r["detector_factor"]
                             for r in correlation_table(experiment)}
        self._layers = stack(experiment)
        self._e_mid = float(beam_energy_at(experiment, 0, [self._layers[0].thickness / 2], self._layers)[0])
        #: The γ rays of the [reaction] nucleus (all of them when the events hold one nucleus).
        self._own = (gammas["nucleus"] == 0) if "nucleus" in gammas.columns else np.ones(len(gammas), dtype=bool)
        self._emitters = emitters(experiment)
        self._lines = transition_lines(experiment)
        self._other: Optional[dict] = None
        #: The events have one particle per collision (elastic scattering is what the excitation leaves), so the
        #: particles in the rings are all of Rutherford scattering. A run taken before that counted the excited
        #: particles on top of it: there the normalisation is to the particles alone, as it was.
        self._whole = run is None or bool(run.summary.get("one_particle"))

    # -- the steps --------------------------------------------------------------------------------------------------

    def _gate(self, s: Settings) -> tuple:
        """Mask over the γ rays whose particle passes the gate, and the number of particles in the same
        rings (for the Rutherford normalisation)."""
        g, ev = self.gammas, self.events
        c = ev.columns
        rows = g["particle"]
        dets = list(ev.detectors) if s.detectors is None else list(s.detectors)
        keep = np.isin(c["detector"][rows], [ev.detectors.index(d) for d in dets]) & g["counted"]
        elastic = np.isin(c["detector"], [ev.detectors.index(d) for d in dets]) & c["counted"]
        if s.rings is not None:
            keep &= np.isin(c["segment_i"][rows], list(s.rings))
            elastic &= np.isin(c["segment_i"], list(s.rings))
        if s.particle_energy is not None:
            lo, hi = sorted(float(x) for x in s.particle_energy)
            keep &= (c["measured"][rows] >= lo) & (c["measured"][rows] <= hi)
        exc = self.experiment.excitation
        if not self._whole:
            excited = [i for i, lab in enumerate(ev.channels) if "excited" in lab]
            elastic &= ~np.isin(c["channel"], excited)
        if s.particle_gate == "inelastic":
            ion = beam_ion(self.experiment)
            target = data.nuclide(max(self._layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
            tb_el = TwoBody(ion, target, self._e_mid)
            tb_in = TwoBody(ion, target, self._e_mid, excitation_mev=self.e0,
                            excite="recoil" if exc.excite == "target" else "ejectile")
            th = c["theta"][rows]
            recoil = c["recoil"][rows]
            pred = {}
            for name, tb in (("elastic", tb_el), ("inelastic", tb_in)):
                e = np.full(len(rows), np.nan)
                for particle, m in (("ejectile", ~recoil), ("recoil", recoil)):
                    if m.any():
                        first = tb.at_lab(th[m], particle)[0]
                        e[m] = np.nan_to_num(np.asarray(first.energy, dtype=float))
                pred[name] = e
            measured = c["measured"][rows]
            # What each group would measure: the kinematic energy, less the losses on the way out of the target
            # and in the detector's dead layer, as the generator computes them.
            expect = {name: self._measured_energy(rows, pred[name]) for name in pred}
            fwhm = np.array([_q(d.resolution).to("MeV") if d.resolution is not None else 0.05
                             for d in self.experiment.detectors])[c["detector"][rows]]
            # The reaction's depth in the target spreads each group by the energy lost on the way out from the
            # back: the gate is that wide, plus the resolution.
            spread = np.abs(self._measured_energy(rows, pred["inelastic"], depth="back") - expect["inelastic"])
            near = np.abs(measured - expect["inelastic"]) <= s.gate_width * np.maximum(fwhm, 0.02) + spread / 2
            closer = np.abs(measured - expect["inelastic"]) < np.abs(measured - expect["elastic"])
            # Whether the two groups can be told apart at all: through a thick target they spread over the
            # energies of every depth, and may overlap (backlog: Allmond 2014, fix 7).
            base = keep & self._own
            overlap = np.zeros(len(rows), dtype=bool)
            if base.any():
                overlap[base] = self._groups_overlap(rows[base], fwhm[base])
            w_base = g["weight"][base].sum()
            self._gate_overlap = float(g["weight"][base & overlap].sum() / w_base) if w_base > 0 else 0.0
            keep &= near & closer
        # What the gate does, known here because the simulation knows the truth: the share of the excitations it
        # keeps, and the particles it lets through (they carry no γ ray, so they only add to the randoms).
        in_dets = self._in_rings(s)
        w_all = g["weight"][in_dets & self._own].sum()
        self._gate_efficiency = float(g["weight"][keep & self._own].sum() / w_all) if w_all > 0 else 1.0
        return keep, float(np.sum(c["weight"][elastic]) * ev.beam_time_s * g.live_fraction)

    def _groups_overlap(self, rows: np.ndarray, fwhm: np.ndarray) -> np.ndarray:
        """For each particle row: whether the elastic and the inelastic group overlap at its lab direction. Each
        group spans the energies it is measured with from every depth of the target (the kinematic energy at the
        beam energy there, less the losses on the way out from there), widened by half the detector's resolution
        on each side; the two overlap when the elastic span covers more than :data:`GATE_GROUP_COVER` of the
        inelastic one."""
        exp, exc = self.experiment, self.experiment.excitation
        lay = self._layers[0]
        ion = beam_ion(exp)
        target = data.nuclide(max(lay.nuclides.items(), key=lambda kv: kv[1])[0]).name
        m1, m2 = data.nuclide(ion).atomic_mass_u, data.nuclide(target).atomic_mass_u
        threshold = 1.05 * self.e0 * (m1 + m2) / m2
        depths = {"front": 0.0, "middle": 0.5 * lay.thickness, "back": 0.999 * lay.thickness}
        energies = {k: float(beam_energy_at(exp, 0, [z], self._layers)[0]) for k, z in depths.items()}
        th, recoil = self.events.columns["theta"][rows], self.events.columns["recoil"][rows]
        span = {}
        for name, kw in (("elastic", {}), ("inelastic", {"excitation_mev": self.e0, "excite": "recoil"
                                                          if exc.excite == "target" else "ejectile"})):
            values = []
            for where, e_beam in energies.items():
                if name == "inelastic":
                    e_beam = max(e_beam, threshold)  # deeper, nothing is excited
                if not e_beam > 0:
                    continue
                tb = TwoBody(ion, target, e_beam, **kw)
                e = np.full(len(rows), np.nan)
                for particle, m in (("ejectile", ~recoil), ("recoil", recoil)):
                    if m.any():
                        e[m] = np.asarray(tb.at_lab(th[m], particle)[0].energy, dtype=float)
                values.append(self._measured_energy(rows, np.nan_to_num(e), depth=where))
            if not values:
                return np.zeros(len(rows), dtype=bool)
            values = np.array(values)
            span[name] = (values.min(axis=0) - fwhm / 2, values.max(axis=0) + fwhm / 2)
        (el_lo, el_hi), (in_lo, in_hi) = span["elastic"], span["inelastic"]
        common = np.clip(np.minimum(el_hi, in_hi) - np.maximum(el_lo, in_lo), 0.0, None)
        cover = common / np.maximum(in_hi - in_lo, 1e-9)
        return cover > GATE_GROUP_COVER

    def _in_rings(self, s: Settings) -> np.ndarray:
        """Mask over the γ rays whose particle was counted in the detectors and rings of the settings, whatever
        its energy."""
        g, ev = self.gammas, self.events
        c = ev.columns
        rows = g["particle"]
        dets = list(ev.detectors) if s.detectors is None else list(s.detectors)
        m = np.isin(c["detector"][rows], [ev.detectors.index(d) for d in dets]) & g["counted"]
        if s.rings is not None:
            m &= np.isin(c["segment_i"][rows], list(s.rings))
        return m

    def _measured_energy(self, rows, energy: np.ndarray, depth: str = "middle") -> np.ndarray:
        """The energy a particle leaving the reaction with ``energy`` along each row's track would measure: after
        the way out of the target (from its front, middle or back) and the detector's dead layer, as the
        generator does."""
        from .detectors import Array

        c = self.events.columns
        th, ph = c["theta"][rows], c["phi"][rows]
        dirs = np.c_[np.sin(np.radians(th)) * np.cos(np.radians(ph)), np.sin(np.radians(th)) * np.sin(np.radians(ph)),
                     np.cos(np.radians(th))]
        out = np.zeros(len(rows))
        ion = beam_ion(self.experiment)
        target = data.nuclide(max(self._layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        array = Array.from_experiment(self.experiment)
        for particle, m in (("ejectile", ~c["recoil"][rows]), ("recoil", c["recoil"][rows])):
            if not m.any():
                continue
            species = ion if particle == "ejectile" else target
            z = self._layers[0].thickness * {"front": 0.001, "middle": 0.5, "back": 0.999}[depth]
            e_out = exit_energy(self.experiment, self._layers, 0, np.full(m.sum(), z), species,
                                np.nan_to_num(energy[m]), dirs[m])
            for i, g in enumerate(array.geometries):
                mm = m & (c["detector"][rows] == i)
                if not mm.any():
                    continue
                resp = array.response(g.name, species)
                e_face = e_out[mm[m]]
                # Incidence on the face, in steps of 2°, to call the stopping tables a few times.
                cos_i = np.clip(-(dirs[mm] @ g.n), 1e-6, 1.0)
                inc = np.round(np.degrees(np.arccos(cos_i)) / 2) * 2
                dep = np.empty(len(e_face))
                for angle in np.unique(inc):
                    sel = inc == angle
                    dep[sel] = np.asarray(resp.deposited(e_face[sel], float(angle)), dtype=float)
                out[mm] = dep
        return out

    def _spectrum(self, s: Settings, keep: np.ndarray) -> tuple:
        """The gated, corrected spectrum summed over the γ-ray detectors, and the random coincidences expected
        per bin in the run."""
        exc = self.experiment.excitation
        which = s.correction if s.correction != "emitter" else ("recoil" if exc.excite == "target" else "projectile")
        key = {"projectile": "corrected_projectile", "recoil": "corrected_recoil"}[which]
        g = self.gammas
        energies = self._energies(s, key)
        sigma = max(float(np.mean([c.response.fwhm(self.e0) for c in g.crystals])) / FWHM_PER_SIGMA, 1e-4)
        width = sigma / 4
        # The window holds the peak however much Doppler width the segments leave: up to ±5% of the energy.
        half = max(12 * sigma, 0.05 * self.e0)
        edges = np.arange(self.e0 - half, self.e0 + half + width, width)
        t = g.events.beam_time_s * g.live_fraction
        counts, _ = np.histogram(energies[keep], bins=edges, weights=g["weight"][keep] * t)
        sigma = min(max(_half_maximum_sigma(counts, edges, self.e0), sigma), 0.03 * self.e0)
        rnd = np.zeros(len(edges) - 1)
        if s.subtract_randoms:
            dets = list(self.events.detectors) if s.detectors is None else list(s.detectors)
            for d in dets:
                for name in g.detector_names():
                    rnd += g.random_spectrum(name, d, edges) * t
            rng = np.random.default_rng(g.seed + 104729)
            counts = counts + rng.poisson(rnd)
        self._variance = self._bin_variance(energies, keep, edges, rnd)
        return counts, edges, rnd, sigma, which

    def _bin_variance(self, energies: np.ndarray, keep: np.ndarray, edges: np.ndarray, rnd: np.ndarray):
        """The variance of each bin of a spectrum of the weighted sample scaled to the beam time: the sum of the
        squared weights of its entries, and the random coincidences drawn into it. None for a run, where each
        event counts once and the variance of a bin is its counts."""
        if self.taken is not None:
            return None
        g = self.gammas
        t = g.events.beam_time_s * g.live_fraction
        sumw2, _ = np.histogram(energies[keep], bins=edges, weights=(g["weight"][keep] * t) ** 2)
        return sumw2 + rnd

    def _energies(self, s: Settings, key: str, offset_mm: Optional[float] = None) -> np.ndarray:
        """The corrected energies with the geometry the analysis assumes: the true one, or the target moved by
        ``Settings.target_offset_mm`` (or ``offset_mm``) along the beam."""
        offset = s.target_offset_mm if offset_mm is None else offset_mm
        if not offset:
            return self.gammas[key]
        from .alignment import with_offset
        from .gamma_events import recorrect

        return recorrect(self.gammas, with_offset(self.experiment, offset))[key]

    def _mean_probability(self, s: Settings, keep: Optional[np.ndarray] = None, energy_scale: float = 1.0,
                          mode: str = "first", exact: bool = False) -> float:
        """⟨P⟩/B: the excitations of the ``[reaction]`` state expected per particle counted in the gate's rings, per
        unit B(E2↑), in first order. The Rutherford normalisation divides the excitations by those particles, so
        this is what B(E2) is read against.

        It comes from the forward model itself, the cross sections the rates and the event generator use
        (:func:`physim.nuclear.rates.reaction_at`), not from a formula of its own: at each particle's lab angle,
        the Coulomb-excitation cross section over that of every particle that reaches the angle (elastic and
        excited together), each integrated over the depth of its layer with the beam energy there
        (:meth:`_forward_ratio`), then averaged over the particles in the rings, weighted as the events are. So
        the excited particle's own CM angle and lab-to-CM Jacobian, the excitation probability through the depth
        of a thick target, and whatever Rutherford factor the forward cross section carries all enter as they do
        in the simulation.

        ``energy_scale`` multiplies the energy the beam enters the target with (the beam-energy systematic), which
        then slows down from there: deep in a thick target that is a larger part of what is left. ``exact`` takes P
        from the first-order Coulomb excitation at each energy rather than from the forward model's interpolation
        through the depth (the two agree to a few parts in 10⁴; the budget compares like with like). ``mode`` =
        "forward" gives instead the forward model's own excitations per particle (every level a level scheme
        populates, to all orders if the events are), not per unit B(E2)."""
        rows = self._elastic_rows(s)
        if not rows.any():
            return 0.0
        c = self.events.columns
        th, rec, w = c["theta"][rows], c["recoil"][rows], c["weight"][rows]
        r = np.zeros(len(th))
        for particle, m in (("ejectile", ~rec), ("recoil", rec)):
            if m.any():
                grid, ratio = self._forward_ratio(particle, energy_scale, mode, exact)
                r[m] = np.interp(th[m], grid, ratio)
        out = float(np.average(r, weights=w)) if w.sum() > 0 else 0.0
        return out / self.experiment.excitation.b_up_e2fm if mode == "first" else out

    def _first_order(self, channel, energy: float, exact: bool):
        """P(θ) of the ``[reaction]`` state alone, in first order, at beam energy ``energy`` for a
        Coulomb-excitation channel: the forward model's own P through the depth when that is first order of the
        one state (``exact`` False, and no level scheme), else the first-order Coulomb excitation at that energy."""
        try:
            dp = depth_probability(self.experiment, channel, self._layers)
        except ValueError:
            return _Zero()
        if not exact and not isinstance(dp.nodes[0], SchemeProbability):
            return dp.at(energy)
        try:
            return coulex_for(self.experiment, channel, self._layers, energy=energy)
        except ValueError:  # below the excitation threshold
            return _Zero()

    def _forward_ratio(self, particle: str, energy_scale: float = 1.0, mode: str = "first",
                       exact: bool = False) -> tuple:
        """(lab angles, excited over all) for ``particle`` ("ejectile" or "recoil"): the forward model's
        Coulomb-excitation cross section in the lab, summed over the channels that excite the ``[reaction]``
        nucleus, over the lab cross section of every channel that sends that particle there (elastic scattering,
        which the excitation reduces by 1 − P, and the excitation itself; for a run taken before the events had
        one particle per collision, elastic scattering alone, as the normalisation counts). Each channel is
        integrated over the depth of its layer the beam crosses (:data:`FORWARD_DEPTH_POINTS`), at the beam
        energy at each depth, through the kinematics of its own reaction (the excited particle's CM angle and
        Jacobian, both branches) as :class:`~physim.nuclear.rates.Rates` integrates it. With ``mode`` "first"
        the excitation is the first-order one of the ``[reaction]`` state (:meth:`_first_order`) inside the same
        cross section. Tabulated every :data:`FORWARD_STEP_DEG` over the lab angles of the counted particles,
        and kept."""
        key = (particle, round(float(energy_scale), 12), mode, bool(exact))
        kept = self.__dict__.setdefault("_forward", {})
        if key in kept:
            return kept[key]
        exp, layers = self.experiment, self._layers
        c = self.events.columns
        th = c["theta"][c["counted"] & (c["recoil"] == (particle == "recoil"))]
        step = FORWARD_STEP_DEG
        lo = max(math.floor(float(th.min()) / step) - 1, 1) * step if len(th) else 0.0
        hi = min(math.ceil(float(th.max()) / step) + 1, round(180.0 / step) - 1) * step if len(th) else 180.0
        grid = np.arange(lo, hi + step / 2, step) if hi > lo else np.array([lo, lo + step])
        num, den = np.zeros(len(grid)), np.zeros(len(grid))
        beam = beam_ion(exp)
        entering = exp if energy_scale == 1.0 else _with_beam_energy(exp, exp.beam.energy_mev * energy_scale)
        x, wq = np.polynomial.legendre.leggauss(FORWARD_DEPTH_POINTS)
        cos_t = math.cos(math.radians(tilt_deg(exp)))
        for ch in channels(layers, exp):
            lay = layers[ch.layer]
            entry = float(beam_energy_at(entering, ch.layer, [0.0], layers)[0])
            if not entry > 0:
                continue
            # The part of the layer the beam crosses, and the share of the layer's nuclei at each depth point.
            reach = min(float(stopping(beam, lay.material).range(entry)) * cos_t, lay.thickness)
            zs = (x + 1) / 2 * reach
            share = wq / 2 * reach / lay.thickness
            es = np.asarray(beam_energy_at(entering, ch.layer, zs, layers), dtype=float)
            excited = ch.excitation is not None
            reacs, ws, firsts = [], [], []
            for e, wz in zip(es, share):
                if not e > 0:
                    continue
                try:
                    reac = reaction_at(exp, ch, float(e), layers)
                except ValueError:  # below the excitation threshold, or no energy left
                    continue
                reacs.append(reac)
                ws.append(wz)
                if excited and mode == "first":
                    # The same forward cross section, with P that of the [reaction] state alone in first order.
                    firsts.append(type(reac)(beam, ch, float(e), self._first_order(ch, float(e), exact)))
            if not reacs:
                continue
            weights = np.asarray(ws)[:, None] * ch.atoms_per_cm2
            total = (weights * sum(sig for sig, _ in self.rates._branches_many(reacs, grid, particle))).sum(axis=0)
            if not excited or self._whole:
                den += total
            if excited:
                if mode == "first":
                    num += (weights * sum(sig for sig, _ in self.rates._branches_many(firsts, grid, particle))
                            ).sum(axis=0)
                else:
                    num += total
        ratio = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
        kept[key] = (grid, ratio)
        return kept[key]

    def _elastic_rows(self, s: Settings) -> np.ndarray:
        c = self.events.columns
        ev = self.events
        dets = list(ev.detectors) if s.detectors is None else list(s.detectors)
        m = np.isin(c["detector"], [ev.detectors.index(d) for d in dets]) & c["counted"]
        if s.rings is not None:
            m &= np.isin(c["segment_i"], list(s.rings))
        if self._whole:  # one particle per collision: every particle in the rings, as a detector counts them
            return m
        excited = [i for i, lab in enumerate(ev.channels) if "excited" in lab]
        return m & ~np.isin(c["channel"], excited)

    # -- the run ----------------------------------------------------------------------------------------------------

    def run(self, settings: Optional[Settings] = None) -> Result:
        s = settings or Settings()
        exc = self.experiment.excitation
        g = self.gammas
        steps = []
        keep, elastic_counts = self._gate(s)
        steps.append({"step": "particle gate", "text": f"{s.particle_gate} group, {s.gate_width:g} FWHM wide, "
                      + ("all rings" if s.rings is None else f"rings {list(s.rings)}") + f": {int(keep.sum())} of "
                      f"{len(keep)} γ rays kept; {elastic_counts:.4g} particles in the same rings"})
        counts, edges, rnd, sigma, which = self._spectrum(s, keep)
        steps.append({"step": "Doppler correction", "text": f"for the {which}"})
        search = self._search(s, keep)
        steps.append({"step": "peak search", "text": self._search_text(search)})
        fit = self._peak(s, 0, self.e0, counts, edges, sigma, self._variance)
        area, area_unc = self._area(s, fit, edges, rnd)
        steps.append({"step": "peak fit", "text": f"failed: {fit.failed}" if fit.failed else
                      f"{fit.shape} within ±{s.fit_half_width:g} FWHM: area "
                      f"{fit.area:.1f} ± {fit.area_unc:.1f} at {1e3 * fit.centroid:.1f} keV, FWHM {fit.fwhm_kev:.1f} "
                      f"keV, χ² {fit.chi2:.0f} for {fit.bins} bins"
                      + ("" if fit.tails is None else f"; {100 * (1 - fit.tails['inside']):.1f} % of the peak lies "
                         "beyond ±3σ, in its tails" + self._simulated_text(fit))})
        steps.append({"step": "random coincidences", "text": f"{fit.randoms:.1f} expected in the peak, subtracted"
                      if s.subtract_randoms else "not subtracted"})
        # The efficiency, branch and correlation of the gated γ rays: weighted over the pairs that contributed.
        rows = g["particle"][keep]
        w = g["weight"][keep]
        eff, corr = self._efficiency_and_correlation(keep)
        branch = 1.0  # one state decaying to the ground state; conversion is below 1e-3 above 1 MeV
        gate_eff = getattr(self, "_gate_efficiency", 1.0) if s.particle_gate == "inelastic" else 1.0
        yield_ = area / (eff * corr * branch * gate_eff)
        steps.append({"step": "yield", "text": f"area / (efficiency {eff:.4%} × correlation {corr:.3f} × branch "
                      f"{branch:g} × gate efficiency {gate_eff:.3f}) = {yield_:.4g} excitations in the gated "
                      "rings"})
        if s.particle_gate == "inelastic":
            notes_gate = (f"The gate keeps {gate_eff:.1%} of the excitations (known from the simulation, as an "
                          "experimentalist would take it from one); the yield is divided by that.")
        p_per_b = self._mean_probability(s, keep)
        notes = []
        norm_var = None
        # The other nucleus's peak, when the events hold it: found, fitted and integrated like the first.
        other = self._other_line(s, elastic_counts, search)
        if s.normalisation == "target" and other is not None and not s.reference:
            norm_counts = other["yield"]
            ok = norm_counts > 0 and p_per_b > 0 and other["expected_per_particle"] > 0
            # The particles the peak stands for: its yield over what the level scheme gives per particle.
            b = yield_ * other["expected_per_particle"] / norm_counts / p_per_b if ok else float("nan")
            ref_unc = 0.03
            steps.append({"step": "normalisation", "text": f"to the measured {other['label']} peak at "
                          f"{other['energy_kev']:.1f} keV (corrected for the {other['correction']}): area "
                          f"{other['area']:.1f} ± {other['area_unc']:.1f}, over its efficiency "
                          f"{other['efficiency']:.4%} and correlation {other['correlation']:.3f} that is "
                          f"{norm_counts:.4g} γ rays; its level scheme gives {other['expected_per_particle']:.4g} "
                          f"per particle, to all orders, so the peak stands for "
                          f"{norm_counts / other['expected_per_particle']:.4g} particles "
                          f"({elastic_counts:.4g} were counted); ⟨P⟩/B(E2) = {p_per_b:.3e} per e²fm⁴"})
            notes.append(f"Normalised to the {other['label']} peak of the same run: the beam current, the target "
                         "thickness, the dead time and the number of particles counted cancel. Its matrix "
                         f"elements are taken as known to {100 * ref_unc:.0f}%.")
            # Both peaks count: the reference peak's own statistics, not those of a number of particles.
            norm_var = (other["area_unc"] / other["area"]) ** 2 if other["area"] > 0 else float("inf")
        elif s.normalisation == "rutherford":
            norm_counts = elastic_counts
            # Excitations / particles in the same rings = ⟨P⟩ in first order.
            b = yield_ / norm_counts / p_per_b if norm_counts > 0 and p_per_b > 0 else float("nan")
            steps.append({"step": "normalisation", "text": f"to the {norm_counts:.4g} particles in the same "
                          f"rings: ⟨P⟩ = {yield_ / norm_counts:.3e}, and ⟨P⟩/B(E2) = {p_per_b:.3e} per "
                          "e²fm⁴ (first order, from the forward cross sections at each particle's lab angle, through "
                          "the depth of the target)"})
            ref_unc = 0.0
        else:
            if not s.reference:
                raise ValueError("normalisation to the other nucleus needs its level scheme in the setup (its "
                                 "peak in the run is then the reference), or Settings.reference: a known "
                                 "transition")
            ref = self._reference(s, keep, w)
            norm_counts = ref["counts"]
            b = yield_ / norm_counts * ref["b_up"] * ref["p_per_b"] / p_per_b if norm_counts > 0 else float("nan")
            ref_unc = float(s.reference.get("unc", 0.0))
            steps.append({"step": "normalisation", "text": f"to the reference transition at "
                          f"{1e3 * ref['energy']:.1f} keV, B(E2↑) = {ref['b_up']:.4g} e²fm⁴ ± {100 * ref_unc:.1f}%: "
                          f"{ref['raw_counts']:.4g} counts in its peak, {norm_counts:.4g} over its efficiency"})
            notes.append("The reference peak is the analytic expectation with Poisson noise: the simulation excites "
                         "one state.")
        # To all orders, when the events are: the coupled equations, with the matrix element as the unknown.
        b_first = b
        full = self._all_orders(s, b, p_per_b, fit) if area > 0 else None
        if full is not None:
            b = full["b"]
            steps.append({"step": "all orders", "text": full["text"]})
        statistical = math.sqrt((area_unc / area) ** 2 + (1 / max(norm_counts, 1.0) if norm_var is None
                                                          else norm_var)) if area > 0 else float("inf")
        mc = self._monte_carlo_uncertainty(s, keep, fit) if area > 0 else float("inf")
        if not area > 0:
            b = float("nan")
            notes.append(f"The peak fit failed: {fit.failed}. There is no area, and no B(E2)." if fit.failed else
                         "No peak: the gate or the rings leave too few γ rays to fit.")
        if fit.tails is not None and fit.tails.get("simulated"):
            notes.append("The counts are too few to fix the tails of the peak: their share and shape come from a "
                         "simulation of this setup, and the line \"peak tails (simulated)\" of the budget is what "
                         "that may be wrong by.")
        fits = [fit] + ([PeakFit.of(other["fit"])] if s.normalisation == "target" and other is not None
                        and not s.reference else [])
        budget = self._budget(s, keep, b, p_per_b, eff, corr, ref_unc, fits)
        if full is not None and full["gain"] > 1.0:
            # To all orders B(E2) moves more than the yield does (see _all_orders): so do its uncertainties.
            statistical *= full["gain"]
            mc *= full["gain"]
            budget = {k: v * full["gain"] for k, v in budget.items()}
        systematic = math.sqrt(sum(v**2 for v in budget.values()))
        truth = exc.b_up_e2fm
        total = b * math.hypot(statistical, systematic) if b == b else float("nan")
        pull = (b - truth) / total if total and total > 0 else float("nan")
        a_mass = data.nuclide(ex_nuclide(self.experiment)).A
        wu = weisskopf_unit(2, a_mass)
        shape = self._shape(b, a_mass)
        t_s = g.events.beam_time_s * g.live_fraction
        per_shift = area / t_s * SHIFT_H * 3600 if t_s > 0 else 0.0
        # Precision 1/√N for the peak alone: N = 1/precision² counts.
        hours = (1 / s.wanted_precision**2) / (area / t_s) / 3600 if area > 0 else float("inf")
        gate_overlap, gate_warning = None, None
        if s.particle_gate == "inelastic":
            notes.append(notes_gate)
            gate_overlap = getattr(self, "_gate_overlap", 0.0)
            if gate_overlap > GATE_OVERLAP_WARN:
                gate_warning = (
                    f"This result used the inelastic particle gate, which cannot separate the groups here: for "
                    f"{gate_overlap:.0%} of the excitations in these rings the energies of the elastic group, "
                    "from every depth of the target, reach into those of the inelastic group. The gate then keeps "
                    f"only part of the excitations ({gate_eff:.1%}), from the depths it happens to pick out, and "
                    "the result leans on the simulation's word for that share. Use the gate \"all\" "
                    "(particle_gate=\"all\"): with the Rutherford normalisation every particle in the rings counts.")
                notes.insert(0, gate_warning)
                steps[0]["text"] += f". Warning: the elastic and inelastic groups overlap ({gate_overlap:.0%})"
        whole = None
        source = "the Monte Carlo sample"
        if self.taken is not None:
            # Every event of a run counts once: there are no Monte Carlo weights, and the statistical uncertainty is
            # the run's own, from its real part. The whole run, scaled beyond it, would have scale× the counts.
            mc = 0.0
            scale = self.taken.scale
            whole = statistical / math.sqrt(scale) if scale > 1 else statistical
            source = f"the run's real events ({int(round(area))} in the peak, {self.taken.real_s / 60:.3g} min of beam)"
            notes.append(f"The statistical uncertainty comes from {source}; each event counts once, so there is no "
                         "separate Monte Carlo uncertainty."
                         + (f" The whole {self.taken.duration_s / 3600:.3g} h run, scaled from it with its own Poisson "
                            f"statistics, would give {100 * whole:.2g}% instead of {100 * statistical:.2g}%."
                            if scale > 1 else ""))
        spectrum = {"edges": edges.tolist(), "counts": np.asarray(counts, dtype=float).tolist(),
                    "randoms": rnd.tolist(), "correction": which, "energy": float(self.e0)}
        lines = [self._own_line(s, fit, area, area_unc, eff, corr, gate_eff, yield_, elastic_counts, spectrum,
                                search)]
        if lines[0]["blends"]:
            notes.append("The ±3σ window of the peak also holds " + ", ".join(lines[0]["blends"]) + ": "
                         "those counts are in the area.")
        if other is not None:
            lines.append(other)
            failed = other["fit"].get("failed")
            notes.append(f"The fit of the {other['label']} peak ({other['energy_kev']:.1f} keV, corrected for the "
                         f"{other['correction']}) failed: {failed}." if failed else
                         f"The {other['label']} peak ({other['energy_kev']:.1f} keV, corrected for the "
                         f"{other['correction']}): {other['per_particle']:.4g} γ rays per particle "
                         f"measured, {other['expected_per_particle']:.4g} expected from its level scheme to all "
                         f"orders; measured over expected {other['ratio']:.3f}.")
        result = Result(s.as_dict(), steps, fit, area, eff, corr, branch, s.normalisation, norm_counts, b,
                        statistical, systematic, mc, b * 1e-4, b / wu, budget, truth, pull, shape, per_shift, hours,
                        notes, source, whole, spectrum, search, lines, b_first, s.particle_gate, gate_overlap,
                        gate_warning)
        if self.history:
            changed = {k: v for k, v in s.as_dict().items() if v != self.history[-1].settings.get(k)}
            if changed:
                result.notes.append(f"Changed from the previous run: {changed}; B(E2) was "
                                    f"{self.history[-1].b_e2fm4:.4g} e²fm⁴, now {b:.4g}.")
        self.history.append(result)
        return result

    # -- the peak search and the peaks of both nuclei ----------------------------------------------------------------

    def _detector_fwhm(self, energy) -> np.ndarray:
        """The γ-ray detectors' resolution (FWHM, MeV) at ``energy`` (MeV), averaged over the crystals."""
        e = np.maximum(np.atleast_1d(np.asarray(energy, dtype=float)), 1e-3)
        return np.mean([np.asarray(c.response.fwhm(e), dtype=float) for c in self.gammas.crystals], axis=0)

    def _random_rate(self, s: Settings, edges: np.ndarray) -> np.ndarray:
        """Random coincidences expected per bin and per second with the gate's particle detectors."""
        g = self.gammas
        rnd = np.zeros(len(edges) - 1)
        dets = list(self.events.detectors) if s.detectors is None else list(s.detectors)
        for d in dets:
            for name in g.detector_names():
                rnd += g.random_spectrum(name, d, edges)
        return rnd

    def _search(self, s: Settings, keep: np.ndarray) -> dict:
        """The peak search in both corrected spectra over their whole range. The spectrum corrected for the
        ``[reaction]`` nucleus is that of the particle gate; the other one is that of every particle counted in
        the same rings (the gate is made for the first nucleus). Random coincidences are subtracted: their lines
        are not Doppler shifted and would stand sharp in every spectrum."""
        from .gamma_events import spectrum_top

        g = self.gammas
        exc = self.experiment.excitation
        mine = "recoil" if exc.excite == "target" else "projectile"
        t = g.events.beam_time_s * g.live_fraction
        top = max(spectrum_top(self.experiment), 1.3 * self.e0)
        step = float(min(max(self._detector_fwhm(0.5)[0] / 4, 0.25e-3), 4e-3))
        edges = np.arange(0.0, top + step, step)
        rnd = self._random_rate(s, edges) * t if s.subtract_randoms else np.zeros(len(edges) - 1)
        drawn = np.random.default_rng(g.seed + 15_485_863).poisson(rnd) if s.subtract_randoms else rnd
        wide = self._in_rings(s)
        by_role = {em["correction"]: em["nuclide"] for em in self._emitters}
        out = {}
        for which in ("projectile", "recoil"):
            key = {"projectile": "corrected_projectile", "recoil": "corrected_recoil"}[which]
            m = keep if which == mine else wide
            energies = self._energies(s, key)
            true, _ = np.histogram(energies[m], bins=edges, weights=g["weight"][m] * t)
            counts = true + drawn - rnd
            var = self._bin_variance(energies, m, edges, rnd)
            var = true + rnd if var is None else var
            peaks = []
            for pk in find_peaks(counts, edges, self._detector_fwhm, variance=var):
                try:
                    fit = fit_peak(counts, edges, pk["energy"], pk["sigma"], s.fit_half_width,
                                   max_sigma=0.03 * pk["energy"], variance=var)
                except ValueError:
                    continue
                if fit.failed or not fit.area > 0:
                    continue
                # Two candidates the fit brings to one peak are that peak once (the more significant came first).
                if any(abs(fit.centroid - q["energy"]) < max(fit.sigma, q["sigma"]) for q in peaks):
                    continue
                peaks.append({"energy": fit.centroid, "sigma": fit.sigma, "significance": pk["significance"],
                              "area": fit.area, "area_unc": fit.area_unc})
            identify(peaks, self._lines, which)
            # A bump the fit finds wider than 3% of its energy is no peak: it is what the wrong correction, or a
            # Compton edge, leaves. One on a transition stays.
            peaks = [pk for pk in peaks if pk["label"] or FWHM_PER_SIGMA * pk["sigma"] <= 0.03 * pk["energy"]]
            peaks.sort(key=lambda pk: pk["energy"])
            out[which] = {
                "edges": edges.tolist(), "counts": counts.tolist(),
                "gate": ("the particle gate" if which == mine else "every particle counted in the same rings"),
                "nuclide": by_role.get(which),
                "peaks": [{"energy_kev": 1e3 * pk["energy"], "fwhm_kev": 1e3 * FWHM_PER_SIGMA * pk["sigma"],
                           "area": pk["area"], "area_unc": pk["area_unc"], "significance": pk["significance"],
                           "label": pk["label"], "nuclide": pk["nuclide"], "line_kev": pk["line_kev"],
                           "near": pk["near"]} for pk in peaks]}
        return out

    def _search_text(self, search: dict) -> str:
        parts = []
        for which in ("projectile", "recoil"):
            sp = search[which]
            named = [pk for pk in sp["peaks"] if pk["label"]]
            text = f"corrected for the {which}: {len(sp['peaks'])} peak{'s' if len(sp['peaks']) != 1 else ''} found"
            if named:
                text += ", " + "; ".join(f"{pk['label']} at {pk['energy_kev']:.1f} keV" for pk in named)
            elif sp["nuclide"]:
                text += f", none on a transition of {sp['nuclide']}"
            parts.append(text)
        return ". ".join(parts)

    def _found(self, search: dict, correction: str, energy_kev: float, half_kev: float):
        """The peak of the search nearest ``energy_kev`` in the spectrum of ``correction``, within ``half_kev``."""
        near = [pk for pk in search[correction]["peaks"] if abs(pk["energy_kev"] - energy_kev) <= half_kev]
        return min(near, key=lambda pk: abs(pk["energy_kev"] - energy_kev)) if near else None

    def _blended(self, s: Settings, nucleus: int, energy_kev: float, fit: PeakFit) -> bool:
        """Whether the nucleus has another line in the fit window, which tails would take for their own."""
        half = 1e3 * s.fit_half_width * FWHM_PER_SIGMA * fit.sigma
        return any(ln["nucleus"] == nucleus and abs(ln["energy_kev"] - energy_kev) > 0.05
                   and abs(ln["energy_kev"] - 1e3 * fit.centroid) <= half for ln in self._lines)

    def _peak(self, s: Settings, nucleus: int, energy: float, counts, edges, sigma: float, variance,
              like: Optional[PeakFit] = None) -> PeakFit:
        """The fit of the peak at ``energy`` (MeV) of one nucleus's corrected spectrum, ``sigma`` a first idea of its
        width. A Gaussian first; with its tails where the counts show them (:func:`fit_tails`); in a run whose counts
        cannot, with those of the simulated peak of the same setup (:meth:`simulated_peak`), and as the Gaussian where
        that one shows none either. As the Gaussian, too, when the nucleus has another line where the tails would be.
        A fit that ends at a bound is :attr:`PeakFit.failed`, with no area. ``like``: a fit of the same peak whose
        kind of tails (its own, the simulated ones, none) this one is to have too."""
        top = 0.03 * energy  # no Doppler residue is wider than a few percent
        # Two passes: the peak after the Doppler correction is wider than the resolution by what the segments
        # leave, so the window is set from a first fit's width.
        first = fit_peak(counts, edges, energy, sigma, 3 * s.fit_half_width, max_sigma=top, variance=variance)
        fit = fit_peak(counts, edges, energy, sigma if first.failed else first.sigma, s.fit_half_width, max_sigma=top,
                       variance=variance)
        if like is not None and like.tails is None:
            return fit
        if self._blended(s, nucleus, 1e3 * energy, fit):
            return self._checked(fit)
        simulated = like is not None and like.tails["simulated"]
        if not fit.failed and not simulated:
            own = fit_tails(counts, edges, fit, s.fit_half_width, variance=variance,
                            worth=TAILS_WORTH if like is None else -math.inf)
            if own.tails is not None:
                return own
        if self.taken is None or not float(np.sum(counts)) > 0:
            return self._checked(fit)
        sim = self.simulated_peak(s, nucleus)
        if sim is None:
            return self._checked(fit)
        shape = tail_shape(sim)
        if shape is not None:
            got = fit_tails(counts, edges, fit, s.fit_half_width, shape=shape)
            if got.tails is not None:
                gauss = sim.tails["gaussian_area"]
                got.tails.update(simulated_loss=max(1 - gauss / sim.area, 0.0),
                                 gaussian_area=None if fit.failed else fit.area)
            return self._checked(got)
        # A Gaussian peak: the free fit, unless it is not of a width, or not at a place, that the setup can give it.
        if (fit.failed or not 1 / SHAPE_REACH < fit.sigma / sim.sigma < SHAPE_REACH
                or abs(fit.centroid - sim.centroid) > SHAPE_REACH * sim.sigma):
            fit = fit_peak(counts, edges, sim.centroid, sim.sigma, s.fit_half_width, max_sigma=SHAPE_REACH * sim.sigma)
        return self._checked(fit)

    @staticmethod
    def _checked(fit: PeakFit) -> PeakFit:
        """``fit``, without its area if it failed."""
        return replace(fit, area=0.0, area_unc=0.0, gauss_area=0.0) if fit.failed else fit

    def _area(self, s: Settings, fit: PeakFit, edges, rnd) -> tuple:
        """The fit's area less the random coincidences expected in its ±3σ window, and its uncertainty."""
        window = np.abs((edges[:-1] + edges[1:]) / 2 - fit.centroid) < 3 * fit.sigma
        fit.randoms = float(rnd[window].sum()) if s.subtract_randoms and not fit.failed else 0.0
        return fit.area - fit.randoms, math.hypot(fit.area_unc, math.sqrt(max(fit.randoms, 0.0)))

    @staticmethod
    def _simulated_text(fit: PeakFit) -> str:
        t = fit.tails
        if not t.get("simulated"):
            return ""
        kev = 1e3 * FWHM_PER_SIGMA * t["simulated_sigma"]
        return (", whose share and shape are those of a simulation of this setup (the counts are too few to fix them; "
                f"the simulated peak has a FWHM of {kev:.1f} keV)")

    def simulated_peak(self, settings: Optional[Settings] = None, nucleus: int = 0) -> Optional[PeakFit]:
        """The peak that a simulation of the run's own setup gives under the same gate, rings, crystals, Doppler
        correction and fit window, with far more counts than a run and its own seed, fitted as a weighted sample's
        is: its tails (:func:`tail_shape`) are what a run of few counts takes for its own. Kept, so that a second fit
        does not simulate again. None when nothing of it could be fitted."""
        s = settings or Settings()
        which = (s.detectors, s.rings, nucleus)
        key = json.dumps([which, s.particle_gate, s.particle_energy, s.gate_width, s.correction, s.fit_half_width,
                          s.target_offset_mm], default=_plain)
        kept = self.__dict__.setdefault("_simulated", {})
        if key not in kept:
            kept[key] = None
            g = self._simulated_gammas(s, nucleus)
            if g is not None:
                a = self._of(g)
                quiet = replace(s, subtract_randoms=False)
                try:
                    if nucleus == 0:
                        counts, edges, _, sigma, _ = a._spectrum(quiet, a._gate(quiet)[0])
                        fit = a._peak(quiet, 0, self.e0, counts, edges, sigma, a._variance)
                    else:
                        o = self._other_setup()
                        counts, edges, _, sigma, var = a._other_spectrum(quiet, o)
                        fit = a._peak(quiet, 1, o["energy"], counts, edges, sigma, var)
                    if not fit.failed and fit.area > 0:
                        if fit.tails is not None:
                            fit.tails["gaussian_area"] = self._gaussian_area(quiet, fit, counts, edges)
                        kept[key] = fit
                except ValueError:
                    pass
        return kept[key]

    @staticmethod
    def _gaussian_area(s: Settings, fit: PeakFit, counts, edges) -> float:
        """The area of the Gaussian on a line (:func:`fit_peak`) in the window of ``fit``."""
        return fit_peak(counts, edges, fit.centroid, fit.sigma, s.fit_half_width).area

    def _of(self, gammas) -> Analysis:
        """This analysis with other γ rays of the same setup, a weighted sample."""
        import copy

        a = copy.copy(self)
        a.gammas, a.events, a.taken, a.history, a._whole = gammas, gammas.events, None, [], True
        a._own = (gammas["nucleus"] == 0) if "nucleus" in gammas.columns else np.ones(len(gammas), dtype=bool)
        for name in ("_forward", "_simulated", "_variance", "_gate_overlap", "_gate_efficiency"):
            a.__dict__.pop(name, None)
        return a

    def _simulated_gammas(self, s: Settings, nucleus: int):
        """The γ rays of :meth:`simulated_peak`: the particles that the gate's detectors and rings count, from
        :data:`SHAPE_EVENTS` reactions, drawn in proportion to their weights so that each γ ray counts about the
        same, as many as :data:`SHAPE_COUNTS` full-energy counts need. Kept for the last rings asked for."""
        from .events import simulate
        from .gamma_events import simulate_gammas

        key = json.dumps([s.detectors, s.rings, nucleus], default=_plain)
        kept = self.__dict__.setdefault("_shape_gammas", {})
        if key in kept:
            return kept[key]
        kept.clear()
        g = self.gammas
        seed = SHAPE_SEED + (1 if int(g.seed) == SHAPE_SEED else 0)
        ev = simulate(self.experiment, SHAPE_EVENTS, seed)
        a = self._of(replace(g, events=ev))
        rows = a._elastic_rows(s)
        if nucleus == 0:
            rows &= np.isin(ev.columns["channel"], [i for i, lab in enumerate(ev.channels) if "excited" in lab])
        rows = np.flatnonzero(rows)
        w = np.asarray(ev.columns["weight"], dtype=float)[rows]
        if not w.sum() > 0:
            kept[key] = None
            return None
        energy = self.e0 if nucleus == 0 else self._other_setup()["energy"]
        eff = sum(float(c.response.peak_efficiency(energy, c.element)) for c in g.crystals)
        n = int(min(SHAPE_COUNTS / max(eff, 1e-9), SHAPE_EMISSIONS) / SHAPE_COPIES) + 1
        rng = np.random.default_rng(seed)
        pick = np.sort(rng.choice(rows, size=n, p=w / w.sum()))
        cols = {k: np.asarray(v)[pick] for k, v in ev.columns.items()}
        cols["weight"] = np.full(n, float(w.sum()) / n)
        cols["event"] = np.arange(n, dtype=np.asarray(ev.columns["event"]).dtype)
        kept[key] = simulate_gammas(self.experiment, seed=seed, particle_events=replace(ev, columns=cols),
                                    gammas_per_event=SHAPE_COPIES, rates=self.rates, other=nucleus == 1,
                                    distance_um=g.distance_um)
        return kept[key]

    def _tails_uncertainty(self, fit: PeakFit) -> float:
        """What the simulated share of the peak's tails may be wrong by, relative to the area: the change the share
        makes (the area with it against that of the Gaussian, or, if larger, what the Gaussian loses of the simulated
        peak), times how far the simulated width is from the measured one (and the uncertainty of that), between
        :data:`SHAPE_KNOWN` and 1."""
        t = fit.tails
        if t is None or not t.get("simulated") or not fit.area > 0:
            return 0.0
        change = max(abs(1 - t["gaussian_area"] / fit.area) if t.get("gaussian_area") else 0.0, t["simulated_loss"])
        off = math.hypot(fit.sigma / t["simulated_sigma"] - 1, 1 / math.sqrt(2 * fit.area))
        return change * min(max(off, SHAPE_KNOWN), 1.0)

    def _blends(self, nucleus: int, energy_kev: float, fit: PeakFit) -> list:
        """The other transitions of the same nucleus inside the ±3σ window the area is counted in."""
        half = 3e3 * fit.sigma
        return [f"{ln['label']} ({ln['energy_kev']:.1f} keV)" for ln in self._lines
                if ln["nucleus"] == nucleus and abs(ln["energy_kev"] - energy_kev) > 0.05
                and abs(ln["energy_kev"] - 1e3 * fit.centroid) <= half]

    def _own_line(self, s: Settings, fit: PeakFit, area: float, area_unc: float, eff: float, corr: float,
                  gate_eff: float, yield_: float, elastic_counts: float, spectrum: dict, search: dict) -> dict:
        """The ``[reaction]`` transition's peak as an entry of :attr:`Result.lines`."""
        e_kev = 1e3 * self.e0
        label = next((ln["label"] for ln in self._lines if ln["nucleus"] == 0
                      and abs(ln["energy_kev"] - e_kev) < 1.0), f"{self.gammas.emitter} {e_kev:.1f} keV")
        gate = (f"the {s.particle_gate} group" if s.particle_gate == "inelastic" else "every particle") \
            + (", all rings" if s.rings is None else f", rings {list(s.rings)}")
        return {"nucleus": 0, "nuclide": self.gammas.emitter, "label": label, "energy_kev": e_kev,
                "correction": spectrum["correction"], "gate": gate, "spectrum": spectrum, "fit": asdict(fit),
                "area": area, "area_unc": area_unc, "randoms": fit.randoms, "efficiency": eff, "correlation": corr,
                "gate_efficiency": gate_eff, "yield": yield_,
                "per_particle": yield_ / elastic_counts if elastic_counts > 0 else float("nan"),
                "expected_per_particle": None, "ratio": None,
                "found": self._found(search, spectrum["correction"], 1e3 * fit.centroid,
                                     max(1.5e3 * FWHM_PER_SIGMA * fit.sigma, 2.0)),
                "blends": self._blends(0, e_kev, fit)}

    def _other_setup(self) -> Optional[dict]:
        """What the analysis needs of the other nucleus of the collision, once: its first transition, the
        correction that is right for it, its angular-correlation factors, and the transition's yield per
        collision against the CM angle from its level scheme (to all orders). None when the events do not hold
        its γ rays."""
        if self._other is None:
            em = next((x for x in self._emitters if x["nucleus"] == 1), None)
            # A run taken before the events held both nuclei has no γ ray of the other one: nothing to measure.
            if em is None or self._own.all():
                self._other = {}
            else:
                sol = em["solution"]
                levels = sol.scheme.levels
                lev_f = sol.ground
                lev_i = min((n for n in range(len(levels)) if n != lev_f), key=lambda n: levels[n].energy.value)
                grid, states = sol.table()
                line = next(ln for ln in self._lines if ln["nucleus"] == 1 and ln["initial"] == lev_i
                            and ln["final"] == lev_f)
                key = max(self._layers[0].nuclides.items(), key=lambda kv: kv[1])[0]
                channel = f"{data.nuclide(key).name} ({self._layers[0].name})"
                self._other = {
                    "line": line, "energy": 1e-3 * line["energy_kev"], "correction": em["correction"],
                    "correlation": {(r["particle_detector"], r["gamma_detector"]): r["detector_factor"]
                                    for r in correlation_table(self.experiment, nucleus="other")},
                    "grid": np.r_[0.0, grid],
                    "yield": np.r_[0.0, [st.gamma_yield(lev_i, lev_f) for st in states]],
                    "channels": [i for i, name in enumerate(self.events.channels)
                                 if name == channel or name.startswith(channel + ",")]}
        return self._other or None

    def _other_spectrum(self, s: Settings, o: dict) -> tuple:
        """The spectrum of the other nucleus's line (``o`` from :meth:`_other_setup`), as :meth:`_spectrum` is that of
        the ``[reaction]`` one: (counts, edges, random coincidences expected, a first idea of σ, the variances)."""
        g = self.gammas
        e_line = o["energy"]
        key = {"projectile": "corrected_projectile", "recoil": "corrected_recoil"}[o["correction"]]
        keep = self._in_rings(s)
        sigma = max(float(self._detector_fwhm(e_line)[0]) / FWHM_PER_SIGMA, 1e-4)
        width = sigma / 4
        half = max(12 * sigma, 0.05 * e_line)
        edges = np.arange(e_line - half, e_line + half + width, width)
        t = g.events.beam_time_s * g.live_fraction
        energies = self._energies(s, key)
        counts, _ = np.histogram(energies[keep], bins=edges, weights=g["weight"][keep] * t)
        sigma = min(max(_half_maximum_sigma(counts, edges, e_line), sigma), 0.03 * e_line)
        rnd = np.zeros(len(edges) - 1)
        if s.subtract_randoms:
            rnd = self._random_rate(s, edges) * t
            counts = counts + np.random.default_rng(g.seed + 104729 + 17).poisson(rnd)
        return counts, edges, rnd, sigma, self._bin_variance(energies, keep, edges, rnd)

    def _other_line(self, s: Settings, elastic_counts: float, search: dict) -> Optional[dict]:
        """The first transition of the other nucleus, measured like the ``[reaction]`` one: its peak in the
        spectrum corrected for that nucleus, in coincidence with every particle counted in the gate's rings,
        fitted and integrated, its yield per particle, and what its level scheme gives. An entry of
        :attr:`Result.lines`, or None when the events do not hold that nucleus."""
        o = self._other_setup()
        if o is None:
            return None
        e_line, which = o["energy"], o["correction"]
        keep = self._in_rings(s)
        counts, edges, rnd, sigma, var = self._other_spectrum(s, o)
        try:
            fit = self._peak(s, 1, e_line, counts, edges, sigma, var)
        except ValueError:
            return None
        area, area_unc = self._area(s, fit, edges, rnd)
        eff, corr = self._efficiency_and_correlation(keep, e_line, o["correlation"], ~self._own, which)
        yield_ = area / (eff * corr) if eff > 0 and corr > 0 else float("nan")
        # What the level scheme gives per particle in the same rings: the transition's yield per
        # collision at each particle's CM angle (no excitation on the other nuclides of the target stack).
        elastic = self._elastic_rows(s)
        c = self.events.columns
        y = np.interp(c["theta_cm"][elastic], o["grid"], o["yield"])
        y = np.where(np.isin(c["channel"][elastic], o["channels"]), y, 0.0)
        expected = float(np.average(y, weights=c["weight"][elastic])) if elastic.any() else 0.0
        per_particle = yield_ / elastic_counts if elastic_counts > 0 else float("nan")
        return {"nucleus": 1, "nuclide": o["line"]["nuclide"], "label": o["line"]["label"],
                "energy_kev": o["line"]["energy_kev"], "correction": which,
                "gate": "every particle" + (", all rings" if s.rings is None else f", rings {list(s.rings)}"),
                "spectrum": {"edges": edges.tolist(), "counts": np.asarray(counts, dtype=float).tolist(),
                             "randoms": rnd.tolist(), "correction": which, "energy": float(e_line)},
                "fit": asdict(fit), "area": area, "area_unc": area_unc, "randoms": fit.randoms,
                "efficiency": eff, "correlation": corr, "gate_efficiency": 1.0, "yield": yield_,
                "per_particle": per_particle, "expected_per_particle": expected,
                "ratio": per_particle / expected if expected > 0 else float("nan"),
                "found": self._found(search, which, 1e3 * fit.centroid, max(1.5e3 * FWHM_PER_SIGMA * fit.sigma, 2.0)),
                "blends": self._blends(1, o["line"]["energy_kev"], fit)}

    # -- B(E2) to all orders ---------------------------------------------------------------------------------------

    def _scheme_yield(self, b_up: float, level: int, lo: float, hi: float) -> tuple:
        """The γ rays of every transition per collision against the CM angle, from the level scheme's coupled
        equations with B(E2↑) of the transition from the ground state to ``level`` set to ``b_up`` (e²fm⁴) and
        every other matrix element as the scheme has it: (angles from ``lo`` to ``hi`` in steps of 2°, {(initial,
        final): yields}). Kept: the same question comes back when only the gate changes."""
        import copy

        from .coupled import CoupledChannels
        from .gamma import EVENT_TOLERANCE

        sol = self._emitters[0]["solution"]
        key = (round(b_up, 4), level, round(lo, 3), round(hi, 3))
        kept = self.__dict__.setdefault("_scheme_yields", {})
        if key not in kept:
            scheme = copy.deepcopy(sol.scheme)
            old = scheme.matrix_element(sol.ground, level, "E2")
            spin = float(scheme.levels[sol.ground].spin)
            size = math.sqrt(b_up * (2 * spin + 1))
            scheme.set_matrix_element(sol.ground, level, "E2", math.copysign(size, old.value if old else 1.0),
                                      source="user")
            cc = CoupledChannels(sol.beam, sol.target, sol.beam_energy, scheme, excite=sol.excite,
                                 tolerance=EVENT_TOLERANCE)
            grid = np.arange(max(2.0, 2.0 * math.floor(lo / 2.0)), min(180.0, 2.0 * math.ceil(hi / 2.0)) + 1e-9, 2.0)
            if len(grid) < 2:
                grid = np.array([max(grid[0] - 2.0, 2.0), grid[0]]) if grid[0] > 2.0 else np.array([2.0, 4.0])
            states = [cc._populated(float(t), cc._densities(a)) for t, a in zip(grid, cc.solve_many(grid))]
            if len(kept) > 24:
                kept.clear()
            kept[key] = (grid, {(i, f): np.array([st.gamma_yield(i, f) for st in states])
                                for i, f in cc.transitions()})
        return kept[key]

    def _all_orders(self, s: Settings, b_first: float, p_per_b: float, fit: PeakFit) -> Optional[dict]:
        """B(E2↑) of the ``[reaction]`` transition from the measured yield with the level scheme's coupled
        equations in place of first-order theory: the value for which they give, per particle in the
        gate's rings, the γ rays that first order would need ``b_first`` for, counting every transition of the
        scheme inside the peak's ±3σ window. None when first order is all there is (one state, first-order
        events, ``Settings.order = "first"``) or when ``b_first`` is no number."""
        from .coupled import CoupledChannels
        from .gamma import excitation_order, state_transition

        sol = self._emitters[0]["solution"] if self._emitters else None
        if (s.order == "first" or not isinstance(sol, CoupledChannels) or excitation_order(self.experiment) != "all"
                or not b_first == b_first or b_first <= 0 or p_per_b <= 0):
            return None
        level, ground = state_transition(self.experiment)
        own = next((ln for ln in self._lines if ln["nucleus"] == 0 and (ln["initial"], ln["final"]) == (level, ground)),
                   None)
        if own is None:
            return None
        inside = [ln for ln in self._lines if ln["nucleus"] == 0
                  and abs(ln["energy_kev"] - 1e3 * fit.centroid) <= 3e3 * fit.sigma]
        if own not in inside:
            inside.append(own)
        lines = tuple(sorted((ln["initial"], ln["final"]) for ln in inside))
        rows = self._elastic_rows(s)
        if not rows.any():
            return None
        c = self.events.columns
        th, w = c["theta_cm"][rows], c["weight"][rows]
        measured = b_first * p_per_b  # γ rays per particle: what the yield and the normalisation gave
        lo, hi = float(th.min()), float(th.max())
        # The scheme's solution is at the middle of the target and its yields are taken at the particles' CM
        # angles; the forward model excites through the depth, with the excited particle's own kinematics. The
        # ratio of the two for the excitations of the scheme as it is set (the forward model's excitations per
        # particle over the solution's at the particles' angles) carries the difference over to the yields.
        mid = float(np.average(SchemeProbability(sol).probability(th), weights=w))
        forward = self._mean_probability(s, mode="forward")
        kappa = forward / mid if mid > 0 and forward > 0 else 1.0

        def expected(b_up: float, which: tuple = lines) -> float:
            grid, y = self._scheme_yield(b_up, level, lo, hi)
            return kappa * float(np.average(np.interp(th, grid, sum(y[k] for k in which if k in y)), weights=w))

        # The yield is nearly proportional to B(E2): one pass by proportion, then the straight line through the
        # two points, which is good to a few parts in ten thousand; a third pass only if the two were far apart.
        b0, f0 = b_first, expected(b_first)
        if not f0 > 0:
            return None
        b1 = b0 * measured / f0
        f1 = expected(b1)
        if not f1 > 0:
            return None
        b = b1 + (measured - f1) * (b1 - b0) / (f1 - f0) if f1 != f0 else b1
        if abs(b - b1) > 0.02 * b1:
            b2, f2 = b, expected(b)
            b = b2 + (measured - f2) * (b2 - b1) / (f2 - f1) if f2 != f1 else b2
        # How much B(E2) changes for a change of the yield: more than in proportion, since a strong excitation
        # saturates (the yield grows more slowly than B(E2)). Every uncertainty of the yield grows by it.
        slope = math.log(f1 / f0) / math.log(b1 / b0) if b1 != b0 and f1 != f0 else 1.0
        gain = min(max(1.0 / slope, 1.0), 4.0) if slope > 0 else 1.0
        others = [ln for ln in inside if ln is not own]
        levels = sol.scheme.levels
        fed = any(lev.energy.value > levels[level].energy.value for lev in levels)
        own_share = expected(b1, ((level, ground),)) / f1
        text = (f"first order would read {b_first:.4g} e²fm⁴ from this yield. The level scheme's coupled equations "
                f"give the measured {measured:.4g} γ rays per particle in the window for B(E2↑) = {b:.4g} "
                f"e²fm⁴ ({100 * (b / b_first - 1):+.1f}%)"
                + (": they count the feeding of the state from the levels above" if fed else "")
                + ((", and " if fed else ": they count ") + ", ".join(f"{ln['label']} ({ln['energy_kev']:.1f} keV)" for ln in others)
                   + f", which the window also holds ({100 * (1 - own_share):.1f}% of its γ rays)" if others else "")
                + ". Every other matrix element is taken from the scheme."
                + (f" The scheme's yields, solved at the middle of the target, are scaled by {kappa:.4f}: the "
                   "forward model's excitations per particle (through the depth of the target, with the excited "
                   "particle's kinematics) over the solution's." if abs(kappa - 1) > 1e-3 else "")
                + (f" The excitation is strong enough to saturate: the yield grows as B(E2) to the power "
                   f"{1 / gain:.2f}, so a 1% change of the yield is a {gain:.2f}% change of B(E2), and every "
                   "uncertainty of the yield counts that much more." if gain > 1.02 else ""))
        return {"b": b, "text": text, "lines": lines, "gain": gain}

    def _monte_carlo_uncertainty(self, s: Settings, keep: np.ndarray, fit: PeakFit) -> float:
        """The relative uncertainty of the result that comes from the finite Monte Carlo sample: the weights of
        the γ rays in the peak window, and of the particles it is normalised to."""
        g = self.gammas
        exc = self.experiment.excitation
        which = s.correction if s.correction != "emitter" else ("recoil" if exc.excite == "target" else "projectile")
        x = g[{"projectile": "corrected_projectile", "recoil": "corrected_recoil"}[which]]
        inside = keep & (np.abs(x - fit.centroid) < 3 * fit.sigma)
        w = g["weight"][inside]
        peak = math.sqrt(float(np.sum(w**2))) / float(np.sum(w)) if w.sum() > 0 else float("inf")
        el = self.events.columns["weight"][self._elastic_rows(s)]
        elastic = math.sqrt(float(np.sum(el**2))) / float(np.sum(el)) if el.sum() > 0 else 0.0
        return math.hypot(peak, elastic)

    def _efficiency_and_correlation(self, keep: np.ndarray, energy: Optional[float] = None,
                                    correlation: Optional[dict] = None, own: Optional[np.ndarray] = None,
                                    correction: Optional[str] = None) -> tuple:
        """The full-energy-peak efficiency of all crystals together (the spectrum is summed), each crystal's at the
        energy it sees the line at (Doppler shifted), and the angular-correlation factor averaged over the
        excitations in the gate: for each excitation, the efficiency-weighted mean of its particle detector's
        factors. For the ``[reaction]`` transition, or for another one at ``energy`` (MeV) with its
        ``correlation`` factors, the γ rays ``own`` of its nucleus and the Doppler ``correction`` right for it."""
        g = self.gammas
        exc = self.experiment.excitation
        if correction is None:
            correction = "recoil" if exc.excite == "target" else "projectile"
        energy = self.e0 if energy is None else energy
        correlation = self._correlation if correlation is None else correlation
        mine = keep & (self._own if own is None else own)
        rows = np.unique(g["particle"][mine])
        c = self.events.columns
        # Each crystal sees the line at its Doppler-shifted energy, and its efficiency is that at this energy: the
        # shift is read from the data, as the measured energy over the Doppler-corrected one of the γ rays in the
        # gate that the crystal counted (the correction for the nucleus that emits the line).
        measured = np.asarray(g["measured"], dtype=float)
        fixed = np.asarray(g[{"projectile": "corrected_projectile", "recoil": "corrected_recoil"}[correction]],
                           dtype=float)
        good = mine & g["counted"] & (measured > 0) & (fixed > 0)
        crystal = np.asarray(g["crystal"])
        shift = np.ones(len(g.crystals))
        if good.any():
            wsum = np.bincount(crystal[good], weights=g["weight"][good], minlength=len(g.crystals))
            ratio = np.bincount(crystal[good], weights=g["weight"][good] * measured[good] / fixed[good],
                                minlength=len(g.crystals))
            shift = np.where(wsum > 0, ratio / np.where(wsum > 0, wsum, 1.0), 1.0)
        effs = {}
        for k, cr in enumerate(g.crystals):
            name = cr.name.rsplit(" ", 1)[0] if cr.element is not None else cr.name
            effs[name] = effs.get(name, 0.0) + float(cr.response.peak_efficiency(energy * shift[k], cr.element))
        total_eff = sum(effs.values())
        per_det = {d: sum(e * correlation.get((d, n), 1.0) for n, e in effs.items()) / total_eff
                   for d in self.events.detectors}
        factors = np.array([per_det[self.events.detectors[d]] for d in c["detector"][rows]])
        w = c["weight"][rows]
        return float(total_eff), float(np.average(factors, weights=w)) if len(w) else 1.0

    def _reference(self, s: Settings, keep: np.ndarray, w: np.ndarray) -> dict:
        """The reference transition's expected peak counts (analytic, with Poisson noise) and its ⟨P⟩/B."""
        ref = s.reference
        energy = _q(ref["energy"]).to("MeV")
        from .coulex import parse_b

        b_val = parse_b(ref["b_up"], 2)
        exc = self.experiment.excitation
        ion = beam_ion(self.experiment)
        target = data.nuclide(max(self._layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        other = "projectile" if exc.excite == "target" else "target"
        cx = Coulex(ion, target, self._e_mid, excite=other, energy=energy, multipolarity="E2", b_up=b_val)
        elastic_rows = self._elastic_rows(s)
        th = self.events.columns["theta_cm"][elastic_rows]
        p = float(np.average(np.asarray(cx.probability(th)), weights=self.events.columns["weight"][elastic_rows]))
        _, elastic = self._gate(s)
        eff = sum(float(c.response.peak_efficiency(energy, c.element)) for c in self.gammas.crystals)
        expected = elastic * p * eff
        rng = np.random.default_rng(self.gammas.seed + 7)
        counts = float(rng.poisson(expected)) if expected < 1e9 else expected
        # The reference's yield, like the unknown's, is its peak over its efficiency (its angular correlation is
        # taken as isotropic here, since the simulation does not excite it).
        return {"energy": energy, "b_up": b_val, "p_per_b": p / b_val if b_val else 0.0, "counts": counts / eff,
                "raw_counts": counts}

    def _budget(self, s: Settings, keep: np.ndarray, b: float, p_per_b: float, eff: float, corr: float,
                ref_unc: float, fits: list) -> dict:
        """Relative systematic contributions, each by changing its input by one standard deviation. ``fits``: the
        fit of the peak, and that of the peak it is normalised to."""
        if not b == b or b <= 0:
            return {}
        out = {"efficiency": s.efficiency_unc, "angular correlation": s.correlation_unc}
        # Beam energy: ⟨P⟩/B changes with the orbit. The energy of the beam that enters the target is raised by its
        # uncertainty, and it slows down from there; both sides with the first-order P at each energy (the forward
        # model's P is only tabulated through the energies the nominal beam has).
        base = self._mean_probability(s, keep, exact=True)
        shifted = self._mean_probability(s, keep, energy_scale=1 + s.beam_energy_unc, exact=True)
        out["beam energy"] = abs(shifted / base - 1) if base > 0 else 0.0
        # The target's place: the Doppler correction with the target moved by its uncertainty shifts and broadens
        # the peak, and the γ-ray efficiency goes as 1/d²; the peak area is found again with the moved target.
        d_mm = float(np.mean([gd.distance_mm() for gd in self.experiment.gamma_detectors]))
        out["detector positions"] = math.hypot(2 * s.position_unc_mm / d_mm,
                                               self._area_change(s, keep, s.position_unc_mm, fits[0]))
        tails = [self._tails_uncertainty(f) for f in fits]
        if any(tails):
            out["peak tails (simulated)"] = math.sqrt(sum(x**2 for x in tails))
        if s.normalisation == "target":
            out["reference B(E2)"] = ref_unc  # or, for a measured reference peak, its matrix elements
        out["matrix elements assumed"] = 0.0  # first order: no other matrix element enters
        return out

    def _area_change(self, s: Settings, keep: np.ndarray, delta_mm: float, like: PeakFit) -> float:
        """The relative change of the peak area when the target assumed in the correction moves by ±delta, each area
        from a fit of the kind ``like`` is (with its own tails, the simulated ones, or none)."""
        exc = self.experiment.excitation
        which = s.correction if s.correction != "emitter" else ("recoil" if exc.excite == "target" else "projectile")
        key = {"projectile": "corrected_projectile", "recoil": "corrected_recoil"}[which]
        g = self.gammas
        t = g.events.beam_time_s * g.live_fraction
        sigma = max(float(np.mean([c.response.fwhm(self.e0) for c in g.crystals])) / FWHM_PER_SIGMA, 1e-4)
        areas = []
        for offset in (s.target_offset_mm, s.target_offset_mm + delta_mm, s.target_offset_mm - delta_mm):
            x = self._energies(s, key, offset)
            half = max(12 * sigma, 0.05 * self.e0)
            edges = np.arange(self.e0 - half, self.e0 + half + sigma / 4, sigma / 4)
            counts, _ = np.histogram(x[keep], bins=edges, weights=g["weight"][keep] * t)
            guess = min(max(_half_maximum_sigma(counts, edges, self.e0), sigma), 0.03 * self.e0)
            var = self._bin_variance(x, keep, edges, np.zeros(len(edges) - 1))
            try:
                fit = self._peak(s, 0, self.e0, counts, edges, guess, var, like=like)
                areas.append(fit.area if not fit.failed else float("nan"))
            except ValueError:
                areas.append(float("nan"))
        if not areas[0] or not all(a == a for a in areas):
            return 0.0
        return float(max(abs(areas[1] - areas[0]), abs(areas[2] - areas[0])) / areas[0])

    def _shape(self, b: float, a_mass: int) -> dict:
        z = data.nuclide(ex_nuclide(self.experiment)).Z
        if not b == b or b <= 0:
            return {}
        r0 = 1.2 * a_mass ** (1 / 3)
        beta2 = 4 * math.pi / (3 * z * r0**2) * math.sqrt(b)
        q0 = math.sqrt(16 * math.pi * b / 5)
        out = {"beta2": beta2, "q0_efm2": q0, "qs_2plus_efm2": -2 / 7 * q0,
               "note": "β₂, Q₀ and Q_s(2⁺) assume a rigid axially symmetric rotor; the sign of Q_s is not "
                       "measured here (item 58)."}
        scheme = self.experiment.levels.get("target" if self.experiment.excitation.excite == "target" else "beam")
        if scheme is not None:
            twos = [lev for lev in scheme.levels if lev.spin == 2 and lev.energy.value > 0]
            fours = [lev for lev in scheme.levels if lev.spin == 4]
            if twos and fours:
                out["e4_over_e2"] = fours[0].energy.value / twos[0].energy.value
        return out


def ex_nuclide(experiment) -> str:
    """The nucleus the reaction excites."""
    exc = experiment.excitation
    if exc.excite == "projectile":
        return beam_ion(experiment)
    layers = stack(experiment)
    return data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name


__all__ = ["Analysis", "PeakFit", "RESULT_FILE", "Result", "SHIFT_H", "Settings", "ex_nuclide", "fit_peak", "fit_tails",
           "load_result", "save_result", "tail_shape"]
