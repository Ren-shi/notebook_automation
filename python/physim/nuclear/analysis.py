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
   ``gate_width`` of the first. Rings (or strips) can be chosen.
2. **Doppler correction** for the emitting nucleus (or the other one, to see what happens).
3. **Peak search:** each corrected spectrum (for the projectile, and for the recoil) is searched for peaks over
   its whole range (:func:`find_peaks`), each peak found is fitted, and those that sit on a transition of the
   nucleus that correction is right for are named (:func:`identify`). The other nucleus's lines are smeared by
   the wrong correction: each nucleus shows its peaks in its own spectrum.
4. **Peak fit:** a Gaussian on a straight line, by weighted least squares within ``fit_half_width`` resolutions
   of the transition energy, with the random coincidences expected in the window subtracted.
5. **Yield:** the area over the full-energy-peak efficiency, the γ-ray branch (internal conversion) and the
   angular-correlation factor.
6. **Normalisation**, chosen by the user: ``"rutherford"``, to the particles counted in the same rings,
   so the beam current and the target thickness cancel; or ``"target"``, to a transition of the other nucleus of
   the collision whose strength is known. When the events hold that nucleus's γ rays (the setup has its level
   scheme), its peak is found, fitted and integrated in its own corrected spectrum like the first one, and its
   yield is compared with what its level scheme gives to all orders (:attr:`Result.lines`); otherwise the
   reference is ``Settings.reference``, its peak taken from its analytic rate with Poisson noise (noted in the
   result).
7. **B(E2)** from the first-order proportionality of the excitation probability to B(E2), in e²fm⁴, e²b² and
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
from dataclasses import asdict, dataclass, field, fields
from typing import Optional

import numpy as np

from . import data
from .coulex import Coulex
from .gamma import correlation_table, emitters, transition_lines
from .kinematics import TwoBody
from .levels import weisskopf_unit
from .quantity import Quantity
from .rates import Rates, _q, beam_energy_at, beam_ion, exit_energy, stack
from .response import FWHM_PER_SIGMA, Response

#: Hours in a shift, for the counts per shift.
SHIFT_H = 8.0


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

    def curve(self, energy) -> np.ndarray:
        """The fitted Gaussian on its line at ``energy`` (MeV), counts per bin."""
        e = np.asarray(energy, dtype=float)
        return (self.height * np.exp(-0.5 * ((e - self.centroid) / self.sigma) ** 2) + self.level
                + self.slope * (e - self.origin))

    @property
    def fwhm_kev(self) -> float:
        return 1e3 * FWHM_PER_SIGMA * self.sigma


def fit_peak(counts: np.ndarray, edges: np.ndarray, energy: float, sigma0: float, half_width: float = 4.0,
             iterations: int = 40, max_sigma: Optional[float] = None, variance=None) -> PeakFit:
    """Fit a Gaussian on a straight line to ``counts`` within ``half_width`` × FWHM of ``energy``, by weighted
    least squares (Gauss–Newton; the weights are the counts, at least one per bin), then count the peak within
    ±3σ of its centre over the line. Returns the area in counts.

    ``variance`` is the variance of each bin when it is not its counts: a spectrum of weighted events scaled to
    a beam time holds far more counts than it has entries, and a fit that took the counts for their own variance
    would trust every bin too much and follow the scatter of the few entries. Only the fit's weights use it; the
    uncertainty of the area stays that of the counts, which is what the beam time gives."""
    centres = (edges[:-1] + edges[1:]) / 2
    width = float(edges[1] - edges[0])
    half = half_width * FWHM_PER_SIGMA * sigma0
    m = np.abs(centres - energy) <= half
    x, y = centres[m] - energy, counts[m].astype(float)
    if m.sum() < 6:
        raise ValueError("the fit window holds too few bins; widen it or use finer bins")
    if variance is None:
        w = 1 / np.maximum(y, 1.0)
    else:
        var = np.asarray(variance, dtype=float)[m]
        # A bin without entries is as uncertain as one entry (as a bin without counts is taken as one count).
        filled = var[var > 0]
        w = 1 / np.maximum(var, float(np.median(filled)) if filled.size else 1.0)
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
        p[2] = min(max(abs(p[2]), 0.3 * sigma0), max_sigma if max_sigma else 4.0 * sigma0)
        p[1] = min(max(p[1], -2 * sigma0), 2 * sigma0)
        p[0] = max(p[0], 0.0)
        if np.max(np.abs(step[:3]) / (np.abs(p[:3]) + 1e-12)) < 1e-7:
            break
    r = y - model(p)
    chi2 = float(np.sum(w * r * r))
    gauss_area = p[0] * p[2] * math.sqrt(2 * math.pi) / width
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
                   origin=float(energy))


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

    @property
    def total_unc(self) -> float:
        return math.hypot(self.statistical, self.systematic)

    def to_dict(self) -> dict:
        """The result as plain data, to keep with a run (:func:`save_result`)."""
        return json.loads(json.dumps(asdict(self), default=_plain))

    @classmethod
    def from_dict(cls, d: dict) -> Result:
        d = {k: v for k, v in d.items() if k in {f.name for f in fields(cls)}}
        fit = d.get("fit") or {}
        d["fit"] = PeakFit(**{k: v for k, v in fit.items() if k in {f.name for f in fields(PeakFit)}})
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
            keep &= near & closer
        # What the gate does, known here because the simulation knows the truth: the share of the excitations it
        # keeps, and the particles it lets through (they carry no γ ray, so they only add to the randoms).
        in_dets = self._in_rings(s)
        w_all = g["weight"][in_dets & self._own].sum()
        self._gate_efficiency = float(g["weight"][keep & self._own].sum() / w_all) if w_all > 0 else 1.0
        return keep, float(np.sum(c["weight"][elastic]) * ev.beam_time_s * g.live_fraction)

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
        the way out of the target (from its middle, or from its back) and the detector's dead layer, as the
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
            z = self._layers[0].thickness * (0.5 if depth == "middle" else 0.999)
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

    def _mean_probability(self, s: Settings, keep: np.ndarray, beam_energy: Optional[float] = None,
                          energy: Optional[float] = None, b_up: Optional[float] = None) -> float:
        """The first-order excitation probability per unit B(E2↑), averaged over the particles in the gate
        (weighted as the events are), at the given beam energy and state."""
        exc = self.experiment.excitation
        e_beam = beam_energy if beam_energy is not None else self._e_mid
        e_star = energy if energy is not None else self.e0
        b = b_up if b_up is not None else exc.b_up_e2fm
        ion = beam_ion(self.experiment)
        target = data.nuclide(max(self._layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name
        cx = Coulex(ion, target, e_beam, excite=exc.excite, energy=e_star, multipolarity=exc.multipolarity, b_up=b)
        # Averaged over the particles in the gate's rings, as the Rutherford normalisation divides by them.
        elastic = self._elastic_rows(s)
        th = self.events.columns["theta_cm"][elastic]
        w = self.events.columns["weight"][elastic]
        p = np.asarray(cx.probability(th))
        return float(np.average(p, weights=w) / b) if elastic.any() else 0.0

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
        # Two passes: the peak after the Doppler correction is wider than the resolution by what the segments
        # leave, so the window is set from a first fit's width.
        top = 0.03 * self.e0  # no Doppler residue is wider than a few percent
        var = self._variance
        first = fit_peak(counts, edges, self.e0, sigma, 3 * s.fit_half_width, max_sigma=top, variance=var)
        fit = fit_peak(counts, edges, self.e0, first.sigma, s.fit_half_width, max_sigma=top, variance=var)
        window = np.abs((edges[:-1] + edges[1:]) / 2 - fit.centroid) < 3 * fit.sigma
        fit.randoms = float(rnd[window].sum()) if s.subtract_randoms else 0.0
        area = fit.area - fit.randoms
        area_unc = math.hypot(fit.area_unc, math.sqrt(max(fit.randoms, 0.0)))
        steps.append({"step": "peak fit", "text": f"Gaussian on a line within ±{s.fit_half_width:g} FWHM: area "
                      f"{fit.area:.1f} ± {fit.area_unc:.1f} at {1e3 * fit.centroid:.1f} keV, FWHM {fit.fwhm_kev:.1f} "
                      f"keV, χ² {fit.chi2:.0f} for {fit.bins} bins"})
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
                          f"rings: ⟨P⟩ = {yield_ / norm_counts:.3e}, and ⟨P⟩/B(E2) = {p_per_b:.3e} per e²fm⁴"})
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
            notes.append("No peak: the gate or the rings leave too few γ rays to fit.")
        budget = self._budget(s, keep, b, p_per_b, eff, corr, ref_unc)
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
        if s.particle_gate == "inelastic":
            notes.append(notes_gate)
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
            notes.append(f"The {other['label']} peak ({other['energy_kev']:.1f} keV, corrected for the "
                         f"{other['correction']}): {other['per_particle']:.4g} γ rays per particle "
                         f"measured, {other['expected_per_particle']:.4g} expected from its level scheme to all "
                         f"orders; measured over expected {other['ratio']:.3f}.")
        result = Result(s.as_dict(), steps, fit, area, eff, corr, branch, s.normalisation, norm_counts, b,
                        statistical, systematic, mc, b * 1e-4, b / wu, budget, truth, pull, shape, per_shift, hours,
                        notes, source, whole, spectrum, search, lines, b_first)
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
                if not fit.area > 0:
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

    def _other_line(self, s: Settings, elastic_counts: float, search: dict) -> Optional[dict]:
        """The first transition of the other nucleus, measured like the ``[reaction]`` one: its peak in the
        spectrum corrected for that nucleus, in coincidence with every particle counted in the gate's rings,
        fitted and integrated, its yield per particle, and what its level scheme gives. An entry of
        :attr:`Result.lines`, or None when the events do not hold that nucleus."""
        o = self._other_setup()
        if o is None:
            return None
        g = self.gammas
        e_line, which = o["energy"], o["correction"]
        key = {"projectile": "corrected_projectile", "recoil": "corrected_recoil"}[which]
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
        top = 0.03 * e_line
        var = self._bin_variance(energies, keep, edges, rnd)
        try:
            first = fit_peak(counts, edges, e_line, sigma, 3 * s.fit_half_width, max_sigma=top, variance=var)
            fit = fit_peak(counts, edges, e_line, first.sigma, s.fit_half_width, max_sigma=top, variance=var)
        except ValueError:
            return None
        window = np.abs((edges[:-1] + edges[1:]) / 2 - fit.centroid) < 3 * fit.sigma
        fit.randoms = float(rnd[window].sum()) if s.subtract_randoms else 0.0
        area = fit.area - fit.randoms
        area_unc = math.hypot(fit.area_unc, math.sqrt(max(fit.randoms, 0.0)))
        eff, corr = self._efficiency_and_correlation(keep, e_line, o["correlation"], ~self._own)
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

        def expected(b_up: float, which: tuple = lines) -> float:
            grid, y = self._scheme_yield(b_up, level, lo, hi)
            return float(np.average(np.interp(th, grid, sum(y[k] for k in which if k in y)), weights=w))

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
        own_share = expected(b1, ((level, ground),)) / f1
        text = (f"first order would read {b_first:.4g} e²fm⁴ from this yield. The level scheme's coupled equations "
                f"give the measured {measured:.4g} γ rays per particle in the window for B(E2↑) = {b:.4g} "
                f"e²fm⁴ ({100 * (b / b_first - 1):+.1f}%): they count the feeding of the state from the levels above"
                + (", and " + ", ".join(f"{ln['label']} ({ln['energy_kev']:.1f} keV)" for ln in others)
                   + f", which the window also holds ({100 * (1 - own_share):.1f}% of its γ rays)" if others else "")
                + ". Every other matrix element is taken from the scheme."
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
                                    correlation: Optional[dict] = None, own: Optional[np.ndarray] = None) -> tuple:
        """The full-energy-peak efficiency of all crystals together (the spectrum is summed) and the
        angular-correlation factor averaged over the excitations in the gate: for each excitation, the
        efficiency-weighted mean of its particle detector's factors. For the ``[reaction]`` transition, or for
        another one at ``energy`` (MeV) with its ``correlation`` factors and the γ rays ``own`` of its nucleus."""
        g = self.gammas
        energy = self.e0 if energy is None else energy
        correlation = self._correlation if correlation is None else correlation
        rows = np.unique(g["particle"][keep & (self._own if own is None else own)])
        c = self.events.columns
        effs = {}
        for cr in g.crystals:
            name = cr.name.rsplit(" ", 1)[0] if cr.element is not None else cr.name
            effs[name] = effs.get(name, 0.0) + float(cr.response.peak_efficiency(energy, cr.element))
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
                ref_unc: float) -> dict:
        """Relative systematic contributions, each by changing its input by one standard deviation."""
        if not b == b or b <= 0:
            return {}
        out = {"efficiency": s.efficiency_unc, "angular correlation": s.correlation_unc}
        # Beam energy: ⟨P⟩/B changes with the orbit.
        shifted = self._mean_probability(s, keep, beam_energy=self._e_mid * (1 + s.beam_energy_unc))
        out["beam energy"] = abs(shifted / p_per_b - 1) if p_per_b > 0 else 0.0
        # The target's place: the Doppler correction with the target moved by its uncertainty shifts and broadens
        # the peak, and the γ-ray efficiency goes as 1/d²; the peak area is found again with the moved target.
        d_mm = float(np.mean([gd.distance_mm() for gd in self.experiment.gamma_detectors]))
        out["detector positions"] = math.hypot(2 * s.position_unc_mm / d_mm,
                                               self._area_change(s, keep, s.position_unc_mm))
        if s.normalisation == "target":
            out["reference B(E2)"] = ref_unc  # or, for a measured reference peak, its matrix elements
        out["matrix elements assumed"] = 0.0  # first order: no other matrix element enters
        return out

    def _area_change(self, s: Settings, keep: np.ndarray, delta_mm: float) -> float:
        """The relative change of the peak area when the target assumed in the correction moves by ±delta."""
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
                first = fit_peak(counts, edges, self.e0, guess, 3 * s.fit_half_width, max_sigma=0.03 * self.e0,
                                 variance=var)
                areas.append(fit_peak(counts, edges, self.e0, first.sigma, s.fit_half_width,
                                      max_sigma=0.03 * self.e0, variance=var).area)
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


__all__ = ["Analysis", "PeakFit", "RESULT_FILE", "Result", "SHIFT_H", "Settings", "ex_nuclide", "fit_peak",
           "load_result", "save_result"]
