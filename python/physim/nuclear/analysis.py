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
3. **Peak fit:** a Gaussian on a straight line, by weighted least squares within ``fit_half_width`` resolutions
   of the transition energy, with the random coincidences expected in the window subtracted.
4. **Yield:** the area over the full-energy-peak efficiency, the γ-ray branch (internal conversion) and the
   angular-correlation factor.
5. **Normalisation**, chosen by the user: ``"rutherford"``, to the elastic particles counted in the same rings,
   so the beam current and the target thickness cancel; or ``"target"``, to a reference transition of known
   B(E2) (``Settings.reference``), whose peak is taken from its analytic rate with Poisson noise, since the
   simulation excites one state (noted in the result).
6. **B(E2)** from the first-order proportionality of the excitation probability to B(E2), in e²fm⁴, e²b² and
   Weisskopf units, with the statistical uncertainty of the fit and a budget of systematic contributions, each
   found by changing its input by one standard deviation and running again.
7. **Shape:** β₂, the strength in Weisskopf units, E(4⁺)/E(2⁺) if a level scheme is there, and the quadrupole
   moment a rigid rotor would have, with the warning that these readings depend on the rotor model.

The result also holds the comparison with the value that was put in, and the beam time for a chosen precision.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, fields
from typing import Optional

import numpy as np

from . import data
from .coulex import Coulex
from .gamma import correlation_table
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

    @property
    def fwhm_kev(self) -> float:
        return 1e3 * FWHM_PER_SIGMA * self.sigma


def fit_peak(counts: np.ndarray, edges: np.ndarray, energy: float, sigma0: float, half_width: float = 4.0,
             iterations: int = 40, max_sigma: Optional[float] = None) -> PeakFit:
    """Fit a Gaussian on a straight line to ``counts`` within ``half_width`` × FWHM of ``energy``, by weighted
    least squares (Gauss–Newton; the weights are the counts, at least one per bin), then count the peak within
    ±3σ of its centre over the line. Returns the area in counts."""
    centres = (edges[:-1] + edges[1:]) / 2
    width = float(edges[1] - edges[0])
    half = half_width * FWHM_PER_SIGMA * sigma0
    m = np.abs(centres - energy) <= half
    x, y = centres[m] - energy, counts[m].astype(float)
    if m.sum() < 6:
        raise ValueError("the fit window holds too few bins; widen it or use finer bins")
    w = 1 / np.maximum(y, 1.0)
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
                   gauss_area=float(gauss_area))


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

    @property
    def total_unc(self) -> float:
        return math.hypot(self.statistical, self.systematic)


class Analysis:
    """The analysis of one simulated experiment (``gammas`` from
    :func:`physim.nuclear.gamma_events.simulate_gammas`). Every :meth:`run` is kept in :attr:`history` with
    its settings, so a change and its effect are on record."""

    def __init__(self, experiment, gammas):
        self.experiment = experiment
        self.gammas = gammas
        self.events = gammas.events
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

    # -- the steps --------------------------------------------------------------------------------------------------

    def _gate(self, s: Settings) -> tuple:
        """Mask over the γ rays whose particle passes the gate, and the number of elastic particles in the same
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
        # keeps, and the elastic particles it lets through (they carry no γ ray, so they only add to the randoms).
        in_dets = np.isin(c["detector"][rows], [ev.detectors.index(d) for d in dets]) & g["counted"]
        if s.rings is not None:
            in_dets &= np.isin(c["segment_i"][rows], list(s.rings))
        w_all = g["weight"][in_dets].sum()
        self._gate_efficiency = float(g["weight"][keep].sum() / w_all) if w_all > 0 else 1.0
        return keep, float(np.sum(c["weight"][elastic]) * ev.beam_time_s * g.live_fraction)

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
        return counts, edges, rnd, sigma, which

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
        # Averaged over the elastic particles in the gate's rings, as the Rutherford normalisation divides by them.
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
                      f"{len(keep)} γ rays kept; {elastic_counts:.4g} elastic particles in the same rings"})
        counts, edges, rnd, sigma, which = self._spectrum(s, keep)
        steps.append({"step": "Doppler correction", "text": f"for the {which}"})
        # Two passes: the peak after the Doppler correction is wider than the resolution by what the segments
        # leave, so the window is set from a first fit's width.
        top = 0.03 * self.e0  # no Doppler residue is wider than a few percent
        first = fit_peak(counts, edges, self.e0, sigma, 3 * s.fit_half_width, max_sigma=top)
        fit = fit_peak(counts, edges, self.e0, first.sigma, s.fit_half_width, max_sigma=top)
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
        if s.normalisation == "rutherford":
            norm_counts = elastic_counts
            # Excitations / elastic particles in the same rings = ⟨P⟩ in first order.
            b = yield_ / norm_counts / p_per_b if norm_counts > 0 and p_per_b > 0 else float("nan")
            steps.append({"step": "normalisation", "text": f"to the {norm_counts:.4g} elastic particles in the same "
                          f"rings: ⟨P⟩ = {yield_ / norm_counts:.3e}, and ⟨P⟩/B(E2) = {p_per_b:.3e} per e²fm⁴"})
            ref_unc = 0.0
        else:
            if not s.reference:
                raise ValueError("normalisation to the target needs Settings.reference: the known transition")
            ref = self._reference(s, keep, w)
            norm_counts = ref["counts"]
            b = yield_ / norm_counts * ref["b_up"] * ref["p_per_b"] / p_per_b if norm_counts > 0 else float("nan")
            ref_unc = float(s.reference.get("unc", 0.0))
            steps.append({"step": "normalisation", "text": f"to the reference transition at "
                          f"{1e3 * ref['energy']:.1f} keV, B(E2↑) = {ref['b_up']:.4g} e²fm⁴ ± {100 * ref_unc:.1f}%: "
                          f"{ref['raw_counts']:.4g} counts in its peak, {norm_counts:.4g} over its efficiency"})
            notes.append("The reference peak is the analytic expectation with Poisson noise: the simulation excites "
                         "one state.")
        statistical = math.sqrt((area_unc / area) ** 2 + 1 / max(norm_counts, 1.0)) if area > 0 else float("inf")
        mc = self._monte_carlo_uncertainty(s, keep, fit) if area > 0 else float("inf")
        if not area > 0:
            b = float("nan")
            notes.append("No peak: the gate or the rings leave too few γ rays to fit.")
        budget = self._budget(s, keep, b, p_per_b, eff, corr, ref_unc)
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
        result = Result(s.as_dict(), steps, fit, area, eff, corr, branch, s.normalisation, norm_counts, b,
                        statistical, systematic, mc, b * 1e-4, b / wu, budget, truth, pull, shape, per_shift, hours,
                        notes)
        if self.history:
            changed = {k: v for k, v in s.as_dict().items() if v != self.history[-1].settings.get(k)}
            if changed:
                result.notes.append(f"Changed from the previous run: {changed}; B(E2) was "
                                    f"{self.history[-1].b_e2fm4:.4g} e²fm⁴, now {b:.4g}.")
        self.history.append(result)
        return result

    def _monte_carlo_uncertainty(self, s: Settings, keep: np.ndarray, fit: PeakFit) -> float:
        """The relative uncertainty of the result that comes from the finite Monte Carlo sample: the weights of
        the γ rays in the peak window, and of the elastic particles it is normalised to."""
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

    def _efficiency_and_correlation(self, keep: np.ndarray) -> tuple:
        """The full-energy-peak efficiency of all crystals together (the spectrum is summed) and the
        angular-correlation factor averaged over the excitations in the gate: for each excitation, the
        efficiency-weighted mean of its particle detector's factors."""
        g = self.gammas
        rows = np.unique(g["particle"][keep])
        c = self.events.columns
        effs = {}
        for cr in g.crystals:
            name = cr.name.rsplit(" ", 1)[0] if cr.element is not None else cr.name
            effs[name] = effs.get(name, 0.0) + float(cr.response.peak_efficiency(self.e0, cr.element))
        total_eff = sum(effs.values())
        per_det = {d: sum(e * self._correlation.get((d, n), 1.0) for n, e in effs.items()) / total_eff
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
            out["reference B(E2)"] = ref_unc
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
            try:
                first = fit_peak(counts, edges, self.e0, guess, 3 * s.fit_half_width, max_sigma=0.03 * self.e0)
                areas.append(fit_peak(counts, edges, self.e0, first.sigma, s.fit_half_width,
                                      max_sigma=0.03 * self.e0).area)
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


__all__ = ["Analysis", "PeakFit", "Result", "SHIFT_H", "Settings", "ex_nuclide", "fit_peak"]
