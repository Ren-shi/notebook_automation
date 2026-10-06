"""Finding a misplaced target from the data, as one does with real data: the Doppler-corrected peak's centroid
against ring number, crystal by crystal, is flat when the geometry assumed in the correction is right and sloped
when it is not.

::

    from physim.nuclear.alignment import diagnostic, fit_offset, overlay

    rows = diagnostic(gammas, experiment_assumed)      # centroid ± error per (crystal, ring)
    fit = fit_offset(gammas, experiment)                # the target offset along the beam that flattens it
    both = overlay(gammas, experiment, experiment_assumed)   # the corrected peak with each geometry

``gammas`` are the simulated γ rays (:func:`~physim.nuclear.gamma_events.simulate_gammas`), made with the true
geometry; ``experiment_assumed`` differs from it in the target's ``position``. Everything here uses only what an
experimentalist has: the recorded energies, the segment and crystal that fired, and the geometry assumed.
"""

from __future__ import annotations

import copy
import math
from typing import Optional

import numpy as np

from .experiment import Experiment
from .gamma_events import recorrect
from .response import FWHM_PER_SIGMA


def _corrected(gammas, experiment_assumed) -> np.ndarray:
    exc = experiment_assumed.excitation
    key = "corrected_recoil" if exc.excite == "target" else "corrected_projectile"
    return recorrect(gammas, experiment_assumed)[key]


def _centroid(x: np.ndarray, w: np.ndarray, e0: float, sigma: float) -> tuple:
    """Weighted mean and its error of the peak's energies within ±2.5σ of its centre, found by two passes (a
    window about E₀, then about the mean), with the continuum under the peak taken as the level of the bands
    beside it and subtracted from the moments, as one does with real data."""
    centre = e0
    m = np.zeros(len(x), dtype=bool)
    for _ in range(3):
        m = np.abs(x - centre) < 2.5 * sigma
        side = (np.abs(x - centre) >= 2.5 * sigma) & (np.abs(x - centre) < 6 * sigma)
        if m.sum() < 3 or w[m].sum() <= 0:
            return float("nan"), float("nan"), 0
        # Background per unit energy from the bands, taken off the peak window's moments.
        level = float(w[side].sum()) / (7 * sigma) if side.any() else 0.0
        under = level * 5 * sigma
        total = float(w[m].sum()) - under
        if total <= 0:
            return float("nan"), float("nan"), 0
        centre = (float(np.sum(w[m] * x[m])) - under * centre) / total
    sd = math.sqrt(max(float(np.average((x[m] - centre) ** 2, weights=w[m])), 1e-12))
    err = sd * math.sqrt(float(np.sum(w[m] ** 2))) / max(float(np.sum(w[m])) - under, 1e-12)
    return centre, err, int(m.sum())


def diagnostic(gammas, experiment_assumed, least: int = 10) -> list:
    """The corrected peak's centroid per crystal and ring (or strip): rows {"crystal", "detector", "ring",
    "centroid_kev", "error_kev", "n"} with at least ``least`` γ rays, over the run's statistics."""
    x = _corrected(gammas, experiment_assumed)
    e0 = gammas.energy_mev
    c = gammas.events.columns
    rows_p = gammas["particle"]
    t = gammas.events.beam_time_s * gammas.live_fraction
    w = gammas["weight"] * t
    counted = gammas["counted"]
    out = []
    for k, cr in enumerate(gammas.crystals):
        sigma = max(float(cr.response.fwhm(e0)) / FWHM_PER_SIGMA, 1e-4)
        in_crystal = (gammas["crystal"] == k) & counted
        for d, name in enumerate(gammas.events.detectors):
            m_det = in_crystal & (c["detector"][rows_p] == d)
            for ring in np.unique(c["segment_i"][rows_p][m_det]):
                m = m_det & (c["segment_i"][rows_p] == ring)
                if m.sum() < least:
                    continue
                centre, err, n = _centroid(x[m], w[m], e0, 3 * sigma)
                if not math.isfinite(centre):
                    continue
                # Counting statistics of the run, not of the sample: σ/√N with N the counts in the run.
                counts = float(w[m].sum())
                sd_run = err * math.sqrt(float(np.sum(w[m] ** 2))) / float(np.sum(w[m])) if counts > 0 else err
                out.append({"crystal": cr.name, "detector": name, "ring": int(ring), "centroid_kev": 1e3 * centre,
                            "error_kev": 1e3 * max(err, sd_run), "n": n, "counts": counts})
    return out


def flatness(rows: list, e0_kev: Optional[float] = None, reference: Optional[list] = None) -> float:
    """χ² of the centroids about each crystal's own mean over the rings: flat lines make it small whatever the
    energy calibration, a misplaced target slopes them. With ``e0_kev`` the lines are measured against the
    transition energy instead.

    With ``reference`` (the diagnostic of a simulation made with the geometry assumed, corrected with it), the
    centroids are measured against that simulation's own pattern, ring by ring: large crystals and the energy
    loss leave the corrected centroids a little ring-dependent even with the right geometry, and the simulation
    knows that pattern."""
    total = 0.0
    ref_by = {(r["crystal"], r["detector"], r["ring"]): r for r in reference} if reference else {}
    for crystal in {r["crystal"] for r in rows}:
        mine = [r for r in rows if r["crystal"] == crystal and r["error_kev"] > 0]
        if reference:
            pairs = [(r, ref_by[(r["crystal"], r["detector"], r["ring"])]) for r in mine
                     if (r["crystal"], r["detector"], r["ring"]) in ref_by]
            if len(pairs) < 2:
                continue
            diff = np.array([r["centroid_kev"] - q["centroid_kev"] for r, q in pairs])
            err = np.array([math.hypot(r["error_kev"], q["error_kev"]) for r, q in pairs])
            wts = 1 / err**2
            common = float(np.sum(wts * diff) / wts.sum())
            total += float(np.sum(((diff - common) / err) ** 2))
            continue
        if len(mine) < 2 and e0_kev is None:
            continue
        if e0_kev is None:
            wts = np.array([1 / r["error_kev"] ** 2 for r in mine])
            ref = float(np.sum(wts * np.array([r["centroid_kev"] for r in mine])) / wts.sum())
        else:
            ref = e0_kev
        total += sum(((r["centroid_kev"] - ref) / r["error_kev"]) ** 2 for r in mine)
    return float(total)


def with_offset(experiment, offset_mm: float) -> Experiment:
    """The experiment with the target moved by ``offset_mm`` along the beam from where it is."""
    d = experiment.to_dict()
    d["target"]["position"] = f"{experiment.target.position_mm + offset_mm:.6g} mm"
    if abs(experiment.target.position_mm + offset_mm) < 1e-9:
        d["target"].pop("position", None)
    return Experiment.from_dict(d)


def fit_offset(gammas, experiment, span_mm: float = 8.0, steps: int = 17, reference: bool = True,
               seed: int = 2) -> dict:
    """The target offset along the beam (from the position ``experiment`` assumes) that flattens the diagnostic
    plot: χ² of the centroids is scanned over ±``span_mm`` and a parabola through its lowest points gives the
    minimum and its uncertainty (Δχ² = 1). Returns {"offset_mm", "uncertainty_mm", "chi2", "scan"}.

    With ``reference`` (the default) the pattern is judged against that of a simulation made with the assumed
    geometry and corrected with it (the same number of events, another seed), which holds the ring dependence
    that the finite crystals and the energy loss leave even when the geometry is right; without it, against
    flat lines."""
    from .gamma_events import simulate_gammas

    ref = None
    if reference:
        sim = simulate_gammas(experiment, gammas.events.n_events, seed=seed)
        ref = diagnostic(sim, experiment)
    grid = np.linspace(-span_mm, span_mm, steps)
    chi2 = np.array([flatness(diagnostic(gammas, with_offset(experiment, z)), reference=ref) for z in grid])
    k = int(np.argmin(chi2))
    lo, hi = max(k - 2, 0), min(k + 2, len(grid) - 1)
    coeff = np.polyfit(grid[lo:hi + 1], chi2[lo:hi + 1], 2)
    if coeff[0] > 0:
        z_min = -coeff[1] / (2 * coeff[0])
        unc = 1 / math.sqrt(coeff[0])
        best = float(np.polyval(coeff, z_min))
    else:
        z_min, unc, best = float(grid[k]), float(grid[1] - grid[0]), float(chi2[k])
    at_edge = k in (0, len(grid) - 1)
    return {"offset_mm": float(z_min), "uncertainty_mm": float(unc), "chi2": best, "at_edge": at_edge,
            "scan": [{"offset_mm": float(z), "chi2": float(c)} for z, c in zip(grid, chi2)],
            "degrees_of_freedom": len(diagnostic(gammas, experiment)), "against": "simulation" if reference else "flat"}


def overlay(gammas, experiment_true, experiment_assumed, bins: int = 200) -> dict:
    """The corrected peak with the true geometry and with the assumed one, as spectra (counts in the run per
    bin) with the centroid shift and the change in width: {"edges", "true", "assumed", "centroid_true_kev",
    "centroid_assumed_kev", "fwhm_true_kev", "fwhm_assumed_kev"}."""
    e0 = gammas.energy_mev
    t = gammas.events.beam_time_s * gammas.live_fraction
    w = gammas["weight"] * t
    m = gammas["counted"]
    sigma = max(float(np.mean([c.response.fwhm(e0) for c in gammas.crystals])) / FWHM_PER_SIGMA, 1e-4)
    edges = np.linspace(e0 - 20 * sigma, e0 + 20 * sigma, bins + 1)
    out = {"edges": edges}
    for label, exp in (("true", experiment_true), ("assumed", experiment_assumed)):
        x = _corrected(gammas, exp)
        h, _ = np.histogram(x[m], bins=edges, weights=w[m])
        centre, err, _ = _centroid(x[m], w[m], e0, 3 * sigma)
        inside = np.abs(x[m] - centre) < 3 * sigma
        sd = math.sqrt(float(np.average((x[m][inside] - centre) ** 2, weights=w[m][inside]))) if inside.any() else 0.0
        out[label] = h
        out[f"centroid_{label}_kev"] = 1e3 * centre
        out[f"fwhm_{label}_kev"] = 1e3 * FWHM_PER_SIGMA * sd
    out["shift_kev"] = out["centroid_assumed_kev"] - out["centroid_true_kev"]
    out["broadening_kev"] = out["fwhm_assumed_kev"] - out["fwhm_true_kev"]
    return out


__all__ = ["diagnostic", "fit_offset", "flatness", "overlay", "with_offset"]
