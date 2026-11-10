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


#: The most events one simulation of the reference of :func:`fit_offset` generates.
REFERENCE_EVENTS = 400_000

#: The most such simulations the reference is made of.
REFERENCE_CHUNKS = 4


def _corrected(gammas, experiment_assumed) -> np.ndarray:
    exc = experiment_assumed.excitation
    key = "corrected_recoil" if exc.excite == "target" else "corrected_projectile"
    return recorrect(gammas, experiment_assumed)[key]


def _centroid(x: np.ndarray, w: np.ndarray, e0: float, sigma: float) -> tuple:
    """Weighted mean and its error of the peak's energies within ±2.5σ of its centre, found by two passes (a
    window about E₀, then about the mean), with the continuum under the peak taken as the level of the bands
    beside it and subtracted from the moments, as one does with real data. ``sigma`` is the width of the peak
    as it stands in the corrected spectrum (:func:`peak_widths`), not the detector's resolution alone."""
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


def _peak_sigma(x: np.ndarray, w: np.ndarray, e0: float, sigma_det: float) -> float:
    """σ of the peak at ``e0`` as it stands in a spectrum: the width at half maximum of the smoothed histogram
    (bins of the detector's σ, five-bin running mean) over its median level, the crossings interpolated. NaN
    with fewer than 30 γ rays above the half maximum."""
    half = 0.06 * e0
    edges = np.arange(e0 - half, e0 + half + sigma_det, sigma_det)
    if len(x) < 30 or len(edges) < 12:
        return float("nan")
    h, _ = np.histogram(x, bins=edges, weights=w)
    smooth = np.convolve(h, np.ones(5) / 5, mode="same") - float(np.median(h))
    k = int(np.argmax(smooth))
    top = float(smooth[k])
    if top <= 0:
        return float("nan")
    lo = hi = k
    while lo > 0 and smooth[lo - 1] > top / 2:
        lo -= 1
    while hi < len(smooth) - 1 and smooth[hi + 1] > top / 2:
        hi += 1
    if int(np.sum((x >= edges[lo]) & (x < edges[hi + 1]))) < 30:
        return float("nan")
    left = lo - (smooth[lo] - top / 2) / (smooth[lo] - smooth[lo - 1]) if lo > 0 else float(lo)
    right = hi + (smooth[hi] - top / 2) / (smooth[hi] - smooth[hi + 1]) if hi < len(smooth) - 1 else float(hi)
    return float(right - left) * sigma_det / FWHM_PER_SIGMA


def _table(gammas, experiment_assumed) -> dict:
    """What the diagnostic needs of the counted γ rays, corrected with the assumed geometry: the corrected
    energy, the counts each stands for in the run, the crystal, and the particle's detector and ring."""
    x = _corrected(gammas, experiment_assumed)
    c = gammas.events.columns
    keep = gammas["counted"]
    particle = gammas["particle"][keep]
    t = gammas.events.beam_time_s * gammas.live_fraction
    return {"x": x[keep], "w": (gammas["weight"] * t)[keep], "crystal": gammas["crystal"][keep],
            "detector": c["detector"][particle], "ring": c["segment_i"][particle]}


def _sigma_det(crystal, e0: float) -> float:
    return max(float(crystal.response.fwhm(e0)) / FWHM_PER_SIGMA, 1e-4)


def _widths(table: dict, gammas) -> dict:
    e0 = gammas.energy_mev
    out = {}
    for k, cr in enumerate(gammas.crystals):
        m = table["crystal"] == k
        sigma = _sigma_det(cr, e0)
        seen = _peak_sigma(table["x"][m], table["w"][m], e0, sigma)
        out[cr.name] = max(3 * sigma, seen) if math.isfinite(seen) else 3 * sigma
    return out


def peak_widths(gammas, experiment_assumed) -> dict:
    """σ (MeV) of the corrected peak in each crystal, {crystal name: σ}, measured on the spectrum of all its γ
    rays corrected with the assumed geometry. A crystal sees the emitter over its whole opening angle but the
    correction uses its centre, so a fast emitter leaves a peak many times wider than the detector's resolution
    (62Ni at 9% of c in a crystal of 35 mm radius at 120 mm: 20 keV FWHM against 2.5 keV). Never less than three
    times the detector's σ, which is also the width taken when the crystal has too few γ rays to measure it."""
    return _widths(_table(gammas, experiment_assumed), gammas)


def _rows(table: dict, gammas, widths: dict, least: int) -> list:
    e0 = gammas.energy_mev
    out = []
    for k, cr in enumerate(gammas.crystals):
        sigma = widths.get(cr.name, 3 * _sigma_det(cr, e0))
        in_crystal = table["crystal"] == k
        for d, name in enumerate(gammas.events.detectors):
            m_det = in_crystal & (table["detector"] == d)
            for ring in np.unique(table["ring"][m_det]):
                m = m_det & (table["ring"] == ring)
                if m.sum() < least:
                    continue
                x, w = table["x"][m], table["w"][m]
                centre, err, n = _centroid(x, w, e0, sigma)
                if not math.isfinite(centre) or n < least:
                    continue
                out.append({"crystal": cr.name, "detector": name, "ring": int(ring), "centroid_kev": 1e3 * centre,
                            "error_kev": 1e3 * err, "n": n, "counts": float(w.sum())})
    return out


def diagnostic(gammas, experiment_assumed, least: int = 10, widths: Optional[dict] = None) -> list:
    """The corrected peak's centroid per crystal and ring (or strip): rows {"crystal", "detector", "ring",
    "centroid_kev", "error_kev", "n", "counts"} with at least ``least`` γ rays in the peak's window, over the
    run's statistics. The window is ±2.5σ with σ the peak's width in that crystal (:func:`peak_widths`, measured
    here unless ``widths`` gives them: a scan over geometries keeps one set, so that every step and the
    reference take the same window)."""
    table = _table(gammas, experiment_assumed)
    return _rows(table, gammas, widths if widths is not None else _widths(table, gammas), least)


def _flatness(rows: list, e0_kev: Optional[float] = None, reference: Optional[list] = None) -> tuple:
    """(χ², points compared, means taken out) of :func:`flatness`."""
    total, points, means = 0.0, 0, 0
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
            points, means = points + len(pairs), means + 1
            continue
        if len(mine) < 2 and e0_kev is None:
            continue
        if e0_kev is None:
            wts = np.array([1 / r["error_kev"] ** 2 for r in mine])
            ref = float(np.sum(wts * np.array([r["centroid_kev"] for r in mine])) / wts.sum())
            means += 1
        else:
            ref = e0_kev
        total += sum(((r["centroid_kev"] - ref) / r["error_kev"]) ** 2 for r in mine)
        points += len(mine)
    return float(total), points, means


def flatness(rows: list, e0_kev: Optional[float] = None, reference: Optional[list] = None) -> float:
    """χ² of the centroids about each crystal's own mean over the rings: flat lines make it small whatever the
    energy calibration, a misplaced target slopes them. With ``e0_kev`` the lines are measured against the
    transition energy instead.

    With ``reference`` (the diagnostic of a simulation made with the geometry assumed, corrected with it), the
    centroids are measured against that simulation's own pattern, ring by ring: large crystals and the energy
    loss leave the corrected centroids a little ring-dependent even with the right geometry, and the simulation
    knows that pattern."""
    return _flatness(rows, e0_kev, reference)[0]


def with_offset(experiment, offset_mm: float) -> Experiment:
    """The experiment with the target moved by ``offset_mm`` along the beam from where it is."""
    d = experiment.to_dict()
    d["target"]["position"] = f"{experiment.target.position_mm + offset_mm:.6g} mm"
    if abs(experiment.target.position_mm + offset_mm) < 1e-9:
        d["target"].pop("position", None)
    return Experiment.from_dict(d)


def _reference(experiment, gammas, wanted: int, widths: dict, seed: int) -> tuple:
    """The diagnostic of a simulation made with ``experiment`` and corrected with it, holding about ``wanted``
    counted γ rays: (rows, γ rays, events simulated).

    A run counts every reaction once, tens of millions in an hour of beam, and a simulation of that many
    weighted events does not fit in memory (over 14 GB for 63 million). What sets the reference's precision is
    its γ rays, a few per hundred events, so it is made of simulations of at most :data:`REFERENCE_EVENTS`
    events each, one after the other (only the γ rays of each are kept), until it has as many γ rays as the
    data, and of at most :data:`REFERENCE_CHUNKS` of them. A first simulation with hardly any γ rays (see
    :func:`fit_offset`) is left alone: more of the same would give no pattern either."""
    from .gamma_events import simulate_gammas

    events = max(min(int(gammas.events.n_events), REFERENCE_EVENTS), 1)
    tables, chunks, first = [], 1, None
    k = 0
    while k < chunks:
        sim = simulate_gammas(experiment, events, seed=seed + k)
        tables.append(_table(sim, experiment))
        if first is None:
            first = sim
            got = max(len(tables[0]["x"]), 1)
            need = wanted / got
            chunks = 1 if need > 4 * REFERENCE_CHUNKS else min(max(int(round(need)), 1), REFERENCE_CHUNKS)
        del sim
        k += 1
    table = {key: np.concatenate([t[key] for t in tables]) for key in tables[0]}
    table["w"] = table["w"] / len(tables)
    return _rows(table, first, widths, 10), len(table["x"]), events * len(tables)


def _covers(reference: list, rows: list) -> bool:
    """Whether the reference has a centroid for at least half of the points of ``rows``."""
    have = {(r["crystal"], r["detector"], r["ring"]) for r in reference}
    return 2 * sum((r["crystal"], r["detector"], r["ring"]) in have for r in rows) >= max(len(rows), 1)


def _minimum(grid: np.ndarray, chi2: np.ndarray) -> tuple:
    """The minimum of a scan by a parabola through its five lowest points: (offset, uncertainty for Δχ² = 1,
    χ² there, index of the lowest point)."""
    k = int(np.argmin(chi2))
    lo, hi = max(k - 2, 0), min(k + 2, len(grid) - 1)
    coeff = np.polyfit(grid[lo:hi + 1], chi2[lo:hi + 1], 2) if hi - lo >= 2 else (0.0, 0.0, 0.0)
    z_min = -coeff[1] / (2 * coeff[0]) if coeff[0] > 0 else float("nan")
    if grid[lo] <= z_min <= grid[hi]:  # a parabola whose minimum lies outside its own points says nothing
        return float(z_min), 1 / math.sqrt(coeff[0]), float(np.polyval(coeff, z_min)), k
    return float(grid[k]), float(grid[1] - grid[0]) if len(grid) > 1 else float("nan"), float(chi2[k]), k


def fit_offset(gammas, experiment, span_mm: float = 8.0, steps: int = 17, reference: bool = True,
               seed: int = 2, refine: int = 2) -> dict:
    """The target offset along the beam (from the position ``experiment`` assumes) that flattens the diagnostic
    plot: χ² of the centroids is scanned over ±``span_mm`` and a parabola through its lowest points gives the
    minimum and its uncertainty (Δχ² = 1, scaled by √(χ²/ν) when the centroids scatter by more than their
    errors). Returns {"offset_mm", "uncertainty_mm", "chi2", "degrees_of_freedom", "at_edge", "scan", "against",
    "passes", "reference_at_mm", "reference_gammas", "reference_events"}.

    With ``reference`` (the default) the pattern is judged against that of a simulation made with the assumed
    geometry and corrected with it, which holds the ring dependence that the finite crystals and the energy loss
    leave even when the geometry is right; without it, against flat lines. The reference holds as many γ rays as
    ``gammas`` where that fits (see :func:`_reference`): its errors are in χ² and in the uncertainty.

    That pattern belongs to the geometry the reference was made with. When the fit finds the target more than
    two uncertainties from where it was assumed, the reference is made again with the target where the fit put
    it and the scan judged against that, up to ``refine`` times and until the offset moves by less than its
    uncertainty: with the reference left at an assumed geometry 5 mm off, the fit overshoots by about 1 mm.
    ``"reference_at_mm"`` is the offset the reference in use was made with, and ``"passes"`` how many were made.

    A simulation may give too few γ rays for a pattern: with a detector's edge moved to where it sees the slow
    recoils of grazing collisions, nearly every generated event is one of those. Such a reference is not used:
    ``"against"`` is then ``"flat"`` (flat lines) if it was the first, and a later one ends the refinement."""
    widths = peak_widths(gammas, experiment)
    grid = np.linspace(-span_mm, span_mm, steps)
    scanned = [diagnostic(gammas, with_offset(experiment, z), widths=widths) for z in grid]
    wanted = int(np.sum(gammas["counted"]))

    def solve(ref):
        parts = [_flatness(rows, reference=ref) for rows in scanned]
        chi2 = np.array([p[0] for p in parts])
        z_min, unc, best, k = _minimum(grid, chi2)
        dof = max(parts[k][1] - parts[k][2] - 1, 1)
        if best > dof:
            unc *= math.sqrt(best / dof)
        return {"offset_mm": z_min, "uncertainty_mm": float(unc), "chi2": best, "degrees_of_freedom": int(dof),
                "at_edge": k in (0, len(grid) - 1),
                "scan": [{"offset_mm": float(z), "chi2": float(c)} for z, c in zip(grid, chi2)]}

    used = {"against": "flat", "passes": 0, "reference_at_mm": None, "reference_gammas": 0, "reference_events": 0}
    if not reference:
        return {**solve(None), **used}
    middle = scanned[len(scanned) // 2]
    ref, n_gammas, n_events = _reference(experiment, gammas, wanted, widths, seed)
    passes = 1
    if _covers(ref, middle):
        used = {"against": "simulation", "reference_at_mm": 0.0, "reference_gammas": n_gammas,
                "reference_events": n_events}
        out = solve(ref)
    else:
        out = solve(None)
    while passes <= refine and not out["at_edge"] and abs(out["offset_mm"]) > 2 * out["uncertainty_mm"]:
        there = with_offset(experiment, out["offset_mm"])
        ref, n_gammas, n_events = _reference(there, gammas, wanted, widths, seed)
        passes += 1
        if not _covers(ref, middle):
            break
        used = {"against": "simulation", "reference_at_mm": out["offset_mm"], "reference_gammas": n_gammas,
                "reference_events": n_events}
        before, out = out["offset_mm"], solve(ref)
        if abs(out["offset_mm"] - before) < out["uncertainty_mm"]:
            break
    return {**out, **used, "passes": passes}


def overlay(gammas, experiment_true, experiment_assumed, bins: int = 200) -> dict:
    """The corrected peak with the true geometry and with the assumed one, as spectra (counts in the run per
    bin) with the centroid shift and the change in width: {"edges", "true", "assumed", "centroid_true_kev",
    "centroid_assumed_kev", "fwhm_true_kev", "fwhm_assumed_kev"}. The spectra span ±8σ of the peak as it stands
    with the true geometry (all crystals together, :func:`_peak_sigma`) and each width is the peak's second
    moment within ±3σ of its centroid."""
    e0 = gammas.energy_mev
    t = gammas.events.beam_time_s * gammas.live_fraction
    w = gammas["weight"] * t
    m = gammas["counted"]
    sigma_det = max(float(np.mean([c.response.fwhm(e0) for c in gammas.crystals])) / FWHM_PER_SIGMA, 1e-4)
    seen = _peak_sigma(_corrected(gammas, experiment_true)[m], w[m], e0, sigma_det)
    sigma = max(3 * sigma_det, seen) if math.isfinite(seen) else 3 * sigma_det
    edges = np.linspace(e0 - 8 * sigma, e0 + 8 * sigma, bins + 1)
    out = {"edges": edges}
    for label, exp in (("true", experiment_true), ("assumed", experiment_assumed)):
        x = _corrected(gammas, exp)
        h, _ = np.histogram(x[m], bins=edges, weights=w[m])
        centre, err, _ = _centroid(x[m], w[m], e0, sigma)
        inside = np.abs(x[m] - centre) < 3 * sigma
        sd = math.sqrt(float(np.average((x[m][inside] - centre) ** 2, weights=w[m][inside]))) if inside.any() else 0.0
        out[label] = h
        out[f"centroid_{label}_kev"] = 1e3 * centre
        out[f"fwhm_{label}_kev"] = 1e3 * FWHM_PER_SIGMA * sd
    out["shift_kev"] = out["centroid_assumed_kev"] - out["centroid_true_kev"]
    out["broadening_kev"] = out["fwhm_assumed_kev"] - out["fwhm_true_kev"]
    return out


__all__ = ["REFERENCE_CHUNKS", "REFERENCE_EVENTS", "diagnostic", "fit_offset", "flatness", "overlay", "peak_widths",
           "with_offset"]
