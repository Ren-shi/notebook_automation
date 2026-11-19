"""Lifetime analysis for plunger (RDDS) and Doppler-shift attenuation (DSAM) measurements: from spectra to τ.

Everything here reads spectra (arrays of counts and bin edges) or numbers taken from them (peak areas against
distance), never a setup or a run, so it works the same on simulated and on measured data (``planning/lifetimes.md``,
wave 2C; the signatures are the frozen interfaces of wave 0)::

    from physim.nuclear import lifetime_analysis as la

    areas = la.peak_areas(spectra, 652.0, 0.045, 36.5, 2.5)            # per distance: (shifted, unshifted, errors)
    s, u, e = (np.array(x) for x in zip(*areas))
    la.decay_curve(distances, u / (u + s), ..., 0.045)                  # fit of the unshifted share (Bateman)
    la.ddcm(distances, u, s, (e[:, 1], e[:, 0]), 0.045)                 # differential decay-curve method
    la.fit_lineshape(counts, edges, 1193.1, kernel, 52.8, fwhm_kev=2.8) # DSAM

Units: time ps, distance µm, γ-ray energy keV, speed as β. τ is the proper (rest-frame) mean life: plunger flight
times are divided by γ = 1/√(1 − β²) unless ``proper=False`` (DSAM kernels are in lab time, and time dilation, below
0.2% at β < 6%, is left out there, as in :mod:`~physim.nuclear.lifetimes`).

**Results** are :class:`TauResult`: τ, its statistical uncertainty and a ``budget`` that always names
``"statistics"``, ``"distance"`` (the offset of the zero distance), ``"velocity"``, ``"feeding"`` and ``"stopping"``
(each the change of τ when that input moves by its uncertainty; 0 where it does not enter, with a note). The
stopping entry of a DSAM fit scales the kernel's time axis by ± ``stopping_uncertainty`` (default 7%: stopping powers
of heavy ions at these speeds are known to about 5–10%, and a uniform scaling of the stopping is a uniform scaling of
the slowing-down time).

**Feeding** of a decay curve or a line shape is a :class:`Feeding`: the cascade in :mod:`~physim.nuclear.lifetimes`'
convention (every level's τ, the feeding matrix row = from, column = to, the populations), the level being measured,
and which other τ are fitted too (a side-feeding time, say). The DDCM takes the feeder's intensities instead.

**No scipy:** the fits are small Levenberg–Marquardt and golden-section searches written here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

from . import lifetimes as lt
from .analysis import FWHM_PER_SIGMA, fit_peak

#: The budget's entries, always present.
BUDGET_KEYS = ("statistics", "distance", "velocity", "feeding", "stopping")
#: The default relative uncertainty of the stopping powers for a DSAM fit (the kernel's time axis is scaled by ± it).
STOPPING_UNCERTAINTY = 0.07


# ---------------------------------------------------------------------------------------------------------------
# Results and feeding


@dataclass
class TauResult:
    """A measured mean life.

    ``tau_ps`` and its statistical uncertainty ``stat_ps``; ``budget`` {source: ps} with at least
    :data:`BUDGET_KEYS` (``budget["statistics"]`` is ``stat_ps``); ``method`` ("decay curve", "DDCM", "DDCM gated",
    "DDCM three-foil", "line shape"); ``per_distance`` (a list of dicts, one per distance, for plots: the τ(x) of the
    DDCM, the data and the model of a decay curve) or None; ``notes``, plain sentences; ``params``, other fitted
    numbers (χ², side-feeding time, asymmetric errors, ...)."""

    tau_ps: float
    stat_ps: float
    budget: dict
    method: str
    per_distance: Optional[list] = None
    notes: list = field(default_factory=list)
    params: dict = field(default_factory=dict)

    @property
    def total_ps(self) -> float:
        """All the budget's entries added in quadrature."""
        return float(math.sqrt(sum(float(v) ** 2 for v in self.budget.values())))

    def as_dict(self) -> dict:
        return {"tau_ps": self.tau_ps, "stat_ps": self.stat_ps, "total_ps": self.total_ps,
                "budget": dict(self.budget), "method": self.method, "per_distance": self.per_distance,
                "notes": list(self.notes), "params": dict(self.params)}


def _budget(stat: float, **others) -> dict:
    out = {k: 0.0 for k in BUDGET_KEYS}
    out["statistics"] = float(stat)
    for k, v in others.items():
        out[k] = float(abs(v))
    return out


@dataclass
class Feeding:
    """The cascade that feeds the measured level, in :mod:`~physim.nuclear.lifetimes`' convention.

    ``taus_ps``: every level's mean life (ps); the measured level's entry, and those in ``free``, are starting
    guesses. ``feeding``: (levels × levels) share of the row level's decays that go to the column level, or None.
    ``populations``: the nuclei put directly into each level. ``level``: the index of the measured level.
    ``free``: indices of other levels whose τ is fitted too (a side-feeding time). ``tau_errors_ps``: one
    uncertainty per level (0 or None where exact or fitted), for the budget's "feeding" entry."""

    taus_ps: Sequence[float]
    feeding: Optional[np.ndarray]
    populations: Sequence[float]
    level: int
    free: tuple = ()
    tau_errors_ps: Optional[Sequence[Optional[float]]] = None

    @classmethod
    def side(cls, side_ps: float, side_share: float, *, free: bool = False, side_error_ps: float = 0.0,
             tau_guess_ps: float = 1.0) -> "Feeding":
        """The usual DSAM case: the measured level (index 1), a share ``side_share`` of whose population comes
        through a side-feeding cascade with the effective mean life ``side_ps`` (index 0), fitted when ``free``."""
        if not 0 <= side_share <= 1:
            raise ValueError("the side-feeding share must be in [0, 1]")
        return cls([float(side_ps), float(tau_guess_ps)], np.array([[0.0, 1.0], [0.0, 0.0]]),
                   [side_share, 1.0 - side_share], level=1, free=(0,) if free else (),
                   tau_errors_ps=[float(side_error_ps), 0.0])

    def checked(self) -> "Feeding":
        n = len(self.taus_ps)
        if not 0 <= int(self.level) < n:
            raise ValueError("the measured level must be one of the cascade's levels")
        if int(self.level) in tuple(self.free):
            raise ValueError("the measured level's τ is always fitted: leave it out of 'free'")
        if any(not 0 <= int(i) < n for i in self.free):
            raise ValueError("'free' lists indices of the cascade's levels")
        lt._cascade(self.taus_ps, self.feeding, self.populations)  # raises on a bad cascade
        return self

    def errors(self) -> list:
        e = self.tau_errors_ps
        if e is None:
            return [0.0] * len(self.taus_ps)
        if len(e) != len(self.taus_ps):
            raise ValueError("give one τ uncertainty per level of the cascade")
        return [0.0 if x is None else float(x) for x in e]


def _one_level() -> Feeding:
    return Feeding([1.0], None, [1.0], 0)


# ---------------------------------------------------------------------------------------------------------------
# Small numerical tools (no scipy)


def _jacobian(fun: Callable, p: np.ndarray, r0: np.ndarray, lower, upper) -> np.ndarray:
    jac = np.empty((len(r0), len(p)))
    for j in range(len(p)):
        h = 1e-6 * max(abs(p[j]), 1e-3)
        hi, lo = p.copy(), p.copy()
        hi[j] = min(p[j] + h, upper[j])
        lo[j] = max(p[j] - h, lower[j])
        if hi[j] == lo[j]:
            jac[:, j] = 0.0
            continue
        jac[:, j] = (fun(hi) - fun(lo)) / (hi[j] - lo[j])
    return jac


def _least_squares(fun: Callable, p0, lower=None, upper=None, max_iter: int = 200) -> tuple:
    """Levenberg–Marquardt on weighted residuals ``fun(p)``, with box bounds by clipping. Returns (p, covariance
    (JᵀJ)⁻¹ at the solution, χ² = Σ r²)."""
    p = np.array(p0, dtype=float)
    n = len(p)
    lower = np.full(n, -np.inf) if lower is None else np.asarray(lower, dtype=float)
    upper = np.full(n, np.inf) if upper is None else np.asarray(upper, dtype=float)
    p = np.clip(p, lower, upper)
    r = fun(p)
    chi2 = float(r @ r)
    lam = 1e-3
    for _ in range(max_iter):
        jac = _jacobian(fun, p, r, lower, upper)
        a = jac.T @ jac
        g = jac.T @ r
        improved = False
        while lam < 1e12:
            try:
                step = np.linalg.solve(a + lam * np.diag(np.diag(a) + 1e-30), -g)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            q = np.clip(p + step, lower, upper)
            rq = fun(q)
            cq = float(rq @ rq)
            if np.isfinite(cq) and cq <= chi2:
                improved = True
                done = chi2 - cq <= 1e-12 * max(chi2, 1e-300) and np.all(np.abs(q - p) <= 1e-9 * (np.abs(p) + 1e-9))
                p, r, chi2 = q, rq, cq
                lam = max(lam / 10, 1e-12)
                break
            lam *= 10
        if not improved or done:
            break
    jac = _jacobian(fun, p, r, lower, upper)
    cov = np.linalg.pinv(jac.T @ jac)
    return p, cov, chi2


def _golden(f: Callable, a: float, b: float, tol: float = 1e-6) -> tuple:
    """Minimum of a unimodal f on [a, b]: (x, f(x))."""
    g = (math.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while abs(b - a) > tol * (abs(a) + abs(b) + 1e-12):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (c, fc) if fc < fd else (d, fd)


def _minimise_1d(f: Callable, lo: float, hi: float, n_scan: int = 25, tol: float = 1e-6) -> tuple:
    """Scan f on n_scan points of [lo, hi], then golden-section between the best point's neighbours."""
    xs = np.linspace(lo, hi, n_scan)
    fs = np.array([f(x) for x in xs])
    i = int(np.nanargmin(fs))
    a, b = xs[max(i - 1, 0)], xs[min(i + 1, n_scan - 1)]
    x, fx = _golden(f, a, b, tol)
    if fs[i] < fx:
        x, fx = xs[i], fs[i]
    return float(x), float(fx), i in (0, n_scan - 1)


def _crossing(f: Callable, x0: float, f_target: float, direction: int, limit: float) -> Optional[float]:
    """The x beyond x0 (towards ``direction``) where f rises through f_target, by stepping out then bisection; None
    when f stays below it up to ``limit``."""
    step = 0.05
    a, x = x0, x0
    while True:
        x = x0 + direction * step
        if (direction > 0 and x >= limit) or (direction < 0 and x <= limit):
            x = limit
            if f(x) < f_target:
                return None
            break
        if f(x) >= f_target:
            break
        a = x
        step *= 1.6
    b = x
    for _ in range(60):
        m = 0.5 * (a + b)
        if f(m) < f_target:
            a = m
        else:
            b = m
        if abs(b - a) < 1e-7:
            break
    return 0.5 * (a + b)


def _nelder_mead(f: Callable, x0, step: float, tol: float = 1e-7, max_iter: int = 400) -> tuple:
    """Minimum of f by the Nelder–Mead simplex: (x, f(x))."""
    n = len(x0)
    pts = [np.array(x0, dtype=float)]
    for i in range(n):
        p = pts[0].copy()
        p[i] += step
        pts.append(p)
    vals = [f(p) for p in pts]
    for _ in range(max_iter):
        order = np.argsort(vals)
        pts, vals = [pts[i] for i in order], [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) <= tol * (abs(vals[0]) + 1e-12) and \
                max(np.max(np.abs(p - pts[0])) for p in pts[1:]) < 1e-5:
            break
        centre = np.mean(pts[:-1], axis=0)
        xr = centre + (centre - pts[-1])
        fr = f(xr)
        if fr < vals[0]:
            xe = centre + 2 * (centre - pts[-1])
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = centre + 0.5 * (pts[-1] - centre)
            fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = pts[0] + 0.5 * (pts[i] - pts[0])
                    vals[i] = f(pts[i])
    i = int(np.argmin(vals))
    return pts[i], float(vals[i])


def _hessian(f: Callable, x, h: float) -> np.ndarray:
    """Central-difference Hessian of f at x."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    f0 = f(x)
    hess = np.empty((n, n))
    for i in range(n):
        e_i = np.zeros(n)
        e_i[i] = h
        hess[i, i] = (f(x + e_i) - 2 * f0 + f(x - e_i)) / h ** 2
        for j in range(i):
            e_j = np.zeros(n)
            e_j[j] = h
            hess[i, j] = hess[j, i] = (f(x + e_i + e_j) - f(x + e_i - e_j) - f(x - e_i + e_j)
                                       + f(x - e_i - e_j)) / (4 * h * h)
    return hess


def _phi(x) -> np.ndarray:
    """The standard normal distribution function (Abramowitz & Stegun 7.1.26 for erf; |error| < 1.5·10⁻⁷)."""
    x = np.asarray(x, dtype=float) / math.sqrt(2.0)
    z = np.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * z)
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    erf = 1.0 - poly * np.exp(-z * z)
    return 0.5 * (1.0 + np.sign(x) * erf)


def _flight(distances_um, beta, offset_um: float, cos_recoil: float, proper: bool) -> np.ndarray:
    d = np.asarray(distances_um, dtype=float) + float(offset_um)
    if np.any(d < 0):
        raise ValueError("a distance plus the offset is negative")
    return np.asarray(lt.flight_times(d, beta, cos_recoil, proper=proper), dtype=float)


def _arr(x, n: int, name: str) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    if a.shape != (n,):
        raise ValueError(f"{name} needs one value per distance ({n})")
    return a


# ---------------------------------------------------------------------------------------------------------------
# Peak areas


def _spectrum(item, edges_kev) -> tuple:
    if isinstance(item, (tuple, list)) and len(item) == 2 and np.ndim(item[0]) == 1 and np.ndim(item[1]) == 1:
        counts, edges = item
    else:
        if edges_kev is None:
            raise ValueError("give each spectrum as (counts, edges_kev), or the shared edges_kev")
        counts, edges = item, edges_kev
    counts = np.asarray(counts, dtype=float)
    edges = np.asarray(edges, dtype=float)
    if edges.shape != (len(counts) + 1,) or np.any(np.diff(edges) <= 0):
        raise ValueError("a spectrum needs increasing bin edges, one more than its bins")
    return counts, edges


def _doublet(counts, edges, e_u, e_s, sigma0, half_width) -> tuple:
    """Two Gaussians (unshifted at e_u, shifted at e_s, their own widths) on a straight line, fitted by weighted
    least squares (weights 1/max(counts, 1)). Returns ((area_s, area_u), (err_s, err_u))."""
    centres = 0.5 * (edges[1:] + edges[:-1])
    fwhm = FWHM_PER_SIGMA * sigma0
    lo, hi = min(e_u, e_s) - half_width * fwhm, max(e_u, e_s) + half_width * fwhm
    m = (centres >= lo) & (centres <= hi)
    if m.sum() < 10:
        raise ValueError("the fit window holds too few bins; widen it or use finer bins")
    x, y, w = centres[m], counts[m], np.diff(edges)[m]
    origin = 0.5 * (e_u + e_s)
    wt = 1.0 / np.sqrt(np.maximum(y, 1.0))
    edge_level = float(np.median(np.concatenate([y[:3], y[-3:]])))

    def near(e):
        k = np.abs(x - e) <= fwhm
        return max(float(np.sum(y[k] - edge_level)), 1.0)

    def model(p):
        a_u, a_s, mu_u, mu_s, s_u, s_s, b0, b1 = p
        g_u = np.exp(-0.5 * ((x - mu_u) / s_u) ** 2) / (s_u * math.sqrt(2 * math.pi))
        g_s = np.exp(-0.5 * ((x - mu_s) / s_s) ** 2) / (s_s * math.sqrt(2 * math.pi))
        return w * (a_u * g_u + a_s * g_s) + b0 + b1 * (x - origin)

    def resid(p):
        return (y - model(p)) * wt

    shift_room = max(fwhm, 0.15 * abs(e_s - e_u))
    p0 = [near(e_u), near(e_s), e_u, e_s, sigma0, sigma0, edge_level, 0.0]
    lower = [0.0, 0.0, e_u - fwhm, e_s - shift_room, 0.3 * sigma0, 0.3 * sigma0, -np.inf, -np.inf]
    upper = [np.inf, np.inf, e_u + fwhm, e_s + shift_room, 4 * sigma0, 6 * sigma0, np.inf, np.inf]
    p, cov, _ = _least_squares(resid, p0, lower, upper)
    err = np.sqrt(np.maximum(np.diag(cov), 0.0))
    return (float(p[1]), float(p[0])), (float(err[1]), float(err[0]))


def peak_areas(spectra_by_distance, e0_kev, beta, theta_deg, fwhm_kev, *, edges_kev=None,
               half_width: float = 4.0, fwhm_shifted_kev: Optional[float] = None) -> list:
    """The shifted and unshifted peak areas (counts) in each spectrum: ``[(shifted, unshifted, (err_shifted,
    err_unshifted)), ...]``, one per distance.

    ``spectra_by_distance``: each a (counts, edges_kev) pair, or counts with the shared ``edges_kev``. The unshifted
    peak is at ``e0_kev``, the shifted one at E₀ + :func:`lifetimes.peak_separation` (β, the detector at
    ``theta_deg`` to the recoil direction). Peaks further apart than 2 × ``half_width`` × FWHM are fitted one by one
    with :func:`physim.nuclear.analysis.fit_peak` (a Gaussian on a straight line, the area counted within ±3σ);
    closer ones together, two Gaussians with their own widths on one straight line (the shifted peak is broadened by
    the spread of velocities and angles: ``fwhm_shifted_kev`` is its starting width)."""
    if not fwhm_kev > 0:
        raise ValueError("the detector's FWHM must be positive")
    sep = float(lt.peak_separation(float(e0_kev), float(beta), float(theta_deg)))
    e_u, e_s = float(e0_kev), float(e0_kev) + sep
    sigma0 = fwhm_kev / FWHM_PER_SIGMA
    sigma_s = (fwhm_shifted_kev or fwhm_kev) / FWHM_PER_SIGMA
    out = []
    for item in spectra_by_distance:
        counts, edges = _spectrum(item, edges_kev)
        if abs(sep) >= 2 * half_width * fwhm_kev:
            fu = fit_peak(counts, edges, e_u, sigma0, half_width)
            fs = fit_peak(counts, edges, e_s, sigma_s, half_width, max_sigma=6 * sigma_s)
            out.append((fs.area, fu.area, (fs.area_unc, fu.area_unc)))
        else:
            (a_s, a_u), (s_s, s_u) = _doublet(counts, edges, e_u, e_s, max(sigma0, 1e-9), half_width)
            out.append((a_s, a_u, (s_s, s_u)))
    return out


# ---------------------------------------------------------------------------------------------------------------
# Decay curves


def _taus_with(feed: Feeding, values: dict) -> list:
    taus = [float(x) for x in feed.taus_ps]
    for i, v in values.items():
        taus[i] = float(v)
    return taus


def _decay_fit(t, r, sr, feed: Feeding, free_norm: bool, guess: Optional[float]) -> tuple:
    """(τ, σ_τ, χ², ndf, fitted {index: τ}, norm, model) for the unshifted share r(t) of the measured level."""
    lev = int(feed.level)
    idx = [lev] + [int(i) for i in feed.free]
    tau0 = guess
    if tau0 is None:  # from the data: where the share falls to 1/e, roughly
        good = (r > 0) & (t > 0) & (r < 0.99)
        tau0 = float(np.median(-t[good] / np.log(np.clip(r[good], 1e-6, 0.99)))) if good.any() else 10.0
        if not np.isfinite(tau0) or tau0 <= 0:
            tau0 = 10.0
    p0 = [math.log(tau0)] + [math.log(max(float(feed.taus_ps[i]), 1e-6)) for i in feed.free]
    if free_norm:
        p0.append(1.0)

    def model(p):
        taus = _taus_with(feed, {i: math.exp(v) for i, v in zip(idx, p[:len(idx)])})
        norm = p[len(idx)] if free_norm else 1.0
        return norm * lt.unshifted(taus, feed.feeding, feed.populations, t)[lev]

    def resid(p):
        return (r - model(p)) / sr

    lower = [math.log(1e-6)] * len(idx) + ([0.0] if free_norm else [])
    upper = [math.log(1e9)] * len(idx) + ([np.inf] if free_norm else [])
    p, cov, chi2 = _least_squares(resid, p0, lower, upper)
    ndf = len(r) - len(p)
    tau = math.exp(p[0])
    sigma = tau * math.sqrt(max(cov[0, 0], 0.0))
    fitted = {i: math.exp(v) for i, v in zip(idx, p[:len(idx)])}
    return tau, sigma, chi2, ndf, fitted, (p[len(idx)] if free_norm else 1.0), model(p)


def decay_curve(distances_um, unshifted, errors, beta, feeding: Optional[Feeding] = None, *, offset_um: float = 0.0,
                offset_error_um: float = 0.0, beta_error: float = 0.0, free_norm: bool = False, proper: bool = True,
                cos_recoil: float = 1.0) -> TauResult:
    """τ from a fit of the unshifted share against distance: R(d) = I_u/(I_u + I_s), with ``errors`` its uncertainty
    at each distance. Without ``feeding``, R = exp(−t/τ); with a :class:`Feeding`, R is the measured level's
    :func:`lifetimes.unshifted` for the cascade (Bateman's equations), the feeders' τ fixed except those in
    ``feeding.free``. t = (d + ``offset_um``)/(βc cos θ), divided by γ when ``proper``.

    ``free_norm`` fits a normalisation too (R = a·exp(−t/τ) as Radeck et al. fit), for shares that do not reach 1 at
    d = 0. The statistical uncertainty comes from the fit's covariance, scaled by √(χ²/ndf) when that is above 1.
    The budget: ``offset_error_um`` (distance), ``beta_error`` (absolute σ of β: velocity), the feeders'
    ``tau_errors_ps`` (feeding); stopping does not enter a flight time (0)."""
    d = np.asarray(distances_um, dtype=float).ravel()
    n = len(d)
    r = _arr(unshifted, n, "unshifted")
    sr = _arr(errors, n, "errors")
    if np.any(~(sr > 0)):
        raise ValueError("every error must be positive")
    feed = (feeding or _one_level()).checked()
    if n <= 1 + len(feed.free) + int(free_norm):
        raise ValueError("too few distances for the parameters fitted")

    def run(off=offset_um, b=beta, f=feed, guess=None):
        t = _flight(d, b, off, cos_recoil, proper)
        return _decay_fit(t, r, sr, f, free_norm, guess)

    tau, sigma, chi2, ndf, fitted, norm, model = run()
    notes = []
    red = chi2 / ndf if ndf > 0 else float("nan")
    if ndf > 0 and red > 1:
        sigma *= math.sqrt(red)
        notes.append(f"χ²/ndf = {red:.2f}: the statistical uncertainty is scaled by its square root.")
    distance = velocity = 0.0
    if offset_error_um > 0:
        distance = 0.5 * abs(run(off=offset_um + offset_error_um, guess=tau)[0]
                             - run(off=offset_um - offset_error_um, guess=tau)[0])
    if beta_error > 0:
        velocity = 0.5 * abs(run(b=beta + beta_error, guess=tau)[0] - run(b=beta - beta_error, guess=tau)[0])
    feeding_err = 0.0
    for i, e in enumerate(feed.errors()):
        if e > 0 and i != feed.level and i not in feed.free:
            shifts = []
            for sgn in (1, -1):
                taus = list(feed.taus_ps)
                taus[i] = max(taus[i] + sgn * e, 1e-9)
                f2 = Feeding(taus, feed.feeding, feed.populations, feed.level, feed.free, None)
                shifts.append(run(f=f2, guess=tau)[0])
            feeding_err = math.hypot(feeding_err, 0.5 * abs(shifts[0] - shifts[1]))
    notes.append("Stopping does not enter the flight times of a plunger: its budget entry is 0.")
    t = _flight(d, beta, offset_um, cos_recoil, proper)
    per = [{"distance_um": float(d[k]), "flight_time_ps": float(t[k]), "unshifted": float(r[k]),
            "error": float(sr[k]), "model": float(model[k])} for k in range(n)]
    params = {"chi2": float(chi2), "ndf": int(ndf), "norm": float(norm)}
    for i, v in fitted.items():
        if i != feed.level:
            params[f"tau_{i}_ps"] = float(v)
    return TauResult(float(tau), float(sigma), _budget(sigma, distance=distance, velocity=velocity,
                                                        feeding=feeding_err), "decay curve", per, notes, params)


# ---------------------------------------------------------------------------------------------------------------
# DDCM


def _curve_design(t: np.ndarray, polynomial, knots_t: np.ndarray, scale: float) -> tuple:
    """(X, dX/dt) of the curve fitted to the shifted intensity against flight time: a single polynomial in t of
    the given degree (as NaPyTau), quadratic pieces joined with a continuous slope at the knots ("piecewise", as
    NAPATAU), or a sum of decays with fixed time constants, log-spaced over the flight times ("smooth": the decay
    curve of any cascade is a sum of exponentials, Bateman)."""
    x = t / scale
    if isinstance(polynomial, (int, np.integer)) and not isinstance(polynomial, bool):
        k = int(polynomial)
        if k < 1:
            raise ValueError("the polynomial's degree must be at least 1")
        cols = [x ** j for j in range(k + 1)]
        dcols = [j * x ** (j - 1) / scale if j else np.zeros_like(x) for j in range(k + 1)]
    elif polynomial == "piecewise":
        cols = [np.ones_like(x), x, x * x]
        dcols = [np.zeros_like(x), np.ones_like(x) / scale, 2 * x / scale]
        for kn in knots_t / scale:
            cols.append(np.where(x > kn, (x - kn) ** 2, 0.0))
            dcols.append(np.where(x > kn, 2 * (x - kn), 0.0) / scale)
    elif polynomial == "smooth":
        pos = t[t > 0]
        if len(pos) < 2:
            raise ValueError("the smooth curve needs at least two distances above zero")
        m = int(min(5, max(2, len(t) - 2)))
        cols, dcols = [np.ones_like(t)], [np.zeros_like(t)]
        for theta in np.geomspace(float(pos.min()), float(t.max()), m):
            cols.append(-np.expm1(-t / theta))
            dcols.append(np.exp(-t / theta) / theta)
    else:
        raise ValueError("polynomial is 'piecewise', 'smooth' or a degree (an integer)")
    return np.column_stack(cols), np.column_stack(dcols)


def _default_knots(t: np.ndarray, per_piece: int = 2) -> np.ndarray:
    """Knots between the points (at their geometric mean), so that each quadratic piece holds about ``per_piece``
    of them: NAPATAU's short pieces. Two per piece follow a decay curve over two decades of distance to a few
    0.1% when the unshifted points constrain the slope too (the joint fit)."""
    ts = np.sort(np.unique(t))
    n_pieces = max(1, len(ts) // per_piece)
    if n_pieces == 1:
        return np.array([])
    out = []
    for c in np.linspace(0, len(ts), n_pieces + 1)[1:-1]:
        i = min(max(int(round(c)), 1), len(ts) - 1)
        out.append(math.sqrt(ts[i - 1] * ts[i]) if ts[i - 1] > 0 else 0.5 * ts[i])
    return np.array(sorted(set(out)))


def _tau_points(num, s_num, slope, dx_mat, cov) -> tuple:
    """τ(x) = num/slope at each point and its covariance (the numerator's error and the curve's)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        tau_i = num / slope
        jac = -(tau_i / slope)[:, None] * dx_mat
        c = np.diag((s_num / slope) ** 2) + jac @ cov @ jac.T
    err_i = np.sqrt(np.maximum(np.diag(c), 0.0))
    return tau_i, err_i, c


def _ddcm_core(t, num, s_num, der, s_der, polynomial, knots_t, region, weighted: bool, fit: str = "joint") -> dict:
    """The DDCM on normalised intensities: ``num`` the decays after each point (unshifted, feeding subtracted),
    ``der`` the decays before it (shifted), against flight time t.

    ``fit="joint"`` (NAPATAU's way): one τ for the region, and the curve P(t) fitted to the shifted points and,
    as τ·P′(t), to the numerators at once; for each trial τ the curve is a linear fit, so τ comes from a profile
    χ²(τ) and its uncertainty where χ² rises by 1 (by √(χ²/ndf) more when that is above 1). ``fit="shifted"``
    (NaPyTau's way): the curve fitted to the shifted points alone, τ(x) = num/P′ at each point, and τ their mean
    weighted by their errors, with its variance from their full covariance (the curve correlates them).
    Either way the τ(x) of every point are returned, from the fitted curve."""
    sel = np.ones(len(t), dtype=bool)
    if region is not None:
        sel = (t >= region[0]) & (t <= region[1])
    if sel.sum() < 3:
        raise ValueError("fewer than three distances in the region")
    ts = t[sel]
    scale = float(np.max(ts)) if np.max(ts) > 0 else 1.0
    if polynomial == "piecewise":
        knots = _default_knots(ts) if knots_t is None else np.asarray(knots_t, dtype=float)
    else:
        knots = np.array([])
    x_all, dx_all = _curve_design(t, polynomial, knots, scale)
    x_mat, dx_mat = x_all[sel], dx_all[sel]
    p = x_mat.shape[1]
    if fit == "shifted":
        if p >= len(ts):
            raise ValueError("too few distances for the curve: use fewer knots or a lower degree")
        w = 1.0 / s_der[sel] ** 2 if weighted else np.ones(len(ts))
        a_inv = np.linalg.pinv(x_mat.T @ (x_mat * w[:, None]))
        coef = a_inv @ (x_mat.T @ (w * der[sel]))
        chi2 = float(np.sum((der[sel] - x_mat @ coef) ** 2 / s_der[sel] ** 2))
        ndf = len(ts) - p
        # Sandwich covariance of the coefficients (right for the unweighted fit too), scaled when χ²/ndf > 1.
        cov = a_inv @ (x_mat.T @ (x_mat * (w * w * s_der[sel] ** 2)[:, None])) @ a_inv
        if ndf > 0 and chi2 / ndf > 1:
            cov = cov * (chi2 / ndf)
        slope = dx_all @ coef
        tau_i, err_i, c = _tau_points(num, s_num, slope, dx_all, cov)
        use = sel & np.isfinite(tau_i) & (tau_i > 0) & (slope > 0) & np.isfinite(err_i) & (err_i > 0)
        if not use.any():
            raise ValueError("no distance gives a usable τ(x): check the intensities and the region")
        k = np.flatnonzero(use)
        wi = 1.0 / err_i[k] ** 2
        wi /= wi.sum()
        tau = float(wi @ tau_i[k])
        stat = math.sqrt(max(float(wi @ c[np.ix_(k, k)] @ wi), 0.0))
        plus = minus = stat
    elif fit == "joint":
        if 2 * len(ts) <= p + 1:
            raise ValueError("too few distances for the curve: use fewer knots or a lower degree")
        ws, wu = 1.0 / s_der[sel], 1.0 / s_num[sel]
        b = np.concatenate([der[sel] * ws, num[sel] * wu])

        def solve(tau):
            a = np.vstack([x_mat * ws[:, None], tau * dx_mat * wu[:, None]])
            coef = np.linalg.lstsq(a, b, rcond=None)[0]
            r = a @ coef - b
            return coef, float(r @ r), a

        def chi2_of(ln_tau):
            return solve(math.exp(ln_tau))[1]

        with np.errstate(divide="ignore", invalid="ignore"):
            rough = num[sel] / np.gradient(der[sel], ts) if len(ts) > 1 else np.array([1.0])
        rough = rough[np.isfinite(rough) & (rough > 0)]
        g = float(np.median(rough)) if len(rough) else float(np.median(ts))
        x0, f0, _ = _minimise_1d(chi2_of, math.log(g) - 4.0, math.log(g) + 4.0, n_scan=41, tol=1e-9)
        ndf = 2 * len(ts) - p - 1
        red = f0 / ndf if ndf > 0 else 1.0
        target = f0 + max(red, 1.0)  # Δχ² = 1, scaled by χ²/ndf when that is above 1
        hi = _crossing(chi2_of, x0, target, +1, x0 + 8.0)
        lo = _crossing(chi2_of, x0, target, -1, x0 - 8.0)
        tau = math.exp(x0)
        plus = math.exp(hi) - tau if hi is not None else float("inf")
        minus = tau - math.exp(lo) if lo is not None else tau
        stat = 0.5 * (plus + minus) if math.isfinite(plus) else minus
        coef, chi2, a = solve(tau)
        cov = np.linalg.pinv(a.T @ a) * max(red, 1.0)
        slope = dx_all @ coef
        tau_i, err_i, c = _tau_points(num, s_num, slope, dx_all, cov)
        use = sel & np.isfinite(tau_i) & (tau_i > 0) & (slope > 0) & np.isfinite(err_i) & (err_i > 0)
    else:
        raise ValueError("fit is 'joint' or 'shifted'")
    k = np.flatnonzero(use)
    spread = float(np.sum((tau_i[k] - tau) ** 2 / err_i[k] ** 2) / max(len(k) - 1, 1)) if len(k) else float("nan")
    return {"tau": tau, "stat": stat, "plus": plus, "minus": minus, "tau_i": tau_i, "err_i": err_i, "use": use,
            "slope": slope, "curve": x_all @ coef, "chi2": chi2, "ndf": ndf, "spread": spread, "knots_t": knots,
            "fit": fit}


def _ddcm_inputs(n, unshifted, shifted, errors, feeder_unshifted, feeder_shifted, norm, normalise) -> dict:
    u = _arr(unshifted, n, "unshifted")
    s = _arr(shifted, n, "shifted")
    errs = list(errors) if isinstance(errors, (tuple, list)) else list(np.asarray(errors, dtype=float).T)
    if len(errs) not in (2, 4):
        raise ValueError("errors are (σ_unshifted, σ_shifted), plus (σ_feeder_unshifted, σ_feeder_shifted) "
                         "with a feeder")
    su, ss = _arr(errs[0], n, "errors"), _arr(errs[1], n, "errors")
    out = {"u": u, "s": s, "su": su, "ss": ss}
    if feeder_unshifted is not None:
        bu = _arr(feeder_unshifted, n, "feeder_unshifted")
        sbu = _arr(errs[2], n, "errors") if len(errs) == 4 else np.sqrt(np.maximum(np.abs(bu), 1.0))
        out.update(bu=bu, sbu=sbu)
    if feeder_shifted is not None:
        bs = _arr(feeder_shifted, n, "feeder_shifted")
        sbs = _arr(errs[3], n, "errors") if len(errs) == 4 else np.sqrt(np.maximum(np.abs(bs), 1.0))
        out.update(bs=bs, sbs=sbs)
    if norm is not None:
        nrm = _arr(norm, n, "norm")
    elif normalise:
        nrm = u + s
    else:
        nrm = np.ones(n)
    if np.any(~(nrm > 0)):
        raise ValueError("the normalisation must be positive at every distance")
    for key in ("u", "s", "su", "ss", "bu", "sbu", "bs", "sbs"):
        if key in out:
            out[key] = out[key] / nrm
    return out


def _ddcm_result(t_of, d_eff, data, alpha, gated, polynomial, knots_um, region_um, weighted, fit, *, beta,
                 beta_error, offset_um, offset_error_um, alpha_error, method, extra_velocity=None) -> TauResult:
    """The DDCM with its budget; ``t_of(offset, beta, ...)`` gives the flight times of the points."""
    use_feeder = (not gated) and ("bu" in data)

    def run(off=offset_um, b=beta, a=alpha, **kw):
        t = t_of(off, b, **kw)
        num = data["u"] - (a * data["bu"] if use_feeder else 0.0)
        s_num = np.sqrt(data["su"] ** 2 + ((a * data["sbu"]) ** 2 if use_feeder else 0.0))
        knots_t = None if knots_um is None else t_of(off, b, distances=np.asarray(knots_um, dtype=float))
        region = None if region_um is None else t_of(off, b, distances=np.asarray(region_um, dtype=float))
        return _ddcm_core(t, num, s_num, data["s"], data["ss"], polynomial, knots_t, region, weighted, fit), t

    core, t = run()
    tau, stat = core["tau"], core["stat"]
    notes = []
    if core["spread"] > 1.5:
        notes.append(f"The τ(x) scatter more than their errors allow (χ²/ndf = {core['spread']:.2f} about the "
                     "mean): look at the τ(x) plot and the region of sensitivity.")
    distance = velocity = feeding = 0.0
    if offset_error_um > 0:
        distance = 0.5 * abs(run(off=offset_um + offset_error_um)[0]["tau"]
                             - run(off=offset_um - offset_error_um)[0]["tau"])
    if beta_error > 0:
        velocity = 0.5 * abs(run(b=beta + beta_error)[0]["tau"] - run(b=beta - beta_error)[0]["tau"])
    if extra_velocity is not None:
        velocity = math.hypot(velocity, extra_velocity(run))
    if use_feeder and alpha_error > 0:
        feeding = 0.5 * abs(run(a=alpha + alpha_error)[0]["tau"] - run(a=alpha - alpha_error)[0]["tau"])
    if gated:
        notes.append("Gated on the shifted component of the feeding transition: no feeding correction is needed.")
    elif not use_feeder:
        notes.append("No feeder given: τ(x) is without a feeding correction (right only for a level fed directly).")
    notes.append("Stopping does not enter the flight times of a plunger: its budget entry is 0.")
    per = []
    for k in range(len(t)):
        row = {"distance_um": float(d_eff[k]), "flight_time_ps": float(t[k]), "tau_ps": float(core["tau_i"][k]),
               "error_ps": float(core["err_i"][k]), "used": bool(core["use"][k]),
               "unshifted": float(data["u"][k]), "shifted": float(data["s"][k]), "curve": float(core["curve"][k])}
        if "bu" in data:
            row["feeder_unshifted"] = float(data["bu"][k])
        if "bs" in data:
            row["feeder_shifted"] = float(data["bs"][k])
        per.append(row)
    params = {"chi2": core["chi2"], "ndf": core["ndf"], "spread": core["spread"], "fit": fit,
              "plus_ps": core["plus"], "minus_ps": core["minus"], "polynomial": polynomial,
              "knots_ps": [float(x) for x in core["knots_t"]]}
    if use_feeder:
        params["alpha"] = float(alpha)
    return TauResult(float(tau), float(stat), _budget(stat, distance=distance, velocity=velocity, feeding=feeding),
                     method, per, notes, params)


def ddcm(distances_um, unshifted, shifted, errors, beta, feeder_unshifted=None, feeder_shifted=None, alpha=1.0,
         polynomial="piecewise", *, gated: bool = False, norm=None, normalise: bool = True, alpha_error: float = 0.0,
         offset_um: float = 0.0, offset_error_um: float = 0.0, beta_error: float = 0.0, knots_um=None,
         region_um=None, weighted: bool = True, proper: bool = True, cos_recoil: float = 1.0,
         fit: str = "joint") -> TauResult:
    """The differential decay-curve method (Dewald et al., Z. Phys. A 334, 163 (1989); Böhm et al., NIM A 329, 248
    (1993) for the gated form): at each distance

        τ(x) = [I_u(x) − α I_u^B(x)] / (dI_s/dt)(x)     (singles, B the feeding transition), or
        τ(x) = {B_s, A_u}(x) / (d{B_s, A_s}/dt)(x)      (``gated``: A's components in coincidence with B's shifted one),

    the derivative from a curve P(t) through the shifted intensity against flight time t = (x + ``offset_um``)/(βc):
    ``polynomial="piecewise"`` (quadratic pieces with a continuous slope, as NAPATAU; knots between every two
    points, or ``knots_um``), ``"smooth"`` (a sum of decays with fixed time constants log-spaced over the flight
    times: any cascade's decay curve is one) or an integer (one polynomial in t of that degree, as NaPyTau). The
    shifted intensity carries the derivative and the unshifted one the numerator, so their errors are independent.

    ``fit="joint"`` (default, NAPATAU's way): one τ over the distances in ``region_um`` (two distances; all by
    default), the curve fitted to the shifted points and, as τ·P′, to the numerators together; τ from the profile
    χ²(τ) and its uncertainty where χ² rises by 1 (scaled by χ²/ndf when above 1; asymmetric errors in ``params``).
    ``fit="shifted"`` (NaPyTau's way): the curve through the shifted points alone (``weighted=False``: unweighted,
    as NaPyTau), and τ the mean of the τ(x) weighted by their errors, its uncertainty from their covariance. The
    curve fitted to the shifted points alone follows a decay curve less well, and the mean inherits its bias (a few
    % for quadratic pieces), which the joint fit avoids. Either way ``per_distance`` holds every τ(x).

    ``errors``: (σ_unshifted, σ_shifted), plus (σ_feeder_unshifted, σ_feeder_shifted) with a feeder (√counts
    otherwise). Singles are normalised at each distance by I_u + I_s of the transition (``normalise``), or by
    ``norm`` (one positive number per distance: the beam count, a total of all lines); gated data need ``norm``
    (their I_u + I_s follows the feeder's decay) or are taken as given. ``alpha`` scales the feeder's intensity into
    the transition's: with efficiency-corrected intensities it is the share of the feeder's decays that go through
    the transition (1 when it is the only branch); ``feeder_shifted`` is kept for the plot only. Budget:
    ``offset_error_um``, ``beta_error`` (absolute σ of β), ``alpha_error``."""
    d = np.asarray(distances_um, dtype=float).ravel()
    n = len(d)
    if gated and feeder_unshifted is not None:
        raise ValueError("a gated DDCM takes no feeder intensities: the gate removes the feeding")
    if gated and norm is None and normalise:
        normalise = False
    data = _ddcm_inputs(n, unshifted, shifted, errors, feeder_unshifted, feeder_shifted, norm, normalise)

    def t_of(off, b, distances=None):
        return _flight(d if distances is None else distances, b, off, cos_recoil, proper)

    method = "DDCM gated" if gated else "DDCM"
    return _ddcm_result(t_of, d, data, float(alpha), gated, polynomial, knots_um, region_um, weighted, fit,
                        beta=beta, beta_error=beta_error, offset_um=offset_um, offset_error_um=offset_error_um,
                        alpha_error=alpha_error, method=method)


def ddcm_three_foil(distances_um, fully_shifted, reduced_shifted, unshifted, errors, beta, beta_reduced,
                    feeder_fully_shifted=None, feeder_reduced_shifted=None, feeder_unshifted=None, alpha=1.0,
                    polynomial="piecewise", *, gated: bool = False, norm=None, alpha_error: float = 0.0,
                    offset_um: float = 0.0, offset_error_um: float = 0.0, beta_error: float = 0.0,
                    beta_reduced_error: float = 0.0, region_um=None, weighted: bool = True,
                    proper: bool = True, fit: str = "joint") -> TauResult:
    """The DDCM for a differential (three-foil) plunger: target, a gap x₁, a degrader that slows the recoil from β
    to β′, a gap x₂, a stopper (Beckers et al., NIM A 1042, 167418 (2022), Eqs. 13–14; CoCoDiff, ⁵⁰Cr). Each pair
    of distances gives three peaks: fully shifted (decays in the first gap), reduced shifted (in the second) and
    unshifted (in the stopper), and so two points of one decay curve:

    - at the degrader, t₁ = x₁/(βc): decays after it, I_rs + I_u, against decays before it, I_fs;
    - at the stopper, t₂ = x₁/(βc) + x₂/(β′c): decays after it, I_u, against decays before it, I_fs + I_rs.

    Both are put on one axis by renormalising the second gap by v/v′: x_eff = x₁ + x₂ β/β′, flown at β. Then the
    2n points go through :func:`ddcm` (singles with the feeder's three components and ``alpha``, or ``gated``).
    ``distances_um``: (n, 2) pairs (x₁, x₂). ``errors``: (σ_fs, σ_rs, σ_u), plus the feeder's three. The
    normalisation at each pair is the transition's I_fs + I_rs + I_u (or ``norm``). Budget: ``beta_error`` and
    ``beta_reduced_error`` (velocity), ``offset_error_um`` (on x₁), ``alpha_error``."""
    pairs = np.asarray(distances_um, dtype=float)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError("give the distances as (n, 2) pairs: target–degrader and degrader–stopper")
    n = len(pairs)
    fs = _arr(fully_shifted, n, "fully_shifted")
    rs = _arr(reduced_shifted, n, "reduced_shifted")
    uu = _arr(unshifted, n, "unshifted")
    errs = list(errors)
    if len(errs) not in (3, 6):
        raise ValueError("errors are (σ_fs, σ_rs, σ_u), plus the feeder's three")
    sfs, srs, su = (_arr(e, n, "errors") for e in errs[:3])
    if not 0 < beta_reduced < beta:
        raise ValueError("the degrader slows the recoil: 0 < β′ < β")
    if norm is not None:
        nrm = _arr(norm, n, "norm")
    elif gated:
        nrm = np.ones(n)
    else:
        nrm = fs + rs + uu
    # The 2n points: (after, before) at the degrader, then at the stopper.
    after = np.concatenate([rs + uu, uu])
    before = np.concatenate([fs, fs + rs])
    s_after = np.concatenate([np.hypot(srs, su), su])
    s_before = np.concatenate([sfs, np.hypot(sfs, srs)])
    nrm2 = np.concatenate([nrm, nrm])
    feeder = None
    if feeder_unshifted is not None or feeder_reduced_shifted is not None:
        if gated:
            raise ValueError("a gated DDCM takes no feeder intensities: the gate removes the feeding")
        if feeder_unshifted is None or feeder_reduced_shifted is None:
            raise ValueError("the feeding correction needs the feeder's reduced-shifted and unshifted intensities")
        bu = _arr(feeder_unshifted, n, "feeder_unshifted")
        brs = _arr(feeder_reduced_shifted, n, "feeder_reduced_shifted")
        if len(errs) == 6:
            sbrs, sbu = _arr(errs[4], n, "errors"), _arr(errs[5], n, "errors")
        else:
            sbrs, sbu = np.sqrt(np.maximum(np.abs(brs), 1.0)), np.sqrt(np.maximum(np.abs(bu), 1.0))
        feeder = (np.concatenate([brs + bu, bu]), np.concatenate([np.hypot(sbrs, sbu), sbu]))
    data = {"u": after / nrm2, "s": before / nrm2, "su": s_after / nrm2, "ss": s_before / nrm2}
    if feeder is not None:
        data.update(bu=feeder[0] / nrm2, sbu=feeder[1] / nrm2)
    if feeder_fully_shifted is not None:  # for the plot: the feeder's decays before each foil
        bfs = _arr(feeder_fully_shifted, n, "feeder_fully_shifted")
        brs_plot = (np.zeros(n) if feeder_reduced_shifted is None
                    else _arr(feeder_reduced_shifted, n, "feeder_reduced_shifted"))
        data["bs"] = np.concatenate([bfs, bfs + brs_plot]) / nrm2
    x1, x2 = pairs[:, 0], pairs[:, 1]

    def t_of(off, b, b_red=beta_reduced, distances=None):
        if distances is not None:  # region limits, given as effective distances at β
            return _flight(distances, b, off, 1.0, proper)
        x_eff = np.concatenate([x1, x1 + x2 * b / b_red])
        return _flight(x_eff, b, off, 1.0, proper)

    def reduced_velocity(run) -> float:
        if beta_reduced_error <= 0:
            return 0.0
        hi = run(b_red=beta_reduced + beta_reduced_error)[0]["tau"]
        lo = run(b_red=max(beta_reduced - beta_reduced_error, 1e-6))[0]["tau"]
        return 0.5 * abs(hi - lo)

    d_eff = np.concatenate([x1, x1 + x2 * beta / beta_reduced])
    res = _ddcm_result(t_of, d_eff, data, float(alpha), gated, polynomial, None, region_um, weighted, fit, beta=beta,
                       beta_error=beta_error, offset_um=offset_um, offset_error_um=offset_error_um,
                       alpha_error=alpha_error, method="DDCM three-foil", extra_velocity=reduced_velocity)
    for k, row in enumerate(res.per_distance):
        row["foil"] = "degrader" if k < n else "stopper"
        row["pair"] = [float(x1[k % n]), float(x2[k % n])]
    res.notes.append("Three-foil DDCM: the second gap is renormalised by v/v′ (x_eff = x₁ + x₂ β/β′), so the "
                     "degrader and stopper points lie on one decay curve.")
    return res


# ---------------------------------------------------------------------------------------------------------------
# DSAM line shapes


class _ScaledKernel:
    """A kernel with its time axis and its speeds scaled (stopping and velocity budgets): the same recoils, slowing
    ``time_scale`` times slower, moving ``beta_scale`` times faster."""

    def __init__(self, kernel, time_scale: float = 1.0, beta_scale: float = 1.0):
        self.times_ps = np.asarray(kernel.times_ps, dtype=float) * time_scale
        self.weights = getattr(kernel, "weights", None)
        if lt._directional(kernel):
            self.beta = np.asarray(kernel.beta, dtype=float) * beta_scale
            self._kernel, self._scale = kernel, beta_scale
            self.beta_cos = self._beta_cos
        else:
            self.beta_cos = np.asarray(kernel.beta_cos, dtype=float) * beta_scale

    def _beta_cos(self, theta_deg, phi_deg=0.0):
        return np.asarray(self._kernel.beta_cos(theta_deg, phi_deg), dtype=float) * self._scale


def _fwhm_at(fwhm_kev, energy) -> np.ndarray:
    if fwhm_kev is None:
        return np.zeros_like(np.asarray(energy, dtype=float))
    if callable(fwhm_kev):
        return np.asarray(fwhm_kev(np.asarray(energy, dtype=float)), dtype=float)
    return np.full_like(np.asarray(energy, dtype=float), float(fwhm_kev))


def _response(edges: np.ndarray, e0: float, fwhm_kev, line_lo: float, line_hi: float) -> tuple:
    """A fine grid for the line (bins ≤ σ/3, covering the data and the Doppler range, plus 5σ) and the matrix that
    folds it with the detector's Gaussian and puts it in the data's bins."""
    widths = np.diff(edges)
    sigma = float(np.max(_fwhm_at(fwhm_kev, np.array([e0])))) / FWHM_PER_SIGMA
    step = float(np.min(widths)) / 4.0
    if sigma > 0:
        step = min(step, sigma / 3.0)
    lo = min(edges[0], line_lo) - 5 * sigma - step
    hi = max(edges[-1], line_hi) + 5 * sigma + step
    nf = int(math.ceil((hi - lo) / step))
    fine = lo + step * np.arange(nf + 1)
    centres = 0.5 * (fine[1:] + fine[:-1])
    if sigma > 0:
        sig = _fwhm_at(fwhm_kev, centres) / FWHM_PER_SIGMA
        z_hi = (edges[1:, None] - centres[None, :]) / sig[None, :]
        z_lo = (edges[:-1, None] - centres[None, :]) / sig[None, :]
        mat = _phi(z_hi) - _phi(z_lo)
    else:
        idx = np.searchsorted(edges, centres, side="right") - 1
        mat = np.zeros((len(edges) - 1, len(centres)))
        ok = (idx >= 0) & (idx < len(edges) - 1)
        mat[idx[ok], np.flatnonzero(ok)] = 1.0
    return fine, mat


def _background_design(edges: np.ndarray, e0: float, background) -> tuple:
    """(columns of the fitted background (bins, k), a fixed background (bins,))."""
    centres = 0.5 * (edges[1:] + edges[:-1])
    nb = len(centres)
    if isinstance(background, str):
        if background == "linear":
            span = max(float(edges[-1] - edges[0]), 1e-9)
            return np.column_stack([np.ones(nb), (centres - e0) / span]), np.zeros(nb)
        if background == "constant":
            return np.ones((nb, 1)), np.zeros(nb)
        if background == "none":
            return np.zeros((nb, 0)), np.zeros(nb)
        raise ValueError("background is 'linear', 'constant', 'none', or a fixed level (a number or one per bin)")
    fixed = np.broadcast_to(np.asarray(background, dtype=float), (nb,)).copy()
    return np.zeros((nb, 0)), fixed


def _linear_fit(counts, design, fixed, statistic: str) -> tuple:
    """The linear parameters (line area, background) for one shape: weighted least squares, iterated with weights
    1/model for Poisson statistics (maximum likelihood for a linear model). Returns (model, statistic value)."""
    y = counts - fixed
    if statistic == "neyman":
        w = 1.0 / np.maximum(counts, 1.0)
        coef = np.linalg.lstsq(design * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
        m = design @ coef + fixed
        return m, float(np.sum((counts - m) ** 2 * w))
    w = 1.0 / np.maximum(counts, 1.0)
    m = None
    for _ in range(12):
        coef = np.linalg.lstsq(design * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
        m_new = design @ coef + fixed
        if m is not None and np.max(np.abs(m_new - m)) < 1e-9 * (np.max(np.abs(m_new)) + 1.0):
            m = m_new
            break
        m = m_new
        w = 1.0 / np.maximum(m, 1e-3)
    mm = np.maximum(m, 1e-12)
    pos = counts > 0
    dev = 2.0 * float(np.sum(mm - counts) + np.sum(counts[pos] * np.log(counts[pos] / mm[pos])))
    return m, dev


def _lineshape_fit(counts, edges, e0, kernel, theta, feed: Feeding, background, fwhm_kev, statistic, tau_range,
                   start=None) -> dict:
    """Profile fit: the line area and background are solved for each trial τ (and side-feeding τ); τ by a scan
    and golden-section search in ln τ; its uncertainty where the statistic rises by 1."""
    t, beta, proj, _ = lt._kernel_arrays(kernel, theta)
    ratio = lt._ratio(beta, proj)
    line_lo, line_hi = e0 * float(np.min(ratio)), e0 * float(np.max(ratio))
    fine, fold = _response(edges, e0, fwhm_kev, line_lo, line_hi)
    bg_cols, fixed = _background_design(edges, e0, background)
    lev = int(feed.level)
    free = [int(i) for i in feed.free]

    def shape(tau, extra):
        taus = _taus_with(feed, {lev: tau, **extra})
        row = lt.lineshape(e0, taus, feed.feeding, feed.populations, kernel, theta, fine)[lev]
        return fold @ row

    def value(ln_tau, ln_extra=()):
        extra = {i: math.exp(v) for i, v in zip(free, ln_extra)}
        s = shape(math.exp(ln_tau), extra)
        design = np.column_stack([s, bg_cols]) if bg_cols.size else s[:, None]
        return _linear_fit(counts, design, fixed, statistic)

    lo, hi = math.log(tau_range[0]), math.log(tau_range[1])
    if start is not None:
        lo_s, hi_s = max(lo, math.log(start) - 1.0), min(hi, math.log(start) + 1.0)
    else:
        lo_s, hi_s = lo, hi
    if not free:
        def profile(ln_tau) -> float:
            return value(ln_tau)[1]

        x, fmin, edge = _minimise_1d(profile, lo_s, hi_s, n_scan=9 if start is not None else 25, tol=1e-6)
        upper = _crossing(profile, x, fmin + 1.0, +1, hi)
        lower = _crossing(profile, x, fmin + 1.0, -1, lo)
        model, stat_value = value(x)
        tau = math.exp(x)
        return {"tau": tau, "plus": None if upper is None else math.exp(upper) - tau,
                "minus": None if lower is None else tau - math.exp(lower), "statistic": stat_value, "model": model,
                "edge": edge, "side": {}}
    # Free side-feeding times: a coarse grid, then Nelder–Mead in (ln τ, ln τ_sf, ...); the uncertainty from the
    # curvature (statistic rising by 1 along the profile: σ² = 2 (H⁻¹)₀₀).
    side_lo, side_hi = math.log(max(tau_range[0], 1e-5)), math.log(tau_range[1])

    def fun(p) -> float:
        return value(p[0], p[1:])[1]

    if start is not None:
        side0 = [math.log(max(float(feed.taus_ps[i]), 1e-5)) for i in free]
        axes = [np.linspace(lo_s, hi_s, 7)] + [np.linspace(s - 1.5, s + 1.5, 5) for s in side0]
    else:
        axes = [np.linspace(lo, hi, 13)] + [np.linspace(side_lo, side_hi, 7) for _ in free]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, len(axes))
    values = np.array([fun(p) for p in grid])
    p0 = grid[int(np.argmin(values))]
    edge = p0[0] in (axes[0][0], axes[0][-1])
    p, fmin = _nelder_mead(fun, p0, 0.25)
    hess = _hessian(fun, p, 0.03)
    try:
        cov = 2.0 * np.linalg.inv(hess)
        sigma = math.sqrt(cov[0, 0]) if cov[0, 0] > 0 else float("nan")
    except np.linalg.LinAlgError:
        sigma = float("nan")
    model, stat_value = value(p[0], p[1:])
    tau = math.exp(p[0])
    ok = math.isfinite(sigma)
    return {"tau": tau, "plus": tau * math.expm1(sigma) if ok else None,
            "minus": tau * -math.expm1(-sigma) if ok else None, "statistic": stat_value, "model": model,
            "edge": edge, "side": {i: math.exp(v) for i, v in zip(free, p[1:])}}


def fit_lineshape(counts, edges_kev, e0_kev, kernel, theta_deg, feeding: Optional[Feeding] = None,
                  background="linear", *, fwhm_kev=None, statistic: str = "poisson", tau_range_ps=None,
                  stopping_uncertainty: float = STOPPING_UNCERTAINTY, beta_scale_error: float = 0.0) -> TauResult:
    """τ from a DSAM line shape: the counts in bins ``edges_kev`` against :func:`lifetimes.lineshape` of the level
    for the slowing-down ``kernel`` (``transport.kernel``, or any :class:`lifetimes.Kernel`; the detector at
    ``theta_deg`` to the beam), folded with the detector's Gaussian response (``fwhm_kev``: a number, or a function
    of energy in keV; None for none), plus a background ("linear", "constant", "none", or a fixed level: a number or
    one per bin). The kernel's binning does not depend on τ and is kept by :mod:`lifetimes`, so each trial τ is a
    re-weighting of it.

    ``feeding``: a :class:`Feeding` (``Feeding.side(...)`` for a side-feeding time, fixed or ``free=True``; the
    measured level's τ always fitted), or None for a level fed directly. The line area and the background are
    solved exactly at each trial τ (profiled); τ is found by a scan over ``tau_range_ps`` (default: from a twentieth
    of the kernel's first time step to ten times its last time) and a golden-section search, and its uncertainty is
    where the statistic rises by 1 (``params["plus_ps"]``, ``params["minus_ps"]``; ``stat_ps`` is their mean). With
    free side-feeding times: a grid, then Nelder–Mead, and the uncertainty from the curvature at the minimum (the
    profile's rise by 1 for a quadratic statistic), about a second per fit.
    ``statistic``: "poisson" (maximum likelihood) or "neyman" (χ² with the counts as variances, as AhKin).

    Budget: ``stopping_uncertainty`` scales the kernel's time axis by ± that share and refits; ``beta_scale_error``
    scales its speeds; the feeders' ``tau_errors_ps`` (fixed side-feeding times) give "feeding"."""
    counts = np.asarray(counts, dtype=float).ravel()
    edges = np.asarray(edges_kev, dtype=float).ravel()
    if edges.shape != (len(counts) + 1,) or np.any(np.diff(edges) <= 0):
        raise ValueError("the bin edges must increase, one more than the bins")
    if statistic not in ("poisson", "neyman"):
        raise ValueError("statistic is 'poisson' or 'neyman'")
    feed = (feeding or _one_level()).checked()
    t = np.asarray(kernel.times_ps, dtype=float)
    if tau_range_ps is None:
        first = float(t[1]) if len(t) > 1 else 1e-3
        tau_range_ps = (max(first / 20.0, 1e-5), 10.0 * float(t[-1]))
    e0, theta = float(e0_kev), float(theta_deg)

    def fit(k=kernel, f=feed, start=None):
        return _lineshape_fit(counts, edges, e0, k, theta, f, background, fwhm_kev, statistic, tau_range_ps, start)

    best = fit()
    tau = best["tau"]
    notes = []
    if best["edge"]:
        notes.append("The best τ is at the end of the range searched: widen tau_range_ps, or the line does not "
                     "constrain τ (a level much slower or faster than the slowing down).")
    plus, minus = best["plus"], best["minus"]
    if plus is None or minus is None:
        notes.append("The statistic does not rise by 1 on one side within the range: τ is a limit there.")
    errs = [e for e in (plus, minus) if e is not None]
    stat = float(np.mean(errs)) if errs else float("nan")
    stopping = velocity = feeding_err = 0.0
    if stopping_uncertainty > 0:
        hi = fit(_ScaledKernel(kernel, 1.0 + stopping_uncertainty), start=tau)["tau"]
        lo = fit(_ScaledKernel(kernel, 1.0 - stopping_uncertainty), start=tau)["tau"]
        stopping = 0.5 * abs(hi - lo)
    if beta_scale_error > 0:
        hi = fit(_ScaledKernel(kernel, 1.0, 1.0 + beta_scale_error), start=tau)["tau"]
        lo = fit(_ScaledKernel(kernel, 1.0, 1.0 - beta_scale_error), start=tau)["tau"]
        velocity = 0.5 * abs(hi - lo)
    for i, e in enumerate(feed.errors()):
        if e > 0 and i != feed.level and i not in feed.free:
            shifts = []
            for sgn in (1, -1):
                taus = list(feed.taus_ps)
                taus[i] = max(taus[i] + sgn * e, 1e-9)
                f2 = Feeding(taus, feed.feeding, feed.populations, feed.level, feed.free, None)
                shifts.append(fit(f=f2, start=tau)["tau"])
            feeding_err = math.hypot(feeding_err, 0.5 * abs(shifts[0] - shifts[1]))
    notes.append(f"Stopping: the kernel's time axis scaled by ±{100 * stopping_uncertainty:.0f}% (stopping powers "
                 "are known to about 5–10%).")
    notes.append("A line shape has no plunger distance: its budget entry is 0.")
    params = {"plus_ps": plus, "minus_ps": minus, "statistic": statistic, "statistic_min": best["statistic"],
              "bins": int(len(counts)), "model": best["model"].tolist()}
    for i, v in best["side"].items():
        params[f"tau_{i}_ps"] = float(v)
    return TauResult(float(tau), stat, _budget(stat, stopping=stopping, velocity=velocity, feeding=feeding_err),
                     "line shape", None, notes, params)


# ---------------------------------------------------------------------------------------------------------------
# B(Eλ) from a result


def b_from_result(result: TauResult, scheme, level, *, final=None, lam: Optional[str] = None) -> tuple:
    """(B(λ↓), its uncertainty) of a transition of the measured ``level`` (a name as
    :func:`physim.nuclear.gamma.level_names` gives, or an index) from the result's τ, through
    :func:`lifetimes.b_from_tau` with the level scheme's branching (γ-ray intensities and conversion), conversion and
    mixing. The transition is the one to ``final`` (name or index), or the level's strongest. ``lam`` is the
    multipole of B (default: the transition's multipolarity; E2 for an M1+E2 transition with a mixing ratio).

    The uncertainty adds τ's total (all of the budget) and the branching's from the intensities' uncertainties."""
    from .gamma import level_names

    names = level_names(scheme.levels)

    def index(x) -> int:
        if isinstance(x, (int, np.integer)):
            if not 0 <= int(x) < len(scheme.levels):
                raise ValueError(f"no level {x} in the scheme")
            return int(x)
        if x in names:
            return names.index(x)
        raise ValueError(f"no level '{x}' in the scheme (levels: {', '.join(names)})")

    i = index(level)
    branches = [t for t in scheme.transitions if t.initial == i and t.energy is not None and t.energy.value > 0]
    if not branches:
        raise ValueError(f"level {names[i]} has no γ-ray transitions in the scheme")

    def alpha(t) -> float:
        return float(t.conversion.value) if t.conversion is not None else 0.0

    if final is not None:
        f = index(final)
        chosen = [t for t in branches if t.final == f]
        if not chosen:
            raise ValueError(f"no transition from {names[i]} to {names[f]} in the scheme")
        tr = chosen[0]
    elif len(branches) == 1:
        tr = branches[0]
    else:
        if any(t.intensity is None for t in branches):
            raise ValueError(f"level {names[i]} has several branches without intensities: name the final level")
        tr = max(branches, key=lambda t: t.intensity.value * (1 + alpha(t)))
    branching, s_branch = 1.0, 0.0
    if len(branches) > 1:
        if any(t.intensity is None for t in branches):
            raise ValueError(f"the branching of {names[i]} needs the γ-ray intensities of all its transitions")
        tot = [t.intensity.value * (1 + alpha(t)) for t in branches]
        total = sum(tot)
        mine = tr.intensity.value * (1 + alpha(tr))
        branching = mine / total
        var = 0.0
        for t, part in zip(branches, tot):
            unc = t.intensity.unc
            if unc:
                d_part = unc * (1 + alpha(t))
                deriv = ((total - mine) / total ** 2) if t is tr else (-mine / total ** 2)
                var += (deriv * d_part) ** 2
        s_branch = math.sqrt(var)
    mult = (tr.multipolarity or "").replace("(", "").replace(")", "").upper()
    mixing = None
    if lam is None:
        if mult in ("", "E2"):
            lam = "E2"
        elif "+" in mult:
            lam = "E2"
            if tr.mixing is None:
                raise ValueError(f"the {mult} transition needs a mixing ratio, or name lam for a pure multipole")
            mixing = float(tr.mixing.value)
        else:
            lam = mult
    elif "+" in mult and tr.mixing is not None:
        mixing = float(tr.mixing.value)
    tau = float(result.tau_ps)
    b = float(lt.b_from_tau(tr.energy.value, lam, tau, branching, alpha(tr), mixing))
    rel = math.hypot(result.total_ps / tau, s_branch / branching)
    return b, b * rel


# ---------------------------------------------------------------------------------------------------------------
# Deorientation in vacuum


def _legendre(k: int, x) -> np.ndarray:
    c = np.zeros(k + 1)
    c[k] = 1.0
    return np.polynomial.legendre.legval(np.asarray(x, dtype=float), c)


def _w(theta_deg, a_k: dict, g: dict):
    x = np.cos(np.radians(np.asarray(theta_deg, dtype=float)))
    out = 1.0
    for k, a in a_k.items():
        out = out + float(a) * np.asarray(g.get(k, 1.0), dtype=float) * _legendre(int(k), x)
    return out


def deorientation_correction(theta_deg, a_k: dict, g_unshifted: dict, g_shifted: dict):
    """W_U/W_S: the angular distribution of the unshifted γ rays over that of the shifted ones at ``theta_deg``,
    W(θ) = 1 + Σ_k a_k G_k P_k(cos θ), with ``a_k`` the unattenuated coefficients {k: a_k} (A_k/A_0) and the
    attenuations ``g_unshifted`` = G_k(t) at the flight time to the stopper, ``g_shifted`` = G̃_k averaged over
    the decays in flight (:func:`lifetimes.deorientation`, :func:`lifetimes.deorientation_shifted`; arrays per
    distance broadcast). The measured unshifted intensity divided by this factor is what it would be with the
    shifted component's alignment, so the corrected share is R/(R + W_U/W_S (1 − R)) for a measured share R."""
    w_u = _w(theta_deg, a_k, g_unshifted)
    w_s = _w(theta_deg, a_k, g_shifted)
    return w_u / w_s


def _deorientation_inputs(component, theta, n) -> Optional[tuple]:
    if component is None:
        return None
    if theta is None:
        if not isinstance(component, dict):
            raise ValueError("without theta_deg, give each component as {k: (values, errors)} per distance")
        return {int(k): (_arr(v[0], n, f"a_{k}"), _arr(v[1], n, f"a_{k} errors")) for k, v in component.items()}
    values, errors = (np.asarray(x, dtype=float) for x in component)
    if values.shape != (n, len(theta)) or errors.shape != values.shape:
        raise ValueError("with theta_deg, give each component as (values, errors), each (distances, rings)")
    return values, errors


def fit_deorientation(distances_um, angular_unshifted, angular_shifted, a_k0: dict, tau_ps: float, beta: float, *,
                      theta_deg=None, proper: bool = True, offset_um: float = 0.0, guess: Optional[dict] = None
                      ) -> tuple:
    """α_k and Γ_k of the static deorientation model, G_k(t) = α_k + (1 − α_k) e^(−Γ_k t) (Radeck et al., PRC 85,
    014301 (2012)), from the angular distributions of the unshifted and the shifted γ rays against distance.

    The unshifted component carries G_k at the flight time t = (d + ``offset_um``)/(βc) (alignment frozen in the
    stopper), the shifted one G̃_k averaged over the decays in flight (Radeck Eq. 7), for the level's ``tau_ps``;
    ``a_k0`` = {k: a_k} the unattenuated coefficients. The angular distributions are either

    - ``theta_deg`` given: intensities per ring, (values, errors) each (distances, rings), any normalisation per
      distance (fitted): N_d (1 + Σ a_k G_k P_k(cos θ)); or
    - ``theta_deg`` None: the attenuated coefficients a_k G_k fitted per distance, {k: (values, errors)}.

    Either component may be None. Returns ({k: α_k}, {k: Γ_k in 1/ps}, covariance) with the covariance over
    (α_k for k sorted, then Γ_k for k sorted)."""
    d = np.asarray(distances_um, dtype=float).ravel()
    n = len(d)
    ks = sorted(int(k) for k in a_k0)
    if not ks:
        raise ValueError("give the unattenuated coefficients a_k0 = {k: a_k}")
    theta = None if theta_deg is None else np.asarray(theta_deg, dtype=float).ravel()
    unsh = _deorientation_inputs(angular_unshifted, theta, n)
    shif = _deorientation_inputs(angular_shifted, theta, n)
    if unsh is None and shif is None:
        raise ValueError("give the unshifted or the shifted angular distributions (or both)")
    t = _flight(d, beta, offset_um, 1.0, proper)
    guess = guess or {}
    p0 = [float(guess.get(("alpha", k), 0.3)) for k in ks] + [float(guess.get(("gamma", k), 0.2)) for k in ks]

    def attenuations(p):
        alpha = {k: p[j] for j, k in enumerate(ks)}
        gamma = {k: p[len(ks) + j] for j, k in enumerate(ks)}
        return lt.deorientation(t, alpha, gamma), lt.deorientation_shifted(t, tau_ps, alpha, gamma)

    def resid(p):
        g_u, g_s = attenuations(p)
        parts = []
        for comp, g in ((unsh, g_u), (shif, g_s)):
            if comp is None:
                continue
            if theta is None:
                for k in ks:
                    if k in comp:
                        v, e = comp[k]
                        parts.append((v - float(a_k0[k]) * np.asarray(g[k])) / e)
            else:
                v, e = comp
                shape = _w(theta[None, :], a_k0, {k: np.asarray(g[k])[:, None] for k in ks})  # (n, rings)
                w = 1.0 / e ** 2
                norm = np.sum(w * v * shape, axis=1) / np.sum(w * shape * shape, axis=1)  # profiled per distance
                parts.append(((v - norm[:, None] * shape) / e).ravel())
        return np.concatenate(parts)

    lower = [0.0] * len(ks) + [0.0] * len(ks)
    upper = [1.0] * len(ks) + [np.inf] * len(ks)
    p, cov, chi2 = _least_squares(resid, p0, lower, upper)
    ndf = len(resid(p)) - len(p)
    if ndf > 0 and chi2 / ndf > 1:
        cov = cov * (chi2 / ndf)
    alpha = {k: float(p[j]) for j, k in enumerate(ks)}
    gamma = {k: float(p[len(ks) + j]) for j, k in enumerate(ks)}
    return alpha, gamma, cov


__all__ = ["BUDGET_KEYS", "Feeding", "STOPPING_UNCERTAINTY", "TauResult", "b_from_result", "ddcm", "ddcm_three_foil",
           "decay_curve", "deorientation_correction", "fit_deorientation", "fit_lineshape", "peak_areas"]
