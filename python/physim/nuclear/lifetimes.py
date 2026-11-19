"""Lifetime physics for plunger (RDDS) and Doppler-shift attenuation (DSAM) measurements: pure functions on numbers.

Nothing here reads a setup: the planner, the event generators and the analysis call these with numbers they took
from the setup (``planning/lifetimes.md``, wave 1B; the signatures are the frozen interfaces of wave 0)::

    from physim.nuclear import lifetimes as lt

    tau = lt.tau_from_b(652.0, "E2", 1300.0)               # ps, from B(E2↓) in e² fm⁴
    t = lt.flight_times([10, 50, 200], beta=0.045)         # lab ps at each distance (µm)
    lt.unshifted([tau], None, [1.0], t)                    # share of the γ rays emitted after the stopper
    lt.recommend_distances(tau, 0.045)                     # µm, log-spaced over the sensitive range

Units: time ps, distance µm, γ-ray energy keV, speed as β = v/c, B(Eλ) in e² fm^{2λ}, B(Mλ) in µ_N² fm^{2λ−2}.

**Frames of time.** A lifetime τ is the proper (rest-frame) mean life. :func:`flight_times` gives lab times, as the
frozen contract says; with ``proper=True`` it divides by γ = 1/√(1 − β²) (:func:`lorentz_gamma`), which is what the
decay law needs. :func:`bateman`, :func:`unshifted`, :func:`deorientation` and :func:`deorientation_shifted` take
times in the frame of the τ they are given, so the caller chooses: proper flight times with proper τ (time dilation
included), or lab times (dilation left out; it is below 0.2% for β < 6%).

**A cascade** is given by its levels' mean lives ``taus_ps``, a feeding matrix ``feeding`` (levels × levels, row =
from, column = to: the share of the row level's decays, γ rays and conversion together, that go to the column level;
rows sum to at most 1, the rest leaves the listed levels) as in ``residues.ResidueCascade``, and ``populations``, the
number of nuclei put directly into each level at t = 0. A side-feeding cascade with its own time is a level of its
own in the same matrix. ``feeding=None`` means no feeding.

**The DSAM kernel** (:class:`Kernel`) is what ``transport.kernel`` returns: at a grid of lab times after the state
was formed, β·cos θ of a sample of recoils, θ measured from the beam axis, the same recoil in each column. It does
not depend on τ, so a line-shape fit re-weights it for each trial τ (AhKin's idea).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Protocol, runtime_checkable

import numpy as np

from .gamma import doppler_energy, level_names
from .levels import ELECTRIC, HBARC_MEV_FM, LN2, fractions, rate_per_b

#: Speed of light in µm/ps (= mm/ns).
C_UM_PER_PS = 299.792458
#: Proton mass in MeV (CODATA 2018), for the nuclear magneton.
PROTON_MASS_MEV = 938.27208816
#: µ_N / e = ħc / (2 m_p c²), in fm: B(Mλ) in µ_N² fm^{2λ−2} is B in e² fm^{2λ} divided by this squared.
NUCLEAR_MAGNETON_FM = HBARC_MEV_FM / (2.0 * PROTON_MASS_MEV)


def _out(x):
    """A 0-d array as a float, anything else unchanged."""
    a = np.asarray(x)
    return float(a) if a.ndim == 0 else a


# ---------------------------------------------------------------------------------------------------------------
# B(Eλ) and τ


def _multipole(lam: str) -> tuple:
    text = str(lam).strip().upper()
    if len(text) != 2 or text[0] not in "EM" or text[1] not in "123456":
        raise ValueError(f"write the multipole as E1, E2, M1, ...; got '{lam}'")
    return text[0], int(text[1])


def rate_per_unit_b(lam: str, energy_kev: float) -> float:
    """γ-ray emission rate per unit B(λ↓), in 1/s: :func:`levels.rate_per_b` for Eλ (per e² fm^{2λ}); for Mλ the
    same times (µ_N/e)², per µ_N² fm^{2λ−2}."""
    kind, order = _multipole(lam)
    rate = rate_per_b(order, energy_kev)
    return rate * NUCLEAR_MAGNETON_FM**2 if kind == "M" else rate


def _share(lam: str, mixing: Optional[float]) -> float:
    """The multipole's share of the transition's γ-ray intensity. ``mixing`` is the E2/M1 mixing ratio δ."""
    if mixing is None:
        return 1.0
    kind, order = _multipole(lam)
    name = f"{kind}{order}"
    if name not in ("M1", "E2"):
        raise ValueError(f"a mixing ratio here is the E2/M1 ratio δ, for a transition given as M1 or E2, not {name}")
    return fractions("M1+E2", float(mixing))[name]


def _check_b_args(energy_kev, branching, conversion) -> None:
    if np.any(np.asarray(energy_kev, dtype=float) <= 0):
        raise ValueError("the γ-ray energy must be positive")
    b = np.asarray(branching, dtype=float)
    if np.any(b <= 0) or np.any(b > 1 + 1e-12):
        raise ValueError("the branching (the transition's share of the level's decays) must be in (0, 1]")
    if np.any(np.asarray(conversion, dtype=float) < 0):
        raise ValueError("the conversion coefficient α cannot be negative")


def tau_from_b(energy_kev, lam: str, b_down, branching=1.0, conversion=0.0, mixing=None) -> float:
    """The level's mean life (ps) from B(λ↓) of one of its transitions, the inverse of ``LevelScheme._b_down``.

    ``branching`` is the transition's share of all the level's decays, γ rays and conversion electrons together:
    I_γ (1 + α) / Σ_k I_γ,k (1 + α_k). ``conversion`` is the transition's total conversion coefficient α.
    ``mixing`` is the E2/M1 mixing ratio δ (``lam`` then names the multipole B belongs to, M1 or E2): the multipole
    carries the share 1/(1 + δ²) (M1) or δ²/(1 + δ²) (E2) of the γ rays. So

        1/τ = T(λ) B(λ↓) (1 + α) / (share × branching),   T(λ) = :func:`rate_per_unit_b`.
    """
    _check_b_args(energy_kev, branching, conversion)
    b_down = np.asarray(b_down, dtype=float)
    if np.any(b_down <= 0):
        raise ValueError("B(λ↓) must be positive")
    share = _share(lam, mixing)
    if share <= 0:
        raise ValueError(f"with this mixing ratio the transition carries no {lam}: B({lam}) gives no lifetime")
    total_rate = rate_per_unit_b(lam, energy_kev) * b_down / share * (1.0 + np.asarray(conversion)) / branching
    return _out(1e12 / total_rate)


def b_from_tau(energy_kev, lam: str, tau_ps, branching=1.0, conversion=0.0, mixing=None) -> float:
    """B(λ↓) (e² fm^{2λ}, or µ_N² fm^{2λ−2} for Mλ) from the level's mean life in ps; the inverse of
    :func:`tau_from_b`, with the same meaning of ``branching``, ``conversion`` and ``mixing``."""
    _check_b_args(energy_kev, branching, conversion)
    tau_ps = np.asarray(tau_ps, dtype=float)
    if np.any(tau_ps <= 0):
        raise ValueError("the mean life must be positive")
    gamma_rate = 1e12 / tau_ps * branching / (1.0 + np.asarray(conversion))
    return _out(gamma_rate * _share(lam, mixing) / rate_per_unit_b(lam, energy_kev))


def level_taus(scheme) -> dict:
    """{level name: (τ in ps, source)} for each level of a :class:`~physim.nuclear.levels.LevelScheme` that has γ
    rays and a known τ. Names are :func:`physim.nuclear.gamma.level_names` ("2+₁", "4+").

    τ is the level's half-life / ln 2 when it has one (source ``"half-life (ensdf)"``, ``"half-life (user)"``, ...).
    Otherwise it comes from the B(Eλ↓) of its transitions (source ``"B(E2) (user)"``, ``"B(E2) (shell model)"``,
    ...: the matrix element's own source), through :func:`tau_from_b`: the strongest transition with a matrix element
    and a known multipole share gives the total rate, with the branching from the γ-ray intensities and conversion.
    When the intensities are not all known and there are several branches, the rates of the transitions with a
    matrix element are added, which makes τ an upper limit (the source then ends in ", other branches unknown").
    Levels with neither are left out.
    """
    names = level_names(scheme.levels)
    out: dict = {}
    for i, level in enumerate(scheme.levels):
        branches = [t for t in scheme.transitions if t.initial == i and t.energy is not None and t.energy.value > 0]
        if not branches:
            continue
        hl = level.half_life
        if hl is not None and hl.value > 0:
            out[names[i]] = (hl.value / LN2 * 1e12, f"half-life ({hl.source})")
            continue
        found = _taus_from_matrix_elements(scheme, i, branches)
        if found is not None:
            out[names[i]] = found
    return out


def _taus_from_matrix_elements(scheme, i: int, branches: list) -> Optional[tuple]:
    def alpha(t) -> float:
        return float(t.conversion.value) if t.conversion is not None else 0.0

    known = []  # (transition, total decay rate through it in 1/s, multipole, source of its matrix element)
    for t in branches:
        for mult in ELECTRIC:
            b = scheme.b(i, t.final, mult)
            if b is None or b <= 0:
                continue
            if t.multipolarity:
                share = fractions(t.multipolarity, None if t.mixing is None else t.mixing.value).get(mult)
            else:
                share = 1.0  # a typed transition without a multipolarity: the matrix element says what it is
            if not share:
                continue
            rate = rate_per_unit_b(mult, t.energy.value) * b / share * (1.0 + alpha(t))
            known.append((t, rate, mult, scheme.matrix_element(i, t.final, mult).source))
            break
    if not known:
        return None
    if len(branches) == 1:
        t, rate, mult, source = known[0]
        return 1e12 / rate, f"B({mult}) ({source})"
    if all(t.intensity is not None and t.intensity.value >= 0 for t in branches):
        total = sum(t.intensity.value * (1.0 + alpha(t)) for t in branches)
        usable = [k for k in known if k[0].intensity.value > 0]
        if total > 0 and usable:
            t, rate, mult, source = max(usable, key=lambda k: k[0].intensity.value * (1.0 + alpha(k[0])))
            branching = t.intensity.value * (1.0 + alpha(t)) / total
            return branching * 1e12 / rate, f"B({mult}) ({source})"
    rate = sum(k[1] for k in known)
    return 1e12 / rate, f"B({known[0][2]}) ({known[0][3]}), other branches unknown"


# ---------------------------------------------------------------------------------------------------------------
# Decay curves


def _cascade(taus_ps, feeding, populations) -> tuple:
    """(decay constants 1/ps, feeding matrix, populations) checked."""
    taus = np.atleast_1d(np.asarray(taus_ps, dtype=float))
    if taus.ndim != 1:
        raise ValueError("give one mean life per level")
    n = len(taus)
    if np.any(~(taus > 0)):
        raise ValueError("mean lives must be positive (infinity for a level that does not decay)")
    lam = np.where(np.isinf(taus), 0.0, 1.0 / taus)
    f = np.zeros((n, n)) if feeding is None else np.asarray(feeding, dtype=float)
    if f.shape != (n, n):
        raise ValueError(f"the feeding matrix must be {n} × {n} (levels × levels, row = from, column = to)")
    if np.any(f < 0) or np.any(np.diag(f) != 0) or np.any(f.sum(axis=1) > 1 + 1e-9):
        raise ValueError("feeding shares must be non-negative, zero on the diagonal, and sum to at most 1 per row")
    p = np.atleast_1d(np.asarray(populations, dtype=float))
    if p.shape != (n,) or np.any(p < 0):
        raise ValueError(f"give {n} non-negative populations, one per level")
    return lam, f, p


def _expm_metzler(a: np.ndarray, t: np.ndarray, integral: bool = False) -> tuple:
    """exp(a t) for each t ≥ 0 (and ∫₀ᵗ exp(a s) ds when ``integral``), for a matrix whose off-diagonal entries are
    non-negative, as the rate matrix of a cascade is.

    Uniformisation: with c the largest decay constant, a + c·1 has no negative entry, so exp(a h) =
    e^{−ch} Σ_k ((a + c)h)^k / k! is a sum of non-negative terms, and so is each squaring that takes h back up to t.
    Nothing cancels, so every entry keeps its relative accuracy, even far down the tail and for τ₁ ≈ τ₂ where the
    closed form of Bateman's equations subtracts nearly equal numbers. The integral comes from the same exponential
    of the block matrix [[a, σ1], [0, 0]] (Van Loan).
    """
    n = len(a)
    c = max(0.0, float(-np.min(np.diag(a))))
    if integral:
        sigma = c if c > 0 else 1.0
        g = np.zeros((2 * n, 2 * n))
        g[:n, :n] = a
        g[:n, n:] = sigma * np.eye(n)
        bound = c + sigma
    else:
        g, bound = a, c
    size = len(g)
    t = np.asarray(t, dtype=float).ravel()
    if np.any(~(t >= 0)):
        raise ValueError("times must be non-negative")
    if bound == 0:
        e = np.broadcast_to(np.eye(size), (len(t), size, size)).copy()
    else:
        b = g + c * np.eye(size)
        s = np.ceil(np.log2(np.maximum(bound * t, 1e-300) / 0.5)).clip(0, None).astype(int)
        h = t / 2.0**s
        x = b[None] * h[:, None, None]  # column sums ≤ 0.5
        term = np.broadcast_to(np.eye(size), x.shape).copy()
        e = term.copy()
        for k in range(1, 20):
            term = term @ x / k
            e += term
        e *= np.exp(-c * h)[:, None, None]
        for j in range(int(s.max(initial=0))):
            idx = np.flatnonzero(s > j)
            e[idx] = e[idx] @ e[idx]
    if integral:
        return e[:, :n, :n], e[:, :n, n:] / sigma
    return e, None


def _rate_matrix(lam: np.ndarray, f: np.ndarray) -> np.ndarray:
    """dN/dt = a N: a[j, i] = feeding[i, j] λ_i off the diagonal, a[i, i] = −λ_i."""
    return f.T * lam[None, :] - np.diag(lam)


def _passing(f: np.ndarray, x: np.ndarray) -> np.ndarray:
    """The nuclei that will pass through each level from a state x (levels on the last axis): (1 − fᵀ)⁻¹ x."""
    m = np.eye(len(f)) - f.T
    try:
        inv = np.linalg.inv(m)
    except np.linalg.LinAlgError:
        raise ValueError("the feeding matrix has a closed loop: decays must lead out of the levels") from None
    if not np.all(np.isfinite(inv)):
        raise ValueError("the feeding matrix has a closed loop: decays must lead out of the levels")
    return x @ inv.T


def bateman(taus_ps, feeding, populations, t_ps) -> np.ndarray:
    """N_i(t): the nuclei still in each level at time t (Bateman's equations for the cascade; module docstring),
    shape (levels, len(t)), or (levels,) for one time. Times in the frame of τ."""
    lam, f, p = _cascade(taus_ps, feeding, populations)
    t = np.asarray(t_ps, dtype=float)
    if not f.any():
        if np.any(~(t >= 0)):
            raise ValueError("times must be non-negative")
        n = p[None, :] * np.exp(-lam[None, :] * t.ravel()[:, None])
    else:
        e, _ = _expm_metzler(_rate_matrix(lam, f), t)
        n = e @ p  # (times, levels)
    return n.T.reshape((len(lam),) + t.shape)


def _decays_after(lam, f, p, t, integral: bool = False) -> tuple:
    """Per level, the share of its decays after each time (times, levels), and with ``integral`` its integral
    over [0, t] (ps)."""
    if not f.any():  # no feeding: each level on its own, in closed form
        tt = np.asarray(t, dtype=float).ravel()[:, None]
        if np.any(~(tt >= 0)):
            raise ValueError("times must be non-negative")
        reached = np.where(p > 0, 1.0, np.nan)[None, :]
        after = np.exp(-lam[None, :] * tt) * reached
        return after, (_mean_exp(lam[None, :], tt) * reached if integral else None)
    e, s = _expm_metzler(_rate_matrix(lam, f), t, integral)
    total = _passing(f, p)
    with np.errstate(invalid="ignore", divide="ignore"):
        norm = np.where(total > 0, 1.0 / np.where(total > 0, total, 1.0), np.nan)
    after = _passing(f, e @ p) * norm
    return after, (None if s is None else _passing(f, s @ p) * norm)


def unshifted(taus_ps, feeding, populations, flight_time_ps) -> np.ndarray:
    """Per level, the share of its decays (γ rays) that come after ``flight_time_ps``: in a plunger, the unshifted
    (stopped) component, I_u / (I_u + I_s). Shape (levels, len(t)), or (levels,) for one time; ``nan`` for a
    level no nucleus reaches. Times in the frame of τ (``flight_times(..., proper=True)`` for proper τ).

    For one level, exp(−t/τ); with feeding, (1 − fᵀ)⁻¹ N(t) over (1 − fᵀ)⁻¹ N(0): the nuclei that are in the level
    or will still pass through it, over all that ever do."""
    lam, f, p = _cascade(taus_ps, feeding, populations)
    t = np.asarray(flight_time_ps, dtype=float)
    after, _ = _decays_after(lam, f, p, t)
    return after.T.reshape((len(lam),) + t.shape)


def lorentz_gamma(beta) -> float:
    """γ = 1/√(1 − β²)."""
    b = np.asarray(beta, dtype=float)
    if np.any(b < 0) or np.any(b >= 1):
        raise ValueError("β must be in [0, 1)")
    return _out(1.0 / np.sqrt(1.0 - b * b))


def flight_times(distances_um, beta, cos_recoil=1.0, *, proper: bool = False) -> np.ndarray:
    """Time of flight (ps) over each target-to-stopper distance (µm) at speed β, the recoil at angle θ to the beam
    axis (``cos_recoil`` = cos θ; the distance is along the axis): d / (β c cos θ). Lab time, as the frozen
    contract says; ``proper=True`` divides by γ (time dilation), giving the time in the nucleus's frame, which is
    what :func:`unshifted` needs with a proper τ. Arrays broadcast (one β per event, say)."""
    d = np.asarray(distances_um, dtype=float)
    b = np.asarray(beta, dtype=float)
    cos = np.asarray(cos_recoil, dtype=float)
    if np.any(b <= 0) or np.any(b >= 1):
        raise ValueError("β must be in (0, 1)")
    if np.any(cos <= 0) or np.any(cos > 1):
        raise ValueError("the recoil must move towards the stopper: cos θ in (0, 1]")
    if np.any(d < 0):
        raise ValueError("distances cannot be negative")
    t = d / (b * C_UM_PER_PS * cos)
    if proper:
        t = t / lorentz_gamma(b)
    return _out(t)


# ---------------------------------------------------------------------------------------------------------------
# Planning a plunger measurement


def _precision(x: float) -> float:
    """Relative uncertainty of τ per √(counts) from the unshifted share R = e^{−x} at x = t/τ, one level: R has the
    binomial error √(R(1 − R)/N), and τ = t / (−ln R), so σ_τ/τ = √((1 − R)/(R N)) / x = √(e^x − 1) / (x √N)."""
    return math.sqrt(math.expm1(x)) / x


@lru_cache(maxsize=None)
def _best_x() -> float:
    """Where :func:`_precision` is smallest (x ≈ 1.59): golden-section search."""
    lo, hi = 1e-6, 50.0
    for _ in range(200):
        a, b = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
        if _precision(a) < _precision(b):
            hi = b
        else:
            lo = a
    return 0.5 * (lo + hi)


@lru_cache(maxsize=None)
def _sensitive_x(factor: float) -> tuple:
    """(x_min, x_max): where :func:`_precision` is within ``factor`` of its best (at x ≈ 1.59)."""
    def solve(lo, hi, target, rising):
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if (_precision(mid) > target) == rising:
                hi = mid
            else:
                lo = mid
        return math.sqrt(lo * hi)

    best = _best_x()
    target = factor * _precision(best)
    return solve(1e-9, best, target, rising=False), solve(best, 200.0, target, rising=True)


def _decay_length_um(tau_ps: float, beta: float) -> float:
    """γ β c τ: the mean distance (µm) flown before the decay, time dilation included."""
    if not tau_ps > 0:
        raise ValueError("the mean life must be positive")
    if not 0 < beta < 1:
        raise ValueError("β must be in (0, 1)")
    return lorentz_gamma(beta) * beta * C_UM_PER_PS * tau_ps


def sensitive_range(tau_ps: float, beta: float, *, factor: float = 2.0) -> tuple:
    """(d_min, d_max) in µm: the distances at which a measured unshifted share gives τ to within ``factor`` of the
    best relative precision for the same number of counts (one level, recoil along the axis). With the default
    factor 2 (a distance costs at most four times the beam time of the best one for the same σ_τ), this is
    0.18 to 5.1 decay lengths γβcτ, the best at 1.59; the region where dR/dd is large, and roughly what experimenters
    choose (⁹⁸Ru, Radeck et al. 2012: 3–450 µm for vτ = 113 µm)."""
    if not factor > 1:
        raise ValueError("the factor must be above 1")
    x_lo, x_hi = _sensitive_x(float(factor))
    length = _decay_length_um(float(tau_ps), float(beta))
    return x_lo * length, x_hi * length


def _round_sig(x: float, digits: int = 2) -> float:
    if x <= 0:
        return 0.0
    scale = 10.0 ** (math.floor(math.log10(x)) - digits + 1)
    return round(x / scale) * scale


def recommend_distances(tau_ps: float, beta: float, n: int = 8) -> list:
    """``n`` target-to-stopper distances (µm), log-spaced over :func:`sensitive_range` (about 0.18 to 5.1 decay
    lengths γβcτ, around βcτ), rounded to two significant figures; fewer if rounding merges some. One distance is
    the best single one, 1.59 decay lengths."""
    if n < 1:
        raise ValueError("ask for at least one distance")
    d_lo, d_hi = sensitive_range(tau_ps, beta)
    if n == 1:
        return [_round_sig(_best_x() * _decay_length_um(tau_ps, beta))]
    out = sorted({_round_sig(d) for d in np.geomspace(d_lo, d_hi, n)})
    return [float(d) for d in out]


def peak_separation(e0_kev: float, beta, theta_deg) -> float:
    """Energy of the shifted (in-flight) peak minus the unshifted one (keV) for a detector at ``theta_deg`` to the
    recoil direction: E₀ √(1 − β²)/(1 − β cos θ) − E₀ (relativistic, as :func:`physim.nuclear.gamma.doppler_energy`).
    Negative at backward angles."""
    return _out(doppler_energy(e0_kev, beta, np.cos(np.radians(theta_deg))) - e0_kev)


# ---------------------------------------------------------------------------------------------------------------
# DSAM


@runtime_checkable
class Kernel(Protocol):
    """The slowing-down kernel of a DSAM measurement, independent of τ (``transport.kernel`` returns one).

    ``times_ps``: (nt,) increasing lab times since the state was formed, from 0, out to where the recoils have
    stopped; after the last time each recoil keeps its last value. ``beta_cos``: (nt, n) β·cos θ of n recoils at each
    time, θ from the beam axis, the same recoil in each column (0 once it has stopped); the recoils are taken as
    moving along the beam. Or, for recoils with their own directions (``transport.Kernel``), a method
    ``beta_cos(theta_deg, phi_deg)`` giving (nt, n) β·cos(angle between the recoil's velocity and a detector at
    (θ, φ)), with an array ``beta`` (nt, n) of their speeds: their sideways motion then counts. An optional attribute
    ``weights`` (n,) weights the recoils (equal weights without it)."""

    times_ps: np.ndarray
    beta_cos: np.ndarray


@dataclass(frozen=True)
class SampledKernel:
    """A :class:`Kernel` from arrays."""

    times_ps: np.ndarray
    beta_cos: np.ndarray
    weights: Optional[np.ndarray] = None


def exponential_kernel(beta0: float, alpha_ps: float, times_ps=None) -> SampledKernel:
    """The kernel of stopping proportional to speed, every recoil along the beam: β(t) = β₀ e^{−t/α}. For it the
    attenuation is exactly F(τ) = α/(α + τ). Default times: 0 to 40α in steps of α/1000."""
    if times_ps is None:
        times_ps = np.linspace(0.0, 40.0 * alpha_ps, 40001)
    t = np.asarray(times_ps, dtype=float)
    return SampledKernel(t, (beta0 * np.exp(-t / alpha_ps))[:, None])


def _directional(kernel) -> bool:
    """Whether the kernel gives each recoil's own direction (a method ``beta_cos(theta, phi)`` and ``beta``)."""
    return callable(kernel.beta_cos)


def _kernel_arrays(kernel, theta_deg: float = 0.0) -> tuple:
    """(times (nt,), speeds β (nt, n), β·cos(angle to a detector at ``theta_deg`` from the beam) (nt, n), weights).
    For an array kernel the recoils move along the beam: β is |β cos θ_recoil| and the projection is
    β cos θ_recoil · cos θ_detector. A directional kernel gives the projection itself and the full speed."""
    t = np.asarray(kernel.times_ps, dtype=float)
    if _directional(kernel):
        proj = np.asarray(kernel.beta_cos(float(theta_deg), 0.0), dtype=float)
        beta = np.asarray(kernel.beta, dtype=float)
        if proj.ndim == 1:
            proj, beta = proj[:, None], beta[:, None]
        if beta.shape != proj.shape:
            raise ValueError("the kernel's beta and beta_cos must have the same shape (nt, recoils)")
    else:
        beta = np.asarray(kernel.beta_cos, dtype=float)
        if beta.ndim == 1:
            beta = beta[:, None]
        proj = beta * math.cos(math.radians(float(theta_deg)))
    if t.ndim != 1 or len(t) < 1 or beta.shape[0] != len(t):
        raise ValueError("the kernel needs times_ps (nt,) and beta_cos (nt, recoils)")
    if np.any(np.diff(t) <= 0) or t[0] < 0:
        raise ValueError("the kernel's times must start at or after 0 and increase")
    w = getattr(kernel, "weights", None)
    w = np.ones(beta.shape[1]) if w is None else np.asarray(w, dtype=float)
    if w.shape != (beta.shape[1],) or np.any(w < 0) or w.sum() <= 0:
        raise ValueError("the kernel's weights must be one non-negative number per recoil")
    return t, beta, proj, w / w.sum()


def _emission_weights(lam, f, p, t) -> tuple:
    """Per level, the share of its decays in each interval of the kernel's times, (levels, nt − 1), and the node
    weights (levels, nt) that integrate a quantity varying linearly between the times exactly against the decay
    law (decays before the first time go to the first node, after the last time to the last node)."""
    after, integral = _decays_after(lam, f, p, t, integral=True)
    after, integral = after.T, integral.T  # (levels, nt)
    h = np.diff(t)
    shares = after[:, :-1] - after[:, 1:]
    late = (np.diff(integral, axis=1) - h * after[:, 1:]) / h  # ∫ (t − t_k) p dt / h over each interval
    late = np.clip(late, 0.0, shares)
    nodes = np.zeros_like(after)
    nodes[:, :-1] += shares - late
    nodes[:, 1:] += late
    nodes[:, 0] += 1.0 - after[:, 0]
    nodes[:, -1] += after[:, -1]
    return shares, nodes, after


def dsam_attenuation(tau_ps: float, kernel) -> float:
    """F(τ): the mean β·cos θ of the recoils at the moment of emission over its value at formation, for one level
    with mean life τ (in the frame of the kernel's lab times; time dilation, below 0.2% for β < 6%, is left out).
    β·cos θ is taken linear between the kernel's times and integrated exactly against the decay law."""
    t, _, bc, w = _kernel_arrays(kernel, 0.0)  # β·cos to the beam axis
    lam, f, p = _cascade([tau_ps], None, [1.0])
    mean = bc @ w
    if mean[0] == 0:
        raise ValueError("the kernel's recoils do not move at formation")
    _, nodes, _ = _emission_weights(lam, f, p, t)
    return float(nodes[0] @ mean / mean[0])


def _spread(lo: np.ndarray, hi: np.ndarray, weight: np.ndarray, edges: np.ndarray, group=None,
            n_groups: int = 1) -> np.ndarray:
    """Bin contents (n_groups, bins) of segments [lo, hi], each carrying ``weight`` spread evenly along it (a
    zero-length segment is a point) and counted in its ``group`` (all in group 0 without). Only positive numbers
    are added, so the result keeps its relative accuracy however many short segments there are; the parts outside
    the edges are dropped."""
    nb = len(edges) - 1
    size = n_groups * nb
    group = np.zeros(len(lo), dtype=np.int64) if group is None else np.asarray(group, dtype=np.int64)
    jlo = np.searchsorted(edges, lo, side="right") - 1
    jhi = np.searchsorted(edges, hi, side="right") - 1
    length = hi - lo
    one = (jlo == jhi) | (length <= 1e-13)  # E/E₀ is about 1: shorter than that is a point
    inside = one & (jlo >= 0) & (jlo < nb)
    out = np.bincount(group[inside] * nb + jlo[inside], weights=weight[inside], minlength=size).astype(float)
    many = np.flatnonzero(~one & (jhi >= 0) & (jlo < nb))
    if len(many):
        lo, length, weight, jlo, span, group = (lo[many], length[many], weight[many], jlo[many], (jhi - jlo)[many],
                                                group[many])
        for o in range(int(span.max()) + 1):
            j = jlo + o
            m = (o <= span) & (j >= 0) & (j < nb)
            if not m.any():
                continue
            jj = j[m]
            left = np.clip((edges[jj] - lo[m]) / length[m], 0.0, 1.0)
            right = np.clip((edges[jj + 1] - lo[m]) / length[m], 0.0, 1.0)
            out += np.bincount(group[m] * nb + jj, weights=weight[m] * (right - left), minlength=size)
    return out.reshape(n_groups, nb)


#: The largest (intervals × bins) matrix :func:`lineshape` keeps; beyond it each level is binned on its own.
_MATRIX_LIMIT = 4_000_000
_MATRICES: dict = {}


def _ratio(beta: np.ndarray, proj: np.ndarray) -> np.ndarray:
    """E/E₀ = √(1 − β²)/(1 − β cos ψ), ψ the angle between the recoil's velocity and the γ ray."""
    return np.sqrt(1.0 - beta * beta) / (1.0 - proj)


def _kernel_matrix(kernel, t, beta, proj, w, theta_deg: float, x_edges: np.ndarray) -> tuple:
    """The kernel binned once, independent of τ (AhKin's idea): (intervals, bins) the share of each interval's
    emissions in each bin of E/E₀, and the bins of the recoils as they are at the first and the last time. The last
    few are kept, so a fit that calls :func:`lineshape` for many τ bins the kernel once (kernels are taken as
    unchanging). Keyed on the detector's cos θ for an array kernel, and on θ for a directional one (whose
    projection is not a function of cos θ alone)."""
    if _directional(kernel):
        key = (id(kernel), id(kernel.times_ps), id(kernel.beta), "theta", float(theta_deg), x_edges.tobytes())
    else:
        key = (id(kernel), id(kernel.times_ps), id(kernel.beta_cos), "cos", math.cos(math.radians(float(theta_deg))),
               x_edges.tobytes())
    hit = _MATRICES.get(key)
    if hit is not None and hit[0] is kernel:
        return hit[1]
    ratio = _ratio(beta, proj)  # (nt, recoils): E/E₀
    lo = np.minimum(ratio[:-1], ratio[1:]).ravel()
    hi = np.maximum(ratio[:-1], ratio[1:]).ravel()
    groups = np.repeat(np.arange(len(t) - 1), beta.shape[1])
    weights = np.tile(w, len(t) - 1)
    matrix = _spread(lo, hi, weights, x_edges, groups, len(t) - 1)
    first = _spread(ratio[0], ratio[0], w, x_edges)[0]
    last = _spread(ratio[-1], ratio[-1], w, x_edges)[0]
    if len(_MATRICES) >= 4:
        _MATRICES.pop(next(iter(_MATRICES)))
    _MATRICES[key] = (kernel, (matrix, first, last))
    return matrix, first, last


def lineshape(e0_kev, taus_ps, feeding, populations, kernel, theta_deg, edges_kev) -> np.ndarray:
    """The DSAM line shape of each level's γ ray at a detector at ``theta_deg`` to the beam axis: (levels, bins),
    each row the share of that level's γ rays in each energy bin (a row sums to 1 when the bins cover the line).
    ``e0_kev`` is one energy, or one per level. No detector response: fold it in afterwards.

    The emission times follow the cascade (module docstring) and the kernel gives the recoils' velocities at each
    time; the γ ray has E₀ √(1 − β²)/(1 − β cos ψ), ψ the angle between the recoil's velocity and the detector. An
    array kernel's recoils move along the beam (β cos ψ = β cos θ_recoil · cos θ; their sideways motion is left out,
    which matters only near 90°); a directional kernel (``transport.Kernel``) gives β cos ψ for the detector at
    (θ, φ = 0) with each recoil's own direction. Between two kernel times each recoil's energy is spread evenly between its values at the two
    times, with the decays of that interval. The kernel's binning does not depend on τ and is kept for the next
    call, so a fit over τ is fast."""
    t, beta, proj, w = _kernel_arrays(kernel, theta_deg)
    lam, f, p = _cascade(taus_ps, feeding, populations)
    edges = np.asarray(edges_kev, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("the bin edges must increase")
    e0 = np.broadcast_to(np.asarray(e0_kev, dtype=float), lam.shape)
    if np.any(e0 <= 0):
        raise ValueError("the γ-ray energy must be positive")
    shares, _, after = _emission_weights(lam, f, p, t)
    out = np.full((len(lam), len(edges) - 1), np.nan)
    for i in range(len(lam)):
        if not np.isfinite(after[i, 0]):
            continue
        x = edges / e0[i]
        early, late = 1.0 - after[i, 0], after[i, -1]  # before the first time, after the last
        if (len(t) - 1) * (len(edges) - 1) <= _MATRIX_LIMIT:
            matrix, first, last = _kernel_matrix(kernel, t, beta, proj, w, float(theta_deg), x)
            out[i] = shares[i] @ matrix + early * first + late * last
            continue
        ratio = _ratio(beta, proj)
        lo = np.minimum(ratio[:-1], ratio[1:]).ravel()
        hi = np.maximum(ratio[:-1], ratio[1:]).ravel()
        out[i] = _spread(lo, hi, (shares[i][:, None] * w[None, :]).ravel(), x)[0]
        out[i] += _spread(ratio[0], ratio[0], early * w, x)[0] + _spread(ratio[-1], ratio[-1], late * w, x)[0]
    return out


# ---------------------------------------------------------------------------------------------------------------
# Deorientation in vacuum


def _deorientation_args(alpha: dict, gamma_per_ps: dict) -> list:
    if set(alpha) != set(gamma_per_ps):
        raise ValueError("give α_k and Γ_k for the same k")
    out = []
    for k in sorted(alpha):
        a, g = float(alpha[k]), float(gamma_per_ps[k])
        if not 0 <= a <= 1:
            raise ValueError(f"α_{k} must be in [0, 1]")
        if g < 0:
            raise ValueError(f"Γ_{k} cannot be negative")
        out.append((k, a, g))
    return out


def deorientation(t_ps, alpha: dict, gamma_per_ps: dict) -> dict:
    """G_k(t) = α_k + (1 − α_k) e^{−Γ_k t}: the attenuation of the k-th order of the γ-ray angular distribution
    after a time t in vacuum (static model of the hyperfine interaction, Radeck et al., PRC 85, 014301 (2012)).
    ``alpha`` and ``gamma_per_ps`` are {k: value}, usually k = 2, 4. Returns {k: array shaped like t}."""
    t = np.asarray(t_ps, dtype=float)
    if np.any(t < 0):
        raise ValueError("times must be non-negative")
    return {k: _out(a + (1.0 - a) * np.exp(-g * t)) for k, a, g in _deorientation_args(alpha, gamma_per_ps)}


def _mean_exp(rate, t):
    """(1 − e^{−rate·t}) / rate, and t where rate·t is 0, without cancellation."""
    rate = np.asarray(rate, dtype=float)
    x = rate * t
    with np.errstate(invalid="ignore", divide="ignore"):
        value = -np.expm1(-x) / np.where(rate > 0, rate, 1.0)
    return np.where(x > 0, value, t)


def deorientation_shifted(flight_time_ps, tau_ps, alpha: dict, gamma_per_ps: dict) -> dict:
    """G̃_k(T): G_k averaged over the decays in flight, between 0 and the flight time T (Radeck et al. Eq. 7):
    ∫₀ᵀ G_k(t) e^{−t/τ} dt / ∫₀ᵀ e^{−t/τ} dt = α_k + (1 − α_k) [(1 − e^{−(Γ_k + 1/τ)T})/(Γ_k + 1/τ)] /
    [τ (1 − e^{−T/τ})], which is 1 at T = 0. This is the shifted component's attenuation; the unshifted one is
    G_k(T) (:func:`deorientation`), the alignment frozen once the nucleus stops. T and τ in the same frame."""
    t = np.asarray(flight_time_ps, dtype=float)
    if np.any(t < 0):
        raise ValueError("flight times must be non-negative")
    if not tau_ps > 0:
        raise ValueError("the mean life must be positive")
    lam = 0.0 if math.isinf(tau_ps) else 1.0 / tau_ps
    den = _mean_exp(lam, t)
    out = {}
    for k, a, g in _deorientation_args(alpha, gamma_per_ps):
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(t > 0, _mean_exp(g + lam, t) / np.where(t > 0, den, 1.0), 1.0)
        out[k] = _out(a + (1.0 - a) * ratio)
    return out


__all__ = ["C_UM_PER_PS", "Kernel", "SampledKernel", "b_from_tau", "bateman", "deorientation",
           "deorientation_shifted", "dsam_attenuation", "exponential_kernel", "flight_times", "level_taus",
           "lineshape", "lorentz_gamma", "peak_separation", "rate_per_unit_b", "recommend_distances",
           "sensitive_range", "tau_from_b", "unshifted"]
