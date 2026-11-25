"""Deorientation in vacuum: the hyperfine interaction between an excited nucleus that has left the target and the
electrons of its highly charged ion, which turns the nucleus's spin before it decays and so evens out its γ rays.

The statistical tensors of a level that has been in vacuum for a time t are G_k(t) times what the collision (or the
decay that fed it) left. The model is the two-state model as GOSIA has it, with GOSIA's default parameters, so that a
GOSIA run compares like with like::

    from physim.nuclear.deorientation import TwoState

    level = TwoState("194Pt", spin=2)
    level.average(tau_ps=60.2, beta=0.05)          # {2: G_2, 4: G_4} averaged over the decay, as GOSIA
    level.g(t_ps=[0.0, 10.0, 100.0], beta=0.05)    # {2: array, 4: array} after times t in vacuum

Both give {k: value} as :func:`physim.nuclear.lifetimes.deorientation` does, arrays shaped like the times and speeds
together. Times and mean lives are the nucleus's own, in ps; ``beta`` is v/c of the recoil in vacuum.

Sources: F. Bosch and H. Spehl, Z. Phys. A 280 (1977) 329; R. Brenn et al., Z. Phys. A 281 (1977) 219; the GOSIA
manual (2012), section 4.4; and GOSIA's source, gosia_20110524.14.f, subroutines GKVAC, GKK and XSTATIC and function
ATS, whose arithmetic is followed here. Equation numbers are the manual's.

**The model.** The hyperfine field of the ion is *fluctuating* (the electron shells decay after the ion leaves the
target, reoriented at random every τ_c) until, at the rate Λ*, it becomes *static* (the ion's ground state):

    G_k(t) = e^{−Λ*t} G_k^fl(t) + ∫₀ᵗ Λ* e^{−Λ*t′} G_k^fl(t′) G_k^st(t − t′) dt′                          (4.19)
    G_k^st(t) = ⟨α_k⟩ + (1 − ⟨α_k⟩) e^{−Γt}                                                                (4.24)

with the hard core α_k(I, J) = Σ_F (2F + 1)²/(2J + 1) {F F k; I I J}² (:func:`hard_core`), in the static part averaged
over the ion's charge states (:func:`charge_states`), each with the ground-state spin J of the atom of that many
electrons (:func:`atomic_spin`). The fluctuating part has one mean atomic spin J and the rate

    λ_k = (1 − α_k)/τ_c · [1 − exp(−w_k/(1 − α_k))],   w_k = ⟨ω²⟩τ_c² = k(k + 1) J(J + 1)/3 · (ω τ_c)²,
    ω = μ_N g H/(ħJ),   H = K Z (v/c)^x                                                             (4.22, 4.28, 4.29)

and, from 1 with no slope at t = 0 to (1 − λ_kτ_c) e^{−λ_k t} at long times (4.27),

    G_k^fl(t) = e^{−λ_k t} [(1 − λ_kτ_c) + (λ_k/p + 2λ_kτ_c) e^{−pt} − (λ_k/p + λ_kτ_c) e^{−2pt}],
    p = {[9λ_k² + 8λ_kτ_c(w_k − λ_k²)]^{1/2} − 3λ_k} / (4λ_kτ_c).                                           (4.33)

Its average over the decay, (1/τ)∫₀^∞ G_k(t) e^{−t/τ} dt, is GOSIA's closed form (4.30–4.32), :meth:`TwoState.average`.

**GOSIA's, kept as it is.** In p GOSIA's code has w_k, a number, where the manual has ⟨ω²⟩ in ps⁻² (which would make
the curvature of G_k^fl at t = 0 equal to −⟨ω²⟩). The width of the charge-state distribution is half the manual's.
Charge states beyond three widths, and the neutral atom, are left out, so ⟨α_k⟩ has small steps in β. Γ, Λ* and τ_c do
not depend on g, Z or β. The spins of :func:`atomic_spin` are GOSIA's table.

**Where this departs from GOSIA.**

- g is Z/A with the mass number; GOSIA takes the mass typed in EXPT (about 0.3 % of w_k).
- g = 0 (or K = 0) is no interaction: G_k = 1. GOSIA's formula gives 0/0.
- GOSIA works out ⟨α_k⟩ once in a run, for the first experiment, and keeps it for the others (its flag IBYP); here
  every speed has its own (0.002 in G_2 of a ¹¹⁰Cd level at β = 0.037 after an experiment at 0.064).
- GOSIA prints a G_k for every level and k = 2, 4, 6. A level has tensors up to k = 2I only, and only those are given
  here (none for spin 0), unless ``ranks`` asks for more.
- GOSIA's recoil speed is non-relativistic and one per experiment (that of the second level, its NCM); here it is the
  caller's, one per nucleus if wanted.

As in GOSIA, Λ* = 0 switches deorientation off, a level that decays at once (mean life 0) keeps G_k = 1, and only the
time in vacuum counts: GOSIA starts at the reaction (thin target) and switches it off for a thick one (ITTE). Fields in
solids and the transient field of a nucleus crossing a magnetised foil are not included.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, Optional, Union

import numpy as np

from . import angular, data

#: μ_N/ħ in rad ps⁻¹ T⁻¹: GOSIA's 4789 rad s⁻¹ G⁻¹.
MU_N_PER_HBAR = 4789.0 * 1e4 * 1e-12
#: How many nuclei are worked on at once (the charge states of each are held in memory).
_CHUNK = 50_000


@dataclass(frozen=True)
class Parameters:
    """The parameters of the two-state model, with GOSIA's defaults. The numbers are those of GOSIA's ``VAC,`` input:
    1 J, 2 Γ, 3 Λ*, 4 τ_c, 5 g, 6 K (in 10⁸ T: 6e-6), 7 x."""

    #: J, the mean atomic spin of the fluctuating state (a whole or half-odd number above 0).
    atomic_spin: float = 3.0
    #: Γ, the width of the Larmor frequencies in the static state, 1/ps.
    width_per_ps: float = 0.02
    #: Λ*, the rate of going from the fluctuating to the static state, 1/ps. 0 switches deorientation off, as in GOSIA.
    rate_per_ps: float = 0.0345
    #: τ_c, the correlation time of the fluctuating field, ps. From 1 ps: below, GOSIA's p can be negative.
    correlation_ps: float = 3.5
    #: K of the field H = K Z (v/c)^x, tesla.
    field_tesla: float = 600.0
    #: x of the same.
    field_power: float = 0.6
    #: The nuclear g factor; None for Z/A, as GOSIA. Only its size matters.
    g_factor: Optional[float] = None

    def __post_init__(self):
        two_j = 2 * self.atomic_spin
        if not (two_j > 0 and abs(two_j - round(two_j)) < 1e-9):
            raise ValueError("the atomic spin J must be a whole or half-odd number above 0")
        if not (self.width_per_ps >= 0 and self.rate_per_ps >= 0 and self.field_tesla >= 0):
            raise ValueError("Γ, Λ* and K cannot be negative")
        if not (self.correlation_ps >= 1.0 and math.isfinite(self.correlation_ps)):
            raise ValueError("the correlation time τ_c must be 1 ps or more")
        if not math.isfinite(self.field_power) or (self.g_factor is not None and not math.isfinite(self.g_factor)):
            raise ValueError("x and g must be numbers")


DEFAULTS = Parameters()

# The atomic ground-state spin by number of electrons (GOSIA's function ATS): {J: numbers of electrons}; every other
# number, and more than 96, gives 0.
_ATS = {
    0.5: (1, 3, 5, 11, 13, 19, 29, 31, 37, 47, 49, 55, 61, 69, 73, 79, 81, 87),
    1.5: (7, 9, 15, 17, 21, 33, 35, 39, 51, 53, 71, 83, 85, 89),
    2.5: (23, 27, 41, 45, 63, 77),
    3.5: (57, 59, 95),
    4.5: (25, 43, 75),
    5.5: (91, 93),
    6.5: (67,),
    7.5: (65,),
    1.0: (20, 72),
    2.0: (8, 16, 34, 52, 56, 60, 64, 84, 90, 96),
    3.0: (24, 42, 74),
    4.0: (26, 44, 58, 62, 68, 76, 94),
    6.0: (92,),
    8.0: (66,),
}
_ATS_BY_N = {n: j for j, ns in _ATS.items() for n in ns}


def atomic_spin(electrons: int) -> float:
    """The ground-state spin J of an ion with this many electrons, as GOSIA takes it: the neutral atom's of as many
    electrons, mostly; 0 for a bare ion and for more than 96 electrons."""
    return _ATS_BY_N.get(int(electrons), 0.0)


@lru_cache(maxsize=None)
def hard_core(spin: float, j: float, k: int) -> float:
    """α_k(I, J) = Σ_F (2F + 1)²/(2J + 1) {F F k; I I J}²: what a static hyperfine coupling to an atomic spin J leaves
    of the rank-k orientation of a nucleus of spin I, averaged over time (4.23). 1 for J = 0; 0 for k above 2I."""
    total, f = 0.0, abs(spin - j)
    while f <= spin + j + 1e-9:
        total += ((2 * f + 1) * angular.wigner_6j(f, f, k, spin, spin, j)) ** 2 / (2 * j + 1)
        f += 1.0
    return float(total)


def charge_states(z: int, beta) -> tuple:
    """(charge states Q, their shares) of ions of atomic number ``z`` at speeds ``beta`` (above 0): (n × states), the
    shares 0 outside each ion's own range. Gaussian about Q₀ = Z h^0.6, h = 1/(1 + (0.012008 Z^0.45/β)^{5/3})
    (Nikolaev and Dmitriev), of width √(Q₀(1 − h))/2, cut at three widths and at 1 ≤ Q ≤ Z, as GOSIA's XSTATIC. An ion
    too slow to be charged by that rule gets Q = 1."""
    b = np.atleast_1d(np.asarray(beta, dtype=float))
    h = 1.0 / (1.0 + (z**0.45 * 0.012008 / b) ** 1.666667)
    q0 = z * h**0.6
    width = np.sqrt(q0 * (1.0 - h)) / 2.0
    up = np.maximum(np.minimum(np.trunc(q0 + 3.0 * width + 0.5), z), 1.0)
    down = np.maximum(np.trunc(q0 - 3.0 * width - 0.5), 1.0)
    q = down[:, None] + np.arange(int(np.max(up - down, initial=0.0)) + 1)
    # Relative to the likeliest state, so that nothing underflows; a vanishing width leaves the state nearest Q₀.
    x = ((q - q0[:, None]) / np.maximum(width, 1e-9)[:, None]) ** 2 / 2.0
    inside = q <= up[:, None]
    x = np.where(inside, x, np.inf)
    share = np.where(inside, np.exp(-(x - x.min(axis=1, keepdims=True))), 0.0)
    return np.where(inside, q, 0.0).astype(int), share / share.sum(axis=1, keepdims=True)


def _out(x: np.ndarray, shape: tuple):
    return float(x[0]) if shape == () else x.reshape(shape)


class TwoState:
    """The two-state model of a level of spin ``spin`` of ``nuclide`` ("194Pt") with ``parameters`` (GOSIA's defaults
    when None). ``ranks``: the k to work out; None for those of 2, 4, 6 up to 2I, which is what the γ rays' angular
    distribution needs; "all" for 1 … 2I, odd ranks by the same formulas, which is what a density matrix needs (one
    whose other ranks are left as they were is no density matrix any more: eigenvalues to −0.4). With all ranks it
    stays positive up to spin 2 (to 1e-6); above, the model itself, whose ranks each have their own rate, makes some
    pure states negative by up to 0.005."""

    def __init__(self, nuclide, spin: float, parameters: Optional[Parameters] = None,
                 ranks: Union[None, str, Iterable[int]] = None):
        nuc = data.nuclide(nuclide)
        self.z, self.a, self.spin = int(nuc.Z), int(nuc.A), float(spin)
        two_i = 2 * self.spin
        if not (two_i >= 0 and abs(two_i - round(two_i)) < 1e-9):
            raise ValueError("the level's spin must be a whole or half-odd number")
        self.parameters = p = parameters or DEFAULTS
        self.g_factor = float(p.g_factor) if p.g_factor is not None else self.z / self.a
        top = int(round(two_i))
        if ranks is None:
            ranks = [k for k in (2, 4, 6) if k <= top]
        elif isinstance(ranks, str):
            if ranks != "all":
                raise ValueError("ranks are numbers, or 'all'")
            ranks = range(1, top + 1)
        self.ranks = tuple(int(k) for k in ranks)
        if any(k < 1 for k in self.ranks):
            raise ValueError("ranks start at 1: rank 0 is the population, which nothing attenuates")
        #: α_k in the fluctuating state, by rank.
        self._fluctuating = np.array([hard_core(self.spin, p.atomic_spin, k) for k in self.ranks])
        #: α_k in the static state by charge state (0 … Z) and rank.
        self._by_charge = np.array([[hard_core(self.spin, atomic_spin(self.z - q), k) for k in self.ranks]
                                    for q in range(self.z + 1)]).reshape(self.z + 1, len(self.ranks))

    # -- what the model is made of -----------------------------------------------------------------------------------

    @property
    def off(self) -> bool:
        """Whether nothing is attenuated at any speed: Λ* = 0 (GOSIA's switch), or no field or moment."""
        return self.parameters.rate_per_ps < 1e-9 or self.g_factor == 0 or self.parameters.field_tesla == 0

    def alpha_static(self, beta) -> np.ndarray:
        """⟨α_k⟩ of the static state at each speed, averaged over the charge states: (speeds × ranks)."""
        b = np.atleast_1d(np.asarray(beta, dtype=float))
        out = np.ones((len(b), len(self.ranks)))
        for at in range(0, len(b), _CHUNK):
            part = slice(at, at + _CHUNK)
            moving = np.flatnonzero(b[part] > 0)
            if len(moving):
                q, share = charge_states(self.z, b[part][moving])
                out[at + moving] = np.einsum("nq,nqk->nk", share, self._by_charge[q])
        return out

    def rates(self, beta) -> dict:
        """The model's numbers at each speed, each (speeds × ranks): "alpha" (⟨α_k⟩ of the static state), "w"
        (⟨ω²⟩τ_c²), "lam" (λ_k, 1/ps) and "p" (1/ps)."""
        par = self.parameters
        b = np.atleast_1d(np.asarray(beta, dtype=float))
        j, tc = par.atomic_spin, par.correlation_ps
        field = par.field_tesla * self.z * np.power(b, par.field_power, where=b > 0, out=np.zeros_like(b))
        omega = MU_N_PER_HBAR * self.g_factor * field / j                         # rad/ps
        k = np.array(self.ranks, dtype=float)[None, :]
        w = ((omega * tc) ** 2 * j * (j + 1) / 3.0)[:, None] * k * (k + 1)
        free = np.maximum(1.0 - self._fluctuating, 0.0)[None, :]
        live = (free > 1e-12) & (w > 0)
        lam = np.where(live, free * -np.expm1(-w / np.where(free > 1e-12, free, 1.0)) / tc, 0.0)
        one = np.where(live, lam, 1.0)
        root = np.sqrt(9 * one**2 + 8 * one * tc * np.where(live, w - one**2, 0.0))
        # p = (root − 3λ)/(4λτ_c), the difference taken without cancellation; positive, as w ≥ λ² for τ_c ≥ 1 ps.
        p = np.where(live, 2 * (w - one**2) / (root + 3 * one), 1.0)
        return {"alpha": self.alpha_static(b), "w": w, "lam": lam, "p": p}

    # -- attenuation -------------------------------------------------------------------------------------------------

    def _work(self, times, beta, check: str, form) -> dict:
        """{k: ``form``(times, rates)} over the times and speeds broadcast together, 1 where nothing attenuates."""
        t, b = np.broadcast_arrays(np.asarray(times, dtype=float), np.asarray(beta, dtype=float))
        shape, t, b = t.shape, t.ravel(), b.ravel()
        if not np.all(t >= 0):
            raise ValueError(f"{check} cannot be negative")
        if not np.all((b >= 0) & (b < 1)):
            raise ValueError("speeds are v/c, from 0 to 1")
        out = np.ones((len(t), len(self.ranks)))
        if not self.off and self.ranks:
            rows = np.flatnonzero((t > 0) & (b > 0))
            for at in range(0, len(rows), _CHUNK):
                part = rows[at:at + _CHUNK]
                out[part] = form(t[part][:, None], self.rates(b[part]))
        return {k: _out(out[:, i], shape) for i, k in enumerate(self.ranks)}

    def _coefficients(self, r: dict) -> list:
        """G_k^fl(t) e^{−Λ*t} = Σ c e^{−μt}: the three pairs (c, μ)."""
        lam, p, tc = r["lam"], r["p"], self.parameters.correlation_ps
        ratio = lam / p
        return [(c, self.parameters.rate_per_ps + lam + n * p)
                for n, c in enumerate((1.0 - lam * tc, ratio + 2 * lam * tc, -(ratio + lam * tc)))]

    def g(self, t_ps, beta) -> dict:
        """{k: G_k(t)} of nuclei that have been in vacuum for ``t_ps`` (ps, their own time) at speed ``beta``: numbers
        or arrays that broadcast together. 1 at t = 0; at t = inf what the static field leaves."""
        par = self.parameters
        big, gam = par.rate_per_ps, par.width_per_ps

        def form(t, r):
            t = np.minimum(t, 1e15)
            alpha, total = r["alpha"], 0.0
            for c, mu in self._coefficients(r):
                # (e^{−Γt} − e^{−μt})/(μ − Γ), also where μ is close to Γ.
                gap = np.abs(mu - gam)
                with np.errstate(divide="ignore", invalid="ignore"):
                    spread = np.where(gap * t > 1e-12, -np.expm1(-gap * t) / gap, t) * np.exp(-np.minimum(mu, gam) * t)
                e = np.exp(-mu * t)
                total = total + c * (e + big * alpha * -np.expm1(-mu * t) / mu + big * (1.0 - alpha) * spread)
            return total

        return self._work(t_ps, beta, "times", form)

    def average(self, tau_ps, beta) -> dict:
        """{k: G_k} averaged over the decay of a level of mean life ``tau_ps`` (ps) that is in vacuum from the start,
        at speed ``beta``: GOSIA's coefficients (4.30–4.32). 1 for a mean life of 0; :meth:`g` at infinity for inf."""
        par = self.parameters
        big, gam, tc = par.rate_per_ps, par.width_per_ps, par.correlation_ps

        def form(tau, r):
            lam, p, alpha = r["lam"], r["p"], r["alpha"]
            x = 1.0 / tau
            # GOSIA's expressions divided through by τ (and τ²), which also holds for τ = inf.
            with np.errstate(divide="ignore", invalid="ignore"):
                static = np.where(x + gam > 0, (x + alpha * gam) / np.where(x + gam > 0, x + gam, 1.0), 1.0)
            full = big + lam + x
            return (x + big * static) / full * (1.0 + lam * (full - 2 * p * p * tc) / ((full + p) * (full + 2 * p)))

        return self._work(tau_ps, beta, "mean lives", form)


@lru_cache(maxsize=256)
def model(nuclide: str, spin: float, parameters: Optional[Parameters] = None, ranks=None) -> TwoState:
    """The :class:`TwoState` of a level, made once (``ranks`` as a tuple, "all" or None)."""
    return TwoState(nuclide, float(spin), parameters, ranks)


__all__ = ["DEFAULTS", "MU_N_PER_HBAR", "Parameters", "TwoState", "atomic_spin", "charge_states", "hard_core", "model"]
