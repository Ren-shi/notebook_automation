"""Deorientation in vacuum: the hyperfine interaction between an excited nucleus that has left the target and the
electrons of its highly charged, highly excited ion, which attenuates the nucleus's alignment before it decays.

The statistical tensors of a level that has spent a time t in vacuum become G_k(t) · t_kq, so its γ rays come out
more evenly than the collision left them. The model is the two-state model of Brenn and Spehl, in the form and
with the default parameters of GOSIA, so that a GOSIA run compares like with like:

- F. Bosch and H. Spehl, Z. Phys. A 280 (1977) 329; R. Brenn, H. Spehl, A. Weckherlin, H. A. Doubt and
  G. van Middelkoop, Z. Phys. A 281 (1977) 219;
- the GOSIA manual (T. Czosnyka, D. Cline, C. Y. Wu et al., version of 10 May 2012), section 4.4 and table 4.1,
  and the VAC switch of the CONT suboption (section 7.3);
- the GOSIA source, subroutines GKVAC, GKK and XSTATIC and the function ATS (gosia_20110524.14.f), whose
  arithmetic this module follows where the manual leaves something open.

**The model.** The hyperfine field of the ion is either *fluctuating* (the electron shells decaying after the ion
leaves the target, reoriented at random every τ_c) or *static* (the ion's ground state); the ion goes from the
first to the second at the rate Λ*. Then

    G_k(t) = e^{−Λ*t} G_k^fl(t) + ∫₀ᵗ Λ* G_k^fl(t′) e^{−Λ*t′} G_k^st(t − t′) dt′          (manual 4.19)

with the static part

    G_k^st(t) = ⟨α_k⟩ + (1 − ⟨α_k⟩) e^{−Γt},    α_k(I, J) = Σ_F (2F + 1)² / (2J + 1) {F F k; I I J}²   (4.23, 4.24)

averaged over the charge states Q of the ion, whose electronic ground-state spin J is that of the atom with Z − Q
electrons (:func:`atomic_spin`); the charge states are Gaussian about Q₀ = Z h^0.6 with
h = 1 / (1 + (0.012008 Z^0.45 / β)^{5/3}) (Nikolaev and Dmitriev, 4.25–4.26). The fluctuating part has the rate

    λ_k = (1 − ⟨α_k⟩_J) / τ_c · (1 − exp(−⟨ω²⟩τ_c² / (1 − ⟨α_k⟩_J))),
    ⟨ω²⟩ = k(k + 1) J(J + 1) ω₀² / 3,   ω₀ = μ_N g H / (ħ J),   H = K Z (v/c)^x            (4.28, 4.29, 4.22)

for one mean atomic spin J (default 3). At long times G_k^fl = (1 − λ_kτ_c) e^{−λ_k t} (4.27); GOSIA joins it to
t = 0 by two more exponentials, which gives its time-integral coefficient in closed form (4.30–4.33):

    G_k(τ) = λ_I ∫₀^∞ G_k(t) e^{−λ_I t} dt = G^BS [1 + λ_kτ (r − 2p²ττ_c) / ((r + pτ)(r + 2pτ))],
    G^BS = [1 + Λ*τ (1 + ⟨α_k⟩Γτ) / (1 + Γτ)] / r,    r = (Λ* + λ_k)τ + 1,
    p = [(9λ_k² + 8λ_kτ_c(w − λ_k²))^{1/2} − 3λ_k] / (4λ_kτ_c).

τ is the level's mean life. The joining function this corresponds to is

    G_k^fl(t) = e^{−λ_k t} [(1 − λ_kτ_c) + (λ_k/p + 2λ_kτ_c) e^{−pt} − (λ_k/p + λ_kτ_c) e^{−2pt}],

which is 1 with zero slope at t = 0 and has the asymptote of 4.27; with it the convolution of 4.19 has a closed
form (:meth:`TwoState.g`), whose average over the decay time is GOSIA's G_k(τ) exactly (:meth:`TwoState.integral`).
In the manual w is ⟨ω²⟩, which makes the second derivative of G_k^fl at t = 0 equal to −⟨ω²⟩; GOSIA's code uses
w = ⟨ω²⟩τ_c² there (the same quantity as in λ_k), and so does this module, to give GOSIA's numbers.

The defaults (manual table 4.1, GOSIA's VAC switch): J = 3, Γ = 0.02 ps⁻¹, Λ* = 0.0345 ps⁻¹, τ_c = 3.5 ps,
g = Z/A, K = 600 T (6·10⁻⁶ in GOSIA's units), x = 0.6. Odd ranks, which no γ ray's distribution holds but a
density matrix may, follow the same formulas.

**Where it applies.** Only in vacuum: a nucleus decaying while it is still in the target or its backing, or after
it has stopped there, is not deoriented by this (the hyperfine fields of a solid, and the transient field of a
nucleus moving through a magnetised one, are not included). A level without a half-life decays at once and is not
deoriented. GOSIA applies the time integral from the moment of the reaction, as for a thin target, and switches
deorientation off for a thick target (its ``ITTE`` flag).

    from physim.nuclear.deorientation import TwoState

    model = TwoState("194Pt", spin=2.0)
    model.integral(tau_ps=60.2, beta=0.05)        # {k: G_k} averaged over the decay, as GOSIA
    model.g(t_ps=[0.0, 10.0, 100.0], beta=0.05)   # G_k(t) for a nucleus in vacuum for t
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Optional

import numpy as np

from . import angular, data

#: μ_N / ħ in rad s⁻¹ G⁻¹, as GOSIA writes it (4789): with the field in tesla and the time in picoseconds,
#: 4789 × 10⁴ G/T × 10⁻¹² s/ps.
MU_N_PER_HBAR = 4789.0 * 1e4 * 1e-12


@dataclass(frozen=True)
class Parameters:
    """The parameters of the two-state model, with GOSIA's defaults (manual table 4.1)."""

    #: Mean atomic spin J of the fluctuating state.
    atomic_spin: float = 3.0
    #: Γ, width of the Larmor-frequency distribution of the static state, ps⁻¹.
    width_per_ps: float = 0.02
    #: Λ*, rate of transition from the fluctuating to the static state, ps⁻¹.
    rate_per_ps: float = 0.0345
    #: τ_c, correlation time of the fluctuating field, ps.
    correlation_ps: float = 3.5
    #: K of the mean field H = K Z (v/c)^x, tesla.
    field_tesla: float = 600.0
    #: x of the same.
    field_power: float = 0.6
    #: The nuclear g factor; None for Z/A, as GOSIA.
    g_factor: Optional[float] = None


DEFAULTS = Parameters()

# The atomic ground-state spin by number of electrons, as GOSIA's function ATS has it (gosia_20110524.14.f):
# {J: numbers of electrons}; every other number up to 96 gives 0.
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
    """The electronic ground-state spin J of an ion with this many electrons, as GOSIA takes it (its function
    ATS): mostly the neutral atoms' ground terms; 0 for a bare ion and beyond 96 electrons."""
    n = int(electrons)
    if n <= 0 or n > 96:
        return 0.0
    return _ATS_BY_N.get(n, 0.0)


def charge_states(z: int, beta: float) -> tuple:
    """(charge states Q, their weights) of an ion of atomic number ``z`` leaving a foil at speed ``beta``: the
    Nikolaev–Dmitriev mean Q₀ = Z h^0.6 with h = 1 / (1 + (0.012008 Z^0.45 / β)^{5/3}), and a Gaussian of width
    √(Q₀(1 − h)) / 2 over the whole charges within three widths of Q₀ (1 ≤ Q ≤ Z), as GOSIA's XSTATIC. (The
    manual writes the width without the factor 1/2; the code has it.)"""
    beta = max(float(beta), 1e-4)
    h = 1.0 / (1.0 + (z**0.45 * 0.012008 / beta) ** 1.666667)
    q0 = z * h**0.6
    dq = max(math.sqrt(q0 * (1.0 - h)) / 2.0, 1e-6)
    up = min(math.trunc(q0 + 3.0 * dq + 0.5), z)
    down = max(math.trunc(q0 - 3.0 * dq - 0.5), 1)
    q = np.arange(down, up + 1)
    if not len(q):
        q = np.array([max(min(int(round(q0)), z), 1)])
    w = np.exp(-(((q0 - q) / dq) ** 2) / 2.0)
    return q, w / w.sum()


@lru_cache(maxsize=None)
def hard_core(spin: float, j: float, k: int) -> float:
    """α_k(I, J) = Σ_F (2F + 1)² / (2J + 1) {F F k; I I J}²: what of the rank-k orientation of a nucleus of spin
    I a static hyperfine coupling to an electronic spin J leaves, averaged over time (manual 4.23). 1 for J = 0
    or k = 0."""
    total = 0.0
    f = abs(spin - j)
    while f <= spin + j + 1e-9:
        total += ((2 * f + 1) * angular.wigner_6j(f, f, k, spin, spin, j)) ** 2 / (2 * j + 1)
        f += 1.0
    return float(total)


class TwoState:
    """The two-state model for a level of spin ``spin`` of ``nuclide``, with ``parameters`` (GOSIA's defaults if
    None). Ranks k = 0 … 2I; the attenuation of rank 0 is 1."""

    #: The grid of speeds on which the static part is tabulated (it changes slowly with β).
    BETA_GRID = np.linspace(0.002, 0.25, 497)

    def __init__(self, nuclide: str, spin: float, parameters: Optional[Parameters] = None):
        nuc = data.nuclide(nuclide)
        self.z, self.a, self.spin = int(nuc.Z), int(nuc.A), float(spin)
        self.parameters = parameters or DEFAULTS
        p = self.parameters
        self.g_factor = p.g_factor if p.g_factor is not None else self.z / self.a
        self.ranks = np.arange(0, int(round(2 * self.spin)) + 1)
        #: ⟨α_k⟩ of the fluctuating state (the one mean atomic spin), by rank.
        self.alpha_fluctuating = np.array([hard_core(self.spin, p.atomic_spin, int(k)) for k in self.ranks])
        self._static = None

    # -- the parts of the model ---------------------------------------------------------------------------------------

    def alpha_static(self, beta) -> np.ndarray:
        """⟨α_k⟩ of the static state, averaged over the charge states at each speed: (speeds × ranks)."""
        if self._static is None:
            rows = []
            for b in self.BETA_GRID:
                q, w = charge_states(self.z, float(b))
                rows.append([float(np.sum(w * [hard_core(self.spin, atomic_spin(self.z - int(x)), int(k))
                                               for x in q])) for k in self.ranks])
            self._static = np.array(rows)
        b = np.atleast_1d(np.asarray(beta, dtype=float))
        grid = self.BETA_GRID
        pos = np.clip((b - grid[0]) / (grid[1] - grid[0]), 0.0, len(grid) - 1.0)
        i = np.minimum(pos.astype(int), len(grid) - 2)
        f = (pos - i)[:, None]
        return (1 - f) * self._static[i] + f * self._static[i + 1]

    def rates(self, beta) -> dict:
        """The model's quantities at each speed, each (speeds × ranks): "alpha" (static ⟨α_k⟩), "w" (⟨ω²⟩τ_c²),
        "lam" (λ_k, ps⁻¹) and "p" (ps⁻¹)."""
        par = self.parameters
        b = np.atleast_1d(np.asarray(beta, dtype=float))
        j, tc = par.atomic_spin, par.correlation_ps
        field = par.field_tesla * self.z * np.maximum(b, 0.0) ** par.field_power
        omega = MU_N_PER_HBAR * self.g_factor * field / j                     # rad/ps
        k = self.ranks[None, :].astype(float)
        w = ((omega * tc) ** 2 * j * (j + 1) / 3.0)[:, None] * k * (k + 1)
        af = self.alpha_fluctuating[None, :]
        free = np.maximum(1.0 - af, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            lam = np.where(free > 1e-12, free * (1.0 - np.exp(-w / np.where(free > 1e-12, free, 1.0))) / tc, 0.0)
            root = np.sqrt(np.maximum(9 * lam**2 + 8 * lam * tc * (w - lam**2), 0.0))
            p = np.where(lam > 1e-12, (root - 3 * lam) / (4 * np.where(lam > 1e-12, lam, 1.0) * tc), 1.0)
        return {"alpha": self.alpha_static(b), "w": w, "lam": lam, "p": p}

    # -- attenuation ----------------------------------------------------------------------------------------------------

    def g(self, t_ps, beta) -> np.ndarray:
        """G_k(t) for nuclei in vacuum for ``t_ps`` picoseconds at speed ``beta`` (arrays of one length, or
        scalars): (nuclei × ranks), 1 at t = 0."""
        t = np.atleast_1d(np.asarray(t_ps, dtype=float))
        b = np.broadcast_to(np.atleast_1d(np.asarray(beta, dtype=float)), t.shape)
        par = self.parameters
        r = self.rates(b)
        lam, p, alpha = r["lam"], r["p"], r["alpha"]
        tc, big, gam = par.correlation_ps, par.rate_per_ps, par.width_per_ps
        live = lam > 1e-12
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(live, lam / np.where(live, p, 1.0), 0.0)
        coeffs = (1.0 - lam * tc, ratio + 2 * lam * tc, -(ratio + lam * tc))
        tt = np.maximum(t, 0.0)[:, None]
        out = np.zeros_like(lam)
        for j, c in enumerate(coeffs):
            mu = big + lam + j * p
            e_mu = np.exp(-mu * tt)
            static = alpha * (1.0 - e_mu) / mu if big > 0 else 0.0
            gap = mu - gam
            near = np.abs(gap) < 1e-9
            with np.errstate(divide="ignore", invalid="ignore"):
                spread = np.where(near, tt * np.exp(-gam * tt), (np.exp(-gam * tt) - e_mu) / np.where(near, 1.0, gap))
            out = out + c * (e_mu + big * static + big * (1.0 - alpha) * spread)
        out[:, 0] = 1.0
        return out

    def integral(self, tau_ps, beta) -> np.ndarray:
        """G_k averaged over the decay of a level of mean life ``tau_ps`` (ps) that is in vacuum from its
        formation, at speed ``beta``: GOSIA's time-integral coefficients (manual 4.30–4.33), (levels × ranks)."""
        tau = np.atleast_1d(np.asarray(tau_ps, dtype=float))
        b = np.broadcast_to(np.atleast_1d(np.asarray(beta, dtype=float)), tau.shape)
        par = self.parameters
        r = self.rates(b)
        lam, p, alpha = r["lam"], r["p"], r["alpha"]
        tc, big, gam = par.correlation_ps, par.rate_per_ps, par.width_per_ps
        t = tau[:, None]
        up = 1.0 + big * t * (gam * t * alpha + 1.0) / (t * gam + 1.0)
        down = t * (lam + big) + 1.0
        out = up / down * (1.0 + lam * t * (down - 2 * p * p * t * tc) / ((down + p * t) * (down + 2 * p * t)))
        out[:, 0] = 1.0
        return out

    def hard_limit(self, beta) -> np.ndarray:
        """G_k(t → ∞): what remains of the orientation once the static field has had all the time it wants,
        Λ* ⟨α_k⟩ ∫₀^∞ G_k^fl(t) e^{−Λ*t} dt (speeds × ranks)."""
        par = self.parameters
        r = self.rates(beta)
        lam, p, alpha = r["lam"], r["p"], r["alpha"]
        tc, big = par.correlation_ps, par.rate_per_ps
        live = lam > 1e-12
        ratio = np.where(live, lam / np.where(live, p, 1.0), 0.0)
        coeffs = (1.0 - lam * tc, ratio + 2 * lam * tc, -(ratio + lam * tc))
        out = sum(c * big * alpha / (big + lam + j * p) for j, c in enumerate(coeffs))
        out = np.asarray(out, dtype=float)
        out[:, 0] = 1.0
        return out


@lru_cache(maxsize=None)
def rank_projectors(spin: float) -> np.ndarray:
    """The projectors onto the statistical tensors of each rank, as maps on density matrices flattened row by
    row: P[k] is (2I + 1)² × (2I + 1)², and Σ_k P[k] is the identity. Built from the orthonormal tensor operators
    ⟨I m|T_kq|I m′⟩ = (−1)^{I − m} √(2k + 1) (I k I; −m q m′)."""
    ms = [-spin + i for i in range(int(round(2 * spin)) + 1)]
    d = len(ms)
    out = np.zeros((d, d * d, d * d), dtype=complex)
    for k in range(d):
        for q in range(-k, k + 1):
            t = np.zeros((d, d))
            for a, m in enumerate(ms):
                for b, mp in enumerate(ms):
                    if abs(-m + q + mp) < 1e-9:
                        t[a, b] = ((-1) ** int(round(spin - m)) * math.sqrt(2 * k + 1)
                                   * angular.wigner_3j(spin, k, spin, -m, q, mp))
            v = t.reshape(-1)
            out[k] += np.outer(v, v)
    return out


@lru_cache(maxsize=None)
def _tensor_basis(spin: float) -> tuple:
    """(U, ranks): the rows of U are the conjugated tensor operators T_kq flattened, so U ρ is the statistical
    tensors of ρ in the normalisation Tr(T_kq† ρ) and ρ = U† (U ρ); ``ranks`` is k of each row."""
    ms = [-spin + i for i in range(int(round(2 * spin)) + 1)]
    d = len(ms)
    rows, ranks = [], []
    for k in range(d):
        for q in range(-k, k + 1):
            t = np.zeros((d, d))
            for a, m in enumerate(ms):
                for b, mp in enumerate(ms):
                    if abs(-m + q + mp) < 1e-9:
                        t[a, b] = ((-1) ** int(round(spin - m)) * math.sqrt(2 * k + 1)
                                   * angular.wigner_3j(spin, k, spin, -m, q, mp))
            rows.append(t.reshape(-1))
            ranks.append(k)
    return np.array(rows, dtype=complex), np.array(ranks)


def attenuate(rho: np.ndarray, spin: float, factors: np.ndarray) -> np.ndarray:
    """Density matrices (n × (2I + 1)²) with the rank-k part of each multiplied by ``factors[:, k]`` (n × ranks,
    at least 2I + 1 columns): G_k t_kq for each. The trace (rank 0) is kept when ``factors[:, 0]`` is 1. Two
    matrix products, to the tensors and back."""
    rho = np.asarray(rho)
    n, d = rho.shape[0], rho.shape[1]
    if not n:
        return rho
    u, ranks = _tensor_basis(float(spin))
    t = rho.reshape(n, d * d) @ u.T
    t *= np.asarray(factors)[:, ranks]
    return (t @ u.conj()).reshape(n, d, d)


@lru_cache(maxsize=64)
def model(nuclide: str, spin: float, parameters: Optional[Parameters] = None) -> TwoState:
    """The :class:`TwoState` of a level, made once (its table of the static part is kept with it)."""
    return TwoState(nuclide, float(spin), parameters)


class Faded:
    """The states of a solution's table (:meth:`~physim.nuclear.orientation.Excitation.table`) with every level
    deoriented before it decays: at each grid angle, each level's density matrix (as the collision left it, plus
    what the levels above feed it) has its rank-k part multiplied by ``factors[level][angle, k]`` before it emits
    and before it feeds the levels below, as GOSIA treats the feeding (manual section 4.4). Populations are
    unchanged. :meth:`near` interpolates as :meth:`~physim.nuclear.orientation.Excitation.near` does."""

    def __init__(self, solution, factors: dict, step: float = 2.0):
        from .orientation import Populated, tensors_of

        self.solution, self.step = solution, step
        grid, states = solution.table(step)
        self.grid = grid
        levels = solution.scheme.levels
        self.states = []
        for i, st in enumerate(states):
            rho = _faded_cascade(solution, st.density, {n: f[i] for n, f in factors.items()})
            tensors = {n: tensors_of(m, float(levels[n].spin)) for n, m in rho.items()}
            self.states.append(Populated(solution, float(grid[i]), dict(st.direct), tensors, st.isotropic,
                                         dict(st.density)))

    def near(self, theta_cm: float):
        """The deoriented state at ``theta_cm``, interpolated linearly in the tensors between grid angles."""
        from .orientation import Populated

        grid, step = self.grid, self.step
        t = min(max(float(theta_cm), grid[0]), grid[-1])
        i = min(int((t - grid[0]) / step), len(grid) - 2)
        f = (t - grid[i]) / step
        lo, hi = self.states[i], self.states[i + 1]
        tensors = {n: {key: (1 - f) * v + f * hi.tensors[n][key] for key, v in lo.tensors[n].items()}
                   for n in lo.tensors}
        direct = {n: (1 - f) * v + f * hi.direct[n] for n, v in lo.direct.items()}
        return Populated(self.solution, float(theta_cm), direct, tensors, self.solution.isotropic)


def _faded_cascade(solution, direct: dict, factors: dict) -> dict:
    """:meth:`~physim.nuclear.orientation.Excitation._cascade` with each level's density matrix attenuated by
    ``factors[level]`` (one G_k per rank) once everything above has fed it, before it feeds the levels below."""
    from .correlation import _feeding_map

    levels = solution.scheme.levels
    decays = solution._decays
    rho = {n: m.copy() for n, m in direct.items()}
    order = sorted(rho.keys() | {d.final for d in decays.values()}, key=lambda n: -levels[n].energy.value)
    done = set()
    while order:
        n = order.pop(0)
        if n in done or n not in rho:
            done.add(n)
            continue
        done.add(n)
        if n in factors and levels[n].spin is not None:
            rho[n] = attenuate(rho[n][None], float(levels[n].spin), np.asarray(factors[n])[None])[0]
        for (i, f), d in decays.items():
            if i != n or levels[f].spin is None or not d.l1:
                continue
            ji, jf = float(levels[n].spin), float(levels[f].spin)
            share = [(d.l1, 1.0)] if not d.l2 else [(d.l1, 1 / (1 + d.mixing**2)),
                                                    (d.l2, d.mixing**2 / (1 + d.mixing**2))]
            parts = [w * np.einsum("abcd,cd->ab", _feeding_map(ji, jf, int(L)), rho[n]) for L, w in share
                     if angular.triangle(ji, L, jf)]
            if parts:
                fed = sum(parts)
                rho[f] = rho.get(f, np.zeros_like(fed)) + d.branch * fed
                if f not in done and f not in order:
                    order.append(f)
                    order.sort(key=lambda m: -levels[m].energy.value)
    return rho


__all__ = ["DEFAULTS", "Faded", "MU_N_PER_HBAR", "Parameters", "TwoState", "atomic_spin", "attenuate",
           "charge_states", "hard_core", "model", "rank_projectors"]
