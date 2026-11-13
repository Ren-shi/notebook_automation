"""The angular correlation between successive γ rays of one cascade (backlog item 71): a γ ray seen along a
direction leaves the level it feeds oriented by that direction, and the next γ ray follows that orientation.

Each excited nucleus is carried down its cascade as a density matrix ρ in the orbit's frame (see
:mod:`physim.nuclear.orientation`). A transition from spin I to I_f by multipoles L (and L′ with mixing ratio δ)
emits a photon along the unit vector n with the probability

    W(n) ∝ Σ_σ Tr[A_σ ρ_n A_σ^†],    (A_σ)_{m M} = Σ_L c_L(σ) √(2L + 1) ⟨I_f m L σ | I M⟩,

where ρ_n = D^† ρ D is ρ quantised along n (D the rotation taking z to n), σ = ±1 the photon's helicity along n,
c_L = 1 and c_L′(σ) = −σδ (the electric and the magnetic multipole of a mixture differ in sign between the two
helicities). The level fed is left with

    ρ_f(n) = D_f [Σ_σ A_σ ρ_n A_σ^†] D_f^† / Tr[...],

which is what the next transition starts from. Summed over every direction, ρ_f(n) W(n) gives back the unobserved
feeding of :func:`~physim.nuclear.orientation.fed_density`, and W(n) the distribution of
:meth:`~physim.nuclear.orientation.Populated.w`: the singles are those of a cascade without the correlation, and
only how the γ rays of one cascade share directions changes. For an unoriented 4⁺ → 2⁺ → 0⁺ cascade the two γ
rays follow W(θ) = 1 + 0.1020 P₂(cos θ) + 0.0091 P₄(cos θ).

The rotations are exp(−iβJ_y) by the eigenvectors of J_y, worked out once per spin, so that many events are
rotated at once::

    from physim.nuclear.correlation import Emission

    e = Emission(4.0, 2.0, 2)                 # a 4⁺ → 2⁺ E2 transition
    w = e.weights(rho, directions)            # 4π W along each direction (rho: n × 9 × 9, trace 1)
    rho_2 = e.after(rho, directions)          # the 2⁺ state each γ ray leaves (n × 5 × 5, trace 1)

Not included: the conversion electrons' own correlation (a converted transition feeds as an unobserved one), and
deorientation while a level lives.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from . import angular
from .orientation import _ms

#: Events rotated at once (the rotation matrices of a block are n × (2I + 1)² complex numbers).
BLOCK = 20_000


@lru_cache(maxsize=None)
def _jy_eigen(j: float) -> tuple:
    """Eigenvalues and eigenvectors of J_y for spin j, in the basis M ascending from −j."""
    ms = _ms(j)
    d = len(ms)
    jy = np.zeros((d, d), dtype=complex)
    for a, m in enumerate(ms):
        for b, mp in enumerate(ms):
            if abs(m - mp - 1) < 1e-9:
                jy[a, b] += math.sqrt(j * (j + 1) - mp * (mp + 1)) / 2j
            if abs(m - mp + 1) < 1e-9:
                jy[a, b] -= math.sqrt(j * (j + 1) - mp * (mp - 1)) / 2j
    lam, vec = np.linalg.eigh(jy)
    return lam, vec, np.array(ms)


def rotation(j: float, directions: np.ndarray) -> np.ndarray:
    """D^j_{M′M}(φ, θ, 0) for each unit vector (θ, φ) of ``directions`` (n × 3): n × (2j + 1) × (2j + 1), the
    rotation that takes the z axis to the direction."""
    d = np.asarray(directions, dtype=float)
    theta = np.arccos(np.clip(d[:, 2], -1.0, 1.0))
    phi = np.arctan2(d[:, 1], d[:, 0])
    lam, vec, ms = _jy_eigen(float(j))
    small = ((vec[None, :, :] * np.exp(-1j * np.outer(theta, lam))[:, None, :]) @ vec.conj().T).real
    return np.exp(-1j * np.outer(phi, ms))[:, :, None] * small


def isotropic(j: float) -> np.ndarray:
    """The density matrix of an unoriented level of spin j, trace 1."""
    d = int(round(2 * j + 1))
    return np.eye(d, dtype=complex) / d


class Emission:
    """One transition from spin ``ji`` to ``jf`` by the multipole ``l1`` (mixed with ``l2`` by the mixing ratio
    ``mixing``, the amplitude of ``l2`` over that of ``l1``), as a γ ray seen along a direction."""

    def __init__(self, ji: float, jf: float, l1: int, l2: int = 0, mixing: float = 0.0):
        self.ji, self.jf, self.l1, self.l2, self.mixing = float(ji), float(jf), int(l1), int(l2), float(mixing)
        msi, msf = _ms(self.ji), _ms(self.jf)
        self.a = {}
        for s in (1, -1):
            a = np.zeros((len(msf), len(msi)))
            parts = [(self.l1, 1.0)] + ([(self.l2, -s * self.mixing)] if self.l2 and self.mixing else [])
            for i, m in enumerate(msf):
                big = m + s
                if abs(big) > self.ji + 1e-9:
                    continue
                j = int(round(big + self.ji))
                for lam, c in parts:
                    if angular.triangle(self.jf, lam, self.ji):
                        a[i, j] += c * math.sqrt(2 * lam + 1) * angular.clebsch_gordan(self.jf, m, lam, s, self.ji,
                                                                                       big)
            self.a[s] = a
        # How likely each substate along n is to emit along n: W(n) ∝ Σ_M p[M] (ρ_n)_{MM}.
        p = sum(np.einsum("am,am->m", a, a) for a in self.a.values())
        #: The mean of p over the substates: 4π W(n) = Σ_M p[M] (ρ_n)_{MM} / norm averages to Tr ρ.
        self.norm = float(p.mean())
        self.p = p / self.norm
        #: No 4π W(n) exceeds this, for any ρ of trace 1.
        self.bound = float(self.p.max())

    def weights(self, rho: np.ndarray, directions: np.ndarray) -> np.ndarray:
        """4π W along each direction (n × 3, in the frame of ``rho``) for each density matrix of ``rho``
        (n × (2I + 1)², trace 1): 1 on average over the sphere."""
        n = len(directions)
        out = np.zeros(n)
        for lo in range(0, n, BLOCK):
            sl = slice(lo, lo + BLOCK)
            d = rotation(self.ji, directions[sl])
            diag = (d.conj() * (rho[sl] @ d)).sum(axis=1).real
            out[sl] = diag @ self.p
        return np.maximum(out, 0.0)

    def after(self, rho: np.ndarray, directions: np.ndarray) -> np.ndarray:
        """The density matrix of the level fed (n × (2I_f + 1)², trace 1) when each ρ of ``rho`` emits its γ ray
        along the matching direction."""
        n = len(directions)
        df = int(round(2 * self.jf + 1))
        out = np.zeros((n, df, df), dtype=complex)
        for lo in range(0, n, BLOCK):
            sl = slice(lo, lo + BLOCK)
            d = rotation(self.ji, directions[sl])
            along = np.swapaxes(d.conj(), 1, 2) @ rho[sl] @ d
            fed = sum(a @ along @ a.T for a in self.a.values())
            back = rotation(self.jf, directions[sl])
            out[sl] = back @ fed @ np.swapaxes(back.conj(), 1, 2)
        trace = np.einsum("naa->n", out).real
        good = trace > 0
        out[good] /= trace[good][:, None, None]
        out[~good] = isotropic(self.jf)
        return out

    def sample(self, rho: np.ndarray, rng) -> np.ndarray:
        """A direction for each density matrix of ``rho``, drawn from W(n): unit vectors n × 3 in its frame, by
        acceptance against :attr:`bound`."""
        n = len(rho)
        out = np.zeros((n, 3))
        todo = np.arange(n)
        for _ in range(500):
            if not len(todo):
                break
            u = rng.uniform(-1.0, 1.0, len(todo))
            ph = rng.uniform(0.0, 2 * math.pi, len(todo))
            s = np.sqrt(np.maximum(1 - u * u, 0.0))
            cand = np.c_[s * np.cos(ph), s * np.sin(ph), u]
            w = self.weights(rho[todo], cand)
            ok = rng.uniform(0.0, self.bound, len(todo)) <= w
            out[todo[ok]] = cand[ok]
            todo = todo[~ok]
        return out


@lru_cache(maxsize=None)
def _feeding_map(ji: float, jf: float, lam: int) -> np.ndarray:
    """:func:`~physim.nuclear.orientation.fed_density` as a linear map: out[a, b, c, d] takes ρ[c, d] to
    ρ_f[a, b]."""
    from .orientation import fed_density

    di, df = int(round(2 * ji + 1)), int(round(2 * jf + 1))
    out = np.zeros((df, df, di, di), dtype=complex)
    for c in range(di):
        for d in range(di):
            unit = np.zeros((di, di), dtype=complex)
            unit[c, d] = 1.0
            out[:, :, c, d] = fed_density(unit, ji, jf, lam)
    return out


def unobserved(rho: np.ndarray, ji: float, jf: float, l1: int, l2: int = 0, mixing: float = 0.0) -> np.ndarray:
    """The density matrices a transition leaves when what it emits is not followed (a conversion electron):
    :func:`~physim.nuclear.orientation.fed_density` for each ρ of ``rho``, with the multipoles' shares, trace 1.
    Unoriented when no multipole joins the spins."""
    df = int(round(2 * jf + 1))
    share = [(l1, 1.0 / (1 + mixing**2))] + ([(l2, mixing**2 / (1 + mixing**2))] if l2 and mixing else [])
    out = np.zeros((len(rho), df, df), dtype=complex)
    for lam, w in share:
        if lam and angular.triangle(ji, lam, jf):
            out += w * np.einsum("abcd,ncd->nab", _feeding_map(float(ji), float(jf), int(lam)), rho)
    trace = np.einsum("naa->n", out).real
    good = trace > 0
    out[good] /= trace[good][:, None, None]
    out[~good] = isotropic(jf)
    return out


__all__ = ["BLOCK", "Emission", "isotropic", "rotation", "unobserved"]
