"""Multi-step Coulomb excitation: the coupled equations for the amplitudes of every magnetic substate of every
level, integrated along the orbit. Reorientation (the diagonal E2 matrix elements, that is the quadrupole
moments) and interference between paths follow from it, where first order (:mod:`physim.nuclear.orientation`)
has neither.

It is the same interface as first order, so everything built on :class:`~physim.nuclear.orientation.Excitation`
(populations, tensors, decay, the correlation, the cross sections) works unchanged::

    from physim.nuclear.coupled import CoupledChannels

    cc = CoupledChannels("16O", "194Pt", "60 MeV", scheme, excite="target")
    state = cc.at(150.0)                      # populations and tensors with all orders
    cc.cross_sections()                       # mb per level and per γ ray

**The equations.** In the orbit's frame (see :mod:`physim.nuclear.orientation`), with the orbit written as
r = a(ε cosh w + 1), t = (a/v)(ε sinh w + w),

    da_n/dw = −i Σ_m Σ_λμ K_λ ⟨n|M(Eλ, μ)|m⟩ Y_λμ(π/2, 0) e^{−iμφ(w)} e^{i ξ_nm (ε sinh w + w)} / (ε cosh w + 1)^λ a_m,

K_λ = 4π Z e² / (ħc β (2λ+1) a^λ), ξ_nm = (Z₁Z₂e²/ħc)(1/β_n − 1/β_m) with β_n the CM speed after exciting level n,
and ⟨n M_n|M(Eλ, μ)|m M_m⟩ = (−1)^{I_n − M_n} (I_n λ I_m; −M_n μ M_m) ⟨I_n‖M(Eλ)‖I_m⟩ (Wigner–Eckart). The reduced
matrix elements are those of the level scheme, with ⟨m‖M‖n⟩ = (−1)^{I_n − I_m} ⟨n‖M‖m⟩; a diagonal E2 element is
the quadrupole moment's. The orbit (a, β) is the symmetrised one of the ground state and a reference level, so
with one level and weak coupling the solution is first order exactly.

The equations are integrated by a fourth-order Runge–Kutta scheme in w, with steps that follow the fastest phase,
from where the coupling is negligible on the way in to where it is on the way out. The probabilities of all
levels sum to one at every step (the coupling is Hermitian).

Not included: the E1 polarisation correction (virtual excitation of the giant dipole resonance), deorientation,
lifetimes, mutual excitation of both nuclei, and nuclear interference above the safe energy.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from . import angular
from .coulex import HBARC_MEV_FM, ylm_equator
from .orientation import Excitation, _ms
from .rutherford import E2_MEV_FM

#: Reduced matrix elements the scheme keeps.
ELECTRIC = ("E1", "E2", "E3")


class CoupledChannels(Excitation):
    """Coulomb excitation of a level scheme to all orders. ``reference`` is the level whose symmetrised orbit is
    used (the lowest excited level with a matrix element to the ground state if None); ``tolerance`` sets the
    step of the integration (the relative accuracy of the probabilities, roughly); ``max_levels`` caps the
    scheme (the lowest levels are kept)."""

    def __init__(self, beam: str, target: str, beam_energy, scheme, excite: str = "target", ground: int = 0,
                 isotropic: bool = False, reference: Optional[int] = None, tolerance: float = 1e-5,
                 max_levels: int = 30):
        super().__init__(beam, target, beam_energy, scheme, excite, ground, isotropic)
        self.tolerance = tolerance
        levels = scheme.levels
        keep = [n for n in range(len(levels)) if levels[n].spin is not None]
        keep = sorted(keep, key=lambda n: levels[n].energy.value)[:max_levels]
        if ground not in keep:
            raise ValueError("the ground state has no spin in the level scheme")
        self.included = sorted(keep)
        # The basis: every substate of every included level.
        self.basis = [(n, m) for n in self.included for m in _ms(float(levels[n].spin))]
        self._index = {b: i for i, b in enumerate(self.basis)}
        self._level_slices = {n: [i for i, (k, _) in enumerate(self.basis) if k == n] for n in self.included}
        # The orbit: the symmetrised one of the ground state and the reference level.
        self.reference = reference if reference is not None else (min(self.paths) if self.paths else None)
        if self.reference is None or self.reference not in self.paths:
            raise ValueError("no level is connected to the ground state by a matrix element")
        cx = self.paths[self.reference][0][2]
        self.cx = cx
        self.a, self.beta, self.k = cx.a, cx.beta, cx.k
        self.z_exciting = cx.z_exciting
        mu = cx.mu
        e_cm = cx.e_cm
        self._beta_level = {}
        for n in self.included:
            e_n = levels[n].energy.value * 1e-3 - levels[ground].energy.value * 1e-3
            if e_cm - e_n <= 0:
                self.notes.append(f"Level {n} is above the energy available and is left out of the coupling.")
                continue
            self._beta_level[n] = math.sqrt(2 * (e_cm - e_n) / mu)
        self.included = [n for n in self.included if n in self._beta_level]
        self.basis = [(n, m) for n in self.included for m in _ms(float(levels[n].spin))]
        self._index = {b: i for i, b in enumerate(self.basis)}
        self._level_slices = {n: [i for i, (k, _) in enumerate(self.basis) if k == n] for n in self.included}
        self._couplings = self._build_couplings()
        self._mats = np.array([mat for _, _, mat in self._couplings]) if self._couplings else np.zeros((0, 1, 1))
        n = len(self.basis)
        #: ξ of each basis state, k/(ħc β_n): ξ_nm = ξ_n − ξ_m.
        self.xi_state = np.array([self.k / HBARC_MEV_FM / self._beta_level[ni] for ni, _ in self.basis])
        self.xi = self.xi_state[:, None] - self.xi_state[None, :]
        self._grid = None

    # -- the coupling matrices ------------------------------------------------------------------------------------

    def _build_couplings(self) -> list:
        """[(λ, μ, K_λ Y_λμ(π/2,0) × matrix of Wigner–Eckart factors × reduced matrix elements)] for every
        multipole component with a non-zero matrix."""
        levels = self.scheme.levels
        n = len(self.basis)
        out = []
        for mult in ELECTRIC:
            lam = int(mult[1])
            k_lam = 4 * math.pi * self.z_exciting * E2_MEV_FM / (HBARC_MEV_FM * self.beta * (2 * lam + 1) * self.a**lam)
            reduced = {}
            for a in self.included:
                for b in self.included:
                    if b < a:
                        continue
                    me = self.scheme.matrix_element(a, b, mult)
                    if me is None or me.value == 0:
                        continue
                    ia, ib = float(levels[a].spin), float(levels[b].spin)
                    if not angular.triangle(ia, lam, ib):
                        continue
                    pa, pb = levels[a].parity, levels[b].parity
                    if None not in (pa, pb) and pa * pb != (-1) ** lam:
                        continue
                    # ⟨b‖M‖a⟩ as stored (a ≤ b); the other direction by the symmetry of real matrix elements.
                    reduced[(b, a)] = me.value
                    reduced[(a, b)] = (-1) ** int(round(ia - ib)) * me.value
            if not reduced:
                continue
            for mu in range(-lam, lam + 1):
                y = ylm_equator(lam, mu)
                mat = np.zeros((n, n))
                for i, (ni, mi) in enumerate(self.basis):
                    for j, (nj, mj) in enumerate(self.basis):
                        if (ni, nj) not in reduced or abs(mi - mj - mu) > 1e-9:
                            continue
                        ii = float(levels[ni].spin)
                        jj = float(levels[nj].spin)
                        mat[i, j] = ((-1) ** int(round(ii - mi)) * angular.wigner_3j(ii, lam, jj, -mi, mu, mj)
                                     * reduced[(ni, nj)])
                if np.any(mat):
                    out.append((lam, mu, k_lam * y * mat))
        return out

    # -- the integration ----------------------------------------------------------------------------------------

    def _coupling(self, w: np.ndarray, eps: float) -> np.ndarray:
        """F(w) for an array of w values: shape (len(w), N, N), so that da/dw = −i F a (for checks; the
        integration applies F without forming it, see :meth:`_apply`)."""
        ch, sh = np.cosh(w), np.sinh(w)
        r = eps * ch + 1.0
        e_iphi = ((ch + eps) + 1j * math.sqrt(max(eps * eps - 1.0, 0.0)) * sh) / r
        phase = np.exp(1j * self.xi[None, :, :] * (eps * sh + w)[:, None, None])
        f = np.zeros((len(w), len(self.basis), len(self.basis)), dtype=complex)
        for lam, mu, mat in self._couplings:
            f += (e_iphi ** (-mu) / r**lam)[:, None, None] * mat[None, :, :]
        return f * phase

    def _operator(self, w: float, eps: float) -> tuple:
        """(c, p) at one w: the coefficient c_λμ(w) = e^{−iμφ} / r^λ of each coupling matrix, and the phases
        p_n = exp(i ξ_n (ε sinh w + w)) of the basis states, so that F a = p ⊙ Σ c_λμ mat_λμ (p̄ ⊙ a)."""
        ch, sh = math.cosh(w), math.sinh(w)
        r = eps * ch + 1.0
        e_iphi = ((ch + eps) + 1j * math.sqrt(max(eps * eps - 1.0, 0.0)) * sh) / r
        c = np.array([e_iphi ** (-mu) / r**lam for lam, mu, _ in self._couplings])
        p = np.exp(1j * self.xi_state * (eps * sh + w))
        return c, p

    def _apply(self, op: tuple, a: np.ndarray) -> np.ndarray:
        c, p = op
        x = np.conj(p)[:, None] * a
        # Each coupling matrix is real: K real matrix products on the real and imaginary parts, then the K
        # complex coefficients.
        mx = self._mats @ np.ascontiguousarray(x.real) + 1j * (self._mats @ np.ascontiguousarray(x.imag))
        y = np.tensordot(c, mx, axes=(0, 0))
        return p[:, None] * y

    def _w_range(self, eps: float) -> tuple:
        """Where to start and stop: beyond ±w_max the coupling is below the tolerance (the fastest phase has
        averaged it away, or 1/r^λ is small)."""
        tol = self.tolerance * 1e-1
        lam_min = min(lam for lam, _, _ in self._couplings)
        w_amp = math.acosh(max((1 / tol) ** (1 / lam_min) / eps, 1.0)) + 1.0
        xi_max = float(np.max(np.abs(self.xi)))
        if xi_max > 0 and tol * xi_max < 1:
            w_osc = math.log(2 * (1 / (tol * xi_max)) ** (1 / (lam_min + 1)) / eps)
            return -min(w_amp, max(w_osc, 3.0) + 1.0), min(w_amp, max(w_osc, 3.0) + 1.0)
        return -w_amp, w_amp

    def _steps(self, eps: float) -> np.ndarray:
        """The grid of w: steps small enough for the fastest phase, ξ_max (ε cosh w + 1), and for the shape of
        1/r^λ near the closest approach."""
        lo, hi = self._w_range(eps)
        xi_max = float(np.max(np.abs(self.xi)))
        c = 0.3 * (self.tolerance / 1e-5) ** 0.25
        pts = [0.0]
        while pts[-1] < hi:
            w = pts[-1]
            # Near the closest approach 1/r^λ varies on the scale of 1/ε in w; further out only the phase does.
            h = min(0.05 * c * (eps * math.cosh(w) + 1.0) / (eps + 1.0), 0.3 * c,
                    c / (xi_max * (eps * math.cosh(w) + 1.0) + 1.0))
            pts.append(min(w + h, hi))
        right = np.array(pts)
        return np.concatenate([-right[:0:-1], right])

    def solve(self, theta_cm: float) -> np.ndarray:
        """The amplitudes after the collision, a[basis, initial substate], for the projectile scattered to
        ``theta_cm`` degrees in the CM frame (orbit's frame)."""
        s = math.sin(math.radians(theta_cm) / 2)
        n = len(self.basis)
        g = self._level_slices[self.ground]
        a = np.zeros((n, len(g)), dtype=complex)
        for col, i in enumerate(g):
            a[i, col] = 1.0
        if s <= 0 or not self._couplings:
            return a
        eps = 1 / s
        w = self._steps(eps)
        op_next = self._operator(w[0], eps)
        for k in range(len(w) - 1):
            h = w[k + 1] - w[k]
            op0 = op_next
            op1 = self._operator(0.5 * (w[k] + w[k + 1]), eps)
            op_next = self._operator(w[k + 1], eps)
            k1 = -1j * self._apply(op0, a)
            k2 = -1j * self._apply(op1, a + 0.5 * h * k1)
            k3 = -1j * self._apply(op1, a + 0.5 * h * k2)
            k4 = -1j * self._apply(op_next, a + h * k3)
            a = a + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return a

    def amplitudes(self, theta_cm: float) -> dict:
        """{level: a[M_f, M_i]} for every included level, from the coupled solution (the ground state too)."""
        a = self.solve(theta_cm)
        return {n: a[self._level_slices[n], :] for n in self.included}

    def probability(self, level: int, theta_cm: float) -> float:
        """Probability that the collision leaves the nucleus in ``level`` (all orders), averaged over the ground
        state's substates."""
        a = self.amplitudes(theta_cm).get(level)
        if a is None:
            return 0.0
        return float(np.sum(np.abs(a) ** 2) / a.shape[1])

    def probabilities(self, theta_cm: float) -> dict:
        """{level: probability} for every included level; they sum to one."""
        a = self.solve(theta_cm)
        return {n: float(np.sum(np.abs(a[self._level_slices[n], :]) ** 2) / a.shape[1]) for n in self.included}

    def _direct(self, theta_cm: float) -> dict:
        """{level: density matrix} after the collision, excited levels only (the ground state's population is
        not a yield)."""
        amps = self.amplitudes(theta_cm)
        n0 = len(self._level_slices[self.ground])
        return {n: a @ a.conj().T / n0 for n, a in amps.items() if n != self.ground}

    def second_order_check(self, level: int, theta_cm: float) -> tuple:
        """For the tests: the first- and second-order perturbation amplitudes of ``level`` from the ground state,
        by nested quadrature on the same grid, as (first, second) arrays over the level's substates (initial
        substate M_i = the first). Second order is −∫ dw F(w) ∫_{w' < w} dw' F(w') a₀."""
        s = math.sin(math.radians(theta_cm) / 2)
        eps = 1 / s
        w = self._steps(eps)
        f = self._coupling(w, eps)
        g0 = self._level_slices[self.ground][0]
        a0 = np.zeros(len(self.basis), dtype=complex)
        a0[g0] = 1.0
        # First order: cumulative integral of F a₀.
        integrand = f @ a0
        cum = np.zeros_like(integrand)
        for k in range(1, len(w)):
            cum[k] = cum[k - 1] + 0.5 * (w[k] - w[k - 1]) * (integrand[k] + integrand[k - 1])
        first = -1j * cum[-1]
        inner = cum  # ∫_{−∞}^{w} F a₀
        second_integrand = np.einsum("kij,kj->ki", f, inner)
        second = -np.trapezoid(second_integrand, w, axis=0) if hasattr(np, "trapezoid") else \
            -np.trapz(second_integrand, w, axis=0)
        idx = self._level_slices[level]
        return first[idx], second[idx]


__all__ = ["CoupledChannels", "ELECTRIC"]
