"""How Coulomb excitation leaves a nucleus oriented, and the γ rays that follow: excitation amplitudes for every
magnetic substate, statistical tensors, feeding and decay through the level scheme, and the particle–γ angular
correlation.

First-order theory, as :mod:`physim.nuclear.coulex`, but keeping the amplitudes::

    from physim.nuclear.levels import LevelScheme
    from physim.nuclear.orientation import Excitation, simple_scheme

    scheme = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2")     # 0⁺ → 2⁺ → 0⁺
    ex = Excitation("16O", "58Ni", "30 MeV", scheme, excite="target")
    state = ex.at(150.0)                      # the projectile scattered to 150° in the CM frame
    state.population(1)                       # probability that the 2⁺ state is populated
    state.gamma_yield(1, 0)                   # γ rays of 2⁺ → 0⁺ per collision
    state.w(1, 0, theta, phi)                 # their angular distribution, per steradian, in the orbit's frame

**The orbit's frame.** z is perpendicular to the plane of the orbit, along the angular momentum; x is the
symmetry axis of the hyperbola, pointing from the excited nucleus to the other one at closest approach; y
completes it. The beam direction is (−sin(θ/2), cos(θ/2), 0) and the scattered projectile leaves along
(sin(θ/2), cos(θ/2), 0), θ being the CM scattering angle. :meth:`Populated.axes` gives the frame in the
laboratory.

**The amplitude** of going from the substate M_i of the ground state (spin I_i) to the substate M_f of a level
(spin I_f) by a multipole Eλ is

    a = K_λ ⟨I_f‖M(Eλ)‖I_i⟩ (−1)^{I_f − M_f} (I_f λ I_i; −M_f μ M_i) Y_λμ(π/2, 0) I_λ,−μ(θ, ξ),   μ = M_f − M_i,

with K_λ = 4π Z e² / (ħc (2λ + 1) β a^λ), the orbit integrals of :func:`physim.nuclear.coulex.orbit_integrals`
and Z the charge of the nucleus whose field excites. Averaged over M_i, the amplitudes give the density matrix of
the level and its statistical tensors t_kq, normalised so that t_00 is the population.

**γ rays.** A transition from a level of tensors t_kq has the angular distribution

    W(θ, φ) = (b_γ / 4π) Σ_k A_k Σ_q t_kq √(4π / (2k + 1)) Y_kq(θ, φ),   k = 0, 2, 4, ...

with A_k of :func:`physim.nuclear.angular.distribution_coefficient` and b_γ the share of the level's decays that
give this γ ray (internal conversion takes α / (1 + α) of each transition). A level fed from above takes over the
tensors of its parent times U_k. Everything decays in flight, with the orientation it was given: there is no
deorientation and no lifetime.

Not included: multi-step excitation and reorientation (the excitation of each level is the direct one from the
ground state), which is why the matrix elements enter only through their size.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from . import angular, data
from .coulex import HBARC_MEV_FM, Coulex, orbit_integrals, ylm_equator
from .levels import Level, LevelScheme, Transition, Value, components
from .rutherford import E2_MEV_FM

#: Highest rank of the statistical tensors kept (E3 excitation of a 3⁻ state reaches rank 6).
K_MAX = 6


def simple_scheme(nuclide: str, energy, multipolarity: str = "E2", b_up="1 e2fm4",
                  ground_spin: float = 0.0, ground_parity: int = 1) -> LevelScheme:
    """A two-level scheme: the ground state and one state reached by Eλ, which decays back to it.

    ``b_up`` is B(Eλ↑), as text with a unit or in e² fm^{2λ}. From a 0⁺ ground state the excited state has spin λ
    and parity (−1)^λ; from another ground state it is given the highest spin λ allows."""
    from .coulex import parse_b
    from .quantity import Quantity

    lam = int(multipolarity[1])
    e = Quantity.parse(energy).to("keV") if isinstance(energy, str) else float(energy) * 1e3
    spin = ground_spin + lam
    scheme = LevelScheme(data.nuclide(nuclide).name, reference="one excited state, as in the [reaction] section")
    scheme.levels = [Level(Value(0.0, source="assumed"), ground_spin, ground_parity),
                     Level(Value(e, source="user"), spin, ground_parity * (-1) ** lam)]
    scheme.transitions = [Transition(1, 0, Value(e, source="user"), Value(100.0, source="assumed"), multipolarity)]
    me = math.sqrt(parse_b(b_up, lam) * (2 * ground_spin + 1))
    scheme.set_matrix_element(0, 1, multipolarity, me, source="user")
    return scheme


def _ms(j: float) -> list:
    return [-j + i for i in range(int(round(2 * j)) + 1)]


def tensors_of(rho: np.ndarray, spin: float, k_max: int = K_MAX) -> dict:
    """Statistical tensors {(k, q): t_kq} of a density matrix ρ[M, M′] (M ascending from −I):
    t_kq = √(2I + 1) Σ (−1)^{I − M′} ⟨I M I −M′ | k q⟩ ρ[M, M′], so that t_00 = Tr ρ."""
    ms = _ms(spin)
    out = {}
    for k in range(0, min(k_max, int(round(2 * spin))) + 1):
        for q in range(-k, k + 1):
            total = 0.0 + 0.0j
            for i, m in enumerate(ms):
                mp = m - q
                if abs(mp) > spin + 1e-9:
                    continue
                j = int(round(mp + spin))
                total += ((-1) ** int(round(spin - mp)) * angular.clebsch_gordan(spin, m, spin, -mp, k, q)
                          * rho[i, j])
            out[(k, q)] = math.sqrt(2 * spin + 1) * total
    return out


def fed_density(rho: np.ndarray, spin: float, daughter: float, multipole: int) -> np.ndarray:
    """The density matrix a level of spin ``daughter`` is left with when a level of density matrix ρ decays to
    it by an unobserved radiation of multipolarity L (summed over all directions and polarisations):
    ρ′[m, m′] = Σ ⟨I_f m L μ | I M⟩ ⟨I_f m′ L μ | I M′⟩ ρ[M, M′]. The trace is kept."""
    ms, ds = _ms(spin), _ms(daughter)
    out = np.zeros((len(ds), len(ds)), dtype=complex)
    for a, m in enumerate(ds):
        for b, mp in enumerate(ds):
            for mu in range(-multipole, multipole + 1):
                big, bigp = m + mu, mp + mu
                if abs(big) > spin + 1e-9 or abs(bigp) > spin + 1e-9:
                    continue
                out[a, b] += (angular.clebsch_gordan(daughter, m, multipole, mu, spin, big)
                              * angular.clebsch_gordan(daughter, mp, multipole, mu, spin, bigp)
                              * rho[int(round(big + spin)), int(round(bigp + spin))])
    return out


@dataclass
class _Decay:
    """One transition as the decay uses it."""

    initial: int
    final: int
    #: Share of the level's decays by this transition (γ rays and conversion electrons), and of those that give
    #: its γ ray.
    branch: float
    gamma: float
    #: Multipoles L and L′ (0 if pure) and the mixing ratio δ.
    l1: int
    l2: int
    mixing: float


@dataclass
class Populated:
    """The excited nucleus after one collision at a CM scattering angle: the populations and statistical tensors
    of its levels, in the orbit's frame, with the feeding from higher levels included."""

    excitation: Excitation
    theta_cm: float
    #: Probability of exciting each level directly, by level.
    direct: dict
    #: Statistical tensors with feeding: by level, then by the pair (rank, component).
    tensors: dict
    isotropic: bool = False

    def population(self, level: int) -> float:
        """Probability that the level is populated, directly or by feeding."""
        t = self.tensors.get(level)
        return float(t[(0, 0)].real) if t else 0.0

    def gamma_yield(self, initial: int, final: int) -> float:
        """γ rays of the transition per collision."""
        d = self.excitation.decay(initial, final)
        return self.population(initial) * d.gamma

    def coefficients(self, initial: int, final: int) -> dict:
        """{(k, q): a_kq} of the transition's distribution W = (yield / 4π) Σ a_kq √(4π/(2k+1)) Y_kq, with
        a_00 = 1: the orientation t_kq / t_00 times A_k."""
        d = self.excitation.decay(initial, final)
        t = self.tensors.get(initial)
        if not t or t[(0, 0)].real <= 0:
            return {(0, 0): 1.0 + 0.0j}
        if self.isotropic:
            return {(0, 0): 1.0 + 0.0j}
        ji = self.excitation.scheme.levels[initial].spin
        jf = self.excitation.scheme.levels[final].spin
        out = {}
        for (k, q), value in t.items():
            if k % 2:
                continue
            a_k = angular.distribution_coefficient(k, ji, jf, d.l1, d.l2, d.mixing) if k else 1.0
            if a_k:
                out[(k, q)] = a_k * value / t[(0, 0)].real
        return out

    def w(self, initial: int, final: int, theta, phi) -> np.ndarray:
        """γ rays of the transition per collision and per steradian, in the direction (θ, φ) of the orbit's frame
        (radians), in the rest frame of the emitting nucleus."""
        theta, phi = np.asarray(theta, dtype=float), np.asarray(phi, dtype=float)
        total = np.zeros(np.broadcast(theta, phi).shape)
        for (k, q), a in self.coefficients(initial, final).items():
            y = angular.spherical_harmonic(k, q, theta, phi)
            total = total + (a * math.sqrt(4 * math.pi / (2 * k + 1)) * y).real
        return self.gamma_yield(initial, final) / (4 * math.pi) * total

    def axes(self, phi_particle_deg: float = 0.0) -> np.ndarray:
        """The orbit's frame in the laboratory, for a projectile scattered to this CM angle at azimuth
        ``phi_particle_deg``: a 3 × 3 matrix whose columns are x, y and z of the frame. The beam is along the
        laboratory's z. For projectile excitation x and y are reversed, since the exciting nucleus is then on the
        other side."""
        th, ph = math.radians(self.theta_cm), math.radians(phi_particle_deg)
        beam = np.array([0.0, 0.0, 1.0])
        out = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
        s, c = math.sin(th / 2), math.cos(th / 2)
        # At 0° and 180° the plane is not defined: any plane through the beam at this azimuth serves.
        side = np.array([math.cos(ph), math.sin(ph), 0.0])
        x = (out - beam) / (2 * s) if s > 1e-9 else side
        y = (out + beam) / (2 * c) if c > 1e-9 else side
        z = np.cross(x, y)
        if self.excitation.excite == "projectile":
            x, y = -x, -y
        return np.column_stack([x, y, z])

    def w_lab(self, initial: int, final: int, directions, phi_particle_deg: float = 0.0,
              velocity=None) -> np.ndarray:
        """γ rays of the transition per collision and per steradian in the laboratory, along each unit vector of
        ``directions`` (n × 3).

        ``velocity`` is the emitting nucleus's velocity as a vector in units of c. The γ rays are thrown
        forward: a laboratory direction at angle α to the velocity corresponds to cos α′ = (cos α − β) /
        (1 − β cos α) in the rest frame, and the solid angle shrinks by (1 − β²) / (1 − β cos α)²."""
        d = np.atleast_2d(np.asarray(directions, dtype=float))
        d = d / np.linalg.norm(d, axis=1, keepdims=True)
        factor = np.ones(len(d))
        if velocity is not None:
            v = np.asarray(velocity, dtype=float)
            beta = float(np.linalg.norm(v))
            if beta > 0:
                n = v / beta
                cos_a = d @ n
                gamma = 1 / math.sqrt(1 - beta * beta)
                # The direction in the rest frame: the component along the velocity is boosted back.
                par = (cos_a - beta) / (1 - beta * cos_a)
                perp = d - np.outer(cos_a, n)
                scale = 1 / (gamma * (1 - beta * cos_a))
                d = perp * scale[:, None] + np.outer(par, n)
                factor = (1 - beta * beta) / (1 - beta * cos_a) ** 2
        local = d @ self.axes(phi_particle_deg)
        theta = np.arccos(np.clip(local[:, 2], -1.0, 1.0))
        phi = np.arctan2(local[:, 1], local[:, 0])
        return self.w(initial, final, theta, phi) * factor


class Excitation:
    """First-order Coulomb excitation of every level of a scheme that a matrix element connects to the ground
    state, for a ``beam`` on a ``target`` at ``beam_energy`` (MeV or text). ``excite`` says which of the two the
    scheme belongs to. ``ground`` is the index of the state the nucleus starts in.

    Levels above the CM energy, levels without a spin, and matrix elements that the spins or parities do not
    allow are left out; :attr:`notes` says which."""

    def __init__(self, beam: str, target: str, beam_energy, scheme: LevelScheme, excite: str = "target",
                 ground: int = 0, isotropic: bool = False):
        if excite not in ("target", "projectile"):
            raise ValueError("excite must be 'target' or 'projectile'")
        self.beam, self.target, self.beam_energy = beam, target, beam_energy
        self.scheme, self.excite, self.ground, self.isotropic = scheme, excite, ground, isotropic
        self.notes: list = []
        g = scheme.levels[ground]
        if g.spin is None:
            raise ValueError(f"the ground state of {scheme.nuclide} has no spin in the level scheme")
        self.ground_spin = float(g.spin)
        #: {level: [(λ, ⟨f‖M(Eλ)‖i⟩ in e fm^λ, Coulex of that level)]}: the excitation paths.
        self.paths: dict = {}
        for n, lev in enumerate(scheme.levels):
            if n == ground:
                continue
            for mult in ("E1", "E2", "E3"):
                me = scheme.matrix_element(ground, n, mult)
                if me is None or me.value == 0:
                    continue
                lam = int(mult[1])
                if lev.spin is None:
                    self.notes.append(f"Level {n} ({lev.energy.value:g} keV) has no spin: it is not excited.")
                    continue
                if not angular.triangle(self.ground_spin, lam, lev.spin):
                    self.notes.append(f"Level {n}: {mult} cannot connect spins {g.label} and {lev.label}.")
                    continue
                if None not in (g.parity, lev.parity) and g.parity * lev.parity != (-1) ** lam:
                    self.notes.append(f"Level {n}: {mult} does not connect parities {g.label} and {lev.label}.")
                    continue
                b_up = me.value**2 / (2 * self.ground_spin + 1)
                try:
                    cx = Coulex(beam, target, beam_energy, excite=excite, energy=lev.energy.value * 1e-3,
                                multipolarity=mult, b_up=b_up)
                except ValueError:
                    self.notes.append(f"Level {n} ({lev.energy.value:g} keV) is above the energy available.")
                    continue
                self.paths.setdefault(n, []).append((lam, float(me.value), cx))
        self._decays = self._read_decays()
        self._grid: Optional[tuple] = None

    # -- decay ------------------------------------------------------------------------------------------------------

    def _read_decays(self) -> dict:
        out: dict = {}
        levels = self.scheme.levels
        for n in range(len(levels)):
            branches = [t for t in self.scheme.branches(n) if t.final < n or levels[t.final].energy.value
                        < levels[n].energy.value]
            if not branches:
                continue
            weights = []
            for t in branches:
                i = t.intensity.value if t.intensity is not None else (1.0 if len(branches) == 1 else 0.0)
                alpha = t.conversion.value if t.conversion is not None else 0.0
                weights.append((i, alpha))
            total = sum(i * (1 + a) for i, a in weights)
            if total <= 0:
                self.notes.append(f"Level {n}: its γ rays have no intensities, so its decay is not followed.")
                continue
            for t, (i, alpha) in zip(branches, weights):
                ji, jf = levels[t.initial].spin, levels[t.final].spin
                parts = components(t.multipolarity)
                if ji is None or jf is None:
                    l1, l2, delta = 0, 0, 0.0
                elif parts:
                    l1 = int(parts[0][1])
                    l2 = int(parts[1][1]) if len(parts) == 2 and t.mixing is not None else 0
                    delta = t.mixing.value if l2 else 0.0
                else:
                    # No multipolarity given: the lowest the spins allow.
                    l1, l2, delta = max(int(round(abs(ji - jf))), 1), 0, 0.0
                out[(t.initial, t.final)] = _Decay(t.initial, t.final, i * (1 + alpha) / total, i / total, l1, l2,
                                                   delta)
        return out

    def decay(self, initial: int, final: int) -> _Decay:
        """The transition from ``initial`` to ``final`` as the decay uses it."""
        try:
            return self._decays[(initial, final)]
        except KeyError:
            raise KeyError(f"the level scheme has no transition from level {initial} to level {final}") from None

    def transitions(self) -> list:
        """The (initial, final) pairs of the transitions that are followed."""
        return list(self._decays)

    # -- excitation -------------------------------------------------------------------------------------------------

    def amplitudes(self, theta_cm: float) -> dict:
        """{level: a[M_f, M_i]} in the orbit's frame, M ascending from −I, for the projectile scattered to
        ``theta_cm`` degrees in the CM frame."""
        s = math.sin(math.radians(theta_cm) / 2)
        out = {}
        ms_i = _ms(self.ground_spin)
        for n, paths in self.paths.items():
            jf = float(self.scheme.levels[n].spin)
            ms_f = _ms(jf)
            a = np.zeros((len(ms_f), len(ms_i)), dtype=complex)
            if s > 0:
                for lam, me, cx in paths:
                    integrals = orbit_integrals(lam, 1 / s, cx.xi)       # index μ + λ
                    k = (4 * math.pi * cx.z_exciting * E2_MEV_FM
                         / (HBARC_MEV_FM * (2 * lam + 1) * cx.beta * cx.a**lam))
                    for i, mf in enumerate(ms_f):
                        for j, mi in enumerate(ms_i):
                            mu = int(round(mf - mi))
                            if abs(mu) > lam:
                                continue
                            a[i, j] += (k * me * (-1) ** int(round(jf - mf))
                                        * angular.wigner_3j(jf, lam, self.ground_spin, -mf, mu, mi)
                                        * ylm_equator(lam, mu) * integrals[-mu + lam])
            out[n] = a
        return out

    def _direct(self, theta_cm: float) -> dict:
        """{level: density matrix} of the directly excited levels, averaged over the ground state's substates."""
        return {n: a @ a.conj().T / (2 * self.ground_spin + 1) for n, a in self.amplitudes(theta_cm).items()}

    def _cascade(self, direct: dict) -> dict:
        """Density matrices with feeding: each level, from the top down, hands its density matrix to the levels
        it decays to."""
        rho = {n: m.copy() for n, m in direct.items()}
        levels = self.scheme.levels
        order = sorted(rho.keys() | {d.final for d in self._decays.values()},
                       key=lambda n: -levels[n].energy.value)
        done = set()
        while order:
            n = order.pop(0)
            if n in done or n not in rho:
                done.add(n)
                continue
            done.add(n)
            for (i, f), d in self._decays.items():
                if i != n or levels[f].spin is None or not d.l1:
                    continue
                ji, jf = float(levels[n].spin), float(levels[f].spin)
                share = [(d.l1, 1.0)] if not d.l2 else [(d.l1, 1 / (1 + d.mixing**2)),
                                                        (d.l2, d.mixing**2 / (1 + d.mixing**2))]
                fed = sum(w * fed_density(rho[n], ji, jf, L) for L, w in share if angular.triangle(ji, L, jf))
                if isinstance(fed, np.ndarray):
                    rho[f] = rho.get(f, np.zeros_like(fed)) + d.branch * fed
                    if f not in done and f not in order:
                        order.append(f)
                        order.sort(key=lambda m: -levels[m].energy.value)
        return rho

    def at(self, theta_cm: float) -> Populated:
        """The populations and orientation after a collision at CM angle ``theta_cm`` (degrees), computed for
        this angle."""
        direct = self._direct(theta_cm)
        rho = self._cascade(direct)
        tensors = {n: tensors_of(m, float(self.scheme.levels[n].spin)) for n, m in rho.items()}
        return Populated(self, float(theta_cm), {n: float(np.trace(m).real) for n, m in direct.items()}, tensors,
                         self.isotropic)

    def table(self, step: float = 2.0) -> tuple:
        """(angles, [Populated at each]) on a grid of CM angles from ``step`` to 180°, computed once. With a whole
        number of degrees for ``step`` the orbit integrals are those :class:`~physim.nuclear.coulex.Coulex`
        already keeps for its own table."""
        if self._grid is None or self._grid[0][1] - self._grid[0][0] != step:
            grid = np.arange(step, 180.0 + step / 2, step)
            self._grid = (grid, [self.at(float(t)) for t in grid])
        return self._grid

    def near(self, theta_cm: float, step: float = 2.0) -> Populated:
        """The state at ``theta_cm``, interpolated linearly in the tensors between grid angles ``step`` apart
        (fast for many angles)."""
        grid, states = self.table(step)
        t = min(max(float(theta_cm), grid[0]), grid[-1])
        i = min(int((t - grid[0]) / step), len(grid) - 2)
        f = (t - grid[i]) / step
        lo, hi = states[i], states[i + 1]
        tensors = {n: {key: (1 - f) * v + f * hi.tensors[n][key] for key, v in lo.tensors[n].items()}
                   for n in lo.tensors}
        direct = {n: (1 - f) * v + f * hi.direct[n] for n, v in lo.direct.items()}
        return Populated(self, float(theta_cm), direct, tensors, self.isotropic)

    def cross_sections(self, step: float = 2.0) -> dict:
        """Cross sections integrated over all scattering angles, mb: {"direct": {level: σ}, "populated":
        {level: σ with feeding}, "gamma": {(initial, final): σ of the γ ray}}.

        Each level's excitation probability is multiplied by the Rutherford cross section of its own symmetrised
        orbit and integrated on the grid of :meth:`table`; the decay then shares the cross sections out."""
        grid, states = self.table(step)
        th = np.radians(grid)
        direct = {}
        for n, paths in self.paths.items():
            cx = paths[0][2]
            p = np.array([s.direct.get(n, 0.0) for s in states])
            direct[n] = float(2 * math.pi * np.sum(p * np.asarray(cx.rutherford_cm(grid)) * np.sin(th))
                              * math.radians(step))
        levels = self.scheme.levels
        populated = dict(direct)
        for n in sorted(set(populated) | {d.initial for d in self._decays.values()},
                        key=lambda k: -levels[k].energy.value):
            for (i, f), d in self._decays.items():
                if i == n and populated.get(n):
                    populated[f] = populated.get(f, 0.0) + d.branch * populated[n]
        gamma = {(i, f): populated.get(i, 0.0) * d.gamma for (i, f), d in self._decays.items()}
        return {"direct": direct, "populated": populated, "gamma": gamma}

    def probability(self, level: int, theta_cm: float) -> float:
        """Probability of exciting ``level`` directly at CM angle ``theta_cm``: the amplitudes squared, summed
        over the final substates and averaged over the initial ones."""
        return self.at(theta_cm).direct.get(level, 0.0)


__all__ = ["Excitation", "K_MAX", "Populated", "fed_density", "simple_scheme", "tensors_of"]
