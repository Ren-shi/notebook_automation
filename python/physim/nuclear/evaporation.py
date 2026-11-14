"""Statistical decay of a compound nucleus: the app's own evaporation Monte Carlo (backlog item 75, Stage 3).

A hot compound nucleus boils off neutrons, protons and α particles, emits statistical γ rays and, if it is heavy,
may fission. Each event follows one random chain, as PACE4 does: at every step the open decays compete by their
widths in the Hauser–Feshbach picture, one is chosen, and the emitted particle's energy, orbital angular momentum
and the daughter's spin are drawn. The chain stops when particles have less than a thousandth of the total width,
or when the state sits on the yrast line; the state left after the last particle is the entry state from which the
residue's discrete γ cascade starts. The model and its sources are in ``docs/theory/evaporation.md``.

::

    from physim.nuclear.evaporation import evaporate

    r = evaporate("162Dy", 59.3, sigma_l_mb=sigma_l, events=100_000, seed=1)   # 12C + 150Nd at 70 MeV
    r.channel_table()          # [{"label": "5n", "residue": "157Dy", "share": 0.65, "sigma_mb": 660.0, ...}, ...]
    r.mean_energy_mev("n")     # mean neutron kinetic energy in the compound nucleus's frame: 2.4 MeV
    r.particles(0)             # the emissions of event 0

Input
-----
``evaporate(compound, excitation_mev, *, sigma_l_mb=None, spin=None, spin_weights=None, events=10000, seed=None,
decays=DECAYS, model=None)``

- ``compound``: the compound nucleus, as a name (``"172Yb"``) or ``(Z, A)``.
- ``excitation_mev``: its excitation energy E*, MeV.
- The spin distribution, exactly one of: ``sigma_l_mb``, partial-wave cross sections σ_ℓ in mb indexed by ℓ (what
  ``fusion.py`` returns; their sum is σ_fus and scales the channel cross sections); ``spin``, a single J in ħ;
  ``spin_weights``, relative weights indexed by J (no cross sections then). The compound nucleus's spin is taken
  equal to ℓ: the beam's and target's ground-state spins are neglected.
- ``events``: the number of chains; ``seed``: an int, a ``numpy.random.SeedSequence`` or a ``numpy.random.Generator``.
  The same seed gives bit-identical results.
- ``decays``: the decays allowed, a subset of ``("n", "p", "a", "g", "fission")``. Turning off ``"p"``, ``"a"``,
  ``"g"`` and ``"fission"`` leaves pure neutron evaporation. A decay whose daughter is not in the AME2020 mass
  table is closed.
- ``model``: an :class:`EvaporationModel` with the physics parameters (level-density parameter, grid steps,
  radii, GDR width, fission-barrier scaling). The defaults are documented on the class.

Output
------
An :class:`Evaporation` with NumPy arrays. Per event (length ``events``):

- ``initial_spin`` (J of the compound nucleus, ħ), ``residue_z``, ``residue_a``, ``channel`` (an index into
  ``labels``; ``label_of_event`` gives the strings such as ``"p2n"``, ``"a3n"``, ``"4n"``, ``"g"`` for no particle,
  ``"fission"``), ``fission`` (bool).
- ``entry_energy_mev`` and ``entry_spin``: the state (E*, J) of the residue after the last particle. For a fission
  event, the state of the nucleus that fissioned.
- ``recoil_momentum_mev_c`` (events × 3): the residue's momentum in the compound nucleus's rest frame, MeV/c.

Per emission, in flat arrays, event ``i`` owning ``offsets[i]:offsets[i+1]`` in emission order:

- ``kind``: 0 neutron, 1 proton, 2 α, 3 γ (``KINDS``). Statistical γ rays emitted between particles are listed, so
  that energy and momentum balance exactly; γ rays after the last particle are not (they belong to the cascade
  that starts from the entry state).
- ``energy_mev``: the kinetic energy (the photon's energy for a γ) in the compound nucleus's rest frame, MeV.
- ``momentum_mev_c`` (emissions × 3) and ``direction`` (unit vectors) in the same frame.
- ``event``: the event each emission belongs to.

``lab(beta)`` boosts the emissions to the lab for a compound nucleus moving with β along the beam (+z);
``events_list(beta)`` gives each event as ``{"channel", "particles": [(kind, energy_mev, direction), ...]}``, the
form ``neutrons.simulate_selection`` reads.

Summaries: ``channel_table()`` (residue, share, cross section in mb when σ_fus is known, events, statistical error),
``shares()``, ``mean_energy_mev(kind)``, ``multiplicity(kind)``. Kinematics is relativistic and exact: in every event
the emitted particles and the residue (with its entry excitation) carry the compound nucleus's mass-energy
M_CN + E* and zero total momentum.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Optional, Sequence, Union

import numpy as np

from . import data
from .names import SYMBOLS

#: ħc, MeV fm (CODATA 2018).
HBARC = 197.3269804
#: e²/4πε₀, MeV fm.
E2 = 1.439964548
#: Atomic mass unit, MeV.
U_MEV = data.U_KEV * 1e-3
#: Emission kinds, in the order of the ``kind`` codes.
KINDS = ("n", "p", "a", "g")
#: Decays that can be switched on or off.
DECAYS = ("n", "p", "a", "g", "fission")
#: (Z, A) of the emitted particles.
PARTICLES = {"n": (0, 1), "p": (1, 1), "a": (2, 4)}
#: The α particle's binding energy, MeV (AME2020), for the bound on each nucleus's excitation energy.
_ALPHA_BINDING = 28.2957

NuclideLike = Union[str, tuple]


@dataclass(frozen=True)
class EvaporationModel:
    """Physics parameters of the evaporation Monte Carlo. Energies MeV, lengths fm.

    Level densities (back-shifted Fermi gas on a rigid-rotor yrast line):

    - ``a_divisor``: a = A / a_divisor, MeV⁻¹. 10 is PACE4's default and suits excitation energies of tens of MeV
      (a ≈ A/8 is the low-energy value; a falls towards A/10–A/12 as the nucleus heats up).
    - ``pairing_mev``: Δ = pairing_mev / √A. The back shift is 2Δ for even-even, Δ for odd-A and 0 for odd-odd
      nuclei (Gilbert and Cameron's pairing energies, smoothed).
    - ``r0_inertia``: the rigid-sphere moment of inertia I = (2/5) m A (r0 A^⅓)², which sets the yrast line
      E = ħ²J²/2I and the spin cut-off.
    - ``thermal_floor_mev``: the Fermi-gas prefactor U^(-3/2) is evaluated at U + this, so it stays finite next
      to the yrast line.
    - ``constant_temperature``: below the matching energy U_x = 2.5 + 150/A MeV, Gilbert and Cameron's constant
      temperature form, matched in value and slope (the pure Fermi gas has far too few levels there, which lets
      γ rays win over the last neutron).

    Transmission coefficients:

    - ``r0_neutron``, ``k0_neutron``: the black (strongly absorbing) nucleus of Blatt and Weisskopf for neutrons,
      radius r0 A^⅓ and inner wave number k0.
    - Protons: the real Becchetti–Greenlees potential; α particles: the Huizenga–Igo potential. For each ℓ the
      barrier of Coulomb + nuclear + centrifugal gives T_ℓ: Hill–Wheeler (parabolic) above its top, WKB below.

    γ rays: ``gdr_width_mev``, the width of the giant dipole resonance (Brink–Axel Lorentzian with the
    Thomas–Reiche–Kuhn strength and the Berman–Fultz energy).

    Fission (Bohr–Wheeler): ``fission_barrier_scale`` multiplies the liquid-drop barrier;
    ``saddle_inertia_ratio`` is I_saddle / I_ground, which lowers the barrier with spin; ``af_over_an`` is the
    level-density parameter at the saddle over that at the ground state.

    Grids: ``energy_step_mev`` is the step of the excitation-energy and emission-energy grids;
    ``max_emission_mev`` caps the emission energy (the spectra are negligible above it); ``max_steps`` caps the
    chain length. The chain ends when particles have less than ``stop_ratio`` of the total width (PACE4 has a
    similar discriminator, 0.002), or when no statistical γ ray is possible (the state is on the yrast line): the γ
    rays that follow belong to the cascade from the entry state.

    Particle spins are not coupled: the daughter's spin J' runs over |J − ℓ| … J + ℓ, and the widths of neutrons
    and protons are multiplied by 2s + 1 = 2 for the spin states.
    """

    a_divisor: float = 10.0
    pairing_mev: float = 12.0
    r0_inertia: float = 1.2
    thermal_floor_mev: float = 1.0
    constant_temperature: bool = True
    r0_neutron: float = 1.3
    k0_neutron: float = 1.0
    gdr_width_mev: float = 5.0
    fission_barrier_scale: float = 1.0
    saddle_inertia_ratio: float = 2.0
    af_over_an: float = 1.0
    energy_step_mev: float = 0.5
    max_emission_mev: float = 50.0
    max_steps: int = 400
    stop_ratio: float = 1e-3


# ---------------------------------------------------------------------------------------------------------------
# Ingredients, each usable on its own


def nuclear_mass_mev(z: int, a: int) -> Optional[float]:
    """Nuclear mass from AME2020, MeV, or None when the nuclide is not in the table."""
    try:
        return data.nuclide((z, a)).nuclear_mass_mev
    except ValueError:
        return None


def level_density_parameter(a: int, model: EvaporationModel = EvaporationModel()) -> float:
    """Fermi-gas level-density parameter a = A / a_divisor, MeV⁻¹."""
    return a / model.a_divisor


def pairing_shift_mev(z: int, a: int, model: EvaporationModel = EvaporationModel()) -> float:
    """Back shift δ of the Fermi gas, MeV: 2Δ even-even, Δ odd-A, 0 odd-odd, with Δ = pairing_mev/√A."""
    delta = model.pairing_mev / math.sqrt(a)
    return delta * ((z % 2 == 0) + ((a - z) % 2 == 0))


def rotational_constant_mev(a: int, model: EvaporationModel = EvaporationModel(), ratio: float = 1.0) -> float:
    """ħ²/2I of a rigid sphere, MeV, with I = ratio × (2/5) m_u A (r0 A^⅓)²."""
    inertia = ratio * 0.4 * U_MEV * a * (model.r0_inertia * a ** (1 / 3)) ** 2  # MeV fm² / c²
    return HBARC**2 / (2 * inertia)


def yrast_energy_mev(a: int, spin, model: EvaporationModel = EvaporationModel()):
    """Rigid-rotor yrast energy ħ²J²/2I, MeV."""
    return rotational_constant_mev(a, model) * np.asarray(spin, dtype=float) ** 2


def _fermi(u, a, rot, floor):
    """Fermi-gas levels per MeV with projection M at thermal energy u above the yrast energy of M (0 for u ≤ 0)."""
    u = np.asarray(u, dtype=float)
    up = np.maximum(u, 0.0)
    out = math.sqrt(rot) / 12.0 * (up + floor) ** -1.5 * np.exp(2.0 * np.sqrt(a * up))
    return np.where(u > 0, out, 0.0)


class _LevelDensity:
    """ρ_M(E, M) of one nucleus: Gilbert–Cameron, a constant temperature below the matching energy and the
    back-shifted Fermi gas above it, both measured from the rigid-rotor yrast energy of M."""

    def __init__(self, z, a, model, constant_temperature=True):
        self.lda = level_density_parameter(a, model)
        self.delta = pairing_shift_mev(z, a, model)
        self.rot = rotational_constant_mev(a, model)
        self.floor = model.thermal_floor_mev
        self.ct = constant_temperature and model.constant_temperature
        # Matching at U_x = 2.5 + 150/A MeV (Gilbert and Cameron): equal value and slope.
        self.ux = 2.5 + 150.0 / a
        self.tx = 1.0 / (math.sqrt(self.lda / self.ux) - 1.5 / (self.ux + self.floor))
        self.fx = float(_fermi(self.ux, self.lda, self.rot, self.floor))

    def projection(self, e, m):
        """ρ_M(E, M), levels per MeV with projection M at excitation e (arrays broadcast)."""
        u = np.asarray(e, dtype=float) - self.delta - self.rot * np.asarray(m, dtype=float) ** 2
        fg = _fermi(u, self.lda, self.rot, self.floor)
        if not self.ct:
            return fg
        ct = self.fx * np.exp((np.minimum(u, self.ux) - self.ux) / self.tx)
        return np.where(u >= self.ux, fg, np.where(u > -self.delta, ct, 0.0))

    def levels(self, e, j):
        """ρ(E, J) = ρ_M(E, J) − ρ_M(E, J+1), levels of spin J per MeV."""
        j = np.asarray(j, dtype=float)
        return self.projection(e, j) - self.projection(e, j + 1)


def level_density(z: int, a: int, excitation_mev, spin, model: EvaporationModel = EvaporationModel()):
    """Density of levels of spin J (each counted once), MeV⁻¹: ρ(E, J) = ρ_M(E, J) − ρ_M(E, J+1).

    ρ_M(E, M) is the density of states with angular-momentum projection M. Above the matching energy it is the
    back-shifted Fermi gas, ρ_M = (ħ²/2I)^½ / 12 · (U + floor)^(-3/2) · exp(2√(aU)) with U = E − δ − ħ²M²/2I, whose
    yrast line is the rigid rotor's; for U ≫ ħ²J²/2I this is the familiar ρ(E, J) ∝ (2J+1) exp(−J²/2σ²) ρ(E) with
    σ² = I T/ħ². Below U_x = 2.5 + 150/A MeV it is Gilbert and Cameron's constant temperature, exp(U/T), matched in
    value and slope, down to the yrast line (``EvaporationModel.constant_temperature`` turns this part off).
    """
    return _LevelDensity(z, a, model).levels(excitation_mev, spin)


def neutron_transmission(epsilon_mev, a_daughter: int, lmax: int, model: EvaporationModel = EvaporationModel()):
    """Neutron transmission coefficients T_ℓ(ε) of the black nucleus, shape (lmax+1, len(ε)).

    Blatt and Weisskopf (1952), ch. VIII: a purely ingoing wave inside the radius R with wave number
    K = √(k² + k0²) gives T_ℓ = 4 P_ℓ KR / ((KR + P_ℓ)² + S_ℓ²), with the penetration P_ℓ and shift S_ℓ of the free
    spherical waves at ρ = kR. At low energy T₀ ∝ √ε, the 1/v law.
    """
    eps = np.atleast_1d(np.asarray(epsilon_mev, dtype=float))
    mn = nuclear_mass_mev(0, 1)
    md = a_daughter * U_MEV
    mu = mn * md / (mn + md)
    k = np.sqrt(2 * mu * np.maximum(eps, 0.0)) / HBARC
    radius = model.r0_neutron * a_daughter ** (1 / 3)
    x = np.sqrt(k**2 + model.k0_neutron**2) * radius
    rho = k * radius
    t = np.empty((lmax + 1, eps.size))
    p, s = rho.copy(), np.zeros_like(rho)
    for ell in range(lmax + 1):
        if ell:
            d = (ell - s) ** 2 + p**2
            p, s = rho**2 * p / d, rho**2 * (ell - s) / d - ell
        t[ell] = 4 * p * x / ((x + p) ** 2 + s**2)
    return t


def _optical_real(particle: str, z_d: int, a_d: int):
    """(V0 MeV, R fm, a fm, Rc fm) of the real potential for a charged particle on the daughter (Z, A)."""
    third = a_d ** (1 / 3)
    if particle == "p":  # Becchetti and Greenlees, Phys. Rev. 182, 1190 (1969), at E = 0
        n = a_d - z_d
        return 54.0 + 0.4 * z_d / third + 24.0 * (n - z_d) / a_d, 1.17 * third, 0.75, 1.3 * third
    if particle == "a":  # Huizenga and Igo, Nucl. Phys. 29, 462 (1962)
        return 50.0, 1.17 * third + 1.77, 0.576, 1.17 * third
    raise ValueError(particle)


def _potential(particle: str, z_d: int, a_d: int, r, ell):
    """Coulomb (uniform sphere) + real Woods–Saxon + centrifugal potential, MeV, shape (len(ell), len(r))."""
    zp, ap = PARTICLES[particle]
    v0, rn, an, rc = _optical_real(particle, z_d, a_d)
    mu = ap * a_d / (ap + a_d) * U_MEV
    zz = zp * z_d * E2
    vc = np.where(r >= rc, zz / np.maximum(r, 1e-9), zz * (3 - (r / rc) ** 2) / (2 * rc))
    vn = -v0 / (1 + np.exp((r - rn) / an))
    ell = np.asarray(ell, dtype=float)[:, None]
    return vc + vn + HBARC**2 * ell * (ell + 1) / (2 * mu * r**2), mu, zz, rn


def barriers(particle: str, z_d: int, a_d: int, lmax: int):
    """Barrier heights V_ℓ (MeV), curvatures ħω_ℓ (MeV) and radii (fm) for a charged particle, ℓ = 0..lmax.

    The potential is Coulomb (uniform sphere of radius Rc) + the real Woods–Saxon + ħ²ℓ(ℓ+1)/2μr², searched outside
    the Woods–Saxon radius R (inside it the particle is absorbed). Where the pocket has gone at high ℓ, the top is
    at R and the curvature is that of the last ℓ with a pocket.
    """
    v0, rn, an, rc = _optical_real(particle, z_d, a_d)
    r = np.linspace(rn, rn + 25.0, 2501)
    h = r[1] - r[0]
    v, mu, _, _ = _potential(particle, z_d, a_d, r, np.arange(lmax + 1))
    i = np.argmax(v, axis=1)
    vb = v[np.arange(lmax + 1), i]
    rb = r[i]
    hw = np.empty(lmax + 1)
    last = None
    for ell in range(lmax + 1):
        k = i[ell]
        if 0 < k < r.size - 1:
            curv = (v[ell, k + 1] - 2 * v[ell, k] + v[ell, k - 1]) / h**2
            hw[ell] = HBARC * math.sqrt(max(-curv, 1e-6) / mu)
            last = ell
        else:
            if last is None:
                raise ValueError(f"no barrier for {particle} on Z={z_d}, A={a_d}")
            hw[ell] = hw[last]
    return vb, hw, rb


def charged_transmission(particle: str, z_d: int, a_d: int, epsilon_mev, lmax: int):
    """Transmission coefficients T_ℓ(ε) of a proton or α particle, shape (lmax+1, len(ε)).

    Above the barrier top V_ℓ: the parabolic barrier of Hill and Wheeler, T = 1/(1 + exp(2π(V_ℓ − ε)/ħω_ℓ)). Below
    it: the WKB form of Kemble, T = 1/(1 + exp(2θ)), θ = ∫ √(2μ(V − ε))/ħ dr across the barrier, which keeps the
    long Coulomb tail that a parabola misses. Both give ½ at the top.
    """
    eps = tuple(np.atleast_1d(np.asarray(epsilon_mev, dtype=float)))
    return _charged_transmission(particle, z_d, a_d, eps, lmax).copy()


@lru_cache(maxsize=4096)
def _charged_transmission(particle: str, z_d: int, a_d: int, epsilon_mev: tuple, lmax: int):
    eps = np.asarray(epsilon_mev)
    vb, hw, rb = barriers(particle, z_d, a_d, lmax)
    rn = _optical_real(particle, z_d, a_d)[1]
    r_end = 400.0
    r = np.geomspace(rn, r_end, 700)
    dr = np.gradient(r)
    v, mu, zz, _ = _potential(particle, z_d, a_d, r, np.arange(lmax + 1))
    k2 = 2 * mu / HBARC**2
    t = np.empty((lmax + 1, eps.size))
    e = np.maximum(eps, 1e-6)
    # Coulomb tail beyond r_end, to the outer turning point zz/ε: ∫ √(zz/r − ε) dr in closed form.
    rt = zz / e
    x = np.clip(r_end / rt, 0.0, 1.0)
    tail = np.sqrt(k2 * e) * rt * (np.arccos(np.sqrt(x)) - np.sqrt(x * (1 - x)))
    for ell in range(lmax + 1):
        # Outside R the potential rises to its top and falls beyond it, so the region where V > ε is the barrier.
        kappa = np.sqrt(np.maximum(v[ell][None, :] - e[:, None], 0.0) * k2)
        theta = kappa @ dr + np.where(v[ell][-1] > e, tail, 0.0)
        below = 1.0 / (1.0 + np.exp(np.minimum(2 * theta, 700.0)))
        above = 0.5 * (1 + np.tanh(np.pi * (eps - vb[ell]) / hw[ell]))
        t[ell] = np.where(eps >= vb[ell], above, below)
    t[:, eps <= 0] = 0.0
    return t


def gdr_parameters(z: int, a: int, model: EvaporationModel = EvaporationModel()):
    """Giant-dipole energy E_G (Berman–Fultz, MeV), width Γ_G (MeV) and peak σ₀ (Thomas–Reiche–Kuhn, mb)."""
    e_g = 31.2 * a ** (-1 / 3) + 20.6 * a ** (-1 / 6)
    gamma = model.gdr_width_mev
    sigma0 = 120.0 * z * (a - z) / (math.pi * a * gamma)  # ∫σ dE = 60 NZ/A mb MeV, Lorentzian ∫ = π σ₀ Γ / 2
    return e_g, gamma, sigma0


def e1_transmission(epsilon_mev, z: int, a: int, model: EvaporationModel = EvaporationModel()):
    """E1 transmission T(ε) = 2π ε³ f_E1(ε), f_E1 the Brink–Axel Lorentzian strength function (MeV⁻³)."""
    eps = np.asarray(epsilon_mev, dtype=float)
    e_g, gamma, sigma0 = gdr_parameters(z, a, model)
    f = 8.674e-8 * sigma0 * gamma**2 * eps / ((eps**2 - e_g**2) ** 2 + eps**2 * gamma**2)
    return 2 * np.pi * eps**3 * f


def fission_barrier_mev(z: int, a: int, model: EvaporationModel = EvaporationModel()) -> float:
    """Liquid-drop fission barrier at zero spin, MeV.

    Fissility x = (Z²/A) / (50.88 (1 − 1.7826 I²)), surface energy E_s = 17.9439 (1 − 1.7826 I²) A^⅔ MeV, I = (N−Z)/A
    (Myers and Swiatecki's constants), and the barrier B = 0.38 (0.75 − x) E_s for x ≤ 2/3, 0.83 (1 − x)³ E_s for
    2/3 < x < 1, a fit to the Cohen–Swiatecki saddle energies. No shell correction. ``fission_barrier_scale`` scales it.
    """
    i = (a - 2 * z) / a
    k = 1 - 1.7826 * i * i
    x = z * z / a / (50.88 * k)
    es = 17.9439 * k * a ** (2 / 3)
    if x >= 1:
        b = 0.0
    elif x > 2 / 3:
        b = 0.83 * (1 - x) ** 3 * es
    else:
        b = 0.38 * (0.75 - x) * es
    return b * model.fission_barrier_scale


# ---------------------------------------------------------------------------------------------------------------
# The Monte Carlo


def _name(z: int, a: int) -> str:
    return f"{a}{SYMBOLS[z]}"


def channel_label(n_n: int, n_p: int, n_a: int) -> str:
    """``"a p 2n"`` written ``"ap2n"``: α first, then protons, then neutrons; ``"g"`` when nothing was emitted."""
    def part(count, s):
        return "" if count == 0 else (s if count == 1 else f"{count}{s}")
    return (part(n_a, "a") + part(n_p, "p") + part(n_n, "n")) or "g"


class _Nucleus:
    """Level-density and fission data of one nuclide."""

    def __init__(self, z, a, model):
        self.z, self.a = z, a
        self.mass = nuclear_mass_mev(z, a)
        self.rho = _LevelDensity(z, a, model)
        self.lda, self.delta, self.rot, self.floor = self.rho.lda, self.rho.delta, self.rho.rot, self.rho.floor
        self.bf = fission_barrier_mev(z, a, model)
        self.rot_s = rotational_constant_mev(a, model, model.saddle_inertia_ratio)
        self.lda_s = self.lda * model.af_over_an
        self._g = None

    def fermi_m(self, e, m):
        """ρ_M(E, M) for excitation e and projection m (arrays broadcast)."""
        return self.rho.projection(e, m)

    def fission_width(self, e, j):
        """∫ ρ_saddle(E − B_f − ε, J) dε (Bohr–Wheeler, without the common 1/2πρ_CN)."""
        if self._g is None:
            umax = 400.0
            u = np.linspace(0.0, umax, 20001)
            f = (u + self.floor) ** -1.5 * np.exp(2 * np.sqrt(self.lda_s * u))
            g = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(u))])
            self._g = (u, g * math.sqrt(self.rot_s) / 12.0)
        u, g = self._g
        x0 = e - self.bf - self.delta - self.rot_s * j**2
        x1 = e - self.bf - self.delta - self.rot_s * (j + 1.0) ** 2
        return np.interp(x0, u, g, left=0.0) - np.interp(x1, u, g, left=0.0)


class _Simulation:
    """Tables and the random walk of one call to :func:`evaporate`."""

    def __init__(self, z, a, e_max, j_max, decays, model, rng):
        self.model, self.rng, self.decays = model, rng, decays
        self.cn = (z, a)
        self.e_max = e_max
        self.step = model.energy_step_mev
        self.nk = int(round(model.max_emission_mev / self.step)) + 1  # emission-energy nodes 0, Δ, 2Δ, ...
        self.eps_nodes = np.arange(self.nk) * self.step
        self.nuclei: dict = {}
        self.parents: list = []
        self.parent_of: dict = {}
        self.row_len = 4 * (self.nk - 1) + 1
        self.rows = np.zeros((1024, self.row_len))
        self.nrows = 0
        self.particle_mass = {k: nuclear_mass_mev(*PARTICLES[k]) for k in PARTICLES}
        self._tcache: dict = {}
        self.j_cap = j_max

    # --- nuclei and their tables -------------------------------------------------------------------------------

    def nucleus(self, z, a):
        key = (z, a)
        if key not in self.nuclei:
            self.nuclei[key] = _Nucleus(z, a, self.model) if z >= 0 and a > z and nuclear_mass_mev(z, a) else None
        return self.nuclei[key]

    def transmission(self, kind, zd, ad, eps):
        """T[ℓ, ε] for emission of ``kind`` leaving the daughter (zd, ad), and the ℓ range kept."""
        key = (kind, zd, ad)
        if key not in self._tcache:
            if kind == "n":
                t = neutron_transmission(self.eps_nodes, ad, 60, self.model)
                params = None
            else:
                vb, _, _ = barriers(kind, zd, ad, 120)
                lmax = int(np.nonzero(vb < self.eps_nodes[-1] + 15.0)[0].max())
                t = charged_transmission(kind, zd, ad, self.eps_nodes, lmax)
                params = np.log(np.maximum(t, 1e-300))
            keep = np.nonzero(t.max(axis=1) > 1e-7)[0]
            lmax = int(keep.max()) if keep.size else 0
            self._tcache[key] = (t[: lmax + 1], params, lmax)
        t, params, lmax = self._tcache[key]
        if eps is None:
            return t, lmax
        if kind == "n":
            return neutron_transmission(eps, ad, lmax, self.model)
        # Log-linear interpolation between the nodes.
        x = np.clip(eps / self.step, 0.0, self.nk - 1 - 1e-9)
        k = x.astype(np.int64)
        w = x - k
        logt = params[: lmax + 1]
        return np.exp((1 - w) * logt[:, k] + w * logt[:, np.minimum(k + 1, self.nk - 1)])

    def parent(self, z, a):
        key = (z, a)
        if key in self.parent_of:
            return self.parents[self.parent_of[key]]
        nuc = self.nucleus(z, a)
        zc, ac = self.cn
        dz, dn = zc - z, (ac - a) - (zc - z)
        cn = self.nucleus(zc, ac)
        # Highest excitation this nucleus can have: emitting the missing nucleons in the cheapest way.
        mp, mn = self.particle_mass["p"], self.particle_mass["n"]
        emax = self.e_max + cn.mass - nuc.mass - (dz * mp + dn * mn) + _ALPHA_BINDING * (min(dz, dn) // 2) + 1.0
        emax = min(max(emax, 0.0), self.e_max + 1.0)
        ni = int(emax / self.step) + 2
        jg = min(int(math.sqrt(max(emax, 0.0) / nuc.rot)) + 2, self.j_cap)
        p = {"nuc": nuc, "z": z, "a": a, "ni": ni, "jg": jg, "map": np.full((ni, jg + 1), -1, np.int64),
             "ch": []}
        for kind in ("n", "p", "a"):
            zp, ap = PARTICLES[kind]
            d = self.nucleus(z - zp, a - ap) if kind in self.decays else None
            if d is None or (kind != "n" and z - zp < 1):
                p["ch"].append(None)
                continue
            sep = d.mass + self.particle_mass[kind] - nuc.mass
            t, lmax = self.transmission(kind, d.z, d.a, None)
            ms = np.arange(jg + lmax + 3)
            m_idx = np.arange(ni + self.nk) - self.nk  # P row r holds parent bin m = r − nk (negative → empty)
            e_d = m_idx[:, None] * self.step - sep
            ptab = d.fermi_m(e_d, ms[None, :])
            g = 2.0 if kind in ("n", "p") else 1.0  # 2s + 1
            p["ch"].append({"kind": kind, "d": d, "sep": 0.0 + sep, "t": t, "ells": np.arange(lmax + 1),
                            "P": ptab, "g": g, "mass": self.particle_mass[kind]})
        if "g" in self.decays:
            ells = np.array([1])
            t = e1_transmission(self.eps_nodes, z, a, self.model)[None, :]
            ms = np.arange(jg + 5)
            m_idx = np.arange(ni + self.nk) - self.nk
            ptab = nuc.fermi_m(m_idx[:, None] * self.step, ms[None, :])
            p["ch"].append({"kind": "g", "d": nuc, "sep": 0.0, "t": t, "ells": ells, "P": ptab, "g": 1.0,
                            "mass": 0.0})
        else:
            p["ch"].append(None)
        p["fission"] = "fission" in self.decays and nuc.bf < emax
        self.parent_of[key] = len(self.parents)
        self.parents.append(p)
        return p

    def _grow(self, n):
        while self.nrows + n > self.rows.shape[0]:
            self.rows = np.concatenate([self.rows, np.zeros_like(self.rows)])

    def build_rows(self, p, i, j):
        """Cumulative decay weights for the states (bin i, spin j) of parent p; returns their row numbers."""
        nk, step = self.nk, self.step
        out = np.empty(i.size, np.int64)
        chunk = max(1, 200_000 // (nk * 40))
        k = np.arange(nk)
        for s0 in range(0, i.size, chunk):
            ii, jj = i[s0:s0 + chunk], j[s0:s0 + chunk]
            n = ii.size
            w = np.zeros((n, self.row_len))
            mrow = (ii[:, None] - k[None, :]) + nk  # row of P for each node
            for c, ch in enumerate(p["ch"]):
                if ch is None:
                    continue
                ell = ch["ells"]
                lo = np.abs(jj[:, None] - ell[None, :])  # (n, L)
                hi = jj[:, None] + ell[None, :] + 1
                ptab = ch["P"]
                coupled = ptab[mrow[:, :, None], lo[:, None, :]] - ptab[mrow[:, :, None], hi[:, None, :]]
                f = ch["g"] * np.einsum("nkl,lk->nk", coupled, ch["t"])
                w[:, c * (nk - 1):(c + 1) * (nk - 1)] = 0.5 * step * (f[:, 1:] + f[:, :-1])
            if p["fission"]:
                w[:, -1] = np.maximum(p["nuc"].fission_width(ii * step, jj.astype(float)), 0.0)
            np.maximum(w, 0.0, out=w)
            self._grow(n)
            self.rows[self.nrows:self.nrows + n] = np.cumsum(w, axis=1)
            out[s0:s0 + n] = np.arange(self.nrows, self.nrows + n)
            self.nrows += n
        return out

    def rows_for(self, pid, e, j):
        p = self.parents[pid]
        i = np.clip((e / self.step).astype(np.int64), 0, p["ni"] - 1)
        jj = np.minimum(j, p["jg"])
        r = p["map"][i, jj]
        new = r < 0
        if new.any():
            keys = i[new] * (p["jg"] + 1) + jj[new]
            uk, inv = np.unique(keys, return_inverse=True)
            built = self.build_rows(p, uk // (p["jg"] + 1), uk % (p["jg"] + 1))
            p["map"][uk // (p["jg"] + 1), uk % (p["jg"] + 1)] = built
            r[new] = built[inv.ravel()]
        return r

    # --- the walk ----------------------------------------------------------------------------------------------

    def run(self, e0, j0):
        rng, nk = self.rng, self.nk
        n = e0.size
        zc, ac = self.cn
        root = self.parent(zc, ac)
        cn_mass = root["nuc"].mass
        pid = np.zeros(n, np.int64)
        e = e0.astype(float).copy()
        j = j0.astype(np.int64).copy()
        mom = np.zeros((n, 4))
        mom[:, 0] = cn_mass + e
        done = np.zeros(n, bool)
        fission = np.zeros(n, bool)
        counts = np.zeros((n, 3), np.int64)
        # State after the last particle (the entry state).
        snap_pid, snap_e, snap_j, snap_mom = pid.copy(), e.copy(), j.copy(), mom.copy()
        snap_step = np.full(n, -1, np.int64)
        emissions = []  # (event, step, kind, 4-momentum)
        steps = 0
        while not done.all():
            if steps >= self.model.max_steps:
                break
            act = np.nonzero(~done)[0]
            rows = np.empty(act.size, np.int64)
            for q in np.unique(pid[act]):
                sel = pid[act] == q
                rows[sel] = self.rows_for(q, e[act[sel]], j[act[sel]])
            part_total = self.rows[rows, 3 * (nk - 1) - 1]
            stop = part_total <= self.model.stop_ratio * self.rows[rows, -1]
            if "g" in self.decays:
                # No statistical γ ray possible: the state is on the yrast line and decays by discrete γ rays.
                stop |= self.rows[rows, 4 * (nk - 1) - 1] <= part_total
            done[act[stop]] = True
            go = ~stop
            act, rows = act[go], rows[go]
            u = rng.random(act.size) * self.rows[rows, -1]
            # Binary search for the first cumulative weight ≥ u in each event's row.
            lo = np.zeros(act.size, np.int64)
            hi = np.full(act.size, self.row_len - 1, np.int64)
            while True:
                open_ = lo < hi
                if not open_.any():
                    break
                mid = (lo + hi) // 2
                below = self.rows[rows, mid] < u
                lo = np.where(open_ & below, mid + 1, lo)
                hi = np.where(open_ & ~below, mid, hi)
            idx = lo
            channel = np.minimum(idx // (nk - 1), 4)
            kbin = idx % (nk - 1)
            fis = channel == 4
            fission[act[fis]] = True
            done[act[fis]] = True
            frac = rng.random(act.size)
            iso = rng.random((act.size, 2))
            lsel = rng.random(act.size)
            jsel = rng.random(act.size)
            for c in range(4):
                sel = np.nonzero(channel == c)[0]
                if sel.size == 0:
                    continue
                ev = act[sel]
                groups = [(q, sel[pid[ev] == q]) for q in np.unique(pid[ev])]  # before any event moves
                for q, sub in groups:
                    self._emit(c, q, act[sub], kbin[sub], frac[sub], iso[sub], lsel[sub], jsel[sub],
                               pid, e, j, mom, steps, emissions)
                if c < 3:
                    counts[ev, c] += 1
                    snap_pid[ev], snap_e[ev], snap_j[ev] = pid[ev], e[ev], j[ev]
                    snap_mom[ev] = mom[ev]
                    snap_step[ev] = steps
            steps += 1
        truncated = int((~done).sum())
        # Fission: the nucleus that fissioned is the "residue"; otherwise the entry state.
        keep = ~fission
        for arr, snap in ((pid, snap_pid), (e, snap_e), (j, snap_j)):
            arr[keep] = snap[keep]
        mom[keep] = snap_mom[keep]
        last = np.where(fission, steps + 1, snap_step)
        if emissions:
            ev = np.concatenate([x[0] for x in emissions])
            st = np.concatenate([x[1] for x in emissions])
            kd = np.concatenate([x[2] for x in emissions])
            p4 = np.concatenate([x[3] for x in emissions])
            ok = st <= last[ev]
            ev, st, kd, p4 = ev[ok], st[ok], kd[ok], p4[ok]
            order = np.lexsort((st, ev))
            ev, kd, p4 = ev[order], kd[order], p4[order]
        else:
            ev, kd, p4 = np.zeros(0, np.int64), np.zeros(0, np.int8), np.zeros((0, 4))
        zs = np.array([p["z"] for p in self.parents])[pid]
        as_ = np.array([p["a"] for p in self.parents])[pid]
        return dict(z=zs, a=as_, e=e, j=j, mom=mom, fission=fission, counts=counts, ev=ev, kind=kd, p4=p4,
                    truncated=truncated)

    def _emit(self, c, q, ev, kbin, frac, iso, lsel, jsel, pid, e, j, mom, steps, emissions):
        """Emit from the events ``ev`` (all in parent q) by channel c (0 n, 1 p, 2 α, 3 γ)."""
        p = self.parents[q]
        nuc = p["nuc"]
        step = self.step
        eps = (kbin + frac) * step
        jj = j[ev]
        ch = p["ch"][c]
        d, sep, m1, ells = ch["d"], ch["sep"], ch["mass"], ch["ells"]
        lo = np.abs(jj[None, :] - ells[:, None])
        hi = jj[None, :] + ells[:, None] + 1

        def weights(eps, sel):
            ed = e[ev[sel]] - sep - eps
            if c == 3:
                tl = e1_transmission(eps, nuc.z, nuc.a, self.model)[None, :]
            else:
                tl = self.transmission(ch["kind"], d.z, d.a, eps)
            return tl * (d.fermi_m(ed[None, :], lo[:, sel]) - d.fermi_m(ed[None, :], hi[:, sel]))

        # The energy inside the chosen grid bin: rejection against the exact weight, under an envelope of 1.25 ×
        # the larger of the two edge values (the weight is smooth on the scale of one bin).
        lo_e, hi_e = kbin * step, (kbin + 1) * step
        env = 1.25 * np.maximum(weights(lo_e, slice(None)).sum(axis=0), weights(hi_e, slice(None)).sum(axis=0))
        wl = np.zeros((ells.size, ev.size))
        pend = np.nonzero(env > 0)[0]
        for _ in range(200):
            if pend.size == 0:
                break
            x = lo_e[pend] + self.rng.random(pend.size) * step
            wx = weights(x, pend)
            acc = self.rng.random(pend.size) * env[pend] <= wx.sum(axis=0)
            eps[pend[acc]] = x[acc]
            wl[:, pend[acc]] = wx[:, acc]
            pend = pend[~acc]
        # Rare: the state at the bin's energies falls under the yrast line; take the bin's low edge then.
        bad = wl.sum(axis=0) <= 0
        if bad.any():
            eps[bad] = lo_e[bad]
            wl[:, bad] = weights(eps[bad], bad)
        ed = e[ev] - sep - eps
        wl = np.maximum(wl, 0.0)
        cw = np.cumsum(wl, axis=0)
        ell = ells[np.minimum((cw < lsel[None, :] * cw[-1][None, :]).sum(axis=0), ells.size - 1)]
        # Daughter spin J' in |J−ℓ| .. J+ℓ, weighted by ρ(E', J').
        jlo = np.abs(jj - ell)
        width = 2 * int(ell.max()) + 1
        cand = jlo[None, :] + np.arange(width)[:, None]
        valid = cand <= (jj + ell)[None, :]
        wj = np.where(valid, d.fermi_m(ed[None, :], cand) - d.fermi_m(ed[None, :], cand + 1), 0.0)
        wj = np.maximum(wj, 0.0)
        cj = np.cumsum(wj, axis=0)
        jn = jlo + np.minimum((cj < jsel[None, :] * cj[-1][None, :]).sum(axis=0), width - 1)
        # Relativistic two-body decay in the emitter's frame, isotropic, then boosted to the CN frame.
        m2 = d.mass + ed
        mstar = m1 + m2 + eps
        pstar = np.sqrt(eps * (eps + 2 * m1) * (eps + 2 * m2) * (eps + 2 * m1 + 2 * m2)) / (2 * mstar)
        cth = 2 * iso[:, 0] - 1
        sth = np.sqrt(np.maximum(1 - cth**2, 0.0))
        phi = 2 * np.pi * iso[:, 1]
        dirs = np.stack([sth * np.cos(phi), sth * np.sin(phi), cth], axis=1)
        p1 = np.empty((ev.size, 4))
        p1[:, 0] = np.sqrt(pstar**2 + m1**2)
        p1[:, 1:] = pstar[:, None] * dirs
        parent4 = mom[ev]
        p1 = _boost(p1, parent4)
        mom[ev] = parent4 - p1
        # Keep the daughter exactly on its mass shell (rounding only).
        e[ev] = ed
        j[ev] = jn
        if c != 3:
            pid[ev] = self._pid(d.z, d.a)
        emissions.append((ev.copy(), np.full(ev.size, steps, np.int64), np.full(ev.size, c, np.int8), p1))

    def _pid(self, z, a):
        self.parent(z, a)
        return self.parent_of[(z, a)]


def _boost(p, frame):
    """Boost four-vectors p (n×4, E first) from the rest frame of ``frame`` (n×4) to the frame ``frame`` is in."""
    m = np.sqrt(np.maximum(frame[:, 0] ** 2 - (frame[:, 1:] ** 2).sum(axis=1), 0.0))
    b = frame[:, 1:] / frame[:, :1]
    g = frame[:, 0] / m
    bp = (b * p[:, 1:]).sum(axis=1)
    out = np.empty_like(p)
    out[:, 0] = g * (p[:, 0] + bp)
    out[:, 1:] = p[:, 1:] + b * ((g * g / (g + 1) * bp) + g * p[:, 0])[:, None]
    return out


@dataclass
class Evaporation:
    """Result of :func:`evaporate`. Energies MeV, momenta MeV/c, spins ħ, cross sections mb."""

    compound: tuple
    excitation_mev: float
    sigma_fusion_mb: Optional[float]
    model: EvaporationModel
    initial_spin: np.ndarray
    residue_z: np.ndarray
    residue_a: np.ndarray
    channel: np.ndarray
    labels: list
    fission: np.ndarray
    entry_energy_mev: np.ndarray
    entry_spin: np.ndarray
    recoil_momentum_mev_c: np.ndarray
    offsets: np.ndarray
    event: np.ndarray
    kind: np.ndarray
    energy_mev: np.ndarray
    momentum_mev_c: np.ndarray
    counts: np.ndarray = field(repr=False)
    truncated: int = 0
    spin_clamped: int = 0

    @property
    def events(self) -> int:
        return int(self.initial_spin.size)

    @property
    def direction(self) -> np.ndarray:
        """Unit vectors of the emissions, in the compound nucleus's frame."""
        norm = np.linalg.norm(self.momentum_mev_c, axis=1)
        return self.momentum_mev_c / np.where(norm > 0, norm, 1.0)[:, None]

    @property
    def label_of_event(self) -> np.ndarray:
        return np.array(self.labels, dtype=object)[self.channel]

    def particles(self, i: int) -> dict:
        """The emissions of event i: kinds (as letters), energies (MeV) and momenta (MeV/c)."""
        s = slice(self.offsets[i], self.offsets[i + 1])
        return {"kind": [KINDS[k] for k in self.kind[s]], "energy_mev": self.energy_mev[s],
                "momentum_mev_c": self.momentum_mev_c[s]}

    @property
    def mass_mev(self) -> np.ndarray:
        """Rest mass of each emission (0 for a γ ray), MeV."""
        masses = np.array([nuclear_mass_mev(*PARTICLES[k]) for k in KINDS[:3]] + [0.0])
        return masses[self.kind]

    def lab(self, beta: float) -> tuple:
        """Kinetic energies (MeV) and momenta (MeV/c, emissions × 3) in the lab, for a compound nucleus moving with
        velocity β (units of c) along +z, the beam direction."""
        g = 1.0 / math.sqrt(1.0 - beta * beta)
        m = self.mass_mev
        e = self.energy_mev + m
        pz = self.momentum_mev_c[:, 2]
        mom = self.momentum_mev_c.copy()
        mom[:, 2] = g * (pz + beta * e)
        return g * (e + beta * pz) - m, mom

    def events_list(self, beta: float = 0.0) -> list:
        """The events as ``{"channel": label, "particles": [(kind, energy_mev, direction), ...]}``, energies and
        directions in the lab for a compound nucleus moving with β along +z (β = 0: its own frame)."""
        energy, mom = self.lab(beta)
        norm = np.linalg.norm(mom, axis=1)
        dirs = mom / np.where(norm > 0, norm, 1.0)[:, None]
        labels = self.label_of_event
        out = []
        for i in range(self.events):
            s = slice(self.offsets[i], self.offsets[i + 1])
            out.append({"channel": labels[i],
                        "particles": [(KINDS[k], float(en), d) for k, en, d in zip(self.kind[s], energy[s], dirs[s])]})
        return out

    def shares(self) -> dict:
        """Share of events in each channel (labels as keys), summing to one."""
        counts = np.bincount(self.channel, minlength=len(self.labels))
        return {lab: counts[c] / self.events for c, lab in enumerate(self.labels) if counts[c]}

    def channel_table(self) -> list:
        """One row per exit channel, strongest first: label, residue, Z, A, events, share, its statistical error,
        and cross section σ_fus × share in mb (None when σ_fus is unknown)."""
        counts = np.bincount(self.channel, minlength=len(self.labels))
        rows = []
        for c in np.argsort(-counts, kind="stable"):
            if not counts[c]:
                continue
            first = int(np.nonzero(self.channel == c)[0][0])
            z, a = int(self.residue_z[first]), int(self.residue_a[first])
            lab = self.labels[c]
            share = counts[c] / self.events
            rows.append({
                "label": lab,
                "residue": "fission" if lab == "fission" else _name(z, a),
                "Z": None if lab == "fission" else z, "A": None if lab == "fission" else a,
                "events": int(counts[c]), "share": share, "share_error": math.sqrt(counts[c]) / self.events,
                "sigma_mb": None if self.sigma_fusion_mb is None else share * self.sigma_fusion_mb,
            })
        return rows

    def mean_energy_mev(self, kind: str) -> float:
        """Mean kinetic energy of the emitted ``kind`` ("n", "p", "a", "g") in the compound nucleus's frame, MeV."""
        sel = self.kind == KINDS.index(kind)
        return float(self.energy_mev[sel].mean()) if sel.any() else float("nan")

    def multiplicity(self, kind: str) -> float:
        """Mean number of ``kind`` per event (γ: only those between particles)."""
        return float((self.kind == KINDS.index(kind)).sum() / self.events)

    def summary(self) -> dict:
        return {
            "compound": _name(*self.compound), "excitation_mev": self.excitation_mev,
            "sigma_fusion_mb": self.sigma_fusion_mb, "events": self.events,
            "channels": self.channel_table(),
            "mean_energy_mev": {k: self.mean_energy_mev(k) for k in KINDS[:3]},
            "multiplicity": {k: self.multiplicity(k) for k in KINDS[:3]},
            "fission_share": float(self.fission.mean()),
        }


def _spins(n, rng, sigma_l_mb, spin, spin_weights):
    given = [x is not None for x in (sigma_l_mb, spin, spin_weights)]
    if sum(given) != 1:
        raise ValueError("give exactly one of sigma_l_mb, spin, spin_weights")
    if spin is not None:
        if int(spin) != spin or spin < 0:
            raise ValueError("spin must be a whole number of ħ, at least 0")
        return np.full(n, int(spin), np.int64), None
    w = np.asarray(sigma_l_mb if sigma_l_mb is not None else spin_weights, dtype=float)
    if w.ndim != 1 or w.size == 0 or (w < 0).any() or w.sum() <= 0:
        raise ValueError("the spin distribution must be a 1-D array of non-negative numbers, not all zero")
    cdf = np.cumsum(w) / w.sum()
    j = np.minimum(np.searchsorted(cdf, rng.random(n), side="right"), w.size - 1)
    return j.astype(np.int64), (float(w.sum()) if sigma_l_mb is not None else None)


def evaporate(compound: NuclideLike, excitation_mev: float, *, sigma_l_mb: Optional[Sequence[float]] = None,
              spin: Optional[int] = None, spin_weights: Optional[Sequence[float]] = None, events: int = 10000,
              seed=None, decays: Sequence[str] = DECAYS, model: Optional[EvaporationModel] = None) -> Evaporation:
    """Run the evaporation Monte Carlo; see the module docstring for the inputs and outputs."""
    model = model or EvaporationModel()
    if isinstance(compound, tuple):
        z, a = compound
    else:
        nuc = data.nuclide(compound)
        z, a = nuc.Z, nuc.A
    if nuclear_mass_mev(z, a) is None:
        raise ValueError(f"{_name(z, a)} is not in the AME2020 mass table")
    bad = set(decays) - set(DECAYS)
    if bad:
        raise ValueError(f"unknown decays {sorted(bad)}; choose from {DECAYS}")
    if excitation_mev <= 0:
        raise ValueError("the excitation energy must be positive")
    events = int(events)
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    j0, sigma_fus = _spins(events, rng, sigma_l_mb, spin, spin_weights)
    rot = rotational_constant_mev(a, model)
    delta = pairing_shift_mev(z, a, model)
    j_yrast = int(math.floor(math.sqrt(max(excitation_mev - delta, 0.0) / rot)))
    clamped = int((j0 > j_yrast).sum())
    j0 = np.minimum(j0, max(j_yrast - 1, 0))
    sim = _Simulation(z, a, float(excitation_mev), j_yrast + 2, tuple(decays), model, rng)
    out = sim.run(np.full(events, float(excitation_mev)), j0)
    cnt = out["counts"]
    code = np.where(out["fission"], -1, cnt[:, 0] + 1000 * cnt[:, 1] + 1_000_000 * cnt[:, 2])
    ucode, channel = np.unique(code, return_inverse=True)
    channel = channel.ravel().astype(np.int64)
    labels = ["fission" if c < 0 else channel_label(c % 1000, (c // 1000) % 1000, c // 1_000_000) for c in ucode]
    ev, kd, p4 = out["ev"], out["kind"], out["p4"]
    masses = np.array([sim.particle_mass["n"], sim.particle_mass["p"], sim.particle_mass["a"], 0.0])
    offsets = np.concatenate([[0], np.cumsum(np.bincount(ev, minlength=events))]).astype(np.int64)
    return Evaporation(
        compound=(z, a), excitation_mev=float(excitation_mev), sigma_fusion_mb=sigma_fus, model=model,
        initial_spin=j0, residue_z=out["z"].astype(np.int64), residue_a=out["a"].astype(np.int64),
        channel=channel, labels=labels, fission=out["fission"], entry_energy_mev=out["e"],
        entry_spin=out["j"], recoil_momentum_mev_c=out["mom"][:, 1:].copy(), offsets=offsets, event=ev,
        kind=kd.astype(np.int8), energy_mev=p4[:, 0] - masses[kd], momentum_mev_c=p4[:, 1:].copy(),
        counts=cnt, truncated=out["truncated"], spin_clamped=clamped,
    )
