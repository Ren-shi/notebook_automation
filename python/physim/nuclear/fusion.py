"""Fusion of two nuclei: the Bass barrier, Wong's capture cross section and its partial waves, the compound
nucleus they make, and the excitation function averaged over a target's thickness.

::

    from physim.nuclear.fusion import Fusion

    f = Fusion("16O", "208Pb")
    f.barrier()                           # Bass barrier for l = 0: height (MeV, CM), radius (fm), ħω (MeV)
    f.cross_section("90 MeV")             # σ_fus at a lab beam energy, mb
    f.cross_section_cm(80.0)              # at a CM kinetic energy, mb
    f.wong_cm(80.0)                       # Wong's closed form with the l = 0 barrier, mb
    cn = f.compound("90 MeV")             # the compound nucleus: E*, σ_l by l, spin distribution
    f.thick_target("90 MeV", "1 mg/cm2")  # σ averaged over the beam's energy loss in the target
    f.warnings("90 MeV")                  # plain-language validity warnings

Energies are lab kinetic energies of the beam in MeV (a float, a ``Quantity`` or text such as ``"90 MeV"``) unless
a name says cm; CM energies are floats in MeV. Cross sections are in mb, lengths in fm, angular momenta in units
of ħ. Arrays of energies work wherever one energy does, except in :meth:`Fusion.compound` and
:meth:`Fusion.thick_target`, which describe one beam energy.

The model: the nucleus–nucleus potential is Coulomb (point charges) plus the Bass nuclear potential, plus the
centrifugal term ħ²l(l+1)/(2μr²). For each partial wave l the barrier of that potential is replaced by a parabola
of the same height V_l, position R_l and curvature ħω_l, whose transmission is the Hill–Wheeler
T_l = 1/(1 + exp[2π(V_l − E)/ħω_l]). The capture cross section is σ_fus = πƛ² Σ (2l + 1) T_l with ƛ² = ħ²/(2μE).
Partial waves whose potential has no pocket (above the critical angular momentum) do not fuse. Capture is taken
as fusion: no quasi-fission, no incomplete fusion. Details and limits are in the theory note ``theory/fusion``.

Interface for the evaporation code
----------------------------------
:meth:`Fusion.compound` returns a :class:`CompoundNucleus` for one beam energy:

- ``nuclide``: the compound nucleus (Z₁ + Z₂, A₁ + A₂), a :class:`physim.nuclear.data.Nuclide`;
- ``excitation_mev``: its excitation energy E* = √s − M_CN, MeV, exact relativistic kinematics with AME2020 nuclear
  masses (equal to E_cm + Q to well below 1 keV at these energies);
- ``e_cm_mev``: the CM kinetic energy, MeV;
- ``beta``: the compound nucleus's velocity along the beam, units of c;
- ``sigma_l_mb``: a 1-D ``numpy`` array of σ_l in mb indexed by l = 0, 1, 2, ... (its length is the highest l
  that contributes, plus one);
- ``sigma_mb``: their sum, mb; ``spin_distribution()``: σ_l / σ_fus, the same indexing, summing to 1.

With spin-0 beam and target ground states the compound nucleus's spin J equals l. For nuclei with ground-state
spin the array is the orbital part only: this module does not couple the ground-state spins.
:meth:`Fusion.thick_target` returns the same σ_l averaged over the target's depth, and E* at each depth.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Union

import numpy as np

from . import data
from .quantity import Quantity
from .rutherford import E2_MEV_FM

ArrayLike = Union[float, np.ndarray]
EnergyLike = Union[float, str, Quantity, np.ndarray]

#: ħc, MeV fm (CODATA 2018).
HBARC_MEV_FM = 197.3269804

#: Parameters of the Bass universal function g(s) = [A exp(s/d₁) + B exp(s/d₂)]⁻¹ (MeV/fm): A and B in fm/MeV,
#: d₁ and d₂ in fm. "bass80": R. Bass, *Nuclear Reactions with Heavy Ions* (Springer 1980); "bass77": R. Bass,
#: Phys. Rev. Lett. 39, 265 (1977). Both as written out by I. Dutt and R. K. Puri, Phys. Rev. C 81, 064609 (2010),
#: whose Table II also lists Bass 1980 barriers that the tests compare against.
BASS_PARAMETERS = {
    "bass80": {"A": 0.033, "B": 0.007, "d1": 3.5, "d2": 0.65},
    "bass77": {"A": 0.0300, "B": 0.0061, "d1": 3.30, "d2": 0.65},
}

#: Above this many MeV per nucleon of beam, incomplete fusion (break-up of the projectile) takes a sizeable share.
INCOMPLETE_FUSION_MEV_PER_U = 10.0
#: Above this Z₁Z₂, quasi-fission competes with fusion and capture overestimates fusion (an empirical threshold, see
#: K.-H. Schmidt and W. Morawek, Rep. Prog. Phys. 54, 949 (1991), and W. J. Swiatecki's extra-push model,
#: Nucl. Phys. A 376, 275 (1982)).
QUASI_FISSION_Z1Z2 = 1600
#: Highest partial wave searched for a barrier; far above any critical angular momentum of the Bass potential.
L_LIMIT = 400


def _energy_mev(value: EnergyLike) -> ArrayLike:
    if isinstance(value, str):
        value = Quantity.parse(value)
    if isinstance(value, Quantity):
        if value.kind != "energy":
            raise ValueError(f"beam energy '{value}' must be a total energy such as '90 MeV'")
        return value.to("MeV")
    if np.ndim(value):
        return np.asarray(value, dtype=float)
    return float(value)


def _out(x, like):
    return float(x) if np.ndim(like) == 0 else np.asarray(x, dtype=float)


def _transmission(v: ArrayLike, hw: ArrayLike, e: ArrayLike) -> np.ndarray:
    """Hill–Wheeler transmission 1/(1 + exp[2π(V − E)/ħω]) of a parabolic barrier, without overflow."""
    x = 2 * np.pi * (np.asarray(v) - np.asarray(e)) / np.asarray(hw)
    return 0.5 * (1.0 - np.tanh(0.5 * x))


@dataclass(frozen=True)
class Barrier:
    """A fusion barrier for one partial wave, in the CM frame."""

    #: Partial wave (units of ħ).
    l: int
    #: Height V_l of the barrier, MeV.
    height_mev: float
    #: Position R_l of the barrier's top, fm.
    radius_fm: float
    #: Curvature ħω_l = ħ √(−V''(R_l)/μ), MeV.
    hbar_omega_mev: float


@dataclass
class CompoundNucleus:
    """The compound nucleus made at one beam energy, with its spin distribution (see the module docstring)."""

    nuclide: data.Nuclide
    #: Lab kinetic energy of the beam, MeV.
    beam_energy_mev: float
    #: CM kinetic energy, MeV.
    e_cm_mev: float
    #: Excitation energy E* = √s − M_CN, MeV.
    excitation_mev: float
    #: Velocity of the compound nucleus along the beam, units of c.
    beta: float
    #: σ_l in mb, indexed by l.
    sigma_l_mb: np.ndarray

    @property
    def sigma_mb(self) -> float:
        """σ_fus = Σ σ_l, mb."""
        return float(self.sigma_l_mb.sum())

    def spin_distribution(self) -> np.ndarray:
        """σ_l / σ_fus, indexed by l (all zeros when nothing fuses)."""
        total = self.sigma_mb
        return self.sigma_l_mb / total if total > 0 else np.zeros_like(self.sigma_l_mb)

    @property
    def mean_l(self) -> float:
        """⟨l⟩ weighted by σ_l, ħ."""
        p = self.spin_distribution()
        return float(np.dot(np.arange(p.size), p))

    @property
    def max_l(self) -> int:
        """Highest l with σ_l above 10⁻⁶ of the largest σ_l (0 when nothing fuses)."""
        s = self.sigma_l_mb
        if s.size == 0 or s.max() <= 0:
            return 0
        return int(np.nonzero(s > 1e-6 * s.max())[0][-1])


@dataclass
class ThickTarget:
    """Fusion averaged over the beam's path through a target, σ weighted by depth."""

    #: Lab beam energy entering and leaving the target, MeV (0 if the beam stops).
    e_in_mev: float
    e_out_mev: float
    #: Target thickness, mg/cm².
    thickness_mg_cm2: float
    #: Depth-averaged σ_fus = (1/t) ∫₀ᵗ σ(E(x)) dx, mb.
    sigma_mb: float
    #: Depth-averaged σ_l, mb, indexed by l.
    sigma_l_mb: np.ndarray
    #: Depths x (mg/cm²), lab beam energy E(x) (MeV), CM energy, E* (MeV) and σ_fus (mb) on the integration grid.
    depth_mg_cm2: np.ndarray = field(repr=False)
    beam_energy_mev: np.ndarray = field(repr=False)
    e_cm_mev: np.ndarray = field(repr=False)
    excitation_mev: np.ndarray = field(repr=False)
    sigma_depth_mb: np.ndarray = field(repr=False)


class Fusion:
    """Fusion of ``beam`` on ``target`` (nuclide names, e.g. ``"16O"``, ``"208Pb"``) through the Bass potential.

    ``potential`` picks the Bass parameter set: ``"bass80"`` (the default) or ``"bass77"``.
    """

    def __init__(self, beam: str, target: str, potential: str = "bass80"):
        if potential not in BASS_PARAMETERS:
            raise ValueError(f"potential must be one of {', '.join(BASS_PARAMETERS)}, not '{potential}'")
        self.beam, self.target = data.nuclide(beam), data.nuclide(target)
        if self.beam.Z < 1 or self.target.Z < 1:
            raise ValueError("fusion needs two charged nuclei")
        self.potential = potential
        self.compound_nuclide = data.nuclide((self.beam.Z + self.target.Z, self.beam.A + self.target.A))
        self.m1, self.m2 = self.beam.nuclear_mass_mev, self.target.nuclear_mass_mev
        #: Reduced mass, MeV.
        self.mu = self.m1 * self.m2 / (self.m1 + self.m2)
        #: Z₁Z₂e², MeV fm.
        self.k = self.beam.Z * self.target.Z * E2_MEV_FM
        #: Radii R₁, R₂ of the Bass potential, fm.
        self.r1, self.r2 = self._radius(self.beam.A), self._radius(self.target.A)
        self._p = BASS_PARAMETERS[potential]
        self._barriers: Optional[list] = None

    # -- the potential ------------------------------------------------------------------------------------------

    def _radius(self, a: int) -> float:
        if self.potential == "bass77":
            return 1.16 * a ** (1 / 3) - 1.39 * a ** (-1 / 3)
        rs = 1.28 * a ** (1 / 3) - 0.76 + 0.8 * a ** (-1 / 3)
        return rs * (1 - 0.98 / rs**2)

    def _nuclear(self, r: np.ndarray, order: int = 0) -> np.ndarray:
        """Bass nuclear potential −R̄ g(s) (MeV) or its first or second derivative in r."""
        p = self._p
        s = r - self.r1 - self.r2
        e1, e2 = p["A"] * np.exp(s / p["d1"]), p["B"] * np.exp(s / p["d2"])
        d = e1 + e2
        rbar = self.r1 * self.r2 / (self.r1 + self.r2)
        if order == 0:
            return -rbar / d
        d1 = e1 / p["d1"] + e2 / p["d2"]
        if order == 1:
            return rbar * d1 / d**2
        d2 = e1 / p["d1"] ** 2 + e2 / p["d2"] ** 2
        return rbar * (d2 * d - 2 * d1**2) / d**3

    def _centrifugal_coefficient(self, l: ArrayLike) -> ArrayLike:
        """ħ²l(l+1)/(2μ), MeV fm²."""
        l = np.asarray(l, dtype=float)
        return HBARC_MEV_FM**2 * l * (l + 1) / (2 * self.mu)

    def potential_mev(self, r: ArrayLike, l: int = 0) -> ArrayLike:
        """V_l(r) = Z₁Z₂e²/r + V_Bass(r) + ħ²l(l+1)/(2μr²), MeV, at separation ``r`` (fm)."""
        rr = np.asarray(r, dtype=float)
        out = self.k / rr + self._nuclear(rr) + self._centrifugal_coefficient(l) / rr**2
        return _out(out, r)

    def _dv(self, r: np.ndarray, c: np.ndarray) -> np.ndarray:
        return -self.k / r**2 + self._nuclear(r, 1) - 2 * c / r**3

    def _d2v(self, r: np.ndarray, c: np.ndarray) -> np.ndarray:
        return 2 * self.k / r**3 + self._nuclear(r, 2) + 6 * c / r**4

    # -- barriers -------------------------------------------------------------------------------------------------

    def barriers(self) -> list:
        """The barrier of every partial wave that has one, l = 0, 1, ... up to the critical angular momentum (the
        last l whose potential still has a pocket). A list of :class:`Barrier`, indexed by l."""
        if self._barriers is None:
            self._barriers = self._find_barriers()
        return self._barriers

    def _find_barriers(self) -> list:
        rt = self.r1 + self.r2
        r = np.linspace(max(0.3 * rt, 1.0), rt + 15.0, 3000)
        ls = np.arange(L_LIMIT + 1)
        c = self._centrifugal_coefficient(ls)[:, None]
        dv = self._dv(r[None, :], c)
        out = []
        for i, l in enumerate(ls):
            # The barrier's top: the outermost place where dV/dr goes from + (inside) to − (outside).
            sign = np.nonzero((dv[i, :-1] > 0) & (dv[i, 1:] <= 0))[0]
            if sign.size == 0:
                break  # no pocket: this and all higher partial waves do not fuse
            j = sign[-1]
            lo, hi = r[j], r[j + 1]
            cl = c[i, 0]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if self._dv(np.array(mid), cl) > 0:
                    lo = mid
                else:
                    hi = mid
            rb = 0.5 * (lo + hi)
            curv = -float(self._d2v(np.array(rb), cl))
            if curv <= 0:
                break
            vb = float(self.potential_mev(rb, int(l)))
            out.append(Barrier(int(l), vb, float(rb), HBARC_MEV_FM * math.sqrt(curv / self.mu)))
        if not out:
            raise ValueError(f"{self.beam} + {self.target}: the Bass potential has no barrier")
        return out

    def barrier(self, l: int = 0) -> Barrier:
        """The barrier of partial wave ``l`` (default the s-wave: the fusion barrier)."""
        b = self.barriers()
        if not 0 <= l < len(b):
            raise ValueError(f"partial wave l = {l} has no barrier (the critical angular momentum is "
                             f"{self.critical_l})")
        return b[l]

    @property
    def critical_l(self) -> int:
        """The highest partial wave whose potential has a pocket, ħ. Higher ones do not fuse."""
        return len(self.barriers()) - 1

    @property
    def barrier_lab_mev(self) -> float:
        """The s-wave barrier as a lab beam energy, MeV (non-relativistic, as is customary)."""
        return self.barrier().height_mev * (self.m1 + self.m2) / self.m2

    # -- kinematics -----------------------------------------------------------------------------------------------

    def _sqrt_s(self, e_lab: ArrayLike) -> ArrayLike:
        return np.sqrt((self.m1 + self.m2) ** 2 + 2 * self.m2 * np.asarray(e_lab, dtype=float))

    def e_cm(self, beam_energy: EnergyLike) -> ArrayLike:
        """CM kinetic energy for a lab beam energy, MeV: √s − m₁ − m₂."""
        e = _energy_mev(beam_energy)
        return _out(self._sqrt_s(e) - self.m1 - self.m2, e)

    def beam_energy_for_cm(self, e_cm: ArrayLike) -> ArrayLike:
        """The lab beam energy that gives CM kinetic energy ``e_cm`` (MeV), MeV."""
        s = (np.asarray(e_cm, dtype=float) + self.m1 + self.m2) ** 2
        return _out((s - (self.m1 + self.m2) ** 2) / (2 * self.m2), e_cm)

    @property
    def q_value_mev(self) -> float:
        """Q-value of fusion m₁ + m₂ − M_CN, MeV (nuclear masses)."""
        return self.m1 + self.m2 - self.compound_nuclide.nuclear_mass_mev

    def excitation_energy(self, beam_energy: EnergyLike) -> ArrayLike:
        """Excitation energy of the compound nucleus, MeV: E* = √s − M_CN."""
        e = _energy_mev(beam_energy)
        return _out(self._sqrt_s(e) - self.compound_nuclide.nuclear_mass_mev, e)

    # -- cross sections -------------------------------------------------------------------------------------------

    def _lambda2_fm2(self, e_cm: np.ndarray) -> np.ndarray:
        """ƛ² = ħ²/(2μE), fm²."""
        with np.errstate(divide="ignore"):
            return HBARC_MEV_FM**2 / (2 * self.mu * e_cm)

    def partial_waves_cm(self, e_cm: float, l_dependent: bool = True, l_max: Optional[int] = None) -> np.ndarray:
        """σ_l = πƛ²(2l + 1)T_l at CM energy ``e_cm`` (MeV), mb, indexed by l.

        With ``l_dependent`` (the default) each l has its own barrier V_l, R_l, ħω_l, and partial waves above the
        critical angular momentum do not fuse. Without it, every l sees the s-wave parabola raised by
        ħ²l(l+1)/(2μR_B²), as in Wong's derivation; the sum then runs to ``l_max`` (default: until T_l < 10⁻¹⁶).
        """
        e = float(e_cm)
        if e <= 0:
            return np.zeros(1)
        lam2 = float(self._lambda2_fm2(np.array(e)))
        if l_dependent:
            b = self.barriers()
            v = np.array([x.height_mev for x in b])
            hw = np.array([x.hbar_omega_mev for x in b])
            ls = np.arange(len(b))
        else:
            b0 = self.barrier()
            if l_max is None:
                # T_l < 1e-16 once V_l − E > 6 ħω.
                top = max(e - b0.height_mev + 6 * b0.hbar_omega_mev, 0.0)
                l_max = int(math.sqrt(2 * self.mu * top) * b0.radius_fm / HBARC_MEV_FM) + 2
            ls = np.arange(l_max + 1)
            v = b0.height_mev + self._centrifugal_coefficient(ls) / b0.radius_fm**2
            hw = np.full(ls.size, b0.hbar_omega_mev)
        sigma = math.pi * lam2 * (2 * ls + 1) * _transmission(v, hw, e) * 10.0  # fm² → mb
        return sigma

    def cross_section_cm(self, e_cm: ArrayLike, l_dependent: bool = True) -> ArrayLike:
        """σ_fus at CM energy ``e_cm`` (MeV), mb: the sum of :meth:`partial_waves_cm`."""
        es = np.atleast_1d(np.asarray(e_cm, dtype=float))
        out = np.array([self.partial_waves_cm(e, l_dependent).sum() for e in es])
        return _out(out if np.ndim(e_cm) else out[0], e_cm)

    def wong_cm(self, e_cm: ArrayLike) -> ArrayLike:
        """Wong's closed form, mb: σ = (ħω R_B²/2E) ln{1 + exp[2π(E − V_B)/ħω]} with the s-wave barrier."""
        b = self.barrier()
        e = np.asarray(e_cm, dtype=float)
        x = 2 * np.pi * (e - b.height_mev) / b.hbar_omega_mev
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where(e > 0, b.hbar_omega_mev * b.radius_fm**2 / (2 * e) * np.logaddexp(0.0, x), 0.0) * 10.0
        return _out(out, e_cm)

    def classical_cm(self, e_cm: ArrayLike) -> ArrayLike:
        """The sharp cut-off limit πR_B²(1 − V_B/E) (0 below the barrier), mb."""
        b = self.barrier()
        e = np.asarray(e_cm, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = np.where(e > b.height_mev, math.pi * b.radius_fm**2 * (1 - b.height_mev / e), 0.0) * 10.0
        return _out(out, e_cm)

    def cross_section(self, beam_energy: EnergyLike) -> ArrayLike:
        """σ_fus at a lab beam energy, mb (the excitation function at those energies)."""
        return self.cross_section_cm(self.e_cm(beam_energy))

    def excitation_function(self, beam_energies: EnergyLike, thickness=None, material=None,
                            density=None) -> np.ndarray:
        """σ_fus (mb) at each lab beam energy (MeV), or with ``thickness`` its average over the target (see
        :meth:`thick_target`), the beam energy being that entering the target."""
        es = np.atleast_1d(_energy_mev(beam_energies))
        if thickness is None:
            return np.atleast_1d(self.cross_section(es))
        return np.array([self.thick_target(e, thickness, material, density).sigma_mb for e in es])

    def compound(self, beam_energy: EnergyLike) -> CompoundNucleus:
        """The compound nucleus made at one lab beam energy, with σ_l (see the module docstring)."""
        e = _energy_mev(beam_energy)
        if np.ndim(e):
            raise ValueError("compound() takes one beam energy")
        sqrt_s = float(self._sqrt_s(e))
        p = math.sqrt(e * (e + 2 * self.m1))
        e_cm = sqrt_s - self.m1 - self.m2
        sigma_l = self.partial_waves_cm(e_cm)
        nz = np.nonzero(sigma_l > 0)[0]
        sigma_l = sigma_l[: nz[-1] + 1] if nz.size else np.zeros(1)
        return CompoundNucleus(nuclide=self.compound_nuclide, beam_energy_mev=float(e), e_cm_mev=e_cm,
                               excitation_mev=sqrt_s - self.compound_nuclide.nuclear_mass_mev,
                               beta=p / (e + self.m1 + self.m2), sigma_l_mb=sigma_l)

    # -- thick targets --------------------------------------------------------------------------------------------

    def _stopping(self, material, density):
        from .stopping import Stopping  # imported here to keep this module light

        return Stopping(self.beam.name, self.target.name if material is None else material, density)

    def thick_target(self, beam_energy: EnergyLike, thickness, material=None, density=None,
                     points: int = 201) -> ThickTarget:
        """σ_fus averaged over the target's depth: (1/t) ∫₀ᵗ σ(E(x)) dx, with E(x) the beam energy after depth x.

        ``thickness`` is mg/cm² or a ``Quantity``/text (a length needs the material's density). ``material`` is
        what :func:`physim.nuclear.data.material` accepts; by default the target nuclide itself (an enriched
        target). The energy loss comes from :mod:`physim.nuclear.stopping`; the integral is Simpson's rule on
        ``points`` depths (odd). Beyond the beam's range σ is zero but the depth still counts.
        """
        e_in = _energy_mev(beam_energy)
        if np.ndim(e_in):
            raise ValueError("thick_target() takes one beam energy")
        stop = self._stopping(material, density)
        if isinstance(thickness, (int, float, np.floating)):
            t = float(thickness)
        else:
            t = stop.material.areal_density_mg_cm2(thickness)
        if t <= 0:
            raise ValueError("the target thickness must be positive")
        n = points if points % 2 else points + 1
        x = np.linspace(0.0, t, n)
        e = np.array([e_in] + [stop.energy_after(e_in, float(xi)) for xi in x[1:]], dtype=float)
        e_cm = np.where(e > 0, self._sqrt_s(e) - self.m1 - self.m2, 0.0)
        waves = [self.partial_waves_cm(ec) for ec in e_cm]
        width = max(w.size for w in waves)
        table = np.zeros((n, width))
        for i, w in enumerate(waves):
            table[i, : w.size] = w
        weights = np.ones(n)
        weights[1:-1:2], weights[2:-1:2] = 4.0, 2.0
        weights /= weights.sum()  # Simpson's rule divided by the thickness: the depth average
        sigma_l = weights @ table
        nz = np.nonzero(sigma_l > 0)[0]
        sigma_l = sigma_l[: nz[-1] + 1] if nz.size else np.zeros(1)
        sigma_x = table.sum(axis=1)
        return ThickTarget(e_in_mev=float(e_in), e_out_mev=float(e[-1]), thickness_mg_cm2=float(t),
                           sigma_mb=float(weights @ sigma_x), sigma_l_mb=sigma_l, depth_mg_cm2=x,
                           beam_energy_mev=e, e_cm_mev=e_cm,
                           excitation_mev=np.where(e > 0, self._sqrt_s(e) - self.compound_nuclide.nuclear_mass_mev,
                                                   np.nan),
                           sigma_depth_mb=sigma_x)

    # -- validity -------------------------------------------------------------------------------------------------

    def warnings(self, beam_energy: EnergyLike, thickness=None, material=None, density=None) -> list[str]:
        """Plain-language warnings about where this model stops being reliable, at one lab beam energy."""
        e = _energy_mev(beam_energy)
        if np.ndim(e):
            raise ValueError("warnings() takes one beam energy")
        out = []
        b = self.barrier()
        e_cm = float(self.e_cm(e))
        z1z2 = self.beam.Z * self.target.Z
        if z1z2 > QUASI_FISSION_Z1Z2:
            out.append(f"Z₁Z₂ = {z1z2} is above about {QUASI_FISSION_Z1Z2}: quasi-fission competes with fusion, "
                       f"and the capture cross section given here overestimates fusion, often by orders of "
                       f"magnitude.")
        if e_cm < b.height_mev:
            out.append(f"The beam ({e_cm:.1f} MeV CM) is below the fusion barrier ({b.height_mev:.1f} MeV CM): "
                       f"fusion falls by orders of magnitude per few MeV, coupling to the nuclei's excited states "
                       f"raises it well above this estimate (coupled channels), and the parabolic barrier is "
                       f"rough so far down.")
        elif e_cm < b.height_mev + b.hbar_omega_mev:
            out.append(f"The beam ({e_cm:.1f} MeV CM) is near the fusion barrier ({b.height_mev:.1f} MeV CM): "
                       f"coupled-channels effects can change the cross section noticeably here.")
        per_u = e / self.beam.A
        if per_u > INCOMPLETE_FUSION_MEV_PER_U:
            out.append(f"At {per_u:.1f} MeV/u, above about {INCOMPLETE_FUSION_MEV_PER_U:g} MeV/u, incomplete "
                       f"fusion and pre-equilibrium emission take a sizeable share: complete fusion is "
                       f"overestimated.")
        lc = self.critical_l
        cn = self.compound(e)
        if cn.sigma_mb > 0 and cn.max_l >= lc:
            out.append(f"Partial waves up to the critical angular momentum ({lc} ħ) fuse: σ_fus is limited by "
                       f"the pocket of the potential, not by the barrier, and the highest spins may fission "
                       f"instead.")
        if thickness is not None:
            tt = self.thick_target(e, thickness, material, density)
            if tt.e_out_mev <= 0:
                out.append("The beam stops in the target: the depth-averaged cross section counts the depth "
                           "beyond its range as giving no fusion.")
            elif self.e_cm(tt.e_out_mev) < b.height_mev < e_cm:
                out.append(f"The beam drops below the fusion barrier inside the target (it leaves with "
                           f"{tt.e_out_mev:.1f} MeV): the back of the target adds little.")
        return out
