"""Relativistic two-body reaction kinematics: target(beam, ejectile)recoil, for example 208Pb(16O, 16O)208Pb.

The beam hits a target at rest. In the centre-of-mass (CM) frame the ejectile leaves at angle θ* with momentum
p*; boosting back to the lab gives its energy and angle there. Everything is exact special relativity, with
nuclear (bare) masses from AME2020 and any excitation energy E* added to the recoil (or, with ``excite="ejectile"``,
to the ejectile)::

    from physim.nuclear.kinematics import TwoBody

    r = TwoBody("16O", "208Pb", beam_energy="64 MeV")       # elastic scattering
    r.at_cm(90.0)                                           # ejectile lab angle, energy, Jacobian, dE/dθ
    fwd, bwd = r.at_lab(30.0)                               # both solutions where double-valued
    TwoBody("208Pb", "1H", beam_energy="1 GeV").max_angle() # inverse kinematics: 0.278 deg

Angles are in degrees and energies in MeV throughout. Arrays work wherever a single angle does.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Union

import numpy as np

from . import data
from .quantity import Quantity

ArrayLike = Union[float, np.ndarray]

# g = 1 holds exactly for elastic recoils and equal-mass elastic scattering; allow for round-off.
_G_TOL = 1e-12


@dataclass(frozen=True)
class LabPoint:
    """One kinematic solution, for the ejectile or the recoil. Each field is a float or an array."""

    #: CM angle of this particle, degrees.
    theta_cm: ArrayLike
    #: Lab angle, degrees.
    theta_lab: ArrayLike
    #: Lab kinetic energy, MeV.
    energy: ArrayLike
    #: dΩ_cm/dΩ_lab: multiply a CM cross section by this to get the lab cross section.
    jacobian: ArrayLike
    #: dE/dθ_lab, MeV per degree: how fast the energy changes across a detector (kinematic broadening).
    de_dtheta: ArrayLike


def _energy(value: Union[float, str, Quantity]) -> float:
    if isinstance(value, str):
        value = Quantity.parse(value)
    if isinstance(value, Quantity):
        if value.kind != "energy":
            raise ValueError(f"beam energy '{value}' must be a total energy such as '64 MeV'")
        return value.to("MeV")
    return float(value)


class TwoBody:
    """Kinematics of ``target(beam, ejectile)recoil`` at one beam energy.

    ``ejectile`` defaults to the beam (elastic or inelastic scattering); ``recoil`` is whatever conserves charge
    and mass number. ``excitation_mev`` is left in the recoil, or in the ejectile with ``excite="ejectile"``.
    """

    def __init__(self, beam: str, target: str, beam_energy: Union[float, str, Quantity],
                 ejectile: Optional[str] = None, recoil: Optional[str] = None, excitation_mev: float = 0.0,
                 excite: str = "recoil"):
        n1, n2 = data.nuclide(beam), data.nuclide(target)
        n3 = data.nuclide(ejectile) if ejectile is not None else n1
        if recoil is None:
            z4, a4 = n1.Z + n2.Z - n3.Z, n1.A + n2.A - n3.A
            if z4 < 0 or a4 < 1 or a4 < z4:
                raise ValueError(f"no recoil conserves charge and mass number for {n2}({n1}, {n3})")
            n4 = data.nuclide((z4, a4))
        else:
            n4 = data.nuclide(recoil)
            if (n1.Z + n2.Z, n1.A + n2.A) != (n3.Z + n4.Z, n3.A + n4.A):
                raise ValueError(f"{n2}({n1}, {n3}){n4} does not conserve charge and mass number")
        if excite not in ("recoil", "ejectile"):
            raise ValueError("excite must be 'recoil' or 'ejectile'")
        if excitation_mev < 0:
            raise ValueError("excitation_mev must not be negative")
        self.beam, self.target, self.ejectile, self.recoil = n1, n2, n3, n4
        self.excitation_mev = float(excitation_mev)
        self.excite = excite
        self.beam_energy_mev = _energy(beam_energy)
        if self.beam_energy_mev <= 0:
            raise ValueError("the beam energy must be positive")

        m1, m2 = n1.nuclear_mass_mev, n2.nuclear_mass_mev
        m3, m4 = n3.nuclear_mass_mev, n4.nuclear_mass_mev
        if excite == "recoil":
            m4 += self.excitation_mev
        else:
            m3 += self.excitation_mev
        self.m1, self.m2, self.m3, self.m4 = m1, m2, m3, m4
        t = self.beam_energy_mev
        p1 = math.sqrt(t * (t + 2 * m1))
        e_total = t + m1 + m2
        self.s = (m1 + m2) ** 2 + 2 * m2 * t
        self.sqrt_s = math.sqrt(self.s)
        #: Velocity (units of c) and Lorentz factor of the CM frame in the lab.
        self.beta_cm = p1 / e_total
        self.gamma_cm = e_total / self.sqrt_s
        if self.sqrt_s < m3 + m4:
            raise ValueError(f"{self.label}: {self.beam_energy_mev:g} MeV is below the threshold of "
                             f"{self.threshold_mev:.4f} MeV")
        # s − (m3 ± m4)², written so nothing large cancels: at low energy s and (m3 + m4)² agree to many digits.
        q = (m1 + m2) - (m3 + m4)
        plus = q * (m1 + m2 + m3 + m4) + 2 * m2 * t
        minus = (m1 + m2 - m3 + m4) * (m1 + m2 + m3 - m4) + 2 * m2 * t
        #: CM momentum of the outgoing pair, MeV/c.
        self.p_cm = math.sqrt(plus * minus) / (2 * self.sqrt_s)
        self._e_cm = {3: math.hypot(self.p_cm, m3), 4: math.hypot(self.p_cm, m4)}

    # -- reaction properties ----------------------------------------------------------------------------------

    @property
    def label(self) -> str:
        star = "*" if self.excitation_mev else ""
        c = f"{self.ejectile}{star if self.excite == 'ejectile' else ''}"
        d = f"{self.recoil}{star if self.excite == 'recoil' else ''}"
        return f"{self.target}({self.beam}, {c}){d}"

    @property
    def q_value_mev(self) -> float:
        """Q-value including the excitation energy, MeV (nuclear masses)."""
        return self.m1 + self.m2 - self.m3 - self.m4

    @property
    def threshold_mev(self) -> float:
        """Lowest beam kinetic energy at which the reaction can happen, MeV (0 when Q ≥ 0)."""
        need = (self.m3 + self.m4) ** 2 - (self.m1 + self.m2) ** 2
        return max(0.0, need / (2 * self.m2))

    @property
    def cm_energy_mev(self) -> float:
        """Kinetic energy available in the CM frame before the reaction, MeV."""
        return self.sqrt_s - self.m1 - self.m2

    def _g(self, k: int) -> float:
        """β_cm divided by the CM velocity of particle k: above 1, its lab angle has a maximum."""
        return self.beta_cm * self._e_cm[k] / self.p_cm

    def _index(self, particle: str) -> int:
        if particle not in ("ejectile", "recoil"):
            raise ValueError("particle must be 'ejectile' or 'recoil'")
        return 3 if particle == "ejectile" else 4

    def double_valued(self, particle: str = "ejectile") -> bool:
        """True if two CM angles give the same lab angle (so a detector sees two energies)."""
        return self._g(self._index(particle)) > 1.0 + _G_TOL

    def max_angle(self, particle: str = "ejectile") -> float:
        """Largest lab angle the particle can reach, degrees (180 if it can go anywhere)."""
        k = self._index(particle)
        g = self._g(k)
        if abs(g - 1.0) <= _G_TOL:  # e.g. elastic scattering of equal masses, or any elastic recoil
            return 90.0
        if g < 1.0:
            return 180.0
        return math.degrees(math.atan(1.0 / (self.gamma_cm * math.sqrt(g * g - 1.0))))

    # -- CM angle -> lab ---------------------------------------------------------------------------------------

    def at_cm(self, theta_cm: ArrayLike, particle: str = "ejectile") -> LabPoint:
        """Lab quantities of the ejectile (or recoil) leaving at CM angle ``theta_cm`` (degrees)."""
        k = self._index(particle)
        m = self.m3 if k == 3 else self.m4
        th = np.radians(np.asarray(theta_cm, dtype=float))
        c, s = np.cos(th), np.sin(th)
        g, gam, p, e = self._g(k), self.gamma_cm, self.p_cm, self._e_cm[k]
        p_par = gam * p * (c + g)  # along the beam
        p_perp = p * s
        theta_lab = np.degrees(np.arctan2(p_perp, p_par))
        energy = gam * (e + self.beta_cm * p * c) - m
        # dθ_lab/dθ_cm from tan θ_lab = sin θ* / (γ (cos θ* + g)).
        denom = gam**2 * (c + g) ** 2 + s**2
        dlab_dcm = gam * (1.0 + g * c) / denom
        with np.errstate(divide="ignore", invalid="ignore"):
            sin_lab = np.sin(np.radians(theta_lab))
            # At θ* = 0 or 180° both sines vanish; the ratio tends to (dθ*/dθ_lab)².
            jac = np.where(np.abs(s) > 1e-9, np.abs(s / (sin_lab * dlab_dcm)), 1.0 / dlab_dcm**2)
            de_dtheta = -gam * self.beta_cm * p * s / dlab_dcm * (math.pi / 180.0)  # per degree
        return LabPoint(*(_out(x, theta_cm) for x in (np.degrees(th), theta_lab, energy, jac, de_dtheta)))

    # -- lab angle -> CM ---------------------------------------------------------------------------------------

    def theta_cm(self, theta_lab: ArrayLike, particle: str = "ejectile") -> tuple:
        """CM angles (degrees) giving lab angle ``theta_lab``: (first, second).

        ``first`` is the solution with the smaller CM angle (the higher-energy one). ``second`` exists only where
        the kinematics are double-valued and is NaN otherwise; both are NaN beyond the maximum lab angle.
        """
        k = self._index(particle)
        g, gam = self._g(k), self.gamma_cm
        tl = np.radians(np.asarray(theta_lab, dtype=float))
        a, b = np.cos(tl), gam * np.sin(tl)
        # sin θ* cos θ_lab − γ sin θ_lab cos θ* = γ g sin θ_lab, i.e. R sin(θ* − δ) = b g.
        r = np.hypot(a, b)
        delta = np.arctan2(b, a)
        with np.errstate(invalid="ignore"):
            x = np.arcsin(np.clip(b * g / r, -1.0, 1.0))
            inside = np.abs(b * g / r) <= 1.0 + 1e-12
            cands = [delta + x, delta + np.pi - x]
            sols = []
            for cand in cands:
                cand = np.where(inside, cand, np.nan)
                ok = (cand >= -1e-12) & (cand <= np.pi + 1e-12)
                # The squared equation also admits θ_lab + π; keep solutions whose momentum points the right way.
                ok &= (np.sign(np.cos(cand) + g) == np.sign(a)) | (np.abs(a) < 1e-12)
                sols.append(np.where(ok, np.clip(cand, 0.0, np.pi), np.nan))
        first = np.fmin(sols[0], sols[1])
        second = np.where(np.isclose(sols[0], sols[1], rtol=0, atol=1e-12) | np.isnan(sols[0]) | np.isnan(sols[1]),
                          np.nan, np.fmax(sols[0], sols[1]))
        if g <= 1.0 + _G_TOL:
            second = np.full_like(first, np.nan)
        return _out(np.degrees(first), theta_lab), _out(np.degrees(second), theta_lab)

    def at_lab(self, theta_lab: ArrayLike, particle: str = "ejectile") -> tuple:
        """Lab quantities at lab angle ``theta_lab`` (degrees): (first, second) :class:`LabPoint` solutions.

        ``second`` is ``None`` when the particle's kinematics are single-valued.
        """
        first, second = self.theta_cm(theta_lab, particle)
        p1 = self.at_cm(first, particle)
        p2 = self.at_cm(second, particle) if self.double_valued(particle) else None
        return p1, p2

    def recoil_for(self, theta_cm_ejectile: ArrayLike) -> LabPoint:
        """The recoil partnering an ejectile at CM angle ``theta_cm_ejectile`` (degrees)."""
        return self.at_cm(180.0 - np.asarray(theta_cm_ejectile, dtype=float), "recoil")

    def __repr__(self) -> str:
        return f"TwoBody({self.label} at {self.beam_energy_mev:g} MeV)"


def _out(x: np.ndarray, like: ArrayLike) -> ArrayLike:
    return float(x) if np.ndim(like) == 0 else np.asarray(x, dtype=float)


def kinematic_factor(beam_mass: float, target_mass: float, theta_lab: ArrayLike) -> ArrayLike:
    """Non-relativistic elastic kinematic factor K = E_scattered / E_beam at lab angle θ (degrees).

    K = [(m₁ cos θ + √(m₂² − m₁² sin² θ)) / (m₁ + m₂)]² (the forward solution).
    """
    th = np.radians(np.asarray(theta_lab, dtype=float))
    m1, m2 = beam_mass, target_mass
    k = ((m1 * np.cos(th) + np.sqrt(m2**2 - (m1 * np.sin(th)) ** 2)) / (m1 + m2)) ** 2
    return _out(k, theta_lab)


def elastic(experiment) -> dict:
    """Elastic kinematics for every target nuclide of an :class:`~physim.nuclear.Experiment`.

    Returns ``{(Z, A): (TwoBody, fraction of target atoms)}``; a natural element gives one entry per isotope.
    """
    beam = experiment.beam
    mat = experiment.target.material_data()
    atoms = mat.atoms_per_cm2("1 mg/cm2")
    total = sum(atoms.values())
    out = {}
    for (z, a), n in atoms.items():
        out[(z, a)] = (TwoBody(data.nuclide((beam.Z, beam.A)).name, data.nuclide((z, a)).name, beam.energy_mev),
                       n / total)
    return out
