"""Rutherford scattering: cross sections, distance of closest approach, impact parameters, validity checks, and
classical Coulomb trajectories integrated with physim's engine.

::

    from physim.nuclear.rutherford import Rutherford

    r = Rutherford("4He", "197Au", beam_energy="20 MeV")
    r.cross_section_cm(50.0)          # dσ/dΩ in the CM frame, mb/sr
    r.cross_section_lab(49.115)       # in the lab frame
    r.closest_approach(50.0)          # fm
    r.grazing_angle()                 # CM angle where nuclear forces set in (180 if never)
    r.warnings()                      # plain-language validity warnings

Angles are CM angles of the ejectile in degrees unless a name says lab. Rutherford's formula is classical and
non-relativistic in the CM kinetic energy, which is what experiments and other codes (LISE++) use.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Union

import numpy as np

from . import data
from .kinematics import TwoBody
from .quantity import Quantity

#: e²/(4πε₀) in MeV fm (CODATA 2018).
E2_MEV_FM = 1.43996448
#: Fine-structure constant.
ALPHA = 7.2973525693e-3

ArrayLike = Union[float, np.ndarray]


def _out(x, like):
    return float(x) if np.ndim(like) == 0 else np.asarray(x, dtype=float)


class Rutherford:
    """Elastic Coulomb scattering of ``beam`` on ``target`` at a lab kinetic energy (MeV, or text with a unit).

    ``r0`` (fm) and ``margin`` (fm) set the interaction radius R = r0 (A₁^⅓ + A₂^⅓) + margin at which nuclear forces
    are assumed to start; the defaults (1.2 fm, 2 fm) put R just outside the touching nuclear surfaces.
    """

    def __init__(self, beam: str, target: str, beam_energy: Union[float, str, Quantity], r0: float = 1.2,
                 margin: float = 2.0):
        self.kinematics = TwoBody(beam, target, beam_energy)
        self.beam, self.target = self.kinematics.beam, self.kinematics.target
        self.r0, self.margin = float(r0), float(margin)
        #: Z₁Z₂e², MeV fm.
        self.k = self.beam.Z * self.target.Z * E2_MEV_FM
        #: Kinetic energy in the CM frame, MeV.
        self.e_cm = self.kinematics.cm_energy_mev
        #: Head-on distance of closest approach d₀ = Z₁Z₂e²/E_cm, fm.
        self.d0 = self.k / self.e_cm

    @property
    def beam_energy_mev(self) -> float:
        return self.kinematics.beam_energy_mev

    # -- cross sections -----------------------------------------------------------------------------------------

    def cross_section_cm(self, theta_cm: ArrayLike) -> ArrayLike:
        """dσ/dΩ in the CM frame, mb/sr: (Z₁Z₂e² / 4E_cm)² / sin⁴(θ/2)."""
        s = np.sin(np.radians(np.asarray(theta_cm, dtype=float)) / 2)
        with np.errstate(divide="ignore"):
            out = (self.k / (4 * self.e_cm)) ** 2 / s**4 * 10.0  # fm² → mb
        return _out(out, theta_cm)

    def cross_section_lab(self, theta_lab: ArrayLike, particle: str = "ejectile") -> tuple:
        """dσ/dΩ in the lab at lab angle ``theta_lab`` for the ejectile (or the recoil), mb/sr.

        Returns (first, second): ``second`` is the other kinematic solution where the kinematics are
        double-valued, else ``None``. For the recoil, the cross section is that of the ejectile partner at
        180° − θ*.
        """
        out = []
        for point in self.kinematics.at_lab(theta_lab, particle):
            if point is None:
                out.append(None)
                continue
            th_ej = point.theta_cm if particle == "ejectile" else 180.0 - np.asarray(point.theta_cm)
            out.append(_out(self.cross_section_cm(th_ej) * np.asarray(point.jacobian), theta_lab))
        return tuple(out)

    def integrated(self, theta_cm_min: float, theta_cm_max: float, dphi_deg: float = 360.0) -> float:
        """∫ dσ/dΩ dΩ between two CM angles over an azimuthal range, mb (analytic)."""
        a, b = (math.radians(x) / 2 for x in (theta_cm_min, theta_cm_max))
        if a <= 0:
            raise ValueError("the Rutherford cross section integrated down to 0° is infinite")
        # With u = sin(θ/2): sinθ dθ = 4u du, so ∫ sinθ dθ / sin⁴(θ/2) = [−2/u²] between the limits.
        val = (self.k / (4 * self.e_cm)) ** 2 * 2 * (1 / math.sin(a) ** 2 - 1 / math.sin(b) ** 2) * 10.0
        return val * math.radians(dphi_deg)

    # -- geometry of the orbit ----------------------------------------------------------------------------------

    def closest_approach(self, theta_cm: ArrayLike) -> ArrayLike:
        """Distance of closest approach for scattering to ``theta_cm``, fm: (d₀/2)(1 + 1/sin(θ/2))."""
        s = np.sin(np.radians(np.asarray(theta_cm, dtype=float)) / 2)
        with np.errstate(divide="ignore"):
            return _out(self.d0 / 2 * (1 + 1 / s), theta_cm)

    def closest_approach_for_impact(self, b: ArrayLike) -> ArrayLike:
        """Distance of closest approach for impact parameter ``b`` (fm): d₀/2 + √((d₀/2)² + b²)."""
        b = np.asarray(b, dtype=float)
        return _out(self.d0 / 2 + np.sqrt((self.d0 / 2) ** 2 + b**2), b)

    def impact_parameter(self, theta_cm: ArrayLike) -> ArrayLike:
        """Impact parameter that scatters to ``theta_cm``, fm: b = (d₀/2) cot(θ/2)."""
        t = np.radians(np.asarray(theta_cm, dtype=float)) / 2
        with np.errstate(divide="ignore"):
            return _out(self.d0 / 2 / np.tan(t), theta_cm)

    def angle_for_impact(self, b: ArrayLike) -> ArrayLike:
        """CM scattering angle for impact parameter ``b`` (fm), degrees: θ = 2 arctan(d₀ / 2b)."""
        b = np.asarray(b, dtype=float)
        return _out(np.degrees(2 * np.arctan2(self.d0 / 2, b)), b)

    # -- where Rutherford stops being valid ---------------------------------------------------------------------

    @property
    def interaction_radius(self) -> float:
        """R = r0 (A₁^⅓ + A₂^⅓) + margin, fm: closer than this, nuclear forces act."""
        return self.r0 * (self.beam.A ** (1 / 3) + self.target.A ** (1 / 3)) + self.margin

    @property
    def coulomb_barrier_cm(self) -> float:
        """Coulomb barrier at the interaction radius, CM frame, MeV: Z₁Z₂e²/R."""
        return self.k / self.interaction_radius

    @property
    def coulomb_barrier_lab(self) -> float:
        """The same barrier as a lab beam energy, MeV (non-relativistic conversion)."""
        m1, m2 = self.beam.atomic_mass_u, self.target.atomic_mass_u
        return self.coulomb_barrier_cm * (m1 + m2) / m2

    def grazing_angle(self) -> float:
        """CM angle at which the closest approach reaches the interaction radius, degrees. Rutherford holds at
        smaller angles; 180 means it holds everywhere (the beam is below the barrier)."""
        x = 2 * self.interaction_radius / self.d0 - 1
        if x <= 1:
            return 180.0
        return math.degrees(2 * math.asin(1 / x))

    @property
    def sommerfeld(self) -> float:
        """Sommerfeld parameter η = Z₁Z₂e²/(ħv), v the relative velocity (the beam's, with the target at rest):
        ≫ 1 means classical orbits."""
        t, m = self.kinematics.beam_energy_mev, self.kinematics.m1
        beta = math.sqrt(t * (t + 2 * m)) / (t + m)
        return self.beam.Z * self.target.Z * ALPHA / beta

    def screening_correction(self) -> float:
        """Electron screening factor on the cross section at large angles: 1 − 0.049 Z₁ Z₂^(4/3) / E_cm[keV]
        (L'Ecuyer, Davies and Matsunami 1979)."""
        return 1.0 - 0.049 * self.beam.Z * self.target.Z ** (4 / 3) / (self.e_cm * 1e3)

    def identical(self) -> bool:
        """True for identical particles, where Mott scattering (interference) applies instead."""
        return (self.beam.Z, self.beam.A) == (self.target.Z, self.target.A)

    def warnings(self, theta_cm_max: Optional[float] = None) -> list[str]:
        """Plain-language warnings about where this setup leaves Rutherford's formula."""
        out = []
        graze = self.grazing_angle()
        limit = 180.0 if theta_cm_max is None else theta_cm_max
        if graze < limit:
            out.append(f"Above {graze:.1f}° (CM) the nuclei come within {self.interaction_radius:.1f} fm of each "
                       f"other and nuclear forces act: the cross section there is not Rutherford. The beam "
                       f"({self.e_cm:.1f} MeV CM) is above the Coulomb barrier ({self.coulomb_barrier_cm:.1f} MeV).")
        if self.sommerfeld < 5:  # η of a few or more: semiclassical orbits are fine
            out.append(f"The Sommerfeld parameter is η = {self.sommerfeld:.1f}, not ≫ 1: classical orbits are only "
                       f"approximate.")
        scr = 1 - self.screening_correction()
        if scr > 0.01:
            out.append(f"Electron screening lowers the cross section by about {100 * scr:.1f}% at large angles "
                       f"(more at small angles).")
        if self.identical():
            out.append("Beam and target are identical nuclei: Mott scattering (interference) applies, not "
                       "Rutherford's formula.")
        return out

    # -- trajectories ---------------------------------------------------------------------------------------------

    def trajectories(self, impact_parameters: ArrayLike, distance: Optional[float] = None,
                     rtol: float = 1e-11) -> list:
        """Integrate classical Coulomb orbits in the CM frame with physim's engine.

        The relative motion is one particle of reduced mass μ in the field of a fixed charge Z₁Z₂e², started at
        ``distance`` (default 2000 d₀) on the incoming side of the orbit with impact parameter b, and followed to
        about the same distance on the way out. Units: fm, MeV, time in fm/c.
        Returns a list of :class:`Orbit`.
        """
        import physim as ps  # the engine; imported here to keep this module light

        m1 = self.beam.nuclear_mass_mev
        m2 = self.target.nuclear_mass_mev
        mu = m1 * m2 / (m1 + m2)
        v_inf = math.sqrt(2 * self.e_cm / mu)  # non-relativistic, in units of c
        start = 2000 * self.d0 if distance is None else float(distance)
        # At a finite distance the particle already has potential energy k/r: start slower, so that the total
        # energy is E_cm, and offset sideways so that the angular momentum is μ v∞ b. That is the orbit with
        # impact parameter b, picked up part-way along its incoming branch.
        v = math.sqrt(2 * (self.e_cm - self.k / start) / mu)
        orbits = []
        for b in np.atleast_1d(np.asarray(impact_parameters, dtype=float)):
            w = ps.World(integrator="dopri5")
            y0 = b * v_inf / v
            x0 = -math.sqrt(start**2 - y0**2)
            w.add_particle([x0, y0, 0.0], vel=[v, 0.0, 0.0], mass=mu, charge=1.0)
            w.add_particle([0.0, 0.0, 0.0], mass=1e30, charge=self.k / E2_MEV_FM / 1.0)
            w.pin(1)
            w.add_force(ps.Coulomb(k=E2_MEV_FM))
            t_end = 2 * start / v_inf
            traj = w.run_adaptive(t_end, rtol=rtol, atol=1e-12 * start, energies=False,
                                  times=np.linspace(0, t_end, 801))
            orbits.append(Orbit(b=float(b), mu=mu, k=self.k, pos=np.asarray(traj.pos)[:, 0, :2],
                                vel=np.asarray(traj.vel)[:, 0, :2]))
        return orbits


@dataclass
class Orbit:
    """A Coulomb orbit in the CM frame: positions (fm) and velocities (c) of the relative motion."""

    b: float
    mu: float
    k: float
    pos: np.ndarray
    vel: np.ndarray

    def _elements(self, i: int) -> tuple:
        """(e, ψ, φ, L) of the hyperbola through point i: e cos ψ = 1 + L²/(μkr), e sin ψ = L ṙ / k."""
        x, y = self.pos[i]
        vx, vy = self.vel[i]
        r = math.hypot(x, y)
        L = self.mu * (x * vy - y * vx)
        rdot = (x * vx + y * vy) / r
        ec, es = 1 + L**2 / (self.mu * self.k * r), L * rdot / self.k
        return math.hypot(ec, es), math.atan2(es, ec), math.atan2(y, x), L

    def deflection(self) -> float:
        """Scattering angle in degrees, from the asymptotes of the orbit through the first and last points.

        Each end point fixes the hyperbola through it; the incoming and outgoing asymptotes then follow exactly,
        so the result does not depend on how far out the integration starts and stops.
        """
        e0, psi0, phi0, L = self._elements(0)
        e1, psi1, phi1, _ = self._elements(-1)
        sign = math.copysign(1.0, L)  # φ grows with time for L > 0
        # ψ = φ − φ_p (φ_p: direction of closest approach); the orbit runs from ψ = −s·arccos(1/e) to +s·arccos(1/e).
        incoming = (phi0 - psi0) - sign * math.acos(1 / e0)  # direction the particle comes from
        outgoing = (phi1 - psi1) + sign * math.acos(1 / e1)  # direction it leaves towards
        # Velocity in: pointing from `incoming` towards the centre (incoming + π); velocity out: along `outgoing`.
        turn = (outgoing - incoming - math.pi + math.pi) % (2 * math.pi) - math.pi
        return math.degrees(abs(turn))

    def closest_approach(self) -> float:
        """Smallest distance reached along the recorded orbit, fm."""
        return float(np.min(np.hypot(self.pos[:, 0], self.pos[:, 1])))
