"""Coulomb excitation in first-order semiclassical theory (Alder and Winther).

Below the Coulomb barrier the nuclei never touch: the electric field of one nucleus excites the other as they pass.
The relative motion is a classical Rutherford orbit, and the probability of exciting a state of energy E* with
multipolarity λ (from a 0⁺ ground state) is first-order perturbation theory along that orbit::

    from physim.nuclear.coulex import Coulex

    c = Coulex("16O", "58Ni", "48 MeV", excite="target", energy="1.454 MeV", multipolarity="E2",
               b_up="0.0695 e2b2")
    c.probability(120.0)              # excitation probability at a CM angle (degrees)
    c.cross_section_cm(120.0)         # dσ/dΩ for exciting the state, mb/sr
    c.total()                         # integrated over all angles, mb
    c.xi                              # adiabaticity ξ: excitation is suppressed when ξ ≳ 1
    c.safe(120.0)                     # Cline's safe-distance criterion at that angle

Units: angles in degrees (CM unless a name says lab), energies MeV, lengths fm, B(Eλ) in e² fm^(2λ).
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Optional, Union

import numpy as np

from . import data
from .kinematics import TwoBody
from .quantity import Quantity
from .rutherford import E2_MEV_FM

#: ħc, MeV fm (CODATA 2018).
HBARC_MEV_FM = 197.3269804
#: Cline's safe distance: closest approach at least 1.25 (A₁^⅓ + A₂^⅓) + this many fm.
SAFE_MARGIN_FM = 5.0

#: Units of B(Eλ↑) and their value in e² fm^(2λ), by λ.
B_UNITS = {1: {"e2fm2": 1.0, "e2b": 100.0}, 2: {"e2fm4": 1.0, "e2b2": 1e4}, 3: {"e2fm6": 1.0, "e2b3": 1e6}}

ArrayLike = Union[float, np.ndarray]


def parse_b(text: Union[str, float], lam: int) -> float:
    """B(Eλ↑) from text such as ``"0.0695 e2b2"`` or ``"695 e2fm4"``, in e² fm^(2λ)."""
    if isinstance(text, (int, float)):
        return float(text)
    parts = str(text).split()
    if len(parts) != 2:
        raise ValueError(f"write B(E{lam}) with its unit, e.g. '0.05 {list(B_UNITS[lam])[-1]}'")
    value, unit = float(parts[0]), parts[1]
    if unit not in B_UNITS[lam]:
        raise ValueError(f"unit '{unit}' is not a unit of B(E{lam}); use {' or '.join(B_UNITS[lam])}")
    if value < 0:
        raise ValueError("B(Eλ) must not be negative")
    return value * B_UNITS[lam][unit]


@lru_cache(maxsize=None)
def ylm_equator(lam: int, mu: int) -> float:
    """Y_λμ(θ = π/2, φ = 0), the spherical harmonic on the orbit plane (Condon–Shortley phase)."""
    m = abs(mu)
    leg = np.polynomial.legendre.Legendre.basis(lam).deriv(m)
    plm = (-1) ** m * float(leg(0.0))
    norm = math.sqrt((2 * lam + 1) / (4 * math.pi) * math.factorial(lam - m) / math.factorial(lam + m))
    y = norm * plm
    return y if mu >= 0 else (-1) ** m * y


def orbit_integrals(lam: int, epsilon: float, xi: float, tol: float = 1e-9) -> np.ndarray:
    """I_μ = ∫ e^{iξ(ε sinh w + w)} e^{iμφ(w)} / (ε cosh w + 1)^λ dw over the hyperbolic orbit, for μ = −λ…λ.

    The orbit is r = a(ε cosh w + 1), t = (a/v)(ε sinh w + w), with the position x = a(cosh w + ε),
    y = a √(ε² − 1) sinh w measured from the scattering centre. The integral is done in Gauss–Legendre panels,
    each spanning about a radian of the phase, out to where the integrand has fallen below ``tol``.
    """
    eps_ = max(epsilon, 1.0)
    # Out to where (ε cosh w)^−λ < tol, or (for ξ > 0) where the fast oscillation has averaged it away.
    w_amp = math.acosh(max((1 / tol) ** (1 / lam) / eps_, 1.0)) + 1.0
    w_max = w_amp
    if xi > 0:
        # Beyond w the integrand is ~ e^{−λw}, and oscillates with dφ/dw ≈ ξ ε e^w/2: what is left is below
        # (ε e^w/2)^−λ / (ξ ε e^w/2).
        w_osc = math.log(2 * (1 / (tol * xi)) ** (1 / (lam + 1)) / eps_) if tol * xi < 1 else 0.0
        w_max = min(w_amp, max(w_osc, 3.0) + 1.0)
    edges = [0.0]
    while edges[-1] < w_max:
        w = edges[-1]
        rate = xi * (eps_ * math.cosh(w) + 1.0)
        edges.append(min(w + min(0.5, 2.0 / (rate + 1e-300)), w_max))
    edges = np.array(edges)
    x, wts = np.polynomial.legendre.leggauss(16)
    lo, hi = edges[:-1], edges[1:]
    w = ((hi - lo)[:, None] * (x[None, :] + 1) / 2 + lo[:, None]).ravel()
    weight = ((hi - lo)[:, None] / 2 * wts[None, :]).ravel()
    w = np.concatenate([-w[::-1], w])
    weight = np.concatenate([weight[::-1], weight])
    ch, sh = np.cosh(w), np.sinh(w)
    r = eps_ * ch + 1.0
    phase = np.exp(1j * xi * (eps_ * sh + w))
    e_iphi = ((ch + eps_) + 1j * math.sqrt(max(eps_ * eps_ - 1.0, 0.0)) * sh) / r
    base = weight * phase / r**lam
    return np.array([np.sum(base * e_iphi**mu) for mu in range(-lam, lam + 1)])


class Coulex:
    """First-order Coulomb excitation of one state, of the target or of the projectile.

    ``energy`` is E* (MeV or text), ``multipolarity`` "E1", "E2" or "E3" (from a 0⁺ ground state), ``b_up``
    B(Eλ↑) in e² fm^(2λ) or text with a unit (``"0.0695 e2b2"``). The orbit is symmetrised between the initial and
    final velocities (Alder and Winther): a = Z₁Z₂e²/(μ v_i v_f), ξ = (Z₁Z₂e²/ħ)(1/v_f − 1/v_i).
    """

    def __init__(self, beam: str, target: str, beam_energy, *, excite: str = "target", energy=None,
                 multipolarity: str = "E2", b_up=None, r0: float = 1.25):
        if excite not in ("target", "projectile"):
            raise ValueError("excite must be 'target' or 'projectile'")
        if multipolarity not in ("E1", "E2", "E3"):
            raise ValueError("multipolarity must be E1, E2 or E3")
        self.lam = int(multipolarity[1])
        self.excite = excite
        q = Quantity.parse(energy) if isinstance(energy, str) else energy
        self.e_star = q.to("MeV") if isinstance(q, Quantity) else float(q)
        if self.e_star <= 0:
            raise ValueError("the excitation energy must be positive")
        self.b_up = parse_b(b_up, self.lam)
        self.kinematics = TwoBody(beam, target, beam_energy, excitation_mev=self.e_star,
                                  excite="recoil" if excite == "target" else "ejectile")
        self.beam, self.target = self.kinematics.beam, self.kinematics.target
        self.r0 = r0
        z1, z2 = self.beam.Z, self.target.Z
        #: Charge of the nucleus whose field does the exciting.
        self.z_exciting = z1 if excite == "target" else z2
        m1, m2 = self.beam.nuclear_mass_mev, self.target.nuclear_mass_mev
        self.mu = m1 * m2 / (m1 + m2)
        self.k = z1 * z2 * E2_MEV_FM
        self.e_cm = self.kinematics.cm_energy_mev
        if self.e_cm <= self.e_star:
            raise ValueError("the CM energy is below the excitation energy")
        self.beta_i = math.sqrt(2 * self.e_cm / self.mu)
        self.beta_f = math.sqrt(2 * (self.e_cm - self.e_star) / self.mu)
        #: Symmetrised velocity (units of c), half-distance of closest approach (fm) and adiabaticity.
        self.beta = math.sqrt(self.beta_i * self.beta_f)
        self.a = self.k / (self.mu * self.beta_i * self.beta_f)
        self.xi = self.k / HBARC_MEV_FM * (1 / self.beta_f - 1 / self.beta_i)
        #: Sommerfeld parameter (symmetrised): the semiclassical picture needs η ≫ 1.
        self.eta = self.k / (HBARC_MEV_FM * self.beta)

    # -- excitation probability -------------------------------------------------------------------------------

    def _p(self, theta_cm: float) -> float:
        lam = self.lam
        s = math.sin(math.radians(theta_cm) / 2)
        if s <= 0:
            return 0.0
        integrals = orbit_integrals(lam, 1 / s, self.xi)
        ys = np.array([ylm_equator(lam, mu) for mu in range(-lam, lam + 1)])
        # P = (4π Z e² / (ħc (2λ+1)))² · B/(2λ+1) · Σ_μ |Y_λμ(π/2,0) I_μ / (β a^λ)|²
        pref = (4 * math.pi * self.z_exciting * E2_MEV_FM / (HBARC_MEV_FM * (2 * lam + 1))) ** 2
        amp2 = np.sum(np.abs(ys * integrals) ** 2) / (self.beta * self.a**lam) ** 2
        return pref * self.b_up / (2 * lam + 1) * amp2

    def probability(self, theta_cm: ArrayLike, exact: bool = False) -> ArrayLike:
        """First-order excitation probability for scattering to CM angle ``theta_cm`` (degrees).

        By default a four-point cubic interpolation in exact values every 1° (accurate to ~1e-5 relative, and fast
        for many angles); ``exact=True`` integrates the orbit for each angle."""
        th = np.atleast_1d(np.asarray(theta_cm, dtype=float))
        if exact:
            out = np.array([self._p(t) for t in th])
        else:
            grid, values = self._table()
            h = grid[1] - grid[0]
            t = np.clip(np.nan_to_num(th), 0.0, 180.0)
            i = np.clip(np.floor(t / h).astype(int) - 1, 0, len(grid) - 4)
            x = (t - grid[i]) / h  # position within the four points i … i+3, in grid steps
            y0, y1, y2, y3 = (values[i + k] for k in range(4))
            out = (-y0 * (x - 1) * (x - 2) * (x - 3) / 6 + y1 * x * (x - 2) * (x - 3) / 2
                   - y2 * x * (x - 1) * (x - 3) / 2 + y3 * x * (x - 1) * (x - 2) / 6)
            out = np.where((th <= 0) | np.isnan(th), 0.0, np.maximum(out, 0.0))
        return float(out[0]) if np.ndim(theta_cm) == 0 else out

    def _table(self) -> tuple:
        if getattr(self, "_grid", None) is None:
            grid = np.linspace(0.0, 180.0, 181)
            self._grid = (grid, np.array([self._p(t) for t in grid]))
        return self._grid

    def rutherford_cm(self, theta_cm: ArrayLike) -> ArrayLike:
        """The symmetrised Rutherford cross section (a²/4) / sin⁴(θ/2), mb/sr."""
        s = np.sin(np.radians(np.asarray(theta_cm, dtype=float)) / 2)
        with np.errstate(divide="ignore"):
            out = self.a**2 / 4 / s**4 * 10.0
        return float(out) if np.ndim(theta_cm) == 0 else out

    def cross_section_cm(self, theta_cm: ArrayLike) -> ArrayLike:
        """dσ/dΩ for exciting the state, CM frame, mb/sr: P(θ) × dσ_R/dΩ."""
        return np.asarray(self.probability(theta_cm)) * np.asarray(self.rutherford_cm(theta_cm)) \
            if np.ndim(theta_cm) else self.probability(theta_cm) * self.rutherford_cm(theta_cm)

    def cross_section_lab(self, theta_lab: ArrayLike, particle: str = "ejectile") -> tuple:
        """dσ/dΩ in the lab at the lab angle of the scattered beam (or the recoil), mb/sr: (first, second)
        solution, as :meth:`physim.nuclear.rutherford.Rutherford.cross_section_lab`."""
        out = []
        for point in self.kinematics.at_lab(theta_lab, particle):
            if point is None:
                out.append(None)
                continue
            th_ej = np.asarray(point.theta_cm) if particle == "ejectile" else 180.0 - np.asarray(point.theta_cm)
            th_ej = np.atleast_1d(th_ej)
            ok = ~np.isnan(th_ej)
            sig = np.full(th_ej.shape, np.nan)
            sig[ok] = np.atleast_1d(self.cross_section_cm(th_ej[ok]))
            val = sig * np.atleast_1d(point.jacobian)
            out.append(float(val[0]) if np.ndim(theta_lab) == 0 else val)
        return tuple(out)

    def total(self, theta_min: float = 1.0, theta_max: float = 180.0, points: int = 400) -> float:
        """∫ dσ/dΩ dΩ between two CM angles, mb (excitation dies off at small angles, so the default lower limit
        hardly matters)."""
        u = np.linspace(math.cos(math.radians(theta_min)), math.cos(math.radians(theta_max)), points)
        th = np.degrees(np.arccos(np.clip(u, -1, 1)))
        y = self.cross_section_cm(th)
        return float(2 * math.pi * abs(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(u))))

    # -- validity -------------------------------------------------------------------------------------------------

    @property
    def safe_distance(self) -> float:
        """Cline's safe distance 1.25 (A₁^⅓ + A₂^⅓) + 5 fm."""
        return self.r0 * (self.beam.A ** (1 / 3) + self.target.A ** (1 / 3)) + SAFE_MARGIN_FM

    def closest_approach(self, theta_cm: ArrayLike) -> ArrayLike:
        """Distance of closest approach on the symmetrised orbit, fm: a (1 + 1/sin(θ/2))."""
        s = np.sin(np.radians(np.asarray(theta_cm, dtype=float)) / 2)
        with np.errstate(divide="ignore"):
            out = self.a * (1 + 1 / s)
        return float(out) if np.ndim(theta_cm) == 0 else out

    def safe(self, theta_cm: ArrayLike):
        """True where the closest approach is at least the safe distance (pure Coulomb excitation)."""
        return np.asarray(self.closest_approach(theta_cm)) >= self.safe_distance

    def max_safe_angle(self) -> float:
        """Largest CM angle that still satisfies Cline's criterion, degrees (180 if all do)."""
        x = self.safe_distance / self.a - 1
        if x <= 1:
            return 180.0
        return math.degrees(2 * math.asin(1 / x))

    def warnings(self, theta_cm_max: Optional[float] = None) -> list:
        out = []
        limit = 180.0 if theta_cm_max is None else theta_cm_max
        safe = self.max_safe_angle()
        if safe < limit:
            out.append(f"Coulomb excitation: above {safe:.1f}° (CM) the nuclei come closer than the safe distance "
                       f"({self.safe_distance:.1f} fm): nuclear excitation interferes and first-order Coulomb "
                       f"excitation is not reliable there.")
        if self.eta < 5:
            out.append(f"Coulomb excitation: η = {self.eta:.1f} is not ≫ 1; the semiclassical approximation is "
                       f"rough.")
        p = self.probability(min(limit, 179.0))
        if p > 0.1:
            out.append(f"Coulomb excitation: the excitation probability reaches {p:.2f}; first order is not enough "
                       f"(multi-step excitation and reorientation matter: use GOSIA).")
        return out

    def __repr__(self) -> str:
        who = self.target if self.excite == "target" else self.beam
        return f"Coulex({who} {self.e_star:g} MeV E{self.lam}, {self.kinematics.label})"


__all__ = ["B_UNITS", "Coulex", "HBARC_MEV_FM", "orbit_integrals", "parse_b", "ylm_equator"]
