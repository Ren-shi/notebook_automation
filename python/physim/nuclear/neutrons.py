"""Neutron detectors and channel selection for fusion-evaporation.

A neutron detector here is a cylinder of organic scintillator (liquid or plastic) whose axis points at the target.
A neutron is seen when it scatters on hydrogen and the recoil proton makes more light than the threshold. The
module gives its efficiency against energy, the time of flight, the misidentification between neutrons and γ rays,
and the chance that one neutron fires two neighbouring cells (cross-talk)::

    from physim.nuclear import neutrons

    cell = neutrons.NeutronDetector()                 # a typical EJ-301 cell of a neutron wall
    cell.intrinsic_efficiency(2.0)                    # fraction of 2 MeV neutrons entering the face that are seen
    cell.efficiency(2.0)                              # the same, times the solid angle / 4π
    cell.time_of_flight(2.0)                          # ns from the target
    cell.crosstalk([130.0, 0.0, 0.0], 5.0)            # P(neighbour 130 mm away also fires | this cell fired)

    wall = neutrons.NeutronWall(cell, cells=50, pitch="130 mm")
    det = neutrons.Detection(neutron_efficiency=wall.efficiency(2.0), proton_efficiency=0.6,
                             alpha_efficiency=0.4, crosstalk=wall.crosstalk(2.0))
    channels = neutrons.channels_from({"2n": (2, 0, 0), "p2n": (2, 1, 0), "3n": (3, 0, 0)})
    m = neutrons.selection_matrix(channels, det)      # share of each true channel assigned to each channel
    m.share("3n", "2n")
    mc = neutrons.simulate_selection(events, channels, det, seed=1)   # the same from per-event particle lists

**The efficiency.** For a neutron of energy E entering the face along the axis of a cell of depth D,

    ε(E) = (Σ_H / Σ) (1 − exp(−Σ D)) (1 − E_p,th / E),    Σ = n_H σ_H(E) + n_C σ_C(E),

the probability of a first collision on hydrogen inside the cell, times the share of recoil protons above the
threshold (n–p scattering is isotropic in the CM below about 20 MeV, so the recoil energy is uniform between 0 and
E). E_p,th is the proton energy whose light equals the threshold (in MeVee), from the light-output curve of
Cecil et al. σ_H is Gammel's closed form; σ_C comes from ``data/neutron_scintillators.toml``. Carbon only
attenuates; double scattering, edge effects and the light of α particles from ¹²C(n,α) are left out.

**Channel selection.** Each true channel emits a fixed number of neutrons, protons, α particles (and γ rays). Each
particle is detected independently; its outcome is one of a few counts (a neutron: 0, 1, or 2 with cross-talk;
a proton: none, as a proton, or as an α). The detected counts are the sum of the particles' outcomes, so their
distribution is a convolution of per-particle distributions (multinomial). An event is assigned to the channel
whose multiplicities equal the detected ones, or left unassigned. The Monte Carlo draws the same outcomes particle
by particle from event lists, so it agrees with the analytic matrix within statistics when the same mean
probabilities are used (:meth:`Detection.averaged`).

The theory page is ``docs/theory/neutrons.md``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable, Mapping, Optional, Sequence, Union

import numpy as np

from . import _toml
from .quantity import Quantity

DATA = Path(__file__).resolve().parent / "data"

#: Neutron mass, MeV/c² (CODATA 2018).
NEUTRON_MEV = 939.56542052
#: Speed of light, mm/ns.
C_MM_PER_NS = 299.792458
#: 1 barn in cm².
BARN_CM2 = 1e-24
FWHM_PER_SIGMA = 2.3548200450309493

#: Particle kinds the event lists may use, and the name this module uses for each.
KINDS = {"n": "n", "neutron": "n", "p": "p", "proton": "p", "1H": "p", "a": "a", "alpha": "a", "α": "a",
         "4He": "a", "g": "g", "gamma": "g", "γ": "g"}


# -- data ---------------------------------------------------------------------------------------------------------


@lru_cache(maxsize=None)
def _table() -> dict:
    return _toml.loads((DATA / "neutron_scintillators.toml").read_text(encoding="utf-8"))


def materials() -> dict:
    """The scintillators in the data file, by name (``"EJ-301"``, ``"EJ-309"``, ``"EJ-200"``)."""
    return dict(_table()["materials"])


def material(name: str) -> dict:
    """One scintillator's entry. ``NE213`` and ``BC-501A`` read as EJ-301, ``NE102`` and ``BC-408`` as EJ-200."""
    alias = {"NE213": "EJ-301", "NE-213": "EJ-301", "BC-501A": "EJ-301", "BC501A": "EJ-301",
             "NE102": "EJ-200", "NE-102": "EJ-200", "BC-408": "EJ-200", "BC408": "EJ-200"}
    table = materials()
    key = alias.get(name, name)
    if key not in table:
        raise ValueError(f"unknown scintillator '{name}'; known: {', '.join(table)} (and NE213, BC-501A, NE102)")
    return table[key]


def hydrogen_cross_section(energy_mev):
    """n–p total (elastic) cross section in barns, Gammel's closed form, E in MeV (lab), good to ~1% up to ~40 MeV.

    J. L. Gammel, in *Fast Neutron Physics*, Part II, ed. J. B. Marion and J. L. Fowler (Interscience, 1963),
    p. 2185; quoted in G. F. Knoll, *Radiation Detection and Measurement*, 4th ed. (Wiley, 2010), chapter 15.
    """
    e = np.asarray(energy_mev, dtype=float)
    first = 3 * math.pi / (1.206 * e + (-1.86 + 0.09415 * e + 0.0001306 * e * e) ** 2)
    second = math.pi / (1.206 * e + (0.4223 + 0.13 * e) ** 2)
    return first + second


def carbon_cross_section(energy_mev):
    """Total neutron cross section of natural carbon in barns, interpolated in the data file (typical values)."""
    t = _table()["carbon_total"]
    e = np.asarray(energy_mev, dtype=float)
    return np.interp(np.log(np.clip(e, 1e-6, None)), np.log(t["energy_mev"]), t["sigma_b"])


# -- light output -------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class LightOutput:
    """Proton light output L(E) = a1 E − a2 (1 − exp(−a3 E^a4)), E in MeV, L in MeVee (Cecil et al. 1979)."""

    a1: float
    a2: float
    a3: float
    a4: float

    def __call__(self, proton_mev):
        e = np.asarray(proton_mev, dtype=float)
        return self.a1 * e - self.a2 * (1 - np.exp(-self.a3 * np.clip(e, 0, None) ** self.a4))

    def proton_energy(self, light_mevee: float) -> float:
        """The proton energy, MeV, whose light is ``light_mevee`` (L rises monotonically, so bisection)."""
        if light_mevee <= 0:
            return 0.0
        lo, hi = 0.0, 1.0
        while float(self(hi)) < light_mevee:
            hi *= 2
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if float(self(mid)) < light_mevee else (lo, mid)
        return 0.5 * (lo + hi)


# -- one detector -------------------------------------------------------------------------------------------------


def _quantity(text, unit: str) -> float:
    if isinstance(text, (int, float)):
        return float(text)
    q = text if isinstance(text, Quantity) else Quantity.parse(text)
    return q.to(unit)


def _fraction(x) -> float:
    """A probability given as a number in [0, 1] or as text with %, e.g. ``"2 %"``."""
    v = float(x) if isinstance(x, (int, float)) else _quantity(x, "%") / 100.0
    if not 0.0 <= v <= 1.0:
        raise ValueError(f"a probability must lie between 0 and 1 (0 % and 100 %), got {x!r}")
    return v


#: The fields of :class:`NeutronDetector` whose defaults are typical values, not from a data sheet.
TYPICAL_FIELDS = ("radius", "depth", "distance", "threshold", "time_resolution", "neutron_as_gamma",
                  "gamma_as_neutron", "gamma_efficiency")


@dataclass(frozen=True)
class NeutronDetector:
    """A cylindrical cell of organic scintillator facing the target along its axis.

    Fields are text with units, as in a setup file. The defaults describe a typical EJ-301 (NE213) cell of a
    neutron wall: 127 mm across, 150 mm deep, 510 mm from the target, a threshold of 100 keVee. All defaults
    except ``material`` are typical (:data:`TYPICAL_FIELDS`), to be replaced by the values of the detector used.

    ``threshold`` is the light threshold in electron-equivalent energy: ``"0.1 MeV"`` means 100 keVee.
    ``time_resolution`` is the FWHM of the time-of-flight measurement. ``neutron_as_gamma`` is the probability
    that a detected neutron is called a γ ray by pulse-shape discrimination, ``gamma_as_neutron`` the reverse;
    ``gamma_efficiency`` the probability that a γ ray entering the face gives light above the threshold.
    """

    name: str = "N1"
    material: str = "EJ-301"
    radius: str = "63.5 mm"
    depth: str = "150 mm"
    distance: str = "510 mm"
    threshold: str = "0.1 MeV"
    time_resolution: str = "1.5 ns"
    neutron_as_gamma: str = "2 %"
    gamma_as_neutron: str = "1 %"
    gamma_efficiency: str = "40 %"

    # -- plain numbers in the canonical units ---------------------------------------------------------------------

    @property
    def radius_mm(self) -> float:
        return _quantity(self.radius, "mm")

    @property
    def depth_mm(self) -> float:
        return _quantity(self.depth, "mm")

    @property
    def distance_mm(self) -> float:
        """Distance from the target to the front face, mm."""
        return _quantity(self.distance, "mm")

    @property
    def threshold_mevee(self) -> float:
        return _quantity(self.threshold, "MeV")

    @property
    def time_resolution_ns(self) -> float:
        """FWHM, ns."""
        return _quantity(self.time_resolution, "ns")

    @property
    def misidentification(self) -> tuple:
        """(P(neutron called γ), P(γ called neutron))."""
        return _fraction(self.neutron_as_gamma), _fraction(self.gamma_as_neutron)

    @property
    def scintillator(self) -> dict:
        return material(self.material)

    @property
    def light_output(self) -> LightOutput:
        return LightOutput(**self.scintillator["light_output"])

    def check(self) -> list:
        """Problems with the fields, as readable sentences (empty when the detector is valid)."""
        problems = []
        for name in ("radius", "depth", "distance", "time_resolution"):
            try:
                v = getattr(self, f"{name}_mm" if name != "time_resolution" else "time_resolution_ns")
                if v <= 0:
                    problems.append(f"{self.name}: {name} must be positive, got {getattr(self, name)!r}")
            except ValueError as err:
                problems.append(f"{self.name}: {name}: {err}")
        try:
            if self.threshold_mevee < 0:
                problems.append(f"{self.name}: the threshold cannot be negative")
        except ValueError as err:
            problems.append(f"{self.name}: threshold: {err}")
        for name in ("neutron_as_gamma", "gamma_as_neutron", "gamma_efficiency"):
            try:
                _fraction(getattr(self, name))
            except ValueError as err:
                problems.append(f"{self.name}: {name}: {err}")
        try:
            self.scintillator
        except ValueError as err:
            problems.append(f"{self.name}: {err}")
        return problems

    # -- efficiency -----------------------------------------------------------------------------------------------

    def proton_threshold(self) -> float:
        """The recoil-proton energy, MeV, whose light equals the threshold."""
        return self.light_output.proton_energy(self.threshold_mevee)

    def macroscopic(self, energy_mev) -> tuple:
        """(Σ_H, Σ_C) in 1/mm: number density times cross section for hydrogen and carbon."""
        s = self.scintillator
        to_per_mm = BARN_CM2 * 0.1  # atoms/cm³ × cm² = 1/cm, then 1/cm → 1/mm
        return (s["hydrogen"] * hydrogen_cross_section(energy_mev) * to_per_mm,
                s["carbon"] * carbon_cross_section(energy_mev) * to_per_mm)

    def intrinsic_efficiency(self, energy_mev, path_mm: Optional[float] = None, carbon: bool = True):
        """Probability that a neutron of this energy entering the face along the axis is seen.

        ``path_mm`` overrides the depth (a chord through the cell); ``carbon=False`` drops the carbon attenuation.
        """
        e = np.asarray(energy_mev, dtype=float)
        length = self.depth_mm if path_mm is None else path_mm
        sh, sc = self.macroscopic(e)
        total = sh + (sc if carbon else 0.0)
        e_th = self.proton_threshold()
        above = np.clip(1.0 - e_th / np.where(e > 0, e, np.inf), 0.0, None)
        eff = sh / total * -np.expm1(-total * length) * above
        return float(eff) if eff.ndim == 0 else eff

    def solid_angle(self) -> float:
        """Solid angle of the front face seen from the target, msr (a disc on its axis)."""
        d, r = self.distance_mm, self.radius_mm
        return 2 * math.pi * (1 - d / math.hypot(d, r)) * 1e3

    def efficiency(self, energy_mev):
        """Probability that a neutron emitted isotropically from the target is seen by this cell."""
        return self.intrinsic_efficiency(energy_mev) * self.solid_angle() / (4e3 * math.pi)

    # -- time of flight -------------------------------------------------------------------------------------------

    def flight_path_mm(self) -> float:
        """Mean flight path to the middle of the cell, mm."""
        return self.distance_mm + 0.5 * self.depth_mm

    def time_of_flight(self, energy_mev, path_mm: Optional[float] = None):
        """Flight time from the target, ns, for a neutron of kinetic energy E (relativistic)."""
        return time_of_flight(energy_mev, self.flight_path_mm() if path_mm is None else path_mm)

    def gamma_flash(self) -> float:
        """Flight time of a γ ray to the middle of the cell, ns."""
        return self.flight_path_mm() / C_MM_PER_NS

    def energy_from_tof(self, tof_ns, path_mm: Optional[float] = None):
        return energy_from_tof(tof_ns, self.flight_path_mm() if path_mm is None else path_mm)

    def energy_resolution(self, energy_mev):
        """FWHM of the energy measured by time of flight, MeV.

        σ_E / E = γ(γ+1) √((σ_t/t)² + (σ_L/L)²), with σ_t from ``time_resolution`` and σ_L = depth/√12 (the
        unknown depth of the interaction).
        """
        e = np.asarray(energy_mev, dtype=float)
        g = 1 + e / NEUTRON_MEV
        t = self.time_of_flight(e)
        sigma_t = self.time_resolution_ns / FWHM_PER_SIGMA
        sigma_l = self.depth_mm / math.sqrt(12)
        rel = g * (g + 1) * np.sqrt((sigma_t / t) ** 2 + (sigma_l / self.flight_path_mm()) ** 2)
        out = FWHM_PER_SIGMA * rel * e
        return float(out) if np.ndim(out) == 0 else out

    # -- cross-talk -----------------------------------------------------------------------------------------------

    def crosstalk(self, offset_mm, energy_mev: float, samples: int = 200_000, seed: int = 1) -> float:
        """P(a neighbouring cell also fires | this cell fired), for a neutron of this energy along the axis.

        ``offset_mm`` is the neighbour's centre relative to this cell's centre, in the cell's frame (z along the
        axis, away from the target); the neighbour is an identical, parallel cell. A seeded Monte Carlo: the first
        collision is an n–p scatter at a depth drawn from the attenuation, at a random point of the face; the
        neutron leaves at lab angle θ with cos²θ uniform (isotropic CM), keeps E cos²θ, and gives the proton
        E sin²θ, which must be above the threshold. The scattered neutron escapes this cell with
        exp(−Σ ℓ₁) and fires the neighbour with the single-scatter efficiency along its chord ℓ₂ there.
        Scattering on carbon and second scatters are left out.
        """
        return self.crosstalk_all([offset_mm], energy_mev, samples, seed)

    def crosstalk_all(self, offsets_mm: Sequence, energy_mev: float, samples: int = 200_000, seed: int = 1) -> float:
        """As :meth:`crosstalk`, summed over several neighbours (each scattered neutron reaches at most one)."""
        e = float(energy_mev)
        e_th = self.proton_threshold()
        if e <= 2 * e_th:
            return 0.0  # the proton and the scattered neutron cannot both be above threshold
        rng = np.random.default_rng(seed)
        r_cell, depth = self.radius_mm, self.depth_mm
        sh, sc = self.macroscopic(e)
        total = sh + sc
        # first collision: depth from the truncated exponential, position uniform on the face
        u = rng.random(samples)
        z = -np.log1p(-u * -np.expm1(-total * depth)) / total - depth / 2
        rho, phi = r_cell * np.sqrt(rng.random(samples)), 2 * np.pi * rng.random(samples)
        p0 = np.stack([rho * np.cos(phi), rho * np.sin(phi), z], axis=1)
        # scatter: cos²θ = x uniform; the proton energy is e (1 − x)
        x = rng.random(samples) * (1 - e_th / e)  # conditioned on the proton being above threshold
        cos_t, sin_t = np.sqrt(x), np.sqrt(1 - x)
        az = 2 * np.pi * rng.random(samples)
        d = np.stack([sin_t * np.cos(az), sin_t * np.sin(az), cos_t], axis=1)
        e_out = e * x
        sh_out, sc_out = self.macroscopic(e_out)
        tot_out = sh_out + sc_out
        _, exit1 = _cylinder_chord(p0, d, np.zeros(3), r_cell, depth)
        escape = np.exp(-tot_out * exit1)
        fired = np.zeros(samples)
        above = np.clip(1 - e_th / np.where(e_out > 0, e_out, np.inf), 0, None)
        for off in offsets_mm:
            t_in, t_out = _cylinder_chord(p0, d, np.asarray(off, dtype=float), r_cell, depth)
            chord = np.clip(t_out - np.maximum(t_in, 0.0), 0.0, None)
            fired += sh_out / tot_out * -np.expm1(-tot_out * chord) * above
        return float(np.mean(escape * fired))


def _cylinder_chord(p, d, centre, radius, depth) -> tuple:
    """Entry and exit distances along rays p + t d through a cylinder (axis ∥ z, given centre, radius, depth).

    Rays that miss get t_in = t_out = 0. For a ray starting inside, t_in < 0 and t_out is the distance to the exit.
    """
    q = p - centre
    a = d[:, 0] ** 2 + d[:, 1] ** 2
    b = 2 * (q[:, 0] * d[:, 0] + q[:, 1] * d[:, 1])
    c = q[:, 0] ** 2 + q[:, 1] ** 2 - radius ** 2
    disc = b * b - 4 * a * c
    ok = (disc > 0) & (a > 1e-15)
    sq = np.sqrt(np.where(ok, disc, 0.0))
    safe_a = np.where(ok, a, 1.0)
    t1 = np.where(ok, (-b - sq) / (2 * safe_a), -np.inf)
    t2 = np.where(ok, (-b + sq) / (2 * safe_a), np.inf)
    # rays parallel to the axis: inside the circle for all t, or never
    parallel_inside = (a <= 1e-15) & (c < 0)
    t1 = np.where(parallel_inside, -np.inf, t1)
    t2 = np.where(parallel_inside, np.inf, t2)
    never = ~ok & ~parallel_inside
    dz = d[:, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        za = (-depth / 2 - q[:, 2]) / dz
        zb = (depth / 2 - q[:, 2]) / dz
    flat = np.abs(dz) < 1e-15
    inside_slab = np.abs(q[:, 2]) <= depth / 2
    z1 = np.where(flat, np.where(inside_slab, -np.inf, np.inf), np.minimum(za, zb))
    z2 = np.where(flat, np.where(inside_slab, np.inf, -np.inf), np.maximum(za, zb))
    t_in, t_out = np.maximum(t1, z1), np.minimum(t2, z2)
    miss = never | (t_out <= t_in) | (t_out <= 0)
    return np.where(miss, 0.0, t_in), np.where(miss, 0.0, t_out)


def time_of_flight(energy_mev, path_mm: float):
    """Flight time, ns, of a neutron of kinetic energy E over ``path_mm`` (relativistic)."""
    e = np.asarray(energy_mev, dtype=float)
    g = 1 + e / NEUTRON_MEV
    beta = np.sqrt(1 - 1 / (g * g))
    out = path_mm / (beta * C_MM_PER_NS)
    return float(out) if out.ndim == 0 else out


def energy_from_tof(tof_ns, path_mm: float):
    """Kinetic energy, MeV, of a neutron that takes ``tof_ns`` over ``path_mm``."""
    beta = path_mm / (np.asarray(tof_ns, dtype=float) * C_MM_PER_NS)
    if np.any(beta >= 1):
        raise ValueError("a flight time shorter than light's over the same path is not a neutron")
    out = NEUTRON_MEV * (1 / np.sqrt(1 - beta * beta) - 1)
    return float(out) if np.ndim(out) == 0 else out


def misidentification_from_fom(fom: float) -> float:
    """The probability each way of calling a neutron a γ ray (and back) for a pulse-shape figure of merit.

    FOM = separation / (FWHM_n + FWHM_γ). For two Gaussian peaks of equal width cut halfway between them, each
    peak's tail beyond the cut is ½ erfc(FOM × 2.3548 / √2).
    """
    return 0.5 * math.erfc(fom * FWHM_PER_SIGMA / math.sqrt(2))


# -- an array of cells --------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class NeutronWall:
    """Identical cells at the same distance, packed side by side on a sphere around the target.

    ``layout`` "hexagonal" gives each inner cell 6 neighbours at ``pitch``; "square" gives 4 at ``pitch`` and 4
    at √2 × pitch. Cells on the edge of the wall have fewer; ``edge_share`` (typical: 0) is the share of
    neighbours missing on average. For event lists the cells are taken to fill a cone around the beam from
    ``theta_min`` to the :meth:`theta_max` whose solid angle equals the cells' total, with no gaps.
    """

    cell: NeutronDetector = field(default_factory=NeutronDetector)
    cells: int = 50
    pitch: str = "130 mm"
    layout: str = "hexagonal"
    edge_share: float = 0.0
    theta_min: str = "0 deg"

    @property
    def pitch_mm(self) -> float:
        return _quantity(self.pitch, "mm")

    def solid_angle(self) -> float:
        """Total solid angle of the cells, msr."""
        return self.cells * self.cell.solid_angle()

    def geometric_efficiency(self) -> float:
        return self.solid_angle() / (4e3 * math.pi)

    def theta_max(self) -> float:
        """Polar angle, deg, to which a cone from ``theta_min`` must reach to hold the cells' solid angle."""
        c0 = math.cos(math.radians(_quantity(self.theta_min, "deg")))
        c1 = c0 - self.solid_angle() / (2e3 * math.pi)
        if c1 < -1:
            raise ValueError("the cells cover more than the whole sphere beyond theta_min")
        return math.degrees(math.acos(c1))

    def efficiency(self, energy_mev):
        """Probability that a neutron emitted isotropically is seen by one of the cells."""
        return self.cell.intrinsic_efficiency(energy_mev) * self.geometric_efficiency()

    def neighbour_offsets(self) -> list:
        """Centres of an inner cell's neighbours, mm, in the cell's frame (z towards the back)."""
        p = self.pitch_mm
        if self.layout == "hexagonal":
            return [[p * math.cos(k * math.pi / 3), p * math.sin(k * math.pi / 3), 0.0] for k in range(6)]
        if self.layout == "square":
            return [[p, 0, 0], [-p, 0, 0], [0, p, 0], [0, -p, 0], [p, p, 0], [p, -p, 0], [-p, p, 0], [-p, -p, 0]]
        raise ValueError(f"unknown layout '{self.layout}'; use 'hexagonal' or 'square'")

    def crosstalk(self, energy_mev: float, samples: int = 200_000, seed: int = 1) -> float:
        """P(a second cell fires | one cell fired), for a neutron of this energy (inner cell, less the edge share)."""
        return (1 - self.edge_share) * self.cell.crosstalk_all(self.neighbour_offsets(), energy_mev, samples, seed)

    def detection_probability(self, energy_mev, direction) -> np.ndarray:
        """Per-particle probability for event lists: intrinsic efficiency × packing inside the cone, 0 outside."""
        d = np.atleast_2d(np.asarray(direction, dtype=float))
        cos_t = d[:, 2] / np.linalg.norm(d, axis=1)
        c0 = math.cos(math.radians(_quantity(self.theta_min, "deg")))
        c1 = math.cos(math.radians(self.theta_max()))
        inside = (cos_t <= c0) & (cos_t >= c1)
        return np.where(inside, np.asarray(self.cell.intrinsic_efficiency(energy_mev)), 0.0)


# -- channels -----------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Channel:
    """An exit channel: how many neutrons, protons, α particles and γ rays each event emits."""

    name: str
    neutrons: int = 0
    protons: int = 0
    alphas: int = 0
    gammas: int = 0

    @property
    def multiplicities(self) -> tuple:
        return (self.neutrons, self.protons, self.alphas)


def channel_name(neutrons: int, protons: int, alphas: int) -> str:
    """The usual short name: α first, then p, then n, e.g. ``"αp2n"``; ``"γ"`` for none."""
    def part(k, s):
        return "" if k == 0 else (s if k == 1 else f"{k}{s}")
    return (part(alphas, "α") + part(protons, "p") + part(neutrons, "n")) or "γ"


def channels_from(spec: Union[Mapping, Iterable], gammas: Union[int, Mapping, None] = None) -> list:
    """Channels from ``{"p2n": (2, 1, 0), ...}`` (neutrons, protons, alphas), names, or :class:`Channel` objects.

    A bare name like ``"αp2n"`` is read back into its multiplicities. ``gammas`` gives every channel (an int) or
    each channel (a mapping) its γ-ray multiplicity.
    """
    out = []
    items = spec.items() if isinstance(spec, Mapping) else ((s, None) for s in spec)
    for key, value in items:
        if isinstance(key, Channel):
            ch = key
        elif value is None:
            ch = Channel(key, *parse_channel_name(key))
        else:
            ch = Channel(key, *value)
        g = gammas.get(ch.name, 0) if isinstance(gammas, Mapping) else gammas
        out.append(replace(ch, gammas=int(g)) if g is not None else ch)
    names = [c.name for c in out]
    if len(set(names)) != len(names):
        raise ValueError("two channels have the same name")
    mult = [c.multiplicities for c in out]
    if len(set(mult)) != len(mult):
        raise ValueError("two channels emit the same particles; detected counts cannot tell them apart")
    return out


def parse_channel_name(name: str) -> tuple:
    """(neutrons, protons, alphas) from a name like ``"2pn"``, ``"αp2n"``, ``"a3n"``."""
    counts = {"n": 0, "p": 0, "a": 0}
    i, s = 0, name.replace("α", "a")
    if s == "γ":
        return (0, 0, 0)
    while i < len(s):
        j = i
        while j < len(s) and s[j].isdigit():
            j += 1
        k = int(s[i:j]) if j > i else 1
        if j >= len(s) or s[j] not in counts:
            raise ValueError(f"cannot read the channel name '{name}'; write it like 'p2n' or 'αpn'")
        counts[s[j]] += k
        i = j + 1
    return counts["n"], counts["p"], counts["a"]


# -- detection probabilities --------------------------------------------------------------------------------------

Prob = Union[float, Callable]


@dataclass(frozen=True)
class Detection:
    """Per-particle detection and identification probabilities.

    ``neutron_efficiency``, ``proton_efficiency`` and ``alpha_efficiency``: the probability that one emitted
    particle gives a signal, either a number or a function ``f(energy_mev, direction) -> probability`` applied to
    each particle of an event list (for example :meth:`NeutronWall.detection_probability`). The charged-particle
    numbers come from the silicon or CsI array. ``crosstalk``: the probability that a detected neutron fires a
    second cell, a number or ``f(energy_mev)``. ``gamma_efficiency``: the probability that one γ ray of the event
    fires a neutron cell. Misidentification: ``neutron_as_gamma`` (lost neutron), ``gamma_as_neutron`` (a γ ray
    counted as a neutron), ``proton_as_alpha`` and ``alpha_as_proton``.
    """

    neutron_efficiency: Prob = 0.0
    proton_efficiency: Prob = 0.0
    alpha_efficiency: Prob = 0.0
    crosstalk: Prob = 0.0
    gamma_efficiency: float = 0.0
    neutron_as_gamma: float = 0.0
    gamma_as_neutron: float = 0.0
    proton_as_alpha: float = 0.0
    alpha_as_proton: float = 0.0

    @classmethod
    def from_detector(cls, detector: NeutronDetector, neutron_efficiency: Prob, crosstalk: Prob = 0.0,
                      **charged) -> Detection:
        """With the neutron–γ misidentification of a :class:`NeutronDetector`. ``gamma_efficiency`` (the chance
        that one γ ray of the event fires any cell) must be given in ``charged``, for example the cell's intrinsic
        γ efficiency times the wall's geometric efficiency."""
        n_as_g, g_as_n = detector.misidentification
        args = dict(neutron_efficiency=neutron_efficiency, crosstalk=crosstalk, neutron_as_gamma=n_as_g,
                    gamma_as_neutron=g_as_n)
        args.update(charged)
        return cls(**args)

    def is_numeric(self) -> bool:
        return not any(callable(getattr(self, f)) for f in
                       ("neutron_efficiency", "proton_efficiency", "alpha_efficiency", "crosstalk"))

    # per-particle outcome probabilities, used by both the analytic matrix and the Monte Carlo

    def neutron_outcomes(self, energy=None, direction=None) -> np.ndarray:
        """Probabilities of 0, 1 and 2 counted neutrons for one emitted neutron (rows for arrays of particles)."""
        eps = _evaluate(self.neutron_efficiency, energy, direction)
        c = _evaluate(self.crosstalk, energy)
        q, c = np.broadcast_arrays(eps * (1 - self.neutron_as_gamma), c)
        return np.stack([1 - q, q * (1 - c), q * c], axis=-1)

    def proton_outcomes(self, energy=None, direction=None) -> np.ndarray:
        """Probabilities of (nothing, a proton, an α) for one emitted proton."""
        eps = _evaluate(self.proton_efficiency, energy, direction)
        return np.stack([1 - eps, eps * (1 - self.proton_as_alpha), eps * self.proton_as_alpha], axis=-1)

    def alpha_outcomes(self, energy=None, direction=None) -> np.ndarray:
        """Probabilities of (nothing, a proton, an α) for one emitted α particle."""
        eps = _evaluate(self.alpha_efficiency, energy, direction)
        return np.stack([1 - eps, eps * self.alpha_as_proton, eps * (1 - self.alpha_as_proton)], axis=-1)

    def _outcomes(self, kind: str, energy, direction) -> np.ndarray:
        """Outcome probabilities, one row per particle, for arrays of particles of one kind."""
        method = {"n": self.neutron_outcomes, "p": self.proton_outcomes, "a": self.alpha_outcomes}[kind]
        return np.broadcast_to(method(energy, direction), (len(energy), 3))

    def gamma_fake(self) -> float:
        """Probability that one γ ray of the event is counted as a neutron."""
        return self.gamma_efficiency * self.gamma_as_neutron

    def averaged(self, events: Iterable, channels: Sequence[Channel]) -> dict:
        """Per-channel numeric :class:`Detection` with the mean outcome probabilities of the events' particles.

        With these the analytic matrix is exactly what the Monte Carlo over the same events converges to, as long
        as particles in one event are detected independently.
        """
        channels = list(channels)
        flat = _flatten(events, channels)
        mean = {}
        for kind in ("n", "p", "a"):
            sel = flat.kinds == kind
            probs = self._outcomes(kind, flat.energies[sel], flat.directions[sel])
            rows = flat.truth[flat.event[sel]]
            sums = np.zeros((len(channels), 3))
            np.add.at(sums, rows, probs)
            counts = np.bincount(rows, minlength=len(channels))[:, None]
            mean[kind] = np.where(counts > 0, sums / np.maximum(counts, 1), [1.0, 0.0, 0.0])
        out = {}
        for i, c in enumerate(channels):
            n, p, a = mean["n"][i], mean["p"][i], mean["a"][i]
            q = n[1] + n[2]
            eps_n = q / (1 - self.neutron_as_gamma) if self.neutron_as_gamma < 1 else 0.0
            eps_p, eps_a = p[1] + p[2], a[1] + a[2]
            out[c.name] = replace(
                self, neutron_efficiency=float(eps_n), crosstalk=float(n[2] / q) if q > 0 else 0.0,
                proton_efficiency=float(eps_p), alpha_efficiency=float(eps_a),
                proton_as_alpha=float(p[2] / eps_p) if eps_p > 0 else self.proton_as_alpha,
                alpha_as_proton=float(a[1] / eps_a) if eps_a > 0 else self.alpha_as_proton)
        return out


def _evaluate(prob, energy=None, direction=None):
    if callable(prob):
        try:
            v = prob(energy, direction)
        except TypeError:
            v = prob(energy)
        return np.asarray(v, dtype=float)
    return np.asarray(float(prob))


# -- the selection matrix -----------------------------------------------------------------------------------------

UNASSIGNED = "unassigned"


@dataclass
class SelectionMatrix:
    """Rows: true channels. Columns: the channel each event is assigned to, then ``"unassigned"``.

    ``matrix[i, j]`` is the share of events of true channel i assigned to channel j; each row sums to one.
    ``events`` holds the number of events per row for a Monte Carlo matrix (None for the analytic one).
    """

    true: list
    assigned: list
    matrix: np.ndarray
    events: Optional[np.ndarray] = None

    def share(self, true: str, assigned: str) -> float:
        return float(self.matrix[self.true.index(true), self.assigned.index(assigned)])

    def row(self, true: str) -> dict:
        return dict(zip(self.assigned, self.matrix[self.true.index(true)].tolist()))

    def uncertainty(self) -> np.ndarray:
        """Binomial standard error of each share (Monte Carlo matrices only)."""
        if self.events is None:
            raise ValueError("an analytic matrix has no statistical uncertainty")
        n = np.maximum(self.events, 1)[:, None]
        return np.sqrt(self.matrix * (1 - self.matrix) / n)

    def purity(self, cross_sections: Mapping) -> dict:
        """For each assigned channel, the share of its events that truly belong to it, given the true channels'
        relative cross sections (any unit, missing channels count as zero)."""
        w = np.array([float(cross_sections.get(t, 0.0)) for t in self.true])
        counts = w @ self.matrix
        out = {}
        for j, name in enumerate(self.assigned[:-1]):
            i = self.true.index(name) if name in self.true else None
            out[name] = float(w[i] * self.matrix[i, j] / counts[j]) if i is not None and counts[j] > 0 else 0.0
        return out

    def table(self, digits: int = 3) -> str:
        """Plain-text table of the shares."""
        width = max(10, max(len(a) for a in self.assigned) + 2)
        head = "true \\ assigned".ljust(16) + "".join(a.rjust(width) for a in self.assigned)
        lines = [head]
        for t, row in zip(self.true, self.matrix):
            lines.append(t.ljust(16) + "".join(f"{v:.{digits}f}".rjust(width) for v in row))
        return "\n".join(lines)


def _convolve_power(dist: np.ndarray, k: int) -> np.ndarray:
    out = np.array([1.0])
    for _ in range(k):
        out = np.convolve(out, dist)
    return out


def _charged_joint(proton: np.ndarray, alpha: np.ndarray, n_p: int, n_a: int) -> np.ndarray:
    """Joint distribution of (detected protons, detected α) as a 2D array, from per-particle outcomes."""
    def as_2d(o):
        a = np.zeros((2, 2))
        a[0, 0], a[1, 0], a[0, 1] = o
        return a
    joint = np.ones((1, 1))
    for o, k in ((proton, n_p), (alpha, n_a)):
        kernel = as_2d(o)
        for _ in range(k):
            new = np.zeros((joint.shape[0] + 1, joint.shape[1] + 1))
            for (i, j), w in np.ndenumerate(kernel):
                if w:
                    new[i:i + joint.shape[0], j:j + joint.shape[1]] += w * joint
            joint = new
    return joint


def selection_matrix(channels: Sequence[Channel], detection: Union[Detection, Mapping]) -> SelectionMatrix:
    """The analytic channel-selection matrix.

    ``detection`` is one numeric :class:`Detection` for all channels, or a mapping from channel name to one (for
    example from :meth:`Detection.averaged`). Each true neutron gives 0, 1 or 2 counted neutrons, each γ ray a
    fake neutron with :meth:`Detection.gamma_fake`, each proton and α particle (none, p, α); the counts are
    independent between particles, so their distributions convolve.
    """
    channels = list(channels)
    index = _channel_index(channels)
    names = [c.name for c in channels]
    columns = names + [UNASSIGNED]
    m = np.zeros((len(channels), len(columns)))
    for i, ch in enumerate(channels):
        det = detection[ch.name] if isinstance(detection, Mapping) else detection
        if not det.is_numeric():
            raise ValueError("the analytic matrix needs numbers; use Detection.averaged(events, channels) first")
        n_dist = _convolve_power(det.neutron_outcomes(), ch.neutrons)
        n_dist = np.convolve(n_dist, _convolve_power(np.array([1 - det.gamma_fake(), det.gamma_fake()]),
                                                     ch.gammas))
        pa = _charged_joint(det.proton_outcomes(), det.alpha_outcomes(), ch.protons, ch.alphas)
        for k, pn in enumerate(n_dist):
            if pn == 0:
                continue
            for (jp, ja), w in np.ndenumerate(pa):
                if w == 0:
                    continue
                target = index.get((k, jp, ja))
                m[i, columns.index(target) if target is not None else -1] += pn * w
    return SelectionMatrix(names, columns, m)


# -- the Monte Carlo over event lists -----------------------------------------------------------------------------


def _particle_fields(p) -> tuple:
    """(kind, energy MeV, direction) from a tuple, a dict or an object with attributes."""
    def get(obj, keys, default=None):
        for k in keys:
            if isinstance(obj, Mapping) and k in obj:
                return obj[k]
            if not isinstance(obj, Mapping) and hasattr(obj, k):
                return getattr(obj, k)
        return default
    if isinstance(p, (tuple, list)):
        kind, energy, direction = (list(p) + [None, None])[:3]
    else:
        kind = get(p, ("kind", "type", "particle"))
        energy = get(p, ("energy_mev", "energy", "e_lab", "lab_energy", "e"))
        direction = get(p, ("direction", "lab_direction", "dir"))
    if kind not in KINDS:
        raise ValueError(f"unknown particle kind {kind!r}; use one of {', '.join(sorted(set(KINDS)))}")
    return KINDS[kind], (None if energy is None else float(energy)), direction


def _event_fields(ev) -> tuple:
    """(particles as (kind, energy, direction), extra γ-ray count) from one event."""
    if isinstance(ev, Mapping) or hasattr(ev, "particles"):
        parts = ev["particles"] if isinstance(ev, Mapping) else ev.particles
        g = (ev.get("gammas", 0) if isinstance(ev, Mapping) else getattr(ev, "gammas", 0)) or 0
    else:
        parts, g = ev, 0
    fields = [_particle_fields(p) for p in parts]
    if isinstance(g, (list, tuple)):
        g = len(g)
    return fields, int(g)


def _counts(parts) -> tuple:
    kinds = [k for k, _, _ in parts]
    return kinds.count("n"), kinds.count("p"), kinds.count("a"), kinds.count("g")


def _channel_index(channels: Sequence[Channel]) -> dict:
    return {c.multiplicities: c.name for c in channels}


@dataclass
class _Flat:
    """Event lists flattened to one row per charged particle or neutron (γ rays only counted per event)."""

    event: np.ndarray       # event number of each particle
    kinds: np.ndarray       # "n", "p" or "a"
    energies: np.ndarray    # MeV (NaN when not given)
    directions: np.ndarray  # (N, 3) (NaN when not given)
    truth: np.ndarray       # per event: row of the true channel
    gammas: np.ndarray      # per event: γ-ray count


def _flatten(events: Iterable, channels: Sequence[Channel]) -> _Flat:
    index = {c.multiplicities: i for i, c in enumerate(channels)}
    ev_rows, kinds, energies, dirs, truth, gammas = [], [], [], [], [], []
    for ev in events:
        parts, g_extra = _event_fields(ev)
        counts = _counts(parts)
        row = index.get(counts[:3])
        if row is None:
            continue
        number = len(truth)
        for kind, energy, direction in parts:
            if kind == "g":
                continue
            ev_rows.append(number)
            kinds.append(kind)
            energies.append(np.nan if energy is None else energy)
            dirs.append((np.nan,) * 3 if direction is None else tuple(np.asarray(direction, dtype=float)))
        truth.append(row)
        gammas.append(counts[3] + g_extra)
    return _Flat(np.asarray(ev_rows, dtype=int), np.asarray(kinds, dtype="<U1"), np.asarray(energies, dtype=float),
                 np.asarray(dirs, dtype=float).reshape(-1, 3), np.asarray(truth, dtype=int),
                 np.asarray(gammas, dtype=int))


def simulate_selection(events: Iterable, channels: Sequence[Channel], detection: Detection,
                       seed: int = 1, rng: Optional[np.random.Generator] = None) -> SelectionMatrix:
    """The channel-selection matrix from a seeded Monte Carlo over event lists.

    Each event is a sequence of particles, or a mapping/object with ``particles`` (and optionally ``gammas``, a
    count or a list). A particle is ``(kind, energy_mev, direction)``, a mapping, or an object with ``kind``
    (``"n"``, ``"p"``, ``"alpha"``, ``"gamma"``, ...), ``energy_mev`` (lab) and ``direction`` (lab, a 3-vector).
    The event's true channel is the one with its neutron, proton and α counts; events matching no channel are
    skipped. γ rays in the list add to the event's γ count. Channels' own ``gammas`` are not used here.
    """
    rng = rng if rng is not None else np.random.default_rng(seed)
    channels = list(channels)
    names = [c.name for c in channels]
    columns = names + [UNASSIGNED]
    flat = _flatten(events, channels)
    n_events, truth, gammas = len(flat.truth), flat.truth, flat.gammas
    if n_events == 0:
        return SelectionMatrix(names, columns, np.zeros((len(names), len(columns))), np.zeros(len(names)))

    n_det = np.zeros(n_events, dtype=int)
    p_det = np.zeros(n_events, dtype=int)
    a_det = np.zeros(n_events, dtype=int)
    for kind in ("n", "p", "a"):
        sel = flat.kinds == kind
        if not np.any(sel):
            continue
        probs = detection._outcomes(kind, flat.energies[sel], flat.directions[sel])
        u = rng.random(int(sel.sum()))
        cum = np.cumsum(probs, axis=1)
        outcome = (u[:, None] >= cum[:, :2]).sum(axis=1)  # 0, 1 or 2
        rows = flat.event[sel]
        if kind == "n":
            np.add.at(n_det, rows, outcome)
        else:
            np.add.at(p_det, rows, outcome == 1)
            np.add.at(a_det, rows, outcome == 2)
    if detection.gamma_fake() > 0 and gammas.any():
        n_det += rng.binomial(gammas, detection.gamma_fake())

    col_of = {c.multiplicities: columns.index(c.name) for c in channels}
    assigned = np.array([col_of.get((int(n), int(p), int(a)), len(columns) - 1)
                         for n, p, a in zip(n_det, p_det, a_det)])
    counts = np.zeros((len(names), len(columns)))
    np.add.at(counts, (truth, assigned), 1)
    per_row = counts.sum(axis=1)
    m = counts / np.maximum(per_row, 1)[:, None]
    return SelectionMatrix(names, columns, m, per_row)


def isotropic_events(channels: Sequence[Channel], weights: Optional[Sequence[float]] = None, events: int = 10_000,
                     temperature_mev: float = 1.5, seed: int = 1) -> list:
    """Simple stand-in event lists for tests and examples: channels drawn by weight, particles isotropic with a
    Maxwellian energy spectrum E exp(−E/T) (charged particles shifted up by a typical 4 MeV Coulomb barrier for
    protons and 10 MeV for α). Not an evaporation calculation: the true lists come from the evaporation module.
    """
    rng = np.random.default_rng(seed)
    channels = list(channels)
    w = np.ones(len(channels)) if weights is None else np.asarray(weights, dtype=float)
    picks = rng.choice(len(channels), size=events, p=w / w.sum())
    shift = {"n": 0.0, "p": 4.0, "a": 10.0}
    out = []
    for i in picks:
        ch = channels[i]
        parts = []
        for kind, k in (("n", ch.neutrons), ("p", ch.protons), ("alpha", ch.alphas)):
            for _ in range(k):
                e = rng.gamma(2.0, temperature_mev) + shift[KINDS[kind]]
                c, ph = 2 * rng.random() - 1, 2 * math.pi * rng.random()
                s = math.sqrt(1 - c * c)
                parts.append((kind, e, (s * math.cos(ph), s * math.sin(ph), c)))
        out.append({"particles": parts, "gammas": ch.gammas})
    return out
