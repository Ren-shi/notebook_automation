"""Stopping powers, ranges, energy loss and straggling of ions in matter.

Protons and α particles use the NIST PSTAR and ASTAR tables (ICRU Reports 49) directly; heavier ions scale the
proton stopping at the same velocity by an effective charge. Elements NIST does not tabulate are interpolated
between neighbours, and compounds not in the tables use Bragg additivity. Details and validation are in the theory
note ``theory/stopping`` and the physics register.

::

    from physim.nuclear.stopping import Stopping

    s = Stopping("4He", "Au")
    s.stopping_power(5.5)            # MeV/(mg/cm²), electronic + nuclear
    s.range(5.5)                     # mg/cm²
    s.energy_after(5.5, "1 um")      # MeV left after 1 µm of gold
    s.straggling(5.5, "1 um")        # σ of the energy left, MeV

Energies are total kinetic energies in MeV; thicknesses are areal densities in mg/cm² or a ``Quantity``/text such
as ``"1 um"`` (converted with the material's density).
"""

from __future__ import annotations

import csv
import math
from functools import lru_cache
from typing import Optional, Union

import numpy as np

from . import data
from .quantity import Quantity

ArrayLike = Union[float, np.ndarray]

#: Proton and α masses in u, for velocity scaling (AME2020 nuclear masses).
_U_MEV = data.U_KEV * 1e-3
#: Bohr velocity over c (the fine-structure constant).
ALPHA = 7.2973525693e-3
#: e²/(4πε₀), MeV cm, and the Bohr radius, cm (CODATA 2018).
E2_MEV_CM = 1.43996448e-13
A0_CM = 0.529177210903e-8
#: Electron mass, MeV.
ME_MEV = data.ELECTRON_KEV * 1e-3
#: Bohr straggling constant 4π r_e² (m_e c²)² N_A, MeV² cm²/g per unit (z² Z/A).
BOHR_K = 4 * math.pi * (2.8179403262e-13) ** 2 * ME_MEV**2 * data.AVOGADRO

# Short names in physim → NIST STAR material names (compounds and mixtures).
STAR_ALIASES = {
    "Mylar": "Polyethylene Terephthalate (Mylar)",
    "Kapton": "Kapton Polyimide Film",
    "Polyethylene": "Polyethylene",
    "Polystyrene": "Polystyrene",
    "PMMA": "Polymethyl Methacralate (Lucite, Perspex)",  # NIST's spelling
    "Teflon": "Polytetrafluoroethylene (Teflon",  # as in NIST's list
    "PVC": "Polyvinyl Chloride",
    "CsI": "Cesium Iodide",
    "NaI": "Sodium Iodide",
    "LiF": "Lithium Fluoride",
    "SiO2": "Silicon Dioxide",
    "Water": "Water, Liquid",
    "Air": "Air, Dry (near sea level)",
    "Water, Liquid": "Water, Liquid",
    "Air, Dry (near sea level)": "Air, Dry (near sea level)",
}
@lru_cache(maxsize=None)
def _star(prog: str) -> dict:
    """{material name: (energies MeV, electronic, nuclear, csda range)} with stopping in MeV cm²/g."""
    raw: dict = {}
    with open(data.DATA / f"nist_{prog}.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            raw.setdefault((r["material_no"], r["material"]), []).append(
                (float(r["energy_mev"]), float(r["electronic_mev_cm2_g"]), float(r["nuclear_mev_cm2_g"]),
                 float(r["csda_range_g_cm2"])))
    out = {}
    for key, rows in raw.items():
        arr = np.array(rows)
        out[key] = (arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3])
    return out


@lru_cache(maxsize=None)
def _elements_in_tables() -> dict:
    """{Z: material key} for NIST's elemental entries, named "Z: Element". Carbon is graphite (1.7 g/cm³, the
    density physim uses for carbon) rather than amorphous carbon."""
    out = {}
    for no, name in _star("pstar"):
        head = name.split(":", 1)[0].strip()
        if head.isdigit() and not (head == "6" and "Graphite" not in name):
            out[int(head)] = (no, name)
    return out


def _table(prog: str, key: tuple) -> tuple:
    return _star(prog)[key]


def _star_compound(name: str) -> Optional[tuple]:
    full = STAR_ALIASES.get(name, name)
    for key in _star("pstar"):
        if key[1] == full:
            return key
    return None


def _loglog(x: ArrayLike, xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """Log-log interpolation; below the table S ∝ √E (velocity-proportional stopping), above it constant slope."""
    x = np.asarray(x, dtype=float)
    lx, lxs, lys = np.log(np.maximum(x, 1e-300)), np.log(xs), np.log(ys)
    y = np.interp(lx, lxs, lys)
    lo = lx < lxs[0]
    y = np.where(lo, lys[0] + 0.5 * (lx - lxs[0]), y)
    hi = lx > lxs[-1]
    slope = (lys[-1] - lys[-2]) / (lxs[-1] - lxs[-2])
    y = np.where(hi, lys[-1] + slope * (lx - lxs[-1]), y)
    return np.exp(y)


# ---------------------------------------------------------------------------------------------------------------
# Per-element electronic stopping of protons and α particles at a given energy per nucleon


def _per_electron(prog: str, z: int, e_table: ArrayLike, exclude: frozenset = frozenset()) -> np.ndarray:
    """Electronic stopping per electron per atom (MeV cm²) for a tabulated or interpolated element.

    ``exclude`` drops tabulated elements, to test the interpolation on elements NIST does tabulate.
    """
    tab = {z2: key for z2, key in _elements_in_tables().items() if z2 not in exclude}
    if z in tab:
        e, se, _, _ = _table(prog, tab[z])
        a = data.element(z).molar_mass
        return _loglog(e_table, e, se) * a / data.AVOGADRO / z
    # Between the nearest tabulated neighbours, linearly in ln I: exact for Bethe's formula at high energy.
    zs = sorted(tab)
    below = max(x for x in zs if x < z)
    above = min(x for x in zs if x > z)
    i, i0, i1 = (math.log(data.element(x).I_eV) for x in (z, below, above))
    w = (i - i0) / (i1 - i0) if i1 != i0 else (z - below) / (above - below)
    w = min(max(w, -0.5), 1.5)
    f0, f1 = _per_electron(prog, below, e_table, exclude), _per_electron(prog, above, e_table, exclude)
    return np.exp((1 - w) * np.log(f0) + w * np.log(f1))


def _nuclear_zbl(z1: int, m1: float, z2: int, m2: float, e_mev: ArrayLike) -> np.ndarray:
    """ZBL universal nuclear stopping per atom, MeV cm² (Ziegler, Biersack and Littmark 1985)."""
    e_kev = np.asarray(e_mev, dtype=float) * 1e3
    screen = z1**0.23 + z2**0.23
    eps = 32.53 * m2 * e_kev / (z1 * z2 * (m1 + m2) * screen)
    with np.errstate(divide="ignore", invalid="ignore"):
        sn = np.where(eps <= 30.0,
                      np.log1p(1.1383 * eps) / (2 * (eps + 0.01321 * eps**0.21226 + 0.19593 * np.sqrt(eps))),
                      np.log(eps) / (2 * eps))
    # eV / (1e15 atoms/cm²) → MeV cm²
    return 8.462 * z1 * z2 * m1 * sn / ((m1 + m2) * screen) * 1e-15 * 1e-6


def effective_charge(z1: int, e_per_u: ArrayLike, v_fermi: float = 1.0) -> np.ndarray:
    """Fractional effective charge ζ of a heavy ion (Brandt–Kitagawa, with Ziegler's fit for the ionised fraction).

    q is the ionised fraction from the ion's speed relative to the target electrons (Ziegler, Biersack and
    Littmark 1985); the screening term adds the partial charge seen by close collisions. ``v_fermi`` is the target's
    Fermi velocity in units of v₀; 1 is used for all targets (the per-element values are not tabulated here).
    """
    vf = v_fermi
    v = _beta(np.asarray(e_per_u, dtype=float)) / ALPHA
    with np.errstate(divide="ignore", invalid="ignore"):
        vr = np.where(v >= vf, v * (1 + vf**2 / (5 * v**2)), 0.75 * vf * (1 + 2 * v**2 / (3 * vf**2) - v**4 / (15 * vf**4)))
    yr = np.maximum(0.13, vr / z1 ** (2.0 / 3.0))
    a = -0.803 * yr**0.3 + 1.3167 * yr**0.6 + 0.38157 * yr + 0.008983 * yr**2
    q = np.clip(1.0 - np.exp(-np.minimum(a, 50.0)), 0.0, 1.0)
    lam = 2 * 0.24005 * (1 - q) ** (2.0 / 3.0) / (z1 ** (1.0 / 3.0) * (1 - (1 - q) / 7))
    return q + (1 - q) / (2 * vf**2) * np.log1p((4 * lam * vf / 1.919) ** 2)


def lindhard_scharff(z1: int, z2: int, e_per_u: ArrayLike) -> np.ndarray:
    """Lindhard–Scharff electronic stopping per atom at low velocity, MeV cm²:
    8π e² a₀ Z₁^(7/6) Z₂ / (Z₁^(2/3) + Z₂^(2/3))^(3/2) × v/v₀."""
    v = _beta(np.asarray(e_per_u, dtype=float)) / ALPHA
    return 8 * math.pi * E2_MEV_CM * A0_CM * z1 ** (7 / 6) * z2 / (z1 ** (2 / 3) + z2 ** (2 / 3)) ** 1.5 * v


def _beta(e_per_u: np.ndarray) -> np.ndarray:
    g = 1.0 + e_per_u / _U_MEV
    return np.sqrt(1.0 - 1.0 / g**2)


# ---------------------------------------------------------------------------------------------------------------


class Stopping:
    """Energy loss of one ion in one material, tabulated once on a fine grid for fast, vectorised use.

    ``ion`` is a nuclide name (``"1H"``, ``"4He"``, ``"16O"``); ``material`` anything
    :func:`physim.nuclear.data.material` accepts, or a :class:`~physim.nuclear.data.Material`.
    """

    #: Energies per nucleon of the internal grid, MeV/u.
    GRID = np.geomspace(1e-4, 1e4, 801)

    def __init__(self, ion: str, material, density: Union[Quantity, str, float, None] = None):
        self.ion = data.nuclide(ion)
        self.material = material if isinstance(material, data.Material) else data.material(material, density)
        if density is not None and isinstance(material, data.Material):
            self.material = self.material.with_density(density)
        self.mass_u = self.ion.nuclear_mass_mev / _U_MEV
        e_u = self.GRID
        self._se, self._sn = self._compute(e_u)  # MeV cm²/g
        # Range by integrating dE/S in ln E; below the grid S ∝ √E, so R(E₀) = 2E₀/S(E₀).
        e = e_u * self.mass_u
        s = (self._se + self._sn) * 1e-3  # MeV/(mg/cm²)
        integrand = e / s
        ln_e = np.log(e)
        r = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(ln_e))])
        self._range = r + 2.0 * e[0] / s[0]
        self._e = e
        anchor = self._nist_range_anchor()
        if anchor is not None:
            # Protons or α particles in a NIST material: match NIST's own CSDA range at its lowest energy (1 keV),
            # so ranges agree with the tables below ~100 keV too, where the extrapolation below 1 keV matters.
            e0, r0 = anchor
            ours = math.exp(np.interp(math.log(e0), np.log(e), np.log(self._range)))
            below = e < e0
            self._range = np.where(below, self._range * r0 / ours, self._range + (r0 - ours))

    # -- model ------------------------------------------------------------------------------------------------

    def _compute(self, e_u: np.ndarray) -> tuple:
        z1, m1 = self.ion.Z, self.mass_u
        compound = self._nist_key()
        if self.ion.Z <= 2 and compound is not None:
            prog = "pstar" if z1 == 1 else "astar"
            ref_mass = 1.007276467 if z1 == 1 else 4.001506179
            e, se, sn, _ = _table(prog, compound)
            e_ref = e_u * ref_mass  # the tabulated particle at the same velocity
            se_out = _loglog(e_ref, e, se)
            sn_out = _loglog(e_ref, e, sn) if self.ion.A == round(ref_mass) else self._nuclear(e_u)
            return se_out, sn_out
        per_atom = np.zeros_like(e_u)
        for (z2, _), n in self._atoms_per_unit().items():
            per_atom += n * z2 * self._electronic_per_electron(z2, e_u)
        se = per_atom * data.AVOGADRO / self.material.molar_mass
        return se, self._nuclear(e_u)

    def _nist_key(self) -> Optional[tuple]:
        """The NIST table of this material: a tabulated compound, or a natural tabulated element."""
        key = _star_compound(self.material.name)
        if key is None and self.material.name not in STAR_ALIASES:
            try:
                parts = data.parse_material(self.material.name)
            except ValueError:
                return None
            if len(parts) == 1 and parts[0][1] is None and parts[0][0] in _elements_in_tables():
                key = _elements_in_tables()[parts[0][0]]
        return key

    def _nist_range_anchor(self) -> Optional[tuple]:
        """(lowest tabulated energy MeV, CSDA range mg/cm²) when this is a proton or α in a NIST material."""
        if (self.ion.Z, self.ion.A) not in ((1, 1), (2, 4)):
            return None
        key = self._nist_key()
        if key is None:
            return None
        e, _, _, r = _table("pstar" if self.ion.Z == 1 else "astar", key)
        return float(e[0]), float(r[0]) * 1e3

    def _electronic_per_electron(self, z2: int, e_u: np.ndarray) -> np.ndarray:
        z1 = self.ion.Z
        if z1 == 1:
            return _per_electron("pstar", z2, e_u * 1.007276467)
        if z1 == 2:
            return _per_electron("astar", z2, e_u * 4.001506179)
        # Effective-charge scaling of the proton at the same velocity; at low velocity, where that fails, the larger
        # of it and Lindhard–Scharff joined to full-charge stopping (1/S = 1/S_LS + 1/(Z₁² S_p)).
        sp = _per_electron("pstar", z2, e_u * 1.007276467)
        scaled = (z1 * effective_charge(z1, e_u)) ** 2 * sp
        slow = 1.0 / (z2 / lindhard_scharff(z1, z2, e_u) + 1.0 / (z1**2 * sp))
        return np.maximum(scaled, slow)

    def _atoms_per_unit(self) -> dict:
        out: dict = {}
        for z, a, n in self.material.atoms:
            out[(z, a)] = out.get((z, a), 0.0) + n
        return out

    def _nuclear(self, e_u: np.ndarray) -> np.ndarray:
        z1, m1 = self.ion.Z, self.mass_u
        total = np.zeros_like(e_u)
        for (z2, a2), n in self._atoms_per_unit().items():
            m2 = data.nuclide((z2, a2)).atomic_mass_u
            total += n * _nuclear_zbl(z1, m1, z2, m2, e_u * m1)
        return total * data.AVOGADRO / self.material.molar_mass

    # -- results ----------------------------------------------------------------------------------------------

    def _areal(self, thickness) -> float:
        if isinstance(thickness, (int, float, np.floating)):
            return float(thickness)
        return self.material.areal_density_mg_cm2(thickness)

    def stopping_power(self, energy: ArrayLike, part: str = "total") -> ArrayLike:
        """Stopping power at kinetic energy ``energy`` (MeV), MeV/(mg/cm²): ``"total"``, ``"electronic"`` or
        ``"nuclear"``."""
        if part not in ("total", "electronic", "nuclear"):
            raise ValueError("part must be 'total', 'electronic' or 'nuclear'")
        ys = {"electronic": self._se, "nuclear": self._sn, "total": self._se + self._sn}[part] * 1e-3
        out = _loglog(np.asarray(energy, dtype=float), self._e, np.maximum(ys, 1e-300))
        return float(out) if np.ndim(energy) == 0 else out

    def range(self, energy: ArrayLike) -> ArrayLike:
        """Path length to stop from ``energy`` (MeV), mg/cm² (CSDA range)."""
        e = np.asarray(energy, dtype=float)
        out = np.exp(np.interp(np.log(np.maximum(e, 1e-300)), np.log(self._e), np.log(self._range)))
        out = np.where(e <= 0, 0.0, out)
        return float(out) if np.ndim(energy) == 0 else out

    def energy_after(self, energy: ArrayLike, thickness, tilt_deg: float = 0.0) -> ArrayLike:
        """Kinetic energy (MeV) left after crossing ``thickness`` at ``tilt_deg`` from the layer normal; 0 if the
        ion stops. Uses the range: R(E_out) = R(E_in) − path."""
        path = self._areal(thickness) / math.cos(math.radians(tilt_deg))
        left = self.range(energy) - path
        e = np.exp(np.interp(np.log(np.maximum(left, 1e-300)), np.log(self._range), np.log(self._e)))
        out = np.where(left > 0, e, 0.0)
        return float(out) if np.ndim(energy) == 0 else out

    def energy_loss(self, energy: ArrayLike, thickness, tilt_deg: float = 0.0) -> ArrayLike:
        """Energy lost in the layer, MeV."""
        return np.asarray(energy) - self.energy_after(energy, thickness, tilt_deg) if np.ndim(energy) else (
            energy - self.energy_after(energy, thickness, tilt_deg))

    def thickness_to_stop(self, energy: ArrayLike) -> ArrayLike:
        """Thickness that stops the ion (its range), mg/cm²."""
        return self.range(energy)

    def straggling(self, energy: float, thickness, tilt_deg: float = 0.0, steps: int = 200,
                   propagate: bool = True) -> float:
        """Standard deviation (MeV) of the energy after the layer, from Bohr straggling carried through it.

        Integrates dσ²/dx = 2 (dS/dE) σ² + dΩ²/dx along the path, so a thick layer's growing or shrinking
        spread follows the slope of the stopping power. Ω² uses the effective charge, the relativistic factor
        (1 − β²/2)/(1 − β²) and Lindhard and Scharff's reduction for slow ions. The FWHM is 2.355 σ.
        ``propagate=False`` just adds up Ω² along the path (no growth or shrinking), for comparison.
        """
        path = self._areal(thickness) / math.cos(math.radians(tilt_deg))
        z_over_a = sum(n * z for (z, _), n in self._atoms_per_unit().items()) / self.material.molar_mass
        e, var = float(energy), 0.0
        dx = path / steps
        for _ in range(steps):
            if e <= 0:
                return 0.0
            e_mid = max(self.energy_after(e, dx / 2), 1e-12)
            s = self.stopping_power(e_mid)
            ds_de = (self.stopping_power(e_mid * 1.001) - self.stopping_power(e_mid * 0.999)) / (0.002 * e_mid)
            z_eff2 = self._z_eff2(e_mid / self.mass_u)
            b2 = float(_beta(np.array(e_mid / self.mass_u))) ** 2
            omega2 = (BOHR_K * z_eff2 * z_over_a * dx * 1e-3 * (1 - b2 / 2) / (1 - b2)
                      * self._lindhard_scharff_factor(b2))
            # Two ions δ apart in energy drift as dδ/dx = −S′(E) δ, so the spread is carried by exp(−2 S′ dx):
            # it grows above the Bragg peak (S′ < 0) and shrinks below it.
            var = var * math.exp(-2 * ds_de * dx) * (1.0 if propagate else 0.0) + omega2 if propagate else var + omega2
            e = self.energy_after(e, dx)
        return math.sqrt(var)

    def _lindhard_scharff_factor(self, beta2: float) -> float:
        """Reduction of Bohr straggling for slow ions: L(χ)/2 with L = 1.36 χ^½ − 0.016 χ^(3/2) for χ < 3, where
        χ = v²/(Z₂ v₀²) (Lindhard and Scharff 1953), averaged over the material's elements by electrons."""
        total, weight = 0.0, 0.0
        for (z2, _), n in self._atoms_per_unit().items():
            chi = beta2 / ALPHA**2 / z2
            f = 0.5 * (1.36 * math.sqrt(chi) - 0.016 * chi**1.5) if chi < 3.0 else 1.0
            total += n * z2 * min(f, 1.0)
            weight += n * z2
        return total / weight

    def _z_eff2(self, e_u: float) -> float:
        z1 = self.ion.Z
        if z1 == 1:
            return 1.0
        p = Stopping._proton_cache(self.material)
        return float(self.stopping_power(e_u * self.mass_u, "electronic") / p.stopping_power(
            e_u * p.mass_u, "electronic"))

    @staticmethod
    @lru_cache(maxsize=64)
    def _proton_cache(material) -> Stopping:
        return Stopping("1H", material)

    def __repr__(self) -> str:
        return f"Stopping({self.ion} in {self.material.name})"


def radiation_length_g_cm2(material: data.Material) -> float:
    """Radiation length X₀, g/cm² (Dahl's fit, 716.4 A / (Z(Z+1) ln(287/√Z)), combined by mass fractions)."""
    inv = 0.0
    m_total = material.molar_mass
    for z, a, n in material.atoms:
        m = n * data.nuclide((z, a)).atomic_mass_u if a is not None else n * data.element(z).molar_mass
        a_mass = data.nuclide((z, a)).atomic_mass_u
        x0 = 716.4 * a_mass / (z * (z + 1) * math.log(287 / math.sqrt(z)))
        inv += (m / m_total) / x0
    return 1.0 / inv


def angular_straggling_mrad(ion: str, energy: float, material, thickness, density=None) -> float:
    """Width (σ of the projected angle, mrad) of multiple scattering in a layer: Highland's formula.

    Accurate (≈11%) for fast particles and layers between 1e-3 and 100 radiation lengths; for slow heavy ions it
    is only indicative.
    """
    mat = material if isinstance(material, data.Material) else data.material(material, density)
    n = data.nuclide(ion)
    t = mat.areal_density_mg_cm2(thickness) if not isinstance(thickness, (int, float)) else float(thickness)
    x = t * 1e-3 / radiation_length_g_cm2(mat)
    m = n.nuclear_mass_mev
    p = math.sqrt(energy * (energy + 2 * m))
    beta = p / (energy + m)
    return 1e3 * 13.6 / (beta * p) * n.Z * math.sqrt(x) * (1 + 0.038 * math.log(x * n.Z**2 / beta**2))
