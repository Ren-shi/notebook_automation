"""Nuclear and material data: atomic masses (AME2020), natural isotopic compositions, element and compound
densities and mean excitation energies (NIST). Sources and licences are in ``data/SOURCES.md``.

::

    from physim.nuclear import data

    data.nuclide("208Pb").mass_excess_kev       # -21748.519 keV
    data.q_value(["2H", "3H"], ["4He", "n"])     # 17.589 MeV
    gold = data.material("Au")
    gold.areal_density_mg_cm2("1 um")            # 1.932 mg/cm²
    gold.atoms_per_cm2("1 mg/cm2")               # {(79, 197): 3.06e18}
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Optional, Union

from .names import NAMED_MATERIALS, SYMBOLS, Z_OF, parse_material, parse_nuclide
from .quantity import Quantity

DATA = Path(__file__).resolve().parent / "data"

#: Atomic mass unit, keV (CODATA 2018, as used by AME2020).
U_KEV = 931494.10242
#: Electron mass, keV (CODATA 2018).
ELECTRON_KEV = 510.99895000
#: Avogadro constant, 1/mol (exact).
AVOGADRO = 6.02214076e23

NuclideLike = Union[str, tuple]


def electron_binding_kev(z: int) -> float:
    """Total binding energy of the Z electrons of a neutral atom, keV.

    The fit B = 14.4381 Z^2.39 + 1.55468e-6 Z^5.35 eV of D. Lunney, J. M. Pearson and C. Thibault,
    Rev. Mod. Phys. 75, 1021 (2003), eq. (A4), good to a few tens of eV across the table; the same correction AME
    uses to relate atomic and nuclear masses.
    """
    return (14.4381 * z**2.39 + 1.55468e-6 * z**5.35) * 1e-3


# ---------------------------------------------------------------------------------------------------------------
# Nuclides


@dataclass(frozen=True)
class Nuclide:
    """One nuclide from AME2020. Masses are of the neutral atom unless the name says nuclear."""

    Z: int
    A: int
    mass_excess_kev: float
    mass_excess_unc_kev: float
    #: True when AME marks the mass as estimated from systematics rather than measured (``#`` in the table).
    estimated: bool

    @property
    def N(self) -> int:
        return self.A - self.Z

    @property
    def symbol(self) -> str:
        return SYMBOLS[self.Z]

    @property
    def name(self) -> str:
        return "n" if self.Z == 0 else f"{self.A}{self.symbol}"

    @property
    def atomic_mass_u(self) -> float:
        """Atomic mass, u (for the neutron, its mass)."""
        return self.A + self.mass_excess_kev / U_KEV

    @property
    def atomic_mass_mev(self) -> float:
        return (self.A * U_KEV + self.mass_excess_kev) * 1e-3

    @property
    def nuclear_mass_mev(self) -> float:
        """Mass of the bare nucleus: the atomic mass less Z electron masses, plus the electrons' binding energy."""
        return (self.A * U_KEV + self.mass_excess_kev - self.Z * ELECTRON_KEV + electron_binding_kev(self.Z)) * 1e-3

    def __str__(self) -> str:
        return self.name


def _ame_number(text: str) -> tuple[float, bool]:
    text = text.strip()
    return float(text.replace("#", ".")), "#" in text


@lru_cache(maxsize=None)
def _ame() -> dict[tuple[int, int], Nuclide]:
    """Parse ``mass_1.mas20.txt`` (Fortran format a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,...)."""
    table: dict[tuple[int, int], Nuclide] = {}
    with open(DATA / "mass_1.mas20.txt", encoding="ascii") as f:
        for line in f:
            if len(line) < 54:
                continue
            try:
                z, a = int(line[9:14]), int(line[14:19])
                int(line[4:9])  # N, to confirm this is a data line
            except ValueError:
                continue
            excess, est = _ame_number(line[28:42])
            unc, _ = _ame_number(line[42:54])
            table[(z, a)] = Nuclide(z, a, excess, unc, est)
    return table


def nuclide(name: NuclideLike) -> Nuclide:
    """Look up a nuclide by name (``"208Pb"``, ``"alpha"``, ``"n"``) or as ``(Z, A)``."""
    z, a = name if isinstance(name, tuple) else parse_nuclide(name)
    try:
        return _ame()[(z, a)]
    except KeyError:
        raise ValueError(f"{a}{SYMBOLS[z]} is not in the AME2020 mass table") from None


def nuclides() -> list[Nuclide]:
    """Every nuclide in the AME2020 table."""
    return list(_ame().values())


def q_value(entrance: Iterable[NuclideLike], exit: Iterable[NuclideLike], excitation_mev: float = 0.0) -> float:
    """Q-value in MeV of a reaction ``entrance → exit`` (e.g. ``["2H", "3H"] → ["4He", "n"]``).

    Computed from atomic mass excesses, where electron masses cancel because charge is conserved (electron
    binding energies, at most tens of keV in total, cancel to within a few eV for light ions). ``excitation_mev``
    is the excitation energy left in the exit channel, which lowers Q.
    """
    ins, outs = [nuclide(x) for x in entrance], [nuclide(x) for x in exit]
    for what, f in (("charge", lambda n: n.Z), ("mass number", lambda n: n.A)):
        if sum(map(f, ins)) != sum(map(f, outs)):
            raise ValueError(f"the reaction does not conserve {what}: "
                             f"{' + '.join(map(str, ins))} → {' + '.join(map(str, outs))}")
    return (sum(n.mass_excess_kev for n in ins) - sum(n.mass_excess_kev for n in outs)) * 1e-3 - excitation_mev


# ---------------------------------------------------------------------------------------------------------------
# Elements


def _value(text: str) -> Optional[float]:
    """A NIST number such as ``196.966569(5)`` or ``1.932E+01``, without its uncertainty."""
    m = re.match(r"^\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)", text or "")
    return float(m.group(1)) if m else None


@dataclass(frozen=True)
class Element:
    """An element: natural isotopic composition (NIST), density and mean excitation energy (NIST, Z ≤ 92)."""

    Z: int
    symbol: str
    name: Optional[str]
    #: As published by NIST: a value with uncertainty, an interval ``[a,b]``, or ``[A]`` for no stable isotope.
    standard_atomic_weight: Optional[str]
    #: Natural isotopes as (A, atom fraction).
    isotopes: tuple
    density_g_cm3: Optional[float]
    #: Mean excitation energy, eV (used for stopping powers).
    I_eV: Optional[float]

    @property
    def molar_mass(self) -> float:
        """Mean atomic mass of the natural element, g/mol: abundances (NIST) × AME2020 atomic masses."""
        if not self.isotopes:
            raise ValueError(f"{self.symbol} has no natural isotopic composition; name an isotope, e.g. "
                             f"'{_longest_known(self.Z)}{self.symbol}'")
        total = sum(f for _, f in self.isotopes)
        return sum(f * nuclide((self.Z, a)).atomic_mass_u for a, f in self.isotopes) / total


def _longest_known(z: int) -> int:
    weight = _elements_raw().get(z, {}).get("standard_atomic_weight", "")
    m = re.match(r"^\[(\d+)\]$", weight or "")
    return int(m.group(1)) if m else 2 * z


@lru_cache(maxsize=None)
def _elements_raw() -> dict[int, dict]:
    with open(DATA / "nist_elements.csv", encoding="utf-8") as f:
        return {int(r["Z"]): r for r in csv.DictReader(f)}


@lru_cache(maxsize=None)
def _isotopes_raw() -> dict[int, tuple]:
    out: dict[int, list] = {}
    with open(DATA / "nist_isotopes.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out.setdefault(int(r["Z"]), []).append((int(r["A"]), _value(r["isotopic_composition"])))
    return {z: tuple(v) for z, v in out.items()}


@lru_cache(maxsize=None)
def element(symbol: Union[str, int]) -> Element:
    """Look up an element by symbol (``"Au"``) or atomic number."""
    if isinstance(symbol, str):
        if symbol not in Z_OF:
            raise ValueError(f"'{symbol}' is not an element symbol")
        z = Z_OF[symbol]
    else:
        z = int(symbol)
        if not 1 <= z < len(SYMBOLS):
            raise ValueError(f"no element has Z = {z}")
    row = _elements_raw().get(z, {})
    return Element(Z=z, symbol=SYMBOLS[z], name=row.get("name"),
                   standard_atomic_weight=row.get("standard_atomic_weight"),
                   isotopes=_isotopes_raw().get(z, ()), density_g_cm3=_value(row.get("density_g_cm3")),
                   I_eV=_value(row.get("I_eV")))


# ---------------------------------------------------------------------------------------------------------------
# Materials


@lru_cache(maxsize=None)
def _compounds_raw() -> dict[str, dict]:
    with open(DATA / "nist_compounds.csv", encoding="utf-8") as f:
        return {r["name"]: r for r in csv.DictReader(f)}


#: Short names for NIST compounds (the full NIST names are accepted too, see :func:`compound_names`).
COMPOUND_ALIASES = {
    "Mylar": "Polyethylene Terephthalate, (Mylar)",
    "Polyethylene": "Polyethylene",
    "Polystyrene": "Polystyrene",
    "PMMA": "Polymethyl Methacrylate",
    "Teflon": "Polytetrafluoroethylene, (Teflon)",
    "PVC": "Polyvinyl Chloride",
    "CsI": "Cesium Iodide",
    "LiF": "Lithium Fluride",  # NIST's spelling
    "GaAs": "Gallium Arsenide",
    "CdTe": "Cadmium Telluride",
    "Water": "Water, Liquid",
    "Air": "Air, Dry (near sea level)",
}


def compound_names() -> list[str]:
    """Names of the NIST compounds and mixtures, and their short aliases."""
    return list(COMPOUND_ALIASES) + list(_compounds_raw())


@dataclass(frozen=True)
class Material:
    """A material as nuclides and how many of each make one formula unit.

    ``atoms`` lists (Z, A, atoms per formula unit) with natural elements expanded into their isotopes, so
    ``molar_mass`` is the mass of one formula unit in g/mol. ``density_g_cm3`` is ``None`` when it is not known;
    thicknesses can then only be given as areal densities.
    """

    name: str
    atoms: tuple
    density_g_cm3: Optional[float] = None
    #: Mean excitation energy, eV, when NIST gives one for this material (used for stopping powers).
    I_eV: Optional[float] = None

    @property
    def molar_mass(self) -> float:
        """g/mol per formula unit."""
        return sum(n * nuclide((z, a)).atomic_mass_u for z, a, n in self.atoms)

    @property
    def elements(self) -> dict[int, float]:
        """Atoms of each element (by Z) per formula unit."""
        out: dict[int, float] = {}
        for z, _, n in self.atoms:
            out[z] = out.get(z, 0.0) + n
        return out

    def with_density(self, density: Union[Quantity, str, float, None]) -> Material:
        """The same material with a density (``"1.06 g/cm3"`` or a number in g/cm³)."""
        if density is None:
            return self
        rho = Quantity.parse(density).to("g/cm3") if isinstance(density, str) else (
            density.to("g/cm3") if isinstance(density, Quantity) else float(density))
        return Material(self.name, self.atoms, rho, self.I_eV)

    def areal_density_mg_cm2(self, thickness: Union[Quantity, str]) -> float:
        """Thickness as an areal density, mg/cm²; a length needs the density."""
        q = Quantity.parse(thickness) if isinstance(thickness, str) else thickness
        if q.kind == "areal_density":
            return q.to("mg/cm2")
        if q.kind != "length":
            raise ValueError(f"'{q}' is not a thickness")
        if self.density_g_cm3 is None:
            raise ValueError(f"the density of {self.name} is not known: give the thickness in mg/cm2, "
                             "or give a density")
        return self.density_g_cm3 * q.to("cm") * 1e3

    def thickness_um(self, thickness: Union[Quantity, str]) -> float:
        """Thickness as a length, µm (needs the density)."""
        if self.density_g_cm3 is None:
            raise ValueError(f"the density of {self.name} is not known, so its thickness in um cannot be found")
        return self.areal_density_mg_cm2(thickness) * 1e-3 / self.density_g_cm3 * 1e4

    def atoms_per_cm2(self, thickness: Union[Quantity, str]) -> dict[tuple, float]:
        """Atoms per cm² of each nuclide (Z, A) in a layer of this thickness."""
        units = self.areal_density_mg_cm2(thickness) * 1e-3 * AVOGADRO / self.molar_mass
        out: dict[tuple, float] = {}
        for z, a, n in self.atoms:
            out[(z, a)] = out.get((z, a), 0.0) + n * units
        return out

    @classmethod
    def enriched(cls, symbol: str, fractions: dict[int, float], name: Optional[str] = None) -> Material:
        """An isotopically enriched element, e.g. ``Material.enriched("C", {13: 0.99, 12: 0.01})`` (atom fractions).

        The density is the natural element's, scaled by the ratio of mean atomic masses (same number of atoms per
        volume).
        """
        el = element(symbol)
        total = sum(fractions.values())
        if total <= 0 or any(f < 0 for f in fractions.values()):
            raise ValueError("isotope fractions must be non-negative and not all zero")
        atoms = tuple((el.Z, a, f / total) for a, f in sorted(fractions.items()))
        for _, a, _ in atoms:
            nuclide((el.Z, a))
        m = cls(name or f"enriched {symbol}", atoms, None, el.I_eV)
        if el.density_g_cm3 is not None and el.isotopes:
            m = m.with_density(el.density_g_cm3 * m.molar_mass / el.molar_mass)
        return m


def _natural(z: int, count: float) -> list[tuple]:
    el = element(z)
    if not el.isotopes:
        raise ValueError(f"{el.symbol} has no natural isotopic composition; name an isotope, e.g. "
                         f"'{_longest_known(z)}{el.symbol}'")
    total = sum(f for _, f in el.isotopes)
    return [(z, a, count * f / total) for a, f in el.isotopes]


@lru_cache(maxsize=None)
def _material(name: str) -> Material:
    full = COMPOUND_ALIASES.get(name, name)
    if full in _compounds_raw():
        r = _compounds_raw()[full]
        # Mass fractions → atoms per formula unit, with the formula unit chosen to weigh 1 g/mol.
        atoms: list[tuple] = []
        for part in r["mass_fractions"].split():
            z_text, w_text = part.split(":")
            z, w = int(z_text), float(w_text)
            atoms += _natural(z, w / element(z).molar_mass)
        return Material(name, tuple(atoms), _value(r["density_g_cm3"]), _value(r["I_eV"]))
    parts = parse_material(name)
    atoms = []
    for z, a, n in parts:
        atoms += [(z, a, n)] if a is not None else _natural(z, n)
    density, i_ev = None, None
    if len(parts) == 1:
        z, a, _ = parts[0]
        el = element(z)
        i_ev = el.I_eV
        if el.density_g_cm3 is not None:
            if a is None:
                density = el.density_g_cm3
            elif el.isotopes:  # a pure isotope: same atoms per volume as the natural element
                density = el.density_g_cm3 * nuclide((z, a)).atomic_mass_u / el.molar_mass
    return Material(name, tuple(atoms), density, i_ev)


def material(name: str, density: Union[Quantity, str, float, None] = None) -> Material:
    """A material by name: an element (``"Au"``, natural abundances), an isotope (``"208Pb"``), a formula
    (``"CD2"``) or a NIST compound (``"Mylar"``, ``"CsI"``, ``"Water, Liquid"``; see :func:`compound_names`).

    Elements and NIST compounds come with NIST densities; formulas need ``density`` if thicknesses are given as
    lengths. A ``density`` given here overrides the tabulated one.
    """
    if name in NAMED_MATERIALS and name not in COMPOUND_ALIASES:
        m = _material(NAMED_MATERIALS[name])
        m = Material(name, m.atoms, m.density_g_cm3, m.I_eV)
    else:
        m = _material(name)
    return m.with_density(density)
