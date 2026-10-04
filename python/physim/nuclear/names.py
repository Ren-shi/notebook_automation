"""Reading nuclide and material names: ``"4He"``, ``"alpha"``, ``"197Au"``, ``"CD2"``, ``"Mylar"``.

Only the names are checked here. Masses, densities and isotopic compositions come from the nuclear and material
data (backlog item 34).
"""

from __future__ import annotations

import re

#: Element symbols indexed by atomic number (index 0 is the neutron, written ``n``).
SYMBOLS = (
    "n H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr "
    "Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv "
    "Ts Og"
).split()
Z_OF = {s: z for z, s in enumerate(SYMBOLS) if z > 0}

#: Common names for light nuclides, as (Z, A).
NUCLIDE_ALIASES = {
    "p": (1, 1), "proton": (1, 1), "d": (1, 2), "deuteron": (1, 2), "t": (1, 3), "triton": (1, 3),
    "alpha": (2, 4), "a": (2, 4), "D": (1, 2), "T": (1, 3), "n": (0, 1), "neutron": (0, 1),
}

#: Named materials and the chemical formula each stands for.
NAMED_MATERIALS = {
    "Mylar": "C10H8O4",
    "Kapton": "C22H10N2O5",
    "Polyethylene": "CH2",
}

_NUCLIDE = re.compile(r"^(\d+)([A-Z][a-z]?)$|^([A-Z][a-z]?)-?(\d+)$")
_FORMULA_PART = re.compile(r"(\d*)([A-Z][a-z]?)(\d*\.?\d*)")


def parse_nuclide(name: str) -> tuple[int, int]:
    """(Z, A) of a nuclide written ``"4He"``, ``"He-4"``, ``"He4"`` or an alias such as ``"alpha"``."""
    if not isinstance(name, str):
        raise ValueError(f"a nuclide is written as text such as '4He', got {name!r}")
    if name in NUCLIDE_ALIASES:
        return NUCLIDE_ALIASES[name]
    m = _NUCLIDE.match(name.strip())
    if not m:
        raise ValueError(f"cannot read '{name}' as a nuclide; write the mass number and symbol, e.g. '4He'")
    a_text, symbol = (m.group(1), m.group(2)) if m.group(1) else (m.group(4), m.group(3))
    if symbol not in Z_OF:
        raise ValueError(f"'{symbol}' in '{name}' is not an element symbol")
    z, a = Z_OF[symbol], int(a_text)
    if a < z or a > 3 * z + 20:
        raise ValueError(f"'{name}': mass number {a} is not possible for {symbol} (Z = {z})")
    return z, a


def parse_material(name: str) -> list[tuple[int, int | None, float]]:
    """Composition of a material as (Z, A or None for natural abundance, atoms per formula unit).

    Accepts an element (``"Au"``), an isotope (``"197Au"``, also ``"D"`` for deuterium), a formula
    (``"CD2"``, ``"C10H8O4"``, ``"13CH4"``) or a named material (``"Mylar"``).
    """
    if not isinstance(name, str) or not name.strip():
        raise ValueError(f"a material is written as text such as 'Au', '12C' or 'CD2', got {name!r}")
    text = name.strip()
    if text in NAMED_MATERIALS:
        text = NAMED_MATERIALS[text]
    parts: list[tuple[int, int | None, float]] = []
    pos = 0
    for m in _FORMULA_PART.finditer(text):
        if m.start() != pos:
            break
        pos = m.end()
        a_text, symbol, count_text = m.groups()
        if symbol in ("D", "T"):  # deuterium and tritium
            if a_text:
                raise ValueError(f"material '{name}': write '{symbol}' or '{'2' if symbol == 'D' else '3'}H', not both")
            z, a = 1, 2 if symbol == "D" else 3
        elif symbol in Z_OF:
            z = Z_OF[symbol]
            a = int(a_text) if a_text else None
            if a is not None and (a < z or a > 3 * z + 20):
                raise ValueError(f"material '{name}': mass number {a} is not possible for {symbol}")
        else:
            raise ValueError(f"'{symbol}' in material '{name}' is not an element symbol")
        count = float(count_text) if count_text else 1.0
        if count <= 0:
            raise ValueError(f"material '{name}': the count of {symbol} must be positive")
        parts.append((z, a, count))
    if pos != len(text) or not parts:
        known = ", ".join(NAMED_MATERIALS)
        raise ValueError(f"cannot read '{name}' as a material; use an element ('Au'), an isotope ('197Au'), "
                         f"a formula ('CD2') or a named material ({known})")
    return parts
