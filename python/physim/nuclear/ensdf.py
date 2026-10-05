"""Reading adopted levels and γ rays from a local copy of ENSDF (the Evaluated Nuclear Structure Data File).

ENSDF is not shipped with physim: permission to redistribute it has not been confirmed (backlog item 50). A copy
is downloaded once with ``scripts/fetch_ensdf.py`` into :func:`folder` and read from there, zipped or unzipped::

    from physim.nuclear import ensdf

    ensdf.available()                  # is there a local copy?
    data = ensdf.adopted("194Pt")      # the "ADOPTED LEVELS, GAMMAS" dataset
    data.levels[1].energy              # (328.464, 0.012) keV
    data.levels[1].gammas[0].final     # index of the level the γ ray feeds

The file format is ENSDF's 80-column card format (J. K. Tuli, "Evaluated Nuclear Structure Data File: a manual for
preparation of data sets", BNL-NCS-51655-01/02-Rev). Only what a level scheme needs is read: level energies, spins,
parities, half-lives and moments, and γ-ray energies, intensities, multipolarities, mixing ratios, conversion
coefficients and reduced transition probabilities in Weisskopf units.
"""

from __future__ import annotations

import os
import re
import zipfile
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterator, Optional

from .names import SYMBOLS, parse_nuclide

#: Environment variable naming the folder that holds the local ENSDF copy.
ENV = "PHYSIM_ENSDF"

#: ħ in eV s (CODATA 2018), to turn a level width into a half-life.
HBAR_EV_S = 6.582119569e-16
LN2 = 0.6931471805599453

#: ENSDF time units in seconds (a year is 365.2422 days, as ENSDF defines it).
TIME_UNITS = {"Y": 365.2422 * 86400.0, "D": 86400.0, "H": 3600.0, "M": 60.0, "S": 1.0, "MS": 1e-3, "US": 1e-6,
              "NS": 1e-9, "PS": 1e-12, "FS": 1e-15, "AS": 1e-18}
WIDTH_UNITS = {"EV": 1.0, "KEV": 1e3, "MEV": 1e6}

Number = Optional[tuple]  # (value, uncertainty or None)


class EnsdfMissing(LookupError):
    """There is no local copy of ENSDF, or it has no adopted levels for the nuclide."""


def folder() -> Path:
    """Where the local copy lives: ``$PHYSIM_ENSDF``, or ``~/.physim/ensdf``."""
    return Path(os.environ[ENV]) if os.environ.get(ENV) else Path.home() / ".physim" / "ensdf"


def available() -> bool:
    """Whether a local copy of ENSDF is present."""
    d = folder()
    return d.is_dir() and any(p.is_file() and p.name != "README.md" for p in d.iterdir())


# ---------------------------------------------------------------------------------------------------------------
# Fields


def number(text: str, unc: str = "") -> Number:
    """An ENSDF number and its uncertainty, which is written in units of the value's last digit: ``"328.464"`` with
    ``"12"`` is 328.464 ± 0.012. Returns ``None`` when the field holds no plain number. Asymmetric uncertainties
    (``"+5-3"``) give the larger side; limits and approximations (``LT``, ``AP``, ...) give no uncertainty."""
    text = text.strip()
    try:
        value = float(text)
    except ValueError:
        return None
    mantissa = re.split(r"[eE]", text)[0]
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    exponent = int(re.split(r"[eE]", text)[1]) if re.search(r"[eE]", text) else 0
    step = 10.0 ** (exponent - decimals)
    unc = unc.strip()
    if unc.isdigit():
        return value, int(unc) * step
    m = re.fullmatch(r"\+(\d+)-(\d+)", unc) or re.fullmatch(r"-(\d+)\+(\d+)", unc)
    if m:
        return value, max(int(m.group(1)), int(m.group(2))) * step
    return value, None


def half_life(text: str, unc: str = "") -> Number:
    """A half-life in seconds from ENSDF's T field (``"41.7 PS"``, ``"2.5 EV"`` for a width). ``None`` for
    ``STABLE``, an empty field or a limit that cannot be read."""
    parts = text.split()
    if len(parts) != 2:
        return None
    n = number(parts[0], unc)
    if n is None or n[0] <= 0:
        return None
    value, err = n
    unit = parts[1].upper()
    if unit in TIME_UNITS:
        return value * TIME_UNITS[unit], None if err is None else err * TIME_UNITS[unit]
    if unit in WIDTH_UNITS:  # a width Γ: T½ = ħ ln2 / Γ, with the same relative uncertainty
        t = HBAR_EV_S * LN2 / (value * WIDTH_UNITS[unit])
        return t, None if err is None else t * err / value
    return None


def spin_parity(text: str) -> tuple:
    """(J, parity, tentative) from ENSDF's J field. J and parity are ``None`` unless the field names one value:
    ``"2+"`` → (2.0, +1, False); ``"(3/2)-"`` → (1.5, −1, True); ``"2+,3+"`` → (None, None, True)."""
    t = text.strip()
    tentative = "(" in t or "[" in t
    m = re.fullmatch(r"(\d+)(/2)?([+-])?", re.sub(r"[()\[\]]", "", t))
    if not m:
        return None, None, tentative or bool(t)
    j = int(m.group(1)) / (2.0 if m.group(2) else 1.0)
    parity = {"+": 1, "-": -1, None: None}[m.group(3)]
    return j, parity, tentative


# ---------------------------------------------------------------------------------------------------------------
# Records


@dataclass
class Gamma:
    energy: tuple
    intensity: Number = None
    multipolarity: str = ""
    mixing: Number = None
    conversion: Number = None
    #: Reduced transition probabilities in Weisskopf units, by multipolarity: {"E2": (49.2, 1.4)}.
    weisskopf: dict = field(default_factory=dict)
    #: Energy of the final level if the file names it (``FL=``), keV.
    final_energy: Optional[float] = None
    #: Index of the level this γ ray feeds, or ``None`` when no level matches.
    final: Optional[int] = None


@dataclass
class Level:
    energy: tuple
    jpi: str = ""
    half_life: Number = None
    stable: bool = False
    #: Static moments from the level's continuation records: {"E2": (−0.48, 0.14)} in barns, {"M1": ...} in μN.
    moments: dict = field(default_factory=dict)
    #: Moments written without a sign, so only their size is known.
    unsigned: set = field(default_factory=set)
    gammas: list = field(default_factory=list)


@dataclass
class Adopted:
    """The adopted levels and γ rays of one nuclide."""

    nuclide: str
    #: The dataset's reference and date fields, for the citation.
    reference: str
    date: str
    levels: list
    #: γ rays left out because no final level matched.
    unplaced: int = 0


_ITEM = re.compile(r"^\s*(B[EM]\d)W\s*=\s*([-+.\dEe]+)\s*(\S*)|^\s*(MOM[EM]\d)\s*=\s*([-+.\dEe]+)\s*(\S*)"
                   r"|^\s*(FL)\s*=\s*([-+.\dEe]+)|^\s*(CC)\s*=\s*([-+.\dEe]+)\s*(\S*)")


def _continuation(text: str, level: Optional[Level], gamma: Optional[Gamma]) -> None:
    for item in text.split("$"):
        m = _ITEM.match(item)
        if not m:
            continue
        if m.group(1) and gamma is not None:
            n = number(m.group(2), m.group(3))
            if n is not None:
                gamma.weisskopf[m.group(1)[1:]] = n
        elif m.group(4) and level is not None and gamma is None:
            n = number(m.group(5), m.group(6))
            if n is not None:
                level.moments[m.group(4)[3:]] = n
                if m.group(5)[0] not in "+-":
                    level.unsigned.add(m.group(4)[3:])
        elif m.group(7) and gamma is not None:
            try:
                gamma.final_energy = float(m.group(8))
            except ValueError:  # a final level such as "495.5+X" cannot be placed by energy
                pass
        elif m.group(9) and gamma is not None and gamma.conversion is None:
            gamma.conversion = number(m.group(10), m.group(11))


def nucid(z: int, a: int) -> str:
    """ENSDF's five-character nuclide identifier: ``"194PT"``, ``" 58NI"``."""
    return f"{a:>3d}{SYMBOLS[z].upper():<2s}"


def parse(lines: Iterator[str], z: int, a: int, mass_mev: Optional[float] = None) -> Optional[Adopted]:
    """The adopted dataset of nuclide (Z, A) from the lines of an ENSDF file, or ``None`` if it is not there.

    ``mass_mev`` (the nuclear mass) lets the recoil energy be subtracted when a γ ray is matched to its final
    level; without it A × 931.494 MeV is used.
    """
    want = nucid(z, a)
    mass_kev = (mass_mev or a * 931.494) * 1e3
    data: Optional[Adopted] = None
    inside = False
    level: Optional[Level] = None
    gamma: Optional[Gamma] = None
    for raw in lines:
        line = raw.rstrip("\r\n").ljust(80)
        if not line.strip():  # a blank line ends a dataset
            if inside:
                break
            continue
        if not inside:
            if line[:5] == want and line[5:9] == "    " and line[9:39].strip().upper().startswith("ADOPTED LEVELS"):
                inside = True
                data = Adopted(f"{a}{SYMBOLS[z]}", line[39:65].strip(), line[74:80].strip(), [])
            continue
        if line[:5] != want or line[6] != " ":  # another nuclide's record, or a comment
            continue
        kind, primary = line[7], line[5] in " 1"
        if kind == "H":  # the history record names the publication: CIT=NDS 177, 1 (2021)
            cit = re.search(r"CIT=([^$]+)", line[9:])
            if cit and not data.reference:
                data.reference = cit.group(1).strip()
            continue
        if line[8] != " ":  # a particle record (delayed particles), not a level or γ ray
            continue
        if kind == "L" and primary:
            gamma = None
            e = number(line[9:19], line[19:21])
            if e is None:  # an energy such as "1229.5+X": the level cannot be placed
                level = None
                continue
            t_text = line[39:49].strip()
            level = Level(e, line[21:39].strip(), half_life(t_text, line[49:55]), t_text.upper() == "STABLE")
            data.levels.append(level)
        elif kind == "G" and primary:
            gamma = None
            e = number(line[9:19], line[19:21])
            if level is None or e is None:
                continue
            gamma = Gamma(e, number(line[21:29], line[29:31]), line[31:41].strip(),
                          number(line[41:49], line[49:55]), number(line[55:62], line[62:64]))
            level.gammas.append(gamma)
        elif kind in "LG" and not primary:
            _continuation(line[9:], level, gamma if kind == "G" else None)
        elif primary:
            gamma = None
    if data is None:
        return None
    _place(data, mass_kev)
    return data


def _place(data: Adopted, mass_kev: float) -> None:
    """Match each γ ray to the level it feeds: the one named by ``FL=``, else the one nearest to the level energy
    less the γ-ray energy and the recoil energy E²/(2Mc²)."""
    energies = [lv.energy[0] for lv in data.levels]
    for i, lv in enumerate(data.levels):
        kept = []
        for g in lv.gammas:
            e = g.energy[0]
            target = g.final_energy if g.final_energy is not None else lv.energy[0] - e - e * e / (2 * mass_kev)
            below = [k for k in range(len(energies)) if k != i and energies[k] < lv.energy[0]]
            if below:
                k = min(below, key=lambda n: abs(energies[n] - target))
                tolerance = max(1.0, 3 * ((g.energy[1] or 0.0) + (lv.energy[1] or 0.0)
                                          + (data.levels[k].energy[1] or 0.0)))
                if abs(energies[k] - target) <= tolerance:
                    g.final = k
                    kept.append(g)
                    continue
            data.unplaced += 1
        lv.gammas = kept


# ---------------------------------------------------------------------------------------------------------------
# The local copy


def _sources(a: int) -> Iterator[Iterator[str]]:
    """Line iterators over the files of the local copy that may hold mass A: files and zip members whose name ends
    in the mass number (``ensdf.194``) first, then every other file."""
    d = folder()
    suffix = f".{a:03d}"
    paths = sorted(p for p in d.rglob("*") if p.is_file() and p.name != "README.md")
    plain = [p for p in paths if p.suffix.lower() != ".zip"]
    zips = [p for p in paths if p.suffix.lower() == ".zip"]

    def read_path(p: Path) -> Iterator[str]:
        with open(p, encoding="latin-1") as f:
            yield from f

    def read_member(zp: Path, name: str) -> Iterator[str]:
        with zipfile.ZipFile(zp) as zf, zf.open(name) as f:
            for raw in f:
                yield raw.decode("latin-1")

    members = [(zp, n) for zp in zips for n in zipfile.ZipFile(zp).namelist() if not n.endswith("/")]
    for exact in (True, False):
        for p in plain:
            if p.name.endswith(suffix) == exact:
                yield read_path(p)
        for zp, n in members:
            if n.endswith(suffix) == exact and not n.lower().endswith(".zip"):
                yield read_member(zp, n)


@lru_cache(maxsize=64)
def _adopted(z: int, a: int, where: str) -> Adopted:
    from . import data as _data

    try:
        mass = _data.nuclide((z, a)).nuclear_mass_mev
    except ValueError:
        mass = None
    for lines in _sources(a):
        found = parse(lines, z, a, mass)
        if found is not None:
            return found
    raise EnsdfMissing(f"the local ENSDF copy in {where} has no adopted levels for {a}{SYMBOLS[z]}")


def adopted(nuclide: str) -> Adopted:
    """The adopted levels and γ rays of a nuclide from the local copy. Raises :class:`EnsdfMissing`."""
    z, a = parse_nuclide(nuclide)
    if not available():
        raise EnsdfMissing(f"there is no local copy of ENSDF in {folder()}; download one with "
                           f"'python scripts/fetch_ensdf.py', or type the level scheme in the setup file")
    return _adopted(z, a, str(folder()))


__all__ = ["ENV", "Adopted", "EnsdfMissing", "Gamma", "Level", "adopted", "available", "folder", "half_life",
           "nucid", "number", "parse", "spin_parity"]
