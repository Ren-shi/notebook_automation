"""Level schemes: the levels, γ-ray transitions and reduced matrix elements of a nucleus (backlog item 50).

A scheme is read from the local ENSDF copy (:mod:`physim.nuclear.ensdf`) or typed in a setup file, and every value
records where it came from::

    from physim.nuclear.levels import LevelScheme

    pt = LevelScheme.from_ensdf("194Pt", max_energy_kev=1500)
    pt.levels[1].energy                       # Value(328.464, 0.012, "ensdf")
    pt.matrix_element(0, 1, "E2")             # Value(..., "derived"): ⟨2⁺‖E2‖0⁺⟩ in e fm²
    pt.b(0, 1, "E2")                          # B(E2; 0⁺ → 2⁺) in e² fm⁴
    pt.set_matrix_element(0, 1, "E2", "1.28 eb")   # the user's value takes precedence from now on

Sources: ``"ensdf"`` (as evaluated), ``"derived"`` (computed here from ENSDF values), ``"assumed"`` (a choice made
because the data do not say, such as a sign) and ``"user"``.

Conventions: energies in keV, half-lives in seconds, reduced matrix elements ⟨f‖M(Eλ)‖i⟩ in e fm^λ, with
B(Eλ; i → f) = ⟨f‖M(Eλ)‖i⟩² / (2J_i + 1). ENSDF gives no signs: derived matrix elements are positive, and the
sign is marked as assumed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Optional, Union

from . import ensdf
from .names import parse_nuclide
from .quantity import Quantity

SOURCES = ("ensdf", "derived", "assumed", "user")
#: Electric multipolarities that can carry a matrix element.
ELECTRIC = ("E1", "E2", "E3")

#: e² in MeV fm, ħ in MeV s and ħc in MeV fm (CODATA 2018).
E2_MEV_FM = 1.439964547
HBAR_MEV_S = 6.582119569e-22
HBARC_MEV_FM = 197.3269804
LN2 = math.log(2.0)

#: Units of a reduced matrix element and their value in e fm^λ, by λ.
ME_UNITS = {1: {"efm": 1.0, "eb0.5": 10.0}, 2: {"efm2": 1.0, "eb": 100.0}, 3: {"efm3": 1.0, "eb1.5": 1000.0}}


@dataclass(frozen=True)
class Value:
    """A number with its uncertainty (or ``None``), where it came from, and a note on how."""

    value: float
    unc: Optional[float] = None
    source: str = "ensdf"
    note: str = ""

    def __post_init__(self):
        if self.source not in SOURCES:
            raise ValueError(f"source must be one of {', '.join(SOURCES)}, got '{self.source}'")


@dataclass
class Level:
    energy: Value
    #: Spin and parity (+1 or −1), or ``None`` where ENSDF does not give a single value.
    spin: Optional[float] = None
    parity: Optional[int] = None
    #: The assignment as written in ENSDF, e.g. ``"(2+)"``.
    jpi: str = ""
    half_life: Optional[Value] = None

    @property
    def label(self) -> str:
        if self.spin is None:
            return self.jpi or "?"
        j = str(int(self.spin)) if float(self.spin).is_integer() else f"{int(2 * self.spin)}/2"
        return j + {1: "+", -1: "-", None: ""}[self.parity]


@dataclass
class Transition:
    """A γ-ray transition from level ``initial`` down to level ``final`` (indices into the scheme's levels)."""

    initial: int
    final: int
    energy: Value
    #: Relative photon intensity among the γ rays of the initial level.
    intensity: Optional[Value] = None
    multipolarity: str = ""
    mixing: Optional[Value] = None
    #: Total internal conversion coefficient α.
    conversion: Optional[Value] = None
    #: Reduced transition probabilities in Weisskopf units as ENSDF gives them, by multipolarity.
    weisskopf: dict = field(default_factory=dict)


@dataclass
class MatrixElement:
    """⟨b‖M(Eλ)‖a⟩ between levels ``a`` and ``b`` (a ≤ b; a = b is a static moment), in e fm^λ."""

    a: int
    b: int
    multipolarity: str
    value: Value


# ---------------------------------------------------------------------------------------------------------------
# Formulas


def weisskopf_unit(lam: int, a: int) -> float:
    """The Weisskopf single-particle estimate B_W(Eλ) = (1/4π) (3/(λ+3))² (1.2 A^⅓)^{2λ}, in e² fm^{2λ}."""
    return (3.0 / (lam + 3)) ** 2 * (1.2 * a ** (1 / 3)) ** (2 * lam) / (4 * math.pi)


def rate_per_b(lam: int, energy_kev: float) -> float:
    """γ-ray emission rate for unit B(Eλ↓), in 1/s per e² fm^{2λ}:
    T = 8π(λ+1) / (λ [(2λ+1)!!]²) · (e²/ħ) · (E_γ/ħc)^{2λ+1} · B(Eλ)."""
    double_factorial = math.prod(range(2 * lam + 1, 0, -2))
    return (8 * math.pi * (lam + 1) / (lam * double_factorial**2) * E2_MEV_FM / HBAR_MEV_S
            * (energy_kev * 1e-3 / HBARC_MEV_FM) ** (2 * lam + 1))


def quadrupole_factor(j: float) -> float:
    """Q_s / ⟨J‖M(E2)‖J⟩: the spectroscopic quadrupole moment (e fm²) per unit diagonal matrix element,
    √(16π/5) √(J(2J−1) / ((2J+1)(J+1)(2J+3))). For J = 2 it is 0.7579."""
    if j < 1:
        return 0.0
    return math.sqrt(16 * math.pi / 5) * math.sqrt(j * (2 * j - 1) / ((2 * j + 1) * (j + 1) * (2 * j + 3)))


def components(multipolarity: str) -> Optional[list]:
    """The multipoles of an ENSDF multipolarity, in the order written: ``"M1+E2"`` → ["M1", "E2"]; ``"(E2)"`` →
    ["E2"]. ``None`` when it does not name them (``"D+Q"``, ``"M1,E2"``, empty)."""
    text = multipolarity.replace("(", "").replace(")", "").replace("[", "").replace("]", "").strip()
    parts = text.split("+") if text else []
    if not 1 <= len(parts) <= 2 or not all(len(p) == 2 and p[0] in "EM" and p[1].isdigit() for p in parts):
        return None
    return parts


def fractions(multipolarity: str, mixing: Optional[float]) -> dict:
    """Share of the γ-ray intensity carried by each multipole: {"M1": 1/(1+δ²), "E2": δ²/(1+δ²)}. A mixed
    transition without a mixing ratio gives an empty result, as the shares are then unknown."""
    parts = components(multipolarity)
    if parts is None:
        return {}
    if len(parts) == 1:
        return {parts[0]: 1.0}
    if mixing is None:
        return {}
    d2 = mixing * mixing
    return {parts[0]: 1.0 / (1.0 + d2), parts[1]: d2 / (1.0 + d2)}


def parse_matrix_element(text: Union[str, float], lam: int) -> float:
    """A reduced matrix element from text with a unit (``"1.28 eb"``, ``"128 efm2"``), in e fm^λ."""
    if isinstance(text, (int, float)) and not isinstance(text, bool):
        return float(text)
    parts = str(text).split()
    units = ME_UNITS[lam]
    if len(parts) != 2 or parts[1] not in units:
        raise ValueError(f"write an E{lam} matrix element with its unit ({' or '.join(units)}), "
                         f"e.g. '1 {list(units)[-1]}'; got '{text}'")
    try:
        return float(parts[0]) * units[parts[1]]
    except ValueError:
        raise ValueError(f"cannot read '{parts[0]}' as a number") from None


# ---------------------------------------------------------------------------------------------------------------
# The scheme


@dataclass
class LevelScheme:
    nuclide: str
    levels: list = field(default_factory=list)
    transitions: list = field(default_factory=list)
    matrix_elements: list = field(default_factory=list)
    #: Where the scheme came from, for the citation ("ENSDF, adopted levels of 194Pt, NDS 177, 1 (2021)").
    reference: str = ""
    #: What could not be used, in plain words, for the user to read.
    notes: list = field(default_factory=list)

    @property
    def A(self) -> int:  # noqa: N802
        return parse_nuclide(self.nuclide)[1]

    # -- reading from ENSDF -----------------------------------------------------------------------------------------

    @classmethod
    def from_ensdf(cls, nuclide: str, max_energy_kev: Optional[float] = None,
                   connected_only: bool = True) -> LevelScheme:
        """The adopted levels of ``nuclide`` from the local ENSDF copy, with matrix elements derived, cut to the
        levels below ``max_energy_kev`` that Coulomb excitation can reach (see :meth:`select`). Raises
        :class:`~physim.nuclear.ensdf.EnsdfMissing` without a local copy."""
        scheme = cls.from_adopted(ensdf.adopted(nuclide))
        if not scheme.levels:
            raise ensdf.EnsdfMissing(f"ENSDF has no level of {scheme.nuclide} with a known energy")
        scheme.derive()
        if max_energy_kev is not None or connected_only:
            scheme = scheme.select(max_energy_kev, connected_only)
        return scheme

    @classmethod
    def from_adopted(cls, data: ensdf.Adopted) -> LevelScheme:
        """The scheme as ENSDF gives it, before anything is derived."""

        def val(n) -> Optional[Value]:
            return None if n is None else Value(n[0], n[1], "ensdf")

        scheme = cls(data.nuclide, reference=("ENSDF, adopted levels and γ rays of " + data.nuclide
                                              + (f", {data.reference}" if data.reference else "")
                                              + (f" (file dated {data.date})" if data.date else "")))
        for i, lv in enumerate(data.levels):
            j, parity, _ = ensdf.spin_parity(lv.jpi)
            scheme.levels.append(Level(val(lv.energy), j, parity, lv.jpi, val(lv.half_life)))
            for g in lv.gammas:
                scheme.transitions.append(Transition(
                    i, g.final, val(g.energy), val(g.intensity), g.multipolarity, val(g.mixing), val(g.conversion),
                    {k: val(v) for k, v in g.weisskopf.items()}))
            q = lv.moments.get("E2")
            if q is not None and j is not None and quadrupole_factor(j) > 0:
                f = quadrupole_factor(j)  # Q in barns → matrix element in e fm²
                scheme.matrix_elements.append(MatrixElement(i, i, "E2", Value(
                    q[0] * 100.0 / f, None if q[1] is None else q[1] * 100.0 / f, "derived",
                    f"from the quadrupole moment Q = {q[0]:g} b in ENSDF"
                    + ("; ENSDF gives no sign, positive assumed" if "E2" in lv.unsigned else ""))))
        if data.unplaced:
            scheme.notes.append(f"{data.unplaced} γ ray(s) were left out because no final level matched")
        return scheme

    # -- matrix elements --------------------------------------------------------------------------------------------

    def _find(self, a: int, b: int, multipolarity: str) -> Optional[MatrixElement]:
        a, b = min(a, b), max(a, b)
        for me in self.matrix_elements:
            if (me.a, me.b, me.multipolarity) == (a, b, multipolarity):
                return me
        return None

    def matrix_element(self, a: int, b: int, multipolarity: str = "E2") -> Optional[Value]:
        """⟨b‖M(Eλ)‖a⟩ in e fm^λ, or ``None`` if the scheme has none for this pair. Its size does not depend on
        the order of a and b."""
        me = self._find(a, b, multipolarity)
        return None if me is None else me.value

    def set_matrix_element(self, a: int, b: int, multipolarity: str, value: Union[str, float],
                           unc: Optional[float] = None, source: str = "user", note: str = "") -> None:
        """Set a matrix element (text with a unit, or a number in e fm^λ). A value set here is kept by
        :meth:`derive`."""
        if multipolarity not in ELECTRIC:
            raise ValueError(f"matrix elements are kept for {', '.join(ELECTRIC)}, not '{multipolarity}'")
        for k in (a, b):
            if not 0 <= k < len(self.levels):
                raise ValueError(f"there is no level {k}; the scheme has {len(self.levels)} levels")
        v = Value(parse_matrix_element(value, int(multipolarity[1])), unc, source, note)
        me = self._find(a, b, multipolarity)
        if me is None:
            self.matrix_elements.append(MatrixElement(min(a, b), max(a, b), multipolarity, v))
        else:
            me.value = v

    def b(self, initial: int, final: int, multipolarity: str = "E2") -> Optional[float]:
        """B(Eλ; initial → final) = ⟨f‖M‖i⟩² / (2J_i + 1), in e² fm^{2λ}."""
        v = self.matrix_element(initial, final, multipolarity)
        j = self.levels[initial].spin
        if v is None or j is None:
            return None
        return v.value**2 / (2 * j + 1)

    def b_weisskopf(self, initial: int, final: int, multipolarity: str = "E2") -> Optional[float]:
        """B(Eλ; initial → final) in Weisskopf units."""
        b = self.b(initial, final, multipolarity)
        return None if b is None else b / weisskopf_unit(int(multipolarity[1]), self.A)

    def quadrupole_moment(self, level: int) -> Optional[float]:
        """Spectroscopic quadrupole moment of a level, e fm², from its diagonal E2 matrix element."""
        v = self.matrix_element(level, level, "E2")
        j = self.levels[level].spin
        return None if v is None or j is None else v.value * quadrupole_factor(j)

    def branches(self, level: int) -> list:
        """The γ-ray transitions leaving a level."""
        return [t for t in self.transitions if t.initial == level]

    def derive(self) -> None:
        """Fill in the size of each E1, E2 and E3 matrix element that ENSDF's values determine. Values set by
        the user are kept. Each transition uses, in order of preference:

        1. ENSDF's B(Eλ) in Weisskopf units: B↓ = B_W × the Weisskopf unit;
        2. the level's half-life with its branching: the γ-ray rate of the transition is
           (ln 2 / T½) × I_γ / Σ_k I_k (1 + α_k), of which the multipole carries the share its mixing ratio gives.

        Then ⟨f‖M‖i⟩ = √((2J_i + 1) B↓), taken positive: the sign is not in the data.
        """
        for t in self.transitions:
            up, low = self.levels[t.initial], self.levels[t.final]
            if up.spin is None or low.spin is None:
                continue
            for mult in ELECTRIC:
                existing = self._find(t.initial, t.final, mult)
                if existing is not None and existing.value.source == "user":
                    continue
                found = self._b_down(t, mult)
                if found is None:
                    continue
                b_down, rel, how = found
                m = math.sqrt((2 * up.spin + 1) * b_down)
                value = Value(m, None if rel is None else m * rel / 2, "derived", how + "; sign assumed positive")
                if existing is None:
                    self.matrix_elements.append(MatrixElement(min(t.initial, t.final), max(t.initial, t.final),
                                                              mult, value))
                else:
                    existing.value = value

    def _b_down(self, t: Transition, mult: str) -> Optional[tuple]:
        """(B(Eλ↓) in e² fm^{2λ}, its relative uncertainty or None, how it was found) for one transition."""
        lam = int(mult[1])
        w = t.weisskopf.get(mult)
        if w is not None and w.value > 0:
            return (w.value * weisskopf_unit(lam, self.A), None if w.unc is None else w.unc / w.value,
                    f"from B({mult}) = {w.value:g} W.u. in ENSDF")
        share = fractions(t.multipolarity, None if t.mixing is None else t.mixing.value).get(mult)
        t_half = self.levels[t.initial].half_life
        branches = self.branches(t.initial)
        if not share or t_half is None or t.intensity is None or t.intensity.value <= 0:
            return None
        if any(k.intensity is None for k in branches):
            return None
        total = sum(k.intensity.value * (1 + (k.conversion.value if k.conversion else 0.0)) for k in branches)
        rate = LN2 / t_half.value * t.intensity.value / total * share
        rel = None if t_half.unc is None else t_half.unc / t_half.value
        return (rate / rate_per_b(lam, t.energy.value), rel,
                f"from the half-life {t_half.value:.3g} s and the γ-ray branching in ENSDF")

    # -- truncation -------------------------------------------------------------------------------------------------

    def connected(self) -> set:
        """Levels that Coulomb excitation can reach from the ground state: joined to it by a chain of E1, E2 or E3
        matrix elements."""
        reach, grew = {0}, True
        while grew:
            grew = False
            for me in self.matrix_elements:
                if me.a != me.b and (me.a in reach) != (me.b in reach):
                    reach |= {me.a, me.b}
                    grew = True
        return reach

    def select(self, max_energy_kev: Optional[float] = None, connected_only: bool = True) -> LevelScheme:
        """A scheme cut down to the levels at or below ``max_energy_kev`` that excitation can reach, together with
        every level their γ rays feed (so branching ratios stay complete). Levels are renumbered."""
        keep = set(range(len(self.levels)))
        if connected_only:
            keep &= self.connected()
        if max_energy_kev is not None:
            keep = {i for i in keep if self.levels[i].energy.value <= max_energy_kev}
        keep.add(0)
        grew = True
        while grew:  # levels fed by the γ rays of kept levels
            fed = {t.final for t in self.transitions if t.initial in keep} - keep
            grew = bool(fed)
            keep |= fed
        order = sorted(keep)
        new = {old: k for k, old in enumerate(order)}
        notes = list(self.notes)
        if len(order) < len(self.levels):
            notes.append(f"{len(order)} of {len(self.levels)} levels kept"
                         + (f", up to {max_energy_kev:g} keV" if max_energy_kev is not None else ""))
        return LevelScheme(
            self.nuclide, [self.levels[i] for i in order],
            [replace(t, initial=new[t.initial], final=new[t.final]) for t in self.transitions
             if t.initial in keep and t.final in keep],
            [MatrixElement(new[m.a], new[m.b], m.multipolarity, m.value) for m in self.matrix_elements
             if m.a in keep and m.b in keep],
            self.reference, notes)

    # -- the setup file ---------------------------------------------------------------------------------------------

    def to_dict(self) -> dict:
        """The scheme as a ``[levels.beam]`` or ``[levels.target]`` section of a setup file."""

        def put(d: dict, key: str, v: Optional[Value], unit: str = "", scale: float = 1.0) -> None:
            if v is None:
                return
            d[key] = f"{_num(v.value * scale)} {unit}" if unit else v.value * scale
            if v.unc is not None:
                d[key + "_unc"] = v.unc * scale

        out: dict[str, Any] = {"nuclide": self.nuclide}
        if self.reference:
            out["reference"] = self.reference
        out["level"] = []
        for lv in self.levels:
            d: dict[str, Any] = {}
            put(d, "energy", lv.energy, "keV")
            if lv.jpi or lv.spin is not None:
                d["jpi"] = lv.jpi or lv.label
            put(d, "half_life", lv.half_life, "s")
            d["source"] = lv.energy.source
            out["level"].append(d)
        if self.transitions:
            out["transition"] = []
        for t in self.transitions:
            d = {"from": t.initial, "to": t.final}
            put(d, "energy", t.energy, "keV")
            put(d, "intensity", t.intensity)
            if t.multipolarity:
                d["multipolarity"] = t.multipolarity
            put(d, "mixing", t.mixing)
            put(d, "conversion", t.conversion)
            for mult, w in t.weisskopf.items():
                put(d, f"b_{mult.lower()}_wu", w)
            d["source"] = t.energy.source
            out["transition"].append(d)
        if self.matrix_elements:
            out["matrix_element"] = []
        for m in self.matrix_elements:
            lam = int(m.multipolarity[1])
            d = {"from": m.a, "to": m.b, "multipolarity": m.multipolarity}
            put(d, "value", m.value, list(ME_UNITS[lam])[0])
            d["source"] = m.value.source
            if m.value.note:
                d["note"] = m.value.note
            out["matrix_element"].append(d)
        return out

    @classmethod
    def from_dict(cls, raw: Any, where: str, problems: list) -> Optional[LevelScheme]:
        """Read a ``[levels.*]`` section. Problems are appended, each naming the field; returns ``None`` if any
        were found."""
        before = len(problems)
        if not isinstance(raw, dict):
            problems.append(f"{where}: must be a table with a nuclide and [[{where}.level]] sections")
            return None
        nuclide = raw.get("nuclide")
        try:
            parse_nuclide(nuclide)
        except ValueError as e:
            problems.append(f"{where}: nuclide {e}")
            nuclide = None

        def value(d: dict, key: str, label: str, kind: Optional[str] = None, unit: str = "",
                  required: bool = False) -> Optional[Value]:
            if key not in d:
                if required:
                    problems.append(f"{label}: '{key}' is missing")
                return None
            x, source = d[key], d.get("source", "user")
            if source not in SOURCES:
                problems.append(f"{label}: source must be one of {', '.join(SOURCES)}, got '{source}'")
                source = "user"
            try:
                if kind is not None:
                    q = Quantity.parse(x)
                    if q.kind != kind:
                        raise ValueError(f"'{x}' is not written in {unit} or a unit like it")
                    number = q.to(unit)
                elif isinstance(x, bool) or not isinstance(x, (int, float)):
                    raise ValueError(f"must be a number, got {x!r}")
                else:
                    number = float(x)
            except ValueError as e:
                problems.append(f"{label}: {key} {e}")
                return None
            unc = d.get(key + "_unc")
            if unc is not None and (isinstance(unc, bool) or not isinstance(unc, (int, float)) or unc < 0):
                problems.append(f"{label}: {key}_unc must be a number that is not negative")
                unc = None
            return Value(number, unc, source, d.get("note", "") if key == "value" else "")

        def rows(key: str) -> list:
            items = raw.get(key, [])
            if not isinstance(items, list) or not all(isinstance(x, dict) for x in items):
                problems.append(f"{where}: {key} must be written as [[{where}.{key}]] sections")
                return []
            return items

        scheme = cls(nuclide or "", reference=raw.get("reference", ""))
        for i, d in enumerate(rows("level")):
            label = f"{where} level {i}"
            energy = value(d, "energy", label, "energy", "keV", required=True)
            jpi = d.get("jpi", "")
            j, parity, _ = ensdf.spin_parity(jpi) if isinstance(jpi, str) else (None, None, False)
            if not isinstance(jpi, str):
                problems.append(f"{label}: jpi must be text such as \"2+\"")
            t_half = value(d, "half_life", label, "time", "s")
            if energy is not None:
                if scheme.levels and energy.value <= scheme.levels[-1].energy.value:
                    problems.append(f"{label}: levels must be listed in order of increasing energy")
                scheme.levels.append(Level(energy, j, parity, jpi if isinstance(jpi, str) else "", t_half))
        if not scheme.levels:
            problems.append(f"{where}: no levels; add [[{where}.level]] sections, the ground state first")
        elif scheme.levels[0].energy.value != 0:
            problems.append(f"{where}: the first level must be the ground state, with energy \"0 keV\"")
        n = len(scheme.levels)

        def pair(d: dict, label: str) -> Optional[tuple]:
            a, b = d.get("from"), d.get("to")
            if not all(isinstance(k, int) and not isinstance(k, bool) and 0 <= k < n for k in (a, b)):
                problems.append(f"{label}: 'from' and 'to' must be level numbers from 0 to {n - 1}")
                return None
            return a, b

        for i, d in enumerate(rows("transition")):
            label = f"{where} transition {i}"
            ends = pair(d, label)
            energy = value(d, "energy", label, "energy", "keV", required=True)
            if ends is None or energy is None:
                continue
            if ends[0] <= ends[1]:
                problems.append(f"{label}: a transition goes from a higher level to a lower one")
                continue
            scheme.transitions.append(Transition(
                ends[0], ends[1], energy, value(d, "intensity", label), str(d.get("multipolarity", "")),
                value(d, "mixing", label), value(d, "conversion", label),
                {m: v for m in ("E1", "E2", "E3", "M1", "M2") if (v := value(d, f"b_{m.lower()}_wu", label))}))
        for i, d in enumerate(rows("matrix_element")):
            label = f"{where} matrix element {i}"
            ends = pair(d, label)
            mult = d.get("multipolarity", "E2")
            if mult not in ELECTRIC:
                problems.append(f"{label}: multipolarity must be one of {', '.join(ELECTRIC)}")
                continue
            if ends is None:
                continue
            source = d.get("source", "user")
            try:
                number = parse_matrix_element(d.get("value"), int(mult[1]))
            except ValueError as e:
                problems.append(f"{label}: value: {e}")
                continue
            if source not in SOURCES:
                problems.append(f"{label}: source must be one of {', '.join(SOURCES)}, got '{source}'")
                continue
            scheme.matrix_elements.append(MatrixElement(min(ends), max(ends), mult, Value(
                number, d.get("value_unc"), source, str(d.get("note", "")))))
        return scheme if len(problems) == before else None

    # -- summaries --------------------------------------------------------------------------------------------------

    def first_excitation(self, multipolarity: str = "E2") -> Optional[tuple]:
        """The lowest level joined to the ground state by a matrix element of this multipolarity, as
        (level number, energy in keV, B(Eλ↑) in e² fm^{2λ})."""
        for i in range(1, len(self.levels)):
            b = self.b(0, i, multipolarity)
            if b:
                return i, self.levels[i].energy.value, b
        return None

    def table(self) -> dict:
        """The scheme as rows for display: levels, transitions and matrix elements, each with its source."""
        levels = [{"n": i, "energy_kev": lv.energy.value, "jpi": lv.label,
                   "half_life_s": lv.half_life.value if lv.half_life else None, "source": lv.energy.source}
                  for i, lv in enumerate(self.levels)]
        transitions = []
        for t in self.transitions:
            total = sum(k.intensity.value for k in self.branches(t.initial) if k.intensity)
            transitions.append({
                "from": t.initial, "to": t.final, "energy_kev": t.energy.value, "multipolarity": t.multipolarity,
                "branching": t.intensity.value / total if t.intensity and total > 0 else None,
                "mixing": t.mixing.value if t.mixing else None,
                "conversion": t.conversion.value if t.conversion else None, "source": t.energy.source})
        elements = []
        for m in self.matrix_elements:
            lam = int(m.multipolarity[1])
            low, high = self.levels[m.a], self.levels[m.b]
            row = {"from": m.a, "to": m.b, "multipolarity": m.multipolarity, "value_efm": m.value.value,
                   "unc_efm": m.value.unc, "source": m.value.source, "note": m.value.note,
                   "b_up_e2fm": None, "b_up_wu": None, "q_efm2": None}
            if m.a != m.b and low.spin is not None:
                row["b_up_e2fm"] = m.value.value**2 / (2 * low.spin + 1)
                row["b_up_wu"] = row["b_up_e2fm"] / weisskopf_unit(lam, self.A)
            elif m.a == m.b and high.spin is not None and lam == 2:
                row["q_efm2"] = m.value.value * quadrupole_factor(high.spin)
            elements.append(row)
        return {"nuclide": self.nuclide, "reference": self.reference, "notes": list(self.notes),
                "levels": levels, "transitions": transitions, "matrix_elements": elements}


def _num(x: float) -> str:
    return f"{x:.10g}"


__all__ = ["ELECTRIC", "ME_UNITS", "SOURCES", "Level", "LevelScheme", "MatrixElement", "Transition", "Value",
           "components", "fractions", "parse_matrix_element", "quadrupole_factor", "rate_per_b", "weisskopf_unit"]
