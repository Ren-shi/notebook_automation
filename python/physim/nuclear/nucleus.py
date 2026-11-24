"""The nucleus above the experiments: what is known about one nucleus, by every method and across its experiments.

The results ledger (:mod:`physim.nuclear.results`) holds one row per result. This module is the layer between it
and the screens: it brings the evaluation in as rows, takes data points, and gives the nuclei and one nucleus as
plain data (dicts and lists) that a screen, the report and the MCP tools can all show::

    from physim.nuclear import nucleus

    nucleus.snapshot_evaluation("62Ni")              # ENSDF's adopted values into the nucleus's ledger
    nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 906.0, "e2fm4", direction="up", uncertainty=37.0,
                           method="coulex-fit", reference={"authors": "Allmond et al.", "year": 2014,
                                                           "doi": "10.1103/PhysRevC.90.034309"})
    nucleus.nuclei()                                  # the nuclei, their experiments and what is still open
    nucleus.view("62Ni")                              # one nucleus: its experiments and its observables as rows
    nucleus.observable_view("62Ni", "B(E2)", "2+1 > 0+1")   # the comparison, the routes, what the value gives

**The evaluation as rows** (:func:`evaluation`, :func:`snapshot_evaluation`). The adopted values of the local
ENSDF copy become ``evaluated`` rows: level energies, mean lives, quadrupole moments and B(E2) values, each with
its uncertainty and a reference that names ENSDF, the publication and the data set's date. The snapshot writes
them to the nucleus's ledger; taking it again adds nothing, and a value that a newer ENSDF copy changed gets a new
row that supersedes the old one.

**Data points** (:func:`add_data_point`). Real data comes in as points: an observable of a level or transition
with its value and uncertainty, a method, and a reference or "own measurement". The level or transition is given
as one says it ("2+1", "2+ at 1173 keV", "2+1 > 0+1"; :func:`level`, :func:`transition`). Everything is checked
and :class:`PointProblem` says all that is wrong at once. A point is never changed: a correction is a new row that
supersedes it (:func:`correct_data_point`), and so is a retraction (:func:`retract_data_point`).

**The views** (:func:`nuclei`, :func:`view`, :func:`observable_view`). The evidence stays in its two groups:
independent (measured, evaluated, theory, outside program) and projected (a simulated run analysed back, compared
only with the value the simulation was given). For one observable of one level or transition the view holds the
comparison of :mod:`physim.nuclear.compare`, the routes of the registry (:mod:`physim.nuclear.observables`) each
marked done, projected (simulated runs only), in physim or later and tested for this nucleus (:func:`applies`),
the registry's conversions applied to the reference value, the checks against other observables and the
calculations (:func:`checks`: plain ratios of independent rows, and the rotor reading), what a value is used for,
and the history.

**No formula of its own.** Every number here comes from a function that existed before: the level scheme's own
derivation (:class:`~physim.nuclear.levels.LevelScheme`), :func:`physim.nuclear.lifetimes.level_taus`, and the
conversions of the registry and the comparison engine. Three things are left out because the code has no
function for them: the branching of a level with several γ rays (so B(E2) ↔ τ is converted for levels with one γ ray
only, :func:`decay`), the uncertainty of a matrix element taken from a B value (:func:`gives`), and that of a ratio
of two rows (:func:`checks` shows each value with its own).
"""

from __future__ import annotations

import logging
import math
import re
from datetime import date as _date
from functools import lru_cache
from pathlib import Path
from typing import Optional, Union

from . import compare as _compare
from . import ensdf, lifetimes, observables, results, runs
from .gamma import level_names
from .levels import ME_UNITS, LevelScheme, components
from .names import parse_nuclide
from .quantity import Quantity
from .results import METHODS, PROJECTED, LevelKey, Row, TransitionKey

_log = logging.getLogger(__name__)

#: The name the evaluation's rows carry in their origin (``origin["evaluation"]``).
ENSDF = "ENSDF"

#: The evidence a data point can be: independent evidence that someone brings in. The evaluation comes from the
#: snapshot and projected rows from the analysis of a simulated run.
POINT_EVIDENCE = ("measured", "outside program", "theory")
#: The methods each evidence of a data point allows (a key of :data:`physim.nuclear.results.METHODS`).
POINT_METHODS = {
    "measured": tuple(m for m in METHODS if m not in ("shell-model", "evaluation", "converted")),
    "outside program": tuple(m for m in METHODS if m not in ("evaluation", "converted")),
    "theory": ("shell-model", "other"),
}

#: What a route of the registry can be for a nucleus: an independent row by that method holds now; only
#: simulated runs by that method hold (no measurement); physim has that kind of experiment; the kind comes later;
#: or the route is no experiment of the planner (a calculation, the evaluation, a method known from the literature
#: only: its value comes in as a data point).
DONE, ONLY_PROJECTED, IN_PHYSIM, LATER, NO_EXPERIMENT = "done", "projected", "in physim", "later", "no experiment"

#: The observables whose value must be above 0, and those that cannot be negative.
POSITIVE, NOT_NEGATIVE = ("tau",), ("B(E2)", "B(M1)", "E_level")

_SUBSCRIPT = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
#: What people write in a spin and parity or a level's name, as the keys write it.
_SAID = str.maketrans({"⁺": "+", "⁻": "-", "−": "-", "–": "-", "π": "", "₀": "0", "₁": "1", "₂": "2", "₃": "3",
                       "₄": "4", "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9"})
_ENERGY = re.compile(r"([-+]?\d+(?:\.\d+)?)\s*(keV|MeV)\b", re.IGNORECASE)
_GROUND = re.compile(r"(?<![\w.])(?:g\.?\s?s\.?|ground(?:\s+state)?)(?![\w.])", re.IGNORECASE)
_ARROW = re.compile(r"\s*(?:->|→|>|\bto\b)\s*")
_MULTIPOLE = re.compile(r"[EM][1-6]")


# ---------------------------------------------------------------------------------------------------------------
# A level and a transition, as one says them


def _kev(x) -> Optional[float]:
    """An energy in keV from a number (in keV) or text with a unit ("1172.9 keV", "1.173 MeV")."""
    if x is None:
        return None
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return float(x)
    try:
        return float(x)
    except (TypeError, ValueError):
        return float(Quantity.parse(str(x)).to("keV"))


def _of_nuclide(key, name: str):
    if key.nuclide != name:
        raise ValueError(f"'{key}' is of {key.nuclide}, not of {name}")
    return key


def level(nuclide: str, spec, check: bool = True) -> LevelKey:
    """The key of a level of ``nuclide`` from the way it is said.

    ``spec`` is its Jπ and place ("2+1": ENSDF's first 2⁺), its Jπ and energy ("2+ at 1172.9 keV", ``{"energy":
    "1172.9 keV", "jpi": "2+"}``, ``(1172.9, "2+")``), its energy alone ("1172.9 keV", or a number in keV), "gs"
    for the ground state, a key's text ("62Ni:2+1", "62Ni:2+@1173") or a :class:`~physim.nuclear.results.LevelKey`.
    A level given by its energy is anchored to ENSDF's level within 2 keV, as
    :func:`physim.nuclear.results.level_key` does. A level named by its place must be one ENSDF has; with
    ``check`` False (for looking rows up, where the key may have been written on another computer) the name is
    taken as it is when the local ENSDF copy cannot confirm it. Raises ``ValueError`` with what to write instead.
    """
    name = results.nuclide_name(nuclide)
    if isinstance(spec, LevelKey):
        return _of_nuclide(spec, name)
    energy = jpi = None
    if isinstance(spec, dict):
        unknown = sorted(set(spec) - {"energy", "energy_kev", "jpi"})
        if unknown:
            raise ValueError(f"a level is given by its 'energy' and 'jpi'; '{unknown[0]}' is neither")
        energy, jpi = _kev(spec.get("energy", spec.get("energy_kev"))), spec.get("jpi")
    elif isinstance(spec, (tuple, list)) and len(spec) == 2:
        energy, jpi = _kev(spec[0]), spec[1]
    elif isinstance(spec, (int, float)) and not isinstance(spec, bool):
        energy = float(spec)
    else:
        text = str(spec).translate(_SAID).strip()
        if ":" in text:
            head, _, text = text.partition(":")
            if results.nuclide_name(head) != name:
                raise ValueError(f"'{spec}' is a level of {head.strip()}, not of {name}")
        ground = _GROUND.search(text)
        if ground is not None:  # "gs", or "0+ gs"
            energy, jpi = 0.0, (text[:ground.start()] + " " + text[ground.end():]).strip() or None
        else:
            try:
                named = LevelKey.parse(f"{name}:{text.replace(' ', '')}")
            except ValueError:
                named = None
            if named is not None and named.anchored:  # "2+1": Jπ and the place among the levels of that Jπ
                found = [e for e, j, p, n in results.ensdf_levels(name)
                         if (j, p, n) == (named.j2, named.parity, named.ordinal)]
                if found:
                    return LevelKey(name, named.j2, named.parity, named.ordinal, found[0], True)
                if not check:
                    return named
                why = ("ENSDF has no such level" if results.ensdf_levels(name)
                       else "there is no local copy of ENSDF")
                raise ValueError(f"'{text}' of {name} cannot be placed: {why}. Name the level by its energy and "
                                 f"Jπ, as in '2+ at 1173 keV'")
            if named is not None and named.energy_kev is None:  # "2+(1)": a calculation's own place
                return named
            if named is not None:  # "2+@1173"
                energy, jpi = named.energy_kev, (named.jpi if named.j2 is not None else None)
            else:
                m = _ENERGY.search(text)
                if m is None:
                    raise ValueError(f"cannot read '{spec}' as a level of {name}; write '2+1' (the first 2+), "
                                     f"'2+ at 1173 keV', the energy alone ('1173 keV') or 'gs'")
                energy = float(m.group(1)) * (1e3 if m.group(2).lower() == "mev" else 1.0)
                rest = re.sub(r"\b(?:at|level|state|the)\b|,", " ", text[:m.start()] + " " + text[m.end():])
                jpi = rest.strip() or None
    if energy is None:
        raise ValueError(f"a level of {name} needs its energy (or its place, as in '2+1'); got {spec!r}")
    if energy < 0:
        raise ValueError(f"a level's energy cannot be negative, got {energy:g} keV")
    spin = parity = None
    if jpi is not None and str(jpi).strip():
        spin, parity, _ = ensdf.spin_parity(str(jpi).translate(_SAID))
        if spin is None:
            raise ValueError(f"cannot read '{jpi}' as a spin and parity; write it like '2+', '3/2-' or '(2+)'")
    return results.level_key(name, energy, spin, parity)


def transition(nuclide: str, spec, multipolarity: str = "E2", check: bool = True) -> TransitionKey:
    """The key of a transition of ``nuclide`` from the way it is said: two levels, each as :func:`level` takes
    it, joined by an arrow ("2+1 > 0+1", "0+ gs -> 2+ at 1173 keV", "4+1 to 2+1"), given as a pair or as ``{"from":
    ..., "to": ...}``; a key's text ("62Ni:2+1>0+1:E2") or a
    :class:`~physim.nuclear.results.TransitionKey`. The upper level comes first in the key whichever way the
    two are said: the direction of a B value is stated with the value, not here."""
    name = results.nuclide_name(nuclide)
    if isinstance(spec, TransitionKey):
        return _of_nuclide(spec, name)
    if isinstance(spec, dict):
        ends = [spec.get("from", spec.get("upper")), spec.get("to", spec.get("lower"))]
        multipolarity = spec.get("multipolarity", multipolarity)
    elif isinstance(spec, (tuple, list)) and len(spec) == 2:
        ends = list(spec)
    else:
        parts = str(spec).translate(_SAID).split(":")
        if len(parts) > 1 and _MULTIPOLE.fullmatch(parts[-1].strip()):
            multipolarity = parts.pop().strip()
        if len(parts) > 1:
            head = parts.pop(0)
            if results.nuclide_name(head) != name:
                raise ValueError(f"'{spec}' is a transition of {head.strip()}, not of {name}")
        ends = _ARROW.split(":".join(parts).strip())
    if len(ends) != 2 or any(e is None or not str(e).strip() for e in ends):
        raise ValueError(f"a transition is two levels, as in '2+1 > 0+1' or '2+ at 1173 keV > gs'; got {spec!r}")
    a, b = level(name, ends[0], check), level(name, ends[1], check)
    if a.text == b.text:
        raise ValueError(f"a transition joins two levels; both ends are {a.level}")
    return results.transition_key(a, b, str(multipolarity))


def key_of(nuclide: str, observable: str, key, check: bool = True) -> Union[LevelKey, TransitionKey]:
    """The key of the level or transition that ``key`` names as one says it (:func:`level`, :func:`transition`),
    whichever the observable belongs to."""
    o = observables.get(observable)
    if o.belongs_to == "transition":
        multipolarity = o.multipolarity
        found = transition(nuclide, key, multipolarity or "E2", check)
        if multipolarity and found.multipolarity != multipolarity:
            raise ValueError(f"{o.label} belongs to an {multipolarity} transition, not {found.multipolarity}")
        return found
    return level(nuclide, key, check)


def level_label(key: LevelKey) -> str:
    """A level's name for a heading: "2+₁" (ENSDF's first 2⁺), "2+ at 1173 keV" (known by its energy)."""
    if key.anchored:
        return key.jpi + str(key.ordinal).translate(_SUBSCRIPT)
    if key.ordinal is not None:
        return f"{key.jpi} (calculated no. {key.ordinal})"
    return f"{key.jpi} at {key.energy_kev:.0f} keV"


def label(observable: str, key) -> str:
    """An observable of a level or transition for a heading, in the canonical convention: "B(E2; 0+₁ → 2+₁)↑",
    "τ(2+₁)"."""
    o = observables.get(observable)
    if o.belongs_to == "transition":
        k = key if isinstance(key, TransitionKey) else TransitionKey.parse(str(key))
        ends = (k.final, k.initial) if o.direction == observables.UP else (k.initial, k.final)
        return (f"{o.label[:-1]}; {level_label(ends[0])} → {level_label(ends[1])})"
                + observables.DIRECTIONS.get(o.direction, ""))
    return f"{o.label}({level_label(key if isinstance(key, LevelKey) else LevelKey.parse(str(key)))})"


# ---------------------------------------------------------------------------------------------------------------
# The evaluation as rows


@lru_cache(maxsize=16)
def _adopted_scheme(name: str, where: str) -> Optional[tuple]:
    try:
        data = ensdf.adopted(name)
    except ensdf.EnsdfMissing:
        return None
    scheme = LevelScheme.from_adopted(data)
    keys = tuple(results.level_key(name, lv.energy.value, lv.spin, lv.parity) for lv in scheme.levels)
    return data, scheme, keys


def adopted_scheme(nuclide: str) -> Optional[tuple]:
    """(the adopted data set, its level scheme as ENSDF gives it, the key of each level) from the local ENSDF
    copy; None without a copy or without the nuclide. ``where`` in the cache is the ENSDF folder, so a change of
    it is noticed."""
    if not ensdf.available():
        return None
    return _adopted_scheme(results.nuclide_name(nuclide), str(ensdf.folder()))


def evaluation(nuclide: str, max_energy_kev: Optional[float] = None) -> list:
    """ENSDF's adopted values of a nucleus as ``evaluated`` rows (nothing is written:
    :func:`snapshot_evaluation` does that). Raises :class:`~physim.nuclear.ensdf.EnsdfMissing` without a local
    copy or without the nuclide. ``max_energy_kev`` keeps the levels up to that energy.

    - **E_level**: each level's energy and uncertainty, in keV (not the ground state).
    - **tau**: the mean life in ps of each level with γ rays and a half-life, as
      :func:`physim.nuclear.lifetimes.level_taus` gives it; the half-life's relative uncertainty is kept. A
      half-life that ENSDF gives without an uncertainty is left out: that is how a limit ("3.1 PS GT"), an
      approximate or a systematic value arrives here (the reader keeps no mark of which), and it is no value to
      pull a measurement against.
    - **Qs**: the quadrupole moment in e b as ENSDF writes it, for the levels the level scheme takes it for.
    - **mu**: the magnetic moment in μ_N as ENSDF writes it.
    - **B(M1)**: ENSDF's B(M1)↓ in Weisskopf units, for the γ rays it states one for. None is derived.
    - **B(E2)**: for each γ ray between levels of known spin, ENSDF's B(E2)↓ in Weisskopf units as written; where
      ENSDF gives none, the B(E2)↓ in e² fm⁴ that the level scheme derives from the half-life and the branching
      (:meth:`LevelScheme.derive <physim.nuclear.levels.LevelScheme.derive>`), which depends on the τ row.
      None is derived from a half-life without an uncertainty, nor for a level with a γ ray the reader could not
      place (:attr:`ensdf.Level.dropped <physim.nuclear.ensdf.Level.dropped>`): its branching is not complete.

    Each row's reference names ENSDF, the publication and the data set's date. Where ENSDF writes an asymmetric
    uncertainty the reader keeps the larger side (:func:`physim.nuclear.ensdf.number`).
    """
    name = results.nuclide_name(nuclide)
    data, scheme, keys = _adopted(name)
    year = int(data.date[:4]) if data.date[:4].isdigit() else None
    origin = results.reference("ENSDF, adopted levels and γ rays", year, where=data.reference, evaluation=ENSDF,
                               dataset_date=data.date)
    rows, taken = [], set()

    def add(observable: str, key, value: float, unit: str, **more) -> Optional[Row]:
        if (observable, str(key)) in taken:  # two levels without a Jπ within a keV of each other: one name
            return None
        taken.add((observable, str(key)))
        rows.append(Row(name, observable, key, value, unit, "evaluation", "evaluated", origin=origin, **more))
        return rows[-1]

    inside = {i for i, lv in enumerate(scheme.levels)
              if max_energy_kev is None or lv.energy.value <= max_energy_kev}
    taus = lifetimes.level_taus(scheme)
    names = level_names(scheme.levels)
    tau_rows = {}
    for i in sorted(inside):
        lv, key = scheme.levels[i], keys[i]
        if lv.energy.value > 0:
            add("E_level", key, lv.energy.value, "keV", total=lv.energy.unc)
        tau, hl = taus.get(names[i]), lv.half_life
        # A half-life without an uncertainty is a limit or an approximate value: no evaluated row is made of it.
        if tau is not None and hl is not None and hl.unc is not None and tau[1].startswith("half-life"):
            tau_rows[i] = add("tau", key, tau[0], "ps", detail=f"from the half-life {hl.value:.4g} s in ENSDF",
                              total=tau[0] * hl.unc / hl.value)
        q = data.levels[i].moments.get("E2")
        if q is not None and scheme.matrix_element(i, i, "E2") is not None:
            add("Qs", key, q[0], "eb", total=q[1],
                detail="ENSDF gives no sign" if "E2" in data.levels[i].unsigned else "")
        mu = data.levels[i].moments.get("M1")
        if mu is not None:
            add("mu", key, mu[0], "mu_N", total=mu[1],
                detail="ENSDF gives no sign" if "M1" in data.levels[i].unsigned else "")
    for t in scheme.transitions:
        if t.initial not in inside:
            continue
        if scheme.levels[t.initial].spin is None or scheme.levels[t.final].spin is None:
            continue
        m1 = t.weisskopf.get("M1")
        if m1 is not None and m1.value > 0:  # as written; none is derived (the mixing ratio would be needed)
            add("B(M1)", results.transition_key(keys[t.initial], keys[t.final], "M1"), m1.value,
                observables.WEISSKOPF, direction=observables.DOWN, total=m1.unc)
        key = results.transition_key(keys[t.initial], keys[t.final], "E2")
        w = t.weisskopf.get("E2")
        if w is not None and w.value > 0:
            add("B(E2)", key, w.value, observables.WEISSKOPF, direction=observables.DOWN, total=w.unc)
            continue
        found = scheme._b_down(t, "E2")  # the level scheme's own derivation, as LevelScheme.derive uses it
        if found is None:
            continue
        b_down, relative, how = found
        if relative is None or data.levels[t.initial].dropped:
            continue  # from a limit of the half-life, or with a branching that is not complete
        tau_row = tau_rows.get(t.initial)
        add("B(E2)", key, b_down, "e2fm4", direction=observables.DOWN, total=b_down * relative,
            detail=how + "; the uncertainty is the half-life's alone",
            depends_on=(tau_row.id,) if tau_row is not None else ())
    return rows


def _adopted(name: str) -> tuple:
    found = adopted_scheme(name)
    if found is None:
        ensdf.adopted(name)  # raises EnsdfMissing with the reason: no local copy, or no such nuclide in it
        raise ensdf.EnsdfMissing(f"the local ENSDF copy has no adopted levels for {name}")
    return found


#: What makes two rows the same statement: everything but when and by which code the row was made.
_CONTENT = ("observable", "key", "value", "unit", "direction", "stat", "sys", "total", "plus", "minus", "method",
            "evidence", "detail", "origin", "depends_on")


def _same(a: Row, b: Row) -> bool:
    return all(getattr(a, name) == getattr(b, name) for name in _CONTENT)


def is_evaluation(row: Row) -> bool:
    """Whether a row is a snapshot of ENSDF's adopted values."""
    return row.evidence == "evaluated" and row.origin.get("evaluation") == ENSDF


def snapshot_evaluation(nuclide: str, home: Union[str, Path, None] = None,
                        max_energy_kev: Optional[float] = None) -> dict:
    """Take ENSDF's adopted values of a nucleus (:func:`evaluation`) into its ledger.

    A value the ledger already holds unchanged is not written again, so taking the snapshot twice adds nothing.
    A value that differs from the snapshot that holds now (a newer ENSDF copy: another number, uncertainty or
    data set) gets a new row that supersedes the old one. A value the new copy no longer has stays as it is.
    Returns {"nuclide", "reference", "added", "superseded", "unchanged"}: the rows written for the first time,
    those that replaced an earlier one, and how many were there already.
    """
    name = results.nuclide_name(nuclide)
    fresh = evaluation(name, max_energy_kev)
    held = results.read(results.nucleus_folder(name, home) / results.POINTS_FILE)
    have = {(r.observable, r.key): r for r in results.current(held) if is_evaluation(r)}
    stored, added, superseded, unchanged = {}, [], [], 0
    for row in fresh:  # the τ rows come before the B(E2) rows that depend on them
        made = row.id
        if row.depends_on:
            row = Row.from_dict({**row.to_dict(), "id": "",
                                 "depends_on": [stored.get(i, i) for i in row.depends_on]})
        old = have.get((row.observable, row.key))
        if old is not None and _same(old, row):
            stored[made] = old.id
            unchanged += 1
            continue
        if old is not None:
            row = Row.from_dict({**row.to_dict(), "id": "", "supersedes": old.id})
        stored[made] = row.id
        (superseded if old is not None else added).append(row)
    results.add_points(added + superseded, home)
    return {"nuclide": name, "reference": source(fresh[0]) if fresh else "",
            "dataset_date": fresh[0].origin.get("dataset_date") if fresh else None,
            "added": added, "superseded": superseded, "unchanged": unchanged}


# ---------------------------------------------------------------------------------------------------------------
# Data points


class PointProblem(ValueError):
    """A data point that cannot be taken. ``problems`` lists everything that is wrong with it."""

    def __init__(self, problems: list):
        self.problems = list(problems)
        super().__init__("; ".join(self.problems))


def _reference(ref, problems: list) -> Optional[dict]:
    if not isinstance(ref, dict):
        problems.append("a reference is given as {authors, year, doi, where}, where 'where' is the place in the "
                        "paper (\"Table II\")")
        return None
    unknown = sorted(set(ref) - {"authors", "year", "doi", "where"})
    if unknown:
        problems.append(f"a reference has authors, year, doi and where; '{unknown[0]}' is none of them")
    authors = str(ref.get("authors") or "").strip()
    if not authors:
        problems.append("a reference needs its authors (\"Allmond et al.\")")
    year = ref.get("year")
    try:
        year = int(year)
        if not 1900 <= year <= _date.today().year + 1:
            raise ValueError
    except (TypeError, ValueError):
        problems.append(f"a reference needs the year it was published, got {ref.get('year')!r}")
        year = None
    return results.reference(authors, year, str(ref.get("doi") or "").strip(), str(ref.get("where") or "").strip())


def _own(own, problems: list) -> Optional[dict]:
    if not isinstance(own, dict):
        problems.append("an own measurement is given as {date, note}: when it was taken (\"2027-06-14\") and a "
                        "line on where and how")
        return None
    unknown = sorted(set(own) - {"date", "note"})
    if unknown:
        problems.append(f"an own measurement has a date and a note; '{unknown[0]}' is neither")
    when = str(own.get("date") or "").strip()
    try:
        if not re.fullmatch(r"\d{4}(-\d{2}){0,2}", when):
            raise ValueError
        _date.fromisoformat(when + "-01-01"[len(when) - 4:])
    except ValueError:
        problems.append(f"an own measurement needs its date, written as 2027-06-14 (or 2027-06, or 2027); got "
                        f"{own.get('date')!r}")
    return results.own_measurement(when, str(own.get("note") or "").strip())


def _direction(direction, o, problems: list) -> Optional[str]:
    said = {"↑": observables.UP, "↓": observables.DOWN, "excitation": observables.UP, "decay": observables.DOWN}
    direction = said.get(direction, direction)
    if o.direction is None:
        if direction is not None:
            problems.append(f"{o.label} has no direction: leave it out")
        return None
    if direction not in observables.DIRECTIONS:
        problems.append(f"{o.label} needs its direction: 'up' (↑, from the lower level to the upper one, as "
                        f"Coulomb excitation gives it) or 'down' (↓, the decay)")
        return None
    return direction


#: The arguments of :func:`add_data_point` that hold the uncertainty, in its three forms.
_UNCERTAINTY = ("uncertainty", "stat", "sys", "plus", "minus")


def _above_zero(*amounts) -> bool:
    """Whether one of the uncertainties given is above 0. Something that is not a number counts: it is said as
    its own problem."""
    for x in amounts:
        try:
            if x is not None and (isinstance(x, bool) or float(x) > 0):
                return True
        except (TypeError, ValueError):
            return True
    return False


def _data_point(name: str, row_id: str, home) -> Row:
    """The data point of a nucleus's ledger that ``row_id`` names (an id or its first characters), if it holds
    now. Refuses rows that are not data points: those are replaced where they were made."""
    rows = results.read(results.nucleus_folder(name, home) / results.POINTS_FILE)
    found = [r for r in rows if r.id.startswith(str(row_id).strip())] if str(row_id).strip() else []
    if not found:
        raise ValueError(f"the ledger of {name} has no row '{row_id}'")
    if len(found) > 1:
        raise ValueError(f"'{row_id}' fits {len(found)} rows of {name}; give more of the id")
    row = found[0]
    later = next((r for r in rows if r.supersedes == row.id), None)
    if later is not None:
        raise ValueError(f"row {row.id} does not hold any more: it was "
                         + (f"withdrawn by {later.id} ({later.retracted})" if later.retracted
                            else f"replaced by {later.id}"))
    if row.retracted:
        raise ValueError(f"row {row.id} is a retraction, not a value")
    if not is_data_point(row):
        raise ValueError(f"row {row.id} is not a data point ({source(row)}): the evaluation is replaced by a new "
                         f"snapshot and a calculation by a new calculation")
    return row


def is_data_point(row: Row) -> bool:
    """Whether a row is a data point someone brought in: a published value or an own measurement."""
    return (row.evidence in POINT_EVIDENCE and not row.retracted and not is_evaluation(row)
            and row.origin.get("kind") in ("reference", "own measurement"))


def add_data_point(nuclide: str, observable: str, key, value, unit: Optional[str] = None, *,
                   direction: Optional[str] = None, uncertainty: Optional[float] = None,
                   stat: Optional[float] = None, sys: Optional[float] = None, plus: Optional[float] = None,
                   minus: Optional[float] = None, method: Optional[str] = None, evidence: str = "measured",
                   reference: Optional[dict] = None, own: Optional[dict] = None, detail: str = "",
                   corrects: Optional[str] = None, home: Union[str, Path, None] = None) -> Row:
    """Add a data point to a nucleus's ledger and return its row.

    - ``observable``: a key of the registry ("B(E2)", "tau", "Qs", "E_level").
    - ``key``: the level or transition, as one says it (:func:`level`, :func:`transition`).
    - ``value`` and ``unit``: a unit the observable allows (its canonical one when left out). A mean life is
      above 0; a B value and a level's energy are not negative.
    - ``direction``: "up" or "down" for a B value; a value in Weisskopf units is always "down".
    - the uncertainty, absolute and in the same unit, in one form: one total (``uncertainty``), statistical and
      systematic (``stat``, ``sys``), or ``plus`` and ``minus``. A measured point needs one, above 0.
    - ``method``: a key of :data:`physim.nuclear.results.METHODS` that fits the evidence (:data:`POINT_METHODS`).
    - ``evidence``: "measured" (a published value or one's own), "outside program" (GOSIA and the like) or
      "theory".
    - where it is from: ``reference`` = {authors, year, doi, where in the paper}, or ``own`` = {date, note} for
      an own measurement. One of the two.
    - ``detail``: a free line beside the method ("thin target, GOSIA").
    - ``corrects``: the id of the data point this one replaces. The old row is kept and superseded.

    Raises :class:`PointProblem` with everything that is wrong. A point that is in the ledger already, word for
    word, is refused too.
    """
    problems: list = []
    name = o = old = None
    try:
        name = results.nuclide_name(nuclide)
    except (ValueError, TypeError, KeyError) as e:
        problems.append(f"the nuclide: {e}")
    try:
        o = observables.get(observable)
    except ValueError as e:
        problems.append(str(e))
    if name is not None and o is not None:
        try:
            key = key_of(name, observable, key)
        except ValueError as e:
            problems.append(f"the {o.belongs_to}: {e}")
    try:
        if isinstance(value, bool):
            raise ValueError
        number = float(value)
        if not math.isfinite(number):
            raise ValueError
        if o is not None and o.key in POSITIVE and number <= 0:
            problems.append(f"{o.label} must be above 0, got {number:g}")
        if o is not None and o.key in NOT_NEGATIVE and number < 0:
            problems.append(f"{o.label} cannot be negative, got {number:g}")
        value = number
    except (TypeError, ValueError):
        problems.append(f"the value must be a number, got {value!r}")
    if o is not None:
        unit = unit or o.unit
        if unit not in o.units:
            problems.append(f"'{unit}' is not a unit of {o.label}; use {', '.join(o.units)}")
        direction = _direction(direction, o, problems)
        if unit == observables.WEISSKOPF and direction == observables.UP:
            problems.append(f"{o.label} in Weisskopf units is always the ↓ value: give it as 'down'")
    forms = [uncertainty is not None, stat is not None or sys is not None, plus is not None or minus is not None]
    if sum(forms) > 1:
        problems.append("give the uncertainty in one form: one total, or stat and sys, or plus and minus")
    elif (plus is None) != (minus is None):
        problems.append("an asymmetric uncertainty needs both plus and minus")
    for label_, x in zip(_UNCERTAINTY, (uncertainty, stat, sys, plus, minus)):
        try:
            if x is not None and (isinstance(x, bool) or not (math.isfinite(float(x)) and float(x) >= 0)):
                raise ValueError
        except (TypeError, ValueError):
            problems.append(f"the uncertainty '{label_}' must be a number that is not negative, got {x!r}")
    if evidence not in POINT_EVIDENCE:
        why = {"evaluated": ": the evaluation is taken from ENSDF (snapshot_evaluation)",
               PROJECTED: ": a projected row comes from the analysis of a simulated run, in its experiment"}
        problems.append(f"a data point is {', '.join(POINT_EVIDENCE[:-1])} or {POINT_EVIDENCE[-1]}, not "
                        f"'{evidence}'{why.get(evidence, '')}")
    elif evidence == "measured" and not any(forms):
        problems.append("a measured point needs its uncertainty: one total, or stat and sys, or plus and minus")
    elif evidence == "measured" and not _above_zero(uncertainty, stat, sys, plus, minus):
        problems.append("a measured point needs an uncertainty above 0")
    if method not in METHODS:
        problems.append(f"the method is one of {', '.join(METHODS)}; got {method!r}")
    elif evidence in POINT_METHODS and method not in POINT_METHODS[evidence]:
        problems.append(f"'{method}' ({METHODS[method]}) is not a method of a {evidence} point; use one of "
                        f"{', '.join(POINT_METHODS[evidence])}")
    if (reference is None) == (own is None):
        problems.append("say where the point is from: a reference {authors, year, doi, where} or an own "
                        "measurement {date, note}, one of the two")
        origin = None
    else:
        origin = _reference(reference, problems) if reference is not None else _own(own, problems)
    if corrects is not None and name is not None:
        try:
            old = _data_point(name, corrects, home)
        except ValueError as e:
            problems.append(f"the point to correct: {e}")
    if problems:
        raise PointProblem(problems)
    try:
        row = Row(name, observable, key, value, unit, method, evidence, origin=origin, direction=direction,
                  stat=stat, sys=sys, total=uncertainty, plus=plus, minus=minus, detail=str(detail or "").strip(),
                  supersedes=old.id if old is not None else None)
    except ValueError as e:
        raise PointProblem([str(e)]) from None
    held = results.current(results.read(results.nucleus_folder(name, home) / results.POINTS_FILE))
    twin = next((r for r in held if _same(r, row)), None)
    if twin is not None:
        raise PointProblem([f"this point is in the ledger of {name} already (row {twin.id})"])
    if not results.add_point(row, home):  # the same words in the same second as a row since replaced or withdrawn
        raise PointProblem([f"a row that says exactly this was written this very second (row {row.id}, replaced "
                            f"or withdrawn since); add the point again in a moment"])
    return row


def correct_data_point(nuclide: str, row_id: str, home: Union[str, Path, None] = None, **changes) -> Row:
    """Correct a data point: a new row with ``changes`` (in the words of :func:`add_data_point`: ``value``,
    ``unit``, ``key``, ``direction``, the uncertainty, ``method``, ``evidence``, ``reference`` or ``own``,
    ``detail``) that supersedes the old one, which is kept. A new uncertainty replaces the old one whole, in
    whichever form it is given; so does a new reference or own measurement."""
    name = results.nuclide_name(nuclide)
    old = _data_point(name, row_id, home)
    o = observables.get(old.observable)
    args = {"key": TransitionKey.parse(old.key) if o.belongs_to == "transition" else LevelKey.parse(old.key),
            "value": old.value, "unit": old.unit, "direction": old.direction, "uncertainty": old.total,
            "stat": old.stat, "sys": old.sys, "plus": old.plus, "minus": old.minus, "method": old.method,
            "evidence": old.evidence, "detail": old.detail}
    if old.origin.get("kind") == "reference":
        args["reference"] = {k: old.origin.get(k) for k in ("authors", "year", "doi", "where")}
    else:
        args["own"] = {k: old.origin.get(k) for k in ("date", "note")}
    unknown = sorted(set(changes) - set(args) - {"observable", "reference", "own"})
    if unknown:
        raise PointProblem([f"'{unknown[0]}' is not something a data point has"])
    if any(k in changes for k in _UNCERTAINTY):
        args.update(dict.fromkeys(_UNCERTAINTY))
    if "reference" in changes or "own" in changes:
        args.pop("reference", None)
        args.pop("own", None)
    if "key" not in changes and changes.get("observable", old.observable) != old.observable:
        raise PointProblem(["a point moved to another observable needs its level or transition ('key') again"])
    args.update(changes)
    return add_data_point(name, args.pop("observable", old.observable), args.pop("key"), args.pop("value"),
                          args.pop("unit"), corrects=old.id, home=home, **args)


def retract_data_point(nuclide: str, row_id: str, reason: str, home: Union[str, Path, None] = None) -> Row:
    """Withdraw a data point: a row that supersedes it and says why. Nothing is deleted: the point leaves what
    holds now and stays in the history, with the retraction beside it. Returns the retraction's row."""
    name = results.nuclide_name(nuclide)
    if not str(reason or "").strip():
        raise PointProblem(["a retraction needs its reason (\"entered for the wrong level\")"])
    row = _data_point(name, row_id, home).revised(retracted=str(reason).strip())
    results.add_point(row, home)
    return row


# ---------------------------------------------------------------------------------------------------------------
# Rows as a screen shows them


def source(row: Row) -> str:
    """Where a row is from, in a few words: "ENSDF adopted, NDS 204, 1 (2025)", "Allmond et al. 2014", "own
    measurement, 2027-06-14", "62Ni-on-194Pt · 03-beam-24h", "shell model, JUN45"."""
    o = row.origin
    kind = o.get("kind")
    if kind == "reference":
        if o.get("evaluation"):
            return f"{o['evaluation']} adopted" + (f", {o['where']}" if o.get("where") else "")
        return " ".join(str(x) for x in (o.get("authors"), o.get("year")) if x)
    if kind == "own measurement":
        return "own measurement" + (f", {o['date']}" if o.get("date") else "")
    if kind == "run":
        return f"{o.get('experiment')} · {o.get('run')}"
    if kind == "shell model":
        return "shell model" + (f", {o['interaction']}" if o.get("interaction") else "")
    return "converted from another row"


def value_text(value: float, plus: Optional[float] = None, minus: Optional[float] = None) -> str:
    """A value with its uncertainty as text, with the digits the uncertainty warrants (two significant ones):
    "884 ± 13", "2.080 ± 0.030", "915 +59 −52"; six significant digits for a value without an uncertainty."""
    widths = [u for u in (plus, minus) if u is not None and math.isfinite(u) and u > 0]
    if not widths:
        return f"{value:.6g}"
    digits = max(0, 1 - math.floor(math.log10(min(widths))))
    if plus == minus or plus is None or minus is None:
        return f"{value:.{digits}f} ± {widths[0]:.{digits}f}"
    return f"{value:.{digits}f} +{plus:.{digits}f} −{minus:.{digits}f}"


def cell(row: Row, via: Optional[dict] = None) -> dict:
    """One row as a screen shows it, in the observable's canonical convention (B(E2)↑ in e² fm⁴, τ in ps, Q_s
    in e fm², E in keV): {"id", "evidence", "method", "method_label", "source", "detail", "created", "value",
    "plus", "minus", "unit", "direction", "uncertainty_percent", "text", "truth", "converted", "as_given",
    "via"}. "converted" says how the value was brought to the convention, "as_given" is the row's own value, unit
    and direction, and "via" is the row a converted value was made from ({"id", "observable", "method",
    "source", "text"}; None for a row of the ledger). A row that cannot be brought to the convention (a B↓ between
    levels of unknown spin) keeps its own value and unit and says why in "problem"."""
    o = observables.get(row.observable)
    out = {"id": row.id, "observable": row.observable, "key": row.key, "evidence": row.evidence,
           "method": row.method, "method_label": METHODS[row.method], "source": source(row), "detail": row.detail,
           "created": row.created, "as_given": {"value": row.value, "unit": row.unit, "direction": row.direction},
           "via": None, "problem": ""}
    try:
        c = _compare.canonical(row)
        out.update(value=c.value, plus=c.plus, minus=c.minus, unit=c.unit, direction=c.direction,
                   uncertainty_percent=c.percent, truth=c.truth, converted=c.how)
    except (ValueError, TypeError) as e:
        out.update(value=row.value, plus=row.unc("plus"), minus=row.unc("minus"), unit=row.unit,
                   direction=row.direction, truth=row.truth, converted="", problem=str(e),
                   uncertainty_percent=None if row.relative is None else 100 * row.relative)
    out["text"] = (f"{value_text(out['value'], out['plus'], out['minus'])} {out['unit']}"
                   + (observables.DIRECTIONS.get(out["direction"], "") if o.direction else ""))
    made_from = (via or {}).get(row.id)
    if made_from is not None:
        out["via"] = {"id": made_from.id, "observable": made_from.observable, "method": made_from.method,
                      "method_label": METHODS[made_from.method], "source": source(made_from),
                      "text": f"{observables.get(made_from.observable).label} = "
                              f"{value_text(made_from.value, made_from.unc('plus'), made_from.unc('minus'))} "
                              f"{made_from.unit}"}
        out["source"] = out["via"]["source"]
    return out


# ---------------------------------------------------------------------------------------------------------------
# B(E2) from τ and τ from B(E2): through the level's one γ ray


def decay(nuclide: str, level_key: LevelKey) -> tuple:
    """(transition, inputs, why not): the one γ ray by which a level decays, from ENSDF's adopted γ rays, and
    what a B ↔ τ conversion needs besides the value ({"energy_kev", "branching", "conversion", "mixing"}, as
    :func:`physim.nuclear.compare.convert` takes them).

    Only a level with a single γ ray has an answer: its branching is 1. The γ rays are counted as the ENSDF
    file lists them: a level with a γ ray the reader could not place (no final level matched, or an energy it
    could not read) has more than the one it kept, and no answer. For a level with several γ rays the
    branching is the transition's share of all decays; the code works that out in three places but has no
    function for it, so (None, None, the reason) comes back and nothing is converted. The same when ENSDF does
    not have the level, names no multipolarity, or gives a mixed transition without its mixing ratio.
    """
    found = adopted_scheme(nuclide)
    if found is None:
        return None, None, f"there is no local ENSDF copy with the γ rays of {results.nuclide_name(nuclide)}"
    data, scheme, keys = found
    texts = [k.text for k in keys]
    if level_key.text not in texts:
        return None, None, f"ENSDF has no level {level_key.level}"
    i = texts.index(level_key.text)
    dropped = data.levels[i].dropped
    if dropped:
        return None, None, (f"ENSDF lists {dropped} more γ ray{'s' if dropped != 1 else ''} from the level "
                            f"{level_label(keys[i])} that could not be placed: its branching is not known")
    branches = scheme.branches(i)
    if len(branches) != 1:
        return None, None, (f"the level {level_label(keys[i])} has {len(branches)} γ rays: the conversion needs "
                            f"the branching, which the code has no function for yet"
                            if branches else f"ENSDF has no γ ray from the level {level_label(keys[i])}")
    t = branches[0]
    parts = components(t.multipolarity)
    if parts is None:
        return None, None, f"ENSDF does not name the multipolarity of the {t.energy.value:g} keV γ ray"
    mixing = None
    if len(parts) == 2:
        if set(parts) != {"M1", "E2"} or t.mixing is None:
            return None, None, (f"the {t.energy.value:g} keV γ ray is {t.multipolarity} and ENSDF gives no "
                                f"E2/M1 mixing ratio for it")
        mixing, parts = t.mixing.value, ["E2"]
    try:
        key = results.transition_key(keys[i], keys[t.final], parts[0])
    except ValueError as e:
        return None, None, str(e)
    return key, {"energy_kev": t.energy.value, "branching": 1.0,
                 "conversion": t.conversion.value if t.conversion is not None else 0.0, "mixing": mixing}, ""


def with_converted(nuclide: str, rows: list) -> tuple:
    """(rows, via): the rows that hold now, and with them the B(E2) values that the τ rows give and the τ values
    that the B(E2) rows give, through the level's one γ ray (:func:`decay`,
    :func:`physim.nuclear.compare.convert`). ``via`` maps each converted row's id to the row it was made from.

    Converted rows are not written to the ledger: they are made each time, say what they were made from and
    depend on it, so the comparison does not count them as independent of their source. The evaluation is not
    converted (it states its own B(E2) and τ), nor is a row that was itself converted, nor a row whose partner
    is in the ledger already: a lifetime analysis writes its B(E2) beside its τ, and the B(E2) depends on the τ.
    """
    out, via = list(rows), {}
    decays: dict = {}
    observable_of = {r.id: r.observable for r in rows}
    made_from = {(i, r.observable) for r in rows for i in r.depends_on}

    def linked(row: Row, to: str) -> bool:
        """Whether the ledger holds a row of ``to`` made from this row, or this row was made from one."""
        return (row.id, to) in made_from or any(observable_of.get(i) == to for i in row.depends_on)

    def decay_of(level_text: str) -> tuple:
        if level_text not in decays:
            decays[level_text] = decay(nuclide, LevelKey.parse(level_text))
        return decays[level_text]

    for row in rows:
        if row.evidence == "evaluated" or row.method == "converted" or row.retracted:
            continue
        try:
            if row.observable == "tau":
                key, inputs, _ = decay_of(row.key)
                if key is None or key.multipolarity != "E2" or linked(row, "B(E2)"):
                    continue
                new = _compare.convert(row, "B(E2)", key, **inputs)
            elif row.observable == "B(E2)":
                key, inputs, _ = decay_of(TransitionKey.parse(row.key).initial.text)
                if key is None or key.text != row.key or linked(row, "tau"):
                    continue
                new = _compare.convert(row, "tau", **inputs)
            else:
                continue
        except (ValueError, TypeError):
            continue
        out.append(new)
        via[new.id] = row
    return out, via


# ---------------------------------------------------------------------------------------------------------------
# The routes to an observable and what a value gives


_TIME_WORDS = {"Y": "y", "D": "d", "H": "h", "M": "min", "S": "s", "MS": "ms", "US": "µs", "NS": "ns", "PS": "ps",
               "FS": "fs"}


def _time_words(seconds: float) -> str:
    """A time in the largest of ENSDF's units that keeps it at 1 or more: "5.27 y", "2.6e+06 y", "12.7 h"."""
    units = sorted(((size, name) for name, size in ensdf.TIME_UNITS.items() if name in _TIME_WORDS), reverse=True)
    size, name = next(((s, n) for s, n in units if seconds >= s), units[-1])
    return f"{seconds / size:.3g} {_TIME_WORDS[name]}"


def _ps_words(tau_ps: float) -> str:
    """A mean life given in ps in the unit one says it in: "10 fs", "1 ps", "1 ns"."""
    for factor, unit in ((1e3, "ns"), (1.0, "ps")):
        if tau_ps >= factor:
            return f"{tau_ps / factor:g} {unit}"
    return f"{tau_ps * 1e3:g} fs"


def stability(nuclide: str) -> dict:
    """What the data say about a nucleus as a target or a beam: {"abundance" (its atom fraction in the natural
    element, from the NIST table of :func:`physim.nuclear.data.element`; None if it is not in it), "stable"
    (ENSDF's ground state is marked STABLE; None without a local ENSDF copy or without the nuclide),
    "half_life_s" (the ground state's half-life in ENSDF, or None), "words" (the same in a few words: "stable,
    3.63 % of natural Ni", "T½ = 5.27 y in ENSDF, not in the natural element")}."""
    from . import data
    from .names import SYMBOLS

    name = results.nuclide_name(nuclide)
    z, a = parse_nuclide(name)
    symbol = SYMBOLS[z]
    try:
        abundance = dict(data.element(symbol).isotopes).get(a)
    except (KeyError, ValueError):
        abundance = None
    stable = half_life = None
    found = adopted_scheme(name)
    if found is not None and found[0].levels:
        ground = found[0].levels[0]
        stable = bool(ground.stable)
        half_life = ground.half_life[0] if ground.half_life is not None else None
    known = "stable" if stable else f"T½ = {_time_words(half_life)} in ENSDF" if half_life is not None else ""
    share = f"{100 * abundance:.3g} % of natural {symbol}" if abundance is not None else "not in the natural element"
    return {"abundance": abundance, "stable": stable, "half_life_s": half_life,
            "words": f"{known}, {share}" if known else share if abundance is not None else
            f"{share}, and no half-life of its ground state is at hand"}


def _target_test(name: str, beam_too: bool) -> tuple:
    """Whether the nucleus can be a target (and, with ``beam_too``, whether it then matters): True for a nucleus
    of the natural element or one ENSDF marks stable. For another nucleus nothing in the code says how long-lived
    a target must be, nor which beams exist: None, with the half-life."""
    s = stability(name)
    if s["abundance"] is not None or s["stable"]:
        return True, (f"{name} ({s['words']}) can be the target or the beam" if beam_too
                      else f"a target of {name} ({s['words']}) can be made")
    return None, (f"{name}: {s['words']}. " + ("It would be the beam, or a target made for the purpose; whether "
                                               "either can be had is not checked" if beam_too else
                                               "The method needs a target of it; whether one can be made is not "
                                               "checked"))


def _is_ground(level_key: LevelKey) -> Optional[bool]:
    """Whether a level is the ground state: by ENSDF's first level where there is a local copy; else an anchored
    0⁺₁ of an even-even nucleus is, and for another level it is not known (None)."""
    known = results.ensdf_levels(level_key.nuclide)
    if known:
        e, j2, parity, n = known[0]
        return n is not None and level_key.anchored and (level_key.j2, level_key.parity,
                                                         level_key.ordinal) == (j2, parity, n)
    if level_key.energy_kev is not None and not level_key.anchored:
        return level_key.energy_kev == 0.0
    z, a = parse_nuclide(level_key.nuclide)
    if z % 2 == 0 and a % 2 == 0 and level_key.anchored:
        return (level_key.j2, level_key.parity, level_key.ordinal) == (0, 1, 1)
    return None


def _standing(rows: list, observable: str, key: str, evidence: tuple = ("evaluated", "measured",
                                                                      "outside program")) -> Optional[Row]:
    """The value of one observable of one level or transition that a test or a check reads: the evaluated row
    (the newest) where there is one, else the measured row with the smallest relative uncertainty. Never a
    projected row (a simulated run tests nothing) and, unless ``evidence`` names it, never a calculation."""
    found = [r for r in results.select(rows, observable, key) if r.evidence in evidence and not r.retracted]
    for kind in evidence:
        mine = [r for r in found if r.evidence == kind]
        if not mine:
            continue
        if kind == "evaluated":
            return max(mine, key=lambda r: r.created)
        judged = [r for r in mine if r.relative is not None]
        return min(judged, key=lambda r: r.relative) if judged else mine[-1]
    return None


def _lifetime_test(level_key: LevelKey, method: str, rows: Optional[list], via: Optional[dict]) -> tuple:
    """Whether the level's mean life is in the range physim plans a plunger or DSAM for
    (:func:`physim.nuclear.lifetime_plan.method_for`), with the evaluated τ of the ledger or, without one, the
    best measured one (a τ made from a B(E2) row counts). Never a calculated or a projected τ."""
    from . import lifetime_plan

    ranges = {"rdds": ("a plunger", lifetime_plan.RDDS_RANGE_PS), "dsam": ("DSAM", lifetime_plan.DSAM_RANGE_PS)}
    row = _standing(rows or [], "tau", level_key.text)
    if row is None:
        return None, (f"no τ of the level {level_label(level_key)} is entered (take the evaluation or add a data "
                      f"point): whether it is in the range of a plunger or of DSAM is not tested")
    tau = _compare.canonical(row).value
    what, (low, high) = ranges[method]
    made_from = (via or {}).get(row.id)
    known = f"τ = {tau:.3g} ps (" + (source(row) if made_from is None else
                                     f"from the {observables.get(made_from.observable).label} of "
                                     f"{source(made_from)}") + ")"
    span = f"{_ps_words(low)} to {_ps_words(high)}"
    suited = lifetime_plan.method_for(tau)
    if suited == method:
        return True, f"{known} is in the range physim plans {what} for ({span})"
    other = {"rdds": "a plunger", "dsam": "DSAM"}.get(suited)
    return False, (f"{known} is {'above' if tau >= high else 'below'} the range physim plans {what} for ({span})"
                   + (f"; it is in that of {other}" if other else ""))


@lru_cache(maxsize=64)
def _shell_model_fit(name: str, states: tuple, free_gb: int) -> tuple:
    """(True, False or None; the reason) for a shell-model calculation of ``states`` ((2J, parity, how many),
    ...) of a nucleus with ``free_gb`` GB of RAM free: :func:`physim.nuclear.shellmodel_spaces.estimate` for each
    interaction that holds the nucleus (:func:`~physim.nuclear.shellmodel_spaces.options`)."""
    from . import shellmodel as _sm
    from . import shellmodel_spaces as _spaces

    options = _spaces.options(name)
    if not options:
        return False, (f"none of the interactions physim has holds {name} in its valence space; a calculation "
                       f"run elsewhere can be imported")
    fits, too_big, cannot = [], [], []
    for opt in options:
        if any(p not in opt.parities for _, p, _ in states):
            cannot.append(f"{opt.spec.name} makes no state of that parity")
            continue
        try:
            est = _spaces.estimate(_sm.ShellModelRequest(name, opt.spec.name, states),
                                   free_bytes=free_gb * 1024**3)
        except ValueError as e:
            cannot.append(f"{opt.spec.name}: {e}")
            continue
        (fits if est.fits else too_big).append((est.memory_bytes, opt.spec.name, opt.spec.space.name, est.reason))
    if fits:
        _, interaction, space, reason = min(fits)
        others = [n for _, n, _, _ in sorted(fits)[1:]]
        return True, (f"{interaction} ({space} space): {reason}"
                      + (f" Also fitting: {', '.join(others)}." if others else "")
                      + (f" Too big now: {', '.join(n for _, n, _, _ in sorted(too_big))}." if too_big else ""))
    if too_big:
        _, interaction, space, reason = min(too_big)
        return False, f"{interaction} ({space} space), the smallest: {reason}"
    return False, "; ".join(cannot)


def _shell_model_test(name: str, levels: tuple) -> tuple:
    """Whether a shell-model calculation of the levels fits this computer now: the size estimate behind
    ``estimate_shell_model`` for a request of exactly these states (for 2⁺₁ → 0⁺₁, one 0⁺ and one 2⁺)."""
    from . import shellmodel_spaces as _spaces

    states: dict = {}
    for lv in levels:
        if lv.j2 is None or lv.parity is None:
            return None, (f"the spin and parity of the level {level_label(lv)} are not known: there is no "
                          f"calculation to size")
        states[(lv.j2, lv.parity)] = max(states.get((lv.j2, lv.parity), 0), lv.ordinal or 1)
    try:
        free = _spaces.available_memory_bytes()
    except OSError as e:
        return None, f"the memory free on this computer could not be read ({e}): the size is not tested"
    try:
        return _shell_model_fit(name, tuple(sorted((j2, p, n) for (j2, p), n in states.items())),
                                int(free // 1024**3))
    except (ValueError, KeyError, OSError) as e:
        return None, f"the size of a shell-model calculation could not be estimated ({e})"


def applies(nuclide: str, observable: str, key, route, rows: Optional[list] = None,
            via: Optional[dict] = None) -> tuple:
    """Whether a route applies to this nucleus: (True, False, or None for "not tested"; the reason in a line).

    Only tests the code already has are made, and the reason says which:

    - a route that needs the nucleus as a **target or a beam** (Coulomb excitation; photon, neutron and electron
      scattering, which need a target): True for a nucleus of the natural element or one ENSDF marks stable
      (:func:`stability`). For another nucleus nothing in the code says how long-lived a target must be or which
      beams exist: None, with the half-life;
    - the **Coulex yield** (first order) is of the excitation from the ground state, as the setup's
      ``[reaction]`` is: False for a transition between excited states; and a τ by Coulomb excitation needs the
      level's one γ ray (:func:`decay`);
    - a **lifetime** by a plunger or by DSAM: whether the level's τ is in the range physim plans the method for
      (:func:`physim.nuclear.lifetime_plan.method_for`), with the evaluated τ of ``rows`` or the best measured
      one; None without a τ;
    - the **shell model**: whether a calculation of exactly these states fits this computer now
      (:func:`physim.nuclear.shellmodel_spaces.estimate`, half of the RAM free); False where no interaction
      physim has holds the nucleus;
    - the **evaluation**: whether the local ENSDF copy has adopted levels for the nucleus;
    - Q_s of a level of spin 0 or 1/2 does not exist, by any route.

    A kind of experiment physim does not have yet (inelastic scattering, transfer) and a moment measurement have
    no test: None, and the reason says that nothing is checked. ``rows`` are the nucleus's rows that hold now
    (converted ones among them, ``via`` saying what each was made from). :func:`routes` puts the answer in each
    route as "applies" and "reason".
    """
    name = results.nuclide_name(nuclide)
    o = observables.get(observable)
    try:
        k = key if isinstance(key, (LevelKey, TransitionKey)) else key_of(name, observable, key, check=False)
    except ValueError as e:
        return None, f"the {o.belongs_to} could not be read ({e}): nothing is tested"
    transition_ = k if isinstance(k, TransitionKey) else None
    level_ = k.initial if transition_ is not None else k  # the level that decays: its τ is what a lifetime is of
    if observable == "Qs" and level_.spin is not None:
        try:
            observables.diagonal_matrix_element(1.0, level_.spin)
        except ValueError as e:
            return False, str(e)
    if route.method == "evaluation":
        if not ensdf.available():
            return None, "there is no local copy of ENSDF (scripts/fetch_ensdf.py): not tested"
        found = adopted_scheme(name)
        if found is None:
            return False, f"the local ENSDF copy has no adopted levels for {name}"
        return True, (f"the local ENSDF copy has the adopted levels of {name} ({found[0].reference}, data set of "
                      f"{found[0].date}); the snapshot takes them")
    if route.method == "shell-model":
        return _shell_model_test(name, (transition_.initial, transition_.final) if transition_ else (level_,))
    if route.kind == "lifetimes":
        return _lifetime_test(level_, route.method, rows, via)
    if route.kind == "coulex":
        if transition_ is None and observable == "tau":  # a τ from a B(E2): through the level's one γ ray
            transition_, _, why = decay(name, level_)
            if transition_ is None:
                return None, f"τ from a Coulex B(E2) is not converted here: {why}"
        if route.method == "coulex-yield" and transition_ is not None and _is_ground(transition_.final) is False:
            return False, ("the yield analysis is of the excitation from the ground state (the setup's "
                           "[reaction]); a transition between excited states needs the fit to all orders")
        return _target_test(name, beam_too=True)
    if route.method in ("nrf", "inelastic-neutron", "electron-scattering"):
        return _target_test(name, beam_too=False)
    if route.kind is not None and not route.live:
        from .workbench import KINDS

        kind = next(x.label for x in KINDS if x.key == route.kind)
        return None, f"physim has no {kind[0].lower() + kind[1:]} yet: nothing is tested for this route"
    return None, "physim has no test for this route: nothing is checked"


def routes(nuclide: str, observable: str, key: str, rows: list, via: Optional[dict] = None,
           every: Optional[list] = None, tests: bool = True) -> list:
    """The registry's routes to an observable, each with what it is for this level or transition.

    ``rows`` are the rows of the observable and key that hold now (converted ones among them, with ``via`` saying
    what each was made from: a B(E2) made from a DSAM lifetime counts for the DSAM route). ``every`` are all the
    nucleus's rows that hold now, for the tests of :func:`applies` (the level's τ is among them; ``rows`` when
    left out); with ``tests`` False no test is made and every route has "applies" None and no reason (a table of
    many rows needs the status only). Each route is a dict:
    "method" and "label", "kind" and "kind_label" (the kind of experiment that plans it, or None), "note", "status"
    (:data:`DONE`: an independent row by that method holds now; :data:`ONLY_PROJECTED`: only simulated runs
    by that method hold, so nothing is measured yet; :data:`IN_PHYSIM`: physim has that kind of experiment;
    :data:`LATER`: a kind still to come; :data:`NO_EXPERIMENT`: a calculation, the evaluation or a method known
    from the literature only), "live" (whether physim has the kind of experiment, whatever the status), "means"
    for a route that is no experiment ("calculation", "evaluation", "literature"), "rows", "independent" and
    "projected" (how many rows by that method hold now, in the two groups), "ids" (the ids of those rows), "own"
    (how many of them are one's own: from a run of an experiment or an own measurement, not from the literature),
    "open" (physim has the kind and none of the rows is one's own: the route is still worth planning, even where a
    published value by that method exists) and "applies" and "reason" (:func:`applies`; a route by which an
    independent value is here is never said not to apply: where physim's test says so, "applies" is None and the
    reason says both).
    """
    from .workbench import KINDS

    kinds = {k.key: k.label for k in KINDS}
    by_method: dict = {}
    for r in rows:
        by_method.setdefault((via or {}).get(r.id, r).method, []).append(r)
    out = []
    for route in observables.get(observable).routes:
        mine = by_method.get(route.method, [])
        live = route.live
        independent = sum(r.independent for r in mine)
        status = (DONE if independent else ONLY_PROJECTED if mine else IN_PHYSIM if live
                  else LATER if live is False else NO_EXPERIMENT)
        fits, why = (applies(nuclide, observable, key, route, rows if every is None else every, via) if tests
                     else (None, ""))
        if fits is False and independent:  # a value by this method is here: the test is of what physim plans
            fits, why = None, f"a value by this method is here, so the method reaches it; but {why}"
        own = sum((via or {}).get(r.id, r).origin.get("kind") in ("run", "own measurement") for r in mine)
        out.append({
            "method": route.method, "label": METHODS[route.method], "kind": route.kind,
            "kind_label": kinds.get(route.kind), "note": route.note, "status": status, "live": bool(live),
            "means": (None if route.kind is not None else "calculation" if route.method == "shell-model"
                      else "evaluation" if route.method == "evaluation" else "literature"),
            "rows": len(mine), "independent": independent, "projected": len(mine) - independent,
            "ids": [r.id for r in mine], "own": own, "open": bool(live) and own == 0, "applies": fits,
            "reason": why})
    return out


def _spin(key: LevelKey) -> float:
    if key.spin is None:
        raise ValueError(f"the spin of the level {key.level} is not known")
    return key.spin


def gives(nuclide: str, row: Row) -> list:
    """What a value gives: the conversions the registry offers for the row's observable, applied to the row (the
    reference of a comparison, usually).

    One dict per conversion: "to", "function" and "calls" (the registry's words: what it converts to, the function
    of :mod:`physim.nuclear.observables` and the older function behind it), "needs", then "value", "plus",
    "minus", "unit", "direction", "text" and "how"; "key" where the result belongs to another level or transition
    than the row. Where something needed is missing (a spin, the level's γ ray), "value" is None and "missing"
    says what.

    The uncertainty is carried where the conversion is a plain factor (B↓, Weisskopf units, the diagonal matrix
    element) or goes through :func:`physim.nuclear.compare.convert` (B(E2) ↔ τ). A matrix element from a B value
    is a square root: the code has no function for its uncertainty, so it comes without one.
    """
    c = _compare.canonical(row)
    name = results.nuclide_name(nuclide)

    def scaled(factor: float, unit: str, how: str, direction: Optional[str] = None) -> dict:
        return {"value": c.value * factor, "plus": None if c.plus is None else c.plus * abs(factor),
                "minus": None if c.minus is None else c.minus * abs(factor), "unit": unit, "direction": direction,
                "how": how}

    def of_row(new: Row) -> dict:
        n = _compare.canonical(new)
        return {"value": n.value, "plus": n.plus, "minus": n.minus, "unit": n.unit, "direction": n.direction,
                "how": new.detail, "key": new.key}

    def down(t: TransitionKey) -> tuple:
        ratio = observables.b_up_over_down(_spin(t.final), _spin(t.initial), t.multipolarity)
        return ratio, f"B↑ / {ratio:g} (J = {_spin(t.final):g} and {_spin(t.initial):g})"

    def b_up_over_down() -> dict:
        ratio, how = down(TransitionKey.parse(row.key))
        return scaled(1.0 / ratio, c.unit, how, observables.DOWN)

    def weisskopf() -> dict:
        t = TransitionKey.parse(row.key)
        ratio, how = down(t)
        factor = observables.weisskopf(1.0 / ratio, name, t.multipolarity)
        unit = observables.from_weisskopf(1.0, name, t.multipolarity)
        return scaled(factor, observables.WEISSKOPF, f"{how}; 1 W.u. = {unit:.4g} {c.unit} for {name}",
                      observables.DOWN)

    def matrix_element() -> dict:
        t = TransitionKey.parse(row.key)
        return {"value": observables.matrix_element(c.value, _spin(t.final), t.multipolarity), "plus": None,
                "minus": None, "unit": list(ME_UNITS[int(t.multipolarity[1])])[0], "direction": None,
                "how": "its size, from B↑ and the spin of the lower level; no uncertainty is worked out for it"}

    def tau_from_b() -> dict:
        t = TransitionKey.parse(row.key)
        key, inputs, why = decay(name, t.initial)
        if key is None or key.text != t.text:
            raise ValueError(why or f"the level decays by {key}, not by this transition")
        return of_row(_compare.convert(row, "tau", **inputs))

    def b_from_tau() -> dict:
        key, inputs, why = decay(name, LevelKey.parse(row.key))
        if key is None:
            raise ValueError(why)
        if key.multipolarity != "E2":
            raise ValueError(f"the level decays by an {key.multipolarity} γ ray, which has no B(E2)")
        return of_row(_compare.convert(row, "B(E2)", key, **inputs))

    def diagonal_matrix_element() -> dict:
        j = _spin(LevelKey.parse(row.key))
        factor = observables.diagonal_matrix_element(1.0, j)
        return scaled(factor, "efm2", f"Q_s × {factor:.4g} (J = {j:g})")

    apply = {f.__name__: f for f in (b_up_over_down, weisskopf, matrix_element, tau_from_b, b_from_tau,
                                    diagonal_matrix_element)}
    out = []
    for conv in observables.get(row.observable).conversions:
        entry = {"to": conv.to, "function": conv.function, "calls": conv.calls, "needs": list(conv.needs),
                 "value": None, "plus": None, "minus": None, "unit": "", "direction": None, "text": "", "how": "",
                 "missing": ""}
        try:
            entry.update(apply[conv.function]())
            entry["text"] = (f"{value_text(entry['value'], entry['plus'], entry['minus'])} {entry['unit']}"
                             + observables.DIRECTIONS.get(entry["direction"], ""))
        except KeyError:
            entry["missing"] = f"the nucleus view has no way to apply {conv.function} yet"
        except (ValueError, TypeError) as e:
            entry["missing"] = str(e)
        out.append(entry)
    return out


# ---------------------------------------------------------------------------------------------------------------
# Checks between observables


_UNIT_WORDS = {"e2fm4": "e²fm⁴", "efm2": "e fm²", "mu_N2": "μ_N²", "mu_N": "μ_N"}


def _said(row: Row, via: Optional[dict] = None) -> str:
    """A row's value and where it is from: "884 ± 13 e²fm⁴ ↑ (ENSDF adopted, NDS 204, 1 (2025))", with what a
    converted value was made from."""
    c = cell(row, via)
    arrow = observables.DIRECTIONS.get(c["direction"], "")
    return (f"{value_text(c['value'], c['plus'], c['minus'])} {_UNIT_WORDS.get(c['unit'], c['unit'])}"
            + (f" {arrow}" if arrow else "") + f" ({c['source']}"
            + (f", converted from {observables.get(c['via']['observable']).label}" if c["via"] is not None else "")
            + ")")


def _no_sign(row: Row) -> bool:
    """Whether the row is a moment that ENSDF gives without its sign (:func:`evaluation` says so in the detail)."""
    return row.evidence == "evaluated" and "no sign" in row.detail


def _check(label_: str, text: str, how: str, rows: list, value: float, reference: str = "") -> dict:
    return {"label": label_, "text": text, "how": how, "reference": reference,
            "rows": list(dict.fromkeys(r.id for r in rows)), "value": value}


#: How to read each check, and what it does not show.
_READ_E42 = ("2.0 is a harmonic vibrator, 2.5 a γ-soft rotor, 3.33 a rigid axial rotor; below 2 the nucleus is near a "
             "closed shell and few nucleons carry the states. It is a ratio of energies alone: it shows neither the "
             "sign of the deformation nor how collective the transitions are.")
_READ_B42 = ("10/7 = 1.43 is the rigid axial rotor and 2 the harmonic vibrator; below 1 the 4⁺₁ state does not continue "
             "the collective motion of the 2⁺₁ state (near a closed shell, or a state of another structure). Read it "
             "with E(4⁺₁)/E(2⁺₁): alone it does not tell a rotor from a nucleus in between. Each value has its own "
             "uncertainty; none is worked out for the ratio.")
_READ_ROTOR = ("For a rigid axial rotor B(E2) fixes the size of Q_s(2⁺₁), and the sign of the measured Q_s says prolate "
               "(negative) or oblate (positive). A measured Q_s much smaller than the rotor's says that the state is "
               "not that rotor (a harmonic vibrator has Q_s = 0, a soft or triaxial nucleus lies between): β₂ is then "
               "the size of the strength, not a shape. Agreement in size alone does not prove a rotor: look at "
               "E(4⁺₁)/E(2⁺₁) too. The rotor's numbers carry no uncertainty here.")
_READ_THEORY = {
    "B(E2)": ("1 means that the interaction, its valence space and its effective charges together give the strength; "
              "below 1 the calculation lacks collectivity (a space or effective charges too small), above 1 it has too "
              "much. The ratio does not say which of them is off, and a calculation has no uncertainty, so it is not "
              "judged."),
    "B(M1)": ("1 means that the interaction, its valence space and its effective g factors together give the strength. "
              "The ratio does not say which of them is off, and a calculation has no uncertainty, so it is not judged."),
    "tau": ("Above 1 the calculated transitions are too weak, below 1 too strong. It tests the interaction, its valence "
            "space and its effective charges together and does not say which is off; a calculation has no "
            "uncertainty, so it is not judged."),
    "Qs": ("A negative ratio means that the calculated state has the other shape (the other sign); near 1 the state is "
           "as the calculation has it, as far as this moment shows. It does not say whether the wave function or the "
           "effective charges are off, and a calculation has no uncertainty, so it is not judged."),
    "mu": ("A negative ratio means that the calculation has the other sign; near 1 the state holds the orbits the "
           "calculation gives it, as far as this moment shows. It does not say whether the wave function or the "
           "effective g factors are off, and a calculation has no uncertainty, so it is not judged."),
    "E_level": ("Far from 1, the calculation's wave functions are to be trusted less for this nucleus. Near 1 does not "
                "show that its B(E2) values and moments are right: an interaction is usually fitted to energies."),
}


def checks(nuclide: str, observable: str, key: str, rows: list, via: Optional[dict] = None,
           reference: Optional[Row] = None) -> list:
    """What the values that hold say against each other: the checks between the observable and others of the
    nucleus, and the calculations over the reference. Independent rows only: a projected one never takes part.

    ``rows`` are all the nucleus's rows that hold now (converted ones among them, ``via`` saying what each was made
    from) and ``reference`` the row of the comparison. Where a check needs one value of an observable it takes the
    evaluated one, else the measured one with the smallest uncertainty. Each check is a dict: "label", "text"
    (the numbers, each with its source), "how" (how to read it and what it does not show), "reference" (the book or
    paper the reading is from, "" for none), "rows" (the ids of the rows it used) and "value" (the ratio; for the
    rotor reading the prolate rotor's Q_s(2⁺₁) in e fm²).

    - **E(4⁺₁)/E(2⁺₁)** of an even-even nucleus: with its B(E2) and Q_s of the first states and these energies;
    - **B(E2; 4⁺₁ → 2⁺₁)/B(E2; 2⁺₁ → 0⁺₁)**, both ↓: with either B(E2);
    - **the rotor reading** of B(E2; 0⁺₁ → 2⁺₁) beside Q_s(2⁺₁), with either: the numbers are the Coulex analysis's
      (:func:`physim.nuclear.observables.rotor_reading`);
    - **each calculation over the reference** (a shell-model value over the evaluated one).

    The ratios are plain ones of two rows (:func:`physim.nuclear.compare.ratio`) and carry no uncertainty: the code
    has no function for it. A check that cannot be made is left out.
    """
    name = results.nuclide_name(nuclide)
    key = str(key)
    two, four = f"{name}:2+1", f"{name}:4+1"
    first, second = f"{two}>0+1:E2", f"{four}>2+1:E2"
    out: list = []

    def have(observable_: str, key_: str) -> Optional[Row]:
        return _standing(rows, observable_, key_)

    def energies() -> Optional[tuple]:
        a, b = have("E_level", four), have("E_level", two)
        return None if a is None or b is None else (_compare.ratio(a, b), a, b)

    def e42() -> Optional[dict]:
        found = energies()
        if found is None:
            return None
        value, a, b = found
        return _check("E(4⁺₁)/E(2⁺₁)", f"{value:.2f} = {_said(a, via)} over {_said(b, via)}", _READ_E42, [a, b], value,
                      observables.CASTEN)

    def b42() -> Optional[dict]:
        a, b = have("B(E2)", second), have("B(E2)", first)
        if a is None or b is None:
            return None
        # The rows compare as B↑ and the ratio is of the decays: B↓ = B↑ / b_up_over_down.
        value = _compare.ratio(a, b, observables.b_up_over_down(0, 2) / observables.b_up_over_down(2, 4))
        return _check("B(E2; 4⁺₁ → 2⁺₁)/B(E2; 2⁺₁ → 0⁺₁), both ↓",
                      f"{value:.2f}, from B(E2; 2⁺₁ → 4⁺₁) = {_said(a, via)} and B(E2; 0⁺₁ → 2⁺₁) = {_said(b, via)}",
                      _READ_B42, [a, b], value, observables.CASTEN)

    def rotor() -> Optional[dict]:
        b = have("B(E2)", first)
        r = observables.rotor_reading(_compare.canonical(b).value, name) if b is not None else {}
        if not r:
            return None
        used, q, found = [b], have("Qs", two), energies()
        prolate = r["qs_2plus_efm2"]
        text = (f"B(E2; 0⁺₁ → 2⁺₁) = {_said(b, via)} read as a rigid axial rotor: β₂ = {r['beta2']:.3f}, Q₀ = "
                f"{r['q0_efm2']:.0f} e fm² and Q_s(2⁺₁) = {_signed(prolate)} e fm² if prolate, {_signed(-prolate)} if "
                f"oblate. ")
        if q is None:
            text += "No Q_s(2⁺₁) is entered to set beside it."
        else:
            used.append(q)
            text += (f"{results.EVIDENCE[q.evidence]}: Q_s(2⁺₁) = {_said(q, via)}"
                     + (". ENSDF gives no sign for it: only its size can be set beside the rotor's."
                        if _no_sign(q) else "."))
        if found is not None:
            text += f" E(4⁺₁)/E(2⁺₁) = {found[0]:.2f}."
            used += found[1:]
        return _check("The rotor's Q_s(2⁺₁) from B(E2)", text, _READ_ROTOR, used, prolate,
                      observables.BOHR_MOTTELSON_2)

    def over(row: Row) -> dict:
        value = _compare.ratio(row, reference)
        how = _READ_THEORY[observable]
        if observable in ("Qs", "mu") and _no_sign(reference):
            how = "ENSDF gives no sign for the reference, so the sign of the ratio says nothing. " + how
        made_by = (via or {}).get(row.id, row)  # a calculated τ may come from a calculated B(E2)
        who = source(made_by)
        return _check(f"{who[:1].upper()}{who[1:]}, over the reference",
                      f"{value:.2f} = {_said(row, via)} over {_said(reference, via)}", how, [row, reference], value,
                      observables.RICHTER_2008 if made_by.method == "shell-model" else "")

    wanted = []
    z, a_ = parse_nuclide(name)
    if z % 2 == 0 and a_ % 2 == 0:
        if (observable, key) in (("B(E2)", first), ("B(E2)", second), ("Qs", two), ("E_level", two),
                                 ("E_level", four)):
            wanted.append(e42)
        if (observable, key) in (("B(E2)", first), ("B(E2)", second)):
            wanted.append(b42)
        if (observable, key) in (("B(E2)", first), ("Qs", two)):
            wanted.append(rotor)
    if reference is not None and reference.independent and reference.evidence != "theory":
        wanted += [lambda row=row: over(row) for row in results.select(rows, observable, key, "theory")
                   if row.id != reference.id]
    for make in wanted:
        try:
            found = make()
        except (ValueError, TypeError, AttributeError) as e:  # a value of 0, a B↓ between unknown spins, ...
            _log.debug("a check of %s of %s was left out: %s", observable, key, e)
            continue
        if found is not None:
            out.append(found)
    return out


def _signed(x: float) -> str:
    return f"{x:+.1f}".replace("-", "−")


# ---------------------------------------------------------------------------------------------------------------
# The views: the nuclei, one nucleus, one observable


def partner_nuclide(experiment) -> Optional[str]:
    """The other nucleus of a Coulomb-excitation setup: the one the ``[reaction]`` does not excite (the target's
    main nuclide when the beam is excited, the beam otherwise). None for another kind of reaction, and when beam
    and target are the same nuclide."""
    studied = results.reaction_nuclide(experiment)
    if studied is None:
        return None
    from . import data
    from .rates import beam_ion, stack

    try:
        beam = results.nuclide_name(beam_ion(experiment))
        target = results.nuclide_name(
            data.nuclide(max(stack(experiment)[0].nuclides.items(), key=lambda kv: kv[1])[0]).name)
    except (ValueError, KeyError, IndexError):
        return None
    other = target if studied == beam else beam
    return None if other == studied else other


def _experiments(root: Union[str, Path, None]) -> list:
    """Every experiment under ``root`` with what the nucleus views need: its folder and name, its kind, the nucleus
    it studies and its Coulex partner, its runs and which of them are analysed, and its ledger."""
    from .analysis import RESULT_FILE
    from .workbench import kind_of

    out = []
    for e in sorted(runs.list_experiments(root), key=lambda e: Path(e["path"]).name):
        folder = runs.ExperimentFolder(e["path"])
        try:
            setup = folder.setup()
        except (ValueError, OSError, KeyError, IndexError):
            setup = None
        kind, kind_label = kind_of(setup) if setup is not None else (None, "the setup cannot be read")
        ledger = results.of_experiment(folder.path)
        now = results.current(ledger)
        taken = [p.name for p in folder.run_folders()]
        # A run is analysed when the ledger holds a row from it, or its folder holds a kept Coulex analysis that
        # no row was made from yet (results.backfill makes them).
        analysed = {r.origin.get("run") for r in now if r.origin.get("kind") == "run"}
        analysed |= {p.name for p in folder.run_folders() if (p / RESULT_FILE).is_file()}
        key = None
        if setup is not None and setup.excitation is not None and setup.excitation.multipolarity == "E2":
            try:
                key = results.key_of_reaction(setup).text
            except ValueError:
                key = None
        out.append({"folder": folder.path.name, "path": str(folder.path), "name": folder.name, "kind": kind,
                    "kind_label": kind_label,
                    "studied": results.studied_nuclide(setup) if setup is not None else None,
                    "partner": partner_nuclide(setup) if setup is not None else None, "runs": len(taken),
                    "analysed_runs": len(analysed & set(taken)), "key": key, "ledger": ledger, "now": now})
    return out


def _unassigned(e: dict) -> str:
    """Why an experiment belongs to no nucleus."""
    if e["kind"] is None:
        return "the setup cannot be read"
    if e["kind"] == "elastic":
        return "elastic scattering studies no level of a nucleus"
    return ("the setup does not mark the exit channel it studies ([run] channel), so the experiment does not "
            "belong to a nucleus yet")


def _summary(cells: list) -> Optional[dict]:
    """The rows of one kind of evidence in a table row: how many, the range of the values, the best (the smallest
    relative uncertainty; the newest where none has one) and the newest."""
    if not cells:
        return None
    judged = [c for c in cells if c["uncertainty_percent"] is not None]
    newest = max(reversed(cells), key=lambda c: c["created"])  # of two written in one second, the later one
    return {"count": len(cells), "min": min(c["value"] for c in cells), "max": max(c["value"] for c in cells),
            "best": min(judged, key=lambda c: c["uncertainty_percent"]) if judged else newest, "newest": newest}


def table(nuclide: str, rows: list, extra: tuple = ()) -> list:
    """The observables of a nucleus as rows, one per observable and level or transition.

    ``rows`` are the rows that hold now; the B(E2) and τ values they give each other are added
    (:func:`with_converted`). ``extra`` names (observable, key) pairs to list even without a row (the transition an
    experiment studies before any run is analysed). Each row of the table is a dict:

    - "observable", "key", "label" ("B(E2; 0+₁ → 2+₁)↑"), "unit" (the canonical one) and "rows" (how many hold);
    - "evaluated": the newest evaluated row as a :func:`cell`, or None;
    - "measured", "outside_program", "theory", "projected": None, or {"count", "min", "max", "best", "newest"}
      with the best and the newest as cells; for the projected rows the best one's "uncertainty_percent" is the
      precision the setup reaches;
    - "open": the methods of the routes physim has an experiment for and that no row of one's own (a run or an
      own measurement) has taken yet;
    - "headline": whether there is more than the evaluation (or nothing at all yet): what a short list shows.

    The independent evidence and the projected rows are in separate entries and are never merged. Headline rows
    come first, then B(E2), τ, Q_s and energies, each by level energy.
    """
    name = results.nuclide_name(nuclide)
    rows, via = with_converted(name, rows)
    groups: dict = {}
    for r in rows:
        groups.setdefault((r.observable, r.key), []).append(r)
    for observable, key in extra:
        groups.setdefault((observable, str(key)), [])
    energies = {(j, p, n): e for e, j, p, n in results.ensdf_levels(name) if n is not None}

    def energy(observable: str, key: str) -> float:
        k = TransitionKey.parse(key).initial if observables.get(observable).belongs_to == "transition" \
            else LevelKey.parse(key)
        found = energies.get((k.j2, k.parity, k.ordinal)) if k.anchored else k.energy_kev
        return math.inf if found is None else found

    order = list(observables.OBSERVABLES)
    out = []
    for (observable, key), group in groups.items():
        cells = [cell(r, via) for r in group]
        kinds = {ev: [c for c in cells if c["evidence"] == ev] for ev in results.EVIDENCE}
        evaluated = _summary(kinds["evaluated"])
        out.append({
            "observable": observable, "key": key, "label": label(observable, key),
            "unit": observables.get(observable).unit, "rows": len(group),
            "evaluated": None if evaluated is None else evaluated["newest"],
            "measured": _summary(kinds["measured"]),
            "outside_program": _summary(kinds["outside program"]), "theory": _summary(kinds["theory"]),
            "projected": _summary(kinds[PROJECTED]),
            "open": [r["method"] for r in routes(name, observable, key, group, via, tests=False)
                     if r["open"]],
            "headline": not group or any(c["evidence"] != "evaluated" for c in cells)})
    return sorted(out, key=lambda t: (not t["headline"], order.index(t["observable"]),
                                      energy(t["observable"], t["key"]), t["key"]))


def _entry(name: str, experiments: list, home, full: bool) -> dict:
    """One nucleus of :func:`nuclei` (``full`` False: the headline rows of the table only) or :func:`view`."""
    mine = []
    for e in experiments:
        own = [r for r in e["now"] if r.nuclide == name]
        role = ("studied" if e["studied"] == name else "partner" if e["partner"] == name
                else "results" if own else None)
        if role is None:
            continue
        latest = max(reversed(own), key=lambda r: r.created) if own else None
        if latest is not None:
            latest = {**cell(latest), "label": label(latest.observable, latest.key)}
        mine.append({"folder": e["folder"], "path": e["path"], "name": e["name"], "kind": e["kind"],
                     "kind_label": e["kind_label"], "role": role, "runs": e["runs"],
                     "analysed_runs": e["analysed_runs"], "rows": len(own), "latest": latest})
    points = results.read(results.nucleus_folder(name, home) / results.POINTS_FILE)
    every = sorted([r for e in experiments for r in e["ledger"] if r.nuclide == name] + points,
                   key=lambda r: r.created)
    now = results.current(every)
    studied = [x for x in mine if x["role"] == "studied"]
    rows = table(name, now, tuple(("B(E2)", e["key"]) for e in experiments if e["studied"] == name and e["key"]))
    taken = [r for r in now if is_evaluation(r)]
    still = [{"what": "route", "observable": t["observable"], "key": t["key"], "label": t["label"],
              "methods": t["open"], "text": f"{t['label']}: " + ", ".join(METHODS[m] for m in t["open"])}
             for t in rows if t["headline"] and t["open"]]
    waiting = sum(x["runs"] - x["analysed_runs"] for x in studied)
    if waiting:
        still.append({"what": "runs", "count": waiting,
                      "text": f"{waiting} run{'s' if waiting != 1 else ''} not analysed"})
    if not taken and ensdf.available():
        still.append({"what": "evaluation", "text": "the adopted values of ENSDF are not taken yet"})
    z, a = parse_nuclide(name)
    return {"nuclide": name, "z": z, "a": a, "n": a - z, "experiments": mine,
            "runs": sum(x["runs"] for x in studied), "analysed_runs": sum(x["analysed_runs"] for x in studied),
            "rows": len(now), "points": sum(is_data_point(r) for r in now),
            "evaluation": {"rows": len(taken), "reference": source(taken[-1]) if taken else None,
                           "dataset_date": taken[-1].origin.get("dataset_date") if taken else None,
                           "taken": max(r.created for r in taken) if taken else None},
            "table": rows if full else [t for t in rows if t["headline"]], "open": still}


def _names(experiments: list, home) -> list:
    """The nuclei that have an experiment, a row of one, or a ledger of their own, lightest first."""
    found = {x for e in experiments for x in (e["studied"], e["partner"]) if x}
    found |= {r.nuclide for e in experiments for r in e["ledger"]}
    base = Path(home).expanduser() if home is not None else results.nuclei_home()
    for p in sorted(base.iterdir()) if base.is_dir() else []:
        if (p / results.POINTS_FILE).is_file():
            try:
                found.add(results.nuclide_name(p.name))
            except (ValueError, KeyError):
                continue
    return sorted(found, key=parse_nuclide)


def nuclei(root: Union[str, Path, None] = None, home: Union[str, Path, None] = None) -> dict:
    """The nuclei, each with its experiments, what is known and what is still open: the first screen.

    ``root`` is the experiments folder (default :func:`physim.nuclear.runs.home`) and ``home`` the nuclei home
    (default :func:`physim.nuclear.results.nuclei_home`). Returns {"nuclei": [...], "unassigned": [...]}.

    Each nucleus is a dict: "nuclide", "z", "a", "n"; "experiments" (each {"folder", "path", "name", "kind",
    "kind_label", "role", "runs", "analysed_runs", "rows", "latest"}: the role is "studied" (the nucleus the setup
    studies, :func:`physim.nuclear.results.studied_nuclide`), "partner" (the other nucleus of a Coulex
    experiment) or "results" (the experiment only produced rows for it); "latest" is the newest row of this
    nucleus from the experiment, as a :func:`cell` with its "label"); "runs" and "analysed_runs" (of the
    experiments that study it); "rows" (how many hold now) and "points" (the data points among them);
    "evaluation" ({"rows", "reference", "dataset_date", "taken"}); "table" (the headline rows of
    :func:`table`: :func:`view` has them all); and "open": what is still to do, each {"what", "text", ...} with
    "what" one of "route" (with "observable", "key", "label", "methods": routes physim has an experiment for and
    no row of one's own yet), "runs" (with "count": runs taken and not analysed) and "evaluation" (ENSDF's
    values not taken).

    "unassigned" lists the experiments that belong to no nucleus, each {"folder", "path", "name", "kind",
    "kind_label", "runs", "reason"}: elastic scattering, and a fusion-evaporation setup that does not mark the
    exit channel it studies.
    """
    experiments = _experiments(root)
    return {"nuclei": [_entry(name, experiments, home, full=False) for name in _names(experiments, home)],
            "unassigned": [{**{k: e[k] for k in ("folder", "path", "name", "kind", "kind_label", "runs")},
                            "reason": _unassigned(e)} for e in experiments if e["studied"] is None]}


def view(nuclide: str, root: Union[str, Path, None] = None, home: Union[str, Path, None] = None) -> dict:
    """One nucleus: the same dict as one entry of :func:`nuclei`, with every row of the :func:`table` (the
    evaluation's many level energies included; "headline" marks the rows that have more than the evaluation). A
    nucleus without experiments or rows comes back empty, ready for its first data point."""
    return _entry(results.nuclide_name(nuclide), _experiments(root), home, full=True)


def holding(nuclide: str, observable: Optional[str] = None, key=None, root: Union[str, Path, None] = None,
            home: Union[str, Path, None] = None) -> list:
    """The rows of a nucleus that hold now, each as a :func:`cell`: of one observable, and of one level or
    transition of it (``key``, as one says it), if given. The B(E2) and τ values the rows give each other are
    among them (:func:`with_converted`; their "via" says what each was made from). In the order of the table:
    by observable, then by level or transition, the independent evidence before the projected rows."""
    name = results.nuclide_name(nuclide)
    key = None if key is None else key_of(name, observable, key, check=False).text
    rows, via = with_converted(name, results.current(results.of_nuclide(name, root, home)))
    order, evidence = list(observables.OBSERVABLES), list(results.EVIDENCE)
    found = sorted(results.select(rows, observable, key),
                   key=lambda r: (order.index(r.observable), r.key, evidence.index(r.evidence), r.created))
    return [cell(r, via) for r in found]


def history(rows: list, observable: Optional[str] = None, key=None) -> list:
    """The rows of a nucleus (of one observable and level or transition, if given), every one that was ever
    written, the newest first. Each is a dict: "id", "created", "observable", "key", "evidence", "method",
    "method_label", "source", "detail", "value", "unit", "direction", "text" (the value as given), "status"
    ("holds", "superseded" by a newer row, "retracted", or "retraction": the row that withdrew another),
    "supersedes", "superseded_by" and "reason" (a retraction's, on both rows)."""
    later = {r.supersedes: r for r in rows if r.supersedes}
    out = []
    # The newest first; of rows written in the same second, the one written later (further down its file).
    for r in reversed(sorted(results.select(rows, observable, key), key=lambda r: r.created)):
        after = later.get(r.id)
        status = ("retraction" if r.retracted else "holds" if after is None
                  else "retracted" if after.retracted else "superseded")
        out.append({"id": r.id, "created": r.created, "observable": r.observable, "key": r.key,
                    "evidence": r.evidence, "method": r.method, "method_label": METHODS[r.method],
                    "source": source(r), "detail": r.detail, "value": r.value, "unit": r.unit,
                    "direction": r.direction,
                    "text": f"{value_text(r.value, r.unc('plus'), r.unc('minus'))} {r.unit}"
                            + observables.DIRECTIONS.get(r.direction, ""),
                    "status": status, "supersedes": r.supersedes, "superseded_by": after.id if after else None,
                    "reason": r.retracted or (after.retracted if after is not None else "")})
    return out


def observable_view(nuclide: str, observable: str, key, reference: Optional[str] = None,
                    root: Union[str, Path, None] = None, home: Union[str, Path, None] = None,
                    ledger: Optional[list] = None) -> dict:
    """One observable of one level or transition of a nucleus, by every route.

    ``key`` names the level or transition as one says it; ``reference`` is the id (or its first characters) of
    the row to compare against, the evaluated row when left out. ``ledger`` are the nucleus's rows (superseded
    ones among them) where the caller has them already; they are then not read from ``root`` and ``home`` (the
    report passes its own, which may hold a row that an analysis gives and the ledger does not hold yet). Nothing
    is written. Returns a dict:

    - "nuclide", "observable", "key", "label", "unit", "convention";
    - "comparison": the :class:`physim.nuclear.compare.Comparison` (the object: ``.to_dict()`` makes it plain
      data), or None when no row holds. Its ``entries`` are the independent rows, each against the reference;
      its ``projected`` the simulated runs, each against the value its simulation was given;
    - "cells": every row compared, as a :func:`cell`, by id (the source and what a converted row was made from);
    - "reference": the reference's cell, or None;
    - "routes": :func:`routes`, each with "applies" and "reason" for this nucleus (:func:`applies`);
    - "gives": :func:`gives` of the reference ([] without one);
    - "checks": :func:`checks` against other observables and the calculations ([] when there is none to make);
    - "used_for": what the registry says a value of the observable is used for, [{"text", "reference"}];
    - "history": :func:`history` of the rows of the ledgers (converted values are not in it: their sources are).
    """
    name = results.nuclide_name(nuclide)
    o = observables.get(observable)
    key = key_of(name, observable, key, check=False).text
    every = results.of_nuclide(name, root, home) if ledger is None else list(ledger)
    rows, via = with_converted(name, results.current(every))
    group = results.select(rows, observable, key)
    comparison = None
    if group:
        if reference is not None:
            fits = [r.id for r in group if r.id.startswith(str(reference).strip())]
            if len(fits) != 1:
                raise ValueError(f"'{reference}' names {len(fits) or 'no'} row{'s' if len(fits) > 1 else ''} of "
                                 f"{label(observable, key)} that hold now; the rows are "
                                 + ", ".join(f"{r.id} ({source(via.get(r.id, r))})" for r in group))
            reference = fits[0]
        comparison = _compare.compare(group, reference)
    ref = comparison.reference.row if comparison is not None and comparison.reference is not None else None
    return {"nuclide": name, "observable": observable, "key": key, "label": label(observable, key), "unit": o.unit,
            "convention": o.convention, "comparison": comparison, "cells": {r.id: cell(r, via) for r in group},
            "reference": None if ref is None else cell(ref, via),
            "routes": routes(name, observable, key, group, via, rows),
            "gives": gives(name, ref) if ref is not None else [],
            "checks": checks(name, observable, key, rows, via, ref), "used_for": [u.to_dict() for u in o.used_for],
            "history": history(every, observable, key)}


__all__ = ["DONE", "ENSDF", "IN_PHYSIM", "LATER", "NO_EXPERIMENT", "ONLY_PROJECTED", "POINT_EVIDENCE",
           "POINT_METHODS", "PointProblem", "add_data_point", "adopted_scheme", "applies", "cell", "checks",
           "correct_data_point", "decay", "evaluation", "gives", "history", "holding", "is_data_point",
           "is_evaluation", "key_of", "label", "level", "level_label", "nuclei", "observable_view",
           "partner_nuclide", "retract_data_point", "routes", "snapshot_evaluation", "source", "stability",
           "table", "transition", "value_text", "view", "with_converted"]
