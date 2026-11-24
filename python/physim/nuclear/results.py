"""The results ledger: one row per result, and one name for each level and transition.

A **row** (:class:`Row`) is one value of one observable of one nucleus by one method: what was found, with its
uncertainty, how, on what evidence and where it came from. Rows are appended to JSON-lines files and never
rewritten; a new analysis adds a row that names the one it supersedes, and a value that is withdrawn gets a row
that supersedes it and says why (``retracted``), so the files are the history::

    from physim.nuclear import results

    key = results.transition_key(results.level_key("62Ni", 1172.9, 2, +1), results.level_key("62Ni", 0.0, 0, +1))
    row = results.Row("62Ni", "B(E2)", key, 906.0, "e2fm4", "coulex-fit", "measured", direction="up", total=37.0,
                      origin=results.reference("Allmond et al.", 2014, doi="10.1103/PhysRevC.90.034309"))
    results.add_point(row)                      # <nuclei home>/62Ni/points.jsonl
    results.current(results.of_nuclide("62Ni"))  # every experiment's rows and the nucleus's own, superseded dropped

**Evidence** comes in two groups (:data:`EVIDENCE`). *Independent* (:data:`INDEPENDENT`): ``measured`` (a
published point or one's own data point: a value with its uncertainty), ``evaluated`` (ENSDF adopted), ``theory``
(shell model, rotor) and ``outside program`` (GOSIA and the like). *Projected*: a simulated run analysed back. It
shows the precision a setup would reach; it carries the value the simulation was given (``truth``) and is compared
with that and with how well the quantity is already known, never with a measurement as if it tested it
(:mod:`physim.nuclear.compare`).

**Where the rows live.** Rows an experiment produced are in ``results.jsonl`` in the experiment's folder, so the
folder stays complete on its own. Rows that belong to no experiment (literature and own data points, the
evaluation, shell-model results) are in ``<nuclei home>/<nuclide>/points.jsonl``, where the nuclei home is
``$PHYSIM_NUCLEI`` or the folder ``nuclei`` beside the experiments (:func:`nuclei_home`). :func:`of_nuclide` reads
both; there is no central file to go stale.

**Names.** A level's key (:func:`level_key`) is its nuclide, Jπ and ordinal ("62Ni:2+1" for 2⁺₁), anchored to the
level of the local ENSDF copy within 2 keV, so that the ordinal does not depend on how far a level scheme was
cut. A level ENSDF does not have, or any level without a local copy, keeps its energy and Jπ ("62Ni:2+@1173") and
is marked unanchored. A transition's key (:func:`transition_key`) is its upper level, its lower level and the
multipolarity ("62Ni:2+1>0+1:E2"); the direction of a B value is stated in the row, not in the key. The four
conventions in use elsewhere map onto these: a level scheme's index (:func:`key_of_level`), a name of
:func:`physim.nuclear.gamma.level_names` (:func:`key_of_name`), a shell-model state (2J, π, n)
(:func:`key_of_shell_state`) and the ``[reaction]`` energy (:func:`key_of_reaction`).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from dataclasses import dataclass, field, fields
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Optional, Union

from . import ensdf, observables, runs
from .levels import LevelScheme
from .names import SYMBOLS, parse_nuclide

#: The methods a value can come from, with their labels. The last four are for values from the literature only.
METHODS = {
    "coulex-yield": "Coulex from the yield (first order)",
    "coulex-fit": "Coulex fitted to all orders",
    "rdds": "RDDS (decay curve, DDCM)",
    "dsam": "DSAM (line shape)",
    "shell-model": "Shell model",
    "evaluation": "Evaluation",
    "converted": "Converted from another row",
    "nrf": "NRF",
    "inelastic-neutron": "(n,n′γ)",
    "electron-scattering": "Electron scattering",
    "other": "Other",
}
LITERATURE_ONLY = ("nrf", "inelastic-neutron", "electron-scattering", "other")

#: The kinds of evidence, with their labels, in the two groups: independent of physim's simulation, and projected.
EVIDENCE = {"measured": "Measured", "evaluated": "Evaluated", "theory": "Theory",
            "outside program": "Outside program", "projected": "Projected"}
INDEPENDENT = ("measured", "evaluated", "theory", "outside program")
PROJECTED = "projected"

#: Where a row can come from (``origin["kind"]``): a run of an experiment, a reference, a shell-model job, one's
#: own measurement, or the conversion of another row.
ORIGINS = ("run", "reference", "shell model", "own measurement", "conversion")

#: The file of an experiment's rows (in its folder) and of a nucleus's own rows (in its folder of the nuclei home).
RESULTS_FILE = "results.jsonl"
POINTS_FILE = "points.jsonl"
#: Environment variable naming the nuclei home.
ENV = "PHYSIM_NUCLEI"

#: How far a level may lie from the ENSDF level it is anchored to, keV.
ANCHOR_KEV = 2.0


def nuclide_name(nuclide: str) -> str:
    """A nuclide as the ledger writes it: "Ni-62" → "62Ni"."""
    z, a = parse_nuclide(nuclide)
    return f"{a}{SYMBOLS[z]}"


# ---------------------------------------------------------------------------------------------------------------
# Names: the key of a level and of a transition

_LEVEL = re.compile(r"(\?|\d+(?:/2)?)([+-])?(?:(\d+)|\((\d+)\)|@(\d+(?:\.\d+)?))")
_MULTIPOLE = re.compile(r"[EM][1-6]")


@dataclass(frozen=True)
class LevelKey:
    """The name of a level. ``anchored`` keys carry ENSDF's ordinal and energy; others the energy and Jπ they
    were given (or, for a calculated state without an ENSDF partner, the calculation's own ordinal)."""

    nuclide: str
    #: Twice the spin, and the parity (+1 or −1); None where not known.
    j2: Optional[int] = None
    parity: Optional[int] = None
    #: The place among the levels of this Jπ, counted from 1: ENSDF's when anchored.
    ordinal: Optional[int] = None
    energy_kev: Optional[float] = None
    anchored: bool = False

    @property
    def spin(self) -> Optional[float]:
        return None if self.j2 is None else self.j2 / 2.0

    @property
    def jpi(self) -> str:
        """ "2+", "3/2-", "2" (parity not known) or "?"."""
        if self.j2 is None:
            return "?"
        j = str(self.j2 // 2) if self.j2 % 2 == 0 else f"{self.j2}/2"
        return j + {1: "+", -1: "-", None: ""}[self.parity]

    @property
    def level(self) -> str:
        """The key without the nuclide: "2+1" (anchored), "2+(1)" (a calculation's own first 2⁺) or "2+@1173"
        (known by its energy, to the keV)."""
        if self.ordinal is None and self.energy_kev is None:
            raise ValueError(f"a level of {self.nuclide} is named by its ordinal or its energy; this key has neither")
        if self.anchored:
            return f"{self.jpi}{self.ordinal}"
        if self.ordinal is not None:
            return f"{self.jpi}({self.ordinal})"
        return f"{self.jpi}@{self.energy_kev:.0f}"

    @property
    def text(self) -> str:
        return f"{self.nuclide}:{self.level}"

    def __str__(self) -> str:
        return self.text

    @classmethod
    def parse(cls, text: str) -> LevelKey:
        """The key from its text. An anchored key read back has no energy: the text does not hold it."""
        nuclide, sep, level = str(text).partition(":")
        m = _LEVEL.fullmatch(level) if sep else None
        if m is None:
            raise ValueError(f"cannot read '{text}' as a level; a level is written like '62Ni:2+1' or "
                             f"'62Ni:2+@1173'")
        j, sign, ordinal, own, energy = m.groups()
        j2 = None if j == "?" else int(j[:-2]) if j.endswith("/2") else 2 * int(j)
        parity = {"+": 1, "-": -1, None: None}[sign]
        return cls(nuclide_name(nuclide), j2, parity, int(ordinal or own) if (ordinal or own) else None,
                   float(energy) if energy else None, anchored=ordinal is not None)


@dataclass(frozen=True)
class TransitionKey:
    """The name of a transition: the upper level (the one that decays), the lower level and the multipolarity."""

    initial: LevelKey
    final: LevelKey
    multipolarity: str

    @property
    def nuclide(self) -> str:
        return self.initial.nuclide

    @property
    def anchored(self) -> bool:
        return self.initial.anchored and self.final.anchored

    @property
    def text(self) -> str:
        return f"{self.nuclide}:{self.initial.level}>{self.final.level}:{self.multipolarity}"

    def __str__(self) -> str:
        return self.text

    @classmethod
    def parse(cls, text: str) -> TransitionKey:
        parts = str(text).split(":")
        levels = parts[1].split(">") if len(parts) == 3 else []
        if len(levels) != 2 or not _MULTIPOLE.fullmatch(parts[2]):
            raise ValueError(f"cannot read '{text}' as a transition; a transition is written like "
                             f"'62Ni:2+1>0+1:E2'")
        return cls(LevelKey.parse(f"{parts[0]}:{levels[0]}"), LevelKey.parse(f"{parts[0]}:{levels[1]}"), parts[2])


@lru_cache(maxsize=64)
def _ensdf_levels(nuclide: str, where: str) -> tuple:
    """The adopted levels of a nuclide as ((energy keV, 2J, parity, ordinal), ...); 2J, parity and the ordinal
    are None where ENSDF names no single Jπ. Empty without a local copy or without the nuclide. A tentative
    assignment, "(2+)", counts as that Jπ. ``where`` is the ENSDF folder, so a change of it is noticed."""
    try:
        scheme = LevelScheme.from_adopted(ensdf.adopted(nuclide))
    except (ensdf.EnsdfMissing, ValueError):
        return ()
    seen: dict = {}
    out = []
    for lv in sorted((x for x in scheme.levels if x.energy is not None), key=lambda x: x.energy.value):
        if lv.spin is None or lv.parity is None:
            out.append((float(lv.energy.value), None, None, None))
            continue
        jpi = (int(round(2 * lv.spin)), int(lv.parity))
        seen[jpi] = seen.get(jpi, 0) + 1
        out.append((float(lv.energy.value), jpi[0], jpi[1], seen[jpi]))
    return tuple(out)


def ensdf_levels(nuclide: str) -> tuple:
    """((energy keV, 2J, parity, ordinal), ...) of the local ENSDF copy's adopted levels; empty without one."""
    if not ensdf.available():
        return ()
    return _ensdf_levels(nuclide_name(nuclide), str(ensdf.folder()))


def level_key(nuclide: str, energy_kev: float, spin: Optional[float] = None, parity: Optional[int] = None,
              tolerance_kev: float = ANCHOR_KEV) -> LevelKey:
    """The key of the level of ``nuclide`` at ``energy_kev`` with spin ``spin`` and parity ``parity`` (+1, −1).

    With a local ENSDF copy the level is anchored to the nearest adopted level within ``tolerance_kev`` that has
    the same Jπ, and takes its ordinal and energy from there. Given without a Jπ, it is anchored to the nearest
    adopted level and takes its Jπ too. A level ENSDF does not have (or names with another or no single Jπ), and
    any level without a local copy, keeps the energy and Jπ given and is unanchored.
    """
    name = nuclide_name(nuclide)
    j2 = None if spin is None else int(round(2 * float(spin)))
    near = sorted(((abs(e - energy_kev), e, j, p, n) for e, j, p, n in ensdf_levels(name)
                   if abs(e - energy_kev) <= tolerance_kev), key=lambda found: found[0])
    for _, e, j, p, n in near:
        if j2 is None:
            if n is None:  # the nearest level has no single Jπ: no ordinal to anchor to
                break
            return LevelKey(name, j, p, n, e, True)
        if n is not None and j == j2 and (parity is None or p == parity):
            return LevelKey(name, j, p, n, e, True)
    return LevelKey(name, j2, parity, None, float(energy_kev), False)


def transition_key(initial: LevelKey, final: LevelKey, multipolarity: str = "E2") -> TransitionKey:
    """The key of the transition between two levels of one nuclide. The upper level comes first whichever way
    they are given, where both energies are known (give the upper one first otherwise)."""
    if initial.nuclide != final.nuclide:
        raise ValueError(f"a transition joins two levels of one nuclide, not {initial.nuclide} and {final.nuclide}")
    if not _MULTIPOLE.fullmatch(str(multipolarity)):
        raise ValueError(f"write the multipolarity as E1, E2, M1, ...; got '{multipolarity}'")
    if initial.energy_kev is not None and final.energy_kev is not None and initial.energy_kev < final.energy_kev:
        initial, final = final, initial
    return TransitionKey(initial, final, multipolarity)


def key_of_level(scheme, index: int) -> LevelKey:
    """The key of level ``index`` of a :class:`~physim.nuclear.levels.LevelScheme` (the scheme's own numbering,
    which changes when the scheme is cut)."""
    lv = scheme.levels[index]
    return level_key(scheme.nuclide, lv.energy.value, lv.spin, lv.parity)


def key_of_name(scheme, name: str) -> LevelKey:
    """The key of the level a scheme's :func:`~physim.nuclear.gamma.level_names` calls ``name`` ("2+₁", "4+")."""
    from .gamma import level_names

    names = level_names(scheme.levels)
    if name not in names:
        raise ValueError(f"the level scheme of {scheme.nuclide} has no level '{name}'; it has {', '.join(names)}")
    return key_of_level(scheme, names.index(name))


def key_of_transition(scheme, initial: int, final: int, multipolarity: str = "E2") -> TransitionKey:
    """The key of the transition between two levels of a scheme, by their indices."""
    return transition_key(key_of_level(scheme, initial), key_of_level(scheme, final), multipolarity)


def key_of_shell_state(nuclide: str, state: tuple, energy_kev: Optional[float] = None) -> LevelKey:
    """The key of a shell-model state (2J, parity, index), as :mod:`physim.nuclear.shellmodel` labels it: the
    n-th calculated state of a Jπ is the n-th adopted level of that Jπ. Where ENSDF has no such level (or there
    is no local copy) the key keeps the calculation's own ordinal ("62Ni:2+(1)") and is unanchored."""
    name = nuclide_name(nuclide)
    j2, parity, index = (int(x) for x in state)
    for e, j, p, n in ensdf_levels(name):
        if (j, p, n) == (j2, parity, index):
            return LevelKey(name, j, p, n, e, True)
    return LevelKey(name, j2, parity, index, energy_kev, False)


def reaction_nuclide(experiment) -> Optional[str]:
    """The nucleus a Coulomb-excitation setup's ``[reaction]`` excites; None for another kind of reaction."""
    if getattr(experiment, "excitation", None) is None:
        return None
    from .analysis import ex_nuclide

    return nuclide_name(ex_nuclide(experiment))


def key_of_reaction(experiment) -> TransitionKey:
    """The key of the transition a Coulomb-excitation setup's ``[reaction]`` names by its energy and
    multipolarity: the state at that energy, reached from the ground state. The state's Jπ is what an Eλ
    excitation of a 0⁺ ground state gives (2⁺ for E2), as the ``[reaction]`` assumes; where ENSDF's ground state
    is not 0⁺ the state is found by its energy alone."""
    nuclide = reaction_nuclide(experiment)
    if nuclide is None:
        raise ValueError("the setup's [reaction] names no excited state (it is not Coulomb excitation)")
    exc = experiment.excitation
    lam = int(exc.multipolarity[1])
    ground = level_key(nuclide, 0.0)
    if not ground.anchored:
        ground = level_key(nuclide, 0.0, 0.0, +1)
    if (ground.j2, ground.parity) == (0, 1):
        state = level_key(nuclide, 1e3 * exc.energy_mev, float(lam), (-1) ** lam)
    else:
        state = level_key(nuclide, 1e3 * exc.energy_mev)
    return transition_key(state, ground, exc.multipolarity)


# ---------------------------------------------------------------------------------------------------------------
# A row


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@lru_cache(maxsize=1)
def _code() -> tuple:
    from .. import __version__
    from .report import _commit

    return __version__, _commit()


def code() -> dict:
    """The physim that is running: {"version", "commit"} (the commit is "unknown" outside a git checkout)."""
    version, commit = _code()
    return {"version": version, "commit": commit}


def from_run(experiment: str, run: str, settings: Optional[dict] = None, **more) -> dict:
    """The origin of a row an experiment's run produced: the experiment, the run and the analysis settings."""
    return {"kind": "run", "experiment": experiment, "run": run, "settings": settings or {}, **more}


def reference(authors: str, year: int, doi: str = "", where: str = "", **more) -> dict:
    """The origin of a row from the literature: authors, year, DOI and where in the paper the value stands."""
    return {"kind": "reference", "authors": authors, "year": year, "doi": doi, "where": where, **more}


def shell_model_job(interaction: str, space: str = "", effective_charges: Optional[dict] = None, **more) -> dict:
    """The origin of a shell-model row: the interaction, the valence space and the effective charges."""
    return {"kind": "shell model", "interaction": interaction, "space": space,
            "effective_charges": effective_charges or {}, **more}


def own_measurement(date: str, note: str = "", **more) -> dict:
    """The origin of one's own data point: its date and a note."""
    return {"kind": "own measurement", "date": date, "note": note, **more}


def _plain(x):
    """Plain JSON data (tuples become lists), so that a row read back equals the row written."""
    return json.loads(json.dumps(x, ensure_ascii=False, default=str))


def _amount(name: str, x) -> Optional[float]:
    if x is None:
        return None
    x = float(x)
    if not math.isfinite(x) or x < 0:
        raise ValueError(f"the uncertainty '{name}' must be a number that is not negative, got {x!r}")
    return x


@dataclass(frozen=True)
class Row:
    """One value of one observable of one nucleus by one method.

    The uncertainty is absolute, in the row's unit, and comes in one of three forms: ``stat`` and ``sys`` (either
    or both), one ``total``, or an asymmetric ``plus`` and ``minus``; a row with none (a calculation) is allowed.
    ``budget`` holds the terms of the systematic uncertainty by name, absolute too. Everything is checked on
    construction and a ``ValueError`` says what is wrong. ``id`` is a hash of the content; ``created`` and
    ``code`` are filled in when left out. A row is never changed: :meth:`revised` makes the row that supersedes
    it.
    """

    nuclide: str
    #: A key of :data:`physim.nuclear.observables.OBSERVABLES`.
    observable: str
    #: The level or transition (the text of a :class:`LevelKey` or :class:`TransitionKey`, or the key itself).
    key: str
    value: float
    unit: str
    #: A key of :data:`METHODS`, and a free line of detail beside it.
    method: str
    #: A key of :data:`EVIDENCE`.
    evidence: str
    #: Where the row comes from: :func:`from_run`, :func:`reference`, :func:`shell_model_job`,
    #: :func:`own_measurement`, or the conversion of another row. Its "kind" is one of :data:`ORIGINS`.
    origin: dict = field(default_factory=dict)
    #: "up" (B↑) or "down" (B↓): stated for a B value, never implied; None for an observable without one.
    direction: Optional[str] = None
    stat: Optional[float] = None
    sys: Optional[float] = None
    total: Optional[float] = None
    plus: Optional[float] = None
    minus: Optional[float] = None
    budget: Optional[dict] = None
    detail: str = ""
    #: Projected rows only: the value the simulation was given, in the row's unit and direction.
    truth: Optional[float] = None
    #: The ids of rows this one is not independent of (an evaluation and the measurements it adopts; a value
    #: converted from another row).
    depends_on: tuple = ()
    #: The id of the earlier row that this one, a new analysis of the same thing, replaces.
    supersedes: Optional[str] = None
    #: Why the row it supersedes was withdrawn. A retraction keeps that row's content and holds no value of its
    #: own: :func:`current` leaves it out, the history shows it.
    retracted: str = ""
    created: str = ""
    code: Optional[dict] = None
    id: str = ""

    def __post_init__(self):
        def put(name: str, value) -> None:
            object.__setattr__(self, name, value)

        put("nuclide", nuclide_name(self.nuclide))
        o = observables.get(self.observable)
        if self.method not in METHODS:
            raise ValueError(f"'{self.method}' is not a method; the methods are {', '.join(METHODS)}")
        if self.evidence not in EVIDENCE:
            raise ValueError(f"'{self.evidence}' is not a kind of evidence; the kinds are {', '.join(EVIDENCE)}")
        put("value", float(self.value))
        if not math.isfinite(self.value):
            raise ValueError(f"the value must be a number, got {self.value!r}")
        # The unit and the direction: a B says which way; Weisskopf units are always ↓.
        if self.unit not in o.units:
            raise ValueError(f"'{self.unit}' is not a unit of {o.label}; use {', '.join(o.units)}")
        if o.direction is None:
            if self.direction is not None:
                raise ValueError(f"{o.label} has no direction")
        elif self.direction not in observables.DIRECTIONS:
            raise ValueError(f"{o.label} needs its direction stated: 'up' (↑) or 'down' (↓)")
        elif self.unit == observables.WEISSKOPF and self.direction != observables.DOWN:
            raise ValueError(f"{o.label} in Weisskopf units is always the ↓ value")
        # The key names a level or a transition of this nuclide, as the observable needs.
        kind = TransitionKey if o.belongs_to == "transition" else LevelKey
        key = self.key if isinstance(self.key, kind) else kind.parse(str(self.key))
        if key.nuclide != self.nuclide:
            raise ValueError(f"the key '{key}' is not of {self.nuclide}")
        put("key", key.text)
        # The uncertainty, in one form.
        for name in ("stat", "sys", "total", "plus", "minus"):
            put(name, _amount(name, getattr(self, name)))
        forms = [self.stat is not None or self.sys is not None, self.total is not None,
                 self.plus is not None or self.minus is not None]
        if sum(forms) > 1:
            raise ValueError("give the uncertainty in one form: stat and sys, or one total, or plus and minus")
        if (self.plus is None) != (self.minus is None):
            raise ValueError("an asymmetric uncertainty needs both plus and minus")
        if self.budget is not None:
            put("budget", {str(k): _amount(str(k), v) for k, v in dict(self.budget).items()} or None)
        # Projected rows carry the value the simulation was given; no other row does.
        if self.evidence == PROJECTED and self.truth is None:
            raise ValueError("a projected row needs 'truth', the value the simulation was given")
        if self.evidence != PROJECTED and self.truth is not None:
            raise ValueError("'truth' is for projected rows only: an independent value has none")
        if self.truth is not None:
            put("truth", float(self.truth))
        origin = _plain(dict(self.origin or {}))
        if origin.get("kind") not in ORIGINS:
            raise ValueError(f"the origin needs a 'kind', one of {', '.join(ORIGINS)}; got {origin.get('kind')!r}")
        put("origin", origin)
        put("depends_on", tuple(str(x) for x in (self.depends_on or ())))
        put("supersedes", str(self.supersedes) if self.supersedes else None)
        put("retracted", str(self.retracted or "").strip())
        if self.retracted and not self.supersedes:
            raise ValueError("a retraction names the row it withdraws ('supersedes')")
        put("detail", str(self.detail or ""))
        put("created", self.created or _now())
        put("code", _plain(dict(self.code)) if self.code is not None else code())
        if not self.id:
            put("id", _digest(self.to_dict()))

    # -- what the row says ------------------------------------------------------------------------------------------

    @property
    def independent(self) -> bool:
        """Whether the evidence is independent of physim's simulation (anything but a projected row)."""
        return self.evidence in INDEPENDENT

    @property
    def asymmetric(self) -> bool:
        return self.plus is not None

    def unc(self, side: Optional[str] = None) -> Optional[float]:
        """The total absolute uncertainty, or None for a row without one: stat and sys in quadrature, or the one
        total; for an asymmetric uncertainty the ``side`` asked for ("plus" or "minus"), or the mean of the two."""
        if self.plus is not None:
            return {"plus": self.plus, "minus": self.minus}.get(side, (self.plus + self.minus) / 2.0)
        if self.total is not None:
            return self.total
        if self.stat is None and self.sys is None:
            return None
        return math.hypot(self.stat or 0.0, self.sys or 0.0)

    @property
    def relative(self) -> Optional[float]:
        """The total uncertainty over the size of the value; None without an uncertainty or for a value of 0."""
        u = self.unc()
        return None if u is None or self.value == 0 else u / abs(self.value)

    def revised(self, **changes) -> Row:
        """The row a new analysis puts in this one's place: the same row with ``changes``, made now by the
        running code, and superseding this one."""
        d = self.to_dict()
        d.update(id="", created="", code=None, supersedes=self.id)
        d.update(changes)
        return Row.from_dict(d)

    # -- plain data -------------------------------------------------------------------------------------------------

    def to_dict(self) -> dict:
        """The row as plain data, one line of a ledger file. Fields that are not set are left out."""
        d = {"id": self.id, "created": self.created, "code": self.code}
        for f in fields(self):
            value = getattr(self, f.name)
            if f.name in d or value is None or value == "" or value == ():
                continue
            d[f.name] = list(value) if isinstance(value, tuple) else value
        return d

    @classmethod
    def from_dict(cls, d: dict) -> Row:
        """The row from plain data. Fields this version does not know are ignored; the id is kept as written."""
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in known})


def _digest(d: dict) -> str:
    """The content hash of a row: of everything but the id itself."""
    content = {k: v for k, v in d.items() if k != "id"}
    text = json.dumps(content, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


# ---------------------------------------------------------------------------------------------------------------
# The files: append-only JSON lines


def nuclei_home() -> Path:
    """Where the rows that belong to no experiment live: ``$PHYSIM_NUCLEI``, or the folder ``nuclei`` beside the
    experiments (``~/.physim/nuclei``)."""
    return Path(os.environ.get(ENV) or runs.home().parent / "nuclei").expanduser()


def nucleus_folder(nuclide: str, home: Union[str, Path, None] = None) -> Path:
    """The folder of one nucleus in the nuclei home."""
    return (Path(home).expanduser() if home is not None else nuclei_home()) / nuclide_name(nuclide)


def read(path: Union[str, Path], problems: Optional[list] = None) -> list:
    """The rows of a ledger file, in the order written; none if there is no file. A line that cannot be read (a
    row of a later version, a damaged line) is skipped, and said in ``problems`` if a list is given."""
    path = Path(path)
    out = []
    if not path.is_file():
        return out
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            out.append(Row.from_dict(json.loads(line)))
        except (ValueError, TypeError) as e:
            if problems is not None:
                problems.append(f"{path.name}, line {n}: {e}")
    return out


def _append(path: Path, row: Row) -> bool:
    """Add a row at the end of a ledger file; nothing already there is touched. False (and nothing written) when
    the file already has the row."""
    path.parent.mkdir(parents=True, exist_ok=True)
    before = path.read_text(encoding="utf-8") if path.is_file() else ""
    for line in before.splitlines():
        try:
            if line.strip() and json.loads(line).get("id") == row.id:
                return False
        except (ValueError, AttributeError):
            continue
    with path.open("a", encoding="utf-8", newline="\n") as f:
        f.write(("" if not before or before.endswith("\n") else "\n")
                + json.dumps(row.to_dict(), ensure_ascii=False) + "\n")
    return True


def _append_all(path: Path, rows: list) -> list:
    """Add rows at the end of a ledger file, reading it once. Returns the rows written: those the file did not
    hold already."""
    path.parent.mkdir(parents=True, exist_ok=True)
    before = path.read_text(encoding="utf-8") if path.is_file() else ""
    have = set()
    for line in before.splitlines():
        try:
            have.add(json.loads(line).get("id") if line.strip() else None)
        except (ValueError, AttributeError):
            continue
    new = []
    for row in rows:
        if row.id not in have:
            have.add(row.id)
            new.append(row)
    if new:
        with path.open("a", encoding="utf-8", newline="\n") as f:
            f.write(("" if not before or before.endswith("\n") else "\n")
                    + "".join(json.dumps(row.to_dict(), ensure_ascii=False) + "\n" for row in new))
    return new


def append(folder: Union[str, Path], row: Row) -> bool:
    """Add a row to an experiment's ledger (``results.jsonl`` in its folder). False if it is already there."""
    return _append(Path(folder).expanduser() / RESULTS_FILE, row)


def add_point(row: Row, home: Union[str, Path, None] = None) -> bool:
    """Add a row that belongs to no experiment to its nucleus's ledger (``points.jsonl`` in the nucleus's
    folder of the nuclei home). False if it is already there."""
    return _append(nucleus_folder(row.nuclide, home) / POINTS_FILE, row)


def add_points(rows: list, home: Union[str, Path, None] = None) -> list:
    """Add several rows of one nucleus that belong to no experiment to its ledger, as :func:`add_point` does one,
    reading the file once. Returns the rows written."""
    rows = list(rows)
    if len({r.nuclide for r in rows}) > 1:
        raise ValueError("the rows are of several nuclei: " + ", ".join(sorted({r.nuclide for r in rows})))
    return _append_all(nucleus_folder(rows[0].nuclide, home) / POINTS_FILE, rows) if rows else []


def of_experiment(folder: Union[str, Path]) -> list:
    """Every row an experiment produced, superseded ones included (:func:`current` drops them)."""
    return read(Path(folder).expanduser() / RESULTS_FILE)


def _experiment_folders(root: Union[str, Path, None]) -> list:
    root = Path(root).expanduser() if root is not None else runs.home()
    if not root.is_dir():
        return []
    return sorted(p for p in root.iterdir() if (p / "experiment.toml").is_file())


def of_nuclide(nuclide: str, root: Union[str, Path, None] = None, home: Union[str, Path, None] = None) -> list:
    """Every row of a nucleus: from each experiment's folder under ``root`` (default
    :func:`physim.nuclear.runs.home`) and from the nucleus's own folder in the nuclei home, in the order
    created. Superseded rows are included: this is the history, and :func:`current` of it what holds now."""
    name = nuclide_name(nuclide)
    rows = [r for folder in _experiment_folders(root) for r in of_experiment(folder) if r.nuclide == name]
    rows += read(nucleus_folder(name, home) / POINTS_FILE)
    return sorted(rows, key=lambda r: r.created)


def studied_nuclide(experiment) -> Optional[str]:
    """The nucleus a setup studies: the excited one of a Coulomb-excitation setup, and the residue of the exit
    channel a fusion-evaporation setup marks as studied (``[run] channel``). None for a fusion-evaporation setup
    without that mark (its largest channel only sets the beam time: it is not a choice of nucleus) and for
    elastic scattering."""
    name = reaction_nuclide(experiment)
    if name is not None:
        return name
    fusion, run = getattr(experiment, "fusion", None), getattr(experiment, "run", None)
    if fusion is None or not fusion.channels or run is None or run.channel is None:
        return None
    from .experiment import studied_channel_index

    i = studied_channel_index(fusion.channels, run.channel)
    try:
        return None if i is None else nuclide_name(fusion.channels[i].residue_name)
    except ValueError:  # a residue the mass table does not know
        return None


def experiment_nuclide(folder: Union[str, Path]) -> Optional[str]:
    """The nucleus an experiment studies, from its current setup (:func:`studied_nuclide`). None for a setup
    that studies none and for one that cannot be read."""
    from .experiment import Experiment

    try:
        return studied_nuclide(Experiment.load(Path(folder).expanduser() / "experiment.toml"))
    except (ValueError, OSError, KeyError, IndexError):
        return None


def nuclei(root: Union[str, Path, None] = None, home: Union[str, Path, None] = None) -> list:
    """The nuclei that have rows or experiments, lightest first: dicts with "nuclide", "rows" (how many hold now,
    superseded ones dropped) and "experiments" (the names of the folders of the experiments that study it or
    produced rows for it)."""
    found: dict = {}

    def entry(name: str) -> dict:
        return found.setdefault(name, {"nuclide": name, "rows": [], "experiments": []})

    for folder in _experiment_folders(root):
        rows = of_experiment(folder)
        studied = experiment_nuclide(folder)
        for name in sorted({r.nuclide for r in rows} | ({studied} if studied else set())):
            entry(name)["experiments"].append(folder.name)
        for r in rows:
            entry(r.nuclide)["rows"].append(r)
    base = Path(home).expanduser() if home is not None else nuclei_home()
    for p in sorted(base.iterdir()) if base.is_dir() else []:
        for r in read(p / POINTS_FILE):
            entry(r.nuclide)["rows"].append(r)
    out = [{**e, "rows": len(current(e["rows"]))} for e in found.values()]
    return sorted(out, key=lambda e: parse_nuclide(e["nuclide"]))


def current(rows: list) -> list:
    """The rows that hold now: those no other row of the list supersedes, without the retractions (a retraction
    withdraws the row it supersedes and holds no value itself)."""
    replaced = {r.supersedes for r in rows if r.supersedes}
    return [r for r in rows if r.id not in replaced and not r.retracted]


def history(rows: list, row: Row) -> list:
    """A row and the rows it replaced, the newest first: the chain of ``supersedes`` as far as the list holds it."""
    by_id = {r.id: r for r in rows}
    chain = [row]
    while chain[-1].supersedes in by_id and by_id[chain[-1].supersedes] not in chain:
        chain.append(by_id[chain[-1].supersedes])
    return chain


def select(rows: list, observable: Optional[str] = None, key=None, evidence: Optional[str] = None) -> list:
    """The rows of one observable, one level or transition (a key or its text) and one kind of evidence; each
    filter is optional. Keys match by their text, so unanchored keys of one level match only when they were
    written alike."""
    text = None if key is None else str(key)
    return [r for r in rows if (observable is None or r.observable == observable)
            and (text is None or r.key == text) and (evidence is None or r.evidence == evidence)]


# ---------------------------------------------------------------------------------------------------------------
# Backfill: rows from the analyses that runs already hold


def backfill(folder: Union[str, Path]) -> list:
    """Make projected rows from the ``analysis.json`` files an experiment's runs already hold, and append them to
    the experiment's ledger. Returns the rows added; calling it again adds none.

    Each kept Coulomb-excitation analysis (:class:`physim.nuclear.analysis.Result`) gives one row: B(E2)↑ of the
    transition the run's ``[reaction]`` names, with the statistical and systematic uncertainties and the budget
    made absolute, and the value the setup held as ``truth``. The row is dated by when the analysis was saved;
    the code that made it was not recorded and is left unknown. A run analysed again since the last backfill
    gives a row that supersedes the earlier one. Runs of other reactions, runs without an analysis and analyses
    that cannot be read are passed over.
    """
    from .analysis import RESULT_FILE
    from .experiment import Experiment

    folder = Path(folder).expanduser()
    have = of_experiment(folder)
    added = []
    for run in runs.ExperimentFolder(folder).run_folders():
        try:
            d = json.loads((run / RESULT_FILE).read_text(encoding="utf-8"))
            experiment = Experiment.load(run / "setup.toml")
            if experiment.excitation is None or experiment.excitation.multipolarity != "E2":
                continue
            key = key_of_reaction(experiment)
            value = float(d["b_e2fm4"])
            saved = d.get("saved") or datetime.fromtimestamp((run / RESULT_FILE).stat().st_mtime,
                                                              timezone.utc).isoformat(timespec="seconds")
            row = Row(
                key.nuclide, "B(E2)", key, value, "e2fm4", "coulex-yield", PROJECTED, direction=observables.UP,
                stat=abs(value) * float(d["statistical"]), sys=abs(value) * float(d["systematic"]),
                budget={k: abs(value) * float(v) for k, v in (d.get("budget") or {}).items()},
                detail=", ".join(x for x in (f"{d.get('normalisation')} normalisation" if d.get("normalisation")
                                             else "", f"{d.get('gate')} gate" if d.get("gate") else "") if x),
                truth=float(d["truth_e2fm4"]), created=saved, code={"version": None, "commit": None},
                origin=from_run(folder.name, run.name, d.get("settings"), analysed=saved, backfilled=True))
        except (OSError, ValueError, TypeError, KeyError):
            continue
        earlier = [r for r in have if r.origin.get("run") == run.name and r.observable == row.observable
                   and r.key == row.key]
        if any(r.origin.get("analysed") == saved for r in earlier):
            continue
        earlier = current(earlier)
        if earlier:
            row = Row.from_dict({**row.to_dict(), "id": "", "supersedes": earlier[-1].id})
        if append(folder, row):
            have.append(row)
            added.append(row)
    return added


__all__ = ["ANCHOR_KEV", "ENV", "EVIDENCE", "INDEPENDENT", "LITERATURE_ONLY", "METHODS", "ORIGINS", "POINTS_FILE",
           "PROJECTED", "RESULTS_FILE", "LevelKey", "Row", "TransitionKey", "add_point", "add_points", "append",
           "backfill",
           "code", "current", "ensdf_levels", "experiment_nuclide", "from_run", "history", "key_of_level",
           "key_of_name", "key_of_reaction", "key_of_shell_state", "key_of_transition", "level_key", "nuclei",
           "nuclei_home", "nucleus_folder", "nuclide_name", "of_experiment", "of_nuclide", "own_measurement",
           "reaction_nuclide", "read", "reference", "select", "shell_model_job", "studied_nuclide",
           "transition_key"]
