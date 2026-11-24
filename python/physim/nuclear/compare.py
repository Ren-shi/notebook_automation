"""The comparison engine: the rows of one observable of one transition or level, set beside each other.

Given rows of the results ledger (:mod:`physim.nuclear.results`), :func:`compare` converts them to the
observable's canonical convention (:mod:`physim.nuclear.observables`), takes a reference (the evaluated row
unless another is named) and gives, by a rule fixed beforehand, how each row stands to it::

    from physim.nuclear import compare, results

    rows = results.select(results.of_nuclide("62Ni"), "B(E2)", key)
    c = compare.compare(rows)
    [(e.label, e.difference_percent, e.pull, e.verdict) for e in c.entries]
    [(p.closure, p.reach_percent, p.known_percent, p.best_measured_percent) for p in c.projected]

**Independent rows** (measured, evaluated, theory, outside program) each get the difference from the reference
in percent, the ratio to it, and the pull: the difference over the two uncertainties in quadrature, where an
asymmetric uncertainty counts with the side that faces the other value. The verdict follows the pull: below
:data:`AGREE_BELOW` (2) the two agree, up to :data:`DIFFER_ABOVE` (3) there is tension, above it they differ. A
row without an uncertainty gets its ratio and no pull or verdict. A row that is not independent of the reference
(the two share a source through ``depends_on``) is flagged.

**Projected rows** (a simulated run analysed back) are handled apart and never pulled against the reference.
Each is closed against its ``truth``, the value the simulation was given, over its own total uncertainty; and what
it says about the setup is put beside what is known: the precision it reaches, the precision of the reference, and
that of the best single measurement.

The result is plain data (:class:`Comparison`, :meth:`Comparison.to_dict`) that a table, a plot and a report can
be built from. :func:`convert` makes a row of another observable from a row (a B(E2) from a τ, through the γ-ray
energy, the branching and the conversion), so that it can join the comparison.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, fields, is_dataclass
from typing import Optional, Union

from . import observables
from .results import METHODS, PROJECTED, LevelKey, Row, TransitionKey, current

#: The rule, fixed before any comparison is made: |pull| below 2 agree, 2 to 3 tension, above 3 differ.
AGREE_BELOW = 2.0
DIFFER_ABOVE = 3.0
VERDICTS = {"agree": "Agree", "tension": "Tension", "differ": "Differ"}


def verdict(pull: Optional[float]) -> Optional[str]:
    """ "agree", "tension" or "differ" for a pull, by the rule; None without a pull."""
    if pull is None or not math.isfinite(pull):
        return None
    size = abs(pull)
    return "agree" if size < AGREE_BELOW else "tension" if size <= DIFFER_ABOVE else "differ"


# ---------------------------------------------------------------------------------------------------------------
# The canonical convention


@dataclass(frozen=True)
class Canonical:
    """A row's value and uncertainty in its observable's canonical unit and direction."""

    row: Row
    value: float
    unit: str
    direction: Optional[str]
    #: The upper and lower uncertainty (equal for a symmetric one), or None for a row without one.
    plus: Optional[float]
    minus: Optional[float]
    stat: Optional[float]
    sys: Optional[float]
    truth: Optional[float]
    #: The steps of the conversion in words; empty when the row was in the canonical convention already.
    how: str

    def unc(self, facing: Optional[float] = None) -> Optional[float]:
        """The total uncertainty: towards the value ``facing`` (the upper one if that lies above), or the mean
        of the two sides."""
        if self.plus is None:
            return None
        if facing is None:
            return (self.plus + self.minus) / 2.0
        return self.plus if facing > self.value else self.minus

    @property
    def percent(self) -> Optional[float]:
        """The total uncertainty in percent of the value."""
        u = self.unc()
        return None if u is None or self.value == 0 else 100.0 * u / abs(self.value)


def canonical(row: Row) -> Canonical:
    """``row`` in its observable's canonical convention: B(E2)↑ in e² fm⁴, τ in ps, Q_s in e fm², E in keV. The
    spins that B↓ → B↑ needs come from the row's key, the mass number for Weisskopf units from its nuclide."""
    o = observables.get(row.observable)
    spins = {}
    if o.belongs_to == "transition":
        key = TransitionKey.parse(row.key)
        spins = {"j_upper": key.initial.spin, "j_lower": key.final.spin}
    try:
        factor, how = observables.canonical_factor(row.observable, row.unit, row.direction, nuclide=row.nuclide,
                                                   **spins)
    except ValueError as e:
        raise ValueError(f"row {row.id} ({row.key}): {e}") from None

    def scaled(x: Optional[float]) -> Optional[float]:
        return None if x is None else x * abs(factor)

    plus = row.unc("plus") if factor > 0 else row.unc("minus")
    minus = row.unc("minus") if factor > 0 else row.unc("plus")
    return Canonical(row, row.value * factor, o.unit, o.direction, scaled(plus), scaled(minus), scaled(row.stat),
                     scaled(row.sys), None if row.truth is None else row.truth * factor, how)


def convert(row: Row, to: str, key: Union[TransitionKey, LevelKey, str, None] = None, *, energy_kev: float,
            branching: float = 1.0, conversion: float = 0.0, mixing: Optional[float] = None) -> Row:
    """A row of the observable ``to`` made from ``row``: B(E2) from a level's τ, or τ from a B(E2).

    ``energy_kev`` is the γ-ray energy of the transition, ``branching`` its share of all the level's decays,
    ``conversion`` its conversion coefficient α and ``mixing`` the E2/M1 mixing ratio δ, as
    :func:`physim.nuclear.lifetimes.b_from_tau` takes them, which does the work. ``key`` is the transition the
    B(E2) belongs to (its upper level must be the τ row's level); for τ from B(E2) the level is the transition's
    upper one.

    The new row is in the canonical convention, has the method "converted", the evidence of ``row`` and says
    what it was made from; it depends on ``row`` and on everything ``row`` depends on. B and τ are inversely
    proportional, so the relative uncertainties carry over (to first order) and an asymmetric one swaps sides.
    """
    c = canonical(row)
    inputs = {"energy_kev": float(energy_kev), "branching": float(branching), "conversion": float(conversion),
              "mixing": mixing}
    if (row.observable, to) == ("tau", "B(E2)"):
        if key is None:
            raise ValueError("a B(E2) from a τ needs the key of the transition it belongs to")
        new_key = key if isinstance(key, TransitionKey) else TransitionKey.parse(str(key))
        if new_key.initial.text != row.key:
            raise ValueError(f"the transition {new_key} does not start at the level {row.key} whose τ this is")
        b_down = observables.b_from_tau(c.value, energy_kev, branching, conversion, mixing,
                                        new_key.multipolarity)
        ratio = observables.b_up_over_down(new_key.final.spin, new_key.initial.spin, new_key.multipolarity)
        value, direction, function = b_down * ratio, observables.UP, "lifetimes.b_from_tau"
        truth = None if c.truth is None else ratio * observables.b_from_tau(
            c.truth, energy_kev, branching, conversion, mixing, new_key.multipolarity)
    elif (row.observable, to) == ("B(E2)", "tau"):
        old = TransitionKey.parse(row.key)
        new_key = old.initial
        ratio = observables.b_up_over_down(old.final.spin, old.initial.spin, old.multipolarity)
        value = observables.tau_from_b(c.value / ratio, energy_kev, branching, conversion, mixing,
                                       old.multipolarity)
        direction, function = None, "lifetimes.tau_from_b"
        truth = None if c.truth is None else observables.tau_from_b(
            c.truth / ratio, energy_kev, branching, conversion, mixing, old.multipolarity)
    else:
        raise ValueError(f"there is no conversion from {row.observable} to {to}")
    scale = abs(value / c.value)
    unc = {}
    if row.asymmetric:  # the sides swap: B falls when τ rises
        unc = {"plus": c.minus * scale, "minus": c.plus * scale}
    elif row.total is not None:
        unc = {"total": c.plus * scale}
    else:
        unc = {name: x * scale for name, x in (("stat", c.stat), ("sys", c.sys)) if x is not None}
    detail = (f"from {observables.get(row.observable).label} = {c.value:.5g} {c.unit} with Eγ = {energy_kev:g} keV, "
              f"branching {branching:g}, α = {conversion:g}" + (f", δ = {mixing:g}" if mixing is not None else ""))
    return Row(row.nuclide, to, new_key, value, observables.get(to).unit, "converted", row.evidence,
               direction=direction, detail=detail, truth=truth, depends_on=(row.id, *row.depends_on),
               origin={"kind": "conversion", "from": row.id, "function": function, **inputs}, **unc)


# ---------------------------------------------------------------------------------------------------------------
# The comparison


@dataclass(frozen=True)
class Entry:
    """One independent row beside the reference. Values are in the canonical convention."""

    row: Row
    id: str
    evidence: str
    method: str
    #: The method's label, and the row's own line of detail.
    label: str
    detail: str
    value: float
    plus: Optional[float]
    minus: Optional[float]
    #: The total uncertainty in percent of the value.
    uncertainty_percent: Optional[float]
    #: How the row was brought to the canonical convention; empty if it was in it.
    converted: str
    is_reference: bool = False
    #: The value minus the reference, in the unit and in percent of the reference, and the value over it.
    difference: Optional[float] = None
    difference_percent: Optional[float] = None
    ratio: Optional[float] = None
    #: The difference over the two uncertainties in quadrature, and what the rule makes of it. None for a row
    #: without an uncertainty, which has its ratio only.
    pull: Optional[float] = None
    verdict: Optional[str] = None
    #: False when the row and the reference share a source (``depends_on``): the agreement is then no test.
    independent_of_reference: Optional[bool] = None


@dataclass(frozen=True)
class Projection:
    """One projected row: what a simulated run, analysed back, says about its setup."""

    row: Row
    id: str
    method: str
    label: str
    detail: str
    value: float
    stat: Optional[float]
    sys: Optional[float]
    total: Optional[float]
    #: The precision the setup reaches: the total uncertainty in percent of the value.
    reach_percent: Optional[float]
    #: The value the simulation was given, and the closure: the value minus it, over the total uncertainty.
    truth: float
    closure: Optional[float]
    closure_verdict: Optional[str]
    #: How well the value is known already: the reference's precision and that of the best single measurement.
    known_percent: Optional[float]
    best_measured_percent: Optional[float]
    converted: str = ""


@dataclass(frozen=True)
class Comparison:
    """The rows of one observable of one transition or level, compared."""

    observable: str
    #: The observable's label in the canonical convention ("B(E2)↑"), its unit, and the level or transition.
    label: str
    unit: str
    key: str
    #: The reference (also among ``entries``), or None when there is no evaluated row and none was named.
    reference: Optional[Entry]
    #: The independent rows, the reference among them.
    entries: tuple
    projected: tuple
    #: The measured row with the smallest relative uncertainty, or None.
    best_measured: Optional[Entry]
    #: The rule the verdicts follow.
    agree_below: float = AGREE_BELOW
    differ_above: float = DIFFER_ABOVE

    def to_dict(self) -> dict:
        """The comparison as plain data (each row as :meth:`Row.to_dict`)."""
        return _plain(self)


def _plain(x):
    if isinstance(x, Row):
        return x.to_dict()
    if is_dataclass(x):
        return {f.name: _plain(getattr(x, f.name)) for f in fields(x)}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    return x


def _sources(row: Row, by_id: dict) -> set:
    """The row and everything it depends on, as far as the rows at hand say."""
    found, todo = set(), [row.id]
    while todo:
        i = todo.pop()
        if i in found:
            continue
        found.add(i)
        todo.extend(by_id[i].depends_on if i in by_id else ())
    return found


def _entry(c: Canonical, reference: Optional[Canonical], by_id: dict) -> Entry:
    row = c.row
    e = {"row": row, "id": row.id, "evidence": row.evidence, "method": row.method, "label": METHODS[row.method],
         "detail": row.detail, "value": c.value, "plus": c.plus, "minus": c.minus, "uncertainty_percent": c.percent,
         "converted": c.how}
    if reference is None:
        return Entry(**e)
    if reference.row.id == row.id:
        return Entry(**e, is_reference=True)
    difference = c.value - reference.value
    pull = None
    if c.plus is not None:  # the side of each uncertainty that faces the other value
        width = math.hypot(c.unc(reference.value), reference.unc(c.value) or 0.0)
        pull = difference / width if width > 0 else None
    return Entry(
        **e, difference=difference,
        difference_percent=100.0 * difference / abs(reference.value) if reference.value else None,
        ratio=c.value / reference.value if reference.value else None, pull=pull, verdict=verdict(pull),
        independent_of_reference=not (_sources(row, by_id) & _sources(reference.row, by_id)))


def _projection(c: Canonical, known: Optional[float], best: Optional[float]) -> Projection:
    row = c.row
    width = c.unc(c.truth)
    closure = (c.value - c.truth) / width if width else None
    return Projection(row, row.id, row.method, METHODS[row.method], row.detail, c.value, c.stat, c.sys, c.unc(),
                      c.percent, c.truth, closure, verdict(closure), known, best, c.how)


def compare(rows: list, reference: Union[Row, str, None] = None) -> Comparison:
    """Compare the rows of one observable of one transition or level.

    ``rows`` may hold superseded rows; only those that hold now are compared. ``reference`` is a row or its id;
    left out, the evaluated row is the reference (the newest, if there are several), and without one the rows
    are listed with no differences. A projected row cannot be the reference. Raises ``ValueError`` for rows of
    different observables or of different levels or transitions.
    """
    rows = list(rows)
    if not rows:
        raise ValueError("there are no rows to compare")
    if len({r.observable for r in rows}) > 1 or len({r.key for r in rows}) > 1:
        raise ValueError("the rows are of different observables or of different levels or transitions: "
                         + ", ".join(sorted({f"{r.observable} of {r.key}" for r in rows})))
    by_id = {r.id: r for r in rows}
    now = [canonical(r) for r in current(rows)]
    independent = [c for c in now if c.row.independent]
    if reference is None:
        evaluated = sorted((c for c in independent if c.row.evidence == "evaluated"), key=lambda c: c.row.created)
        ref = evaluated[-1] if evaluated else None
    else:
        wanted = reference.id if isinstance(reference, Row) else str(reference)
        ref = next((c for c in now if c.row.id == wanted), None)
        if ref is None:
            raise ValueError(f"the reference {wanted} is not among the rows that hold now")
        if ref.row.evidence == PROJECTED:
            raise ValueError("a projected row cannot be the reference: it is a simulated run, not a value to test")
    entries = tuple(_entry(c, ref, by_id) for c in independent)
    measured = [e for e in entries if e.evidence == "measured" and e.uncertainty_percent is not None]
    best = min(measured, key=lambda e: e.uncertainty_percent) if measured else None
    o = observables.get(rows[0].observable)
    return Comparison(
        o.key, o.symbol, o.unit, rows[0].key, next((e for e in entries if e.is_reference), None), entries,
        tuple(_projection(c, None if ref is None else ref.percent,
                          None if best is None else best.uncertainty_percent)
              for c in now if c.row.evidence == PROJECTED),
        best)


__all__ = ["AGREE_BELOW", "DIFFER_ABOVE", "VERDICTS", "Canonical", "Comparison", "Entry", "Projection",
           "canonical", "compare", "convert", "verdict"]
