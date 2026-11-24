"""The registry of observables: what each quantity is, how it can be reached and what it converts to.

The known physics of the results pages, written down once as data (:data:`OBSERVABLES`), with thin functions
around the conversions the code already has. The results ledger (:mod:`physim.nuclear.results`) checks its rows
against it, and the comparison engine (:mod:`physim.nuclear.compare`) converts rows to each observable's
canonical convention with it::

    from physim.nuclear import observables as obs

    obs.get("B(E2)").unit                              # "e2fm4": B(E2)↑ in e² fm⁴ is what is compared
    obs.to_canonical("B(E2)", 12.1, "W.u.", "down", nuclide="62Ni", j_upper=2, j_lower=0)
    obs.b_from_tau(2.09, 1172.9)                       # B(E2↓) in e² fm⁴ from τ in ps
    [r.kind for r in obs.get("tau").routes if r.live]  # the kinds of experiment that plan it today

Started with B(E2), τ, Q_s and level energies; B(M1), μ and cross sections follow.

**No formula of its own.** Every conversion here calls a function that existed before this module
(:meth:`LevelScheme.b <physim.nuclear.levels.LevelScheme.b>` for B = ⟨f‖M‖i⟩² / (2J_i + 1),
:func:`~physim.nuclear.levels.weisskopf_unit`, :func:`~physim.nuclear.levels.quadrupole_factor`,
:func:`lifetimes.tau_from_b <physim.nuclear.lifetimes.tau_from_b>` and
:func:`lifetimes.b_from_tau <physim.nuclear.lifetimes.b_from_tau>`) and the unit tables beside them. What the code
holds only inside another class is left out: the rotor reading of a B(E2) (β₂, Q₀ and the rotor's Q_s, in the
Coulex analysis) and half-life ↔ mean life (written inline where it is used).

Conventions: B(E2) is compared as B(E2)↑ in e² fm⁴ (:data:`UP`, the excitation, lower level → upper level); a
value in Weisskopf units is always the ↓ value (the decay) and is labelled so. τ is the level's mean life in ps,
Q_s in e fm², level energies in keV.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from . import ensdf, lifetimes
from .coulex import B_UNITS
from .levels import ME_UNITS, Level, LevelScheme, Value, quadrupole_factor, weisskopf_unit
from .names import parse_nuclide

#: The two directions of a reduced transition probability: ↑ from the lower level to the upper one (the
#: excitation), ↓ from the upper to the lower (the decay).
UP, DOWN = "up", "down"
DIRECTIONS = {UP: "↑", DOWN: "↓"}
#: The unit that is a convention as well: a B in Weisskopf units is always B↓.
WEISSKOPF = "W.u."


@dataclass(frozen=True)
class Route:
    """One way to reach an observable: the method (a key of :data:`physim.nuclear.results.METHODS`) and the kind
    of experiment that plans it (a key of :data:`physim.nuclear.workbench.KINDS`; None for a route that is not an
    experiment of the planner: a calculation, the evaluation, a method known from the literature only)."""

    method: str
    kind: Optional[str] = None
    #: A line on how the route gets the value ("through τ and the branching").
    note: str = ""

    @property
    def live(self) -> Optional[bool]:
        """Whether physim has this kind of experiment today (False: a later kind; None: not an experiment)."""
        if self.kind is None:
            return None
        from .workbench import KINDS

        return next(k.live for k in KINDS if k.key == self.kind)


@dataclass(frozen=True)
class Conversion:
    """What an observable gives: another observable or convention, the function of this module that converts
    to it, the function it calls (which existed before), and what it needs besides the value."""

    to: str
    function: str
    calls: str
    needs: tuple = ()


@dataclass(frozen=True)
class Observable:
    key: str
    label: str
    #: "transition" or "level": what the key of a row of this observable names.
    belongs_to: str
    #: The canonical unit, and every unit a row may come in with its value in the canonical one (None for
    #: Weisskopf units, whose value depends on the mass number).
    unit: str
    units: dict
    #: The canonical direction (B↑), or None for an observable without one.
    direction: Optional[str]
    #: The convention in words, for a column heading or a footnote.
    convention: str
    routes: tuple = ()
    conversions: tuple = ()
    #: What the value is used for. Nearly empty for now: later waves fill it and the owner reviews the text.
    used_for: str = ""

    @property
    def symbol(self) -> str:
        """The label with the canonical direction: "B(E2)↑"."""
        return self.label + (DIRECTIONS[self.direction] if self.direction else "")


_PS = ensdf.TIME_UNITS["PS"]

OBSERVABLES = {o.key: o for o in (
    Observable(
        "B(E2)", "B(E2)", "transition", "e2fm4", {**B_UNITS[2], WEISSKOPF: None}, UP,
        "B(E2)↑ in e² fm⁴, from the lower level to the upper one; in Weisskopf units always the ↓ value",
        routes=(
            Route("coulex-yield", "coulex", "the γ-ray yield per scattered particle, first order"),
            Route("coulex-fit", "coulex", "the matrix elements fitted to the yields, to all orders"),
            Route("rdds", "lifetimes", "through τ, the γ-ray energy, the branching and the conversion"),
            Route("dsam", "lifetimes", "through τ, the γ-ray energy, the branching and the conversion"),
            Route("other", "inelastic", "nuclear excitation above the barrier"),
            Route("shell-model"), Route("evaluation"),
            Route("nrf"), Route("inelastic-neutron", note="through τ (DSAM after (n,n′γ))"),
            Route("electron-scattering"),
        ),
        conversions=(
            Conversion("B(E2)↓", "b_up_over_down", "LevelScheme.b", ("j_lower", "j_upper")),
            Conversion("B(E2)↓ in W.u.", "weisskopf", "levels.weisskopf_unit", ("nuclide",)),
            Conversion("matrix element", "matrix_element", "LevelScheme.b", ("j_initial",)),
            Conversion("tau", "tau_from_b", "lifetimes.tau_from_b", ("energy_kev", "branching", "conversion")),
        )),
    Observable(
        "tau", "τ", "level", "ps", {u.lower(): ensdf.TIME_UNITS[u] / _PS for u in ("FS", "PS", "NS", "S")}, None,
        "the level's mean life in ps (not the half-life)",
        routes=(
            Route("rdds", "lifetimes", "the decay curve over plunger distances"),
            Route("dsam", "lifetimes", "the Doppler-broadened line shape"),
            Route("coulex-yield", "coulex", "from B(E2), the γ-ray energy, the branching and the conversion"),
            Route("coulex-fit", "coulex", "from B(E2), the γ-ray energy, the branching and the conversion"),
            Route("shell-model", note="from the calculated B values"), Route("evaluation"),
            Route("nrf"), Route("inelastic-neutron"),
        ),
        conversions=(
            Conversion("B(E2)", "b_from_tau", "lifetimes.b_from_tau", ("energy_kev", "branching", "conversion")),
        )),
    Observable(
        "Qs", "Q_s", "level", "efm2", dict(ME_UNITS[2]), None,
        "the spectroscopic quadrupole moment in e fm², with its sign",
        routes=(
            Route("coulex-fit", "coulex", "the reorientation effect: the diagonal matrix element"),
            Route("shell-model"), Route("evaluation"),
        ),
        conversions=(
            Conversion("matrix element", "diagonal_matrix_element", "levels.quadrupole_factor", ("j",)),
        )),
    Observable(
        "E_level", "E", "level", "keV", {"keV": 1.0, "MeV": 1e3}, None,
        "the level's excitation energy in keV",
        routes=(
            Route("shell-model"), Route("evaluation"),
            Route("other", "transfer", "the levels a transfer reaction populates"),
            Route("other", "inelastic"),
        )),
)}


def get(key: str) -> Observable:
    """The registry's entry for an observable. Raises ``ValueError`` for one it does not have."""
    try:
        return OBSERVABLES[key]
    except KeyError:
        raise ValueError(f"'{key}' is not an observable of the registry; it has {', '.join(OBSERVABLES)}") from None


def routes(key: str, live: Optional[bool] = None) -> list:
    """The routes to an observable; with ``live`` True those of a kind of experiment physim has today, False
    those of a later kind."""
    found = get(key).routes
    return [r for r in found if live is None or (r.kind is not None and r.live is live)]


# ---------------------------------------------------------------------------------------------------------------
# Conversions: thin functions around what the code already has


def _pair(j_lower: float, j_upper: float, multipolarity: str = "E2") -> LevelScheme:
    """Two levels joined by a unit matrix element, so that the existing :meth:`LevelScheme.b` gives
    B / ⟨f‖M‖i⟩² in each direction (level 0 the lower, level 1 the upper)."""
    scheme = LevelScheme("", [Level(Value(0.0, source="assumed"), float(j_lower)),
                              Level(Value(1.0, source="assumed"), float(j_upper))])
    scheme.set_matrix_element(0, 1, multipolarity, 1.0)
    return scheme


def b_up_over_down(j_lower: float, j_upper: float, multipolarity: str = "E2") -> float:
    """B(Eλ↑) / B(Eλ↓) = (2J_upper + 1) / (2J_lower + 1), from :meth:`LevelScheme.b`: 5 for 0⁺ → 2⁺."""
    scheme = _pair(j_lower, j_upper, multipolarity)
    return scheme.b(0, 1, multipolarity) / scheme.b(1, 0, multipolarity)


def weisskopf(b_down: float, nuclide: str, multipolarity: str = "E2") -> float:
    """B(Eλ↓) in Weisskopf units from its value in e² fm^{2λ} (:func:`levels.weisskopf_unit`)."""
    return b_down / weisskopf_unit(int(multipolarity[1]), parse_nuclide(nuclide)[1])


def from_weisskopf(b_wu: float, nuclide: str, multipolarity: str = "E2") -> float:
    """B(Eλ↓) in e² fm^{2λ} from its value in Weisskopf units."""
    return b_wu * weisskopf_unit(int(multipolarity[1]), parse_nuclide(nuclide)[1])


def matrix_element(b: float, j_initial: float, multipolarity: str = "E2") -> float:
    """The size of ⟨f‖M(Eλ)‖i⟩ in e fm^λ from B(Eλ; i → f) = ⟨f‖M‖i⟩² / (2J_i + 1); the sign is not in a B."""
    return math.sqrt(b / _pair(j_initial, j_initial, multipolarity).b(0, 1, multipolarity))


def tau_from_b(b_down: float, energy_kev: float, branching: float = 1.0, conversion: float = 0.0,
               mixing: Optional[float] = None, multipolarity: str = "E2") -> float:
    """The level's mean life in ps from B(λ↓) of one of its transitions: :func:`lifetimes.tau_from_b`, where
    ``branching``, ``conversion`` and ``mixing`` are explained."""
    return lifetimes.tau_from_b(energy_kev, multipolarity, b_down, branching, conversion, mixing)


def b_from_tau(tau_ps: float, energy_kev: float, branching: float = 1.0, conversion: float = 0.0,
               mixing: Optional[float] = None, multipolarity: str = "E2") -> float:
    """B(λ↓) in e² fm^{2λ} from the level's mean life in ps: :func:`lifetimes.b_from_tau`."""
    return lifetimes.b_from_tau(energy_kev, multipolarity, tau_ps, branching, conversion, mixing)


def diagonal_matrix_element(qs_efm2: float, j: float) -> float:
    """⟨J‖M(E2)‖J⟩ in e fm² from the spectroscopic quadrupole moment (:func:`levels.quadrupole_factor`)."""
    factor = quadrupole_factor(j)
    if factor <= 0:
        raise ValueError(f"a level of spin {j:g} has no spectroscopic quadrupole moment")
    return qs_efm2 / factor


def quadrupole_moment(diagonal_efm2: float, j: float) -> float:
    """Q_s in e fm² from the diagonal matrix element ⟨J‖M(E2)‖J⟩."""
    return diagonal_efm2 * quadrupole_factor(j)


def canonical_factor(key: str, unit: str, direction: Optional[str] = None, *, nuclide: Optional[str] = None,
                     j_lower: Optional[float] = None, j_upper: Optional[float] = None) -> tuple:
    """(factor, how): the number that turns a value of the observable ``key`` given in ``unit`` (and, for a B,
    ``direction``) into its canonical convention, and the steps in words (empty when nothing was converted).

    A value in Weisskopf units needs the ``nuclide``; B↓ needs the spins of the two levels. Raises ``ValueError``
    when the unit is not one of the observable's or something needed is missing.
    """
    o = get(key)
    if unit not in o.units:
        raise ValueError(f"'{unit}' is not a unit of {o.label}; use {', '.join(o.units)}")
    how = []
    if unit == WEISSKOPF:
        if direction != DOWN:
            raise ValueError(f"{o.label} in Weisskopf units is always the ↓ value")
        if nuclide is None:
            raise ValueError("a value in Weisskopf units needs the nuclide: the unit depends on the mass number")
        factor = from_weisskopf(1.0, nuclide)
        how.append(f"1 W.u. = {factor:.4g} {o.unit} for {nuclide}")
    else:
        factor = o.units[unit]
        if factor != 1.0:
            how.append(f"1 {unit} = {factor:g} {o.unit}")
    if o.direction is not None:
        if direction not in DIRECTIONS:
            raise ValueError(f"{o.label} needs its direction stated: '{UP}' (↑) or '{DOWN}' (↓)")
        if direction != o.direction:
            if j_lower is None or j_upper is None:
                raise ValueError(f"turning {o.label}↓ into {o.label}↑ needs the spins of both levels")
            ratio = b_up_over_down(j_lower, j_upper)
            ratio = ratio if o.direction == UP else 1.0 / ratio
            factor *= ratio
            how.append(f"{o.label}{DIRECTIONS[o.direction]} = {ratio:g} × {o.label}{DIRECTIONS[direction]} "
                       f"(J = {j_lower:g} and {j_upper:g})")
    elif direction is not None:
        raise ValueError(f"{o.label} has no direction")
    return factor, "; ".join(how)


def to_canonical(key: str, value: float, unit: str, direction: Optional[str] = None, **needs) -> float:
    """``value`` in the observable's canonical convention (:func:`canonical_factor`)."""
    return value * canonical_factor(key, unit, direction, **needs)[0]


__all__ = ["DIRECTIONS", "DOWN", "OBSERVABLES", "UP", "WEISSKOPF", "Conversion", "Observable", "Route",
           "b_from_tau", "b_up_over_down", "canonical_factor", "diagonal_matrix_element", "from_weisskopf", "get",
           "matrix_element", "quadrupole_moment", "routes", "tau_from_b", "to_canonical", "weisskopf"]
