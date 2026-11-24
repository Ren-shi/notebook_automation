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
    [(u.text, u.reference) for u in obs.get("Qs").used_for]   # what the value tells, each with its reference
    obs.rotor_reading(884.0, "62Ni")                   # β₂, Q₀ and Q_s(2⁺) of a rigid rotor with this B(E2)↑

It holds B(E2), τ, Q_s, B(M1), μ and level energies; cross sections follow.

**No formula of its own.** Every conversion here calls a function that existed before this module
(:meth:`LevelScheme.b <physim.nuclear.levels.LevelScheme.b>` for B = ⟨f‖M‖i⟩² / (2J_i + 1),
:func:`~physim.nuclear.levels.weisskopf_unit` and
:data:`~physim.nuclear.shellmodel_import.WEISSKOPF_M1`, :func:`~physim.nuclear.levels.quadrupole_factor`,
:func:`lifetimes.tau_from_b <physim.nuclear.lifetimes.tau_from_b>` and
:func:`lifetimes.b_from_tau <physim.nuclear.lifetimes.b_from_tau>`) and the unit tables beside them. The rotor
reading of a B(E2) (β₂, Q₀ and the rotor's Q_s) is the Coulex analysis's own
(:meth:`Analysis._shape <physim.nuclear.analysis.Analysis>`), called as it stands (:func:`rotor_reading`). What
the code has no function for is left out: half-life ↔ mean life (written inline where it is used), a B(M1) ↔ τ
conversion in the comparison engine, and everything of backlog item 76.

**The text** (:attr:`Observable.used_for`) is the known physics in a few plain lines per use, each with a
reference: a textbook or a review, most of them cited elsewhere in physim's documentation.

Conventions: a B value is compared as B↑ (:data:`UP`, the excitation, lower level → upper level): B(E2)↑ in
e² fm⁴ and B(M1)↑ in μ_N²; a value in Weisskopf units is always the ↓ value (the decay) and is labelled so. τ is
the level's mean life in ps, Q_s in e fm² and μ in μ_N with their signs, level energies in keV.
"""

from __future__ import annotations

import math
import re
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
class Use:
    """One thing a value is used for: what it tells and where it is used, in a few plain lines, and the
    reference the lines rest on (authors, the journal or book, the year)."""

    text: str
    reference: str

    def to_dict(self) -> dict:
        return {"text": self.text, "reference": self.reference}


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
    #: What the value is used for: :class:`Use` entries, each a few lines with a reference.
    used_for: tuple = ()

    @property
    def symbol(self) -> str:
        """The label with the canonical direction: "B(E2)↑"."""
        return self.label + (DIRECTIONS[self.direction] if self.direction else "")

    @property
    def multipolarity(self) -> Optional[str]:
        """The multipolarity of a transition's observable ("E2" for B(E2), "M1" for B(M1)); None for a level's."""
        m = re.search(r"\(([EM][1-6])\)", self.key)
        return m.group(1) if m else None


_PS = ensdf.TIME_UNITS["PS"]

# The references of the text. The first six are cited elsewhere in physim (docs/theory/evaporation.md twice,
# residues.py, docs/theory/coulex.md, shellmodel_spaces.py and docs/shell-model.md); the others are the standard
# textbook, reviews and tables of their subject. No page numbers: none could be checked against the repository.
BLATT_WEISSKOPF = "J. M. Blatt and V. F. Weisskopf, Theoretical Nuclear Physics (1952)"
BOHR_MOTTELSON_1 = "A. Bohr and B. R. Mottelson, Nuclear Structure, Vol. I (1969)"
BOHR_MOTTELSON_2 = "A. Bohr and B. R. Mottelson, Nuclear Structure, Vol. II (1975)"
ALDER_WINTHER = "K. Alder and A. Winther, Electromagnetic Excitation (1975)"
RICHTER_2008 = "W. A. Richter, S. Mkhize and B. A. Brown, Phys. Rev. C 78, 064302 (2008)"
BROWN_RICHTER_2006 = "B. A. Brown and W. A. Richter, Phys. Rev. C 74, 034315 (2006)"
DEWALD_2012 = "A. Dewald, O. Möller and P. Petkov, Prog. Part. Nucl. Phys. 67 (2012)"
ALEXANDER_FORSTER_1978 = "T. K. Alexander and J. S. Forster, Adv. Nucl. Phys. 10 (1978)"
CASTEN = "R. F. Casten, Nuclear Structure from a Simple Perspective, 2nd ed. (2000)"
RAMAN_2001 = "S. Raman, C. W. Nestor Jr. and P. Tikkanen, At. Data Nucl. Data Tables 78 (2001)"
STONE_2005 = "N. J. Stone, At. Data Nucl. Data Tables 90 (2005)"

_VIA_TAU = "through τ, the γ-ray energy, the branching and the conversion"

OBSERVABLES = {o.key: o for o in (
    Observable(
        "B(E2)", "B(E2)", "transition", "e2fm4", {**B_UNITS[2], WEISSKOPF: None}, UP,
        "B(E2)↑ in e² fm⁴, from the lower level to the upper one; in Weisskopf units always the ↓ value",
        routes=(
            Route("coulex-yield", "coulex", "the γ-ray yield per scattered particle, first order"),
            Route("coulex-fit", "coulex", "the matrix elements fitted to the yields, to all orders"),
            Route("rdds", "lifetimes", _VIA_TAU),
            Route("dsam", "lifetimes", _VIA_TAU),
            Route("inelastic-scattering", "inelastic",
                  "nuclear excitation above the barrier, (p,p′) and the like: it gives the deformation length, "
                  "and a B(E2) only through a model"),
            Route("shell-model"), Route("evaluation"),
            Route("nrf", note="photon scattering from the ground state: the level's width to the ground state"),
            Route("inelastic-neutron", note="through τ (DSAM after (n,n′γ))"),
            Route("electron-scattering", note="the (e,e′) form factor at low momentum transfer"),
        ),
        conversions=(
            Conversion("B(E2)↓", "b_up_over_down", "LevelScheme.b", ("j_lower", "j_upper")),
            Conversion("B(E2)↓ in W.u.", "weisskopf", "levels.weisskopf_unit", ("nuclide",)),
            Conversion("matrix element", "matrix_element", "LevelScheme.b", ("j_initial",)),
            Conversion("tau", "tau_from_b", "lifetimes.tau_from_b", ("energy_kev", "branching", "conversion")),
        ),
        used_for=(
            Use("Collectivity. In Weisskopf units (always the ↓ value) it counts single-particle strengths: about "
                "1 W.u. is one nucleon changing its orbit, ten or more mean many nucleons moving together, and a "
                "well-deformed rotor reaches 100 W.u. and more.", BLATT_WEISSKOPF),
            Use("The lifetime. With the γ-ray energy, the branching and the conversion coefficient it fixes the "
                "level's mean life (for a pure E2 decay τ falls as 1/(E_γ⁵ B(E2))), and so whether the γ ray is "
                "emitted in flight and which lifetime method reaches the level.", BLATT_WEISSKOPF),
            Use("The size of the deformation. Read as a rigid axial rotor, B(E2; 0⁺ → 2⁺₁) gives the intrinsic "
                "quadrupole moment Q₀ and β₂. It gives the size only, never the sign, and the reading holds only "
                "where the level scheme is rotational (E(4⁺₁)/E(2⁺₁) near 3.33).", BOHR_MOTTELSON_2),
            Use("A test of a shell-model interaction and its effective charges. A calculated B(E2) rests on the "
                "wave functions and on the effective charges that stand for the excitations left out of the "
                "valence space; a measured value tests the two together, most sharply over several transitions "
                "or a chain of isotopes.", RICHTER_2008),
            Use("An input to other experiments. A known B(E2) of the collision partner normalises a "
                "Coulomb-excitation yield, and the B(E2) of this nucleus sets the yield and the beam time of "
                "every Coulomb-excitation experiment on it.", ALDER_WINTHER),
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
            Route("nrf", note="from the level's width"), Route("inelastic-neutron", note="DSAM after (n,n′γ)"),
        ),
        conversions=(
            Conversion("B(E2)", "b_from_tau", "lifetimes.b_from_tau", ("energy_kev", "branching", "conversion")),
        ),
        used_for=(
            Use("The transition strengths. With the branching ratios, the conversion coefficients and the mixing "
                "ratios, τ gives the B(E2) and B(M1) of every γ ray from the level with no excitation cross "
                "section in between: a route independent of Coulomb excitation.", DEWALD_2012),
            Use("The choice of the method. From about 1 ps to 1 ns the nucleus decays in flight behind the target "
                "and a plunger measures the decay curve (RDDS); from about 10 fs to 1 ps it decays while it slows "
                "down in a backing and the line shape carries τ (DSAM).",
                f"{DEWALD_2012}; {ALEXANDER_FORSTER_1978}"),
            Use("The Doppler correction. τ against the time the recoil needs to leave the target or to stop "
                "decides whether a γ ray is emitted in flight (fully shifted), at rest, or in between, which is "
                "the tail of the line.", ALEXANDER_FORSTER_1978),
        )),
    Observable(
        "Qs", "Q_s", "level", "efm2", dict(ME_UNITS[2]), None,
        "the spectroscopic quadrupole moment in e fm², with its sign",
        routes=(
            Route("coulex-fit", "coulex", "the reorientation effect: the diagonal matrix element"),
            Route("shell-model"), Route("evaluation"),
            Route("moment", note="for a ground or long-lived state: laser spectroscopy, atomic beams, muonic "
                                 "X rays"),
        ),
        conversions=(
            Conversion("matrix element", "diagonal_matrix_element", "levels.quadrupole_factor", ("j",)),
        ),
        used_for=(
            Use("The shape and its sign. For the 2⁺₁ state of a K = 0 rotor Q_s = −(2/7) Q₀: negative for a "
                "prolate nucleus, positive for an oblate one; a harmonic vibrator has Q_s = 0. B(E2) gives the "
                "size of the deformation, Q_s its sign.", BOHR_MOTTELSON_2),
            Use("How it is measured. In Coulomb excitation the diagonal matrix element changes the yield of the "
                "state with the scattering angle and with the collision partner (the reorientation effect), so "
                "Q_s comes from a fit to all orders, together with the interference through higher states.",
                ALDER_WINTHER),
            Use("A test of the shell model. Q_s rests on the same effective charges as B(E2) but on one state's "
                "own wave function, so its sign and size test that state.", RICHTER_2008),
        )),
    Observable(
        "B(M1)", "B(M1)", "transition", "mu_N2", {"mu_N2": 1.0, WEISSKOPF: None}, UP,
        "B(M1)↑ in μ_N², from the lower level to the upper one; in Weisskopf units always the ↓ value",
        routes=(
            Route("rdds", "lifetimes", _VIA_TAU + ", and the E2/M1 mixing ratio"),
            Route("dsam", "lifetimes", _VIA_TAU + ", and the E2/M1 mixing ratio"),
            Route("shell-model"), Route("evaluation"),
            Route("nrf", note="photon scattering from the ground state: the level's width to the ground state"),
            Route("inelastic-neutron", note="through τ (DSAM after (n,n′γ))"),
            Route("electron-scattering", note="(e,e′) at backward angles"),
        ),
        conversions=(
            Conversion("B(M1)↓", "b_up_over_down", "LevelScheme.b", ("j_lower", "j_upper")),
            Conversion("B(M1)↓ in W.u.", "weisskopf", "shellmodel_import.WEISSKOPF_M1", ()),
        ),
        used_for=(
            Use("The magnetic strength of a transition. It follows the spin and orbital g factors of the few "
                "nucleons that change their orbit, where B(E2) follows the collective motion of the charge. "
                "1 W.u. is 1.79 μ_N²; a strong M1 transition is of that order.", BLATT_WEISSKOPF),
            Use("With the B(E2) of the same γ ray it fixes the size of the E2/M1 mixing ratio; in the shell model "
                "it tests the effective g factors as B(E2) tests the effective charges.", RICHTER_2008),
        )),
    Observable(
        "mu", "μ", "level", "mu_N", {"mu_N": 1.0}, None,
        "the magnetic dipole moment in nuclear magnetons μ_N, with its sign",
        routes=(
            Route("moment", note="for an excited state the transient field, recoil in vacuum or perturbed "
                                 "angular correlations; for a ground or long-lived state laser spectroscopy, "
                                 "NMR or atomic beams"),
            Route("shell-model"), Route("evaluation"),
        ),
        used_for=(
            Use("Which orbits the state is made of. μ = g J μ_N: one nucleon in one orbit has the Schmidt value "
                "of that orbit and a collective state has a g factor near Z/A, so the measured g factor says how "
                "much of each the state holds.", f"{BOHR_MOTTELSON_1}; {BOHR_MOTTELSON_2}"),
            Use("A test of a shell-model wave function and of its effective g factors, which a small admixture "
                "of another orbit changes more than it changes the energy.", RICHTER_2008),
            Use("Where the measured values are: the table of nuclear magnetic dipole and electric quadrupole "
                "moments.", STONE_2005),
        )),
    Observable(
        "E_level", "E", "level", "keV", {"keV": 1.0, "MeV": 1e3}, None,
        "the level's excitation energy in keV",
        routes=(
            Route("shell-model"), Route("evaluation"),
            Route("transfer", "transfer", "the levels a transfer reaction populates"),
            Route("inelastic-scattering", "inelastic", "the levels inelastic scattering excites"),
        ),
        used_for=(
            Use("E(4⁺₁)/E(2⁺₁) of an even-even nucleus names the kind of collective motion: 2.0 for a harmonic "
                "vibrator, 2.5 for a γ-soft rotor, 3.33 for a rigid axial rotor, and below 2 near closed shells, "
                "where few nucleons carry the states.", CASTEN),
            Use("E(2⁺₁) itself. It is high at a closed shell and low in a deformed nucleus, and it falls as "
                "B(E2; 0⁺ → 2⁺₁) rises: the two together are the first sign of a shell closure or of its loss.",
                RAMAN_2001),
            Use("A test of a shell-model interaction. Level energies are what an interaction is fitted to, so "
                "calculated minus measured energies say how far it is to be trusted for this nucleus before its "
                "B(E2) values and moments are read.", BROWN_RICHTER_2006),
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
    """B(λ↑) / B(λ↓) = (2J_upper + 1) / (2J_lower + 1), from :meth:`LevelScheme.b`: 5 for 0⁺ → 2⁺. The ratio
    holds the spins only, so for a magnetic multipole (a level scheme keeps electric matrix elements) it is taken
    from the electric one of the same order."""
    multipolarity = "E" + multipolarity[1:] if multipolarity.startswith("M") else multipolarity
    scheme = _pair(j_lower, j_upper, multipolarity)
    return scheme.b(0, 1, multipolarity) / scheme.b(1, 0, multipolarity)


def weisskopf_unit_of(nuclide: str, multipolarity: str = "E2") -> float:
    """One Weisskopf unit of B(λ↓): for Eλ :func:`levels.weisskopf_unit` in e² fm^{2λ} (it grows with the mass
    number); for M1 the constant :data:`physim.nuclear.shellmodel_import.WEISSKOPF_M1`, 1.790 μ_N². The code has
    no Weisskopf unit for a higher magnetic multipole."""
    if multipolarity == "M1":
        from .shellmodel_import import WEISSKOPF_M1

        return WEISSKOPF_M1
    if not multipolarity.startswith("E"):
        raise ValueError(f"physim has no Weisskopf unit for {multipolarity}")
    return weisskopf_unit(int(multipolarity[1]), parse_nuclide(nuclide)[1])


def weisskopf(b_down: float, nuclide: str, multipolarity: str = "E2") -> float:
    """B(λ↓) in Weisskopf units from its value in e² fm^{2λ} (or μ_N² for M1): :func:`weisskopf_unit_of`."""
    return b_down / weisskopf_unit_of(nuclide, multipolarity)


def from_weisskopf(b_wu: float, nuclide: str, multipolarity: str = "E2") -> float:
    """B(λ↓) in e² fm^{2λ} (or μ_N² for M1) from its value in Weisskopf units."""
    return b_wu * weisskopf_unit_of(nuclide, multipolarity)


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


def rotor_reading(b_up_e2fm4: float, nuclide: str) -> dict:
    """B(E2; 0⁺ → 2⁺₁)↑ of an even-even nucleus read as a rigid axially symmetric rotor: {"beta2", "q0_efm2",
    "qs_2plus_efm2" (the prolate rotor's Q_s(2⁺₁), negative; the oblate one's has the other sign), "note"}.
    Empty for a B that is not above 0.

    The numbers are the Coulomb-excitation analysis's own: its method
    :meth:`Analysis._shape <physim.nuclear.analysis.Analysis>` (β₂ = 4π √B / (3 Z R₀²) with R₀ = 1.2 A^⅓ fm,
    Q₀ = √(16π B / 5), Q_s(2⁺) = −(2/7) Q₀) is called as it stands, with a stand-in for the analysis that holds
    the nucleus alone, so no formula is written a second time and a page shows what the Analysis step shows.
    The method is not public: if it comes to need more of an analysis than the nucleus, this raises and the
    check that uses it is left out (:func:`physim.nuclear.nucleus.checks`).
    """
    from types import SimpleNamespace

    from .analysis import Analysis

    z, a = parse_nuclide(nuclide)
    stand_in = SimpleNamespace(experiment=SimpleNamespace(
        beam=SimpleNamespace(Z=z, A=a), excitation=SimpleNamespace(excite="projectile"), levels={}))
    return Analysis._shape(stand_in, float(b_up_e2fm4), a)


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
        factor = from_weisskopf(1.0, nuclide, o.multipolarity or "E2")
        how.append(f"1 W.u. = {factor:.4g} {o.unit}" + (f" for {nuclide}" if o.multipolarity != "M1" else ""))
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
            ratio = b_up_over_down(j_lower, j_upper, o.multipolarity or "E2")
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


__all__ = ["DIRECTIONS", "DOWN", "OBSERVABLES", "UP", "WEISSKOPF", "Conversion", "Observable", "Route", "Use",
           "b_from_tau", "b_up_over_down", "canonical_factor", "diagonal_matrix_element", "from_weisskopf", "get",
           "matrix_element", "quadrupole_moment", "rotor_reading", "routes", "tau_from_b", "to_canonical",
           "weisskopf", "weisskopf_unit_of"]
