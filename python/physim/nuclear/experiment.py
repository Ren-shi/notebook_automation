"""The experiment setup: beam, target, detectors and run conditions, read from and written to a TOML file.

Coordinates: the beam travels along +z through the target at the origin; x is horizontal and y vertical. A
detector at polar angle θ (from +z) and azimuth φ (from +x towards +y) sits at distance d from the target, facing
it unless ``facing`` says otherwise.

Every field is written with its unit (``"5.5 MeV"``, ``"40 mm"``), so a setup file reads unambiguously. Reading
collects every problem in the file and reports them together, each naming the field::

    exp = Experiment.example("alpha_on_gold")  # or Experiment.load("my_setup.toml")
    exp.beam.energy            # Quantity(5.5, 'MeV')
    exp.beam.energy_mev        # total kinetic energy in MeV
    exp.detectors[0].theta = "60 deg"
    exp.validate()             # raises SetupError listing every problem
    exp.save("my_setup.toml")
"""

from __future__ import annotations

import difflib
import math
from dataclasses import dataclass, field, fields
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Optional, Union

from . import _toml, catalogue
from . import data as _data
from .names import parse_material, parse_nuclide
from .quantity import Quantity

#: Schema identifier written at the top of every setup file.
SCHEMA = "physim.experiment/1"

#: Reactions this version can plan. Others are added slice by slice (backlog item 43 onwards).
REACTIONS = ("elastic", "coulex", "fusion-evaporation")

#: Folder of the example setups shipped with the package.
EXAMPLES = Path(__file__).resolve().parent / "examples"


def example_names() -> list[str]:
    """Names of the example setups, for :meth:`Experiment.example`."""
    return sorted(p.stem for p in EXAMPLES.glob("*.toml"))


#: Detector shapes and the size fields each one uses.
SHAPES = {
    "rectangle": ("width", "height", "strips_x", "strips_y"),
    "annular": ("inner_radius", "outer_radius", "rings", "sectors"),
    "circle": ("radius",),
}

QuantityLike = Union[Quantity, str]


class SetupError(ValueError):
    """A setup is invalid. ``problems`` lists every problem found, each naming the field."""

    def __init__(self, problems: list[str]):
        self.problems = list(problems)
        n = len(self.problems)
        head = "the setup has 1 problem:" if n == 1 else f"the setup has {n} problems:"
        super().__init__("\n".join([head] + [f"  - {p}" for p in self.problems]))


# ---------------------------------------------------------------------------------------------------------------
# Field rules


def _positive(q: Quantity) -> Optional[str]:
    return None if q.value > 0 else "must be positive"


def _non_negative(q: Quantity) -> Optional[str]:
    return None if q.value >= 0 else "must not be negative"


def _polar(q: Quantity) -> Optional[str]:
    return None if 0 <= q.to("deg") <= 180 else "must be between 0 and 180 deg"


def _tilt(q: Quantity) -> Optional[str]:
    return None if abs(q.to("deg")) < 90 else "must be between -90 and 90 deg (exclusive)"


@dataclass(frozen=True)
class _Field:
    name: str
    kind: str  # a quantity kind list "energy|energy_per_nucleon", or "str", "int", "nuclide", "material", ...
    required: bool = False
    check: Optional[Callable[[Any], Optional[str]]] = None


def _read_value(spec: _Field, raw: Any) -> Any:
    """Convert one raw TOML value; raises ValueError with a message about the value."""
    kind = spec.kind
    if kind == "str":
        if not isinstance(raw, str):
            raise ValueError(f"must be text, got {raw!r}")
        return raw
    if kind == "int":
        if isinstance(raw, bool) or not isinstance(raw, int):
            raise ValueError(f"must be a whole number, got {raw!r}")
        return raw
    if kind == "bool":
        if not isinstance(raw, bool):
            raise ValueError(f"must be true or false, got {raw!r}")
        return raw
    if kind == "number":
        if isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ValueError(f"must be a number without a unit, got {raw!r}")
        return float(raw)
    if kind == "nuclide":
        parse_nuclide(raw)
        return raw
    if kind == "material":
        parse_material(raw)
        return raw
    if kind == "direction":
        if (not isinstance(raw, list) or len(raw) != 3
                or not all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in raw)):
            raise ValueError(f"must be three numbers [x, y, z], got {raw!r}")
        if math.hypot(*raw) == 0:
            raise ValueError("must not be the zero vector")
        return [float(x) for x in raw]
    if kind in ("absorbers", "curve", "lines"):
        what = {"absorbers": '[["Pb", "1 mm"], ...]: a material and its thickness',
                "curve": '[["122 keV", "1.2 %"], ...]: a γ-ray energy and the efficiency there',
                "lines": '[["1274 keV", "0.5 /s"], ...]: a γ-ray energy and its rate'}[kind]
        if not isinstance(raw, list) or not all(isinstance(x, list) and len(x) == 2 for x in raw):
            raise ValueError(f"must be a list of pairs {what}; got {raw!r}")
        out = []
        for a, b in raw:
            if kind == "absorbers":
                if not isinstance(a, str):
                    raise ValueError(f"must name a material first in each pair, got {a!r}")
                parse_material(a)
                q = b if isinstance(b, Quantity) else Quantity.parse(b)
                if q.kind not in ("length", "areal_density") or q.value <= 0:
                    raise ValueError(f"'{q}' is not a thickness (a length, or mg/cm2)")
                out.append([a, q])
            else:
                e = a if isinstance(a, Quantity) else Quantity.parse(a)
                f = b if isinstance(b, Quantity) else Quantity.parse(b)
                if e.kind != "energy" or e.value <= 0:
                    raise ValueError(f"'{e}' is not a γ-ray energy")
                if kind == "curve" and (f.kind != "fraction" or not 0 < f.to("%") <= 100):
                    raise ValueError(f"'{f}' is not an efficiency (above 0 and at most 100 %)")
                if kind == "lines" and (f.kind != "rate" or f.value < 0):
                    raise ValueError(f"'{f}' is not a rate (counts per second, e.g. '0.5 /s')")
                out.append([e, f])
        if kind == "curve" and len(out) < 2:
            raise ValueError("needs at least two points")
        return out
    if kind == "lengths":
        if not isinstance(raw, list):
            raise ValueError(f'must be a list of lengths ["10 um", "20 um", ...], got {raw!r}')
        out = []
        for i, x in enumerate(raw, start=1):
            if not isinstance(x, (str, Quantity)):
                raise ValueError(f'must be lengths written with a unit, e.g. "10 um"; number {i} is {x!r}')
            q = x if isinstance(x, Quantity) else Quantity.parse(x)
            if q.kind != "length":
                raise ValueError(f"number {i}, '{q}', is not a length")
            if q.value <= 0:
                raise ValueError(f"must all be positive; number {i} is '{q}'")
            out.append(q)
        return out
    if kind == "position":
        if not isinstance(raw, list) or len(raw) != 3:
            raise ValueError(f"must be three lengths [\"x mm\", \"y mm\", \"z mm\"], got {raw!r}")
        qs = [raw_q if isinstance(raw_q, Quantity) else Quantity.parse(raw_q) for raw_q in raw]
        for q in qs:
            if q.kind != "length":
                raise ValueError(f"'{q}' is not a length")
        if all(q.value == 0 for q in qs):
            raise ValueError("must not be at the target (origin)")
        return qs
    # A quantity of one of the listed kinds.
    kinds = kind.split("|")
    q = raw if isinstance(raw, Quantity) else None
    if isinstance(raw, str):
        try:
            q = Quantity.parse(raw)
        except ValueError as err:  # say it with a unit that suits this field, not a generic length
            units = ", ".join(u for k in kinds for u in _units_of(k))
            if "has no unit" in str(err):
                raise ValueError(f"'{raw.strip()}' has no unit; write it with one, e.g. "
                                 f"'{raw.strip()} {_example_unit(kind)}'") from None
            if "unknown unit" in str(err):
                raise ValueError(f"{err}; use one of {units}") from None
            raise ValueError(f"cannot read '{raw}' as a number with a unit, e.g. '1 {_example_unit(kind)}'") from None
    if q is None:
        if isinstance(raw, (int, float)) and not isinstance(raw, bool):
            raise ValueError(f"{raw!r} has no unit; write it as text with a unit, e.g. \"{raw} {_example_unit(kind)}\"")
        raise ValueError(f"must be a number with a unit, e.g. \"1 {_example_unit(kind)}\", got {raw!r}")
    if q.kind not in kinds:
        units = ", ".join(u for k in kinds for u in _units_of(k))
        raise ValueError(f"'{q}' has the wrong unit; use one of {units}")
    return q


def _units_of(kind: str) -> list[str]:
    from .quantity import UNITS

    return list(UNITS[kind])


def _example_unit(kind: str) -> str:
    return _units_of(kind.split("|")[0])[0]


def _read_section(raw: Any, specs: tuple[_Field, ...], where: str, problems: list[str],
                  extra: tuple[str, ...] = ()) -> dict[str, Any]:
    """Read a table against its field specs. Problems are appended; returns the values that were valid."""
    if not isinstance(raw, dict):
        problems.append(f"{where}: must be a table of fields")
        return {}
    names = [s.name for s in specs]
    out: dict[str, Any] = {}
    for key in raw:
        if key not in names and key not in extra:
            hint = difflib.get_close_matches(key, names, n=1)
            problems.append(f"{where}: unknown field '{key}'" + (f"; did you mean '{hint[0]}'?" if hint else ""))
    for spec in specs:
        if spec.name not in raw:
            if spec.required:
                problems.append(f"{where}: '{spec.name}' is missing")
            continue
        try:
            value = _read_value(spec, raw[spec.name])
        except ValueError as e:
            problems.append(f"{where}: {spec.name} {_lower_first(str(e))}")
            continue
        if spec.check is not None:
            message = spec.check(value)
            if message:
                problems.append(f"{where}: {spec.name} {message}, got '{value}'")
                continue
        out[spec.name] = value
    return out


def _lower_first(s: str) -> str:
    """Lower-case a leading capitalised word ("Cannot ..."), but not a symbol or name ("Tc has ...")."""
    word = s.split(" ", 1)[0]
    return s[:1].lower() + s[1:] if len(word) > 2 and word[1:].isalpha() and word[1:].islower() else s


def _q(x: Optional[QuantityLike]) -> Optional[Quantity]:
    return Quantity.parse(x) if isinstance(x, str) else x


def _out(v: Any) -> Any:
    if isinstance(v, Quantity):
        return str(v)
    if isinstance(v, list):
        return [_out(x) for x in v]
    return v


def _to_dict(obj: Any) -> dict[str, Any]:
    """Fields in the order of the section's specs (the order a setup file is written in), omitting unset ones."""
    names = [s.name for s in obj.SPECS] + [f.name for f in fields(obj) if f.name not in {s.name for s in obj.SPECS}]
    return {n: _out(getattr(obj, n)) for n in names
            if getattr(obj, n) is not None and n != "backing" and not n.startswith("_")}


# ---------------------------------------------------------------------------------------------------------------
# Sections


@dataclass
class Beam:
    """The beam: which nuclide, how fast, how intense."""

    nuclide: str
    energy: QuantityLike
    current: QuantityLike
    charge_state: Optional[int] = None
    energy_spread: Optional[QuantityLike] = None
    spot_size: Optional[QuantityLike] = None
    #: A beam halo: this share of the beam (``"2 %"``) arrives spread evenly over a disc of ``halo_radius`` (the
    #: target frame's aperture) instead of the spot.
    halo_fraction: Optional[QuantityLike] = None
    halo_radius: Optional[QuantityLike] = None

    SPECS = (
        _Field("nuclide", "nuclide", required=True),
        _Field("energy", "energy|energy_per_nucleon", required=True, check=_positive),
        _Field("current", "particle_current|electrical_current", required=True, check=_positive),
        _Field("charge_state", "int"),
        _Field("energy_spread", "energy|fraction", check=_non_negative),
        _Field("spot_size", "length", check=_non_negative),
        _Field("halo_fraction", "fraction", check=lambda q: None if 0 <= q.to("%") <= 100 else
               "must be between 0 and 100 %"),
        _Field("halo_radius", "length", check=_non_negative),
    )

    @property
    def Z(self) -> int:
        return parse_nuclide(self.nuclide)[0]

    @property
    def A(self) -> int:
        return parse_nuclide(self.nuclide)[1]

    @property
    def charge(self) -> int:
        """Charge state; fully stripped (Z) unless given."""
        return self.Z if self.charge_state is None else self.charge_state

    @property
    def energy_mev(self) -> float:
        """Total kinetic energy, MeV (``MeV/u`` is multiplied by the mass number A)."""
        q = _q(self.energy)
        return q.canonical * self.A if q.kind == "energy_per_nucleon" else q.canonical

    @property
    def particle_current_pna(self) -> float:
        """Beam particles per second expressed as particle-nA (electrical current divided by the charge state)."""
        q = _q(self.current)
        return q.canonical / self.charge if q.kind == "electrical_current" else q.canonical

    @property
    def particles_per_second(self) -> float:
        return self.particle_current_pna * 1e-9 / 1.602176634e-19


@dataclass
class Layer:
    """A layer of material, such as a target backing."""

    material: str
    thickness: QuantityLike
    density: Optional[QuantityLike] = None

    SPECS = (
        _Field("material", "material", required=True),
        _Field("thickness", "areal_density|length", required=True, check=_positive),
        _Field("density", "density", check=_positive),
    )

    def material_data(self) -> _data.Material:
        """Composition and density of the layer (see :func:`physim.nuclear.data.material`)."""
        return _data.material(self.material, _q(self.density))


@dataclass
class Target:
    """The target, optionally on a backing. ``tilt`` turns it about the vertical (y) axis."""

    material: str
    thickness: QuantityLike
    tilt: Optional[QuantityLike] = None
    density: Optional[QuantityLike] = None
    backing: Optional[Layer] = None
    #: Where the target sits along the beam, from the chamber's centre (0 if left out). Detectors are placed
    #: from the chamber's centre; every angle and distance in the results is from the target.
    position: Optional[QuantityLike] = None
    #: A target ladder: pairs of a material and a thickness; ``selected`` (from 1) is the one in the beam, and
    #: ``material`` and ``thickness`` are then its values.
    ladder: Optional[list] = None
    selected: Optional[int] = None
    #: Thin layers of contaminants on the target's downstream face (behind the backing, if any): pairs of a material
    #: and a thickness, ``[["C", "10 ug/cm2"], ["O", "3 ug/cm2"]]``. Each is a layer of its own, with its own
    #: elastic scattering.
    contaminants: Optional[list] = None

    SPECS = (
        _Field("material", "material"),
        _Field("thickness", "areal_density|length", check=_positive),
        _Field("density", "density", check=_positive),
        _Field("tilt", "angle", check=_tilt),
        _Field("position", "length"),
        _Field("ladder", "absorbers"),
        _Field("selected", "int", check=lambda n: None if n >= 1 else "counts from 1"),
        _Field("contaminants", "absorbers"),
    )

    @property
    def position_mm(self) -> float:
        return _q(self.position).to("mm") if self.position is not None else 0.0

    def material_data(self) -> _data.Material:
        """Composition and density of the target (see :func:`physim.nuclear.data.material`)."""
        return _data.material(self.material, _q(self.density))


@dataclass
class Detector:
    """One detector. Place it with ``theta``/``phi``/``distance`` or with ``position``; size it by ``shape``."""

    shape: str
    thickness: QuantityLike
    name: Optional[str] = None
    #: A model of the detector catalogue (:mod:`physim.nuclear.catalogue`), whose values fill the fields not given.
    model: Optional[str] = None
    theta: Optional[QuantityLike] = None
    phi: Optional[QuantityLike] = None
    distance: Optional[QuantityLike] = None
    position: Optional[list] = None
    facing: Optional[list] = None
    rotation: Optional[QuantityLike] = None
    width: Optional[QuantityLike] = None
    height: Optional[QuantityLike] = None
    strips_x: Optional[int] = None
    strips_y: Optional[int] = None
    inner_radius: Optional[QuantityLike] = None
    outer_radius: Optional[QuantityLike] = None
    rings: Optional[int] = None
    sectors: Optional[int] = None
    radius: Optional[QuantityLike] = None
    material: Optional[str] = None
    density: Optional[QuantityLike] = None
    dead_layer: Optional[QuantityLike] = None
    resolution: Optional[QuantityLike] = None
    threshold: Optional[QuantityLike] = None
    #: The target's place from the chamber's centre, mm (set by the setup; positions in the file are from the
    #: chamber's centre, results are from the target).
    _origin: Optional[tuple] = None

    SPECS = (
        _Field("name", "str"),
        _Field("model", "str"),
        _Field("shape", "str", required=True),
        _Field("theta", "angle", check=_polar),
        _Field("phi", "angle"),
        _Field("distance", "length", check=_positive),
        _Field("position", "position"),
        _Field("facing", "direction"),
        _Field("rotation", "angle"),
        _Field("width", "length", check=_positive),
        _Field("height", "length", check=_positive),
        _Field("strips_x", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("strips_y", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("inner_radius", "length", check=_non_negative),
        _Field("outer_radius", "length", check=_positive),
        _Field("rings", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("sectors", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("radius", "length", check=_positive),
        _Field("material", "material"),
        _Field("density", "density", check=_positive),
        _Field("thickness", "length|areal_density", required=True, check=_positive),
        _Field("dead_layer", "length|areal_density", check=_non_negative),
        _Field("resolution", "energy", check=_non_negative),
        _Field("threshold", "energy", check=_non_negative),
    )

    @property
    def detector_material(self) -> str:
        """Detector material; silicon unless given."""
        return self.material or "Si"

    def material_data(self) -> _data.Material:
        """Composition and density of the detector (and its dead layer)."""
        return _data.material(self.detector_material, _q(self.density))

    def blocking(self) -> list:
        """Dead material of the detector's model that stops particles: annuli in the plane of the face, as
        (name, inner radius, outer radius) in mm. Empty without a model."""
        if self.model is None:
            return []
        return [(b["name"], _q(b["inner_radius"]).to("mm"), _q(b["outer_radius"]).to("mm"))
                for b in catalogue.model(self.model, "particle").blocking]

    def chamber_position_mm(self) -> tuple[float, float, float]:
        """Centre of the detector face, mm, from the chamber's centre (as the setup file places it)."""
        if self.position is not None:
            x, y, z = (_q(c).to("mm") for c in self.position)
            return (x, y, z)
        th = math.radians(_q(self.theta).to("deg"))
        ph = math.radians(_q(self.phi).to("deg")) if self.phi is not None else 0.0
        d = _q(self.distance).to("mm")
        return (d * math.sin(th) * math.cos(ph), d * math.sin(th) * math.sin(ph), d * math.cos(th))

    def position_mm(self) -> tuple[float, float, float]:
        """Centre of the detector face, mm, from the target (the coordinates of the module docstring)."""
        x, y, z = self.chamber_position_mm()
        if self._origin:
            return (x - self._origin[0], y - self._origin[1], z - self._origin[2])
        return (x, y, z)

    def facing_direction(self) -> tuple[float, float, float]:
        """Unit vector the sensitive face looks along; towards the target unless ``facing`` is given."""
        v = self.facing if self.facing is not None else [-c for c in self.position_mm()]
        n = math.hypot(*v)
        return (v[0] / n, v[1] / n, v[2] / n)


@dataclass
class Run:
    """How long the beam runs, and how many counts per detector the experiment needs."""

    beam_time: QuantityLike
    counts_wanted: Optional[int] = None
    #: Full width of the coincidence window (100 ns if left out), and a non-paralysable dead time per count.
    coincidence_window: Optional[QuantityLike] = None
    dead_time: Optional[QuantityLike] = None
    #: Room background in each γ-ray crystal, counts per second in its lines (⁴⁰K, thorium and uranium series).
    room_background: Optional[QuantityLike] = None
    #: Extra γ-ray lines in every crystal, by hand: pairs of an energy and a rate (``[["1274 keV", "0.5 /s"]]``).
    extra_lines: Optional[list] = None
    #: Shaping time of the γ-ray amplifiers: two signals closer than twice it pile up (none if left out).
    shaping_time: Optional[QuantityLike] = None
    #: Fusion-evaporation: the exit channel the experiment studies, by its particles ("4n") or its residue
    #: ("158Dy"). Its strongest line sets the beam time; the largest channel when left out.
    channel: Optional[str] = None

    SPECS = (
        _Field("beam_time", "time", required=True, check=_positive),
        _Field("counts_wanted", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("coincidence_window", "time", check=_positive),
        _Field("dead_time", "time", check=_non_negative),
        _Field("room_background", "rate", check=_non_negative),
        _Field("extra_lines", "lines"),
        _Field("shaping_time", "time", check=_non_negative),
        _Field("channel", "str"),
    )


@dataclass
class Excitation:
    """The excited state of a Coulomb-excitation reaction (``[reaction] type = "coulex"``).

    One state, reached from a 0⁺ ground state: its energy, multipolarity (E1, E2, E3), the reduced transition
    probability ``b_up`` = B(Eλ↑) as text with a unit (``"0.0695 e2b2"``), and which nucleus is excited.

    With a level scheme of the excited nucleus (``[levels]``) the simulated events populate its levels as the
    coupled equations give them, to all orders; ``order = "first"`` keeps first-order theory (each level in one
    step from the ground state).

    ``e1_polarisation`` is GOSIA's parameter d of the E1 polarisation correction to the excitation
    (:mod:`physim.nuclear.coulex`), not GOSIA's input ``DIP,``, which is 1000 d. Left out it is 0.0038, GOSIA's own
    default; 0 leaves the correction out.
    """

    energy: QuantityLike
    b_up: str
    excite: str = "target"
    multipolarity: str = "E2"
    #: How the γ rays leave the nucleus: "correlated" with the scattered particle (also when left out), or
    #: "isotropic".
    emission: Optional[str] = None
    #: How the levels of a level scheme are populated in the events: "all" orders (also when left out) or
    #: "first" order.
    order: Optional[str] = None
    #: GOSIA's E1 polarisation parameter d; None: 0.0038, GOSIA's default; 0: no correction. The value in effect is
    #: :func:`physim.nuclear.gamma.e1_polarisation`.
    e1_polarisation: Optional[float] = None

    SPECS = (
        _Field("excite", "str", check=lambda s: None if s in ("target", "projectile")
               else "must be 'target' or 'projectile'"),
        _Field("energy", "energy", check=_positive),
        _Field("multipolarity", "str", check=lambda s: None if s in ("E1", "E2", "E3") else "must be E1, E2 or E3"),
        _Field("b_up", "str"),
        _Field("emission", "str", check=lambda s: None if s in ("correlated", "isotropic")
               else "must be 'correlated' or 'isotropic'"),
        _Field("order", "str", check=lambda s: None if s in ("all", "first") else "must be 'all' or 'first'"),
        _Field("e1_polarisation", "number", check=lambda d: None if 0 <= d <= 0.05 else
               "must be between 0 and 0.05: GOSIA's parameter d (0.0038 by default, 0 for none), not its DIP, input "
               "(1000 d)"),
    )

    @property
    def energy_mev(self) -> float:
        return _q(self.energy).to("MeV")

    @property
    def lam(self) -> int:
        return int(self.multipolarity[1])

    @property
    def b_up_e2fm(self) -> float:
        """B(Eλ↑) in e² fm^(2λ)."""
        from .coulex import parse_b

        return parse_b(self.b_up, self.lam)


@dataclass
class ExitChannel:
    """One exit channel of a fusion-evaporation reaction (``[[reaction.channels]]``): the light particles the
    compound nucleus boils off (``"4n"``, ``"p2n"``, ``"a n"``), the residue they leave (worked out from the
    compound nucleus if left out) and the channel's cross section with its unit (``"450 mb"``).

    ``temperature`` sets the evaporation spectrum of the channel's particles (the reaction's, or one from the
    excitation energy, if left out). ``entry_energy`` (excitation energy of the residue, MeV) and ``entry_spin``
    (ħ) are where its γ-ray cascade starts. ``source`` says where the cross section comes from, in words
    ("PACE4", "measured", "illustrative"); "placeholder" marks a value nobody has checked, and the planner warns.
    """

    particles: str
    cross_section: QuantityLike
    residue: Optional[str] = None
    temperature: Optional[QuantityLike] = None
    entry_energy: Optional[QuantityLike] = None
    entry_spin: Optional[float] = None
    source: Optional[str] = None
    #: The residue's name as the setup's compound nucleus gives it (set when the setup is read).
    _residue: Optional[str] = None

    SPECS = (
        _Field("residue", "nuclide"),
        _Field("particles", "str", required=True),
        _Field("cross_section", "cross_section", required=True, check=_positive),
        _Field("temperature", "energy", check=_positive),
        _Field("entry_energy", "energy", check=_non_negative),
        _Field("entry_spin", "number", check=lambda x: None if 0 <= x <= 120 and (2 * x).is_integer()
               else "must be a whole or half-whole number of ħ between 0 and 120"),
        _Field("source", "str"),
    )

    @property
    def cross_section_mb(self) -> float:
        return _q(self.cross_section).to("mb")

    @property
    def residue_name(self) -> str:
        """The residue, as AME names it ("158Er"): from the setup's compound nucleus, or as given."""
        if self._residue is not None:
            return self._residue
        return _data.nuclide(self.residue).name if self.residue else "?"


@dataclass
class FusionEvaporation:
    """The fusion-evaporation reaction (``[reaction] type = "fusion-evaporation"``): the beam fuses with the
    target's main nuclide into a compound nucleus, which boils off light particles; each exit channel of
    ``channels`` has its cross section.

    ``temperature`` is the default nuclear temperature of the evaporation spectra. ``alignment`` is the width
    σ/J of the compound nucleus's distribution of magnetic substates about the beam (its spin is perpendicular to
    the beam; 0.3 if left out). ``side_feeding`` is the share of the γ-ray cascades that enter the known levels
    below the top one (25 % if left out).

    For lifetime measurements, every cascade reaches its first known level through a rotational side-feeding band
    (the LINESHAPE model): ``side_feeding_levels`` E2 transitions (a whole number, 5 if left out; 0 makes the
    feeding prompt) with the quadrupole moment ``side_feeding_q`` (in barns, "6 b" if left out). A longer band, or
    a smaller quadrupole moment, makes the side feeding slower.
    """

    channels: list = field(default_factory=list)
    temperature: Optional[QuantityLike] = None
    alignment: Optional[float] = None
    side_feeding: Optional[QuantityLike] = None
    side_feeding_q: Optional[QuantityLike] = None
    side_feeding_levels: Optional[int] = None

    SPECS = (
        _Field("temperature", "energy", check=_positive),
        _Field("alignment", "number", check=lambda x: None if 0 < x <= 2 else "must be above 0 and at most 2"),
        _Field("side_feeding", "fraction", check=lambda q: None if 0 <= q.to("%") <= 100
               else "must be between 0 and 100 %"),
        # A quadrupole moment is an area: barns, the units of a cross section.
        _Field("side_feeding_q", "cross_section", check=lambda q: None if 0 < q.to("b") <= 100
               else "must be a quadrupole moment above 0 and at most 100 b, e.g. \"6 b\""),
        _Field("side_feeding_levels", "int", check=lambda x: None if 0 <= x <= 50
               else "must be a whole number of transitions from 0 to 50"),
    )


@dataclass
class GammaDetector:
    """A γ-ray detector facing the target: one crystal, or the four crystals of a clover.

    ``distance`` is from the target to the front face of the crystals. Give ``radius`` for a single disc, or a
    ``model`` of the detector catalogue (or the crystal fields themselves) for real crystals: ``crystals`` (1 or 4,
    four being a 2 × 2 clover with centres ``crystal_pitch`` apart), ``crystal_diameter`` and ``crystal_length``,
    in a square housing of side ``housing_side`` whose window is ``window_gap`` in front of the crystals.
    """

    theta: QuantityLike
    distance: QuantityLike
    radius: Optional[QuantityLike] = None
    name: Optional[str] = None
    model: Optional[str] = None
    phi: Optional[QuantityLike] = None
    crystals: Optional[int] = None
    crystal_diameter: Optional[QuantityLike] = None
    crystal_length: Optional[QuantityLike] = None
    crystal_pitch: Optional[QuantityLike] = None
    housing_side: Optional[QuantityLike] = None
    window_gap: Optional[QuantityLike] = None
    #: Resolution (FWHM) at 1332 keV; it varies with energy as :mod:`physim.nuclear.response` describes.
    resolution: Optional[QuantityLike] = None
    #: Full-energy-peak efficiency for the γ ray, as a percentage of all γ rays emitted ("2.5 %"), used at every
    #: energy. Without it the efficiency comes from ``efficiency_curve`` or from the response model.
    efficiency: Optional[QuantityLike] = None
    #: The crystal's material: "Ge" (the default) or "LaBr3".
    material: Optional[str] = None
    #: Material between the target and the detector: pairs of a material and a thickness.
    absorbers: Optional[list] = None
    #: A measured full-energy-peak efficiency: pairs of a γ-ray energy and the efficiency there, for the whole
    #: detector where it stands. It replaces the model's efficiency.
    efficiency_curve: Optional[list] = None
    #: Energy below which the crystal records nothing.
    threshold: Optional[QuantityLike] = None
    #: Add-back for a clover: the energies its crystals record together are summed, so a γ ray that scatters from
    #: one crystal into the next still ends in the full-energy peak.
    addback: Optional[bool] = None
    #: The add-back factor at 1332 keV: the full-energy peak with add-back over the peak without
    #: (:data:`physim.nuclear.response.ADDBACK_FACTOR` if left out).
    addback_factor: Optional[float] = None
    #: A Compton-suppression shield around the crystals ("BGO"): a γ ray that scatters out of the crystals is
    #: seen by the shield and the event is rejected. The shield stops particles, as the housing does.
    shield: Optional[str] = None
    shield_thickness: Optional[QuantityLike] = None
    #: By how much the shield lowers the continuum at 1332 keV; the peak is unchanged
    #: (:data:`physim.nuclear.response.SUPPRESSION_FACTOR` if left out).
    suppression_factor: Optional[float] = None
    #: The target's place from the chamber's centre, mm (set by the setup).
    _origin: Optional[tuple] = None

    SPECS = (
        _Field("name", "str"),
        _Field("model", "str"),
        _Field("theta", "angle", required=True, check=_polar),
        _Field("phi", "angle"),
        _Field("distance", "length", required=True, check=_positive),
        _Field("radius", "length", check=_positive),
        _Field("crystals", "int", check=lambda n: None if n in (1, 4) else "must be 1 or 4"),
        _Field("crystal_diameter", "length", check=_positive),
        _Field("crystal_length", "length", check=_positive),
        _Field("crystal_pitch", "length", check=_positive),
        _Field("housing_side", "length", check=_positive),
        _Field("window_gap", "length", check=_non_negative),
        _Field("resolution", "energy", check=_non_negative),
        _Field("efficiency", "fraction", check=lambda q: None if 0 < q.value <= 100
               else f"must be above 0 and at most 100 %, got '{q}'"),
        _Field("material", "str"),
        _Field("absorbers", "absorbers"),
        _Field("efficiency_curve", "curve"),
        _Field("threshold", "energy", check=_non_negative),
        _Field("addback", "bool"),
        _Field("addback_factor", "number", check=lambda x: None if x >= 1 else "must be at least 1"),
        _Field("shield", "str"),
        _Field("shield_thickness", "length", check=_positive),
        _Field("suppression_factor", "number", check=lambda x: None if x >= 1 else "must be at least 1"),
    )

    def geometric_efficiency(self) -> float:
        """Fraction of all directions the crystal faces cover, seen from the target: Σ (1 − cos α)/2 over the
        crystals, α being the half-angle of each."""
        return sum((1.0 - math.cos(math.atan2(r, math.hypot(*c)))) / 2.0 for _, c, r in self.elements())

    def crystal_radius_mm(self) -> float:
        """Radius of one crystal's front face, mm."""
        if self.crystal_diameter is not None:
            return _q(self.crystal_diameter).to("mm") / 2.0
        return _q(self.radius).to("mm")

    def face_axes(self) -> tuple:
        """Two unit vectors across the front face, perpendicular to :meth:`direction`."""
        u = self.direction()
        ref = (0.0, 0.0, 1.0) if abs(u[2]) < 0.9 else (1.0, 0.0, 0.0)
        a = (u[1] * ref[2] - u[2] * ref[1], u[2] * ref[0] - u[0] * ref[2], u[0] * ref[1] - u[1] * ref[0])
        n = math.hypot(*a)
        a = (a[0] / n, a[1] / n, a[2] / n)
        b = (u[1] * a[2] - u[2] * a[1], u[2] * a[0] - u[0] * a[2], u[0] * a[1] - u[1] * a[0])
        return a, b

    def elements(self) -> list:
        """The crystals, as (label, centre of the front face in mm, radius in mm). One crystal has an empty label;
        the four of a clover are A to D."""
        u, d, r = self.direction(), self.distance_mm(), self.crystal_radius_mm()
        centre = (d * u[0], d * u[1], d * u[2])
        if (self.crystals or 1) == 1:
            return [("", centre, r)]
        a, b = self.face_axes()
        h = _q(self.crystal_pitch).to("mm") / 2.0
        return [(label, tuple(centre[k] + sa * h * a[k] + sb * h * b[k] for k in range(3)), r)
                for label, sa, sb in (("A", 1, 1), ("B", -1, 1), ("C", -1, -1), ("D", 1, -1))]

    def shield_mm(self) -> float:
        """Thickness of the shield's wall, mm (the typical value of :data:`physim.nuclear.response.SHIELDS` when
        the setup gives none); 0 without a shield."""
        if self.shield is None:
            return 0.0
        if self.shield_thickness is not None:
            return _q(self.shield_thickness).to("mm")
        from .response import SHIELDS

        return SHIELDS.get(self.shield, 25.0)

    def housing(self) -> Optional[tuple]:
        """The housing's front window, which stops particles: (centre in mm, side in mm) of a square facing the
        target, ``window_gap`` in front of the crystals; with a shield, the shield's outer side. ``None`` when no
        housing is given."""
        if self.housing_side is None:
            return None
        u = self.direction()
        d = self.distance_mm() - (_q(self.window_gap).to("mm") if self.window_gap is not None else 0.0)
        return (d * u[0], d * u[1], d * u[2]), _q(self.housing_side).to("mm") + 2 * self.shield_mm()

    def peak_efficiency(self, energy_mev: float = 1.332492, experiment=None) -> float:
        """Full-energy-peak efficiency at a γ-ray energy, as a fraction of all γ rays emitted at the target:
        :attr:`efficiency` if given, else :attr:`efficiency_curve`, else the response model
        (:class:`physim.nuclear.response.Response`). ``experiment`` adds the chamber wall's attenuation."""
        from .response import Response

        return float(Response(experiment, self).peak_efficiency(energy_mev))

    def chamber_centre_mm(self) -> tuple:
        """Centre of the front face from the chamber's centre, mm, as the setup file places it."""
        th = math.radians(_q(self.theta).to("deg"))
        ph = math.radians(_q(self.phi).to("deg")) if self.phi is not None else 0.0
        d = _q(self.distance).to("mm")
        return (d * math.sin(th) * math.cos(ph), d * math.sin(th) * math.sin(ph), d * math.cos(th))

    def centre_mm(self) -> tuple:
        """Centre of the front face from the target, mm."""
        x, y, z = self.chamber_centre_mm()
        if self._origin:
            return (x - self._origin[0], y - self._origin[1], z - self._origin[2])
        return (x, y, z)

    def distance_mm(self) -> float:
        """Distance from the target to the front face, mm."""
        return math.hypot(*self.centre_mm())

    def direction(self) -> tuple:
        """Unit vector from the target to the centre of the detector."""
        c = self.centre_mm()
        d = math.hypot(*c)
        return (c[0] / d, c[1] / d, c[2] / d)

    def half_angle_deg(self) -> float:
        """Half the opening angle of the detector seen from the target: of the one crystal, or of the circle
        around the four crystals of a clover."""
        r = self.crystal_radius_mm()
        if (self.crystals or 1) == 4:
            r += _q(self.crystal_pitch).to("mm") / math.sqrt(2.0)
        return math.degrees(math.atan2(r, self.distance_mm()))


@dataclass
class Chamber:
    """The scattering chamber: a sphere of ``radius`` around the target with a wall, and the beam pipe. Particle
    detectors sit inside it and γ-ray detectors outside."""

    radius: QuantityLike
    wall_thickness: Optional[QuantityLike] = None
    wall_material: Optional[str] = None
    beam_pipe_radius: Optional[QuantityLike] = None

    SPECS = (
        _Field("radius", "length", required=True, check=_positive),
        _Field("wall_thickness", "length", check=_positive),
        _Field("wall_material", "material"),
        _Field("beam_pipe_radius", "length", check=_positive),
    )


def _scintillator(name: str) -> Optional[str]:
    from .neutrons import materials

    known = materials()
    return None if name in known else f"is not a known scintillator (known: {', '.join(sorted(known))})"


@dataclass
class NeutronWallSetup:
    """A wall of organic-scintillator cells that counts the evaporated neutrons (``[neutron_wall]``, for a
    fusion-evaporation reaction): ``cells`` identical cylinders ``distance`` from the target, side by side at
    ``pitch``, filling a cone around the beam from ``theta_min``. The fields left out take the typical values of
    :class:`physim.nuclear.neutrons.NeutronDetector` (an EJ-301 cell, 127 mm across and 150 mm deep, 100 keVee).

    Three fields are for the channel selection (:func:`physim.nuclear.exit_channels.channel_selection`):
    ``gamma_efficiency``, the share of γ rays entering a cell that give light above its threshold (typical 40 %);
    ``gamma_multiplicity``, the γ rays per event (default: J/2 from the entry spin J plus
    :data:`~physim.nuclear.exit_channels.STATISTICAL_GAMMAS` statistical ones); and ``particle_threshold``, the
    energy a proton or α particle must deposit in a particle detector that sets no ``threshold`` of its own
    (default :data:`~physim.nuclear.exit_channels.PARTICLE_THRESHOLD_MEV`)."""

    cells: int
    distance: QuantityLike
    pitch: Optional[QuantityLike] = None
    radius: Optional[QuantityLike] = None
    depth: Optional[QuantityLike] = None
    material: Optional[str] = None
    threshold: Optional[QuantityLike] = None
    layout: Optional[str] = None
    theta_min: Optional[QuantityLike] = None
    neutron_as_gamma: Optional[QuantityLike] = None
    gamma_as_neutron: Optional[QuantityLike] = None
    gamma_efficiency: Optional[QuantityLike] = None
    gamma_multiplicity: Optional[float] = None
    particle_threshold: Optional[QuantityLike] = None

    SPECS = (
        _Field("cells", "int", required=True, check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("distance", "length", required=True, check=_positive),
        _Field("pitch", "length", check=_positive),
        _Field("radius", "length", check=_positive),
        _Field("depth", "length", check=_positive),
        _Field("material", "str", check=_scintillator),
        _Field("threshold", "energy", check=_positive),
        _Field("layout", "str", check=lambda s: None if s in ("hexagonal", "square")
               else "must be 'hexagonal' or 'square'"),
        _Field("theta_min", "angle", check=lambda q: None if 0 <= q.to("deg") < 180 else "must be 0 to 180 deg"),
        _Field("neutron_as_gamma", "fraction", check=lambda q: None if 0 <= q.to("%") <= 100 else
               "must be 0 to 100 %"),
        _Field("gamma_as_neutron", "fraction", check=lambda q: None if 0 <= q.to("%") <= 100 else
               "must be 0 to 100 %"),
        _Field("gamma_efficiency", "fraction", check=lambda q: None if 0 <= q.to("%") <= 100 else
               "must be 0 to 100 %"),
        _Field("gamma_multiplicity", "number", check=lambda x: None if 0 <= x <= 100 else "must be 0 to 100"),
        _Field("particle_threshold", "energy", check=_non_negative),
    )

    def model(self):
        """The wall as :class:`physim.nuclear.neutrons.NeutronWall`, with the typical values where the setup gives
        none."""
        from . import neutrons

        def text(q, unit: str) -> str:
            return f"{_q(q).to(unit):.10g} {unit}"

        cell = {"distance": text(self.distance, "mm")}
        for name, unit in (("radius", "mm"), ("depth", "mm"), ("threshold", "MeV"), ("neutron_as_gamma", "%"),
                           ("gamma_as_neutron", "%"), ("gamma_efficiency", "%")):
            if getattr(self, name) is not None:
                cell[name] = text(getattr(self, name), unit)
        if self.material is not None:
            cell["material"] = self.material
        wall = {"cell": neutrons.NeutronDetector(**cell), "cells": self.cells}
        if self.pitch is not None:
            wall["pitch"] = text(self.pitch, "mm")
        else:  # cells touching: the pitch is the cell's diameter
            wall["pitch"] = f"{2 * wall['cell'].radius_mm:.10g} mm"
        if self.layout is not None:
            wall["layout"] = self.layout
        if self.theta_min is not None:
            wall["theta_min"] = text(self.theta_min, "deg")
        return neutrons.NeutronWall(**wall)


def _one_of(*names: str) -> Callable[[str], Optional[str]]:
    return lambda s: None if s in names else "must be " + " or ".join(f"'{n}'" for n in names)


def _percent(q: Quantity) -> Optional[str]:
    return None if 0 <= q.to("%") <= 100 else "must be 0 to 100 %"


@dataclass
class ChargedParticleBallSetup:
    """A ball of particle detectors round the target that counts the evaporated protons and α particles
    (``[charged_particle_ball]``, for a fusion-evaporation reaction): ``elements`` detectors of ``material`` (Si or
    CsI) at ``radius`` from the target, covering θ from ``theta_min`` to ``theta_max`` (all of it if left out).
    Described for the scene and the plan; the channel selection does not use it yet."""

    elements: int
    radius: QuantityLike
    material: Optional[str] = None
    thickness: Optional[QuantityLike] = None
    threshold: Optional[QuantityLike] = None
    theta_min: Optional[QuantityLike] = None
    theta_max: Optional[QuantityLike] = None

    SPECS = (
        _Field("elements", "int", required=True, check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("radius", "length", required=True, check=_positive),
        _Field("material", "str", check=_one_of("Si", "CsI")),
        _Field("thickness", "length", check=_positive),
        _Field("threshold", "energy", check=_non_negative),
        _Field("theta_min", "angle", check=_polar),
        _Field("theta_max", "angle", check=_polar),
    )


@dataclass
class RecoilSeparatorSetup:
    """A recoil separator after the target (``[recoil_separator]``, for a fusion-evaporation reaction), described
    by what it accepts rather than by its magnets: residues within ``theta_acceptance`` of the beam (and within
    ``momentum_acceptance``) enter it, ``distance`` downstream of the target. ``kind`` is "gas-filled", "vacuum"
    or "acceptance". Described for the scene and the plan; its transmission is not computed yet."""

    theta_acceptance: QuantityLike
    kind: Optional[str] = None
    distance: Optional[QuantityLike] = None
    momentum_acceptance: Optional[QuantityLike] = None
    beam_suppression: Optional[float] = None
    focal_plane: Optional[str] = None

    SPECS = (
        _Field("theta_acceptance", "angle", required=True,
               check=lambda q: None if 0 < q.to("deg") < 90 else "must be between 0 and 90 deg"),
        _Field("kind", "str", check=_one_of("gas-filled", "vacuum", "acceptance")),
        _Field("distance", "length", check=_positive),
        _Field("momentum_acceptance", "fraction", check=_percent),
        _Field("beam_suppression", "number", check=lambda x: None if x >= 1 else "must be at least 1"),
        _Field("focal_plane", "str"),
    )


@dataclass
class GammaCalorimeterSetup:
    """A shell of scintillators round the target that counts how many γ rays an event gives
    (``[gamma_calorimeter]``, for a fusion-evaporation reaction): ``elements`` detectors of ``material`` (BGO or
    NaI) from ``radius`` outwards, ``thickness`` deep; events with at least ``fold_min`` of them firing are kept.
    Described for the scene and the plan; the fold selection is not computed yet."""

    elements: int
    radius: QuantityLike
    material: Optional[str] = None
    thickness: Optional[QuantityLike] = None
    efficiency: Optional[QuantityLike] = None
    fold_min: Optional[int] = None

    SPECS = (
        _Field("elements", "int", required=True, check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("radius", "length", required=True, check=_positive),
        _Field("material", "str", check=_one_of("BGO", "NaI")),
        _Field("thickness", "length", check=_positive),
        _Field("efficiency", "fraction", check=_percent),
        _Field("fold_min", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
    )


#: The ancillary detectors of a fusion-evaporation setup beside the neutron wall: (section, class, what it does).
ANCILLARIES = (("charged_particle_ball", ChargedParticleBallSetup, "a charged-particle ball counts evaporated protons "
                "and α particles"),
               ("recoil_separator", RecoilSeparatorSetup, "a recoil separator collects the residues"),
               ("gamma_calorimeter", GammaCalorimeterSetup, "a γ calorimeter counts the γ rays of fusion events"))


#: Methods of a lifetime measurement: "rdds" (recoil distance, a plunger) and "dsam" (Doppler-shift attenuation in
#: the target's backing, or in a thick target).
LIFETIME_METHODS = ("rdds", "dsam")

#: Kinds of a plunger's foil: a "stopper" that stops the recoil, a "degrader" it flies through (slower), or the
#: first foil of a "differential" plunger whose second foil (``[lifetime.degrader2]``) stops it.
STOPPER_KINDS = ("stopper", "degrader", "differential")


@dataclass
class Stopper(Layer):
    """The plunger's foil after the gap (``[lifetime.stopper]``): a layer of material and its ``kind``
    (:data:`STOPPER_KINDS`)."""

    kind: str = "stopper"

    SPECS = (_Field("kind", "str", required=True, check=_one_of(*STOPPER_KINDS)),) + Layer.SPECS


@dataclass
class Degrader2(Layer):
    """The second foil of a differential plunger (``[lifetime.degrader2]``), ``gap`` behind the first; it stops
    the recoil."""

    gap: Optional[QuantityLike] = None

    SPECS = Layer.SPECS + (_Field("gap", "length", required=True, check=_positive),)

    @property
    def gap_um(self) -> float:
        return _q(self.gap).to("um")


def _alpha(x: float) -> Optional[str]:
    return None if 0 <= x <= 1 else "must be between 0 and 1"


@dataclass
class Deorientation:
    """Deorientation in vacuum (``lifetime.deorientation``): the recoil's alignment fades while it flies through
    the plunger's gap, G_k(t) = α_k + (1 − α_k) e^(−Γ_k t) for k = 2, 4 (the static model of Radeck et al.,
    PRC 85, 014301). Give ``alpha2``, ``alpha4`` and either the rates ``gamma2``, ``gamma4`` ("0.142 /ps") or the
    g factor ``g`` with ``c2``, ``c4`` ("3.3 ps"), so that Γ_k = |g| / C_k."""

    alpha2: Optional[float] = None
    gamma2: Optional[QuantityLike] = None
    alpha4: Optional[float] = None
    gamma4: Optional[QuantityLike] = None
    g: Optional[float] = None
    c2: Optional[QuantityLike] = None
    c4: Optional[QuantityLike] = None

    SPECS = (
        _Field("alpha2", "number", check=_alpha),
        _Field("gamma2", "rate", check=_non_negative),
        _Field("alpha4", "number", check=_alpha),
        _Field("gamma4", "rate", check=_non_negative),
        _Field("g", "number"),
        _Field("c2", "time", check=_positive),
        _Field("c4", "time", check=_positive),
    )

    @property
    def alpha(self) -> dict:
        """α_k by k (2 and 4)."""
        return {2: float(self.alpha2), 4: float(self.alpha4)}

    @property
    def gamma_per_ps(self) -> dict:
        """Γ_k by k (2 and 4), per ps: as given, or |g| / C_k."""
        if self.gamma2 is not None:
            return {2: _q(self.gamma2).to("/ps"), 4: _q(self.gamma4).to("/ps")}
        return {k: abs(self.g) / _q(c).to("ps") for k, c in ((2, self.c2), (4, self.c4))}


@dataclass
class Lifetime:
    """A lifetime measurement (``[lifetime]``), with a Coulomb-excitation or fusion-evaporation reaction.

    ``method`` "rdds" is a plunger: the target-to-stopper ``distances`` (each plus ``offset``, the distance at
    electrical contact, 0 if left out), the foil after the gap (:class:`Stopper`) and, for a differential plunger,
    the second foil (:class:`Degrader2`); ``deorientation`` optionally describes how the recoil's alignment fades
    in vacuum. ``method`` "dsam" uses the target's backing (or a thick target) to slow the recoil and has no other
    fields.
    """

    method: str
    distances: Optional[list] = None
    offset: Optional[QuantityLike] = None
    deorientation: Optional[Deorientation] = None
    stopper: Optional[Stopper] = None
    degrader2: Optional[Degrader2] = None

    SPECS = (
        _Field("method", "str", required=True, check=_one_of(*LIFETIME_METHODS)),
        _Field("distances", "lengths"),
        _Field("offset", "length", check=_non_negative),
    )

    @property
    def offset_um(self) -> float:
        """The distance at electrical contact, µm (0 if not given)."""
        return _q(self.offset).to("um") if self.offset is not None else 0.0

    @property
    def distances_um(self) -> list:
        """The distances as given, µm (without the offset)."""
        return [_q(x).to("um") for x in self.distances or []]

    @property
    def flight_distances_um(self) -> list:
        """Target-to-stopper distances the recoil flies, µm: each distance plus the offset."""
        return [d + self.offset_um for d in self.distances_um]

    def to_dict(self) -> dict:
        d = {s.name: _out(getattr(self, s.name)) for s in self.SPECS if getattr(self, s.name) is not None}
        if self.deorientation is not None:
            d["deorientation"] = _to_dict(self.deorientation)
        if self.stopper is not None:
            d["stopper"] = _to_dict(self.stopper)
        if self.degrader2 is not None:
            d["degrader2"] = _to_dict(self.degrader2)
        return d


@dataclass
class Experiment:
    """A complete setup. Build it in Python, or :meth:`load` it from a setup file."""

    title: str
    beam: Beam
    target: Target
    detectors: list[Detector]
    run: Run
    description: Optional[str] = None
    reaction: str = "elastic"
    #: The excited state, for ``reaction = "coulex"``.
    excitation: Optional[Excitation] = None
    gamma_detectors: list = field(default_factory=list)
    #: Level schemes (:class:`~physim.nuclear.levels.LevelScheme`) of the nuclei, by role: "beam" and "target",
    #: and for fusion-evaporation each residue by its name ("158Er").
    levels: dict = field(default_factory=dict)
    chamber: Optional[Chamber] = None
    #: The fusion-evaporation reaction and its exit channels, for ``reaction = "fusion-evaporation"``.
    fusion: Optional[FusionEvaporation] = None
    #: The neutron wall, for a fusion-evaporation reaction.
    neutron_wall: Optional[NeutronWallSetup] = None
    #: The other ancillary detectors of a fusion-evaporation reaction (:data:`ANCILLARIES`).
    charged_particle_ball: Optional[ChargedParticleBallSetup] = None
    recoil_separator: Optional[RecoilSeparatorSetup] = None
    gamma_calorimeter: Optional[GammaCalorimeterSetup] = None
    #: A lifetime measurement (plunger or DSAM), for a Coulomb-excitation or fusion-evaporation reaction.
    lifetime: Optional[Lifetime] = None

    # -- reading ------------------------------------------------------------------------------------------------

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Experiment:
        """Build from the structure of a setup file. Raises :class:`SetupError` listing every problem."""
        problems: list[str] = []
        if not isinstance(data, dict):
            raise SetupError(["the setup must be a table of sections"])
        known = ("schema", "title", "description", "reaction", "beam", "target", "detectors", "run",
                 "gamma_detectors", "levels", "chamber", "neutron_wall", "lifetime") + tuple(a[0] for a in ANCILLARIES)
        for key in data:
            if key not in known:
                hint = difflib.get_close_matches(key, known, n=1)
                problems.append(f"unknown section '{key}'" + (f"; did you mean '{hint[0]}'?" if hint else ""))
        schema = data.get("schema")
        if schema is None:
            problems.append(f"'schema' is missing; the first line should be schema = \"{SCHEMA}\"")
        elif schema != SCHEMA:
            problems.append(f"schema '{schema}' is not supported; this version of physim reads '{SCHEMA}'")
        title = data.get("title")
        if not isinstance(title, str) or not title.strip():
            problems.append("'title' is missing or empty")
        description = data.get("description")
        if description is not None and not isinstance(description, str):
            problems.append("description must be text")

        reaction = "elastic"
        excitation = None
        fusion = None
        raw_channels = None
        if "reaction" in data:
            r = _read_section(data["reaction"], (_Field("type", "str", required=True),) + Excitation.SPECS
                              + FusionEvaporation.SPECS, "reaction", problems, extra=("channels",))
            reaction = r.get("type", reaction)
            if "type" in r and reaction not in REACTIONS:
                problems.append(f"reaction: type '{reaction}' is not available yet; available: {', '.join(REACTIONS)}")
            raw_r = data["reaction"] if isinstance(data["reaction"], dict) else {}
            fusion_names = {s.name for s in FusionEvaporation.SPECS}
            fusion_values = {k: v for k, v in r.items() if k in fusion_names}
            state = {k: v for k, v in r.items() if k != "type" and k not in fusion_names}
            if reaction != "fusion-evaporation":
                given = sorted(set(fusion_values) | ({"channels"} if "channels" in raw_r else set()))
                if given:
                    problems.append(f"reaction: {', '.join(given)} only apply to type = \"fusion-evaporation\"")
            else:
                fusion = FusionEvaporation(**fusion_values)
                raw_channels = raw_r.get("channels")
            if reaction == "fusion-evaporation":
                if state:
                    problems.append(f"reaction: {', '.join(sorted(state))} only apply to type = \"coulex\"")
            elif reaction == "coulex":
                hints = {"energy": "the excited state's energy, e.g. \"1.454 MeV\"",
                         "b_up": "B(Eλ↑) of the excited state, e.g. \"0.0695 e2b2\""}
                for need, hint in hints.items():
                    if need not in raw_r:
                        problems.append(f"reaction: '{need}' is missing: {hint}")
                if "energy" in state and "b_up" in state:
                    excitation = Excitation(**state)
                    try:
                        excitation.b_up_e2fm
                    except ValueError as e:
                        problems.append(f"reaction: b_up {_lower_first(str(e))}")
                        excitation = None
            elif state:
                problems.append(f"reaction: {', '.join(sorted(state))} only apply to type = \"coulex\"")

        beam = None
        if "beam" not in data:
            problems.append("the [beam] section is missing")
        else:
            b = _read_section(data["beam"], Beam.SPECS, "beam", problems)
            if "nuclide" in b and "charge_state" in b:
                z = parse_nuclide(b["nuclide"])[0]
                if not 1 <= b["charge_state"] <= z:
                    problems.append(f"beam: charge_state must be between 1 and {z} for {b['nuclide']} "
                                    f"(Z = {z}), got {b['charge_state']}")
            if "nuclide" in b:
                z, a = parse_nuclide(b["nuclide"])
                if z == 0:
                    problems.append("beam: nuclide must be charged; a neutron beam cannot be planned here")
                    del b["nuclide"]
                else:
                    try:
                        _data.nuclide((z, a))
                    except ValueError as e:
                        problems.append(f"beam: nuclide {e}")
                        del b["nuclide"]
            if all(k in b for k in ("nuclide", "energy", "current")):
                beam = Beam(**b)

        target = None
        if "target" not in data:
            problems.append("the [target] section is missing")
        else:
            t = _read_section(data["target"], Target.SPECS, "target", problems, extra=("backing",))
            backing = None
            raw_t = data["target"] if isinstance(data["target"], dict) else {}
            if "backing" in raw_t:
                lb = _read_section(raw_t["backing"], Layer.SPECS, "target.backing", problems)
                _check_material("target.backing", lb, problems)
                if "material" in lb and "thickness" in lb:
                    backing = Layer(**lb)
            if "ladder" in t:
                k = t.get("selected", 1)
                if not 1 <= k <= len(t["ladder"]):
                    problems.append(f"target: selected = {k}, but the ladder has {len(t['ladder'])} targets")
                else:
                    t["material"], t["thickness"] = t["ladder"][k - 1]
            for name in ("material", "thickness"):
                if name not in t and not any(p.startswith(f"target: {name}") for p in problems):
                    problems.append(f"target: '{name}' is missing" + (" (or a ladder)" if name == "material" else ""))
            _check_material("target", t, problems)
            if "material" in t and "thickness" in t:
                target = Target(**t, backing=backing)

        detectors: list[Detector] = []
        raw_d = data.get("detectors")
        if reaction == "fusion-evaporation" and (raw_d is None or raw_d == []):
            pass  # a fusion-evaporation setup needs no particle detector: the γ rays and ancillaries count
        elif raw_d is None:
            problems.append("no detectors: add at least one [[detectors]] section")
        elif not isinstance(raw_d, list) or not raw_d or not all(isinstance(d, dict) for d in raw_d):
            problems.append("detectors must be written as [[detectors]] sections")
        else:
            names_seen: dict[str, int] = {}
            for i, raw in enumerate(raw_d, start=1):
                label = f"detector {i}" + (f" ({raw['name']})" if isinstance(raw.get("name"), str) else "")
                d = _read_detector(raw, label, problems)
                name = raw.get("name") if isinstance(raw.get("name"), str) else f"D{i}"
                if name in names_seen:
                    problems.append(f"{label}: name '{name}' is already used by detector {names_seen[name]}")
                names_seen.setdefault(name, i)
                if d is not None:
                    detectors.append(d)

        gammas: list = []
        raw_g = data.get("gamma_detectors")
        if raw_g is not None:
            if not isinstance(raw_g, list) or not all(isinstance(g, dict) for g in raw_g):
                problems.append("gamma detectors must be written as [[gamma_detectors]] sections")
            else:
                for i, raw in enumerate(raw_g, start=1):
                    g = _read_gamma_detector(raw, f"gamma detector {i}", problems)
                    if g is not None:
                        gammas.append(g)

        chamber = None
        if "chamber" in data:
            c = _read_section(data["chamber"], Chamber.SPECS, "chamber", problems)
            if "radius" in c:
                chamber = Chamber(**c)
                _check_chamber(chamber, detectors, gammas, problems)

        neutron_wall = None
        if "neutron_wall" in data:
            if reaction != "fusion-evaporation":
                problems.append("neutron_wall: a neutron wall counts evaporated neutrons; it needs [reaction] type = "
                                "\"fusion-evaporation\"")
            w = _read_section(data["neutron_wall"], NeutronWallSetup.SPECS, "neutron_wall", problems)
            if "cells" in w and "distance" in w:
                neutron_wall = NeutronWallSetup(**w)
                try:
                    neutron_wall.model().theta_max()
                except ValueError as e:
                    problems.append(f"neutron_wall: {e}")
        ancillaries = {}
        for name, kind, what in ANCILLARIES:
            if name not in data:
                continue
            if reaction != "fusion-evaporation":
                problems.append(f"{name}: {what}; it needs [reaction] type = \"fusion-evaporation\"")
            a = _read_section(data[name], kind.SPECS, name, problems)
            if all(spec.name in a for spec in kind.SPECS if spec.required):
                ancillaries[name] = kind(**a)
        ball = ancillaries.get("charged_particle_ball")
        if ball is not None and ball.theta_min is not None and ball.theta_max is not None and \
                _q(ball.theta_min).to("deg") >= _q(ball.theta_max).to("deg"):
            problems.append("charged_particle_ball: theta_min must be below theta_max")

        lifetime = None
        if "lifetime" in data:
            lifetime = _read_lifetime(data["lifetime"], reaction, problems)

        if target is not None and target.position is not None:
            origin = (0.0, 0.0, target.position_mm)
            for d in detectors + gammas:
                d._origin = origin
        residues: tuple = ()
        if fusion is not None:
            fusion.channels = _read_channels(raw_channels, beam, target, problems)
            residues = tuple(dict.fromkeys(ch._residue for ch in fusion.channels if ch._residue))
        levels = _read_levels(data.get("levels"), beam, target, problems, residues)

        run = None
        if "run" not in data:
            problems.append("the [run] section is missing")
        else:
            r = _read_section(data["run"], Run.SPECS, "run", problems)
            if "channel" in r:
                _check_studied_channel(r["channel"], reaction, fusion, problems)
            if "beam_time" in r:
                run = Run(**r)

        if lifetime is not None and beam is not None and target is not None:
            channel = run.channel if run is not None else None
            _check_lifetime_foils(lifetime, _recoil(beam, target, excitation, fusion, channel), target, problems)

        if problems:
            raise SetupError(problems)
        return cls(title=title, beam=beam, target=target, detectors=detectors, run=run,
                   description=description, reaction=reaction, excitation=excitation, gamma_detectors=gammas,
                   levels=levels, chamber=chamber, fusion=fusion, neutron_wall=neutron_wall, lifetime=lifetime,
                   **ancillaries)

    @classmethod
    def from_toml(cls, text: str) -> Experiment:
        """Read a setup from TOML text."""
        try:
            data = _toml.loads(text)
        except _toml.TOMLError as e:
            raise SetupError([f"the file is not valid TOML: {e}"]) from None
        return cls.from_dict(data)

    @classmethod
    def load(cls, path: Union[str, Path]) -> Experiment:
        """Read a setup file."""
        return cls.from_toml(Path(path).read_text(encoding="utf-8"))

    @classmethod
    def example(cls, name: str) -> Experiment:
        """One of the example setups shipped with physim (see :func:`example_names`)."""
        path = EXAMPLES / f"{name}.toml"
        if not path.is_file():
            raise ValueError(f"no example setup '{name}'; available: {', '.join(example_names())}")
        return cls.load(path)

    # -- writing ------------------------------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """The structure of a setup file (quantities as text with their units)."""
        d: dict[str, Any] = {"schema": SCHEMA, "title": self.title}
        if self.description is not None:
            d["description"] = self.description
        d["reaction"] = {"type": self.reaction}
        if self.excitation is not None:
            d["reaction"].update(_to_dict(self.excitation))
        if self.fusion is not None:
            d["reaction"].update({s.name: _out(getattr(self.fusion, s.name)) for s in FusionEvaporation.SPECS
                                  if getattr(self.fusion, s.name) is not None})
            d["reaction"]["channels"] = [_to_dict(ch) for ch in self.fusion.channels]
        d["beam"] = _to_dict(self.beam)
        target = _to_dict(self.target)
        if self.target.backing is not None:
            target["backing"] = _to_dict(self.target.backing)
        d["target"] = target
        d["run"] = _to_dict(self.run)
        if self.detectors or self.fusion is None:
            d["detectors"] = [_to_dict(det) for det in self.detectors]
        if self.gamma_detectors:
            d["gamma_detectors"] = [_to_dict(g) for g in self.gamma_detectors]
        if self.chamber is not None:
            d["chamber"] = _to_dict(self.chamber)
        if self.neutron_wall is not None:
            d["neutron_wall"] = _to_dict(self.neutron_wall)
        for name, _, _ in ANCILLARIES:
            if getattr(self, name) is not None:
                d[name] = _to_dict(getattr(self, name))
        if self.lifetime is not None:
            d["lifetime"] = self.lifetime.to_dict()
        if self.levels:
            d["levels"] = {role: scheme.to_dict() for role, scheme in self.levels.items()}
        return d

    def to_toml(self) -> str:
        """The setup as TOML text."""
        return _toml.dumps(self.to_dict())

    def save(self, path: Union[str, Path]) -> None:
        """Validate, then write the setup file."""
        self.validate()
        Path(path).write_text(self.to_toml(), encoding="utf-8")

    def lifetime_recoil(self) -> Optional[tuple]:
        """The nucleus whose lifetime a lifetime measurement measures, as the setup's own numbers give it:
        (nuclide, energy in MeV where it is made, energy in MeV as it leaves the target and its backing; 0 when it
        stops there). It is made in the middle of the target, flying along the beam: for Coulomb excitation in a
        head-on collision (θ_cm = 180°, where excitation is strongest), for fusion-evaporation the studied residue
        at the compound nucleus's speed. None when the emitter does not fly forward (a projectile backscattered in
        normal kinematics) or the numbers cannot be worked out. A single typical recoil: the spread of energies
        and directions is the recoil transport's (:mod:`physim.nuclear.transport`)."""
        channel = self.run.channel if self.run is not None else None
        return _recoil(self.beam, self.target, self.excitation, self.fusion, channel)

    def validate(self) -> None:
        """Check the whole setup, including any edits made in Python. Raises :class:`SetupError`."""
        try:
            d = self.to_dict()
        except (AttributeError, TypeError) as e:
            raise SetupError([f"the setup cannot be written out: {e}"]) from None
        Experiment.from_dict(d)


def _read_gamma_detector(raw: dict[str, Any], label: str, problems: list[str]) -> Optional[GammaDetector]:
    before = len(problems)
    try:
        raw = catalogue.with_defaults(raw, "gamma")
    except ValueError as e:
        problems.append(f"{label}: model {e}")
        return None
    g = _read_section(raw, GammaDetector.SPECS, label, problems)
    if "radius" not in raw and "crystal_diameter" not in raw:
        problems.append(f"{label}: give 'radius', or a 'model' (one of {', '.join(catalogue.names('gamma'))}), or "
                        "'crystal_diameter'")
    if g.get("crystals") == 4 and "crystal_pitch" not in raw:
        problems.append(f"{label}: four crystals need 'crystal_pitch', the distance between their centres")
    if "crystals" in raw and "crystal_diameter" not in raw:
        problems.append(f"{label}: 'crystals' needs 'crystal_diameter'")
    if len(problems) > before or not all(k in g for k in ("theta", "distance")):
        return None
    from . import response

    if g.get("material", "Ge") not in response.CRYSTALS:
        problems.append(f"{label}: material '{g['material']}' has no γ-ray response; use one of "
                        f"{', '.join(response.CRYSTALS)}")
        return None
    for material, thickness in g.get("absorbers", ()):
        try:
            response.transmission(material, thickness, 1.0)
        except ValueError as e:
            problems.append(f"{label}: absorber {material}: {e}")
            return None
    if g.get("addback") and g.get("crystals") != 4:
        problems.append(f"{label}: add-back needs a clover (crystals = 4)")
    if g.get("shield") is not None:
        if g["shield"] not in response.SHIELDS:
            problems.append(f"{label}: shield '{g['shield']}' is not known; use one of {', '.join(response.SHIELDS)}")
        if "housing_side" not in g:
            problems.append(f"{label}: a shield needs 'housing_side', the housing it surrounds")
    if len(problems) > before:
        return None
    return GammaDetector(**g)


def _check_chamber(chamber: Chamber, detectors: list, gammas: list, problems: list[str]) -> None:
    """Particle detectors must fit inside the chamber, and γ-ray detectors must stand outside its wall."""
    radius = _q(chamber.radius).to("mm")
    outer = radius + (_q(chamber.wall_thickness).to("mm") if chamber.wall_thickness is not None else 0.0)
    for i, d in enumerate(detectors, start=1):
        far = math.hypot(*d.chamber_position_mm()) + max(_q(x).to("mm") for x in (
            d.outer_radius, d.radius, d.width, d.height) if x is not None)
        if far > radius:
            problems.append(f"detector {i} ({d.name or f'D{i}'}): reaches {far:.0f} mm from the target, outside "
                            f"the chamber (radius {radius:g} mm)")
    for i, g in enumerate(gammas, start=1):
        near = math.hypot(*g.chamber_centre_mm()) - (_q(g.window_gap).to("mm") if g.window_gap is not None
                                                      else 0.0)
        if near < outer:
            problems.append(f"gamma detector {i} ({g.name or f'G{i}'}): its front is {near:.0f} mm from the "
                            f"target, inside the chamber wall (outer radius {outer:g} mm)")


def _read_detector(raw: dict[str, Any], label: str, problems: list[str]) -> Optional[Detector]:
    before = len(problems)
    try:
        raw = catalogue.with_defaults(raw, "particle")
    except ValueError as e:
        problems.append(f"{label}: model {e}")
        return None
    d = _read_section(raw, Detector.SPECS, label, problems)
    shape = d.get("shape")
    if shape is not None and shape not in SHAPES:
        problems.append(f"{label}: shape '{shape}' is not known; use one of {', '.join(SHAPES)}")
        return None
    if shape is not None:
        size_fields = SHAPES[shape]
        for other, other_fields in SHAPES.items():
            for f in other_fields:
                if f in raw and f not in size_fields:
                    problems.append(f"{label}: '{f}' is not used by a {shape} detector (it belongs to {other})")
        required = {"rectangle": ("width", "height"), "annular": ("inner_radius", "outer_radius"),
                    "circle": ("radius",)}[shape]
        for f in required:
            if f not in raw:
                problems.append(f"{label}: a {shape} detector needs '{f}'")
        if shape == "annular" and "inner_radius" in d and "outer_radius" in d:
            if d["inner_radius"].to("mm") >= d["outer_radius"].to("mm"):
                problems.append(f"{label}: inner_radius must be smaller than outer_radius")
    spherical = [k for k in ("theta", "phi", "distance") if k in raw]
    if "position" in raw and spherical:
        problems.append(f"{label}: give either position or theta/phi/distance, not both")
    elif "position" not in raw:
        for k in ("theta", "distance"):
            if k not in raw:
                problems.append(f"{label}: '{k}' is missing (or give the detector a position)")
    if len(problems) > before:
        return None
    _check_material(label, {"material": d.get("material", "Si"), **{k: d[k] for k in ("density", "thickness",
                                                                                     "dead_layer") if k in d}},
                    problems)
    if len(problems) > before:
        return None
    return Detector(**d)


def main_nuclide(target: Target) -> str:
    """The target's most abundant nuclide (the one a fusion-evaporation reaction fuses with), as AME names it."""
    z, a, _ = max(_data.material(target.material).atoms, key=lambda t: t[2])
    return _data.nuclide((z, a)).name


def _check_studied_channel(which: str, reaction: str, fusion, problems: list[str]) -> None:
    """``[run] channel``: a fusion-evaporation setup's studied channel must be one of its table, named by its
    particles ("4n", in any order) or its residue ("158Dy")."""
    if reaction != "fusion-evaporation":
        problems.append("run: channel is the exit channel a fusion-evaporation experiment studies; it needs "
                        "[reaction] type = \"fusion-evaporation\"")
        return
    if fusion is None or not fusion.channels:
        return  # the channel table's own problems say what is missing
    if studied_channel_index(fusion.channels, which) is None:
        names = ", ".join(ch.particles for ch in fusion.channels)
        problems.append(f"run: channel '{which}' is not in the channel table ({names})")


def studied_channel_index(channels: list, which: Optional[str]) -> Optional[int]:
    """The index of the channel ``which`` names (its particles in any order, or its residue) in a channel table,
    or None. With ``which`` None, the channel with the largest cross section."""
    from .residues import parse_particles

    if not channels:
        return None
    if which is None:
        return max(range(len(channels)), key=lambda i: _q(channels[i].cross_section).canonical)
    try:
        key = tuple(sorted(parse_particles(which)))
    except ValueError:
        key = None
    for i, ch in enumerate(channels):
        if key is not None and tuple(sorted(parse_particles(ch.particles))) == key:
            return i
        if ch._residue is not None and _residue_name(which) == ch._residue:
            return i
    return None


def _residue_name(text: str) -> Optional[str]:
    try:
        return _data.nuclide(parse_nuclide(text)).name
    except ValueError:
        return None


def _read_channels(raw: Any, beam: Optional[Beam], target: Optional[Target], problems: list[str]) -> list:
    """The ``[[reaction.channels]]`` of a fusion-evaporation reaction. Each channel's particles must leave a
    residue of the compound nucleus (beam + the target's main nuclide) that the mass table knows; a residue
    given must be that one; no two channels may emit the same particles."""
    from .residues import parse_particles, residue_of

    if raw is None:
        problems.append("reaction: a fusion-evaporation reaction needs its exit channels: add [[reaction.channels]] "
                        "sections, each with particles (e.g. \"4n\") and a cross_section (e.g. \"300 mb\")")
        return []
    if not isinstance(raw, list) or not all(isinstance(c, dict) for c in raw):
        problems.append("reaction: channels must be written as [[reaction.channels]] sections")
        return []
    if not raw:
        problems.append("reaction: the channel table is empty: add at least one [[reaction.channels]] section")
        return []
    compound = None
    if beam is not None and target is not None:
        try:
            b, t = parse_nuclide(beam.nuclide), parse_nuclide(main_nuclide(target))
            compound = _data.nuclide((b[0] + t[0], b[1] + t[1]))
        except ValueError as e:
            problems.append(f"reaction: the compound nucleus of {beam.nuclide} + {target.material}: {e}")
    out, seen = [], {}
    for i, rc in enumerate(raw, start=1):
        label = f"reaction channel {i}" + (f" ({rc['particles']})" if isinstance(rc.get("particles"), str) else "")
        before = len(problems)
        c = _read_section(rc, ExitChannel.SPECS, label, problems)
        parts = None
        if "particles" in c:
            try:
                parts = parse_particles(c["particles"])
            except ValueError as e:
                problems.append(f"{label}: particles {_lower_first(str(e))}")
        if parts is not None:
            key = tuple(sorted(parts))
            if key in seen:
                problems.append(f"{label}: the same particles as channel {seen[key]}")
            seen.setdefault(key, i)
        residue = None
        if parts is not None and compound is not None:
            try:
                residue = _data.nuclide(residue_of((compound.Z, compound.A), parts)).name
            except ValueError as e:
                problems.append(f"{label}: {compound.name} cannot emit {c['particles']}: {_lower_first(str(e))}")
            if residue is not None and "residue" in c and _data.nuclide(c["residue"]).name != residue:
                problems.append(f"{label}: residue {c['residue']} is not what {compound.name} leaves after emitting "
                                f"{c['particles']} ({residue})")
        if len(problems) > before or not all(k in c for k in ("particles", "cross_section")):
            continue
        ch = ExitChannel(**c)
        ch._residue = residue
        out.append(ch)
    return out


def _read_levels(raw: Any, beam: Optional[Beam], target: Optional[Target], problems: list[str],
                 residues: tuple = ()) -> dict:
    """The ``[levels.beam]`` and ``[levels.target]`` sections, and for fusion-evaporation ``[levels.<residue>]``
    (``[levels.158Er]``): each scheme must belong to a nucleus of its role."""
    from .levels import LevelScheme

    if raw is None:
        return {}
    if not isinstance(raw, dict):
        problems.append("levels must be written as [levels.beam] and [levels.target] sections")
        return {}
    out = {}
    for role, section in raw.items():
        if role not in ("beam", "target") + tuple(residues):
            extra = "".join(f", [levels.{r}]" for r in residues)
            problems.append(f"levels: unknown section '{role}'; use [levels.beam] or [levels.target]{extra}")
            continue
        scheme = LevelScheme.from_dict(section, f"levels.{role}", problems)
        if scheme is None:
            continue
        za = parse_nuclide(scheme.nuclide)
        if role in residues and za != parse_nuclide(role):
            problems.append(f"levels.{role}: nuclide {scheme.nuclide} is not the residue {role}")
        elif role == "beam" and beam is not None and za != parse_nuclide(beam.nuclide):
            problems.append(f"levels.beam: nuclide {scheme.nuclide} is not the beam ({beam.nuclide})")
        elif role == "target" and target is not None:
            try:
                present = {(z, a) for z, a, _ in _data.material(target.material).atoms}
            except ValueError:
                present = {za}  # the material's own problem is reported with the target
            if za not in present:
                problems.append(f"levels.target: nuclide {scheme.nuclide} is not in the target "
                                f"({target.material})")
        out[role] = scheme
    return out


def _read_lifetime(raw: Any, reaction: str, problems: list[str]) -> Optional[Lifetime]:
    """The ``[lifetime]`` section: its method, and for a plunger its distances and foils. Problems for a setup
    that cannot be measured that way (no distances, no stopper, a second foil without a differential plunger)."""
    if reaction not in ("coulex", "fusion-evaporation"):
        problems.append("lifetime: a lifetime measurement (plunger or DSAM) needs [reaction] type = \"coulex\" or "
                        "\"fusion-evaporation\"")
    before = len(problems)
    values = _read_section(raw, Lifetime.SPECS, "lifetime", problems,
                           extra=("stopper", "degrader2", "deorientation"))
    raw_l = raw if isinstance(raw, dict) else {}
    method = values.get("method")
    parts: dict[str, Any] = {}
    for name, kind in (("stopper", Stopper), ("degrader2", Degrader2)):
        if name in raw_l:
            v = _read_section(raw_l[name], kind.SPECS, f"lifetime.{name}", problems)
            _check_material(f"lifetime.{name}", v, problems)
            if all(spec.name in v for spec in kind.SPECS if spec.required):
                parts[name] = kind(**v)
    if "deorientation" in raw_l:
        parts["deorientation"] = _read_deorientation(raw_l["deorientation"], problems)
    if method == "rdds":
        if not values.get("distances") and not any(p.startswith("lifetime: distances") for p in problems):
            problems.append("lifetime: a plunger (RDDS) needs its target-to-stopper distances, e.g. "
                            "distances = [\"10 um\", \"50 um\", \"200 um\"]")
        raw_s = raw_l.get("stopper")
        kind = raw_s.get("kind") if isinstance(raw_s, dict) else None
        if raw_s is None:
            problems.append("lifetime: a plunger (RDDS) needs the foil after the gap: add [lifetime.stopper] with "
                            "kind (\"stopper\", \"degrader\" or \"differential\"), material and thickness")
        if kind == "differential" and "degrader2" not in raw_l:
            problems.append("lifetime: a differential plunger needs its second foil, which stops the recoil: add "
                            "[lifetime.degrader2] with material, thickness and gap")
        elif kind in ("stopper", "degrader") and "degrader2" in raw_l:
            problems.append("lifetime.degrader2: a second foil only applies to a stopper of kind \"differential\"")
    elif method == "dsam":
        given = [n for n in ("distances", "offset", "deorientation", "stopper", "degrader2") if n in raw_l]
        if given:
            problems.append(f"lifetime: {', '.join(given)} only apply to method = \"rdds\" (a plunger); a DSAM "
                            "measurement slows the recoil in the target's backing")
    if len(problems) > before or method is None:
        return None
    return Lifetime(method=method, distances=values.get("distances"), offset=values.get("offset"), **parts)


def _read_deorientation(raw: Any, problems: list[str]) -> Optional[Deorientation]:
    """``lifetime.deorientation``: α₂, α₄ and either the rates Γ₂, Γ₄ or g with C₂, C₄ — one form, complete."""
    before = len(problems)
    v = _read_section(raw, Deorientation.SPECS, "lifetime.deorientation", problems)
    raw_d = raw if isinstance(raw, dict) else {}
    rates = [n for n in ("gamma2", "gamma4") if n in raw_d]
    via_g = [n for n in ("g", "c2", "c4") if n in raw_d]
    forms = "give alpha2 and alpha4 with either gamma2 and gamma4 (\"0.142 /ps\"), or g, c2 and c4 (\"3.3 ps\")"
    if not isinstance(raw, dict):
        return None
    if rates and via_g:
        problems.append(f"lifetime.deorientation: {', '.join(rates + via_g)} mix the two forms; {forms}")
    else:
        need = ["alpha2", "alpha4"] + (["g", "c2", "c4"] if via_g else ["gamma2", "gamma4"])
        missing = [n for n in need if n not in raw_d]
        if missing:
            problems.append(f"lifetime.deorientation: {', '.join(missing)} missing; {forms}")
    if len(problems) > before:
        return None
    return Deorientation(**v)


@lru_cache(maxsize=64)
def _stopping(ion: str, material: str, density: Optional[str]):
    from .stopping import Stopping

    return Stopping(ion, _data.material(material, Quantity.parse(density) if density else None))


def _slowing(ion: str, layer, energy: float, tilt_deg: float = 0.0, share: float = 1.0) -> tuple:
    """The stopping of ``ion`` in ``layer`` (a Target or Layer), the mg/cm² it crosses (``share`` of the layer)
    and the energy left after them."""
    density = str(_q(layer.density)) if layer.density is not None else None
    s = _stopping(ion, layer.material, density)
    areal = share * layer.material_data().areal_density_mg_cm2(_q(layer.thickness))
    return s, areal, float(s.energy_after(energy, areal, tilt_deg))


def _recoil(beam, target, excitation, fusion, channel) -> Optional[tuple]:
    """See :meth:`Experiment.lifetime_recoil`."""
    try:
        tilt = _q(target.tilt).to("deg") if target.tilt is not None else 0.0
        _, _, e_mid = _slowing(beam.nuclide, target, beam.energy_mev, tilt, share=0.5)
        if e_mid <= 0:
            return None
        if excitation is not None:
            from .kinematics import TwoBody

            k = TwoBody(beam.nuclide, main_nuclide(target), e_mid, excitation_mev=excitation.energy_mev,
                        excite="recoil" if excitation.excite == "target" else "ejectile")
            point = k.recoil_for(180.0) if excitation.excite == "target" else k.at_cm(180.0, "ejectile")
            if float(point.theta_lab) >= 90.0:
                return None
            emitter = (k.recoil if excitation.excite == "target" else k.ejectile).name
            energy = float(point.energy)
        elif fusion is not None and fusion.channels:
            i = studied_channel_index(fusion.channels, channel)
            if i is None or fusion.channels[i]._residue is None:
                return None
            emitter = fusion.channels[i]._residue
            a_beam, a_target = beam.A, parse_nuclide(main_nuclide(target))[1]
            # The compound nucleus takes the beam's momentum; the residue keeps its speed.
            energy = e_mid * a_beam / (a_beam + a_target) * parse_nuclide(emitter)[1] / (a_beam + a_target)
        else:
            return None
        _, _, left = _slowing(emitter, target, energy, tilt, share=0.5)
        if target.backing is not None and left > 0:
            _, _, left = _slowing(emitter, target.backing, left, tilt)
        return emitter, energy, left
    except (ValueError, KeyError, ZeroDivisionError):
        return None


def _check_lifetime_foils(lifetime: Lifetime, recoil: Optional[tuple], target: Target, problems: list[str]) -> None:
    """Problems a lifetime setup has with its recoil, from the setup alone: for a plunger, a recoil that never
    leaves the target, a stopper that does not stop it or a degrader that does; for DSAM, nothing to stop it in."""
    if lifetime.method == "dsam":
        if target.backing is None and (recoil is None or recoil[2] > 0):
            problems.append("lifetime: a DSAM measurement needs the recoil to slow down and stop in matter: give the "
                            "target a backing ([target.backing]), or a target thick enough to stop the recoil")
        elif target.backing is not None and recoil is not None and recoil[2] > 0:
            problems.append(f"lifetime: the backing ({target.backing.thickness} of {target.backing.material}) does not "
                            f"stop the recoil ({recoil[0]} made at {recoil[1]:.3g} MeV leaves it at {recoil[2]:.3g} "
                            "MeV); a DSAM backing must stop it: make the backing thicker")
        return
    if recoil is None or lifetime.stopper is None:
        return
    emitter, made, left = recoil
    if left <= 0:
        where = "target and its backing" if target.backing is not None else "target"
        problems.append(f"lifetime: the recoil ({emitter}, made at {made:.3g} MeV) stops in the {where} and never "
                        "reaches the plunger's gap")
        return
    st = lifetime.stopper
    try:
        s, areal, after = _slowing(emitter, st, left)
    except ValueError:
        return
    need = float(s.range(left))
    foil = f"{st.thickness} of {st.material}"
    if st.kind == "stopper" and after > 0:
        problems.append(f"lifetime.stopper: {foil} does not stop the recoil ({emitter} at {left:.3g} MeV after the "
                        f"target needs {need:.3g} mg/cm2); make it thicker")
    elif st.kind in ("degrader", "differential") and after <= 0:
        problems.append(f"lifetime.stopper: the {st.kind} foil ({foil}) stops the recoil ({emitter} at {left:.3g} "
                        f"MeV after the target, range {need:.3g} mg/cm2); it must let the recoil through: make it "
                        "thinner" + (", or set kind = \"stopper\"" if st.kind == "degrader" else ""))
    elif st.kind == "differential" and lifetime.degrader2 is not None:
        try:
            s2, _, after2 = _slowing(emitter, lifetime.degrader2, after)
        except ValueError:
            return
        if after2 > 0:
            d2 = lifetime.degrader2
            problems.append(f"lifetime.degrader2: {d2.thickness} of {d2.material} does not stop the recoil "
                            f"({emitter} at {after:.3g} MeV after the first foil needs {float(s2.range(after)):.3g} "
                            "mg/cm2); make it thicker")


def _check_material(where: str, values: dict[str, Any], problems: list[str]) -> None:
    """The material must exist in the data, and lengths need a density to convert to areal densities."""
    if "material" not in values:
        return
    try:
        m = _data.material(values["material"], values.get("density"))
    except ValueError as e:
        problems.append(f"{where}: material {_lower_first(str(e))}")
        return
    if m.density_g_cm3 is None:
        for key in ("thickness", "dead_layer"):
            q = values.get(key)
            if isinstance(q, Quantity) and q.kind == "length":
                problems.append(f"{where}: {key} '{q}' is a length, but the density of {values['material']} is "
                                "not known; give 'density' or the thickness in mg/cm2")


__all__ = ["SCHEMA", "REACTIONS", "SHAPES", "EXAMPLES", "LIFETIME_METHODS", "STOPPER_KINDS", "example_names",
           "main_nuclide", "Beam", "Chamber", "Degrader2", "Deorientation", "Detector", "Excitation", "ExitChannel",
           "Experiment", "FusionEvaporation", "GammaDetector", "Layer", "Lifetime", "NeutronWallSetup", "Run",
           "SetupError", "Stopper", "Target"]
