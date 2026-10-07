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
from pathlib import Path
from typing import Any, Callable, Optional, Union

from . import _toml, catalogue
from . import data as _data
from .names import parse_material, parse_nuclide
from .quantity import Quantity

#: Schema identifier written at the top of every setup file.
SCHEMA = "physim.experiment/1"

#: Reactions this version can plan. Others are added slice by slice (backlog item 43 onwards).
REACTIONS = ("elastic", "coulex")

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

    SPECS = (
        _Field("nuclide", "nuclide", required=True),
        _Field("energy", "energy|energy_per_nucleon", required=True, check=_positive),
        _Field("current", "particle_current|electrical_current", required=True, check=_positive),
        _Field("charge_state", "int"),
        _Field("energy_spread", "energy|fraction", check=_non_negative),
        _Field("spot_size", "length", check=_non_negative),
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

    SPECS = (
        _Field("material", "material"),
        _Field("thickness", "areal_density|length", check=_positive),
        _Field("density", "density", check=_positive),
        _Field("tilt", "angle", check=_tilt),
        _Field("position", "length"),
        _Field("ladder", "absorbers"),
        _Field("selected", "int", check=lambda n: None if n >= 1 else "counts from 1"),
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

    SPECS = (
        _Field("beam_time", "time", required=True, check=_positive),
        _Field("counts_wanted", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
        _Field("coincidence_window", "time", check=_positive),
        _Field("dead_time", "time", check=_non_negative),
        _Field("room_background", "rate", check=_non_negative),
        _Field("extra_lines", "lines"),
        _Field("shaping_time", "time", check=_non_negative),
    )


@dataclass
class Excitation:
    """The excited state of a Coulomb-excitation reaction (``[reaction] type = "coulex"``).

    One state, reached from a 0⁺ ground state: its energy, multipolarity (E1, E2, E3), the reduced transition
    probability ``b_up`` = B(Eλ↑) as text with a unit (``"0.0695 e2b2"``), and which nucleus is excited.
    """

    energy: QuantityLike
    b_up: str
    excite: str = "target"
    multipolarity: str = "E2"
    #: How the γ rays leave the nucleus: "correlated" with the scattered particle (also when left out), or
    #: "isotropic".
    emission: Optional[str] = None

    SPECS = (
        _Field("excite", "str", check=lambda s: None if s in ("target", "projectile")
               else "must be 'target' or 'projectile'"),
        _Field("energy", "energy", check=_positive),
        _Field("multipolarity", "str", check=lambda s: None if s in ("E1", "E2", "E3") else "must be E1, E2 or E3"),
        _Field("b_up", "str"),
        _Field("emission", "str", check=lambda s: None if s in ("correlated", "isotropic")
               else "must be 'correlated' or 'isotropic'"),
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
    #: Level schemes (:class:`~physim.nuclear.levels.LevelScheme`) of the nuclei, by role: "beam" and "target".
    levels: dict = field(default_factory=dict)
    chamber: Optional[Chamber] = None

    # -- reading ------------------------------------------------------------------------------------------------

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Experiment:
        """Build from the structure of a setup file. Raises :class:`SetupError` listing every problem."""
        problems: list[str] = []
        if not isinstance(data, dict):
            raise SetupError(["the setup must be a table of sections"])
        known = ("schema", "title", "description", "reaction", "beam", "target", "detectors", "run",
                 "gamma_detectors", "levels", "chamber")
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
        if "reaction" in data:
            r = _read_section(data["reaction"], (_Field("type", "str", required=True),) + Excitation.SPECS,
                              "reaction", problems)
            reaction = r.get("type", reaction)
            if "type" in r and reaction not in REACTIONS:
                problems.append(f"reaction: type '{reaction}' is not available yet; available: {', '.join(REACTIONS)}")
            raw_r = data["reaction"] if isinstance(data["reaction"], dict) else {}
            state = {k: v for k, v in r.items() if k != "type"}
            if reaction == "coulex":
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
        if raw_d is None:
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

        if target is not None and target.position is not None:
            origin = (0.0, 0.0, target.position_mm)
            for d in detectors + gammas:
                d._origin = origin
        levels = _read_levels(data.get("levels"), beam, target, problems)

        run = None
        if "run" not in data:
            problems.append("the [run] section is missing")
        else:
            r = _read_section(data["run"], Run.SPECS, "run", problems)
            if "beam_time" in r:
                run = Run(**r)

        if problems:
            raise SetupError(problems)
        return cls(title=title, beam=beam, target=target, detectors=detectors, run=run,
                   description=description, reaction=reaction, excitation=excitation, gamma_detectors=gammas,
                   levels=levels, chamber=chamber)

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
        d["beam"] = _to_dict(self.beam)
        target = _to_dict(self.target)
        if self.target.backing is not None:
            target["backing"] = _to_dict(self.target.backing)
        d["target"] = target
        d["run"] = _to_dict(self.run)
        d["detectors"] = [_to_dict(det) for det in self.detectors]
        if self.gamma_detectors:
            d["gamma_detectors"] = [_to_dict(g) for g in self.gamma_detectors]
        if self.chamber is not None:
            d["chamber"] = _to_dict(self.chamber)
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


def _read_levels(raw: Any, beam: Optional[Beam], target: Optional[Target], problems: list[str]) -> dict:
    """The ``[levels.beam]`` and ``[levels.target]`` sections: each scheme must belong to a nucleus of its role."""
    from .levels import LevelScheme

    if raw is None:
        return {}
    if not isinstance(raw, dict):
        problems.append("levels must be written as [levels.beam] and [levels.target] sections")
        return {}
    out = {}
    for role, section in raw.items():
        if role not in ("beam", "target"):
            problems.append(f"levels: unknown section '{role}'; use [levels.beam] or [levels.target]")
            continue
        scheme = LevelScheme.from_dict(section, f"levels.{role}", problems)
        if scheme is None:
            continue
        za = parse_nuclide(scheme.nuclide)
        if role == "beam" and beam is not None and za != parse_nuclide(beam.nuclide):
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


__all__ = ["SCHEMA", "REACTIONS", "SHAPES", "EXAMPLES", "example_names", "Beam", "Chamber", "Detector", "Excitation",
           "Experiment", "GammaDetector", "Layer", "Run", "SetupError", "Target"]
