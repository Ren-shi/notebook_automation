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

from . import _toml
from .names import parse_material, parse_nuclide
from .quantity import Quantity

#: Schema identifier written at the top of every setup file.
SCHEMA = "physim.experiment/1"

#: Reactions this version can plan. Others are added slice by slice (backlog item 43 onwards).
REACTIONS = ("elastic",)

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
    q = raw if isinstance(raw, Quantity) else Quantity.parse(raw) if isinstance(raw, str) else None
    if q is None:
        if isinstance(raw, (int, float)) and not isinstance(raw, bool):
            raise ValueError(f"{raw!r} has no unit; write it as text with a unit, e.g. \"{raw} {_example_unit(kind)}\"")
        raise ValueError(f"must be a number with a unit, e.g. \"1 {_example_unit(kind)}\", got {raw!r}")
    kinds = kind.split("|")
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
    return s if s.startswith("'") or s[:2].isupper() else s[:1].lower() + s[1:]


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
    return {n: _out(getattr(obj, n)) for n in names if getattr(obj, n) is not None and n != "backing"}


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

    SPECS = (
        _Field("material", "material", required=True),
        _Field("thickness", "areal_density|length", required=True, check=_positive),
    )


@dataclass
class Target:
    """The target, optionally on a backing. ``tilt`` turns it about the vertical (y) axis."""

    material: str
    thickness: QuantityLike
    tilt: Optional[QuantityLike] = None
    backing: Optional[Layer] = None

    SPECS = (
        _Field("material", "material", required=True),
        _Field("thickness", "areal_density|length", required=True, check=_positive),
        _Field("tilt", "angle", check=_tilt),
    )


@dataclass
class Detector:
    """One detector. Place it with ``theta``/``phi``/``distance`` or with ``position``; size it by ``shape``."""

    shape: str
    thickness: QuantityLike
    name: Optional[str] = None
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
    dead_layer: Optional[QuantityLike] = None
    resolution: Optional[QuantityLike] = None
    threshold: Optional[QuantityLike] = None

    SPECS = (
        _Field("name", "str"),
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
        _Field("thickness", "length|areal_density", required=True, check=_positive),
        _Field("dead_layer", "length|areal_density", check=_non_negative),
        _Field("resolution", "energy", check=_non_negative),
        _Field("threshold", "energy", check=_non_negative),
    )

    @property
    def detector_material(self) -> str:
        """Detector material; silicon unless given."""
        return self.material or "Si"

    def position_mm(self) -> tuple[float, float, float]:
        """Centre of the detector face, mm, in the coordinates described in the module docstring."""
        if self.position is not None:
            x, y, z = (_q(c).to("mm") for c in self.position)
            return (x, y, z)
        th = math.radians(_q(self.theta).to("deg"))
        ph = math.radians(_q(self.phi).to("deg")) if self.phi is not None else 0.0
        d = _q(self.distance).to("mm")
        return (d * math.sin(th) * math.cos(ph), d * math.sin(th) * math.sin(ph), d * math.cos(th))

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

    SPECS = (
        _Field("beam_time", "time", required=True, check=_positive),
        _Field("counts_wanted", "int", check=lambda n: None if n >= 1 else "must be at least 1"),
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

    # -- reading ------------------------------------------------------------------------------------------------

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Experiment:
        """Build from the structure of a setup file. Raises :class:`SetupError` listing every problem."""
        problems: list[str] = []
        if not isinstance(data, dict):
            raise SetupError(["the setup must be a table of sections"])
        known = ("schema", "title", "description", "reaction", "beam", "target", "detectors", "run")
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
        if "reaction" in data:
            r = _read_section(data["reaction"], (_Field("type", "str", required=True),), "reaction", problems)
            reaction = r.get("type", reaction)
            if "type" in r and reaction not in REACTIONS:
                problems.append(f"reaction: type '{reaction}' is not available yet; available: {', '.join(REACTIONS)}")

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
                if len(lb) == 2:
                    backing = Layer(**lb)
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
                   description=description, reaction=reaction)

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
        d["beam"] = _to_dict(self.beam)
        target = _to_dict(self.target)
        if self.target.backing is not None:
            target["backing"] = _to_dict(self.target.backing)
        d["target"] = target
        d["run"] = _to_dict(self.run)
        d["detectors"] = [_to_dict(det) for det in self.detectors]
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


def _read_detector(raw: dict[str, Any], label: str, problems: list[str]) -> Optional[Detector]:
    before = len(problems)
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
    return Detector(**d)


__all__ = ["SCHEMA", "REACTIONS", "SHAPES", "EXAMPLES", "example_names", "Beam", "Detector", "Experiment", "Layer", "Run",
           "SetupError", "Target"]
