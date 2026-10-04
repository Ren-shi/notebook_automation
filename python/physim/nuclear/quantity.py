"""Physical quantities written with their unit, as in a setup file: ``"5.5 MeV"``, ``"40 mm"``.

A :class:`Quantity` keeps the number and unit exactly as given, so a setup file saved again reads the same,
and converts on request with :meth:`Quantity.to`. Each kind of quantity has a canonical unit (the first one
listed in :data:`UNITS`), which is what the physics code works in.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

#: Units accepted for each kind of quantity, as the factor that converts to the canonical (first) unit.
UNITS: dict[str, dict[str, float]] = {
    "energy": {"MeV": 1.0, "eV": 1e-6, "keV": 1e-3, "GeV": 1e3},
    "energy_per_nucleon": {"MeV/u": 1.0, "keV/u": 1e-3, "GeV/u": 1e3, "AMeV": 1.0},
    "length": {"mm": 1.0, "nm": 1e-6, "um": 1e-3, "µm": 1e-3, "cm": 10.0, "m": 1e3},
    "areal_density": {"mg/cm2": 1.0, "ug/cm2": 1e-3, "µg/cm2": 1e-3, "g/cm2": 1e3},
    "angle": {"deg": 1.0, "rad": 180.0 / math.pi, "mrad": 0.18 / math.pi},
    "particle_current": {"pnA": 1.0, "ppA": 1e-3, "puA": 1e3, "pµA": 1e3},
    "electrical_current": {"enA": 1.0, "epA": 1e-3, "euA": 1e3, "eµA": 1e3},
    "time": {"h": 1.0, "s": 1.0 / 3600.0, "min": 1.0 / 60.0, "d": 24.0},
    "fraction": {"%": 1.0},
}

_NUMBER_AND_UNIT = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*([^\s\d.+-]\S*)\s*$")


def _kinds_of(unit: str) -> list[str]:
    return [kind for kind, units in UNITS.items() if unit in units]


@dataclass(frozen=True)
class Quantity:
    """A number with a unit, e.g. ``Quantity(5.5, "MeV")`` or ``Quantity.parse("5.5 MeV")``."""

    value: float
    unit: str

    @classmethod
    def parse(cls, text: str) -> Quantity:
        """Parse ``"<number> <unit>"``. Raises ``ValueError`` with a readable message."""
        if not isinstance(text, str):
            raise ValueError(f"expected a number with a unit such as '40 mm', got {text!r}")
        m = _NUMBER_AND_UNIT.match(text)
        if not m:
            if _NUMBER_AND_UNIT.match(text + " x") and text.strip():
                raise ValueError(f"'{text}' has no unit; write it with one, e.g. '{text.strip()} mm'")
            raise ValueError(f"cannot read '{text}' as a number with a unit, e.g. '40 mm'")
        value, unit = float(m.group(1)), m.group(2)
        if not _kinds_of(unit):
            raise ValueError(f"unknown unit '{unit}' in '{text}'")
        return cls(value, unit)

    @property
    def kind(self) -> str:
        """The kind of quantity (a key of :data:`UNITS`)."""
        kinds = _kinds_of(self.unit)
        if not kinds:
            raise ValueError(f"unknown unit '{self.unit}'")
        return kinds[0]

    def to(self, unit: str) -> float:
        """The value in ``unit``, which must be of the same kind."""
        table = UNITS[self.kind]
        if unit not in table:
            raise ValueError(f"cannot convert {self} to '{unit}': a {self.kind.replace('_', ' ')} "
                             f"is given in {', '.join(table)}")
        return self.value * table[self.unit] / table[unit]

    @property
    def canonical(self) -> float:
        """The value in the canonical unit of its kind (MeV, mm, deg, mg/cm², pnA, h, ...)."""
        return self.value * UNITS[self.kind][self.unit]

    def __str__(self) -> str:
        return f"{format_number(self.value)} {self.unit}"


def format_number(x: float) -> str:
    """Shortest text that reads back as the same float, without a trailing ``.0``."""
    if isinstance(x, float) and x.is_integer() and abs(x) < 1e15:
        return str(int(x))
    return repr(x)
