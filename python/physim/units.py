"""Unit systems. The engine is unit-agnostic: pick a unit system, build the world in its
units (using its ``G``), and convert inputs and outputs here.

Three systems are provided:

- :data:`SI`: metre, kilogram, second.
- :data:`ASTRONOMICAL`: astronomical unit, solar mass, and the Gaussian year
  (``2π / k`` days with Gauss's constant ``k = 0.01720209895``), chosen so that
  ``G = 4π²`` exactly: a massless body on a 1 AU circular orbit about one solar mass has a
  period of exactly 1. The solar mass is the one this implies (``k² AU³ / (day² G_SI)``).
- :data:`DIMENSIONLESS`: ``G = 1`` with no physical scale; it cannot be converted.

Quantities are named (``"length"``, ``"velocity"``, ``"energy"``, ...) or given as a
dimension tuple ``(L, M, T)`` of powers of length, mass and time::

    au = ps.units.ASTRONOMICAL
    v = au.from_si(29.78e3, "velocity")          # Earth's orbital speed, in AU/yr (~6.28)
    au.to_si(1.0, "time") / 86400                # 365.2569 days
    ps.units.convert(1.0, "energy", au, ps.units.SI)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: Newton's constant in SI (CODATA 2018), m³ kg⁻¹ s⁻².
G_SI = 6.67430e-11
#: Astronomical unit, m (IAU 2012, exact).
AU = 1.495978707e11
#: Day, s.
DAY = 86400.0
#: Gauss's gravitational constant, rad/day.
GAUSS_K = 0.01720209895
#: Gaussian year, s: the period of a 1 AU orbit about one solar mass.
GAUSSIAN_YEAR = 2.0 * math.pi / GAUSS_K * DAY
#: Solar mass implied by Gauss's constant, kg.
SOLAR_MASS = GAUSS_K**2 * AU**3 / (DAY**2 * G_SI)

DIMENSIONS = {
    "length": (1, 0, 0),
    "mass": (0, 1, 0),
    "time": (0, 0, 1),
    "velocity": (1, 0, -1),
    "acceleration": (1, 0, -2),
    "force": (1, 1, -2),
    "energy": (2, 1, -2),
    "power": (2, 1, -3),
    "momentum": (1, 1, -1),
    "angular_momentum": (2, 1, -1),
    "moment_of_inertia": (2, 1, 0),
    "torque": (2, 1, -2),
    "frequency": (0, 0, -1),
    "angular_velocity": (0, 0, -1),
    "density": (-3, 1, 0),
    "pressure": (-1, 1, -2),
    "spring_constant": (0, 1, -2),
    "gravitational_parameter": (3, 0, -2),
    "G": (3, -1, -2),
}


def _dimension(quantity):
    if isinstance(quantity, str):
        try:
            return DIMENSIONS[quantity]
        except KeyError:
            raise ValueError(
                f"unknown quantity {quantity!r}; use one of {sorted(DIMENSIONS)} "
                "or a (length, mass, time) tuple of powers"
            ) from None
    dims = tuple(quantity)
    if len(dims) != 3:
        raise ValueError(f"a dimension is a (length, mass, time) tuple, got {quantity!r}")
    return dims


@dataclass(frozen=True)
class UnitSystem:
    """A unit system: the size of its units of length, mass and time in SI (``None`` for a
    dimensionless system), and the value of ``G`` in it."""

    name: str
    length: float | None
    mass: float | None
    time: float | None
    G: float

    def scale(self, quantity):
        """SI value of one unit of ``quantity`` in this system."""
        if self.length is None:
            raise ValueError(f"the {self.name} system has no physical scale to convert")
        a, b, c = _dimension(quantity)
        return self.length**a * self.mass**b * self.time**c

    def to_si(self, value, quantity):
        """Converts ``value`` (scalar or array) from this system to SI."""
        return value * self.scale(quantity)

    def from_si(self, value, quantity):
        """Converts ``value`` from SI to this system."""
        return value / self.scale(quantity)

    def __repr__(self):
        return f"UnitSystem({self.name!r}, G={self.G!r})"


SI = UnitSystem("SI", 1.0, 1.0, 1.0, G_SI)
ASTRONOMICAL = UnitSystem("astronomical", AU, SOLAR_MASS, GAUSSIAN_YEAR, 4.0 * math.pi**2)
DIMENSIONLESS = UnitSystem("dimensionless", None, None, None, 1.0)


def convert(value, quantity, source, target):
    """Converts ``value`` of ``quantity`` from unit system ``source`` to ``target``."""
    return target.from_si(source.to_si(value, quantity), quantity)


def scaled(name, length, mass, time):
    """A unit system with the given units (in SI); ``G`` follows from them."""
    G = G_SI * mass * time**2 / length**3
    return UnitSystem(name, float(length), float(mass), float(time), G)
