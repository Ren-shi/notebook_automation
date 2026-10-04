"""Planning nuclear physics experiments: the experiment setup (beam, target, detectors, run conditions).

This is the first part of the nuclear experiment planner (backlog items 33-43). Kinematics, cross sections,
energy loss and count rates build on the :class:`Experiment` described here.
"""

from .experiment import (
    REACTIONS,
    SCHEMA,
    SHAPES,
    Beam,
    Detector,
    Experiment,
    Layer,
    Run,
    SetupError,
    Target,
)
from .names import parse_material, parse_nuclide
from .quantity import UNITS, Quantity

__all__ = [
    "REACTIONS",
    "SCHEMA",
    "SHAPES",
    "UNITS",
    "Beam",
    "Detector",
    "Experiment",
    "Layer",
    "Quantity",
    "Run",
    "SetupError",
    "Target",
    "parse_material",
    "parse_nuclide",
]
