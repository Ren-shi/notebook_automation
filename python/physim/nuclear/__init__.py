"""Planning nuclear physics experiments: the experiment setup (beam, target, detectors, run conditions).

This is the first part of the nuclear experiment planner (backlog items 33-43). Kinematics, cross sections,
energy loss and count rates build on the :class:`Experiment` described here.
"""

from . import (coulex, data, detectors, events, export, gamma, guide, kinematics, planner, rates, report,
               rootio, rutherford, stopping, validation)
from .experiment import (
    EXAMPLES,
    REACTIONS,
    SCHEMA,
    SHAPES,
    Beam,
    Detector,
    Excitation,
    Experiment,
    GammaDetector,
    Layer,
    Run,
    SetupError,
    Target,
    example_names,
)
from .names import parse_material, parse_nuclide
from .quantity import UNITS, Quantity

__all__ = [
    "EXAMPLES",
    "example_names",
    "coulex",
    "data",
    "detectors",
    "events",
    "export",
    "gamma",
    "guide",
    "kinematics",
    "planner",
    "rates",
    "report",
    "rootio",
    "rutherford",
    "stopping",
    "validation",
    "REACTIONS",
    "SCHEMA",
    "SHAPES",
    "UNITS",
    "Beam",
    "Detector",
    "Excitation",
    "Experiment",
    "GammaDetector",
    "Layer",
    "Quantity",
    "Run",
    "SetupError",
    "Target",
    "parse_material",
    "parse_nuclide",
]
