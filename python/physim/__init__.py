"""physim: classical-mechanics simulation with a Rust engine.

Quick start::

    import physim as ps

    w = ps.World(integrator="yoshida4")
    w.add_particle([1, 0, 0], vel=[0, 1, 0], mass=1e-3)
    w.add_particle([0, 0, 0], mass=1.0)
    w.add_force(ps.NewtonianGravity(G=1.0))
    traj = w.run(dt=1e-3, steps=10_000, record_every=10)
    traj.pos.shape  # (frames, particles, 3)
"""

from ._core import (
    INTEGRATORS,
    AnchorSpring,
    CustomForce,
    Event,
    LinearDrag,
    NewtonianGravity,
    QuadraticDrag,
    Spring,
    Trajectory,
    UniformField,
    World,
)
from .analysis import relative_energy_error

__all__ = [
    "INTEGRATORS",
    "AnchorSpring",
    "CustomForce",
    "Event",
    "LinearDrag",
    "NewtonianGravity",
    "QuadraticDrag",
    "Spring",
    "Trajectory",
    "UniformField",
    "World",
    "relative_energy_error",
]
