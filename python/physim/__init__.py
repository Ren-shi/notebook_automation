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
    DampedSpring,
    Event,
    HarmonicTrap,
    HernquistPotential,
    J2Oblateness,
    LinearDrag,
    ModulatedSpring,
    NewtonianGravity,
    PeriodicForce,
    PlummerPotential,
    PostNewtonian,
    PowerLaw,
    QuadraticDrag,
    Spring,
    SpringNetwork,
    Trajectory,
    UniformField,
    World,
    Yukawa,
    __version__,
)
from .analysis import relative_energy_error
from .io import (
    TrajectoryWriter,
    load_checkpoint,
    load_trajectory,
    save_checkpoint,
    save_trajectory,
)

__all__ = [
    "INTEGRATORS",
    "AnchorSpring",
    "CustomForce",
    "DampedSpring",
    "Event",
    "HarmonicTrap",
    "HernquistPotential",
    "J2Oblateness",
    "LinearDrag",
    "ModulatedSpring",
    "NewtonianGravity",
    "PeriodicForce",
    "PlummerPotential",
    "PostNewtonian",
    "PowerLaw",
    "QuadraticDrag",
    "Spring",
    "SpringNetwork",
    "Trajectory",
    "TrajectoryWriter",
    "UniformField",
    "World",
    "Yukawa",
    "load_checkpoint",
    "load_trajectory",
    "relative_energy_error",
    "save_checkpoint",
    "save_trajectory",
    "__version__",
]
