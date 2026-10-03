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
    BodyGravity,
    BodySpring,
    Coulomb,
    CustomForce,
    DampedSpring,
    ElectricField,
    Event,
    FieldForce,
    HarmonicTrap,
    HenonHeiles,
    HernquistPotential,
    J2Oblateness,
    LennardJones,
    LinearDrag,
    MagneticField,
    ModulatedSpring,
    Morse,
    NewtonianGravity,
    PeriodicForce,
    PlummerPotential,
    PostNewtonian,
    PowerLaw,
    QuadraticDrag,
    RigidSystem,
    RigidTrajectory,
    SoftContact,
    Spring,
    SpringNetwork,
    TabulatedPair,
    Trajectory,
    TreeGravity,
    UniformField,
    World,
    Yukawa,
    __version__,
)
from .analysis import (
    poincare_section,
    radial_distribution,
    relative_energy_error,
    tabulate_pair,
)
from .geometry import (
    box_walls,
    inertia_box,
    inertia_cylinder,
    inertia_ellipsoid,
    inertia_sphere,
    quaternion_from_axis_angle,
)
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
    "BodyGravity",
    "BodySpring",
    "Coulomb",
    "CustomForce",
    "DampedSpring",
    "ElectricField",
    "Event",
    "FieldForce",
    "HarmonicTrap",
    "HenonHeiles",
    "HernquistPotential",
    "J2Oblateness",
    "LennardJones",
    "LinearDrag",
    "MagneticField",
    "ModulatedSpring",
    "Morse",
    "NewtonianGravity",
    "PeriodicForce",
    "PlummerPotential",
    "PostNewtonian",
    "PowerLaw",
    "QuadraticDrag",
    "RigidSystem",
    "RigidTrajectory",
    "SoftContact",
    "Spring",
    "SpringNetwork",
    "TabulatedPair",
    "Trajectory",
    "TrajectoryWriter",
    "TreeGravity",
    "UniformField",
    "World",
    "Yukawa",
    "box_walls",
    "inertia_box",
    "inertia_cylinder",
    "inertia_ellipsoid",
    "inertia_sphere",
    "load_checkpoint",
    "load_trajectory",
    "poincare_section",
    "quaternion_from_axis_angle",
    "radial_distribution",
    "relative_energy_error",
    "save_checkpoint",
    "save_trajectory",
    "tabulate_pair",
    "__version__",
]
