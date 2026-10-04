from typing import Any, Callable, Optional, Sequence, TypedDict, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

INTEGRATORS: list[str]
RANDOM_SETUP_COUNTER: int

def random_uniform(seed: int, counter: int, n: int, start: int = 0) -> NDArray[np.float64]: ...
def random_normal(seed: int, counter: int, n: int, start: int = 0) -> NDArray[np.float64]: ...
__version__: str

class Trajectory:
    t: NDArray[np.float64]  # (F,)
    pos: NDArray[np.float64]  # (F, N, 3)
    vel: NDArray[np.float64]  # (F, N, 3)
    kinetic: Optional[NDArray[np.float64]]  # (F,), None if the run skipped energies
    potential: Optional[NDArray[np.float64]]  # (F,)
    energy: Optional[NDArray[np.float64]]  # (F,)
    tension: NDArray[np.float64]  # (F, C) constraint tensions
    n_particles: int
    n_frames: int
    event_t: NDArray[np.float64]  # (K,)
    event_index: NDArray[np.int64]  # (K,) index into the run's events list
    event_pos: NDArray[np.float64]  # (K, N, 3)
    event_vel: NDArray[np.float64]  # (K, N, 3)
    terminated_by: Optional[int]
    metadata: dict[str, Any]
    def __init__(
        self,
        t: ArrayLike,
        pos: ArrayLike,
        vel: ArrayLike,
        kinetic: Optional[ArrayLike] = None,
        potential: Optional[ArrayLike] = None,
        *,
        tension: Optional[ArrayLike] = None,
        event_t: Optional[ArrayLike] = None,
        event_index: Optional[Sequence[int]] = None,
        event_pos: Optional[ArrayLike] = None,
        event_vel: Optional[ArrayLike] = None,
        terminated_by: Optional[int] = None,
        metadata: Optional[dict[str, Any]] = None,
    ) -> None: ...
    def __len__(self) -> int: ...

class UniformField:
    g: NDArray[np.float64]
    def __init__(self, g: ArrayLike) -> None: ...

class NewtonianGravity:
    G: float
    softening: float
    def __init__(self, G: float = 1.0, softening: float = 0.0) -> None: ...

class Spring:
    def __init__(self, i: int, j: int, k: float, rest_length: float) -> None: ...

class AnchorSpring:
    def __init__(self, i: int, anchor: ArrayLike, k: float, rest_length: float = 0.0) -> None: ...

class LinearDrag:
    def __init__(self, gamma: float) -> None: ...

class QuadraticDrag:
    def __init__(self, c: float) -> None: ...

class DampedSpring:
    def __init__(self, i: int, j: int, k: float, rest_length: float, c: float) -> None: ...

class ModulatedSpring:
    def __init__(
        self,
        i: int,
        to: Union[int, ArrayLike],
        k: float,
        depth: float,
        omega: float,
        phase: float = 0.0,
        rest_length: float = 0.0,
    ) -> None: ...

class SpringNetwork:
    def __init__(
        self,
        i: Sequence[int],
        j: Sequence[int],
        k: Union[float, ArrayLike],
        rest_length: Union[float, ArrayLike],
    ) -> None: ...
    def __len__(self) -> int: ...

class PowerLaw:
    def __init__(self, k: float, n: float, center: Optional[ArrayLike] = None) -> None: ...

class Yukawa:
    def __init__(self, k: float, length: float, center: Optional[ArrayLike] = None) -> None: ...

class PlummerPotential:
    def __init__(self, GM: float, a: float, center: Optional[ArrayLike] = None) -> None: ...

class HernquistPotential:
    def __init__(self, GM: float, a: float, center: Optional[ArrayLike] = None) -> None: ...

class HarmonicTrap:
    def __init__(self, omega: Union[float, ArrayLike], center: Optional[ArrayLike] = None) -> None: ...

class PeriodicForce:
    def __init__(self, i: int, amplitude: ArrayLike, omega: float, phase: float = 0.0) -> None: ...

class PostNewtonian:
    def __init__(self, central: int, c: float, G: float = 1.0) -> None: ...

class J2Oblateness:
    def __init__(
        self,
        central: int,
        J2: float,
        radius: float,
        axis: Optional[ArrayLike] = None,
        G: float = 1.0,
    ) -> None: ...

class HenonHeiles:
    def __init__(self, lam: float = 1.0, center: Optional[ArrayLike] = None) -> None: ...

class ElectricField:
    def __init__(self, E: ArrayLike) -> None: ...

class MagneticField:
    def __init__(self, B: ArrayLike) -> None: ...

class Coulomb:
    def __init__(self, k: float = 1.0, softening: float = 0.0) -> None: ...

class FieldForce:
    def __init__(
        self,
        E: Optional[Callable[[float, NDArray[np.float64]], ArrayLike]] = None,
        B: Optional[Callable[[float, NDArray[np.float64]], ArrayLike]] = None,
        name: str = "FieldForce",
    ) -> None: ...

class TreeGravity:
    def __init__(
        self, G: float = 1.0, softening: float = 0.0, theta: float = 0.5, quadrupole: bool = False
    ) -> None: ...

class SoftContact:
    def __init__(self, k: float, damping: float = 0.0, law: str = "linear") -> None: ...

class LennardJones:
    def __init__(
        self,
        epsilon: float = 1.0,
        sigma: float = 1.0,
        cutoff: float = 2.5,
        shift: bool = True,
        box: Optional[Union[float, ArrayLike]] = None,
    ) -> None: ...

class Morse:
    def __init__(
        self,
        depth: float,
        a: float,
        r0: float,
        cutoff: float,
        shift: bool = True,
        box: Optional[Union[float, ArrayLike]] = None,
    ) -> None: ...

class TabulatedPair:
    def __init__(
        self,
        r_min: float,
        dr: float,
        V: Sequence[float],
        dV: Sequence[float],
        cutoff: float,
        shift: bool = True,
        box: Optional[Union[float, ArrayLike]] = None,
    ) -> None: ...

class CustomForce:
    def __init__(
        self,
        acceleration: Callable[[float, NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]], ArrayLike],
        potential: Optional[Callable[[float, NDArray[np.float64], NDArray[np.float64]], float]] = None,
        velocity_dependent: bool = True,
        name: Optional[str] = None,
    ) -> None: ...

class Event:
    name: str
    direction: int
    terminal: bool
    def __init__(
        self,
        function: Callable[[float, NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]], float],
        direction: int = 0,
        terminal: bool = False,
        name: Optional[str] = None,
    ) -> None: ...
    @staticmethod
    def radial_velocity(i: int, j: Optional[int] = None, direction: int = 0, terminal: bool = False) -> Event: ...
    @staticmethod
    def coordinate(i: int, axis: int | str, value: float = 0.0, direction: int = 0, terminal: bool = False) -> Event: ...
    @staticmethod
    def separation(i: int, j: int, distance: float, direction: int = 0, terminal: bool = False) -> Event: ...

Force = Union[
    UniformField,
    NewtonianGravity,
    TreeGravity,
    ParticleMesh,
    SoftContact,
    LennardJones,
    Morse,
    TabulatedPair,
    Spring,
    AnchorSpring,
    LinearDrag,
    QuadraticDrag,
    DampedSpring,
    ModulatedSpring,
    SpringNetwork,
    PowerLaw,
    Yukawa,
    PlummerPotential,
    HernquistPotential,
    HarmonicTrap,
    HenonHeiles,
    PeriodicForce,
    PostNewtonian,
    J2Oblateness,
    ElectricField,
    MagneticField,
    Coulomb,
    FieldForce,
    CustomForce,
]

class IntegratorInfo(TypedDict):
    name: str
    order: int
    symplectic: bool

class World:
    integrator: str
    t: float
    positions: NDArray[np.float64]
    velocities: NDArray[np.float64]
    masses: NDArray[np.float64]
    charges: NDArray[np.float64]
    radii: NDArray[np.float64]
    @property
    def pinned(self) -> NDArray[np.bool_]: ...
    @property
    def n_particles(self) -> int: ...
    @property
    def forces(self) -> dict[int, str]: ...
    @property
    def integrator_info(self) -> IntegratorInfo: ...
    def __init__(self, integrator: str = "verlet") -> None: ...
    def add_particle(
        self, pos: ArrayLike, vel: Optional[ArrayLike] = None, mass: float = 1.0,
        charge: float = 0.0,
        radius: float = 0.0,
    ) -> int: ...
    def remove_particle(self, i: int) -> None: ...
    def pin(self, i: int, pinned: bool = True) -> None: ...
    def add_force(self, force: Force) -> int: ...
    def remove_force(self, id: int) -> None: ...
    def replace_force(self, id: int, force: Force) -> None: ...
    def force_params(self, id: int) -> dict[str, float | NDArray[np.float64]]: ...
    def set_force_params(self, id: int, **params: float | ArrayLike) -> None: ...
    def clear_forces(self) -> None: ...
    def step(self, dt: float, n: int = 1) -> None: ...
    def run(
        self,
        dt: float,
        steps: int,
        record_every: int = 1,
        events: Optional[Sequence[Event]] = None,
        *,
        energies: bool = True,
        sink: Optional[Callable[[Trajectory], Any]] = None,
        chunk_size: int = 1024,
    ) -> Trajectory: ...
    def run_adaptive(
        self,
        t_end: float,
        *,
        rtol: float = 1e-9,
        atol: float = 1e-12,
        times: Optional[ArrayLike] = None,
        events: Optional[Sequence[Event]] = None,
        energies: bool = True,
        first_step: Optional[float] = None,
        max_step: Optional[float] = None,
        max_steps: int = 10_000_000,
        sink: Optional[Callable[[Trajectory], Any]] = None,
        chunk_size: int = 1024,
    ) -> Trajectory: ...
    def lyapunov(
        self,
        dt: float,
        steps: int,
        n: int = 1,
        *,
        renormalize_every: int = 1,
        record_every: int = 100,
        seed: int = 1,
    ) -> dict[str, NDArray[np.float64]]: ...
    def set_collisions(
        self,
        restitution: float = 1.0,
        walls: Optional[Sequence[tuple[ArrayLike, float]]] = None,
        *,
        between_particles: bool = True,
        max_per_step: int = 100_000,
    ) -> None: ...
    def clear_collisions(self) -> None: ...
    @property
    def collisions(self) -> Optional[dict[str, Any]]: ...
    @property
    def collision_count(self) -> int: ...
    def use_langevin(self, temperature: float, friction: float = 1.0, seed: Optional[int] = None) -> None: ...
    seed: int
    def use_nose_hoover(self, temperature: float, tau: float = 1.0) -> None: ...
    def temperature(self) -> float: ...
    def pressure(self, volume: float) -> float: ...
    def thermostat_energy(self) -> float: ...
    def use_composition(self, weights: Sequence[float], order: int, name: str = "composition") -> None: ...
    def use_splitting(
        self, ops: Sequence[tuple[str, float]], order: int, name: str = "splitting"
    ) -> None: ...
    def add_rod(self, i: int, to: Union[int, ArrayLike], length: Optional[float] = None) -> int: ...
    def remove_constraint(self, id: int) -> None: ...
    def clear_constraints(self) -> None: ...
    @property
    def constraints(self) -> dict[int, str]: ...
    def constraint_tensions(self) -> NDArray[np.float64]: ...
    constraint_tolerance: float
    def checkpoint(self) -> dict[str, Any]: ...
    @staticmethod
    def from_checkpoint(
        checkpoint: dict[str, Any], custom_forces: Optional[dict[Any, Any]] = None
    ) -> World: ...
    def accelerations(self) -> NDArray[np.float64]: ...
    def kinetic_energy(self) -> float: ...
    def potential_energy(self) -> float: ...
    def total_energy(self) -> float: ...
    def momentum(self) -> NDArray[np.float64]: ...
    def angular_momentum(self) -> NDArray[np.float64]: ...
    def center_of_mass(self) -> NDArray[np.float64]: ...

class BodyGravity:
    def __init__(self, g: ArrayLike) -> None: ...
    @property
    def g(self) -> NDArray[np.float64]: ...

class BodySpring:
    def __init__(
        self,
        a: Union[tuple[int, ArrayLike], ArrayLike],
        b: Union[tuple[int, ArrayLike], ArrayLike],
        k: float,
        rest_length: float = 0.0,
    ) -> None: ...

class RigidTrajectory:
    t: NDArray[np.float64]
    pos: NDArray[np.float64]
    vel: NDArray[np.float64]
    orientation: NDArray[np.float64]
    rotation: NDArray[np.float64]
    angular_momentum: NDArray[np.float64]
    angular_velocity: NDArray[np.float64]
    kinetic: NDArray[np.float64]
    potential: NDArray[np.float64]
    energy: NDArray[np.float64]
    n_bodies: int
    n_frames: int
    def __len__(self) -> int: ...

class RigidSystem:
    t: float
    def __init__(self, order: int = 2) -> None: ...
    def add_body(
        self,
        mass: float,
        inertia: ArrayLike,
        position: Optional[ArrayLike] = None,
        velocity: Optional[ArrayLike] = None,
        orientation: Optional[ArrayLike] = None,
        angular_velocity: Optional[ArrayLike] = None,
        angular_momentum: Optional[ArrayLike] = None,
        pivot: Optional[ArrayLike] = None,
        center_of_mass: Optional[ArrayLike] = None,
    ) -> int: ...
    def set_body_state(
        self,
        i: int,
        position: Optional[ArrayLike] = None,
        velocity: Optional[ArrayLike] = None,
        orientation: Optional[ArrayLike] = None,
        angular_velocity: Optional[ArrayLike] = None,
        angular_momentum: Optional[ArrayLike] = None,
    ) -> None: ...
    def add_force(
        self,
        force: Union[
            BodyGravity,
            BodySpring,
            Callable[[float, NDArray[np.float64], NDArray[np.float64]], tuple[ArrayLike, ArrayLike]],
        ],
    ) -> None: ...
    def step(self, dt: float) -> None: ...
    def run(self, dt: float, steps: int, record_every: int = 1) -> RigidTrajectory: ...
    def wrenches(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...
    def kinetic_energy(self) -> float: ...
    def potential_energy(self) -> float: ...
    def total_energy(self) -> float: ...
    def momentum(self) -> NDArray[np.float64]: ...
    def angular_momentum(self) -> NDArray[np.float64]: ...
    @property
    def order(self) -> int: ...
    @property
    def n_bodies(self) -> int: ...
    def __len__(self) -> int: ...
    @property
    def masses(self) -> NDArray[np.float64]: ...
    @property
    def inertia(self) -> NDArray[np.float64]: ...
    @property
    def positions(self) -> NDArray[np.float64]: ...
    @property
    def velocities(self) -> NDArray[np.float64]: ...
    @property
    def centers_of_mass(self) -> NDArray[np.float64]: ...
    @property
    def orientations(self) -> NDArray[np.float64]: ...
    @property
    def rotation_matrices(self) -> NDArray[np.float64]: ...
    @property
    def angular_momenta(self) -> NDArray[np.float64]: ...
    @property
    def angular_velocities(self) -> NDArray[np.float64]: ...

class ParticleMesh:
    def __init__(
        self,
        box_size: float,
        cells: int = 64,
        G: float = 1.0,
        center: Optional[ArrayLike] = None,
        periodic: bool = False,
        softening: Optional[float] = None,
    ) -> None: ...
    @property
    def spacing(self) -> float: ...

class WaveEquation:
    u: NDArray[np.float64]
    v: NDArray[np.float64]
    t: float
    def __init__(self, shape: Sequence[int], spacing: float, c: float = 1.0, boundary: str = "dirichlet") -> None: ...
    @property
    def axes(self) -> list[NDArray[np.float64]]: ...
    @property
    def max_stable_dt(self) -> float: ...
    def step(self, dt: float, n: int = 1) -> None: ...
    def energy(self) -> float: ...

class HeatEquation:
    u: NDArray[np.float64]
    t: float
    def __init__(
        self,
        shape: Sequence[int],
        spacing: float,
        diffusivity: float = 1.0,
        boundary: str = "neumann",
        method: str = "crank_nicolson",
    ) -> None: ...
    @property
    def axes(self) -> list[NDArray[np.float64]]: ...
    @property
    def max_explicit_dt(self) -> float: ...
    @property
    def last_iterations(self) -> int: ...
    def step(self, dt: float, n: int = 1) -> None: ...
    def total(self) -> float: ...

def solve_poisson(
    f: ArrayLike, spacing: float, boundary: str = "periodic", phi: Optional[ArrayLike] = None
) -> NDArray[np.float64]: ...

class Schrodinger:
    psi: NDArray[np.complex128]
    potential: NDArray[np.float64]  # assign an array, a function V(t, x[, y[, z]]) or None
    absorber: NDArray[np.float64]
    t: float
    order: int
    def __init__(
        self,
        shape: Sequence[int],
        spacing: float,
        hbar: float = 1.0,
        mass: float = 1.0,
        boundary: str = "periodic",
        method: Optional[str] = None,
        order: int = 2,
        origin: Optional[Union[float, ArrayLike]] = None,
        potential: Optional[Union[ArrayLike, Callable[..., ArrayLike]]] = None,
    ) -> None: ...
    @property
    def density(self) -> NDArray[np.float64]: ...
    @property
    def potential_function(self) -> Optional[Callable[..., ArrayLike]]: ...
    @property
    def method(self) -> str: ...
    @property
    def boundary(self) -> str: ...
    @property
    def shape(self) -> list[int]: ...
    @property
    def spacing(self) -> float: ...
    @property
    def hbar(self) -> float: ...
    @property
    def mass(self) -> float: ...
    @property
    def cell_volume(self) -> float: ...
    @property
    def axes(self) -> list[NDArray[np.float64]]: ...
    @property
    def last_iterations(self) -> int: ...
    def set_gaussian(
        self,
        center: Optional[Union[float, ArrayLike]] = None,
        width: Optional[Union[float, ArrayLike]] = None,
        momentum: Optional[Union[float, ArrayLike]] = None,
    ) -> None: ...
    def absorbing_layer(self, width: float, strength: float) -> None: ...
    def step(self, dt: float, n: int = 1) -> None: ...
    def run(
        self, dt: float, steps: int, record_every: int = 1
    ) -> tuple[NDArray[np.float64], NDArray[np.complex128]]: ...
    def normalize(self) -> None: ...
    def norm(self) -> float: ...
    def position(self) -> NDArray[np.float64]: ...
    def position_variance(self) -> NDArray[np.float64]: ...
    def momentum(self) -> NDArray[np.float64]: ...
    def kinetic_energy(self) -> float: ...
    def potential_energy(self) -> float: ...
    def energy(self) -> float: ...
    def eigenstates(
        self, count: int = 1, tol: float = 1e-9, max_iter: int = 1000
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...
