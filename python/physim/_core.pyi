from typing import Any, Callable, Optional, Sequence, TypedDict, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

INTEGRATORS: list[str]
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

Force = UniformField | NewtonianGravity | Spring | AnchorSpring | LinearDrag | QuadraticDrag | CustomForce

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
    @property
    def pinned(self) -> NDArray[np.bool_]: ...
    @property
    def n_particles(self) -> int: ...
    @property
    def forces(self) -> dict[int, str]: ...
    @property
    def integrator_info(self) -> IntegratorInfo: ...
    def __init__(self, integrator: str = "verlet") -> None: ...
    def add_particle(self, pos: ArrayLike, vel: Optional[ArrayLike] = None, mass: float = 1.0) -> int: ...
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
