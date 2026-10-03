//! `physim`: a classical-mechanics engine for point particles.
//!
//! The pieces are deliberately small and swappable:
//! - [`State`]: positions, velocities, masses and time.
//! - [`Force`]: anything that adds accelerations (and optionally a potential).
//! - [`Integrator`]: a time-stepping scheme.
//! - [`World`]: ties them together and records [`Trajectory`]s, with fixed steps
//!   ([`World::run`]) or adaptive ones ([`World::run_adaptive`]).

pub mod adaptive;
pub mod chaos;
pub mod checkpoint;
pub mod constraints;
pub mod error;
pub mod events;
pub mod forces;
pub mod integrators;
mod parallel;
pub mod state;
pub mod vec3;
pub mod world;

#[cfg(feature = "python")]
mod python;

pub use adaptive::{AdaptiveOptions, AdaptiveOutcome, AdaptiveRun, AdaptiveStats, Output};
pub use chaos::{LyapunovOptions, LyapunovRun};
pub use checkpoint::{Checkpoint, SavedForce};
pub use constraints::{Anchor, ConstraintId, Constraints, Rod};
pub use error::{Result, SimError};
pub use events::{Direction, Event, EventFunction, EventHit};
pub use forces::{
    AnchorSpring, BuiltinForce, ClosureForce, Coulomb, DampedSpring, ElectricField, FieldFunctions,
    Force, ForceId, ForceSet, HarmonicTrap, HenonHeiles, HernquistPotential, J2Oblateness,
    LinearDrag, MagneticField, ModulatedSpring, NewtonianGravity, Param, PeriodicForce,
    PlummerPotential, PostNewtonian, PowerLaw, QuadraticDrag, Spring, SpringNetwork, TreeGravity,
    UniformField, Yukawa,
};
pub use integrators::Integrator;
pub use state::State;
pub use vec3::Vec3;
pub use world::{Frame, Recorder, RunFailure, RunOptions, Trajectory, World};
