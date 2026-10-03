//! `physim`: a classical-mechanics engine for point particles.
//!
//! The pieces are deliberately small and swappable:
//! - [`State`]: positions, velocities, masses and time.
//! - [`Force`]: anything that adds accelerations (and optionally a potential).
//! - [`Integrator`]: a time-stepping scheme.
//! - [`World`]: ties them together and records [`Trajectory`]s.

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

pub use checkpoint::{Checkpoint, SavedForce};
pub use constraints::{Anchor, ConstraintId, Constraints, Rod};
pub use error::{Result, SimError};
pub use events::{Direction, Event, EventFunction, EventHit};
pub use forces::{
    AnchorSpring, BuiltinForce, Force, ForceId, ForceSet, LinearDrag, NewtonianGravity, Param,
    QuadraticDrag, Spring, UniformField,
};
pub use integrators::Integrator;
pub use state::State;
pub use vec3::Vec3;
pub use world::{Frame, Recorder, RunFailure, RunOptions, Trajectory, World};
