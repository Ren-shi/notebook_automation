//! `physim`: a classical-mechanics engine for point particles.
//!
//! The pieces are deliberately small and swappable:
//! - [`State`]: positions, velocities, masses and time.
//! - [`Force`]: anything that adds accelerations (and optionally a potential).
//! - [`Integrator`]: a time-stepping scheme.
//! - [`World`]: ties them together and records [`Trajectory`]s.

pub mod error;
pub mod forces;
pub mod integrators;
pub mod state;
pub mod vec3;
pub mod world;

#[cfg(feature = "python")]
mod python;

pub use error::{Result, SimError};
pub use forces::{
    AnchorSpring, Force, ForceId, ForceSet, LinearDrag, NewtonianGravity, Param, QuadraticDrag,
    Spring, UniformField,
};
pub use integrators::Integrator;
pub use state::State;
pub use vec3::Vec3;
pub use world::{RunFailure, Trajectory, World};
