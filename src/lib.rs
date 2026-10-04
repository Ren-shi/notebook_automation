//! `physim`: a classical-mechanics engine for point particles, with grid solvers for fields
//! and the Schrödinger equation.
//!
//! The pieces are deliberately small and swappable:
//! - [`State`]: positions, velocities, masses and time.
//! - [`Force`]: anything that adds accelerations (and optionally a potential).
//! - [`Integrator`]: a time-stepping scheme.
//! - [`fields`]: grids, the wave, heat, Poisson and Schrödinger equations, and particle-mesh
//!   gravity.
//! - [`RigidSystem`]: rigid bodies with orientation, torques and a symplectic rotation step.
//! - [`World`]: ties them together and records [`Trajectory`]s, with fixed steps
//!   ([`World::run`]) or adaptive ones ([`World::run_adaptive`]).

pub mod adaptive;
mod broadphase;
pub mod chaos;
pub mod checkpoint;
pub mod collisions;
pub mod constraints;
pub mod error;
pub mod events;
pub mod fields;
pub mod forces;
pub mod integrators;
mod parallel;
pub mod rigid;
pub mod rng;
pub mod state;
pub mod vec3;
pub mod world;

#[cfg(feature = "python")]
mod python;

pub use adaptive::{AdaptiveOptions, AdaptiveOutcome, AdaptiveRun, AdaptiveStats, Output};
pub use chaos::{LyapunovOptions, LyapunovRun};
pub use checkpoint::{Checkpoint, SavedForce};
pub use collisions::{Collisions, Wall};
pub use constraints::{Anchor, ConstraintId, Constraints, Rod};
pub use error::{Result, SimError};
pub use events::{Direction, Event, EventFunction, EventHit};
pub use fields::pm::ParticleMesh;
pub use fields::quantum::{PotentialFn, QuantumMethod, Schrodinger};
pub use fields::solvers::{poisson, Heat, HeatMethod, Wave};
pub use fields::{Boundary, Grid};
pub use forces::{
    AnchorSpring, BuiltinForce, ClosureForce, ContactLaw, Coulomb, DampedSpring, ElectricField,
    FieldFunctions, Force, ForceId, ForceSet, HarmonicTrap, HenonHeiles, HernquistPotential,
    J2Oblateness, LinearDrag, MagneticField, ModulatedSpring, NewtonianGravity, PairKind,
    PairPotential, Param, PeriodicForce, PlummerPotential, PostNewtonian, PowerLaw, QuadraticDrag,
    SoftContact, Spring, SpringNetwork, TreeGravity, UniformField, Yukawa,
};
pub use integrators::Integrator;
pub use rigid::{
    Attachment, BodyClosure, BodyForce, BodyGravity, BodySpring, Quat, RigidBody, RigidSystem,
    RigidTrajectory, Wrench,
};
pub use state::State;
pub use vec3::Vec3;
pub use world::{Frame, Recorder, RunFailure, RunOptions, Trajectory, World};
