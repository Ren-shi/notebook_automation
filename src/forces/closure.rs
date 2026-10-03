//! Forces defined by Rust closures, for new physics without a dedicated struct.

use super::Force;
use crate::error::Result;
use crate::parallel;
use crate::vec3::Vec3;

type AccelFn = dyn Fn(f64, &[Vec3], &[Vec3], &[f64], &mut [Vec3]) -> Result<()> + Send + Sync;
type PotentialFn = dyn Fn(f64, &[Vec3], &[f64]) -> Result<f64> + Send + Sync;

/// A force given by closures. [`ClosureForce::new`] takes the whole-system form (add to
/// `acc` like [`Force::accumulate`]); [`ClosureForce::per_particle`] takes the acceleration
/// of one particle from its own state, evaluated in parallel for large systems.
///
/// ```
/// use physim::{ClosureForce, Vec3};
/// // A constant field pulling everything towards -z, with its potential.
/// let f = ClosureForce::per_particle("Down", |_t, _r, _v, _m| Vec3::new(0.0, 0.0, -1.0))
///     .with_potential(|_t, pos, mass| Ok(pos.iter().zip(mass).map(|(r, m)| m * r.z).sum()))
///     .velocity_independent();
/// ```
///
/// Closures must be pure functions of their arguments (see [`Force`]).
pub struct ClosureForce {
    name: String,
    accel: Box<AccelFn>,
    potential: Option<Box<PotentialFn>>,
    velocity_dependent: bool,
}

impl ClosureForce {
    /// `accel(t, pos, vel, mass, acc)` adds every particle's acceleration to `acc`.
    /// Assumed velocity-dependent until [`ClosureForce::velocity_independent`] says otherwise.
    pub fn new(
        name: impl Into<String>,
        accel: impl Fn(f64, &[Vec3], &[Vec3], &[f64], &mut [Vec3]) -> Result<()> + Send + Sync + 'static,
    ) -> Self {
        Self {
            name: name.into(),
            accel: Box::new(accel),
            potential: None,
            velocity_dependent: true,
        }
    }

    /// `accel(t, r, v, m)` returns the acceleration of one particle with position `r`,
    /// velocity `v` and mass `m`.
    pub fn per_particle(
        name: impl Into<String>,
        accel: impl Fn(f64, Vec3, Vec3, f64) -> Vec3 + Send + Sync + 'static,
    ) -> Self {
        Self::new(name, move |t, pos, vel, mass, acc| {
            parallel::for_each_indexed(acc, |i, a| *a += accel(t, pos[i], vel[i], mass[i]));
            Ok(())
        })
    }

    /// Adds a potential energy `potential(t, pos, mass)`, making the force conservative
    /// in energy diagnostics.
    pub fn with_potential(
        mut self,
        potential: impl Fn(f64, &[Vec3], &[f64]) -> Result<f64> + Send + Sync + 'static,
    ) -> Self {
        self.potential = Some(Box::new(potential));
        self
    }

    /// Declares that the acceleration ignores velocities, so symplectic integrators treat
    /// the force exactly and may reuse its evaluations.
    pub fn velocity_independent(mut self) -> Self {
        self.velocity_dependent = false;
        self
    }
}

impl Force for ClosureForce {
    fn accumulate(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        (self.accel)(t, pos, vel, mass, acc)
    }

    fn potential(&self, t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        self.potential.as_ref().map(|p| p(t, pos, mass)).transpose()
    }

    fn velocity_dependent(&self) -> bool {
        self.velocity_dependent
    }

    fn name(&self) -> String {
        self.name.clone()
    }
}
