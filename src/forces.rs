//! Forces act on the whole system at once: each one adds its contribution to the
//! per-particle accelerations. Implement [`Force`] to add new physics in Rust.

use crate::error::{invalid, Result};
use crate::vec3::Vec3;

pub trait Force: Send + Sync {
    /// Add this force's acceleration contribution for every particle to `acc`.
    /// `acc` is pre-sized to the number of particles and must only be added to.
    fn accumulate(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()>;

    /// Potential energy of the configuration, or `None` for non-conservative forces.
    fn potential(&self, _t: f64, _pos: &[Vec3], _mass: &[f64]) -> Result<Option<f64>> {
        Ok(None)
    }

    /// Whether the acceleration depends on velocity. Integrators that split position
    /// and velocity updates (the symplectic ones) are only exact for `false`.
    fn velocity_dependent(&self) -> bool {
        false
    }

    fn name(&self) -> String;
}

/// The collection of forces acting on a system.
#[derive(Default)]
pub struct ForceSet {
    forces: Vec<Box<dyn Force>>,
}

impl ForceSet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, force: Box<dyn Force>) {
        self.forces.push(force);
    }

    pub fn clear(&mut self) {
        self.forces.clear();
    }

    pub fn iter(&self) -> impl Iterator<Item = &dyn Force> {
        self.forces.iter().map(|f| f.as_ref())
    }

    pub fn len(&self) -> usize {
        self.forces.len()
    }

    pub fn is_empty(&self) -> bool {
        self.forces.is_empty()
    }

    /// Overwrites `acc` with the total acceleration of every particle.
    pub fn accelerations(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut Vec<Vec3>,
    ) -> Result<()> {
        acc.clear();
        acc.resize(pos.len(), Vec3::ZERO);
        for f in &self.forces {
            f.accumulate(t, pos, vel, mass, acc)?;
        }
        Ok(())
    }

    /// Sum of the potentials of all conservative forces (non-conservative ones are skipped).
    pub fn potential(&self, t: f64, pos: &[Vec3], mass: &[f64]) -> Result<f64> {
        let mut total = 0.0;
        for f in &self.forces {
            total += f.potential(t, pos, mass)?.unwrap_or(0.0);
        }
        Ok(total)
    }
}

fn check_index(i: usize, n: usize, what: &str) -> Result<()> {
    if i >= n {
        return invalid(format!(
            "{what}: particle index {i} out of range (have {n} particles)"
        ));
    }
    Ok(())
}

/// Constant acceleration field, e.g. surface gravity `g = (0, -9.81, 0)`.
#[derive(Debug, Clone)]
pub struct UniformField {
    pub g: Vec3,
}

impl Force for UniformField {
    fn accumulate(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        for a in acc {
            *a += self.g;
        }
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        Ok(Some(
            pos.iter().zip(mass).map(|(r, m)| -m * self.g.dot(*r)).sum(),
        ))
    }

    fn name(&self) -> String {
        "UniformField".into()
    }
}

/// Pairwise Newtonian gravity with Plummer softening `eps`:
/// `U = -G m_i m_j / sqrt(r^2 + eps^2)`. Direct O(N^2) summation.
#[derive(Debug, Clone)]
pub struct NewtonianGravity {
    pub g: f64,
    pub softening: f64,
}

impl Force for NewtonianGravity {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let eps2 = self.softening * self.softening;
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                let d = pos[j] - pos[i];
                let r2 = d.norm_squared() + eps2;
                let inv_r3 = 1.0 / (r2 * r2.sqrt());
                let s = self.g * inv_r3;
                acc[i] += d * (s * mass[j]);
                acc[j] -= d * (s * mass[i]);
            }
        }
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let eps2 = self.softening * self.softening;
        let mut u = 0.0;
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                let r = ((pos[j] - pos[i]).norm_squared() + eps2).sqrt();
                u -= self.g * mass[i] * mass[j] / r;
            }
        }
        Ok(Some(u))
    }

    fn name(&self) -> String {
        "NewtonianGravity".into()
    }
}

/// Hookean spring between particles `i` and `j`: `U = k/2 (|r_j - r_i| - L)^2`.
#[derive(Debug, Clone)]
pub struct Spring {
    pub i: usize,
    pub j: usize,
    pub k: f64,
    pub rest_length: f64,
}

impl Force for Spring {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        check_index(self.i, pos.len(), "Spring")?;
        check_index(self.j, pos.len(), "Spring")?;
        let d = pos[self.j] - pos[self.i];
        let r = d.norm();
        if r == 0.0 {
            // Direction undefined; only possible to resolve when the spring has zero rest length.
            return Ok(());
        }
        let f = d * (self.k * (r - self.rest_length) / r); // force on i
        acc[self.i] += f / mass[self.i];
        acc[self.j] -= f / mass[self.j];
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], _mass: &[f64]) -> Result<Option<f64>> {
        check_index(self.i, pos.len(), "Spring")?;
        check_index(self.j, pos.len(), "Spring")?;
        let stretch = (pos[self.j] - pos[self.i]).norm() - self.rest_length;
        Ok(Some(0.5 * self.k * stretch * stretch))
    }

    fn name(&self) -> String {
        format!("Spring({}, {})", self.i, self.j)
    }
}

/// Spring tying particle `i` to a fixed point: `U = k/2 (|r_i - anchor| - L)^2`.
/// With `L = 0` this is an isotropic harmonic trap; a stiff spring with `L > 0`
/// approximates a pendulum rod.
#[derive(Debug, Clone)]
pub struct AnchorSpring {
    pub i: usize,
    pub anchor: Vec3,
    pub k: f64,
    pub rest_length: f64,
}

impl Force for AnchorSpring {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        check_index(self.i, pos.len(), "AnchorSpring")?;
        let d = self.anchor - pos[self.i];
        let f = if self.rest_length == 0.0 {
            d * self.k
        } else {
            let r = d.norm();
            if r == 0.0 {
                return Ok(());
            }
            d * (self.k * (r - self.rest_length) / r)
        };
        acc[self.i] += f / mass[self.i];
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], _mass: &[f64]) -> Result<Option<f64>> {
        check_index(self.i, pos.len(), "AnchorSpring")?;
        let stretch = (pos[self.i] - self.anchor).norm() - self.rest_length;
        Ok(Some(0.5 * self.k * stretch * stretch))
    }

    fn name(&self) -> String {
        format!("AnchorSpring({})", self.i)
    }
}

/// Linear (Stokes) drag `a = -gamma v`; `gamma` is a rate (1/time), independent of mass.
#[derive(Debug, Clone)]
pub struct LinearDrag {
    pub gamma: f64,
}

impl Force for LinearDrag {
    fn accumulate(
        &self,
        _t: f64,
        _pos: &[Vec3],
        vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        for (a, v) in acc.iter_mut().zip(vel) {
            *a -= *v * self.gamma;
        }
        Ok(())
    }

    fn velocity_dependent(&self) -> bool {
        true
    }

    fn name(&self) -> String {
        "LinearDrag".into()
    }
}

/// Quadratic (Newtonian) drag force `F = -c |v| v`, so `a = -(c/m) |v| v`.
#[derive(Debug, Clone)]
pub struct QuadraticDrag {
    pub c: f64,
}

impl Force for QuadraticDrag {
    fn accumulate(
        &self,
        _t: f64,
        _pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        for ((a, v), m) in acc.iter_mut().zip(vel).zip(mass) {
            *a -= *v * (self.c * v.norm() / m);
        }
        Ok(())
    }

    fn velocity_dependent(&self) -> bool {
        true
    }

    fn name(&self) -> String {
        "QuadraticDrag".into()
    }
}
