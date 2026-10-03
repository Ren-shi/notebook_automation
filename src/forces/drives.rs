//! Explicitly time-dependent forcing.

use super::{check_index, check_massive, scalar, shift, vector, BuiltinForce, Force, Param};
use crate::error::{invalid, Result};
use crate::vec3::Vec3;

/// Sinusoidal force on particle `i`: `F(t) = amplitude cos(ω t + φ)` (a driven oscillator).
/// Non-conservative: it does work on the system, so no potential is reported.
#[derive(Debug, Clone)]
pub struct PeriodicForce {
    pub i: usize,
    /// Force amplitude (a 3-vector, so the drive has a direction).
    pub amplitude: Vec3,
    pub omega: f64,
    pub phase: f64,
}

impl Force for PeriodicForce {
    fn jacobian_vector(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        _dpos: &[Vec3],
        _dvel: &[Vec3],
        _out: &mut [Vec3],
    ) -> Result<bool> {
        Ok(true) // independent of the state
    }

    fn accumulate(
        &self,
        t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        check_index(self.i, pos.len(), "PeriodicForce")?;
        check_massive(self.i, mass, "PeriodicForce")?;
        acc[self.i] += self.amplitude * ((self.omega * t + self.phase).cos() / mass[self.i]);
        Ok(())
    }

    fn name(&self) -> String {
        format!("PeriodicForce({})", self.i)
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::PeriodicForce(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("amplitude", Param::Vector(self.amplitude)),
            ("omega", Param::Scalar(self.omega)),
            ("phase", Param::Scalar(self.phase)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "amplitude" => self.amplitude = vector(name, value)?,
            "omega" => self.omega = scalar(name, value)?,
            "phase" => self.phase = scalar(name, value)?,
            _ => return invalid(format!("PeriodicForce has no parameter {name:?}")),
        }
        Ok(())
    }

    fn references_particle(&self, i: usize) -> bool {
        self.i == i
    }

    fn particle_removed(&mut self, removed: usize) {
        shift(&mut self.i, removed);
    }
}
