//! Soft contact between particles with a radius: overlapping spheres repel.

use super::{non_negative, BuiltinForce, Force, Param};
use crate::broadphase::candidate_pairs;
use crate::error::{invalid, Result};
use crate::vec3::Vec3;

/// Elastic law of [`SoftContact`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactLaw {
    /// `F = k δ` (a linear spring while overlapping).
    Linear,
    /// `F = k δ^{3/2}` (Hertz contact of elastic spheres).
    Hertz,
}

/// Soft contact: particles whose spheres overlap by `δ = r_i + r_j - |x_j - x_i| > 0` are
/// pushed apart along the line of centres by an elastic force (`ContactLaw`) plus a dashpot
/// `-c u_n` on the normal approach speed `u_n` (the total is never attractive). Elastic
/// potential `k δ²/2` (linear) or `2 k δ^{5/2}/5` (Hertz). Uses each particle's radius; a
/// uniform grid finds the overlapping pairs in O(N). Contacts must be resolved in time:
/// use steps of a small fraction of the contact duration (`π √(μ/k)` for the linear law).
#[derive(Debug, Clone)]
pub struct SoftContact {
    pub k: f64,
    /// Dashpot coefficient (force per unit normal approach speed).
    pub damping: f64,
    pub law: ContactLaw,
}

impl SoftContact {
    fn elastic(&self, delta: f64) -> (f64, f64) {
        match self.law {
            ContactLaw::Linear => (self.k * delta, 0.5 * self.k * delta * delta),
            ContactLaw::Hertz => {
                let s = delta.sqrt();
                (self.k * delta * s, 0.4 * self.k * delta * delta * s)
            }
        }
    }
}

impl Force for SoftContact {
    fn accumulate(
        &self,
        _t: f64,
        _p: &[Vec3],
        _v: &[Vec3],
        _m: &[f64],
        _a: &mut [Vec3],
    ) -> Result<()> {
        invalid("SoftContact needs the particles' radii; evaluate it through a ForceSet")
    }

    fn accumulate_full(
        &self,
        _t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        _charge: &[f64],
        radius: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        for (i, j) in candidate_pairs(pos, radius) {
            let d = pos[j] - pos[i];
            let dist = d.norm();
            let delta = radius[i] + radius[j] - dist;
            if delta <= 0.0 || dist == 0.0 {
                continue;
            }
            if mass[i] == 0.0 || mass[j] == 0.0 {
                return invalid(format!(
                    "SoftContact: particles {i} and {j} touch, but one is massless"
                ));
            }
            let n = d / dist;
            let approach = (vel[j] - vel[i]).dot(n);
            let force = (self.elastic(delta).0 - self.damping * approach).max(0.0);
            acc[j] += n * (force / mass[j]);
            acc[i] -= n * (force / mass[i]);
        }
        Ok(())
    }

    fn potential_full(
        &self,
        _t: f64,
        pos: &[Vec3],
        _mass: &[f64],
        _charge: &[f64],
        radius: &[f64],
    ) -> Result<Option<f64>> {
        let mut u = 0.0;
        for (i, j) in candidate_pairs(pos, radius) {
            let delta = radius[i] + radius[j] - (pos[j] - pos[i]).norm();
            if delta > 0.0 {
                u += self.elastic(delta).1;
            }
        }
        Ok(Some(u))
    }

    fn velocity_dependent(&self) -> bool {
        self.damping != 0.0
    }

    fn name(&self) -> String {
        "SoftContact".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::SoftContact(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("k", Param::Scalar(self.k)),
            ("damping", Param::Scalar(self.damping)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "k" => self.k = non_negative(name, value)?,
            "damping" => self.damping = non_negative(name, value)?,
            _ => return invalid(format!("SoftContact has no parameter {name:?}")),
        }
        Ok(())
    }
}
