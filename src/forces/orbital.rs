//! Perturbations of orbits about a central body: general relativity at first
//! post-Newtonian order, and the oblateness (J2) of the central body.

use super::{check_index, scalar, shift, vector, BuiltinForce, Force, Param};
use crate::error::{invalid, Result};
use crate::parallel;
use crate::vec3::Vec3;

fn positive(name: &str, value: Param) -> Result<f64> {
    let x = scalar(name, value)?;
    if x <= 0.0 {
        return invalid(format!("{name} must be positive, got {x}"));
    }
    Ok(x)
}

/// First post-Newtonian (1PN) correction from the central body `central`, in the
/// test-particle limit (harmonic coordinates). Every other particle, at position `r` and
/// velocity `v` relative to the central body, gets
/// `a = GM / (c² r³) [(4 GM / r - v²) r + 4 (r·v) v]` with `M` the central body's mass.
/// Add it to Newtonian gravity: it produces the relativistic periapsis advance
/// `6π GM / (c² a (1 - e²))` per orbit. No back-reaction on the central body.
#[derive(Debug, Clone)]
pub struct PostNewtonian {
    pub central: usize,
    pub g: f64,
    /// Speed of light, in the simulation's units.
    pub c: f64,
}

impl Force for PostNewtonian {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        check_index(self.central, pos.len(), "PostNewtonian")?;
        positive("c", Param::Scalar(self.c))?;
        let gm = self.g * mass[self.central];
        let (rc, vc) = (pos[self.central], vel[self.central]);
        let c2 = self.c * self.c;
        parallel::for_each_indexed(acc, |p, a| {
            if p == self.central {
                return;
            }
            let r = pos[p] - rc;
            let v = vel[p] - vc;
            let d = r.norm();
            if d == 0.0 {
                return;
            }
            let s = gm / (c2 * d * d * d);
            *a += (r * (4.0 * gm / d - v.norm_squared()) + v * (4.0 * r.dot(v))) * s;
        });
        Ok(())
    }

    fn velocity_dependent(&self) -> bool {
        true
    }

    fn name(&self) -> String {
        format!("PostNewtonian({})", self.central)
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::PostNewtonian(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![("G", Param::Scalar(self.g)), ("c", Param::Scalar(self.c))]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "G" => self.g = scalar(name, value)?,
            "c" => self.c = positive(name, value)?,
            _ => return invalid(format!("PostNewtonian has no parameter {name:?}")),
        }
        Ok(())
    }

    fn references_particle(&self, i: usize) -> bool {
        self.central == i
    }

    fn particle_removed(&mut self, removed: usize) {
        shift(&mut self.central, removed);
    }
}

/// Oblateness of the central body `central`: the J2 term of its gravity field,
/// `Φ = GM J2 R² (3 z²/r² - 1) / (2 r³)` per unit mass, with `z` measured along the unit
/// `axis` (the body's spin axis) and `r` from the body's centre. Add it to Newtonian
/// gravity. The central body feels the equal and opposite reaction, so momentum is conserved.
#[derive(Debug, Clone)]
pub struct J2Oblateness {
    pub central: usize,
    pub g: f64,
    pub j2: f64,
    /// Equatorial radius `R` of the central body.
    pub radius: f64,
    /// Spin axis; normalised when used.
    pub axis: Vec3,
}

impl J2Oblateness {
    fn unit_axis(&self) -> Result<Vec3> {
        let n = self.axis.norm();
        if !(n > 0.0 && n.is_finite()) {
            return invalid("J2Oblateness: axis must be a non-zero, finite vector");
        }
        Ok(self.axis / n)
    }
}

impl Force for J2Oblateness {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        check_index(self.central, pos.len(), "J2Oblateness")?;
        let mc = mass[self.central];
        if mc == 0.0 {
            return invalid(format!(
                "J2Oblateness: central body {} is massless",
                self.central
            ));
        }
        let e = self.unit_axis()?;
        let k = 1.5 * self.g * mc * self.j2 * self.radius * self.radius;
        let rc = pos[self.central];
        let mut reaction = Vec3::ZERO;
        for (p, (a, &m)) in acc.iter_mut().zip(mass).enumerate() {
            if p == self.central {
                continue;
            }
            let d = pos[p] - rc;
            let r2 = d.norm_squared();
            if r2 == 0.0 {
                continue;
            }
            let z = d.dot(e);
            let r5 = r2 * r2 * r2.sqrt();
            let ap = (e * (2.0 * z) + d * (1.0 - 5.0 * z * z / r2)) * (-k / r5);
            *a += ap;
            reaction -= ap * m;
        }
        acc[self.central] += reaction / mc;
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        check_index(self.central, pos.len(), "J2Oblateness")?;
        let e = self.unit_axis()?;
        let k = 0.5 * self.g * mass[self.central] * self.j2 * self.radius * self.radius;
        let rc = pos[self.central];
        let mut u = 0.0;
        for (p, &m) in mass.iter().enumerate() {
            if p == self.central || m == 0.0 {
                continue;
            }
            let d = pos[p] - rc;
            let r2 = d.norm_squared();
            let z = d.dot(e);
            u += m * k * (3.0 * z * z / r2 - 1.0) / (r2 * r2.sqrt());
        }
        Ok(Some(u))
    }

    fn name(&self) -> String {
        format!("J2Oblateness({})", self.central)
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::J2Oblateness(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("G", Param::Scalar(self.g)),
            ("J2", Param::Scalar(self.j2)),
            ("radius", Param::Scalar(self.radius)),
            ("axis", Param::Vector(self.axis)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "G" => self.g = scalar(name, value)?,
            "J2" => self.j2 = scalar(name, value)?,
            "radius" => self.radius = positive(name, value)?,
            "axis" => {
                let axis = vector(name, value)?;
                let old = std::mem::replace(&mut self.axis, axis);
                if let Err(e) = self.unit_axis() {
                    self.axis = old;
                    return Err(e);
                }
            }
            _ => return invalid(format!("J2Oblateness has no parameter {name:?}")),
        }
        Ok(())
    }

    fn references_particle(&self, i: usize) -> bool {
        self.central == i
    }

    fn particle_removed(&mut self, removed: usize) {
        shift(&mut self.central, removed);
    }
}
