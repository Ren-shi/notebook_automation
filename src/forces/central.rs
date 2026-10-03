//! External central potentials `Φ(r)` about a fixed point, acting on every particle.
//! Each gives an acceleration `-∇Φ` independent of the particle's mass (so `Φ` is a
//! potential per unit mass) and a potential energy `Σ m_i Φ(r_i)`.

use super::{non_negative, scalar, vector, BuiltinForce, Force, Param};
use crate::error::{invalid, Result};
use crate::parallel;
use crate::vec3::Vec3;

/// Adds `-Φ'(r) r̂` to every acceleration; `radial(r)` returns `Φ'(r) / r`.
/// A particle exactly at the centre gets no acceleration (the direction is undefined).
fn accumulate_central(
    center: Vec3,
    pos: &[Vec3],
    acc: &mut [Vec3],
    radial: impl Fn(f64) -> f64 + Sync,
) {
    parallel::for_each_indexed(acc, |i, a| {
        let d = pos[i] - center;
        let r = d.norm();
        if r > 0.0 {
            *a -= d * radial(r);
        }
    });
}

fn central_potential(center: Vec3, pos: &[Vec3], mass: &[f64], phi: impl Fn(f64) -> f64) -> f64 {
    pos.iter()
        .zip(mass)
        .filter(|(_, &m)| m != 0.0)
        .map(|(r, m)| m * phi((*r - center).norm()))
        .sum()
}

fn positive(name: &str, value: Param) -> Result<f64> {
    let x = scalar(name, value)?;
    if x <= 0.0 {
        return invalid(format!("{name} must be positive, got {x}"));
    }
    Ok(x)
}

/// Power-law potential `Φ = k r^n` about `center` (`Φ = k ln r` for `n = 0`).
/// `n = 2` is a harmonic trap with `ω² = 2k`; `n = -1, k = -GM` is a fixed Kepler centre;
/// `n = 0` gives a flat rotation curve with circular speed `√k`.
#[derive(Debug, Clone)]
pub struct PowerLaw {
    pub center: Vec3,
    pub k: f64,
    pub n: f64,
}

impl Force for PowerLaw {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let (k, n) = (self.k, self.n);
        if n == 0.0 {
            accumulate_central(self.center, pos, acc, |r| k / (r * r));
        } else {
            accumulate_central(self.center, pos, acc, |r| k * n * r.powf(n - 2.0));
        }
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let (k, n) = (self.k, self.n);
        Ok(Some(central_potential(self.center, pos, mass, |r| {
            if n == 0.0 {
                k * r.ln()
            } else {
                k * r.powf(n)
            }
        })))
    }

    fn name(&self) -> String {
        "PowerLaw".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::PowerLaw(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("center", Param::Vector(self.center)),
            ("k", Param::Scalar(self.k)),
            ("n", Param::Scalar(self.n)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "center" => self.center = vector(name, value)?,
            "k" => self.k = scalar(name, value)?,
            "n" => self.n = scalar(name, value)?,
            _ => return invalid(format!("PowerLaw has no parameter {name:?}")),
        }
        Ok(())
    }
}

/// Yukawa (screened) potential `Φ = -k e^{-r/λ} / r` about `center`. Reduces to a Kepler
/// centre with `GM = k` for `λ → ∞`.
#[derive(Debug, Clone)]
pub struct Yukawa {
    pub center: Vec3,
    pub k: f64,
    /// Screening length λ (positive).
    pub length: f64,
}

impl Force for Yukawa {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        positive("length", Param::Scalar(self.length))?;
        let (k, l) = (self.k, self.length);
        accumulate_central(self.center, pos, acc, |r| {
            k * (-r / l).exp() * (1.0 / r + 1.0 / l) / (r * r)
        });
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let (k, l) = (self.k, self.length);
        Ok(Some(central_potential(self.center, pos, mass, |r| {
            -k * (-r / l).exp() / r
        })))
    }

    fn name(&self) -> String {
        "Yukawa".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::Yukawa(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("center", Param::Vector(self.center)),
            ("k", Param::Scalar(self.k)),
            ("length", Param::Scalar(self.length)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "center" => self.center = vector(name, value)?,
            "k" => self.k = scalar(name, value)?,
            "length" => self.length = positive(name, value)?,
            _ => return invalid(format!("Yukawa has no parameter {name:?}")),
        }
        Ok(())
    }
}

/// Plummer sphere potential `Φ = -GM / √(r² + a²)` about `center` (a cored halo or cluster).
#[derive(Debug, Clone)]
pub struct PlummerPotential {
    pub center: Vec3,
    /// G times the total mass.
    pub gm: f64,
    /// Scale radius `a` (non-negative; 0 is a point mass).
    pub a: f64,
}

impl Force for PlummerPotential {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let (gm, a2) = (self.gm, self.a * self.a);
        accumulate_central(self.center, pos, acc, |r| {
            let s2 = r * r + a2;
            gm / (s2 * s2.sqrt())
        });
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let (gm, a2) = (self.gm, self.a * self.a);
        Ok(Some(central_potential(self.center, pos, mass, |r| {
            -gm / (r * r + a2).sqrt()
        })))
    }

    fn name(&self) -> String {
        "PlummerPotential".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::PlummerPotential(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("center", Param::Vector(self.center)),
            ("GM", Param::Scalar(self.gm)),
            ("a", Param::Scalar(self.a)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "center" => self.center = vector(name, value)?,
            "GM" => self.gm = scalar(name, value)?,
            "a" => self.a = non_negative(name, value)?,
            _ => return invalid(format!("PlummerPotential has no parameter {name:?}")),
        }
        Ok(())
    }
}

/// Hernquist potential `Φ = -GM / (r + a)` about `center` (galactic bulges and halos).
#[derive(Debug, Clone)]
pub struct HernquistPotential {
    pub center: Vec3,
    /// G times the total mass.
    pub gm: f64,
    /// Scale radius `a` (non-negative; 0 is a point mass).
    pub a: f64,
}

impl Force for HernquistPotential {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let (gm, a) = (self.gm, self.a);
        accumulate_central(self.center, pos, acc, |r| gm / (r * (r + a) * (r + a)));
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let (gm, a) = (self.gm, self.a);
        Ok(Some(central_potential(self.center, pos, mass, |r| {
            -gm / (r + a)
        })))
    }

    fn name(&self) -> String {
        "HernquistPotential".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::HernquistPotential(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("center", Param::Vector(self.center)),
            ("GM", Param::Scalar(self.gm)),
            ("a", Param::Scalar(self.a)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "center" => self.center = vector(name, value)?,
            "GM" => self.gm = scalar(name, value)?,
            "a" => self.a = non_negative(name, value)?,
            _ => return invalid(format!("HernquistPotential has no parameter {name:?}")),
        }
        Ok(())
    }
}

/// Harmonic trap `Φ = ½ Σ_k ω_k² (x_k - c_k)²` with per-axis angular frequencies `omega`
/// (isotropic when all three are equal). Acts on every particle regardless of mass.
#[derive(Debug, Clone)]
pub struct HarmonicTrap {
    pub center: Vec3,
    pub omega: Vec3,
}

impl Force for HarmonicTrap {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let w = self.omega;
        let w2 = Vec3::new(w.x * w.x, w.y * w.y, w.z * w.z);
        parallel::for_each_indexed(acc, |i, a| {
            let d = pos[i] - self.center;
            *a -= Vec3::new(w2.x * d.x, w2.y * d.y, w2.z * d.z);
        });
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let w = self.omega;
        Ok(Some(
            pos.iter()
                .zip(mass)
                .map(|(r, m)| {
                    let d = *r - self.center;
                    0.5 * m
                        * (w.x * w.x * d.x * d.x + w.y * w.y * d.y * d.y + w.z * w.z * d.z * d.z)
                })
                .sum(),
        ))
    }

    fn name(&self) -> String {
        "HarmonicTrap".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::HarmonicTrap(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("center", Param::Vector(self.center)),
            ("omega", Param::Vector(self.omega)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "center" => self.center = vector(name, value)?,
            "omega" => self.omega = vector(name, value)?,
            _ => return invalid(format!("HarmonicTrap has no parameter {name:?}")),
        }
        Ok(())
    }
}
