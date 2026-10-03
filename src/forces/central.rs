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

/// Adds the directional derivative of `-g(r) d` (with `d = r - center`) along `dpos`:
/// `-g δd - g'(r) (d·δd / r) d`; `g_and_slope(r)` returns `(g, dg/dr)`.
fn jvp_central(
    center: Vec3,
    pos: &[Vec3],
    dpos: &[Vec3],
    out: &mut [Vec3],
    g_and_slope: impl Fn(f64) -> (f64, f64),
) {
    for ((x, dx), o) in pos.iter().zip(dpos).zip(out.iter_mut()) {
        let d = *x - center;
        let r = d.norm();
        if r > 0.0 {
            let (g, dg) = g_and_slope(r);
            *o -= *dx * g + d * (dg * d.dot(*dx) / r);
        }
    }
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
    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let (k, n) = (self.k, self.n);
        jvp_central(self.center, pos, dpos, out, |r| {
            if n == 0.0 {
                (k / (r * r), -2.0 * k / (r * r * r))
            } else {
                (k * n * r.powf(n - 2.0), k * n * (n - 2.0) * r.powf(n - 3.0))
            }
        });
        Ok(true)
    }

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
    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let (k, l) = (self.k, self.length);
        jvp_central(self.center, pos, dpos, out, |r| {
            let e = (-r / l).exp();
            let inner = 1.0 / (r * r * r) + 1.0 / (l * r * r);
            let d_inner = -3.0 / (r * r * r * r) - 2.0 / (l * r * r * r);
            (k * e * inner, k * e * (d_inner - inner / l))
        });
        Ok(true)
    }

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
    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let (gm, a2) = (self.gm, self.a * self.a);
        jvp_central(self.center, pos, dpos, out, |r| {
            let s2 = r * r + a2;
            let s3 = s2 * s2.sqrt();
            (gm / s3, -3.0 * gm * r / (s3 * s2))
        });
        Ok(true)
    }

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
    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let (gm, a) = (self.gm, self.a);
        jvp_central(self.center, pos, dpos, out, |r| {
            let ra = r + a;
            (
                gm / (r * ra * ra),
                -gm * (3.0 * r + a) / (r * r * ra * ra * ra),
            )
        });
        Ok(true)
    }

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
    fn jacobian_vector(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let w = self.omega;
        for (o, d) in out.iter_mut().zip(dpos) {
            *o -= Vec3::new(w.x * w.x * d.x, w.y * w.y * d.y, w.z * w.z * d.z);
        }
        Ok(true)
    }

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

/// The Hénon-Heiles potential `Φ = ½(x² + y²) + λ (x² y - y³/3)` per unit mass in the
/// xy-plane about `center` (`λ = 1` is the classic system; z is force-free). The standard
/// test bed for chaos: mostly regular below E = 1/12, mostly chaotic near the escape
/// energy E = 1/6 (for λ = 1).
#[derive(Debug, Clone)]
pub struct HenonHeiles {
    pub center: Vec3,
    pub lambda: f64,
}

impl Force for HenonHeiles {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let l = self.lambda;
        parallel::for_each_indexed(acc, |i, a| {
            let d = pos[i] - self.center;
            a.x -= d.x + 2.0 * l * d.x * d.y;
            a.y -= d.y + l * (d.x * d.x - d.y * d.y);
        });
        Ok(())
    }

    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let l = self.lambda;
        for ((x, dx), o) in pos.iter().zip(dpos).zip(out.iter_mut()) {
            let d = *x - self.center;
            o.x -= dx.x + 2.0 * l * (d.y * dx.x + d.x * dx.y);
            o.y -= dx.y + 2.0 * l * (d.x * dx.x - d.y * dx.y);
        }
        Ok(true)
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let l = self.lambda;
        Ok(Some(
            pos.iter()
                .zip(mass)
                .map(|(r, m)| {
                    let d = *r - self.center;
                    m * (0.5 * (d.x * d.x + d.y * d.y)
                        + l * (d.x * d.x * d.y - d.y * d.y * d.y / 3.0))
                })
                .sum(),
        ))
    }

    fn name(&self) -> String {
        "HenonHeiles".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::HenonHeiles(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("center", Param::Vector(self.center)),
            ("lam", Param::Scalar(self.lambda)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "center" => self.center = vector(name, value)?,
            "lam" => self.lambda = scalar(name, value)?,
            _ => return invalid(format!("HenonHeiles has no parameter {name:?}")),
        }
        Ok(())
    }
}
