//! Electromagnetic forces on charged particles: uniform and user-defined electric and
//! magnetic fields (the Lorentz force `F = q (E + v × B)`), and pairwise Coulomb repulsion.
//! Charges are set per particle with `World::set_charges`/`set_charge`.

use super::{non_negative, scalar, vector, BuiltinForce, Force, InverseSquare, Param};
use crate::error::{invalid, Result};
use crate::vec3::Vec3;

fn needs_charges(name: &str) -> Result<()> {
    invalid(format!(
        "{name} acts on charges; evaluate it through a ForceSet"
    ))
}

/// `q_i / m_i` for every particle; fails for a charged massless particle.
fn charge_to_mass(name: &str, mass: &[f64], charge: &[f64]) -> Result<Vec<f64>> {
    mass.iter()
        .zip(charge)
        .enumerate()
        .map(|(i, (&m, &q))| {
            if q == 0.0 {
                Ok(0.0)
            } else if m == 0.0 {
                invalid(format!("{name}: particle {i} is charged but massless"))
            } else {
                Ok(q / m)
            }
        })
        .collect()
}

/// Uniform electric field `E`: `a = (q/m) E`, potential `-Σ q E·r`.
#[derive(Debug, Clone)]
pub struct ElectricField {
    pub e: Vec3,
}

impl Force for ElectricField {
    fn accumulate(
        &self,
        _t: f64,
        _p: &[Vec3],
        _v: &[Vec3],
        _m: &[f64],
        _a: &mut [Vec3],
    ) -> Result<()> {
        needs_charges("ElectricField")
    }

    fn accumulate_charged(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        charge: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let qm = charge_to_mass("ElectricField", mass, charge)?;
        for (a, k) in acc.iter_mut().zip(qm) {
            *a += self.e * k;
        }
        Ok(())
    }

    fn potential_charged(
        &self,
        _t: f64,
        pos: &[Vec3],
        _mass: &[f64],
        charge: &[f64],
    ) -> Result<Option<f64>> {
        Ok(Some(
            pos.iter()
                .zip(charge)
                .map(|(r, q)| -q * self.e.dot(*r))
                .sum(),
        ))
    }

    fn jacobian_vector_charged(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        _charge: &[f64],
        _dpos: &[Vec3],
        _dvel: &[Vec3],
        _out: &mut [Vec3],
    ) -> Result<bool> {
        Ok(true) // independent of the state
    }

    fn name(&self) -> String {
        "ElectricField".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::ElectricField(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![("E", Param::Vector(self.e))]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "E" => self.e = vector(name, value)?,
            _ => return invalid(format!("ElectricField has no parameter {name:?}")),
        }
        Ok(())
    }
}

/// Uniform magnetic field `B`: `a = (q/m) v × B`. Does no work, so it has no potential.
/// The `boris` integrator rotates velocities about it exactly.
#[derive(Debug, Clone)]
pub struct MagneticField {
    pub b: Vec3,
}

impl Force for MagneticField {
    fn accumulate(
        &self,
        _t: f64,
        _p: &[Vec3],
        _v: &[Vec3],
        _m: &[f64],
        _a: &mut [Vec3],
    ) -> Result<()> {
        needs_charges("MagneticField")
    }

    fn accumulate_charged(
        &self,
        _t: f64,
        _pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        charge: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let qm = charge_to_mass("MagneticField", mass, charge)?;
        for ((a, v), k) in acc.iter_mut().zip(vel).zip(qm) {
            *a += v.cross(self.b) * k;
        }
        Ok(())
    }

    fn jacobian_vector_charged(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        charge: &[f64],
        _dpos: &[Vec3],
        dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let qm = charge_to_mass("MagneticField", mass, charge)?;
        for ((o, dv), k) in out.iter_mut().zip(dvel).zip(qm) {
            *o += dv.cross(self.b) * k;
        }
        Ok(true)
    }

    fn magnetic_field(&self, _t: f64, _pos: &[Vec3], out: &mut [Vec3]) -> Result<bool> {
        for b in out.iter_mut() {
            *b += self.b;
        }
        Ok(true)
    }

    fn accumulate_electric(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        _charge: &[f64],
        _acc: &mut [Vec3],
    ) -> Result<()> {
        Ok(())
    }

    fn velocity_dependent(&self) -> bool {
        true
    }

    fn name(&self) -> String {
        "MagneticField".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::MagneticField(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![("B", Param::Vector(self.b))]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "B" => self.b = vector(name, value)?,
            _ => return invalid(format!("MagneticField has no parameter {name:?}")),
        }
        Ok(())
    }
}

/// Pairwise Coulomb interaction `U = k Σ q_i q_j / √(r² + ε²)` (like charges repel). Shares
/// its pair loop, parallel blocking and softening with [`super::NewtonianGravity`].
#[derive(Debug, Clone)]
pub struct Coulomb {
    /// Coulomb constant `k = 1/(4π ε0)` in the simulation's units.
    pub k: f64,
    pub softening: f64,
}

impl Coulomb {
    fn kernel(&self) -> InverseSquare {
        InverseSquare {
            coef: -self.k,
            softening: self.softening,
        }
    }
}

impl Force for Coulomb {
    fn accumulate(
        &self,
        _t: f64,
        _p: &[Vec3],
        _v: &[Vec3],
        _m: &[f64],
        _a: &mut [Vec3],
    ) -> Result<()> {
        needs_charges("Coulomb")
    }

    fn accumulate_charged(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        charge: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let qm = charge_to_mass("Coulomb", mass, charge)?;
        self.kernel().accumulate(pos, charge, |i| qm[i], acc);
        Ok(())
    }

    fn potential_charged(
        &self,
        _t: f64,
        pos: &[Vec3],
        _mass: &[f64],
        charge: &[f64],
    ) -> Result<Option<f64>> {
        Ok(Some(self.kernel().potential(pos, charge)))
    }

    fn jacobian_vector_charged(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        charge: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        let qm = charge_to_mass("Coulomb", mass, charge)?;
        self.kernel()
            .jacobian_vector(pos, charge, |i| qm[i], dpos, out);
        Ok(true)
    }

    fn name(&self) -> String {
        "Coulomb".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::Coulomb(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("k", Param::Scalar(self.k)),
            ("softening", Param::Scalar(self.softening)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "k" => self.k = scalar(name, value)?,
            "softening" => self.softening = non_negative(name, value)?,
            _ => return invalid(format!("Coulomb has no parameter {name:?}")),
        }
        Ok(())
    }
}

type FieldFn = dyn Fn(f64, &[Vec3], &mut [Vec3]) -> Result<()> + Send + Sync;

/// Electric and magnetic fields given by functions of time and position:
/// `e(t, pos, out)` and `b(t, pos, out)` add the field at every particle position to `out`.
/// The Lorentz force `a = (q/m)(E + v × B)` follows; `boris` uses `B` directly.
pub struct FieldFunctions {
    name: String,
    e: Option<Box<FieldFn>>,
    b: Option<Box<FieldFn>>,
}

impl FieldFunctions {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            e: None,
            b: None,
        }
    }

    pub fn electric(
        mut self,
        e: impl Fn(f64, &[Vec3], &mut [Vec3]) -> Result<()> + Send + Sync + 'static,
    ) -> Self {
        self.e = Some(Box::new(e));
        self
    }

    pub fn magnetic(
        mut self,
        b: impl Fn(f64, &[Vec3], &mut [Vec3]) -> Result<()> + Send + Sync + 'static,
    ) -> Self {
        self.b = Some(Box::new(b));
        self
    }

    fn field(f: &Option<Box<FieldFn>>, t: f64, pos: &[Vec3]) -> Result<Option<Vec<Vec3>>> {
        match f {
            None => Ok(None),
            Some(f) => {
                let mut out = vec![Vec3::ZERO; pos.len()];
                f(t, pos, &mut out)?;
                Ok(Some(out))
            }
        }
    }
}

impl Force for FieldFunctions {
    fn accumulate(
        &self,
        _t: f64,
        _p: &[Vec3],
        _v: &[Vec3],
        _m: &[f64],
        _a: &mut [Vec3],
    ) -> Result<()> {
        needs_charges(&self.name)
    }

    fn accumulate_charged(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        charge: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        self.accumulate_electric(t, pos, vel, mass, charge, acc)?;
        if let Some(b) = Self::field(&self.b, t, pos)? {
            let qm = charge_to_mass(&self.name, mass, charge)?;
            for (i, a) in acc.iter_mut().enumerate() {
                *a += vel[i].cross(b[i]) * qm[i];
            }
        }
        Ok(())
    }

    fn accumulate_electric(
        &self,
        t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        charge: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        if let Some(e) = Self::field(&self.e, t, pos)? {
            let qm = charge_to_mass(&self.name, mass, charge)?;
            for (i, a) in acc.iter_mut().enumerate() {
                *a += e[i] * qm[i];
            }
        }
        Ok(())
    }

    fn magnetic_field(&self, t: f64, pos: &[Vec3], out: &mut [Vec3]) -> Result<bool> {
        if let Some(b) = Self::field(&self.b, t, pos)? {
            for (o, b) in out.iter_mut().zip(b) {
                *o += b;
            }
        }
        Ok(true)
    }

    fn velocity_dependent(&self) -> bool {
        self.b.is_some()
    }

    fn name(&self) -> String {
        self.name.clone()
    }
}
