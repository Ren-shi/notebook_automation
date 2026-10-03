//! Springs beyond [`super::Spring`]: with a dashpot, with a time-modulated stiffness, and
//! whole networks of bonds in one force.

use super::{check_index, check_massive, non_negative, scalar, shift, BuiltinForce, Force, Param};
use crate::constraints::Anchor;
use crate::error::{invalid, Result};
use crate::vec3::Vec3;

/// Spring with a dashpot between particles `i` and `j`: the force on `i` along the unit
/// bond vector `n` (from `i` to `j`) is `[k (r - L) + c (v_j - v_i)·n] n`, so the dashpot
/// damps stretching only. Potential `k/2 (r - L)²`; the dashpot dissipates.
#[derive(Debug, Clone)]
pub struct DampedSpring {
    pub i: usize,
    pub j: usize,
    pub k: f64,
    pub rest_length: f64,
    /// Damping coefficient (force per unit stretching speed), non-negative.
    pub c: f64,
}

impl Force for DampedSpring {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        check_index(self.i, pos.len(), "DampedSpring")?;
        check_index(self.j, pos.len(), "DampedSpring")?;
        check_massive(self.i, mass, "DampedSpring")?;
        check_massive(self.j, mass, "DampedSpring")?;
        let d = pos[self.j] - pos[self.i];
        let r = d.norm();
        if r == 0.0 {
            return Ok(());
        }
        let n = d / r;
        let stretch_rate = (vel[self.j] - vel[self.i]).dot(n);
        let f = n * (self.k * (r - self.rest_length) + self.c * stretch_rate); // force on i
        acc[self.i] += f / mass[self.i];
        acc[self.j] -= f / mass[self.j];
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], _mass: &[f64]) -> Result<Option<f64>> {
        check_index(self.i, pos.len(), "DampedSpring")?;
        check_index(self.j, pos.len(), "DampedSpring")?;
        let stretch = (pos[self.j] - pos[self.i]).norm() - self.rest_length;
        Ok(Some(0.5 * self.k * stretch * stretch))
    }

    fn velocity_dependent(&self) -> bool {
        self.c != 0.0
    }

    fn name(&self) -> String {
        format!("DampedSpring({}, {})", self.i, self.j)
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::DampedSpring(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("k", Param::Scalar(self.k)),
            ("rest_length", Param::Scalar(self.rest_length)),
            ("c", Param::Scalar(self.c)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "k" => self.k = scalar(name, value)?,
            "rest_length" => self.rest_length = non_negative(name, value)?,
            "c" => self.c = non_negative(name, value)?,
            _ => return invalid(format!("DampedSpring has no parameter {name:?}")),
        }
        Ok(())
    }

    fn references_particle(&self, i: usize) -> bool {
        self.i == i || self.j == i
    }

    fn particle_removed(&mut self, removed: usize) {
        shift(&mut self.i, removed);
        shift(&mut self.j, removed);
    }
}

/// Spring whose stiffness is modulated in time, `k(t) = k (1 + depth cos(ω t + φ))`,
/// tying particle `i` to another particle or a fixed point: parametric driving (the
/// Mathieu equation). Potential `k(t)/2 (r - L)²` at the current time.
#[derive(Debug, Clone)]
pub struct ModulatedSpring {
    pub i: usize,
    pub to: Anchor,
    pub k: f64,
    pub depth: f64,
    pub omega: f64,
    pub phase: f64,
    pub rest_length: f64,
}

impl ModulatedSpring {
    fn stiffness(&self, t: f64) -> f64 {
        self.k * (1.0 + self.depth * (self.omega * t + self.phase).cos())
    }

    fn other(&self, pos: &[Vec3]) -> Vec3 {
        match self.to {
            Anchor::Particle(j) => pos[j],
            Anchor::Point(p) => p,
        }
    }

    fn check(&self, pos: &[Vec3], mass: &[f64]) -> Result<()> {
        check_index(self.i, pos.len(), "ModulatedSpring")?;
        check_massive(self.i, mass, "ModulatedSpring")?;
        if let Anchor::Particle(j) = self.to {
            check_index(j, pos.len(), "ModulatedSpring")?;
            check_massive(j, mass, "ModulatedSpring")?;
        }
        Ok(())
    }
}

impl Force for ModulatedSpring {
    fn accumulate(
        &self,
        t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        self.check(pos, mass)?;
        let d = self.other(pos) - pos[self.i];
        let k = self.stiffness(t);
        let f = if self.rest_length == 0.0 {
            d * k
        } else {
            let r = d.norm();
            if r == 0.0 {
                return Ok(());
            }
            d * (k * (r - self.rest_length) / r)
        };
        acc[self.i] += f / mass[self.i];
        if let Anchor::Particle(j) = self.to {
            acc[j] -= f / mass[j];
        }
        Ok(())
    }

    fn potential(&self, t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        self.check(pos, mass)?;
        let stretch = (self.other(pos) - pos[self.i]).norm() - self.rest_length;
        Ok(Some(0.5 * self.stiffness(t) * stretch * stretch))
    }

    fn name(&self) -> String {
        match self.to {
            Anchor::Particle(j) => format!("ModulatedSpring({}, {j})", self.i),
            Anchor::Point(p) => format!("ModulatedSpring({}, [{}, {}, {}])", self.i, p.x, p.y, p.z),
        }
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::ModulatedSpring(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("k", Param::Scalar(self.k)),
            ("depth", Param::Scalar(self.depth)),
            ("omega", Param::Scalar(self.omega)),
            ("phase", Param::Scalar(self.phase)),
            ("rest_length", Param::Scalar(self.rest_length)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "k" => self.k = scalar(name, value)?,
            "depth" => self.depth = scalar(name, value)?,
            "omega" => self.omega = scalar(name, value)?,
            "phase" => self.phase = scalar(name, value)?,
            "rest_length" => self.rest_length = non_negative(name, value)?,
            _ => return invalid(format!("ModulatedSpring has no parameter {name:?}")),
        }
        Ok(())
    }

    fn references_particle(&self, i: usize) -> bool {
        self.i == i || self.to == Anchor::Particle(i)
    }

    fn particle_removed(&mut self, removed: usize) {
        shift(&mut self.i, removed);
        if let Anchor::Particle(j) = &mut self.to {
            shift(j, removed);
        }
    }
}

/// Many Hookean springs in one force, stored as arrays: bond `b` joins particles `i[b]`
/// and `j[b]` with stiffness `k[b]` and rest length `rest_length[b]`. Equivalent to one
/// [`super::Spring`] per bond (bit for bit), without per-spring dispatch.
#[derive(Debug, Clone)]
pub struct SpringNetwork {
    i: Vec<usize>,
    j: Vec<usize>,
    k: Vec<f64>,
    rest_length: Vec<f64>,
}

impl SpringNetwork {
    /// Fails if the arrays differ in length, a bond joins a particle to itself, or a
    /// stiffness or rest length is not finite (rest lengths must also be non-negative).
    pub fn new(i: Vec<usize>, j: Vec<usize>, k: Vec<f64>, rest_length: Vec<f64>) -> Result<Self> {
        let n = i.len();
        if j.len() != n || k.len() != n || rest_length.len() != n {
            return invalid(format!(
                "SpringNetwork: i, j, k and rest_length must have equal lengths (got {}, {}, {}, {})",
                n,
                j.len(),
                k.len(),
                rest_length.len()
            ));
        }
        if let Some(b) = (0..n).find(|&b| i[b] == j[b]) {
            return invalid(format!(
                "SpringNetwork: bond {b} joins particle {} to itself",
                i[b]
            ));
        }
        if let Some(b) = (0..n).find(|&b| !k[b].is_finite()) {
            return invalid(format!(
                "SpringNetwork: k[{b}] must be finite, got {}",
                k[b]
            ));
        }
        if let Some(b) = (0..n).find(|&b| !(rest_length[b].is_finite() && rest_length[b] >= 0.0)) {
            return invalid(format!(
                "SpringNetwork: rest_length[{b}] must be finite and non-negative, got {}",
                rest_length[b]
            ));
        }
        Ok(Self {
            i,
            j,
            k,
            rest_length,
        })
    }

    pub fn len(&self) -> usize {
        self.i.len()
    }

    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    pub fn i(&self) -> &[usize] {
        &self.i
    }

    pub fn j(&self) -> &[usize] {
        &self.j
    }

    pub fn k(&self) -> &[f64] {
        &self.k
    }

    pub fn rest_length(&self) -> &[f64] {
        &self.rest_length
    }

    fn check(&self, n: usize, mass: &[f64]) -> Result<()> {
        for (b, (&i, &j)) in self.i.iter().zip(&self.j).enumerate() {
            if i >= n || j >= n {
                return invalid(format!(
                    "SpringNetwork: bond {b} ({i}, {j}) refers to a particle out of range (have {n} particles)"
                ));
            }
            if mass[i] == 0.0 || mass[j] == 0.0 {
                return invalid(format!(
                    "SpringNetwork: bond {b} ({i}, {j}) acts on a massless particle"
                ));
            }
        }
        Ok(())
    }
}

impl Force for SpringNetwork {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let n = pos.len();
        for b in 0..self.i.len() {
            let (i, j) = (self.i[b], self.j[b]);
            // Checked inline: a separate validation pass costs ~20% of the evaluation.
            if i >= n || j >= n || mass[i] == 0.0 || mass[j] == 0.0 {
                return self.check(n, mass);
            }
            let d = pos[j] - pos[i];
            let r = d.norm();
            if r == 0.0 {
                continue;
            }
            let f = d * (self.k[b] * (r - self.rest_length[b]) / r);
            acc[i] += f / mass[i];
            acc[j] -= f / mass[j];
        }
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        self.check(pos.len(), mass)?;
        let mut u = 0.0;
        for b in 0..self.i.len() {
            let stretch = (pos[self.j[b]] - pos[self.i[b]]).norm() - self.rest_length[b];
            u += 0.5 * self.k[b] * stretch * stretch;
        }
        Ok(Some(u))
    }

    fn name(&self) -> String {
        format!("SpringNetwork({} bonds)", self.i.len())
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::SpringNetwork(self.clone()))
    }

    fn references_particle(&self, p: usize) -> bool {
        self.i.contains(&p) || self.j.contains(&p)
    }

    fn particle_removed(&mut self, removed: usize) {
        for x in self.i.iter_mut().chain(self.j.iter_mut()) {
            shift(x, removed);
        }
    }
}
