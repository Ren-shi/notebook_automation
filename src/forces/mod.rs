//! Forces act on the whole system at once: each one adds its contribution to the
//! per-particle accelerations. Implement [`Force`] to add new physics in Rust.

use std::sync::atomic::{AtomicU64, Ordering};

use crate::error::{invalid, Result};
use crate::parallel;
use crate::vec3::Vec3;

mod central;
mod closure;
mod drives;
mod orbital;
mod springs;

pub use central::{
    HarmonicTrap, HenonHeiles, HernquistPotential, PlummerPotential, PowerLaw, Yukawa,
};
pub use closure::ClosureForce;
pub use drives::PeriodicForce;
pub use orbital::{J2Oblateness, PostNewtonian};
pub use springs::{DampedSpring, ModulatedSpring, SpringNetwork};

/// Value of a tunable force parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Param {
    Scalar(f64),
    Vector(Vec3),
}

/// Forces must be pure functions of their arguments (`t`, positions, velocities, masses) and
/// their parameters: integrators may reuse an acceleration computed for identical inputs.
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

    /// Adds the Jacobian-vector product of this force's accelerations to `out`:
    /// `out_i += Σ_j (∂a_i/∂x_j · dpos_j + ∂a_i/∂v_j · dvel_j)`. Used by the variational
    /// equations (Lyapunov exponents). Returns `Ok(false)` if not implemented, in which case
    /// [`ForceSet::jacobian_vector`] uses central finite differences instead.
    #[allow(clippy::too_many_arguments)]
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
        Ok(false)
    }

    /// For the built-in forces, a copy that can be saved and rebuilt (see
    /// [`BuiltinForce`]). Forces defined elsewhere return `None` and must be supplied
    /// again when a checkpoint is loaded.
    fn builtin(&self) -> Option<BuiltinForce> {
        None
    }

    /// Current values of the tunable parameters.
    fn params(&self) -> Vec<(&'static str, Param)> {
        Vec::new()
    }

    /// Changes one parameter. Must leave the force unchanged on error.
    fn set_param(&mut self, name: &str, _value: Param) -> Result<()> {
        invalid(format!("{} has no parameter {name:?}", self.name()))
    }

    /// Whether this force refers to particle `i` by index.
    fn references_particle(&self, _i: usize) -> bool {
        false
    }

    /// Particle `removed` was deleted: shift stored indices above it down by one.
    /// Only called when `references_particle(removed)` is false.
    fn particle_removed(&mut self, _removed: usize) {}
}

/// A built-in force, as plain data: what checkpoints store.
#[derive(Debug, Clone)]
pub enum BuiltinForce {
    UniformField(UniformField),
    NewtonianGravity(NewtonianGravity),
    Spring(Spring),
    AnchorSpring(AnchorSpring),
    LinearDrag(LinearDrag),
    QuadraticDrag(QuadraticDrag),
    DampedSpring(DampedSpring),
    ModulatedSpring(ModulatedSpring),
    SpringNetwork(SpringNetwork),
    PowerLaw(PowerLaw),
    Yukawa(Yukawa),
    PlummerPotential(PlummerPotential),
    HernquistPotential(HernquistPotential),
    HarmonicTrap(HarmonicTrap),
    PeriodicForce(PeriodicForce),
    PostNewtonian(PostNewtonian),
    J2Oblateness(J2Oblateness),
    HenonHeiles(HenonHeiles),
}

impl BuiltinForce {
    pub fn into_force(self) -> Box<dyn Force> {
        match self {
            BuiltinForce::UniformField(f) => Box::new(f),
            BuiltinForce::NewtonianGravity(f) => Box::new(f),
            BuiltinForce::Spring(f) => Box::new(f),
            BuiltinForce::AnchorSpring(f) => Box::new(f),
            BuiltinForce::LinearDrag(f) => Box::new(f),
            BuiltinForce::QuadraticDrag(f) => Box::new(f),
            BuiltinForce::DampedSpring(f) => Box::new(f),
            BuiltinForce::ModulatedSpring(f) => Box::new(f),
            BuiltinForce::SpringNetwork(f) => Box::new(f),
            BuiltinForce::PowerLaw(f) => Box::new(f),
            BuiltinForce::Yukawa(f) => Box::new(f),
            BuiltinForce::PlummerPotential(f) => Box::new(f),
            BuiltinForce::HernquistPotential(f) => Box::new(f),
            BuiltinForce::HarmonicTrap(f) => Box::new(f),
            BuiltinForce::PeriodicForce(f) => Box::new(f),
            BuiltinForce::PostNewtonian(f) => Box::new(f),
            BuiltinForce::J2Oblateness(f) => Box::new(f),
            BuiltinForce::HenonHeiles(f) => Box::new(f),
        }
    }
}

/// Identifies a force within a [`ForceSet`]; stays valid until that force is removed.
pub type ForceId = usize;

/// A value no force set has used before, so a cached acceleration can never be
/// mistaken for one computed by a different set (e.g. after `world.forces` is replaced).
fn fresh_version() -> u64 {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    NEXT.fetch_add(1, Ordering::Relaxed)
}

/// The collection of forces acting on a system.
pub struct ForceSet {
    forces: Vec<(ForceId, Box<dyn Force>)>,
    next_id: ForceId,
    /// Changed on every change to the set or to a force's parameters.
    version: u64,
}

impl Default for ForceSet {
    fn default() -> Self {
        Self {
            forces: Vec::new(),
            next_id: 0,
            version: fresh_version(),
        }
    }
}

impl ForceSet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, force: Box<dyn Force>) -> ForceId {
        let id = self.next_id;
        self.next_id += 1;
        self.version = fresh_version();
        self.forces.push((id, force));
        id
    }

    fn position(&self, id: ForceId) -> Result<usize> {
        match self.forces.iter().position(|(i, _)| *i == id) {
            Some(p) => Ok(p),
            None => invalid(format!("no force with id {id}")),
        }
    }

    pub fn remove(&mut self, id: ForceId) -> Result<Box<dyn Force>> {
        let p = self.position(id)?;
        self.version = fresh_version();
        Ok(self.forces.remove(p).1)
    }

    /// Swaps in a new force under the same id and evaluation order.
    pub fn replace(&mut self, id: ForceId, force: Box<dyn Force>) -> Result<Box<dyn Force>> {
        let p = self.position(id)?;
        self.version = fresh_version();
        Ok(std::mem::replace(&mut self.forces[p].1, force))
    }

    pub fn get(&self, id: ForceId) -> Result<&dyn Force> {
        let p = self.position(id)?;
        Ok(self.forces[p].1.as_ref())
    }

    /// Sets several parameters of one force; on any error, already-applied ones are rolled back.
    pub fn set_params(&mut self, id: ForceId, values: &[(String, Param)]) -> Result<()> {
        let p = self.position(id)?;
        self.version = fresh_version();
        let force = &mut self.forces[p].1;
        let old = force.params();
        for (k, (name, value)) in values.iter().enumerate() {
            if let Err(e) = force.set_param(name, *value) {
                for (prev_name, _) in &values[..k] {
                    if let Some((n, v)) = old.iter().find(|(n, _)| n == prev_name) {
                        let _ = force.set_param(n, *v);
                    }
                }
                return Err(e);
            }
        }
        Ok(())
    }

    /// Rebuilds a set with the given ids, e.g. from a checkpoint. Later [`ForceSet::add`]
    /// calls hand out ids from `next_id`.
    pub fn from_parts(forces: Vec<(ForceId, Box<dyn Force>)>, next_id: ForceId) -> Result<Self> {
        for (k, (id, _)) in forces.iter().enumerate() {
            if *id >= next_id {
                return invalid(format!("force id {id} is not below next_id {next_id}"));
            }
            if forces[..k].iter().any(|(other, _)| other == id) {
                return invalid(format!("duplicate force id {id}"));
            }
        }
        Ok(Self {
            forces,
            next_id,
            version: fresh_version(),
        })
    }

    /// The id the next [`ForceSet::add`] will return.
    pub fn next_id(&self) -> ForceId {
        self.next_id
    }

    pub fn clear(&mut self) {
        self.version = fresh_version();
        self.forces.clear();
    }

    pub fn iter(&self) -> impl Iterator<Item = (ForceId, &dyn Force)> {
        self.forces.iter().map(|(id, f)| (*id, f.as_ref()))
    }

    pub fn len(&self) -> usize {
        self.forces.len()
    }

    pub fn is_empty(&self) -> bool {
        self.forces.is_empty()
    }

    /// Changes whenever forces are added, removed, replaced or re-parameterised.
    pub fn version(&self) -> u64 {
        self.version
    }

    /// Whether any force depends on velocity.
    pub fn velocity_dependent(&self) -> bool {
        self.forces.iter().any(|(_, f)| f.velocity_dependent())
    }

    /// Names of the forces that refer to particle `i` by index.
    pub fn referencing(&self, i: usize) -> Vec<String> {
        self.forces
            .iter()
            .filter(|(_, f)| f.references_particle(i))
            .map(|(id, f)| format!("{} (id {id})", f.name()))
            .collect()
    }

    pub(crate) fn particle_removed(&mut self, removed: usize) {
        self.version = fresh_version();
        for (_, f) in &mut self.forces {
            f.particle_removed(removed);
        }
    }

    /// Overwrites `acc` with the total acceleration of every particle.
    /// Pinned particles get zero acceleration (they still exert forces on others).
    pub fn accelerations(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        pinned: &[bool],
        acc: &mut Vec<Vec3>,
    ) -> Result<()> {
        acc.clear();
        acc.resize(pos.len(), Vec3::ZERO);
        for (_, f) in &self.forces {
            f.accumulate(t, pos, vel, mass, acc)?;
        }
        for (a, &p) in acc.iter_mut().zip(pinned) {
            if p {
                *a = Vec3::ZERO;
            }
        }
        Ok(())
    }

    /// Overwrites `out` with the directional derivative of the total acceleration along
    /// `(dpos, dvel)`: analytic where a force provides it, central finite differences
    /// otherwise (relative accuracy ~1e-10). Pinned particles get zero.
    #[allow(clippy::too_many_arguments)]
    pub fn jacobian_vector(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        pinned: &[bool],
        dpos: &[Vec3],
        dvel: &[Vec3],
        out: &mut Vec<Vec3>,
    ) -> Result<()> {
        let n = pos.len();
        out.clear();
        out.resize(n, Vec3::ZERO);
        let mut fd: Option<FiniteDifference> = None;
        for (_, f) in &self.forces {
            if f.jacobian_vector(t, pos, vel, mass, dpos, dvel, out)? {
                continue;
            }
            let fd = fd.get_or_insert_with(|| FiniteDifference::new(pos, vel, dpos, dvel));
            fd.add(f.as_ref(), t, mass, out)?;
        }
        for (o, &p) in out.iter_mut().zip(pinned) {
            if p {
                *o = Vec3::ZERO;
            }
        }
        Ok(())
    }

    /// Sum of the potentials of all conservative forces (non-conservative ones are skipped).
    pub fn potential(&self, t: f64, pos: &[Vec3], mass: &[f64]) -> Result<f64> {
        let mut total = 0.0;
        for (_, f) in &self.forces {
            total += f.potential(t, pos, mass)?.unwrap_or(0.0);
        }
        Ok(total)
    }
}

/// Central differences of single forces along one direction, with shared perturbed states.
struct FiniteDifference {
    eps: f64,
    plus: (Vec<Vec3>, Vec<Vec3>),
    minus: (Vec<Vec3>, Vec<Vec3>),
    a_plus: Vec<Vec3>,
    a_minus: Vec<Vec3>,
}

impl FiniteDifference {
    fn new(pos: &[Vec3], vel: &[Vec3], dpos: &[Vec3], dvel: &[Vec3]) -> Self {
        let max = |v: &[Vec3]| v.iter().map(|x| x.norm()).fold(0.0, f64::max);
        let size = 1.0 + max(pos).max(max(vel));
        let dir = max(dpos).max(max(dvel));
        // ~cbrt(machine epsilon) relative step, optimal for central differences.
        let eps = if dir > 0.0 { 6e-6 * size / dir } else { 0.0 };
        let shift = |sign: f64| {
            (
                pos.iter()
                    .zip(dpos)
                    .map(|(x, d)| *x + *d * (sign * eps))
                    .collect(),
                vel.iter()
                    .zip(dvel)
                    .map(|(v, d)| *v + *d * (sign * eps))
                    .collect(),
            )
        };
        Self {
            eps,
            plus: shift(1.0),
            minus: shift(-1.0),
            a_plus: Vec::new(),
            a_minus: Vec::new(),
        }
    }

    fn add(&mut self, f: &dyn Force, t: f64, mass: &[f64], out: &mut [Vec3]) -> Result<()> {
        if self.eps == 0.0 {
            return Ok(());
        }
        let n = out.len();
        for (buf, (p, v)) in [
            (&mut self.a_plus, &self.plus),
            (&mut self.a_minus, &self.minus),
        ] {
            buf.clear();
            buf.resize(n, Vec3::ZERO);
            f.accumulate(t, p, v, mass, buf)?;
        }
        let scale = 0.5 / self.eps;
        for ((o, a), b) in out.iter_mut().zip(&self.a_plus).zip(&self.a_minus) {
            *o += (*a - *b) * scale;
        }
        Ok(())
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

/// Forces defined as a force (rather than an acceleration) have no meaning on a massless particle.
fn check_massive(i: usize, mass: &[f64], what: &str) -> Result<()> {
    if mass[i] == 0.0 {
        return invalid(format!("{what}: acts on massless particle {i}"));
    }
    Ok(())
}

fn scalar(name: &str, value: Param) -> Result<f64> {
    match value {
        Param::Scalar(x) if x.is_finite() => Ok(x),
        Param::Scalar(x) => invalid(format!("{name} must be finite, got {x}")),
        Param::Vector(_) => invalid(format!("{name} is a scalar, got a vector")),
    }
}

fn non_negative(name: &str, value: Param) -> Result<f64> {
    let x = scalar(name, value)?;
    if x < 0.0 {
        return invalid(format!("{name} must be non-negative, got {x}"));
    }
    Ok(x)
}

fn vector(name: &str, value: Param) -> Result<Vec3> {
    match value {
        Param::Vector(v) => Ok(v),
        Param::Scalar(_) => invalid(format!("{name} is a 3-vector, got a scalar")),
    }
}

/// Directional derivative of the spring force `k (1 - L/r) d` along `δd`.
fn spring_jvp(k: f64, rest_length: f64, d: Vec3, r: f64, dd: Vec3) -> Vec3 {
    dd * (k * (1.0 - rest_length / r)) + d * (k * rest_length * d.dot(dd) / (r * r * r))
}

/// Shifts index `i` down if it is above a removed particle.
fn shift(i: &mut usize, removed: usize) {
    if *i > removed {
        *i -= 1;
    }
}

/// Constant acceleration field, e.g. surface gravity `g = (0, -9.81, 0)`.
#[derive(Debug, Clone)]
pub struct UniformField {
    pub g: Vec3,
}

impl Force for UniformField {
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
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        parallel::for_each_indexed(acc, |_, a| *a += self.g);
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

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::UniformField(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![("g", Param::Vector(self.g))]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "g" => self.g = vector(name, value)?,
            _ => return invalid(format!("UniformField has no parameter {name:?}")),
        }
        Ok(())
    }
}

/// Pairwise Newtonian gravity with Plummer softening `eps`:
/// `U = -G m_i m_j / sqrt(r^2 + eps^2)`. Direct O(N^2) summation.
/// Massless particles feel gravity but do not source it (test particles).
#[derive(Debug, Clone)]
pub struct NewtonianGravity {
    pub g: f64,
    pub softening: f64,
}

impl NewtonianGravity {
    /// Adds the interactions of pairs `(i, j)` with `lo <= i < hi`, `i < j`.
    /// `acc` holds particles `n - acc.len()..n`: the whole system, or `lo..n` for a block buffer.
    fn accumulate_rows(&self, lo: usize, hi: usize, pos: &[Vec3], mass: &[f64], acc: &mut [Vec3]) {
        let eps2 = self.softening * self.softening;
        let offset = pos.len() - acc.len();
        for i in lo..hi {
            let mut ai = Vec3::ZERO;
            for j in (i + 1)..pos.len() {
                let d = pos[j] - pos[i];
                let r2 = d.norm_squared() + eps2;
                let inv_r3 = 1.0 / (r2 * r2.sqrt());
                let s = self.g * inv_r3;
                ai += d * (s * mass[j]);
                acc[j - offset] -= d * (s * mass[i]);
            }
            acc[i - offset] += ai;
        }
    }

    fn potential_rows(&self, lo: usize, hi: usize, pos: &[Vec3], mass: &[f64]) -> f64 {
        let eps2 = self.softening * self.softening;
        let mut u = 0.0;
        for i in lo..hi {
            for j in (i + 1)..pos.len() {
                let r = ((pos[j] - pos[i]).norm_squared() + eps2).sqrt();
                u -= self.g * mass[i] * mass[j] / r;
            }
        }
        u
    }
}

impl Force for NewtonianGravity {
    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        // δ(d/s³) = δd/s³ - 3 d (d·δd)/s⁵ for d = x_j - x_i, s² = |d|² + ε².
        let eps2 = self.softening * self.softening;
        for i in 0..pos.len() {
            for j in (i + 1)..pos.len() {
                let d = pos[j] - pos[i];
                let dd = dpos[j] - dpos[i];
                let s2 = d.norm_squared() + eps2;
                let inv_s3 = 1.0 / (s2 * s2.sqrt());
                let w = (dd - d * (3.0 * d.dot(dd) / s2)) * (self.g * inv_s3);
                out[i] += w * mass[j];
                out[j] -= w * mass[i];
            }
        }
        Ok(true)
    }

    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let blocks = parallel::pair_blocks(pos.len());
        if blocks.len() == 1 {
            self.accumulate_rows(0, pos.len(), pos, mass, acc);
            return Ok(());
        }
        // Each block of rows writes its own buffer (every pair is computed once, i < j);
        // the buffers are then added in block order, so results do not depend on threads.
        // Rows i >= lo only touch particles >= lo, so a block's buffer covers lo..n.
        let partials = parallel::map_blocks(&blocks, |lo, hi| {
            let mut buf = vec![Vec3::ZERO; pos.len() - lo];
            self.accumulate_rows(lo, hi, pos, mass, &mut buf);
            (lo, buf)
        });
        parallel::add_partials(acc, &partials);
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        let blocks = parallel::pair_blocks(pos.len());
        let partials =
            parallel::map_blocks(&blocks, |lo, hi| self.potential_rows(lo, hi, pos, mass));
        Ok(Some(partials.iter().sum()))
    }

    fn name(&self) -> String {
        "NewtonianGravity".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::NewtonianGravity(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("G", Param::Scalar(self.g)),
            ("softening", Param::Scalar(self.softening)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "G" => self.g = scalar(name, value)?,
            "softening" => self.softening = non_negative(name, value)?,
            _ => return invalid(format!("NewtonianGravity has no parameter {name:?}")),
        }
        Ok(())
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
    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        check_index(self.i, pos.len(), "Spring")?;
        check_index(self.j, pos.len(), "Spring")?;
        let d = pos[self.j] - pos[self.i];
        let r = d.norm();
        if r > 0.0 {
            let dd = dpos[self.j] - dpos[self.i];
            let df = spring_jvp(self.k, self.rest_length, d, r, dd);
            out[self.i] += df / mass[self.i];
            out[self.j] -= df / mass[self.j];
        }
        Ok(true)
    }

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
        check_massive(self.i, mass, "Spring")?;
        check_massive(self.j, mass, "Spring")?;
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

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::Spring(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("k", Param::Scalar(self.k)),
            ("rest_length", Param::Scalar(self.rest_length)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "k" => self.k = scalar(name, value)?,
            "rest_length" => self.rest_length = non_negative(name, value)?,
            _ => return invalid(format!("Spring has no parameter {name:?}")),
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
    fn jacobian_vector(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        dpos: &[Vec3],
        _dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        check_index(self.i, pos.len(), "AnchorSpring")?;
        let d = self.anchor - pos[self.i];
        let dd = -dpos[self.i];
        let df = if self.rest_length == 0.0 {
            dd * self.k
        } else {
            let r = d.norm();
            if r == 0.0 {
                return Ok(true);
            }
            spring_jvp(self.k, self.rest_length, d, r, dd)
        };
        out[self.i] += df / mass[self.i];
        Ok(true)
    }

    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        check_index(self.i, pos.len(), "AnchorSpring")?;
        check_massive(self.i, mass, "AnchorSpring")?;
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

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::AnchorSpring(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("anchor", Param::Vector(self.anchor)),
            ("k", Param::Scalar(self.k)),
            ("rest_length", Param::Scalar(self.rest_length)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "anchor" => self.anchor = vector(name, value)?,
            "k" => self.k = scalar(name, value)?,
            "rest_length" => self.rest_length = non_negative(name, value)?,
            _ => return invalid(format!("AnchorSpring has no parameter {name:?}")),
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

/// Linear (Stokes) drag `a = -gamma v`; `gamma` is a rate (1/time), independent of mass.
#[derive(Debug, Clone)]
pub struct LinearDrag {
    pub gamma: f64,
}

impl Force for LinearDrag {
    fn jacobian_vector(
        &self,
        _t: f64,
        _pos: &[Vec3],
        _vel: &[Vec3],
        _mass: &[f64],
        _dpos: &[Vec3],
        dvel: &[Vec3],
        out: &mut [Vec3],
    ) -> Result<bool> {
        for (o, d) in out.iter_mut().zip(dvel) {
            *o -= *d * self.gamma;
        }
        Ok(true)
    }

    fn accumulate(
        &self,
        _t: f64,
        _pos: &[Vec3],
        vel: &[Vec3],
        _mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        parallel::for_each_indexed(acc, |i, a| *a -= vel[i] * self.gamma);
        Ok(())
    }

    fn velocity_dependent(&self) -> bool {
        true
    }

    fn name(&self) -> String {
        "LinearDrag".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::LinearDrag(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![("gamma", Param::Scalar(self.gamma))]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "gamma" => self.gamma = scalar(name, value)?,
            _ => return invalid(format!("LinearDrag has no parameter {name:?}")),
        }
        Ok(())
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
        for (i, ((a, v), m)) in acc.iter_mut().zip(vel).zip(mass).enumerate() {
            if *m == 0.0 {
                if *v == Vec3::ZERO {
                    continue;
                }
                return invalid(format!(
                    "QuadraticDrag: acts on moving massless particle {i}"
                ));
            }
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

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::QuadraticDrag(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![("c", Param::Scalar(self.c))]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "c" => self.c = scalar(name, value)?,
            _ => return invalid(format!("QuadraticDrag has no parameter {name:?}")),
        }
        Ok(())
    }
}
