//! Time integrators. Each owns its scratch buffers so stepping does not allocate.
//! Implement [`Integrator`] to add a new scheme.

use crate::constraints::Constraints;
use crate::error::{invalid, Result};
use crate::forces::ForceSet;
use crate::state::State;
use crate::vec3::Vec3;

pub trait Integrator: Send + Sync {
    /// Advance `state` (including `state.t`) by `dt`.
    fn step(&mut self, state: &mut State, forces: &ForceSet, dt: f64) -> Result<()>;
    fn name(&self) -> &'static str;
    /// Global order of accuracy.
    fn order(&self) -> u32;
    /// Symplectic for velocity-independent forces (bounded long-time energy error).
    fn symplectic(&self) -> bool;

    /// Advance `state` by `dt` subject to distance `constraints` (see
    /// [`crate::constraints`]), writing each constraint's tension at the end of the step to
    /// `tension`. Only schemes built from velocity Verlet support this.
    fn step_constrained(
        &mut self,
        _state: &mut State,
        _forces: &ForceSet,
        _constraints: &Constraints,
        _tension: &mut Vec<f64>,
        _dt: f64,
    ) -> Result<()> {
        invalid(format!(
            "the {} integrator does not support constraints; use verlet or yoshida4",
            self.name()
        ))
    }
}

/// Names accepted by [`by_name`].
pub const NAMES: &[&str] = &[
    "explicit_euler",
    "symplectic_euler",
    "verlet",
    "yoshida4",
    "rk4",
];

pub fn by_name(name: &str) -> Result<Box<dyn Integrator>> {
    Ok(match name {
        "explicit_euler" | "euler" => Box::new(ExplicitEuler::default()),
        "symplectic_euler" => Box::new(SymplecticEuler::default()),
        "verlet" | "velocity_verlet" | "leapfrog" => Box::new(VelocityVerlet::default()),
        "yoshida4" => Box::new(Yoshida4::default()),
        "rk4" => Box::new(Rk4::default()),
        _ => {
            return invalid(format!(
                "unknown integrator {name:?}; expected one of {NAMES:?}"
            ))
        }
    })
}

/// Forward Euler. First order, not symplectic: energy drifts. Included as a baseline.
#[derive(Default)]
pub struct ExplicitEuler {
    acc: Vec<Vec3>,
}

impl Integrator for ExplicitEuler {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        forces.accelerations(s.t, &s.pos, &s.vel, &s.mass, &s.pinned, &mut self.acc)?;
        for ((x, v), a) in s.pos.iter_mut().zip(s.vel.iter_mut()).zip(&self.acc) {
            *x += *v * dt;
            *v += *a * dt;
        }
        s.t += dt;
        Ok(())
    }
    fn name(&self) -> &'static str {
        "explicit_euler"
    }
    fn order(&self) -> u32 {
        1
    }
    fn symplectic(&self) -> bool {
        false
    }
}

/// Semi-implicit (symplectic) Euler: kick then drift. First order, symplectic.
#[derive(Default)]
pub struct SymplecticEuler {
    acc: Vec<Vec3>,
}

impl Integrator for SymplecticEuler {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        forces.accelerations(s.t, &s.pos, &s.vel, &s.mass, &s.pinned, &mut self.acc)?;
        for ((x, v), a) in s.pos.iter_mut().zip(s.vel.iter_mut()).zip(&self.acc) {
            *v += *a * dt;
            *x += *v * dt;
        }
        s.t += dt;
        Ok(())
    }
    fn name(&self) -> &'static str {
        "symplectic_euler"
    }
    fn order(&self) -> u32 {
        1
    }
    fn symplectic(&self) -> bool {
        true
    }
}

/// Acceleration from the last evaluation, with a snapshot of the inputs that produced it.
///
/// Kick-drift-kick schemes evaluate the acceleration at the end of one (sub)step and again,
/// at the same positions and time, at the start of the next. For velocity-independent forces
/// the cache returns the stored result instead, halving the cost of Verlet. Validity is
/// checked against the actual inputs (bit-for-bit) and the force set's version, so external
/// edits to the state or the forces are always picked up.
#[derive(Default)]
struct AccelCache {
    acc: Vec<Vec3>,
    valid: bool,
    forces_version: u64,
    t: f64,
    pos: Vec<Vec3>,
    mass: Vec<f64>,
    pinned: Vec<bool>,
}

fn same_bits(a: &[Vec3], b: &[Vec3]) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|(p, q)| {
            p.x.to_bits() == q.x.to_bits()
                && p.y.to_bits() == q.y.to_bits()
                && p.z.to_bits() == q.z.to_bits()
        })
}

impl AccelCache {
    fn hit(&self, s: &State, forces: &ForceSet) -> bool {
        self.valid
            && self.forces_version == forces.version()
            && self.t.to_bits() == s.t.to_bits()
            && same_bits(&self.pos, &s.pos)
            && self.pinned == s.pinned
            && self.mass.len() == s.mass.len()
            && self
                .mass
                .iter()
                .zip(&s.mass)
                .all(|(a, b)| a.to_bits() == b.to_bits())
    }

    /// Acceleration for the current state, reused when the inputs are unchanged.
    fn get(&mut self, s: &State, forces: &ForceSet) -> Result<&[Vec3]> {
        if !self.hit(s, forces) {
            self.valid = false;
            forces.accelerations(s.t, &s.pos, &s.vel, &s.mass, &s.pinned, &mut self.acc)?;
            // Velocity-dependent accelerations change with every kick, so never reuse them.
            if !forces.velocity_dependent() {
                self.valid = true;
                self.forces_version = forces.version();
                self.t = s.t;
                self.pos.clone_from(&s.pos);
                self.mass.clone_from(&s.mass);
                self.pinned.clone_from(&s.pinned);
            }
        }
        Ok(&self.acc)
    }
}

/// One kick-drift-kick leapfrog substep of length `h`.
fn kdk(s: &mut State, forces: &ForceSet, cache: &mut AccelCache, h: f64) -> Result<()> {
    let acc = cache.get(s, forces)?;
    for ((x, v), a) in s.pos.iter_mut().zip(s.vel.iter_mut()).zip(acc) {
        *v += *a * (0.5 * h);
        *x += *v * h;
    }
    s.t += h;
    let acc = cache.get(s, forces)?;
    for (v, a) in s.vel.iter_mut().zip(acc) {
        *v += *a * (0.5 * h);
    }
    Ok(())
}

/// One RATTLE substep of length `h`: [`kdk`] with the positions projected onto the
/// constraints after the drift and the velocities after the final kick.
fn rattle(
    s: &mut State,
    forces: &ForceSet,
    constraints: &Constraints,
    cache: &mut AccelCache,
    old: &mut Vec<Vec3>,
    tension: &mut Vec<f64>,
    h: f64,
) -> Result<()> {
    old.clone_from(&s.pos);
    let acc = cache.get(s, forces)?;
    for ((x, v), a) in s.pos.iter_mut().zip(s.vel.iter_mut()).zip(acc) {
        *v += *a * (0.5 * h);
        *x += *v * h;
    }
    constraints.project_positions(s, old, h)?;
    s.t += h;
    let acc = cache.get(s, forces)?;
    for (v, a) in s.vel.iter_mut().zip(acc) {
        *v += *a * (0.5 * h);
    }
    constraints.project_velocities(s, h, tension)
}

/// Velocity Verlet (kick-drift-kick leapfrog). Second order, symplectic, time-reversible.
/// One force evaluation per step for velocity-independent forces (the end-of-step
/// acceleration is reused), two otherwise.
/// With constraints this is RATTLE.
#[derive(Default)]
pub struct VelocityVerlet {
    cache: AccelCache,
    old: Vec<Vec3>,
}

impl Integrator for VelocityVerlet {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        kdk(s, forces, &mut self.cache, dt)
    }
    fn step_constrained(
        &mut self,
        s: &mut State,
        forces: &ForceSet,
        constraints: &Constraints,
        tension: &mut Vec<f64>,
        dt: f64,
    ) -> Result<()> {
        rattle(
            s,
            forces,
            constraints,
            &mut self.cache,
            &mut self.old,
            tension,
            dt,
        )
    }
    fn name(&self) -> &'static str {
        "verlet"
    }
    fn order(&self) -> u32 {
        2
    }
    fn symplectic(&self) -> bool {
        true
    }
}

/// Yoshida (1990) fourth-order symplectic composition of three leapfrog substeps.
/// Three force evaluations per step for velocity-independent forces, six otherwise.
/// With constraints, the same composition of RATTLE substeps.
#[derive(Default)]
pub struct Yoshida4 {
    cache: AccelCache,
    old: Vec<Vec3>,
}

/// Triple-jump weights `[w1, w0, w1]`.
fn yoshida_weights() -> [f64; 3] {
    let cbrt2 = 2f64.cbrt();
    let w1 = 1.0 / (2.0 - cbrt2);
    let w0 = -cbrt2 / (2.0 - cbrt2);
    [w1, w0, w1]
}

impl Integrator for Yoshida4 {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        for w in yoshida_weights() {
            kdk(s, forces, &mut self.cache, w * dt)?;
        }
        Ok(())
    }
    fn step_constrained(
        &mut self,
        s: &mut State,
        forces: &ForceSet,
        constraints: &Constraints,
        tension: &mut Vec<f64>,
        dt: f64,
    ) -> Result<()> {
        for w in yoshida_weights() {
            rattle(
                s,
                forces,
                constraints,
                &mut self.cache,
                &mut self.old,
                tension,
                w * dt,
            )?;
        }
        Ok(())
    }
    fn name(&self) -> &'static str {
        "yoshida4"
    }
    fn order(&self) -> u32 {
        4
    }
    fn symplectic(&self) -> bool {
        true
    }
}

/// Classical fourth-order Runge-Kutta on (x, v). Not symplectic, but treats
/// velocity-dependent and explicitly time-dependent forces consistently.
#[derive(Default)]
pub struct Rk4 {
    x0: Vec<Vec3>,
    v0: Vec<Vec3>,
    xs: Vec<Vec3>,
    vs: Vec<Vec3>,
    acc: Vec<Vec3>,
    dx: Vec<Vec3>,
    dv: Vec<Vec3>,
}

impl Integrator for Rk4 {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        const C: [f64; 4] = [0.0, 0.5, 0.5, 1.0];
        const W: [f64; 4] = [1.0, 2.0, 2.0, 1.0];
        let n = s.len();
        self.x0.clone_from(&s.pos);
        self.v0.clone_from(&s.vel);
        self.xs.clone_from(&s.pos);
        self.vs.clone_from(&s.vel);
        self.dx.clear();
        self.dx.resize(n, Vec3::ZERO);
        self.dv.clear();
        self.dv.resize(n, Vec3::ZERO);

        for stage in 0..4 {
            // (xs, vs) holds the previous stage's slopes (k_x = v, k_v = a) on entry.
            if stage > 0 {
                let h = C[stage] * dt;
                for p in 0..n {
                    let kx = self.vs[p];
                    let kv = self.acc[p];
                    self.xs[p] = self.x0[p] + kx * h;
                    self.vs[p] = self.v0[p] + kv * h;
                }
            }
            forces.accelerations(
                s.t + C[stage] * dt,
                &self.xs,
                &self.vs,
                &s.mass,
                &s.pinned,
                &mut self.acc,
            )?;
            for p in 0..n {
                self.dx[p] += self.vs[p] * W[stage];
                self.dv[p] += self.acc[p] * W[stage];
            }
        }

        for p in 0..n {
            s.pos[p] = self.x0[p] + self.dx[p] * (dt / 6.0);
            s.vel[p] = self.v0[p] + self.dv[p] * (dt / 6.0);
        }
        s.t += dt;
        Ok(())
    }
    fn name(&self) -> &'static str {
        "rk4"
    }
    fn order(&self) -> u32 {
        4
    }
    fn symplectic(&self) -> bool {
        false
    }
}
