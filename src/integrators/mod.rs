//! Time integrators. Each owns its scratch buffers so stepping does not allocate.
//! Implement [`Integrator`] to add a new scheme.

use crate::constraints::Constraints;
use crate::error::{invalid, Result};
use crate::forces::ForceSet;
use crate::state::State;
use crate::vec3::Vec3;

mod boris;
mod gauss;
mod splitting;
mod wisdom_holman;

pub use boris::Boris;
pub use gauss::GaussLegendre;
pub use splitting::{Composition, Op, Scheme, Splitting};
pub use wisdom_holman::WisdomHolman;

pub trait Integrator: Send + Sync {
    /// Advance `state` (including `state.t`) by `dt`.
    fn step(&mut self, state: &mut State, forces: &ForceSet, dt: f64) -> Result<()>;
    fn name(&self) -> &str;
    /// Global order of accuracy.
    fn order(&self) -> u32;
    /// Symplectic for velocity-independent forces (bounded long-time energy error).
    fn symplectic(&self) -> bool;

    /// For user-defined schemes, the coefficients needed to rebuild them (stored in
    /// checkpoints). Built-in integrators are rebuilt from their name and return `None`.
    fn scheme(&self) -> Option<Scheme> {
        None
    }

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
    "dopri5",
    "yoshida6",
    "yoshida8",
    "pefrl",
    "blanes_moan4",
    "gauss2",
    "gauss4",
    "gauss6",
    "wisdom_holman",
    "boris",
];

pub fn by_name(name: &str) -> Result<Box<dyn Integrator>> {
    Ok(match name {
        "explicit_euler" | "euler" => Box::new(ExplicitEuler::default()),
        "symplectic_euler" => Box::new(SymplecticEuler::default()),
        "verlet" | "velocity_verlet" | "leapfrog" => Box::new(VelocityVerlet::default()),
        "yoshida4" | "forest_ruth" => Box::new(Yoshida4::default()),
        "yoshida6" => Box::new(Composition::yoshida6()),
        "yoshida8" => Box::new(Composition::yoshida8()),
        "pefrl" | "omelyan" => Box::new(Splitting::pefrl()),
        "blanes_moan4" => Box::new(Splitting::blanes_moan4()),
        "gauss2" | "implicit_midpoint" => Box::new(GaussLegendre::order2()),
        "gauss4" => Box::new(GaussLegendre::order4()),
        "gauss6" => Box::new(GaussLegendre::order6()),
        "wisdom_holman" => Box::new(WisdomHolman::default()),
        "boris" => Box::new(Boris::default()),
        "rk4" => Box::new(Rk4::default()),
        "dopri5" | "dormand_prince" => Box::new(Dopri5::default()),
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
        forces.accelerations(
            s.t,
            &s.pos,
            &s.vel,
            &s.mass,
            &s.charge,
            &s.radius,
            &s.pinned,
            &mut self.acc,
        )?;
        for ((x, v), a) in s.pos.iter_mut().zip(s.vel.iter_mut()).zip(&self.acc) {
            *x += *v * dt;
            *v += *a * dt;
        }
        s.t += dt;
        Ok(())
    }
    fn name(&self) -> &str {
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
        forces.accelerations(
            s.t,
            &s.pos,
            &s.vel,
            &s.mass,
            &s.charge,
            &s.radius,
            &s.pinned,
            &mut self.acc,
        )?;
        for ((x, v), a) in s.pos.iter_mut().zip(s.vel.iter_mut()).zip(&self.acc) {
            *v += *a * dt;
            *x += *v * dt;
        }
        s.t += dt;
        Ok(())
    }
    fn name(&self) -> &str {
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
/// at the same positions and time, at the start of the next; first-same-as-last Runge-Kutta
/// schemes do the same with the full state. The cache returns the stored result instead.
/// Validity is checked against the actual inputs (bit-for-bit, velocities included when a
/// force depends on them) and the force set's version, so external edits to the state or
/// the forces are always picked up.
#[derive(Default)]
struct AccelCache {
    acc: Vec<Vec3>,
    valid: bool,
    forces_version: u64,
    t: f64,
    pos: Vec<Vec3>,
    /// Only compared when the forces are velocity-dependent.
    vel: Vec<Vec3>,
    mass: Vec<f64>,
    charge: Vec<f64>,
    radius: Vec<f64>,
    pinned: Vec<bool>,
    /// Number of times the forces were actually evaluated.
    evaluations: u64,
}

fn same_scalars(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
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
    fn hit(&self, at: &At<'_>, forces: &ForceSet) -> bool {
        self.valid
            && self.forces_version == forces.version()
            && self.t.to_bits() == at.t.to_bits()
            && same_bits(&self.pos, at.pos)
            && (!forces.velocity_dependent() || same_bits(&self.vel, at.vel))
            && self.pinned == at.pinned
            && self.mass.len() == at.mass.len()
            && same_scalars(&self.mass, at.mass)
            && same_scalars(&self.charge, at.charge)
            && same_scalars(&self.radius, at.radius)
    }

    /// Acceleration for the current state, reused when the inputs are unchanged.
    fn get(&mut self, s: &State, forces: &ForceSet) -> Result<&[Vec3]> {
        self.get_at(
            &At {
                t: s.t,
                pos: &s.pos,
                vel: &s.vel,
                mass: &s.mass,
                charge: &s.charge,
                radius: &s.radius,
                pinned: &s.pinned,
            },
            forces,
        )
    }

    /// Acceleration at `at`, reused when the inputs are unchanged.
    fn get_at(&mut self, at: &At<'_>, forces: &ForceSet) -> Result<&[Vec3]> {
        if !self.hit(at, forces) {
            self.valid = false;
            self.evaluations += 1;
            forces.accelerations(
                at.t,
                at.pos,
                at.vel,
                at.mass,
                at.charge,
                at.radius,
                at.pinned,
                &mut self.acc,
            )?;
            self.valid = true;
            self.forces_version = forces.version();
            self.t = at.t;
            copy_into(&mut self.pos, at.pos);
            if forces.velocity_dependent() {
                copy_into(&mut self.vel, at.vel);
            }
            copy_into(&mut self.mass, at.mass);
            copy_into(&mut self.charge, at.charge);
            copy_into(&mut self.radius, at.radius);
            copy_into(&mut self.pinned, at.pinned);
        }
        Ok(&self.acc)
    }
}

fn copy_into<T: Copy>(dst: &mut Vec<T>, src: &[T]) {
    dst.clear();
    dst.extend_from_slice(src);
}

/// The inputs of one force evaluation.
struct At<'a> {
    t: f64,
    pos: &'a [Vec3],
    vel: &'a [Vec3],
    mass: &'a [f64],
    charge: &'a [f64],
    radius: &'a [f64],
    pinned: &'a [bool],
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
    fn name(&self) -> &str {
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
    fn name(&self) -> &str {
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
                &s.charge,
                &s.radius,
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
    fn name(&self) -> &str {
        "rk4"
    }
    fn order(&self) -> u32 {
        4
    }
    fn symplectic(&self) -> bool {
        false
    }
}

// Dormand-Prince 5(4) tableau (Hairer, Norsett & Wanner, Solving ODEs I, table II.5.2).
const DP_C: [f64; 7] = [0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
const DP_A: [[f64; 6]; 7] = [
    [0.0; 6],
    [0.2, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    // The fifth-order weights: stage 7 is evaluated at the new state (first same as last).
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];
/// Fifth-order minus embedded fourth-order weights: the local error estimate.
const DP_E: [f64; 7] = [
    71.0 / 57600.0,
    0.0,
    -71.0 / 16695.0,
    71.0 / 1920.0,
    -17253.0 / 339200.0,
    22.0 / 525.0,
    -1.0 / 40.0,
];
/// Weights of the fourth-order continuous extension (dense output).
const DP_D: [f64; 7] = [
    -12715105075.0 / 11282082432.0,
    0.0,
    87487479700.0 / 32700410799.0,
    -10690763975.0 / 1880347072.0,
    701980252875.0 / 199316789632.0,
    -1453857185.0 / 822651844.0,
    69997945.0 / 29380423.0,
];

/// Dormand-Prince 5(4): the embedded Runge-Kutta pair behind most adaptive ODE solvers.
/// Fifth order, not symplectic, handles velocity-dependent forces. Six force evaluations
/// per step (the last stage is reused as the first of the next step).
///
/// As a fixed-step integrator it simply takes the fifth-order solution.
/// [`World::run_adaptive`](crate::World::run_adaptive) also uses its embedded error estimate
/// to choose the step and its continuous extension to output at arbitrary times.
#[derive(Default)]
pub struct Dopri5 {
    cache: AccelCache,
    /// Force evaluations made outside the cache.
    direct_evaluations: u64,
    t0: f64,
    h: f64,
    x0: Vec<Vec3>,
    v0: Vec<Vec3>,
    /// Stage slopes: `kx[i]` is the velocity and `kv[i]` the acceleration at stage `i`.
    kx: [Vec<Vec3>; 7],
    kv: [Vec<Vec3>; 7],
    /// Stage state; after [`Dopri5::attempt`], the fifth-order solution.
    xs: Vec<Vec3>,
    vs: Vec<Vec3>,
    /// Coefficients of the continuous extension of the last attempted step.
    dense_x: [Vec<Vec3>; 5],
    dense_v: [Vec<Vec3>; 5],
}

impl Dopri5 {
    /// Force evaluations so far.
    pub(crate) fn evaluations(&self) -> u64 {
        self.cache.evaluations + self.direct_evaluations
    }

    /// Starts a step from `s`, evaluating the forces there unless the last step ended here.
    pub(crate) fn begin(&mut self, s: &State, forces: &ForceSet) -> Result<()> {
        self.t0 = s.t;
        copy_into(&mut self.x0, &s.pos);
        copy_into(&mut self.v0, &s.vel);
        copy_into(&mut self.kx[0], &s.vel);
        let acc = self.cache.get(s, forces)?;
        copy_into(&mut self.kv[0], acc);
        Ok(())
    }

    /// Computes the step of size `h` from the state given to [`Dopri5::begin`], leaving
    /// the fifth-order solution in the stage buffers (see [`Dopri5::commit`]).
    pub(crate) fn attempt(
        &mut self,
        forces: &ForceSet,
        mass: &[f64],
        charge: &[f64],
        radius: &[f64],
        pinned: &[bool],
        h: f64,
    ) -> Result<()> {
        self.h = h;
        let n = self.x0.len();
        for i in 1..7 {
            copy_into(&mut self.xs, &self.x0);
            copy_into(&mut self.vs, &self.v0);
            for (j, &a) in DP_A[i][..i].iter().enumerate() {
                if a == 0.0 {
                    continue;
                }
                let w = a * h;
                for p in 0..n {
                    self.xs[p] += self.kx[j][p] * w;
                    self.vs[p] += self.kv[j][p] * w;
                }
            }
            copy_into(&mut self.kx[i], &self.vs);
            let t = self.t0 + DP_C[i] * h;
            if i < 6 {
                self.direct_evaluations += 1;
                forces.accelerations(
                    t,
                    &self.xs,
                    &self.vs,
                    mass,
                    charge,
                    radius,
                    pinned,
                    &mut self.kv[i],
                )?;
            } else {
                // The new state: cached, so the next step's first stage is free.
                let at = At {
                    t,
                    pos: &self.xs,
                    vel: &self.vs,
                    mass,
                    charge,
                    radius,
                    pinned,
                };
                let acc = self.cache.get_at(&at, forces)?;
                copy_into(&mut self.kv[6], acc);
            }
        }
        Ok(())
    }

    /// Scaled RMS norm of the local error estimate of the last attempt (Hairer's norm):
    /// each component's error divided by `atol + rtol * max(|start|, |end|)`.
    /// Infinite if the step produced non-finite values.
    pub(crate) fn error_norm(&self, rtol: f64, atol: f64) -> f64 {
        let n = self.x0.len();
        if n == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        let mut add = |y0: f64, y1: f64, e: f64| {
            let scale = atol + rtol * y0.abs().max(y1.abs());
            let r = e / scale;
            sum += r * r;
        };
        for p in 0..n {
            let (mut ex, mut ev) = (Vec3::ZERO, Vec3::ZERO);
            for (j, &e) in DP_E.iter().enumerate() {
                if e != 0.0 {
                    ex += self.kx[j][p] * e;
                    ev += self.kv[j][p] * e;
                }
            }
            let (ex, ev) = (ex * self.h, ev * self.h);
            let (x0, x1, v0, v1) = (self.x0[p], self.xs[p], self.v0[p], self.vs[p]);
            add(x0.x, x1.x, ex.x);
            add(x0.y, x1.y, ex.y);
            add(x0.z, x1.z, ex.z);
            add(v0.x, v1.x, ev.x);
            add(v0.y, v1.y, ev.y);
            add(v0.z, v1.z, ev.z);
        }
        let norm = (sum / (6 * n) as f64).sqrt();
        if norm.is_finite() {
            norm
        } else {
            f64::INFINITY
        }
    }

    /// Moves the last attempted step's solution into `s`, with time `t` (normally
    /// `start + h`; pass the exact end time to avoid round-off on the final step).
    pub(crate) fn commit(&self, s: &mut State, t: f64) {
        copy_into(&mut s.pos, &self.xs);
        copy_into(&mut s.vel, &self.vs);
        s.t = t;
    }

    /// Prepares [`Dopri5::interpolate`] for the last attempted step.
    pub(crate) fn prepare_dense(&mut self) {
        let n = self.x0.len();
        let h = self.h;
        for d in self.dense_x.iter_mut().chain(self.dense_v.iter_mut()) {
            d.clear();
            d.resize(n, Vec3::ZERO);
        }
        for p in 0..n {
            for (y0, y1, k, dense) in [
                (self.x0[p], self.xs[p], &self.kx, &mut self.dense_x),
                (self.v0[p], self.vs[p], &self.kv, &mut self.dense_v),
            ] {
                let diff = y1 - y0;
                let b = k[0][p] * h - diff;
                let mut d = Vec3::ZERO;
                for (j, &w) in DP_D.iter().enumerate() {
                    if w != 0.0 {
                        d += k[j][p] * w;
                    }
                }
                dense[0][p] = y0;
                dense[1][p] = diff;
                dense[2][p] = b;
                dense[3][p] = diff - k[6][p] * h - b;
                dense[4][p] = d * h;
            }
        }
    }

    /// Position and velocity at fraction `theta` of the last attempted step, from the
    /// fourth-order continuous extension (call [`Dopri5::prepare_dense`] first).
    pub(crate) fn interpolate(&self, theta: f64, pos: &mut Vec<Vec3>, vel: &mut Vec<Vec3>) {
        let u = 1.0 - theta;
        let eval = |d: &[Vec<Vec3>; 5], p: usize| {
            d[0][p] + (d[1][p] + (d[2][p] + (d[3][p] + d[4][p] * u) * theta) * u) * theta
        };
        let n = self.x0.len();
        pos.clear();
        vel.clear();
        pos.extend((0..n).map(|p| eval(&self.dense_x, p)));
        vel.extend((0..n).map(|p| eval(&self.dense_v, p)));
    }

    /// Initial step guess (Hairer's algorithm): one extra force evaluation. `direction`
    /// is +1 or -1. Call after [`Dopri5::begin`].
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn initial_step(
        &mut self,
        forces: &ForceSet,
        mass: &[f64],
        charge: &[f64],
        radius: &[f64],
        pinned: &[bool],
        direction: f64,
        rtol: f64,
        atol: f64,
    ) -> Result<f64> {
        let n = self.x0.len();
        if n == 0 {
            return Ok(1.0);
        }
        let norm = |f: &dyn Fn(usize) -> [(f64, f64); 6]| {
            let mut sum = 0.0;
            for p in 0..n {
                for (y, v) in f(p) {
                    let r = v / (atol + rtol * y.abs());
                    sum += r * r;
                }
            }
            (sum / (6 * n) as f64).sqrt()
        };
        let comps = |a: Vec3, b: Vec3, c: Vec3, d: Vec3| {
            [
                (a.x, c.x),
                (a.y, c.y),
                (a.z, c.z),
                (b.x, d.x),
                (b.y, d.y),
                (b.z, d.z),
            ]
        };
        let (x0, v0, kx, kv) = (&self.x0, &self.v0, &self.kx[0], &self.kv[0]);
        let d0 = norm(&|p| comps(x0[p], v0[p], x0[p], v0[p]));
        let d1 = norm(&|p| comps(x0[p], v0[p], kx[p], kv[p]));
        let mut h = if d0 < 1e-10 || d1 < 1e-10 {
            1e-6
        } else {
            0.01 * d0 / d1
        };
        // An explicit Euler step to estimate the second derivative.
        self.xs.clear();
        self.vs.clear();
        for p in 0..n {
            self.xs.push(x0[p] + kx[p] * (h * direction));
            self.vs.push(v0[p] + kv[p] * (h * direction));
        }
        self.direct_evaluations += 1;
        forces.accelerations(
            self.t0 + h * direction,
            &self.xs,
            &self.vs,
            mass,
            charge,
            radius,
            pinned,
            &mut self.kv[1],
        )?;
        let (x0, v0, kx, kv, kv1, vs) = (
            &self.x0,
            &self.v0,
            &self.kx[0],
            &self.kv[0],
            &self.kv[1],
            &self.vs,
        );
        let d2 = norm(&|p| comps(x0[p], v0[p], vs[p] - kx[p], kv1[p] - kv[p])) / h;
        let der = d1.max(d2);
        let h1 = if der <= 1e-15 {
            (h * 1e-3).max(1e-6)
        } else {
            (0.01 / der).powf(0.2)
        };
        h = (100.0 * h).min(h1);
        if h.is_finite() && h > 0.0 {
            Ok(h)
        } else {
            invalid("could not choose an initial step: the forces are not finite")
        }
    }
}

impl Integrator for Dopri5 {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        self.begin(s, forces)?;
        self.attempt(forces, &s.mass, &s.charge, &s.radius, &s.pinned, dt)?;
        self.commit(s, s.t + dt);
        Ok(())
    }
    fn name(&self) -> &str {
        "dopri5"
    }
    fn order(&self) -> u32 {
        5
    }
    fn symplectic(&self) -> bool {
        false
    }
}
