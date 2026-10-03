//! Time integrators. Each owns its scratch buffers so stepping does not allocate.
//! Implement [`Integrator`] to add a new scheme.

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
        forces.accelerations(s.t, &s.pos, &s.vel, &s.mass, &mut self.acc)?;
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
        forces.accelerations(s.t, &s.pos, &s.vel, &s.mass, &mut self.acc)?;
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

/// One kick-drift-kick leapfrog substep of length `h`.
fn kdk(s: &mut State, forces: &ForceSet, acc: &mut Vec<Vec3>, h: f64) -> Result<()> {
    forces.accelerations(s.t, &s.pos, &s.vel, &s.mass, acc)?;
    for ((x, v), a) in s.pos.iter_mut().zip(s.vel.iter_mut()).zip(acc.iter()) {
        *v += *a * (0.5 * h);
        *x += *v * h;
    }
    s.t += h;
    forces.accelerations(s.t, &s.pos, &s.vel, &s.mass, acc)?;
    for (v, a) in s.vel.iter_mut().zip(acc.iter()) {
        *v += *a * (0.5 * h);
    }
    Ok(())
}

/// Velocity Verlet (kick-drift-kick leapfrog). Second order, symplectic, time-reversible.
#[derive(Default)]
pub struct VelocityVerlet {
    acc: Vec<Vec3>,
}

impl Integrator for VelocityVerlet {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        kdk(s, forces, &mut self.acc, dt)
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
#[derive(Default)]
pub struct Yoshida4 {
    acc: Vec<Vec3>,
}

impl Integrator for Yoshida4 {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        let cbrt2 = 2f64.cbrt();
        let w1 = 1.0 / (2.0 - cbrt2);
        let w0 = -cbrt2 / (2.0 - cbrt2);
        for w in [w1, w0, w1] {
            kdk(s, forces, &mut self.acc, w * dt)?;
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
