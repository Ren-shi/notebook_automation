//! Gauss-Legendre implicit Runge-Kutta methods: symplectic and symmetric for any
//! Hamiltonian system, including velocity-dependent forces; they also conserve every
//! quadratic invariant exactly (e.g. kinetic energy in a pure magnetic field).

use super::Integrator;
use crate::error::{invalid, Result};
use crate::forces::ForceSet;
use crate::state::State;
use crate::vec3::Vec3;

/// `s`-stage Gauss-Legendre collocation (order `2s`), with stages solved by fixed-point
/// iteration to round-off. Each step costs `s` force evaluations per iteration (typically
/// 5-20 iterations for non-stiff problems); fails if the iteration does not converge
/// (the step is too large).
pub struct GaussLegendre {
    stages: usize,
    a: Vec<Vec<f64>>,
    b: Vec<f64>,
    c: Vec<f64>,
    name: &'static str,
    /// Slopes of each stage: velocity and acceleration.
    kx: Vec<Vec<Vec3>>,
    kv: Vec<Vec<Vec3>>,
    xs: Vec<Vec3>,
    vs: Vec<Vec3>,
    acc: Vec<Vec3>,
}

impl GaussLegendre {
    fn with(name: &'static str, a: Vec<Vec<f64>>, b: Vec<f64>, c: Vec<f64>) -> Self {
        let stages = b.len();
        Self {
            stages,
            a,
            b,
            c,
            name,
            kx: vec![Vec::new(); stages],
            kv: vec![Vec::new(); stages],
            xs: Vec::new(),
            vs: Vec::new(),
            acc: Vec::new(),
        }
    }

    /// Implicit midpoint rule: order 2, one stage.
    pub fn order2() -> Self {
        Self::with("gauss2", vec![vec![0.5]], vec![1.0], vec![0.5])
    }

    /// Two-stage Gauss-Legendre: order 4.
    pub fn order4() -> Self {
        let r = 3f64.sqrt() / 6.0;
        Self::with(
            "gauss4",
            vec![vec![0.25, 0.25 - r], vec![0.25 + r, 0.25]],
            vec![0.5, 0.5],
            vec![0.5 - r, 0.5 + r],
        )
    }

    /// Three-stage Gauss-Legendre: order 6.
    pub fn order6() -> Self {
        let r = 15f64.sqrt();
        Self::with(
            "gauss6",
            vec![
                vec![5.0 / 36.0, 2.0 / 9.0 - r / 15.0, 5.0 / 36.0 - r / 30.0],
                vec![5.0 / 36.0 + r / 24.0, 2.0 / 9.0, 5.0 / 36.0 - r / 24.0],
                vec![5.0 / 36.0 + r / 30.0, 2.0 / 9.0 + r / 15.0, 5.0 / 36.0],
            ],
            vec![5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0],
            vec![0.5 - r / 10.0, 0.5, 0.5 + r / 10.0],
        )
    }
}

const MAX_ITERATIONS: usize = 100;

impl Integrator for GaussLegendre {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        let n = s.len();
        // Initial guess: every stage slope equals the slope at the start.
        forces.accelerations(
            s.t,
            &s.pos,
            &s.vel,
            &s.mass,
            &s.charge,
            &s.pinned,
            &mut self.acc,
        )?;
        for i in 0..self.stages {
            self.kx[i].clone_from(&s.vel);
            self.kv[i].clone_from(&self.acc);
        }
        let scale = |x: &[Vec3]| x.iter().map(|v| v.norm()).fold(0.0, f64::max);
        let mut last_change = f64::INFINITY;
        let mut converged = false;
        for _ in 0..MAX_ITERATIONS {
            let mut change: f64 = 0.0;
            for i in 0..self.stages {
                self.xs.clone_from(&s.pos);
                self.vs.clone_from(&s.vel);
                for j in 0..self.stages {
                    let w = self.a[i][j] * dt;
                    for p in 0..n {
                        self.xs[p] += self.kx[j][p] * w;
                        self.vs[p] += self.kv[j][p] * w;
                    }
                }
                forces.accelerations(
                    s.t + self.c[i] * dt,
                    &self.xs,
                    &self.vs,
                    &s.mass,
                    &s.charge,
                    &s.pinned,
                    &mut self.acc,
                )?;
                for p in 0..n {
                    change = change
                        .max((self.vs[p] - self.kx[i][p]).norm())
                        .max((self.acc[p] - self.kv[i][p]).norm() * dt.abs());
                }
                self.kx[i].clone_from(&self.vs);
                self.kv[i].clone_from(&self.acc);
            }
            let size = scale(&s.vel)
                .max(scale(&self.acc) * dt.abs())
                .max(f64::MIN_POSITIVE);
            // Iterate to round-off: stop once the update vanishes or stops shrinking.
            if change <= 4.0 * f64::EPSILON * size
                || (change >= last_change && change < 1e-10 * size)
            {
                converged = true;
                break;
            }
            last_change = change;
        }
        if !converged {
            return invalid(format!(
                "{}: the implicit stage equations did not converge at t = {} (time step too large?)",
                self.name, s.t
            ));
        }
        for p in 0..n {
            let (mut dx, mut dv) = (Vec3::ZERO, Vec3::ZERO);
            for i in 0..self.stages {
                dx += self.kx[i][p] * self.b[i];
                dv += self.kv[i][p] * self.b[i];
            }
            s.pos[p] += dx * dt;
            s.vel[p] += dv * dt;
        }
        s.t += dt;
        Ok(())
    }
    fn name(&self) -> &str {
        self.name
    }
    fn order(&self) -> u32 {
        2 * self.stages as u32
    }
    fn symplectic(&self) -> bool {
        true
    }
}
