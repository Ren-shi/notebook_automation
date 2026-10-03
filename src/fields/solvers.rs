//! Time-dependent and elliptic solvers on a [`Grid`].
//!
//! - [`Wave`]: `u_tt = c² ∇²u` by leapfrog (kick-drift-kick), second order in space and time,
//!   stable for `c dt / h ≤ 1/√d` in `d` dimensions; it conserves a discrete energy to O(dt²)
//!   without drift.
//! - [`Heat`]: `u_t = D ∇²u`, explicitly (forward Euler, first order in time, stable for
//!   `D dt / h² ≤ 1/(2d)`) or by Crank-Nicolson (second order in time, unconditionally
//!   stable; the implicit system is solved by conjugate gradients).
//! - [`poisson`]: `∇²φ = f`, by FFT on periodic grids (exact for the discrete Laplacian) and by
//!   conjugate gradients otherwise.

use super::fft::{fft_nd, Complex};
use super::{solve_shifted, Boundary, Grid};
use crate::error::{invalid, Result};

/// Tolerance of the implicit solves (relative residual).
const CG_TOL: f64 = 1e-12;

/// The wave equation `u_tt = c² ∇²u`.
#[derive(Debug, Clone)]
pub struct Wave {
    pub grid: Grid,
    pub c: f64,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub t: f64,
    lap: Vec<f64>,
    lap_valid: bool,
}

impl Wave {
    /// A field at rest (`u = v = 0`) with wave speed `c`.
    pub fn new(grid: Grid, c: f64) -> Result<Self> {
        if !(c > 0.0 && c.is_finite()) {
            return invalid(format!("wave speed must be positive, got {c}"));
        }
        let n = grid.len();
        Ok(Self {
            grid,
            c,
            u: vec![0.0; n],
            v: vec![0.0; n],
            t: 0.0,
            lap: vec![0.0; n],
            lap_valid: false,
        })
    }

    /// Sets the displacement and velocity fields.
    pub fn set_state(&mut self, u: Vec<f64>, v: Vec<f64>) -> Result<()> {
        self.grid.check(&u, "u")?;
        self.grid.check(&v, "v")?;
        self.u = u;
        self.v = v;
        self.lap_valid = false;
        Ok(())
    }

    /// Largest stable step, `h / (c √d)`.
    pub fn max_stable_dt(&self) -> f64 {
        self.grid.h / (self.c * (self.grid.dims() as f64).sqrt())
    }

    /// Takes `n` leapfrog steps of `dt`.
    pub fn step(&mut self, dt: f64, n: usize) -> Result<()> {
        if !(dt > 0.0 && dt <= self.max_stable_dt() * (1.0 + 1e-12)) {
            return invalid(format!(
                "wave step must be in (0, h / (c √d)] = (0, {:.6e}] for stability, got {dt}",
                self.max_stable_dt()
            ));
        }
        let c2 = self.c * self.c;
        for _ in 0..n {
            if !self.lap_valid {
                self.grid.laplacian(&self.u, &mut self.lap);
            }
            for (v, l) in self.v.iter_mut().zip(&self.lap) {
                *v += 0.5 * dt * c2 * l;
            }
            for i in 0..self.u.len() {
                if !self.grid.is_fixed(i) {
                    self.u[i] += dt * self.v[i];
                }
            }
            self.grid.laplacian(&self.u, &mut self.lap);
            self.lap_valid = true;
            for (v, l) in self.v.iter_mut().zip(&self.lap) {
                *v += 0.5 * dt * c2 * l;
            }
            self.t += dt;
        }
        Ok(())
    }

    /// Discrete energy `½ Σ v² + ½ c² Σ_edges ((u_a - u_b) / h)²`, times the cell volume `h^d`.
    pub fn energy(&self) -> f64 {
        let g = &self.grid;
        let kinetic: f64 = self
            .v
            .iter()
            .enumerate()
            .filter(|(i, _)| !g.is_fixed(*i))
            .map(|(_, v)| v * v)
            .sum();
        let strides = [g.n[1] * g.n[2], g.n[2], 1];
        let mut strain = 0.0;
        for idx in 0..g.len() {
            for (&n, &stride) in g.n.iter().zip(&strides).take(g.dims()) {
                let i = (idx / stride) % n;
                let next = if i + 1 < n {
                    idx + stride
                } else if g.boundary == Boundary::Periodic {
                    idx + stride - n * stride
                } else {
                    continue;
                };
                let d = (self.u[next] - self.u[idx]) / g.h;
                strain += d * d;
            }
        }
        0.5 * (kinetic + self.c * self.c * strain) * g.h.powi(g.dims() as i32)
    }
}

/// Time integration scheme of [`Heat`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeatMethod {
    Explicit,
    CrankNicolson,
}

/// The heat (diffusion) equation `u_t = D ∇²u`.
#[derive(Debug, Clone)]
pub struct Heat {
    pub grid: Grid,
    pub diffusivity: f64,
    pub method: HeatMethod,
    pub u: Vec<f64>,
    pub t: f64,
    /// Conjugate-gradient iterations used by the last Crank-Nicolson step.
    pub last_iterations: usize,
}

impl Heat {
    pub fn new(grid: Grid, diffusivity: f64, method: HeatMethod) -> Result<Self> {
        if !(diffusivity > 0.0 && diffusivity.is_finite()) {
            return invalid(format!("diffusivity must be positive, got {diffusivity}"));
        }
        Ok(Self {
            grid,
            diffusivity,
            method,
            u: vec![0.0; grid.len()],
            t: 0.0,
            last_iterations: 0,
        })
    }

    pub fn set_u(&mut self, u: Vec<f64>) -> Result<()> {
        self.grid.check(&u, "u")?;
        self.u = u;
        Ok(())
    }

    /// Largest stable explicit step, `h² / (2 d D)`.
    pub fn max_explicit_dt(&self) -> f64 {
        self.grid.h * self.grid.h / (2.0 * self.grid.dims() as f64 * self.diffusivity)
    }

    /// Takes `n` steps of `dt`.
    pub fn step(&mut self, dt: f64, n: usize) -> Result<()> {
        if !(dt > 0.0 && dt.is_finite()) {
            return invalid(format!("dt must be positive, got {dt}"));
        }
        let len = self.u.len();
        let mut lap = vec![0.0; len];
        match self.method {
            HeatMethod::Explicit => {
                if dt > self.max_explicit_dt() * (1.0 + 1e-12) {
                    return invalid(format!(
                        "explicit heat step must be at most h² / (2 d D) = {:.6e} for stability, got {dt} \
                         (use Crank-Nicolson for larger steps)",
                        self.max_explicit_dt()
                    ));
                }
                for _ in 0..n {
                    self.grid.laplacian(&self.u, &mut lap);
                    for (u, l) in self.u.iter_mut().zip(&lap) {
                        *u += dt * self.diffusivity * l;
                    }
                    self.t += dt;
                }
            }
            HeatMethod::CrankNicolson => {
                let a = 0.5 * dt * self.diffusivity;
                let mut rhs = vec![0.0; len];
                for _ in 0..n {
                    self.grid.laplacian(&self.u, &mut lap);
                    for i in 0..len {
                        rhs[i] = self.u[i] + a * lap[i];
                    }
                    self.last_iterations =
                        solve_shifted(&self.grid, 1.0, a, &rhs, &mut self.u, CG_TOL)?;
                    self.t += dt;
                }
            }
        }
        Ok(())
    }

    /// `Σ u h^d`: conserved with Neumann and periodic boundaries.
    pub fn total(&self) -> f64 {
        self.u.iter().sum::<f64>() * self.grid.h.powi(self.grid.dims() as i32)
    }
}

/// Solves `∇²φ = f` on `grid`, writing `φ` into `phi`.
///
/// - Periodic grids (every size a power of two): by FFT, with the discrete Laplacian's
///   eigenvalues, so the result solves the finite-difference equations exactly. The mean of
///   `f` is removed (it has no periodic solution) and `φ` has zero mean.
/// - Dirichlet grids: by conjugate gradients; the boundary nodes of `phi` give the boundary
///   values (set them before calling, e.g. to zero).
/// - Neumann grids: by conjugate gradients after removing the mean of `f`; `φ` has zero mean.
///
/// Returns the number of iterations (0 for the FFT).
pub fn poisson(grid: &Grid, f: &[f64], phi: &mut [f64]) -> Result<usize> {
    grid.check(f, "f")?;
    grid.check(phi, "phi")?;
    let mean = f.iter().sum::<f64>() / f.len() as f64;
    match grid.boundary {
        Boundary::Periodic => {
            let mut data: Vec<Complex> = f.iter().map(|&v| (v - mean, 0.0)).collect();
            fft_nd(&mut data, grid.n, false)?;
            let h2 = grid.h * grid.h;
            let strides = [grid.n[1] * grid.n[2], grid.n[2], 1];
            for (idx, d) in data.iter_mut().enumerate() {
                let mut lam = 0.0;
                for (&n, &stride) in grid.n.iter().zip(&strides).take(grid.dims()) {
                    let k = (idx / stride) % n;
                    let theta = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
                    lam += (2.0 * theta.cos() - 2.0) / h2;
                }
                if lam == 0.0 {
                    *d = (0.0, 0.0);
                } else {
                    *d = (d.0 / lam, d.1 / lam);
                }
            }
            fft_nd(&mut data, grid.n, true)?;
            for (p, d) in phi.iter_mut().zip(&data) {
                *p = d.0;
            }
            Ok(0)
        }
        Boundary::Dirichlet => {
            let b: Vec<f64> = f.iter().map(|v| -v).collect();
            solve_shifted(grid, 0.0, 1.0, &b, phi, CG_TOL)
        }
        Boundary::Neumann => {
            let b: Vec<f64> = f.iter().map(|v| mean - v).collect();
            phi.iter_mut().for_each(|p| *p = 0.0);
            let it = solve_shifted(grid, 0.0, 1.0, &b, phi, CG_TOL)?;
            let m = phi.iter().sum::<f64>() / phi.len() as f64;
            phi.iter_mut().for_each(|p| *p -= m);
            Ok(it)
        }
    }
}
