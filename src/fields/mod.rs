//! Fields on regular grids: finite differences in 1, 2 and 3 dimensions, boundary
//! conditions, and solvers for the wave, heat and Poisson equations ([`solvers`]), plus
//! particle-mesh gravity coupling particles to a grid ([`pm`]).
//!
//! A [`Grid`] has `n` points per axis (unused axes have `n = 1`) spaced `h` apart, stored
//! row-major: point `(i, j, k)` is at index `(i * ny + j) * nz + k`. Boundary conditions:
//!
//! - [`Boundary::Dirichlet`]: the first and last point of each axis are boundary nodes whose
//!   values are held fixed (whatever the field holds there).
//! - [`Boundary::Neumann`]: zero normal derivative. Points are cell centres and the boundary is
//!   half a spacing outside the first and last point (the ghost value equals its neighbour),
//!   which keeps the discrete Laplacian symmetric.
//! - [`Boundary::Periodic`]: the grid wraps around; point `n` is point `0`.
//!
//! The Laplacian is the standard second-order `(u[i+1] - 2 u[i] + u[i-1]) / h²` stencil
//! summed over the active axes.

pub mod fft;
pub mod pm;
pub mod solvers;

use crate::error::{invalid, Result};

/// Boundary condition of a [`Grid`] (the same on every face).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    Dirichlet,
    Neumann,
    Periodic,
}

impl Boundary {
    pub fn parse(name: &str) -> Result<Self> {
        match name {
            "dirichlet" => Ok(Boundary::Dirichlet),
            "neumann" => Ok(Boundary::Neumann),
            "periodic" => Ok(Boundary::Periodic),
            _ => invalid(format!(
                "boundary must be 'dirichlet', 'neumann' or 'periodic', got {name:?}"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Boundary::Dirichlet => "dirichlet",
            Boundary::Neumann => "neumann",
            Boundary::Periodic => "periodic",
        }
    }
}

/// A regular grid of `n[0] × n[1] × n[2]` points spaced `h` apart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid {
    pub n: [usize; 3],
    pub h: f64,
    pub boundary: Boundary,
}

impl Grid {
    /// A grid with `shape` (1 to 3 sizes) and spacing `h`.
    pub fn new(shape: &[usize], h: f64, boundary: Boundary) -> Result<Self> {
        if shape.is_empty() || shape.len() > 3 {
            return invalid(format!("a grid has 1 to 3 axes, got {}", shape.len()));
        }
        let min = if boundary == Boundary::Dirichlet {
            3
        } else {
            2
        };
        if shape.iter().any(|&s| s < min) {
            return invalid(format!(
                "every axis needs at least {min} points with {} boundaries, got {shape:?}",
                boundary.name()
            ));
        }
        if !(h > 0.0 && h.is_finite()) {
            return invalid(format!("grid spacing must be positive, got {h}"));
        }
        let mut n = [1; 3];
        n[..shape.len()].copy_from_slice(shape);
        Ok(Self { n, h, boundary })
    }

    /// Number of axes (1, 2 or 3).
    pub fn dims(&self) -> usize {
        if self.n[2] > 1 {
            3
        } else if self.n[1] > 1 {
            2
        } else {
            1
        }
    }

    pub fn shape(&self) -> Vec<usize> {
        self.n[..self.dims()].to_vec()
    }

    pub fn len(&self) -> usize {
        self.n[0] * self.n[1] * self.n[2]
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn strides(&self) -> [usize; 3] {
        [self.n[1] * self.n[2], self.n[2], 1]
    }

    pub(crate) fn check(&self, u: &[f64], what: &str) -> Result<()> {
        if u.len() != self.len() {
            return invalid(format!(
                "{what}: expected {} values for a grid of shape {:?}, got {}",
                self.len(),
                self.shape(),
                u.len()
            ));
        }
        Ok(())
    }

    /// Whether flat index `idx` is a (fixed) Dirichlet boundary node.
    pub fn is_fixed(&self, idx: usize) -> bool {
        if self.boundary != Boundary::Dirichlet {
            return false;
        }
        let s = self.strides();
        (0..self.dims()).any(|a| {
            let i = (idx / s[a]) % self.n[a];
            i == 0 || i == self.n[a] - 1
        })
    }

    /// `out = ∇²u` (zero at Dirichlet boundary nodes).
    pub fn laplacian(&self, u: &[f64], out: &mut [f64]) {
        let s = self.strides();
        let inv_h2 = 1.0 / (self.h * self.h);
        let dims = self.dims();
        for (idx, o) in out.iter_mut().enumerate() {
            if self.is_fixed(idx) {
                *o = 0.0;
                continue;
            }
            let mut acc = 0.0;
            for (&n, &stride) in self.n.iter().zip(&s).take(dims) {
                let i = (idx / stride) % n;
                let up = if i + 1 < n {
                    idx + stride
                } else {
                    match self.boundary {
                        Boundary::Periodic => idx + stride - n * stride,
                        _ => idx, // Neumann ghost equals the edge value
                    }
                };
                let down = if i > 0 {
                    idx - stride
                } else {
                    match self.boundary {
                        Boundary::Periodic => idx + (n - 1) * stride,
                        _ => idx,
                    }
                };
                acc += u[up] - 2.0 * u[idx] + u[down];
            }
            *o = acc * inv_h2;
        }
    }

    /// Coordinates of point `idx` (Neumann grids are cell-centred: `(i + ½) h`).
    pub fn coordinates(&self, idx: usize) -> [f64; 3] {
        let s = self.strides();
        let shift = if self.boundary == Boundary::Neumann {
            0.5
        } else {
            0.0
        };
        let mut x = [0.0; 3];
        for a in 0..self.dims() {
            x[a] = (((idx / s[a]) % self.n[a]) as f64 + shift) * self.h;
        }
        x
    }
}

/// Solves `(alpha I - beta ∇²) x = b` (`alpha, beta ≥ 0`) by conjugate gradients, holding
/// Dirichlet boundary nodes at their values in `x` (the initial guess). Returns the
/// iterations used. With `alpha = 0` and Neumann or periodic boundaries the operator is
/// singular; `b` must then have zero mean.
pub(crate) fn solve_shifted(
    grid: &Grid,
    alpha: f64,
    beta: f64,
    b: &[f64],
    x: &mut [f64],
    tol: f64,
) -> Result<usize> {
    let n = grid.len();
    let mut lap = vec![0.0; n];
    // Operator on vectors that vanish at fixed nodes (symmetric on the free nodes).
    let apply = |p: &[f64], out: &mut [f64], lap: &mut [f64]| {
        grid.laplacian(p, lap);
        for i in 0..n {
            out[i] = if grid.is_fixed(i) {
                0.0
            } else {
                alpha * p[i] - beta * lap[i]
            };
        }
    };
    // r = b - A x on the free nodes (A x uses the boundary values in x).
    let mut ax = vec![0.0; n];
    grid.laplacian(x, &mut lap);
    for i in 0..n {
        ax[i] = alpha * x[i] - beta * lap[i];
    }
    let mut r: Vec<f64> = (0..n)
        .map(|i| if grid.is_fixed(i) { 0.0 } else { b[i] - ax[i] })
        .collect();
    let norm_b = b
        .iter()
        .enumerate()
        .filter(|(i, _)| !grid.is_fixed(*i))
        .map(|(_, v)| v * v)
        .sum::<f64>()
        .sqrt()
        .max(f64::MIN_POSITIVE);
    let mut p = r.clone();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    let mut ap = vec![0.0; n];
    let max_iter = 10 * n + 100;
    for it in 0..max_iter {
        if rr.sqrt() <= tol * norm_b {
            return Ok(it);
        }
        apply(&p, &mut ap, &mut lap);
        let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
        if pap <= 0.0 {
            return invalid(
                "conjugate gradients broke down (the operator is not positive definite)",
            );
        }
        let step = rr / pap;
        for i in 0..n {
            x[i] += step * p[i];
            r[i] -= step * ap[i];
        }
        let rr_new: f64 = r.iter().map(|v| v * v).sum();
        let mix = rr_new / rr;
        rr = rr_new;
        for i in 0..n {
            p[i] = r[i] + mix * p[i];
        }
    }
    invalid(format!(
        "conjugate gradients did not converge in {max_iter} iterations"
    ))
}
