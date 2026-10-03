//! Particle-mesh (PM) gravity: deposit the particles' mass on a grid (cloud-in-cell), solve
//! for the field with FFTs, and interpolate the accelerations back with the same weights.
//!
//! Two modes:
//!
//! - **isolated** (default): open boundaries. The mass grid is convolved with the softened
//!   Newtonian kernel on a zero-padded grid of twice the size (Hockney & Eastwood), so there
//!   are no periodic images. Forces match direct summation with the same softening to the grid
//!   accuracy at separations of a few cells, at a cost of O(N + M log M) for M grid points,
//!   independent of how the particles are arranged. Every particle must lie inside the box.
//! - **periodic**: the box repeats in every direction (cosmological boxes, plasmas). Poisson's
//!   equation `∇²φ = 4πG (ρ - ρ̄)` is solved with the discrete Laplacian's FFT eigenvalues and
//!   `a = -∇φ` by central differences. Positions are wrapped into the box.
//!
//! Using the same cloud-in-cell weights for deposit and interpolation, with an antisymmetric
//! force kernel (isolated) or a symmetric difference operator (periodic), means a particle
//! exerts no force on itself and total momentum is conserved to round-off.

use std::sync::{Arc, Mutex};

use super::fft::{fft_nd, fft_nd_pruned, Complex};
use crate::error::{invalid, Result};
use crate::forces::{BuiltinForce, Force, Param};
use crate::vec3::Vec3;

/// Particle-mesh gravity (see the module documentation).
pub struct ParticleMesh {
    pub g: f64,
    /// Grid points per axis (a power of two, at least 4).
    pub cells: usize,
    /// Side length of the cubic box.
    pub box_size: f64,
    /// Centre of the box.
    pub center: Vec3,
    pub periodic: bool,
    /// Plummer softening of the isolated kernel; `None` uses one grid spacing.
    pub softening: Option<f64>,
    /// Fourier transforms of the isolated kernels (potential, then x, y, z force), built on
    /// first use. Each kernel's transform is purely real or purely imaginary by symmetry, so
    /// only that part is kept.
    kernels: Mutex<Option<Arc<Vec<Vec<f64>>>>>,
}

impl Clone for ParticleMesh {
    fn clone(&self) -> Self {
        Self::new(
            self.g,
            self.cells,
            self.box_size,
            self.center,
            self.periodic,
            self.softening,
        )
    }
}

impl std::fmt::Debug for ParticleMesh {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ParticleMesh")
            .field("g", &self.g)
            .field("cells", &self.cells)
            .field("box_size", &self.box_size)
            .field("center", &self.center)
            .field("periodic", &self.periodic)
            .field("softening", &self.softening)
            .finish()
    }
}

/// Cloud-in-cell stencil of one particle: 8 (grid index, weight) pairs.
type Stencil = [(usize, f64); 8];

/// Accelerations at the grid nodes, and optionally the potential there.
type NodeFields = (Vec<[f64; 3]>, Option<Vec<f64>>);

impl ParticleMesh {
    pub fn new(
        g: f64,
        cells: usize,
        box_size: f64,
        center: Vec3,
        periodic: bool,
        softening: Option<f64>,
    ) -> Self {
        Self {
            g,
            cells,
            box_size,
            center,
            periodic,
            softening,
            kernels: Mutex::new(None),
        }
    }

    /// Grid spacing `box_size / cells`.
    pub fn spacing(&self) -> f64 {
        self.box_size / self.cells as f64
    }

    fn eps(&self) -> f64 {
        self.softening.unwrap_or_else(|| self.spacing())
    }

    pub fn validate(&self) -> Result<()> {
        if !(self.cells >= 4 && self.cells.is_power_of_two() && self.cells <= 512) {
            return invalid(format!(
                "cells must be a power of two between 4 and 512, got {}",
                self.cells
            ));
        }
        if !(self.box_size > 0.0 && self.box_size.is_finite()) {
            return invalid(format!("box_size must be positive, got {}", self.box_size));
        }
        if !self.g.is_finite() {
            return invalid("G must be finite");
        }
        if let Some(e) = self.softening {
            if !(e >= 0.0 && e.is_finite()) {
                return invalid(format!("softening must be non-negative, got {e}"));
            }
            if e == 0.0 && !self.periodic {
                return invalid("isolated particle-mesh gravity needs a positive softening");
            }
        }
        Ok(())
    }

    /// Cloud-in-cell weights of the particle at `x`.
    fn stencil(&self, x: Vec3, index: usize) -> Result<Stencil> {
        let n = self.cells;
        let h = self.spacing();
        let lo = self.center - Vec3::new(1.0, 1.0, 1.0) * (0.5 * self.box_size);
        let mut base = [0usize; 3];
        let mut frac = [0.0; 3];
        for a in 0..3 {
            let mut s = (x[a] - lo[a]) / h;
            if self.periodic {
                s = s.rem_euclid(n as f64);
                let i = (s.floor() as usize).min(n - 1);
                base[a] = i;
                frac[a] = s - i as f64;
            } else {
                // Nodes sit at lo + i h, i = 0..n-1; the particle needs two of them per axis.
                if !(s >= 0.0 && s <= (n - 1) as f64) {
                    return invalid(format!(
                        "particle {index} at {x:?} is outside the particle-mesh box \
                         (centre {:?}, side {}, nodes up to {} from the centre)",
                        self.center,
                        self.box_size,
                        0.5 * self.box_size - h
                    ));
                }
                let i = (s.floor() as usize).min(n - 2);
                base[a] = i;
                frac[a] = s - i as f64;
            }
        }
        let mut out = [(0usize, 0.0); 8];
        for (c, o) in out.iter_mut().enumerate() {
            let mut idx = 0;
            let mut w = 1.0;
            for a in 0..3 {
                let up = (c >> a) & 1;
                let i = (base[a] + up) % n;
                idx = idx * n + i;
                w *= if up == 1 { frac[a] } else { 1.0 - frac[a] };
            }
            *o = (idx, w);
        }
        Ok(out)
    }

    fn stencils(&self, pos: &[Vec3]) -> Result<Vec<Stencil>> {
        pos.iter()
            .enumerate()
            .map(|(i, &x)| self.stencil(x, i))
            .collect()
    }

    /// Kernel transforms on the padded grid of `2n` per axis.
    fn isolated_kernels(&self) -> Result<Arc<Vec<Vec<f64>>>> {
        let mut cache = self.kernels.lock().unwrap();
        if let Some(k) = cache.as_ref() {
            return Ok(Arc::clone(k));
        }
        let m = 2 * self.cells;
        let h = self.spacing();
        let eps2 = self.eps().powi(2);
        let len = m * m * m;
        let signed = |i: usize| -> f64 {
            if i <= self.cells {
                i as f64
            } else {
                i as f64 - m as f64
            }
        };
        let mut out = Vec::with_capacity(4);
        for kind in 0..4 {
            let mut data: Vec<Complex> = vec![(0.0, 0.0); len];
            for (idx, d) in data.iter_mut().enumerate() {
                let i = idx / (m * m);
                let j = (idx / m) % m;
                let k = idx % m;
                // The displacement exactly half-way across (index n) is ambiguous; both signs
                // give the same distance, and the odd force kernels vanish there.
                let r = Vec3::new(signed(i), signed(j), signed(k)) * h;
                let s2 = r.norm_squared() + eps2;
                d.0 = match kind {
                    0 => -1.0 / s2.sqrt(),
                    a => {
                        let comp = r[a - 1];
                        let on_edge = [i, j, k][a - 1] == self.cells;
                        if on_edge {
                            0.0
                        } else {
                            -comp / (s2 * s2.sqrt())
                        }
                    }
                };
            }
            fft_nd(&mut data, [m, m, m], false)?;
            // Even kernel: real transform. Odd kernels: imaginary transform.
            out.push(
                data.iter()
                    .map(|c| if kind == 0 { c.0 } else { c.1 })
                    .collect(),
            );
        }
        let out = Arc::new(out);
        *cache = Some(Arc::clone(&out));
        Ok(out)
    }

    /// Accelerations on the grid nodes (isolated mode), and optionally the potential.
    fn isolated_fields(
        &self,
        stencils: &[Stencil],
        mass: &[f64],
        want_potential: bool,
    ) -> Result<NodeFields> {
        let n = self.cells;
        let m = 2 * n;
        let kernels = self.isolated_kernels()?;
        // Mass on the padded grid.
        let mut rho: Vec<Complex> = vec![(0.0, 0.0); m * m * m];
        for (st, &mj) in stencils.iter().zip(mass) {
            for &(idx, w) in st {
                let (i, j, k) = (idx / (n * n), (idx / n) % n, idx % n);
                rho[(i * m + j) * m + k].0 += mj * w;
            }
        }
        // Only the first octant holds mass, and only the first octant of the result is used.
        fft_nd_pruned(&mut rho, [m, m, m], false, [n, n, n])?;
        let field = |kind: usize| -> Result<Vec<f64>> {
            let kern = &kernels[kind];
            let mut prod: Vec<Complex> = rho
                .iter()
                .zip(kern)
                .map(|(&(a, b), &k)| {
                    if kind == 0 {
                        (a * k, b * k) // × real
                    } else {
                        (-b * k, a * k) // × i k
                    }
                })
                .collect();
            fft_nd_pruned(&mut prod, [m, m, m], true, [n, n, n])?;
            let mut out = vec![0.0; n * n * n];
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        out[(i * n + j) * n + k] = self.g * prod[(i * m + j) * m + k].0;
                    }
                }
            }
            Ok(out)
        };
        // The three components are independent transforms.
        #[cfg(feature = "parallel")]
        let (ax, (ay, az)) = rayon::join(|| field(1), || rayon::join(|| field(2), || field(3)));
        #[cfg(not(feature = "parallel"))]
        let (ax, ay, az) = (field(1), field(2), field(3));
        let (ax, ay, az) = (ax?, ay?, az?);
        let acc = (0..n * n * n).map(|i| [ax[i], ay[i], az[i]]).collect();
        let pot = if want_potential {
            Some(field(0)?)
        } else {
            None
        };
        Ok((acc, pot))
    }

    /// Accelerations on the grid nodes (periodic mode).
    fn periodic_fields(&self, stencils: &[Stencil], mass: &[f64]) -> Result<Vec<[f64; 3]>> {
        let n = self.cells;
        let h = self.spacing();
        let cell_volume = h * h * h;
        let mut rho: Vec<Complex> = vec![(0.0, 0.0); n * n * n];
        for (st, &mj) in stencils.iter().zip(mass) {
            for &(idx, w) in st {
                rho[idx].0 += mj * w / cell_volume;
            }
        }
        fft_nd(&mut rho, [n, n, n], false)?;
        let two_pi = 2.0 * std::f64::consts::PI;
        for (idx, d) in rho.iter_mut().enumerate() {
            let ks = [idx / (n * n), (idx / n) % n, idx % n];
            let lam: f64 = ks
                .iter()
                .map(|&k| (2.0 * (two_pi * k as f64 / n as f64).cos() - 2.0) / (h * h))
                .sum();
            let scale = if lam == 0.0 {
                0.0 // the mean density: Jeans swindle
            } else {
                4.0 * std::f64::consts::PI * self.g / lam
            };
            *d = (d.0 * scale, d.1 * scale);
        }
        fft_nd(&mut rho, [n, n, n], true)?;
        let phi: Vec<f64> = rho.iter().map(|c| c.0).collect();
        let at = |i: usize, j: usize, k: usize| phi[((i % n) * n + j % n) * n + k % n];
        let mut acc = vec![[0.0; 3]; n * n * n];
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    acc[(i * n + j) * n + k] = [
                        -(at(i + 1, j, k) - at(i + n - 1, j, k)) / (2.0 * h),
                        -(at(i, j + 1, k) - at(i, j + n - 1, k)) / (2.0 * h),
                        -(at(i, j, k + 1) - at(i, j, k + n - 1)) / (2.0 * h),
                    ];
                }
            }
        }
        Ok(acc)
    }
}

impl Force for ParticleMesh {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        self.validate()?;
        let stencils = self.stencils(pos)?;
        let grid_acc = if self.periodic {
            self.periodic_fields(&stencils, mass)?
        } else {
            self.isolated_fields(&stencils, mass, false)?.0
        };
        for (a, st) in acc.iter_mut().zip(&stencils) {
            let mut sum = [0.0; 3];
            for &(idx, w) in st {
                for (s, g) in sum.iter_mut().zip(grid_acc[idx]) {
                    *s += w * g;
                }
            }
            *a += Vec3::from(sum);
        }
        Ok(())
    }

    /// Isolated mode: `½ Σ_j m_j φ(x_j)` with each particle's own (grid) self-energy removed,
    /// an approximation to the pairwise softened energy at the grid accuracy. Periodic mode
    /// reports no potential.
    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        if self.periodic {
            return Ok(None);
        }
        self.validate()?;
        let stencils = self.stencils(pos)?;
        let (_, phi) = self.isolated_fields(&stencils, mass, true)?;
        let phi = phi.expect("potential requested");
        let n = self.cells;
        let h = self.spacing();
        let eps2 = self.eps().powi(2);
        let mut u = 0.0;
        for (st, &mj) in stencils.iter().zip(mass) {
            let mut phi_j = 0.0;
            let mut self_phi = 0.0;
            for &(a, wa) in st {
                phi_j += wa * phi[a];
                for &(b, wb) in st {
                    let d = |x: usize, y: usize| {
                        let (xi, xj, xk) = (x / (n * n), (x / n) % n, x % n);
                        let (yi, yj, yk) = (y / (n * n), (y / n) % n, y % n);
                        Vec3::new(
                            xi as f64 - yi as f64,
                            xj as f64 - yj as f64,
                            xk as f64 - yk as f64,
                        ) * h
                    };
                    let s2 = d(a, b).norm_squared() + eps2;
                    self_phi += wa * wb * mj * (-self.g / s2.sqrt());
                }
            }
            u += 0.5 * mj * (phi_j - self_phi);
        }
        Ok(Some(u))
    }

    fn name(&self) -> String {
        "ParticleMesh".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::ParticleMesh(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![("G", Param::Scalar(self.g))]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match (name, value) {
            ("G", Param::Scalar(v)) if v.is_finite() => {
                self.g = v;
                Ok(())
            }
            _ => invalid(format!(
                "ParticleMesh has a settable scalar parameter G only, got {name}"
            )),
        }
    }
}
