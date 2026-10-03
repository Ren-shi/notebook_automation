//! The time-dependent Schrödinger equation for one particle on a [`Grid`]:
//!
//! `iħ ∂ψ/∂t = Hψ`, with `H = -ħ²/(2m) ∇² + V(x, t) - iW(x)`,
//!
//! where `V` is the potential and `W ≥ 0` an optional absorbing (imaginary) potential that
//! removes outgoing waves near the edges. The wavefunction is normalised so that
//! `Σ |ψ|² h^d = 1`.
//!
//! - [`QuantumMethod::SplitStep`] (periodic grids, every size a power of two): Strang splitting
//!   `e^{-iVτ/2ħ} e^{-iTτ/ħ} e^{-iVτ/2ħ}` with the kinetic factor applied exactly in Fourier
//!   space, so space is resolved spectrally; second order in time, or fourth order with
//!   `order = 4` (Yoshida's composition of three Strang steps). Exactly unitary without an
//!   absorber.
//! - [`QuantumMethod::CrankNicolson`] (any boundary): `(1 + iτH/2ħ) ψ' = (1 - iτH/2ħ) ψ` with
//!   the second-order Laplacian and `H` at the midpoint time; second order in space and time
//!   and unitary, so the norm and (for a static potential) the discrete energy are conserved
//!   up to the linear-solver tolerance. One-dimensional systems are solved directly (Thomas
//!   algorithm, Sherman-Morrison on periodic grids), larger ones by Jacobi-preconditioned
//!   BiCGSTAB. On Dirichlet grids ψ vanishes on the edge nodes (hard walls).
//!
//! [`Schrodinger::eigenstates`] finds the lowest bound states of the same discretisation
//! (spectral for split-step, finite differences for Crank-Nicolson) with LOBPCG (locally
//! optimal block preconditioned conjugate gradients).

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

use super::fft::{fft_nd, Complex};
use super::{Boundary, Grid};
use crate::error::{invalid, Result};
use crate::rng::normal;

/// Relative residual of the iterative linear solves.
const SOLVE_TOL: f64 = 1e-12;

/// A time-dependent potential: fills `v` (one value per grid point) at time `t`.
pub type PotentialFn = Box<dyn FnMut(f64, &mut [f64]) -> Result<()> + Send + Sync>;

/// Time integration scheme of [`Schrodinger`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuantumMethod {
    SplitStep,
    CrankNicolson,
}

impl QuantumMethod {
    pub fn parse(name: &str) -> Result<Self> {
        match name {
            "split_step" => Ok(QuantumMethod::SplitStep),
            "crank_nicolson" => Ok(QuantumMethod::CrankNicolson),
            _ => invalid(format!(
                "method must be 'split_step' or 'crank_nicolson', got {name:?}"
            )),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            QuantumMethod::SplitStep => "split_step",
            QuantumMethod::CrankNicolson => "crank_nicolson",
        }
    }
}

/// A complex number with arithmetic, for the solvers (the FFT works on [`Complex`] tuples).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct C {
    re: f64,
    im: f64,
}

impl C {
    const ZERO: C = C { re: 0.0, im: 0.0 };

    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// `e^{iθ}`.
    fn cis(theta: f64) -> Self {
        let (s, c) = theta.sin_cos();
        Self { re: c, im: s }
    }

    fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    fn norm2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

impl From<Complex> for C {
    fn from((re, im): Complex) -> Self {
        Self { re, im }
    }
}

impl From<C> for Complex {
    fn from(c: C) -> Self {
        (c.re, c.im)
    }
}

impl Add for C {
    type Output = C;
    fn add(self, o: C) -> C {
        C::new(self.re + o.re, self.im + o.im)
    }
}

impl Sub for C {
    type Output = C;
    fn sub(self, o: C) -> C {
        C::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for C {
    type Output = C;
    fn mul(self, o: C) -> C {
        C::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

impl Mul<f64> for C {
    type Output = C;
    fn mul(self, s: f64) -> C {
        C::new(self.re * s, self.im * s)
    }
}

impl Div for C {
    type Output = C;
    fn div(self, o: C) -> C {
        let d = o.norm2();
        let n = self * o.conj();
        C::new(n.re / d, n.im / d)
    }
}

impl Neg for C {
    type Output = C;
    fn neg(self) -> C {
        C::new(-self.re, -self.im)
    }
}

impl AddAssign for C {
    fn add_assign(&mut self, o: C) {
        *self = *self + o;
    }
}

impl SubAssign for C {
    fn sub_assign(&mut self, o: C) {
        *self = *self - o;
    }
}

/// `Σ conj(a_i) b_i`.
fn dot(a: &[C], b: &[C]) -> C {
    a.iter()
        .zip(b)
        .fold(C::ZERO, |acc, (x, y)| acc + x.conj() * *y)
}

fn norm(a: &[C]) -> f64 {
    a.iter().map(|v| v.norm2()).sum::<f64>().sqrt()
}

/// Yoshida's fourth-order composition: Strang steps of `w1 τ`, `w0 τ`, `w1 τ`.
fn yoshida_weights() -> [f64; 3] {
    let cbrt2 = 2f64.cbrt();
    let w1 = 1.0 / (2.0 - cbrt2);
    [w1, -cbrt2 * w1, w1]
}

/// The Schrödinger equation for one particle on a grid. See the [module docs](self).
pub struct Schrodinger {
    pub grid: Grid,
    pub hbar: f64,
    pub mass: f64,
    pub method: QuantumMethod,
    /// Coordinates of grid point `(0, 0, 0)`.
    pub origin: [f64; 3],
    /// The wavefunction at the grid points (row-major, like the grid).
    pub psi: Vec<Complex>,
    pub t: f64,
    /// Linear-solver iterations used by the last Crank-Nicolson step (0 for the direct 1D
    /// solve).
    pub last_iterations: usize,
    order: usize,
    potential: Vec<f64>,
    potential_fn: Option<PotentialFn>,
    /// Time at which `potential_fn` last filled `potential` (NaN: never).
    potential_t: f64,
    absorber: Vec<f64>,
    /// `false` at Dirichlet boundary nodes.
    free: Vec<bool>,
    /// `|k|²` at each point of the Fourier grid (split-step only).
    k2: Vec<f64>,
    /// `e^{-iħk²τ/2m}` per step `τ`.
    kinetic_cache: Vec<(f64, Vec<C>)>,
    /// `e^{-(iV + W)τ/ħ}` per step `τ`, for static potentials.
    kick_cache: Vec<(f64, Vec<C>)>,
    /// ψ at the start of the current step, restored if the step fails.
    backup: Vec<Complex>,
}

impl Schrodinger {
    /// ψ = 0 in zero potential. Split-step needs a periodic grid with power-of-two sizes.
    pub fn new(grid: Grid, hbar: f64, mass: f64, method: QuantumMethod) -> Result<Self> {
        if !(hbar > 0.0 && hbar.is_finite()) {
            return invalid(format!("hbar must be positive, got {hbar}"));
        }
        if !(mass > 0.0 && mass.is_finite()) {
            return invalid(format!("mass must be positive, got {mass}"));
        }
        if method == QuantumMethod::SplitStep {
            if grid.boundary != Boundary::Periodic {
                return invalid(format!(
                    "split-step needs a periodic grid, got {} boundaries (use crank_nicolson)",
                    grid.boundary.name()
                ));
            }
            if grid.n.iter().any(|&s| s > 1 && !s.is_power_of_two()) {
                return invalid(format!(
                    "split-step needs grid sizes that are powers of two, got {:?}",
                    grid.shape()
                ));
            }
        }
        let len = grid.len();
        let free = (0..len).map(|i| !grid.is_fixed(i)).collect();
        let k2 = if method == QuantumMethod::SplitStep {
            let strides = [grid.n[1] * grid.n[2], grid.n[2], 1];
            (0..len)
                .map(|idx| {
                    let mut k2 = 0.0;
                    for (&n, &stride) in grid.n.iter().zip(&strides).take(grid.dims()) {
                        let j = (idx / stride) % n;
                        let j = if j < n / 2 {
                            j as f64
                        } else {
                            j as f64 - n as f64
                        };
                        let k = 2.0 * std::f64::consts::PI * j / (n as f64 * grid.h);
                        k2 += k * k;
                    }
                    k2
                })
                .collect()
        } else {
            Vec::new()
        };
        Ok(Self {
            grid,
            hbar,
            mass,
            method,
            origin: [0.0; 3],
            psi: vec![(0.0, 0.0); len],
            t: 0.0,
            last_iterations: 0,
            order: 2,
            potential: vec![0.0; len],
            potential_fn: None,
            potential_t: f64::NAN,
            absorber: vec![0.0; len],
            free,
            k2,
            kinetic_cache: Vec::new(),
            kick_cache: Vec::new(),
            backup: Vec::new(),
        })
    }

    /// `ħ² / 2m`.
    fn beta(&self) -> f64 {
        self.hbar * self.hbar / (2.0 * self.mass)
    }

    /// The volume element `h^d`.
    pub fn cell_volume(&self) -> f64 {
        self.grid.h.powi(self.grid.dims() as i32)
    }

    /// Coordinates of point `idx`: `origin + grid.coordinates(idx)`.
    pub fn coordinates(&self, idx: usize) -> [f64; 3] {
        let x = self.grid.coordinates(idx);
        [
            self.origin[0] + x[0],
            self.origin[1] + x[1],
            self.origin[2] + x[2],
        ]
    }

    /// Order of the split-step composition (2 or 4).
    pub fn order(&self) -> usize {
        self.order
    }

    pub fn set_order(&mut self, order: usize) -> Result<()> {
        if order != 2 && order != 4 {
            return invalid(format!("order must be 2 or 4, got {order}"));
        }
        if order == 4 && self.method != QuantumMethod::SplitStep {
            return invalid("order 4 is only available with the split-step method");
        }
        self.order = order;
        Ok(())
    }

    /// Sets ψ (its values on Dirichlet boundary nodes are replaced by zero).
    pub fn set_psi(&mut self, psi: Vec<Complex>) -> Result<()> {
        if psi.len() != self.grid.len() {
            return invalid(format!(
                "psi: expected {} values for a grid of shape {:?}, got {}",
                self.grid.len(),
                self.grid.shape(),
                psi.len()
            ));
        }
        if psi.iter().any(|v| !(v.0.is_finite() && v.1.is_finite())) {
            return invalid("psi must be finite");
        }
        self.psi = psi;
        for (v, &free) in self.psi.iter_mut().zip(&self.free) {
            if !free {
                *v = (0.0, 0.0);
            }
        }
        Ok(())
    }

    /// Scales ψ to unit norm.
    pub fn normalize(&mut self) -> Result<()> {
        let n = self.norm();
        if n.is_nan() || n <= 0.0 {
            return invalid("cannot normalise a zero wavefunction");
        }
        let s = 1.0 / n.sqrt();
        self.psi.iter_mut().for_each(|v| *v = (v.0 * s, v.1 * s));
        Ok(())
    }

    /// A normalised Gaussian wave packet `Π exp(-(x - c)²/4σ²) e^{ip·x/ħ}`: centre `center`,
    /// position spread `sigma` and mean momentum `momentum` along each axis.
    pub fn set_gaussian(
        &mut self,
        center: [f64; 3],
        sigma: [f64; 3],
        momentum: [f64; 3],
    ) -> Result<()> {
        let dims = self.grid.dims();
        if sigma[..dims].iter().any(|&s| !(s > 0.0 && s.is_finite())) {
            return invalid(format!(
                "packet width must be positive, got {:?}",
                &sigma[..dims]
            ));
        }
        let psi = (0..self.grid.len())
            .map(|idx| {
                let x = self.coordinates(idx);
                let (mut amp, mut phase) = (0.0, 0.0);
                for a in 0..dims {
                    let d = x[a] - center[a];
                    amp -= d * d / (4.0 * sigma[a] * sigma[a]);
                    phase += momentum[a] * x[a] / self.hbar;
                }
                (C::cis(phase) * amp.exp()).into()
            })
            .collect();
        self.set_psi(psi)?;
        self.normalize()
    }

    /// The potential at the time it was last evaluated (the current time after a step).
    pub fn potential(&self) -> &[f64] {
        &self.potential
    }

    /// Sets a static potential (replacing any potential function).
    pub fn set_potential(&mut self, v: Vec<f64>) -> Result<()> {
        self.grid.check(&v, "potential")?;
        if v.iter().any(|x| !x.is_finite()) {
            return invalid("potential must be finite");
        }
        self.potential = v;
        self.potential_fn = None;
        self.kick_cache.clear();
        Ok(())
    }

    /// Sets a time-dependent potential, evaluated now and then whenever the solver needs it.
    pub fn set_potential_fn(&mut self, f: PotentialFn) -> Result<()> {
        self.potential_fn = Some(f);
        self.potential_t = f64::NAN;
        self.kick_cache.clear();
        self.refresh_potential(self.t)
    }

    /// Sets the clock (re-evaluating a time-dependent potential).
    pub fn set_time(&mut self, t: f64) -> Result<()> {
        if !t.is_finite() {
            return invalid(format!("t must be finite, got {t}"));
        }
        self.t = t;
        self.refresh_potential(t)
    }

    /// Whether the potential is a function of time.
    pub fn time_dependent(&self) -> bool {
        self.potential_fn.is_some()
    }

    fn refresh_potential(&mut self, t: f64) -> Result<()> {
        if let Some(f) = self.potential_fn.as_mut() {
            if self.potential_t.to_bits() != t.to_bits() {
                f(t, &mut self.potential)?;
                if self.potential.iter().any(|x| !x.is_finite()) {
                    self.potential_t = f64::NAN;
                    return invalid(format!("potential at t = {t} is not finite"));
                }
                self.potential_t = t;
            }
        }
        Ok(())
    }

    /// The absorbing potential `W`.
    pub fn absorber(&self) -> &[f64] {
        &self.absorber
    }

    /// Sets the absorbing potential `W ≥ 0` (all zeros: none).
    pub fn set_absorber(&mut self, w: Vec<f64>) -> Result<()> {
        self.grid.check(&w, "absorber")?;
        if w.iter().any(|x| !(x.is_finite() && *x >= 0.0)) {
            return invalid("absorber must be finite and non-negative");
        }
        self.absorber = w;
        self.kick_cache.clear();
        Ok(())
    }

    /// A quadratic absorbing layer `W = strength ((width - d) / width)²` within `width` of
    /// the edges of the grid (`d`: distance to the nearest edge).
    pub fn set_absorbing_layer(&mut self, width: f64, strength: f64) -> Result<()> {
        if !(width > 0.0 && width.is_finite()) {
            return invalid(format!(
                "absorbing layer width must be positive, got {width}"
            ));
        }
        if !(strength >= 0.0 && strength.is_finite()) {
            return invalid(format!(
                "absorbing layer strength must be non-negative, got {strength}"
            ));
        }
        let g = self.grid;
        let strides = [g.n[1] * g.n[2], g.n[2], 1];
        let w = (0..g.len())
            .map(|idx| {
                let d = (0..g.dims())
                    .map(|a| {
                        let i = (idx / strides[a]) % g.n[a];
                        i.min(g.n[a] - 1 - i) as f64 * g.h
                    })
                    .fold(f64::INFINITY, f64::min);
                if d < width {
                    strength * ((width - d) / width).powi(2)
                } else {
                    0.0
                }
            })
            .collect();
        self.set_absorber(w)
    }

    /// Takes `n` steps of `dt`. If a step fails (a potential function raises, or the linear
    /// solver does not converge), ψ and `t` are left at the last completed step.
    pub fn step(&mut self, dt: f64, n: usize) -> Result<()> {
        if !(dt > 0.0 && dt.is_finite()) {
            return invalid(format!("dt must be positive, got {dt}"));
        }
        for _ in 0..n {
            let t0 = self.t;
            self.backup.clone_from(&self.psi);
            let result = match self.method {
                QuantumMethod::SplitStep => self.split_step(dt),
                QuantumMethod::CrankNicolson => self.crank_nicolson(dt),
            };
            if let Err(e) = result {
                std::mem::swap(&mut self.psi, &mut self.backup);
                self.t = t0;
                return Err(e);
            }
            self.t = t0 + dt;
        }
        // Observables use the potential at the current time.
        self.refresh_potential(self.t)
    }

    fn split_step(&mut self, dt: f64) -> Result<()> {
        let weights = if self.order == 4 {
            yoshida_weights().to_vec()
        } else {
            vec![1.0]
        };
        for w in weights {
            let tau = w * dt;
            self.kick(0.5 * tau)?;
            self.drift(tau)?;
            self.t += tau;
            self.kick(0.5 * tau)?;
        }
        Ok(())
    }

    /// ψ ← e^{-(iV + W)τ/ħ} ψ with the potential at the current time.
    fn kick(&mut self, tau: f64) -> Result<()> {
        let c = tau / self.hbar;
        if self.potential_fn.is_some() {
            self.refresh_potential(self.t)?;
            for ((p, &v), &w) in self.psi.iter_mut().zip(&self.potential).zip(&self.absorber) {
                *p = (C::from(*p) * C::cis(-v * c) * (-w * c).exp()).into();
            }
            return Ok(());
        }
        let pos = match self.kick_cache.iter().position(|(s, _)| *s == tau) {
            Some(pos) => pos,
            None => {
                let factors = self
                    .potential
                    .iter()
                    .zip(&self.absorber)
                    .map(|(&v, &w)| C::cis(-v * c) * (-w * c).exp())
                    .collect();
                push_cache(&mut self.kick_cache, tau, factors)
            }
        };
        for (p, f) in self.psi.iter_mut().zip(&self.kick_cache[pos].1) {
            *p = (C::from(*p) * *f).into();
        }
        Ok(())
    }

    /// ψ ← e^{-iTτ/ħ} ψ, exactly in Fourier space.
    fn drift(&mut self, tau: f64) -> Result<()> {
        let pos = match self.kinetic_cache.iter().position(|(s, _)| *s == tau) {
            Some(pos) => pos,
            None => {
                let c = self.hbar * tau / (2.0 * self.mass);
                let factors = self.k2.iter().map(|&k2| C::cis(-k2 * c)).collect();
                push_cache(&mut self.kinetic_cache, tau, factors)
            }
        };
        fft_nd(&mut self.psi, self.grid.n, false)?;
        for (p, f) in self.psi.iter_mut().zip(&self.kinetic_cache[pos].1) {
            *p = (C::from(*p) * *f).into();
        }
        fft_nd(&mut self.psi, self.grid.n, true)
    }

    /// `out = Hψ` with finite differences (zero at fixed nodes).
    fn apply_h_fd(&self, psi: &[C], out: &mut [C]) {
        self.grid.laplacian(psi, out);
        let beta = self.beta();
        for i in 0..psi.len() {
            out[i] = if self.free[i] {
                out[i] * -beta + psi[i] * C::new(self.potential[i], -self.absorber[i])
            } else {
                C::ZERO
            };
        }
    }

    fn crank_nicolson(&mut self, dt: f64) -> Result<()> {
        self.refresh_potential(self.t + 0.5 * dt)?;
        let a = 0.5 * dt / self.hbar;
        let ia = C::new(0.0, a);
        let psi: Vec<C> = self.psi.iter().map(|&v| C::from(v)).collect();
        let mut h_psi = vec![C::ZERO; psi.len()];
        self.apply_h_fd(&psi, &mut h_psi);
        // (1 + iaH) ψ' = (1 - iaH) ψ; the matrix is the identity at fixed nodes.
        let rhs: Vec<C> = (0..psi.len())
            .map(|i| {
                if self.free[i] {
                    psi[i] - ia * h_psi[i]
                } else {
                    C::ZERO
                }
            })
            .collect();
        let beta_h2 = self.beta() / (self.grid.h * self.grid.h);
        let diag: Vec<C> = (0..psi.len())
            .map(|i| {
                if self.free[i] {
                    let h_ii = C::new(
                        2.0 * self.grid.dims() as f64 * beta_h2 + self.potential[i],
                        -self.absorber[i],
                    );
                    C::new(1.0, 0.0) + ia * h_ii
                } else {
                    C::new(1.0, 0.0)
                }
            })
            .collect();
        let mut x = psi;
        if self.grid.dims() == 1 {
            self.solve_tridiagonal(&diag, ia * -beta_h2, &rhs, &mut x)?;
            self.last_iterations = 0;
        } else {
            let mut scratch = vec![C::ZERO; x.len()];
            let apply = |v: &[C], out: &mut [C]| {
                self.apply_h_fd(v, &mut scratch);
                for i in 0..v.len() {
                    out[i] = if self.free[i] {
                        v[i] + ia * scratch[i]
                    } else {
                        v[i]
                    };
                }
            };
            self.last_iterations = bicgstab(apply, &diag, &rhs, &mut x, SOLVE_TOL)?;
        }
        for (p, v) in self.psi.iter_mut().zip(&x) {
            *p = (*v).into();
        }
        Ok(())
    }

    /// Solves the 1D Crank-Nicolson system: `diag` on the diagonal (assuming two neighbours
    /// per point) and `off` off it, with the grid's boundary condition.
    fn solve_tridiagonal(&self, diag: &[C], off: C, rhs: &[C], x: &mut [C]) -> Result<()> {
        let n = diag.len();
        match self.grid.boundary {
            Boundary::Dirichlet => {
                // Unknowns are the interior nodes; the edge nodes stay zero.
                thomas(&diag[1..n - 1], off, &rhs[1..n - 1], &mut x[1..n - 1]);
                x[0] = C::ZERO;
                x[n - 1] = C::ZERO;
            }
            Boundary::Neumann => {
                // The ghost value equals the edge value: one neighbour term folds into the
                // diagonal at each end.
                let mut d = diag.to_vec();
                d[0] += off;
                d[n - 1] += off;
                thomas(&d, off, rhs, x);
            }
            Boundary::Periodic if n < 3 => {
                return invalid("periodic Crank-Nicolson needs at least 3 points")
            }
            Boundary::Periodic => cyclic_thomas(diag, off, rhs, x),
        }
        if x.iter().any(|v| !(v.re.is_finite() && v.im.is_finite())) {
            return invalid("the Crank-Nicolson solve produced non-finite values");
        }
        Ok(())
    }

    /// Fourier transform of ψ (split-step grids only).
    fn spectrum(&self) -> Vec<Complex> {
        let mut data = self.psi.clone();
        // Sizes were checked when the solver was built.
        fft_nd(&mut data, self.grid.n, false).expect("split-step grid sizes are powers of two");
        data
    }

    /// `Σ |ψ|² h^d` (1 for a normalised state; decreases when an absorber removes probability).
    pub fn norm(&self) -> f64 {
        self.psi.iter().map(|&v| C::from(v).norm2()).sum::<f64>() * self.cell_volume()
    }

    fn weight_sum(&self) -> Result<f64> {
        let s: f64 = self.psi.iter().map(|&v| C::from(v).norm2()).sum();
        if s.is_nan() || s <= 0.0 {
            return invalid("the wavefunction is zero");
        }
        Ok(s)
    }

    /// `⟨x⟩` along each axis (zero for unused axes), normalised by the current norm.
    pub fn position(&self) -> Result<[f64; 3]> {
        let total = self.weight_sum()?;
        let mut m = [0.0; 3];
        for (idx, &v) in self.psi.iter().enumerate() {
            let p = C::from(v).norm2();
            let x = self.coordinates(idx);
            for a in 0..self.grid.dims() {
                m[a] += p * x[a];
            }
        }
        Ok(m.map(|v| v / total))
    }

    /// `⟨x²⟩ - ⟨x⟩²` along each axis.
    pub fn position_variance(&self) -> Result<[f64; 3]> {
        let total = self.weight_sum()?;
        let mean = self.position()?;
        let mut var = [0.0; 3];
        for (idx, &v) in self.psi.iter().enumerate() {
            let p = C::from(v).norm2();
            let x = self.coordinates(idx);
            for a in 0..self.grid.dims() {
                var[a] += p * (x[a] - mean[a]).powi(2);
            }
        }
        Ok(var.map(|v| v / total))
    }

    /// `⟨p⟩` along each axis: spectral for split-step, central differences otherwise.
    pub fn momentum(&self) -> Result<[f64; 3]> {
        let total = self.weight_sum()?;
        let g = self.grid;
        let strides = [g.n[1] * g.n[2], g.n[2], 1];
        let mut p = [0.0; 3];
        if self.method == QuantumMethod::SplitStep {
            let spec = self.spectrum();
            let mut weight = 0.0;
            for (idx, &v) in spec.iter().enumerate() {
                let w = C::from(v).norm2();
                weight += w;
                for a in 0..g.dims() {
                    let n = g.n[a];
                    let j = (idx / strides[a]) % n;
                    // The Nyquist mode has no definite sign: it carries no mean momentum.
                    let j = if j < n / 2 {
                        j as f64
                    } else if j == n / 2 {
                        0.0
                    } else {
                        j as f64 - n as f64
                    };
                    p[a] += w * 2.0 * std::f64::consts::PI * j / (n as f64 * g.h);
                }
            }
            return Ok(p.map(|v| self.hbar * v / weight));
        }
        for (idx, &v) in self.psi.iter().enumerate() {
            if !self.free[idx] {
                continue;
            }
            let psi = C::from(v);
            for (a, pa) in p.iter_mut().enumerate().take(g.dims()) {
                let (up, down) = neighbours(&g, idx, a);
                let up = up.map_or(C::ZERO, |j| C::from(self.psi[j]));
                let down = down.map_or(C::ZERO, |j| C::from(self.psi[j]));
                *pa += (psi.conj() * (up - down)).im;
            }
        }
        Ok(p.map(|v| self.hbar * v / (2.0 * g.h * total)))
    }

    /// `⟨T⟩`, with the kinetic operator of the method (spectral or finite differences).
    pub fn kinetic_energy(&self) -> Result<f64> {
        let total = self.weight_sum()?;
        if self.method == QuantumMethod::SplitStep {
            let spec = self.spectrum();
            let (mut num, mut den) = (0.0, 0.0);
            for (&v, &k2) in spec.iter().zip(&self.k2) {
                let w = C::from(v).norm2();
                num += w * k2;
                den += w;
            }
            return Ok(self.beta() * num / den);
        }
        let psi: Vec<C> = self.psi.iter().map(|&v| C::from(v)).collect();
        let mut lap = vec![C::ZERO; psi.len()];
        self.grid.laplacian(&psi, &mut lap);
        Ok(-self.beta() * dot(&psi, &lap).re / total)
    }

    /// `⟨V⟩` at the current time.
    pub fn potential_energy(&self) -> Result<f64> {
        let total = self.weight_sum()?;
        let e: f64 = self
            .psi
            .iter()
            .zip(&self.potential)
            .map(|(&p, &v)| C::from(p).norm2() * v)
            .sum();
        Ok(e / total)
    }

    /// `⟨H⟩ = ⟨T⟩ + ⟨V⟩` (the absorber is not included).
    pub fn energy(&self) -> Result<f64> {
        Ok(self.kinetic_energy()? + self.potential_energy()?)
    }

    /// `out = Hx` for a real vector (no absorber), with the method's discretisation.
    fn apply_h_real(&self, x: &[f64], out: &mut [f64]) {
        let beta = self.beta();
        if self.method == QuantumMethod::SplitStep {
            let mut data: Vec<Complex> = x.iter().map(|&v| (v, 0.0)).collect();
            fft_nd(&mut data, self.grid.n, false).expect("split-step grid sizes are powers of two");
            for (d, &k2) in data.iter_mut().zip(&self.k2) {
                *d = (d.0 * beta * k2, d.1 * beta * k2);
            }
            fft_nd(&mut data, self.grid.n, true).expect("split-step grid sizes are powers of two");
            for i in 0..x.len() {
                out[i] = data[i].0 + self.potential[i] * x[i];
            }
        } else {
            self.grid.laplacian(x, out);
            for i in 0..x.len() {
                out[i] = if self.free[i] {
                    -beta * out[i] + self.potential[i] * x[i]
                } else {
                    0.0
                };
            }
        }
    }

    /// The `count` lowest eigenstates of `H` (with the potential at the current time and
    /// without the absorber), as `(energies, states)`: energies ascending, states real and
    /// normalised (`Σ ψ² h^d = 1`), with their largest value positive.
    ///
    /// LOBPCG (locally optimal block preconditioned conjugate gradients) on a block with extra
    /// guard vectors, until every residual `‖Hψ - Eψ‖` is below `tol (E - σ)`, where `σ` is a
    /// shift just below the potential minimum. Preconditioner: `(βk² + c)⁻¹` in Fourier space
    /// for the spectral operator; an exact tridiagonal solve of `H - σ` for 1D finite
    /// differences; a few Jacobi-preconditioned conjugate-gradient iterations otherwise.
    pub fn eigenstates(
        &mut self,
        count: usize,
        tol: f64,
        max_iter: usize,
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
        self.refresh_potential(self.t)?;
        let len = self.grid.len();
        let n_free = self.free.iter().filter(|&&f| f).count();
        if count == 0 || count > n_free {
            return invalid(format!(
                "count must be between 1 and the number of free grid points ({n_free}), got {count}"
            ));
        }
        if !(tol > 0.0 && tol.is_finite()) {
            return invalid(format!("tol must be positive, got {tol}"));
        }
        // Guard vectors speed up convergence of the highest wanted states.
        let m = (count + (count / 2).max(2)).min(n_free);
        let beta = self.beta();
        let g = self.grid;
        let extent = (0..g.dims())
            .map(|a| g.n[a] as f64 * g.h)
            .fold(0.0, f64::max);
        let v_min = (0..len)
            .filter(|&i| self.free[i])
            .map(|i| self.potential[i])
            .fold(f64::INFINITY, f64::min);
        // H - σ ≥ β / extent² > 0, so the shifted operator is positive definite.
        let sigma = v_min - beta / (extent * extent);
        let apply = |x: &[f64], out: &mut [f64]| self.apply_h_real(x, out);
        let precond = self.eigen_preconditioner(sigma);

        let mut x: Vec<Vec<f64>> = (0..m)
            .map(|j| {
                (0..len)
                    .map(|i| {
                        if self.free[i] {
                            normal(0x5eed, j as u64, i as u64)
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect();
        orthonormalize(&mut x)?;
        let mut hx: Vec<Vec<f64>> = x.iter().map(|v| apply_new(&apply, v)).collect();
        let mut theta = rayleigh_ritz(&mut x, &mut hx, m);
        let mut p: Vec<Vec<f64>> = Vec::new();
        let mut worst = f64::INFINITY;
        for _ in 0..max_iter {
            let r: Vec<Vec<f64>> = (0..m)
                .map(|j| {
                    hx[j]
                        .iter()
                        .zip(&x[j])
                        .map(|(h, v)| h - theta[j] * v)
                        .collect()
                })
                .collect();
            worst = (0..count)
                .map(|j| dot_real(&r[j], &r[j]).sqrt() / (theta[j] - sigma))
                .fold(0.0, f64::max);
            if worst <= tol {
                let dv = self.cell_volume();
                let states = x
                    .into_iter()
                    .take(count)
                    .map(|mut v| {
                        let big = v.iter().fold(0.0f64, |a, b| a.max(b.abs()));
                        let first = v.iter().find(|e| e.abs() > 0.5 * big).copied();
                        let sign = if first.unwrap_or(1.0) < 0.0 {
                            -1.0
                        } else {
                            1.0
                        };
                        let s = sign / (dot_real(&v, &v) * dv).sqrt();
                        v.iter_mut().for_each(|e| *e *= s);
                        v
                    })
                    .collect();
                theta.truncate(count);
                return Ok((theta, states));
            }
            // Search directions: preconditioned residuals and the previous step, made
            // orthonormal to X and to each other (nearly dependent ones are dropped).
            let mut extra: Vec<Vec<f64>> = r
                .iter()
                .map(|rj| {
                    let mut w = vec![0.0; len];
                    precond(rj, &mut w);
                    w
                })
                .collect();
            extra.append(&mut p);
            orthonormalize_against(&x, &mut extra);
            let he: Vec<Vec<f64>> = extra.iter().map(|v| apply_new(&apply, v)).collect();
            // Rayleigh-Ritz in span[X, W, P].
            let k = m + extra.len();
            let mut a = vec![vec![0.0; k]; k];
            for (i, &t) in theta.iter().enumerate() {
                a[i][i] = t;
            }
            for (i, hxi) in hx.iter().enumerate() {
                for (j, e) in extra.iter().enumerate() {
                    let v = dot_real(hxi, e);
                    a[i][m + j] = v;
                    a[m + j][i] = v;
                }
            }
            for i in 0..extra.len() {
                for j in i..extra.len() {
                    let v = 0.5 * (dot_real(&extra[i], &he[j]) + dot_real(&he[i], &extra[j]));
                    a[m + i][m + j] = v;
                    a[m + j][m + i] = v;
                }
            }
            let (evals, c) = jacobi_eigen(a);
            let combine = |basis: &[Vec<f64>], first_row: usize, col: usize| {
                let mut out = vec![0.0; len];
                for (i, b) in basis.iter().enumerate() {
                    let coef = c[first_row + i][col];
                    out.iter_mut().zip(b).for_each(|(o, v)| *o += coef * v);
                }
                out
            };
            let mut new_x = Vec::with_capacity(m);
            let mut new_hx = Vec::with_capacity(m);
            for col in 0..m {
                // The new step P is the part of the update outside span(X).
                let pe = combine(&extra, m, col);
                let hpe = combine(&he, m, col);
                let mut xc = combine(&x, 0, col);
                let mut hxc = combine(&hx, 0, col);
                xc.iter_mut().zip(&pe).for_each(|(a, b)| *a += b);
                hxc.iter_mut().zip(&hpe).for_each(|(a, b)| *a += b);
                new_x.push(xc);
                new_hx.push(hxc);
                p.push(pe);
            }
            x = new_x;
            hx = new_hx;
            theta = evals[..m].to_vec();
            // Rounding slowly erodes orthonormality; restore it (and HX with it) when needed.
            if x.iter().any(|v| (dot_real(v, v) - 1.0).abs() > 1e-10) {
                orthonormalize(&mut x)?;
                hx = x.iter().map(|v| apply_new(&apply, v)).collect();
                theta = rayleigh_ritz(&mut x, &mut hx, m);
            }
        }
        invalid(format!(
            "eigenstates did not converge in {max_iter} iterations (relative residual {worst:.2e}, \
             tol {tol:.0e})"
        ))
    }

    /// `z ≈ (H - σ)⁻¹ r` for [`Self::eigenstates`].
    fn eigen_preconditioner(&self, sigma: f64) -> impl Fn(&[f64], &mut [f64]) + '_ {
        let g = self.grid;
        let len = g.len();
        let beta = self.beta();
        let beta_h2 = beta / (g.h * g.h);
        let shift_mean = (0..len).map(|i| self.potential[i] - sigma).sum::<f64>() / len as f64;
        let diag: Vec<f64> = (0..len)
            .map(|i| {
                if self.free[i] {
                    2.0 * g.dims() as f64 * beta_h2 + self.potential[i] - sigma
                } else {
                    1.0
                }
            })
            .collect();
        move |r: &[f64], z: &mut [f64]| {
            if !self.k2.is_empty() {
                let mut data: Vec<Complex> = r.iter().map(|&v| (v, 0.0)).collect();
                fft_nd(&mut data, g.n, false).expect("split-step grid sizes are powers of two");
                for (d, &k2) in data.iter_mut().zip(&self.k2) {
                    let w = 1.0 / (beta * k2 + shift_mean);
                    *d = (d.0 * w, d.1 * w);
                }
                fft_nd(&mut data, g.n, true).expect("split-step grid sizes are powers of two");
                z.iter_mut().zip(&data).for_each(|(z, d)| *z = d.0);
            } else if g.dims() == 1 {
                let d: Vec<C> = diag.iter().map(|&v| C::new(v, 0.0)).collect();
                let rhs: Vec<C> = r.iter().map(|&v| C::new(v, 0.0)).collect();
                let mut x = vec![C::ZERO; len];
                // H - σ is positive definite, so the tridiagonal solve cannot break down; were
                // it to, LOBPCG would only lose this preconditioned direction.
                if self
                    .solve_tridiagonal(&d, C::new(-beta_h2, 0.0), &rhs, &mut x)
                    .is_err()
                {
                    z.copy_from_slice(r);
                    return;
                }
                z.iter_mut().zip(&x).for_each(|(z, x)| *z = x.re);
            } else {
                // A few conjugate-gradient iterations on H - σ.
                z.iter_mut().for_each(|v| *v = 0.0);
                let apply = |v: &[f64], out: &mut [f64]| {
                    self.apply_h_real(v, out);
                    for i in 0..len {
                        out[i] = if self.free[i] {
                            out[i] - sigma * v[i]
                        } else {
                            v[i]
                        };
                    }
                };
                let jacobi = |v: &[f64], out: &mut [f64]| {
                    out.iter_mut()
                        .zip(v.iter().zip(&diag))
                        .for_each(|(o, (v, d))| *o = v / d);
                };
                pcg(apply, jacobi, r, z, 1e-2, 30);
            }
        }
    }
}

/// Neighbours of `idx` along `axis` (`None` beyond a Dirichlet edge; the point itself
/// beyond a Neumann edge, whose ghost value equals the edge value).
fn neighbours(g: &Grid, idx: usize, axis: usize) -> (Option<usize>, Option<usize>) {
    let strides = [g.n[1] * g.n[2], g.n[2], 1];
    let (n, stride) = (g.n[axis], strides[axis]);
    let i = (idx / stride) % n;
    let edge = |wrapped: usize| match g.boundary {
        Boundary::Periodic => Some(wrapped),
        Boundary::Neumann => Some(idx),
        Boundary::Dirichlet => None,
    };
    let up = if i + 1 < n {
        Some(idx + stride)
    } else {
        edge(idx + stride - n * stride)
    };
    let down = if i > 0 {
        Some(idx - stride)
    } else {
        edge(idx + (n - 1) * stride)
    };
    (up, down)
}

fn push_cache(cache: &mut Vec<(f64, Vec<C>)>, tau: f64, factors: Vec<C>) -> usize {
    // A few step sizes at most are in use at once (Yoshida uses three).
    if cache.len() >= 6 {
        cache.clear();
    }
    cache.push((tau, factors));
    cache.len() - 1
}

/// Solves a tridiagonal system with diagonal `diag` and constant off-diagonal `off`.
fn thomas(diag: &[C], off: C, rhs: &[C], x: &mut [C]) {
    let n = diag.len();
    let mut cp = vec![C::ZERO; n];
    let mut dp = vec![C::ZERO; n];
    cp[0] = off / diag[0];
    dp[0] = rhs[0] / diag[0];
    for i in 1..n {
        let m = diag[i] - off * cp[i - 1];
        cp[i] = off / m;
        dp[i] = (rhs[i] - off * dp[i - 1]) / m;
    }
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
}

/// [`thomas`] with the corners `A[0][n-1] = A[n-1][0] = off` (periodic), by Sherman-Morrison.
fn cyclic_thomas(diag: &[C], off: C, rhs: &[C], x: &mut [C]) {
    let n = diag.len();
    let gamma = -diag[0];
    let mut d = diag.to_vec();
    d[0] -= gamma;
    d[n - 1] -= off * off / gamma;
    thomas(&d, off, rhs, x);
    let mut u = vec![C::ZERO; n];
    u[0] = gamma;
    u[n - 1] = off;
    let mut z = vec![C::ZERO; n];
    thomas(&d, off, &u, &mut z);
    let fact = (x[0] + off * x[n - 1] / gamma) / (C::new(1.0, 0.0) + z[0] + off * z[n - 1] / gamma);
    for (xi, zi) in x.iter_mut().zip(&z) {
        *xi -= fact * *zi;
    }
}

/// Solves `A x = b` by BiCGSTAB with the Jacobi preconditioner `diag(A)`, starting from `x`.
/// Returns the iterations used.
fn bicgstab(
    mut apply: impl FnMut(&[C], &mut [C]),
    diag: &[C],
    b: &[C],
    x: &mut [C],
    tol: f64,
) -> Result<usize> {
    let n = b.len();
    let mut r = vec![C::ZERO; n];
    apply(x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }
    let norm_b = norm(b).max(f64::MIN_POSITIVE);
    if norm(&r) <= tol * norm_b {
        return Ok(0);
    }
    let r_hat = r.clone();
    let (mut rho, mut alpha, mut omega) = (C::new(1.0, 0.0), C::new(1.0, 0.0), C::new(1.0, 0.0));
    let mut v = vec![C::ZERO; n];
    let mut p = vec![C::ZERO; n];
    let mut y = vec![C::ZERO; n];
    let mut z = vec![C::ZERO; n];
    let mut s = vec![C::ZERO; n];
    let mut t = vec![C::ZERO; n];
    let max_iter = 10 * n + 100;
    for it in 1..=max_iter {
        let rho_new = dot(&r_hat, &r);
        if rho_new.norm2() == 0.0 {
            return invalid("BiCGSTAB broke down");
        }
        let beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
            y[i] = p[i] / diag[i];
        }
        apply(&y, &mut v);
        alpha = rho / dot(&r_hat, &v);
        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }
        if norm(&s) <= tol * norm_b {
            for i in 0..n {
                x[i] += alpha * y[i];
            }
            return Ok(it);
        }
        for i in 0..n {
            z[i] = s[i] / diag[i];
        }
        apply(&z, &mut t);
        let tt = dot(&t, &t);
        if tt.norm2() == 0.0 {
            return invalid("BiCGSTAB broke down");
        }
        omega = dot(&t, &s) / tt;
        for i in 0..n {
            x[i] += alpha * y[i] + omega * z[i];
            r[i] = s[i] - omega * t[i];
        }
        if norm(&r) <= tol * norm_b {
            return Ok(it);
        }
    }
    invalid(format!(
        "BiCGSTAB did not converge in {max_iter} iterations"
    ))
}

fn dot_real(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// At most `max_iter` iterations of preconditioned conjugate gradients on the symmetric
/// positive definite system `A x = b` (`precond`: `z ≈ A⁻¹ r`), starting from `x` and stopping
/// once the residual is below `tol ‖b‖`. Returns the iterations used.
fn pcg(
    mut apply: impl FnMut(&[f64], &mut [f64]),
    mut precond: impl FnMut(&[f64], &mut [f64]),
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iter: usize,
) -> usize {
    let n = b.len();
    let mut r = vec![0.0; n];
    apply(x, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
    }
    let norm_b = dot_real(b, b).sqrt().max(f64::MIN_POSITIVE);
    let mut z = vec![0.0; n];
    precond(&r, &mut z);
    let mut p = z.clone();
    let mut rz = dot_real(&r, &z);
    let mut ap = vec![0.0; n];
    for it in 0..max_iter {
        if dot_real(&r, &r).sqrt() <= tol * norm_b {
            return it;
        }
        apply(&p, &mut ap);
        let pap = dot_real(&p, &ap);
        if pap <= 0.0 {
            return it;
        }
        let step = rz / pap;
        for i in 0..n {
            x[i] += step * p[i];
            r[i] -= step * ap[i];
        }
        precond(&r, &mut z);
        let rz_new = dot_real(&r, &z);
        let mix = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + mix * p[i];
        }
    }
    max_iter
}

fn apply_new(apply: &impl Fn(&[f64], &mut [f64]), x: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; x.len()];
    apply(x, &mut out);
    out
}

/// Rotates the orthonormal block `x` (with `hx = H x`) to the lowest `m` Ritz vectors of `H`
/// in its span; returns their Ritz values, ascending.
fn rayleigh_ritz(x: &mut Vec<Vec<f64>>, hx: &mut Vec<Vec<f64>>, m: usize) -> Vec<f64> {
    let k = x.len();
    let mut a = vec![vec![0.0; k]; k];
    for i in 0..k {
        for j in i..k {
            let v = 0.5 * (dot_real(&x[i], &hx[j]) + dot_real(&hx[i], &x[j]));
            a[i][j] = v;
            a[j][i] = v;
        }
    }
    let (evals, c) = jacobi_eigen(a);
    let len = x[0].len();
    let rotate = |basis: &[Vec<f64>]| -> Vec<Vec<f64>> {
        (0..m)
            .map(|col| {
                let mut out = vec![0.0; len];
                for (b, row) in basis.iter().zip(&c) {
                    out.iter_mut().zip(b).for_each(|(o, v)| *o += row[col] * v);
                }
                out
            })
            .collect()
    };
    *x = rotate(x);
    *hx = rotate(hx);
    evals[..m].to_vec()
}

/// Makes `extra` orthonormal to the orthonormal block `x` and to itself (Gram-Schmidt, twice),
/// dropping vectors that are nearly dependent.
fn orthonormalize_against(x: &[Vec<f64>], extra: &mut Vec<Vec<f64>>) {
    let mut kept: Vec<Vec<f64>> = Vec::with_capacity(extra.len());
    for mut v in extra.drain(..) {
        let before = dot_real(&v, &v).sqrt();
        if !(before > 0.0 && before.is_finite()) {
            continue;
        }
        for _ in 0..2 {
            for b in x.iter().chain(kept.iter()) {
                let c = dot_real(b, &v);
                v.iter_mut().zip(b).for_each(|(a, b)| *a -= c * b);
            }
        }
        let after = dot_real(&v, &v).sqrt();
        if after > 1e-8 * before {
            v.iter_mut().for_each(|a| *a /= after);
            kept.push(v);
        }
    }
    *extra = kept;
}

/// Modified Gram-Schmidt, twice for stability.
fn orthonormalize(v: &mut [Vec<f64>]) -> Result<()> {
    for j in 0..v.len() {
        let (done, rest) = v.split_at_mut(j);
        let vj = &mut rest[0];
        for _ in 0..2 {
            for vk in done.iter() {
                let c = dot_real(vk, vj);
                vj.iter_mut().zip(vk).for_each(|(a, b)| *a -= c * b);
            }
        }
        let n = dot_real(vj, vj).sqrt();
        if !(n > 0.0 && n.is_finite()) {
            return invalid("eigenstate basis became linearly dependent");
        }
        vj.iter_mut().for_each(|a| *a /= n);
    }
    Ok(())
}

/// Eigenvalues (ascending) and eigenvectors (columns of `q`) of a small symmetric matrix, by
/// cyclic Jacobi rotations.
fn jacobi_eigen(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let m = a.len();
    let mut q: Vec<Vec<f64>> = (0..m)
        .map(|i| (0..m).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let scale: f64 = a.iter().flatten().map(|v| v * v).sum::<f64>().sqrt();
    for _ in 0..100 {
        let off: f64 = (0..m)
            .flat_map(|i| (0..m).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum::<f64>()
            .sqrt();
        if off <= 1e-15 * scale {
            break;
        }
        for p in 0..m {
            for r in p + 1..m {
                if a[p][r] == 0.0 {
                    continue;
                }
                let theta = (a[r][r] - a[p][p]) / (2.0 * a[p][r]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for row in a.iter_mut() {
                    let (akp, akr) = (row[p], row[r]);
                    row[p] = c * akp - s * akr;
                    row[r] = s * akp + c * akr;
                }
                let (upper, lower) = a.split_at_mut(r);
                for (apk, ark) in upper[p].iter_mut().zip(lower[0].iter_mut()) {
                    let (x, y) = (*apk, *ark);
                    *apk = c * x - s * y;
                    *ark = s * x + c * y;
                }
                for row in q.iter_mut() {
                    let (qp, qr) = (row[p], row[r]);
                    row[p] = c * qp - s * qr;
                    row[r] = s * qp + c * qr;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&i, &j| a[i][i].total_cmp(&a[j][j]));
    let evals = order.iter().map(|&i| a[i][i]).collect();
    let q = q
        .iter()
        .map(|row| order.iter().map(|&i| row[i]).collect())
        .collect();
    (evals, q)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_diagonalises_a_symmetric_matrix() {
        let a = vec![
            vec![4.0, 1.0, 0.5],
            vec![1.0, 3.0, -0.2],
            vec![0.5, -0.2, 1.0],
        ];
        let (e, q) = jacobi_eigen(a.clone());
        assert!(e.windows(2).all(|w| w[0] <= w[1]));
        for (c, &lam) in e.iter().enumerate() {
            for i in 0..3 {
                let av: f64 = (0..3).map(|j| a[i][j] * q[j][c]).sum();
                assert!((av - lam * q[i][c]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn tridiagonal_solvers_match_the_matrix() {
        let n = 7;
        let diag: Vec<C> = (0..n)
            .map(|i| C::new(2.0 + i as f64, 0.3 * i as f64))
            .collect();
        let off = C::new(0.1, -0.7);
        let rhs: Vec<C> = (0..n).map(|i| C::new(i as f64, 1.0 - i as f64)).collect();
        for cyclic in [false, true] {
            let mut x = vec![C::ZERO; n];
            if cyclic {
                cyclic_thomas(&diag, off, &rhs, &mut x);
            } else {
                thomas(&diag, off, &rhs, &mut x);
            }
            for i in 0..n {
                let mut ax = diag[i] * x[i];
                if i > 0 {
                    ax += off * x[i - 1];
                } else if cyclic {
                    ax += off * x[n - 1];
                }
                if i + 1 < n {
                    ax += off * x[i + 1];
                } else if cyclic {
                    ax += off * x[0];
                }
                assert!((ax - rhs[i]).norm2() < 1e-24, "cyclic {cyclic} row {i}");
            }
        }
    }
}
