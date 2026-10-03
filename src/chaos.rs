//! Chaos indicators: Lyapunov exponents (Benettin's method) and MEGNO, from the
//! variational equations integrated alongside the trajectory.
//!
//! Each tangent vector `(δx, δv)` is carried as a block of extra "pseudo-particles" whose
//! acceleration is the Jacobian-vector product `∂a/∂x · δx + ∂a/∂v · δv` (analytic for most
//! built-in forces, finite differences otherwise; see [`ForceSet::jacobian_vector`]). The
//! world's own integrator advances the extended system, so for Verlet-type schemes the
//! tangent vectors follow the exact linearisation of the numerical map. Every
//! `renormalize_every` steps the tangent vectors are re-orthonormalised (Gram-Schmidt / QR),
//! and the logarithms of the stretching factors accumulate into the exponents.

use std::sync::Arc;

use crate::error::{invalid, Result};
use crate::forces::{Force, ForceSet};
use crate::state::State;
use crate::vec3::Vec3;
use crate::world::World;

/// Settings for [`World::lyapunov`].
#[derive(Debug, Clone, Copy)]
pub struct LyapunovOptions {
    /// How many exponents to compute (largest first); at most `6 N`.
    pub n_exponents: usize,
    /// Re-orthonormalise the tangent vectors every this many steps.
    pub renormalize_every: usize,
    /// Record the running estimates every this many steps (and at the end).
    pub record_every: usize,
    /// Seed for the initial tangent vectors.
    pub seed: u64,
}

impl Default for LyapunovOptions {
    fn default() -> Self {
        Self {
            n_exponents: 1,
            renormalize_every: 1,
            record_every: 100,
            seed: 1,
        }
    }
}

/// Result of [`World::lyapunov`].
#[derive(Debug, Clone, Default)]
pub struct LyapunovRun {
    /// Final estimates, largest first (per unit time).
    pub exponents: Vec<f64>,
    /// Times of the recorded estimates.
    pub t: Vec<f64>,
    /// Running estimates of every exponent at each recorded time.
    pub running: Vec<Vec<f64>>,
    /// MEGNO `Y(t)` of the leading tangent vector at each recorded time.
    pub megno: Vec<f64>,
    /// Its time average `⟨Y⟩(t)`: tends to 2 for quasi-periodic motion, 0 for an isochronous
    /// orbit, and grows as `λ t / 2` for chaotic motion.
    pub mean_megno: Vec<f64>,
}

/// The forces of the extended system: the real particles feel the world's forces, and
/// each tangent block feels the Jacobian-vector product along itself.
struct Variational {
    inner: Arc<ForceSet>,
    n: usize,
    blocks: usize,
}

impl Force for Variational {
    fn accumulate(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        let n = self.n;
        let free = vec![false; n]; // pinned entries are zeroed by the outer force set
        let mut a = Vec::new();
        self.inner
            .accelerations(t, &pos[..n], &vel[..n], &mass[..n], &free, &mut a)?;
        acc[..n].copy_from_slice(&a);
        for b in 1..=self.blocks {
            let r = b * n..(b + 1) * n;
            self.inner.jacobian_vector(
                t,
                &pos[..n],
                &vel[..n],
                &mass[..n],
                &free,
                &pos[r.clone()],
                &vel[r.clone()],
                &mut a,
            )?;
            acc[r].copy_from_slice(&a);
        }
        Ok(())
    }

    fn velocity_dependent(&self) -> bool {
        self.inner.velocity_dependent()
    }

    fn name(&self) -> String {
        "Variational".into()
    }
}

/// Tangent vector `b` of the extended state, flattened to `6 N` numbers.
fn tangent(s: &State, n: usize, b: usize) -> Vec<f64> {
    let r = (b + 1) * n..(b + 2) * n;
    s.pos[r.clone()]
        .iter()
        .chain(&s.vel[r])
        .flat_map(|v| v.to_array())
        .collect()
}

fn set_tangent(s: &mut State, n: usize, b: usize, w: &[f64]) {
    for i in 0..n {
        let p = (b + 1) * n + i;
        s.pos[p] = Vec3::new(w[3 * i], w[3 * i + 1], w[3 * i + 2]);
        let q = 3 * (n + i);
        s.vel[p] = Vec3::new(w[q], w[q + 1], w[q + 2]);
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Modified Gram-Schmidt on the tangent vectors; returns the diagonal of R.
fn orthonormalize(s: &mut State, n: usize, blocks: usize) -> Result<Vec<f64>> {
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(blocks);
    let mut diag = Vec::with_capacity(blocks);
    for b in 0..blocks {
        let mut w = tangent(s, n, b);
        for q in &basis {
            let c = dot(&w, q);
            w.iter_mut().zip(q).for_each(|(x, y)| *x -= c * y);
        }
        let norm = dot(&w, &w).sqrt();
        if !(norm > 0.0 && norm.is_finite()) {
            return invalid(format!(
                "tangent vector {b} collapsed (norm {norm}): too many exponents for the free degrees of freedom, or the run blew up"
            ));
        }
        w.iter_mut().for_each(|x| *x /= norm);
        set_tangent(s, n, b, &w);
        diag.push(norm);
        basis.push(w);
    }
    Ok(diag)
}

fn leading_norm(s: &State, n: usize) -> f64 {
    let r = n..2 * n;
    s.pos[r.clone()]
        .iter()
        .chain(&s.vel[r])
        .map(|v| v.norm_squared())
        .sum::<f64>()
        .sqrt()
}

impl World {
    /// Takes `steps` steps of size `dt` with the world's integrator while integrating the
    /// variational equations, and returns Lyapunov exponents and MEGNO. The world advances
    /// as in [`World::run`]. Not available with constraints or `wisdom_holman`.
    pub fn lyapunov(
        &mut self,
        dt: f64,
        steps: usize,
        options: &LyapunovOptions,
    ) -> Result<LyapunovRun> {
        let n = self.state.len();
        let k = options.n_exponents;
        if !self.constraints.is_empty() {
            return invalid("Lyapunov exponents are not available with constraints");
        }
        if self.integrator().name() == "wisdom_holman" {
            return invalid("Lyapunov exponents need an integrator other than wisdom_holman");
        }
        if k == 0 || k > 6 * n {
            return invalid(format!("n_exponents must be between 1 and 6 N = {}", 6 * n));
        }
        if options.renormalize_every == 0 || options.record_every == 0 {
            return invalid("renormalize_every and record_every must be at least 1");
        }
        if !(dt.is_finite() && dt != 0.0) {
            return invalid(format!("dt must be finite and non-zero, got {dt}"));
        }

        // Extended state: real particles, then one block of n pseudo-particles per vector.
        let mut ext = State::new();
        ext.t = self.state.t;
        for b in 0..=k {
            for i in 0..n {
                let (x, v) = if b == 0 {
                    (self.state.pos[i], self.state.vel[i])
                } else {
                    (Vec3::ZERO, Vec3::ZERO)
                };
                let m = if b == 0 { self.state.mass[i] } else { 1.0 };
                ext.add_particle(x, v, m);
                ext.pinned[b * n + i] = self.state.pinned[i];
            }
        }
        // Deterministic pseudo-random initial vectors (xorshift), zero on pinned particles.
        let mut rng = options.seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            (rng >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for b in 0..k {
            let mut w: Vec<f64> = (0..6 * n).map(|_| next()).collect();
            for i in 0..n {
                if self.state.pinned[i] {
                    for c in 0..3 {
                        w[3 * i + c] = 0.0;
                        w[3 * (n + i) + c] = 0.0;
                    }
                }
            }
            set_tangent(&mut ext, n, b, &w);
        }
        orthonormalize(&mut ext, n, k)?;

        let inner = Arc::new(std::mem::take(&mut self.forces));
        let mut forces = ForceSet::new();
        forces.add(Box::new(Variational {
            inner: inner.clone(),
            n,
            blocks: k,
        }));
        let result = self.integrate_tangents(&mut ext, &forces, n, k, dt, steps, options);
        drop(forces);
        self.forces = Arc::try_unwrap(inner).unwrap_or_else(|_| unreachable!("sole owner"));
        let (run, last_good) = result;
        // Leave the world at the last completed step.
        self.state.t = last_good.t;
        self.state.pos.copy_from_slice(&last_good.pos[..n]);
        self.state.vel.copy_from_slice(&last_good.vel[..n]);
        run
    }

    #[allow(clippy::too_many_arguments)]
    fn integrate_tangents(
        &mut self,
        ext: &mut State,
        forces: &ForceSet,
        n: usize,
        k: usize,
        dt: f64,
        steps: usize,
        o: &LyapunovOptions,
    ) -> (Result<LyapunovRun>, State) {
        let t0 = ext.t;
        let mut sums = vec![0.0; k];
        let mut run = LyapunovRun::default();
        let (mut megno_sum, mut mean_integral, mut megno_prev) = (0.0, 0.0, 0.0);
        let mut ln_norm = 0.0; // ln |w_1| since the last renormalisation (which sets it to 1)
        let mut good = ext.clone();
        for step in 1..=steps {
            let t_start = ext.t;
            if let Err(e) = self.integrator_mut().step(ext, forces, dt) {
                return (Err(e), good);
            }
            let finite = ext
                .pos
                .iter()
                .chain(&ext.vel)
                .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite());
            if !finite {
                return (
                    invalid(format!(
                        "state became non-finite while stepping from t = {t_start}"
                    )),
                    good,
                );
            }
            // MEGNO from the growth of the leading vector over this step.
            let ln_new = leading_norm(ext, n).ln();
            let growth = ln_new - ln_norm;
            ln_norm = ln_new;
            let elapsed = ext.t - t0;
            megno_sum += (0.5 * (t_start + ext.t) - t0) * growth;
            let y = 2.0 * megno_sum / elapsed;
            mean_integral += 0.5 * (megno_prev + y) * (ext.t - t_start);
            megno_prev = y;

            let record = step % o.record_every == 0 || step == steps;
            if step % o.renormalize_every == 0 || record {
                match orthonormalize(ext, n, k) {
                    Ok(diag) => {
                        for (s, r) in sums.iter_mut().zip(diag) {
                            *s += r.ln();
                        }
                        ln_norm = 0.0;
                    }
                    Err(e) => return (Err(e), good),
                }
            }
            good.clone_from(ext);
            if record {
                run.t.push(ext.t);
                run.running.push(sums.iter().map(|s| s / elapsed).collect());
                run.megno.push(y);
                run.mean_megno.push(mean_integral / elapsed);
            }
        }
        run.exponents = run.running.last().cloned().unwrap_or_default();
        (Ok(run), good)
    }
}
