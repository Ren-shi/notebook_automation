//! Thermostats: integrators that hold the system at a temperature (k_B = 1).

use super::{AccelCache, Integrator, Scheme};
use crate::error::{invalid, Result};
use crate::forces::ForceSet;
use crate::rng;
use crate::state::State;

fn kick(s: &mut State, acc: &[crate::vec3::Vec3], h: f64) {
    for (v, a) in s.vel.iter_mut().zip(acc) {
        *v += *a * h;
    }
}

fn drift(s: &mut State, h: f64) {
    for (x, v) in s.pos.iter_mut().zip(&s.vel) {
        *x += *v * h;
    }
    s.t += h;
}

/// Langevin dynamics `dv = a dt − γ v dt + √(2 γ T / m) dW` by the BAOAB splitting
/// (Leimkuhler & Matthews): half kick, half drift, exact Ornstein-Uhlenbeck velocity update,
/// half drift, half kick. Samples the canonical distribution with configurational errors of
/// order h²; one force evaluation per step. Random numbers come from a counter-based
/// generator, so runs are reproducible and restarts from checkpoints are exact.
pub struct Langevin {
    pub temperature: f64,
    /// Friction rate γ (1/time).
    pub friction: f64,
    pub seed: u64,
    /// Steps taken (the random-number counter).
    pub counter: u64,
    cache: AccelCache,
}

impl Langevin {
    pub fn new(temperature: f64, friction: f64, seed: u64) -> Result<Self> {
        if !(temperature >= 0.0
            && temperature.is_finite()
            && friction >= 0.0
            && friction.is_finite())
        {
            return invalid(format!(
                "Langevin needs a non-negative temperature and friction, got {temperature}, {friction}"
            ));
        }
        Ok(Self {
            temperature,
            friction,
            seed,
            counter: 0,
            cache: AccelCache::default(),
        })
    }
}

impl Integrator for Langevin {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        let h = 0.5 * dt;
        let acc = self.cache.get(s, forces)?;
        kick(s, acc, h);
        drift(s, h);
        let c1 = (-self.friction * dt.abs()).exp();
        let c2 = (1.0 - c1 * c1).sqrt();
        for i in 0..s.len() {
            if s.pinned[i] {
                continue;
            }
            let m = s.mass[i];
            if m == 0.0 {
                return invalid(format!("Langevin: particle {i} is massless"));
            }
            let sigma = c2 * (self.temperature / m).sqrt();
            let base = 3 * i as u64;
            let xi = crate::vec3::Vec3::new(
                rng::normal(self.seed, self.counter, base),
                rng::normal(self.seed, self.counter, base + 1),
                rng::normal(self.seed, self.counter, base + 2),
            );
            s.vel[i] = s.vel[i] * c1 + xi * sigma;
        }
        self.counter += 1;
        drift(s, h);
        let acc = self.cache.get(s, forces)?;
        kick(s, acc, h);
        Ok(())
    }
    fn name(&self) -> &str {
        "langevin"
    }
    fn order(&self) -> u32 {
        2
    }
    fn symplectic(&self) -> bool {
        false
    }
    fn scheme(&self) -> Option<Scheme> {
        Some(Scheme::Langevin {
            temperature: self.temperature,
            friction: self.friction,
            seed: self.seed,
            counter: self.counter,
        })
    }
}

/// Nosé-Hoover thermostat: `dv/dt = a − ξ v`, `dξ/dt = (Σ m v² − g T) / Q` with
/// `Q = g T τ²` (τ: thermostat relaxation time, g = 3 × free particles). Integrated by the
/// symmetric Trotter splitting of Martyna et al. (thermostat half step, velocity Verlet,
/// thermostat half step). Deterministic and time-reversible; it conserves the extended
/// energy `K + U + Q ξ²/2 + g T η` (`η = ∫ ξ dt`), reported by
/// [`Integrator::thermostat_energy`].
pub struct NoseHoover {
    pub temperature: f64,
    pub tau: f64,
    pub xi: f64,
    pub eta: f64,
    cache: AccelCache,
}

impl NoseHoover {
    pub fn new(temperature: f64, tau: f64) -> Result<Self> {
        if !(temperature > 0.0 && temperature.is_finite() && tau > 0.0 && tau.is_finite()) {
            return invalid(format!(
                "Nosé-Hoover needs a positive temperature and tau, got {temperature}, {tau}"
            ));
        }
        Ok(Self {
            temperature,
            tau,
            xi: 0.0,
            eta: 0.0,
            cache: AccelCache::default(),
        })
    }

    fn free(s: &State) -> Result<(f64, f64)> {
        let (mut twice_k, mut n) = (0.0, 0usize);
        for i in 0..s.len() {
            if s.pinned[i] {
                continue;
            }
            if s.mass[i] == 0.0 {
                return invalid(format!("Nosé-Hoover: particle {i} is massless"));
            }
            twice_k += s.mass[i] * s.vel[i].norm_squared();
            n += 1;
        }
        Ok((twice_k, 3.0 * n as f64))
    }

    /// Thermostat flow for time `h`.
    fn thermostat(&mut self, s: &mut State, h: f64) -> Result<()> {
        let (mut twice_k, g) = Self::free(s)?;
        if g == 0.0 {
            return Ok(());
        }
        let q = g * self.temperature * self.tau * self.tau;
        self.xi += 0.5 * h * (twice_k - g * self.temperature) / q;
        let scale = (-self.xi * h).exp();
        for i in 0..s.len() {
            if !s.pinned[i] {
                s.vel[i] *= scale;
            }
        }
        twice_k *= scale * scale;
        self.eta += self.xi * h;
        self.xi += 0.5 * h * (twice_k - g * self.temperature) / q;
        Ok(())
    }
}

impl Integrator for NoseHoover {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        let h = 0.5 * dt;
        self.thermostat(s, h)?;
        let acc = self.cache.get(s, forces)?;
        kick(s, acc, h);
        drift(s, dt);
        let acc = self.cache.get(s, forces)?;
        kick(s, acc, h);
        self.thermostat(s, h)
    }
    fn name(&self) -> &str {
        "nose_hoover"
    }
    fn order(&self) -> u32 {
        2
    }
    fn symplectic(&self) -> bool {
        false
    }
    fn scheme(&self) -> Option<Scheme> {
        Some(Scheme::NoseHoover {
            temperature: self.temperature,
            tau: self.tau,
            xi: self.xi,
            eta: self.eta,
        })
    }
    fn thermostat_energy(&self, s: &State) -> f64 {
        let g = 3.0 * s.pinned.iter().filter(|p| !**p).count() as f64;
        let q = g * self.temperature * self.tau * self.tau;
        0.5 * q * self.xi * self.xi + g * self.temperature * self.eta
    }
}
