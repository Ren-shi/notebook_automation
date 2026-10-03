//! Wisdom-Holman mixed-variable symplectic integration for systems dominated by one
//! central body (planetary systems), in democratic heliocentric coordinates
//! (Duncan, Levison & Lee 1998).
//!
//! The Hamiltonian splits into Keplerian motion of each body about the central one (solved
//! exactly), the interactions between the other bodies (and any other position-dependent
//! forces), and a small "jump" term from the central body's recoil. The energy error scales as
//! (mass ratio) × h² instead of h², so steps of ~1/20 of the innermost orbital period are
//! enough for planetary systems.

use super::{AccelCache, Integrator};
use crate::error::{invalid, Result};
use crate::forces::{BuiltinForce, ForceSet};
use crate::state::State;
use crate::vec3::Vec3;

/// Stumpff functions `c0..c3` of `z`.
fn stumpff(z: f64) -> [f64; 4] {
    if z.abs() < 0.5 {
        // c2 = Σ (-z)^n / (2n+2)!, c3 = Σ (-z)^n / (2n+3)!, by Horner (error < 1e-19 for
        // |z| < 0.5); then c0 = 1 - z c2 and c1 = 1 - z c3.
        const INV_FACT: [f64; 24] = {
            let mut f = [1.0; 24];
            let mut k = 1;
            while k < 24 {
                f[k] = f[k - 1] / k as f64;
                k += 1;
            }
            f
        };
        let mut c2 = 0.0;
        let mut c3 = 0.0;
        let mut n = 10;
        while n > 0 {
            n -= 1;
            c2 = INV_FACT[2 * n + 2] - z * c2;
            c3 = INV_FACT[2 * n + 3] - z * c3;
        }
        [1.0 - z * c2, 1.0 - z * c3, c2, c3]
    } else if z > 0.0 {
        let s = z.sqrt();
        let (sn, cs) = s.sin_cos();
        [cs, sn / s, (1.0 - cs) / z, (s - sn) / (z * s)]
    } else {
        let s = (-z).sqrt();
        let (sh, ch) = (s.sinh(), s.cosh());
        [ch, sh / s, (ch - 1.0) / (-z), (sh - s) / (-z * s)]
    }
}

/// Advances a Kepler orbit (`gm` = G times the central mass) by time `h`, from position `r`
/// and velocity `v` relative to the central body. Universal variables with Gauss's f and g
/// functions; works for elliptic, parabolic and hyperbolic orbits. `guess` is a starting
/// value for the universal anomaly (e.g. the previous step's); the solution's is returned.
pub(crate) fn kepler_step(
    gm: f64,
    r: Vec3,
    v: Vec3,
    h: f64,
    guess: Option<f64>,
) -> Result<(Vec3, Vec3, f64)> {
    let r0 = r.norm();
    if r0 == 0.0 {
        return invalid("kepler step: body at the position of the central body");
    }
    let eta = r.dot(v);
    let beta = 2.0 * gm / r0 - v.norm_squared();
    let g = |s: f64| {
        let c = stumpff(beta * s * s);
        let (g0, g1, g2, g3) = (c[0], s * c[1], s * s * c[2], s * s * s * c[3]);
        (g0, g1, g2, g3)
    };
    // Initial guess: the caller's; else a second-order expansion of s(h), or for a long
    // elliptic step (where that expansion is poor) the mean-motion estimate ΔE/√β = β h/√GM.
    let mut s = match guess.filter(|g| g.is_finite() && g.signum() == h.signum()) {
        Some(g) => g,
        None => {
            let first = h / r0;
            let second = first - 0.5 * h * h * eta / (r0 * r0 * r0);
            if beta > 0.0 && (second - first).abs() > 0.5 * first.abs() {
                beta * h / gm.sqrt()
            } else {
                second
            }
        }
    };
    if beta > 0.0 {
        // Bound the guess by one orbit's worth of universal anomaly.
        let period_s = 2.0 * std::f64::consts::PI / beta.sqrt();
        s = s.clamp(-period_s, period_s);
    }
    let mut converged = false;
    for iter in 0..60 {
        let (g0, g1, g2, g3) = g(s);
        let f = r0 * g1 + eta * g2 + gm * g3 - h;
        let fp = r0 * g0 + eta * g1 + gm * g2; // the new radius
        let fpp = eta * g0 + (gm - beta * r0) * g1;
        let ds = if iter < 50 {
            // Laguerre-Conway step (robust far from the root), then Newton near it.
            let n = 5.0;
            let disc = ((n - 1.0) * (n - 1.0) * fp * fp - n * (n - 1.0) * f * fpp)
                .abs()
                .sqrt();
            let denom = if fp >= 0.0 { fp + disc } else { fp - disc };
            -n * f / denom
        } else {
            -f / fp
        };
        s += ds;
        // Converged when the step is at round-off, or the residual is at the round-off level
        // of its terms (the iteration can otherwise cycle between neighbouring values).
        let size = (r0 * g1).abs() + (eta * g2).abs() + (gm * g3).abs() + h.abs();
        if ds.abs() <= 4.0 * f64::EPSILON * s.abs().max(f64::MIN_POSITIVE)
            || f.abs() <= 8.0 * f64::EPSILON * size
        {
            converged = true;
            break;
        }
    }
    if !converged || !s.is_finite() {
        return invalid(format!(
            "kepler step did not converge (r = {r0}, h = {h}): orbit too close to the central body?"
        ));
    }
    let (g0, g1, g2, g3) = g(s);
    let rn = r0 * g0 + eta * g1 + gm * g2;
    let f = 1.0 - gm * g2 / r0;
    let gg = h - gm * g3;
    let fd = -gm * g1 / (rn * r0);
    let gd = 1.0 - gm * g2 / rn;
    Ok((r * f + v * gg, r * fd + v * gd, s))
}

/// Wisdom-Holman in democratic heliocentric coordinates. Particle 0 is the central body.
/// Needs exactly one `NewtonianGravity` force without softening (its `G` is used for the
/// Kepler part); other velocity-independent forces that depend only on relative positions
/// (`J2Oblateness`, springs, ...) are treated as perturbations. Second order in time, with
/// error proportional to the planet/star mass ratio. No pinned particles.
#[derive(Default)]
pub struct WisdomHolman {
    cache: AccelCache,
    q: Vec<Vec3>,
    u: Vec<Vec3>,
}

impl WisdomHolman {
    fn gravity_constant(forces: &ForceSet) -> Result<f64> {
        let mut found = None;
        for (_, f) in forces.iter() {
            if let Some(BuiltinForce::NewtonianGravity(g)) = f.builtin() {
                if found.is_some() {
                    return invalid(
                        "wisdom_holman: expected one NewtonianGravity force, found several",
                    );
                }
                if g.softening != 0.0 {
                    return invalid("wisdom_holman: NewtonianGravity must have zero softening");
                }
                found = Some(g.g);
            }
        }
        match found {
            Some(g) => Ok(g),
            None => invalid(
                "wisdom_holman needs a NewtonianGravity force (particle 0 is the central body)",
            ),
        }
    }

    /// Interaction kick for `h`: each body's acceleration minus its Kepler acceleration
    /// towards the central body.
    fn kick(&mut self, s: &State, forces: &ForceSet, gm: f64, h: f64) -> Result<()> {
        let acc = self.cache.get(s, forces)?;
        for ((u, q), a) in self.u.iter_mut().zip(&self.q).zip(acc).skip(1) {
            let r = q.norm();
            let kepler = *q * (-gm / (r * r * r));
            *u += (*a - kepler) * h;
        }
        Ok(())
    }

    /// Central-body recoil: every body moves by `h` times the mean momentum over the central mass.
    fn jump(&mut self, s: &State, h: f64) {
        let p: Vec3 = (1..s.len()).fold(Vec3::ZERO, |p, i| p + self.u[i] * s.mass[i]);
        let shift = p * (h / s.mass[0]);
        for q in &mut self.q[1..] {
            *q += shift;
        }
    }

    /// Writes the physical positions and velocities for heliocentric `q` and barycentric `u`
    /// given the barycentre `x_cm` and its velocity `v_cm`.
    fn to_state(&self, s: &mut State, x_cm: Vec3, v_cm: Vec3) {
        let total: f64 = s.mass.iter().sum();
        let n = s.len();
        let offset = (1..n).fold(Vec3::ZERO, |a, i| a + self.q[i] * s.mass[i]) / total;
        let x0 = x_cm - offset;
        let p = (1..n).fold(Vec3::ZERO, |a, i| a + self.u[i] * s.mass[i]);
        s.pos[0] = x0;
        s.vel[0] = v_cm - p / s.mass[0];
        for i in 1..n {
            s.pos[i] = x0 + self.q[i];
            s.vel[i] = v_cm + self.u[i];
        }
    }
}

impl Integrator for WisdomHolman {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        let n = s.len();
        if n < 2 {
            return Ok(());
        }
        if s.pinned.iter().any(|&p| p) {
            return invalid("wisdom_holman does not support pinned particles");
        }
        if forces.velocity_dependent() {
            return invalid("wisdom_holman needs velocity-independent forces");
        }
        let g = Self::gravity_constant(forces)?;
        let m0 = s.mass[0];
        if m0.is_nan() || m0 <= 0.0 {
            return invalid("wisdom_holman: the central body (particle 0) must have positive mass");
        }
        let gm = g * m0;
        let total: f64 = s.mass.iter().sum();
        let x_cm = s
            .pos
            .iter()
            .zip(&s.mass)
            .fold(Vec3::ZERO, |a, (x, m)| a + *x * *m)
            / total;
        let v_cm = s
            .vel
            .iter()
            .zip(&s.mass)
            .fold(Vec3::ZERO, |a, (v, m)| a + *v * *m)
            / total;
        self.q.clear();
        self.u.clear();
        for i in 0..n {
            self.q.push(s.pos[i] - s.pos[0]);
            self.u.push(s.vel[i] - v_cm);
        }

        self.kick(s, forces, gm, 0.5 * dt)?;
        self.jump(s, 0.5 * dt);
        // No state carried between steps (e.g. Kepler guesses), so restarts are bit-exact.
        for i in 1..n {
            let (q, u, _) = kepler_step(gm, self.q[i], self.u[i], dt, None)?;
            self.q[i] = q;
            self.u[i] = u;
        }
        self.jump(s, 0.5 * dt);
        // The final kick needs the forces at the new positions: write them out first.
        let t0 = s.t;
        s.t = t0 + dt;
        let x_cm_new = x_cm + v_cm * dt;
        self.to_state(s, x_cm_new, v_cm);
        self.kick(s, forces, gm, 0.5 * dt)?;
        self.to_state(s, x_cm_new, v_cm);
        Ok(())
    }
    fn name(&self) -> &str {
        "wisdom_holman"
    }
    fn order(&self) -> u32 {
        2
    }
    fn symplectic(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stumpff_series_matches_closed_forms() {
        for z in [-0.49, -0.3, -1e-3, 1e-8, 0.2, 0.4999] {
            let c = stumpff(z);
            let exact = if z > 0.0 {
                let s: f64 = z.sqrt();
                [
                    s.cos(),
                    s.sin() / s,
                    (1.0 - s.cos()) / z,
                    (s - s.sin()) / (z * s),
                ]
            } else {
                let s: f64 = (-z).sqrt();
                [
                    s.cosh(),
                    s.sinh() / s,
                    (s.cosh() - 1.0) / -z,
                    (s.sinh() - s) / (-z * s),
                ]
            };
            for k in 0..4 {
                // Closed forms lose digits to cancellation at small |z|; compare loosely there.
                let tol = if z.abs() < 0.01 { 1e-6 } else { 1e-14 };
                assert!(
                    (c[k] - exact[k]).abs() < tol * exact[k].abs(),
                    "z = {z}, k = {k}"
                );
            }
        }
    }

    #[test]
    fn kepler_step_returns_to_the_start_after_one_period() {
        for (gm, a, e) in [(1.0, 1.0, 0.0), (2.0, 1.5, 0.6), (0.5, 3.0, 0.95)] {
            let r = Vec3::new(a * (1.0 - e), 0.0, 0.0);
            let v = Vec3::new(0.0, (gm * (1.0 + e) / (a * (1.0 - e))).sqrt(), 0.1 * 0.0);
            let period = 2.0 * std::f64::consts::PI * (a * a * a / gm).sqrt();
            let (mut x, mut u) = (r, v);
            let mut guess = None;
            for _ in 0..7 {
                let (x1, u1, s) = kepler_step(gm, x, u, period / 7.0, guess).unwrap();
                (x, u, guess) = (x1, u1, Some(s));
            }
            assert!((x - r).norm() < 1e-11 * a, "{gm} {a} {e}: {:?}", x - r);
            assert!((u - v).norm() < 1e-11 * v.norm());
        }
        // Hyperbolic: energy and angular momentum are conserved.
        let (gm, r, v) = (1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.3, 1.8, 0.2));
        let (x, u, _) = kepler_step(gm, r, v, 5.0, None).unwrap();
        let energy = |x: Vec3, u: Vec3| 0.5 * u.norm_squared() - gm / x.norm();
        assert!((energy(x, u) - energy(r, v)).abs() < 1e-13);
        assert!((x.cross(u) - r.cross(v)).norm() < 1e-13);
    }
}
