//! The Boris pusher for charged particles in electric and magnetic fields.

use super::Integrator;
use crate::error::{invalid, Result};
use crate::forces::ForceSet;
use crate::state::State;
use crate::vec3::Vec3;

/// Boris (1970): drift half a step, then a half electric kick, an exact-norm rotation about
/// the magnetic field, another half electric kick, and the second half drift. Second order,
/// time-reversible and volume-preserving; in a pure magnetic field it conserves |v| (hence
/// kinetic energy) to round-off for any step. The gyration phase advances by
/// `2 atan(ω h / 2)` per step instead of `ω h`, a relative frequency error of `(ω h)²/12`.
///
/// Magnetic forces (`MagneticField`, `FieldFunctions` with `B`) are rotated; every other
/// force is applied as a kick and must not depend on velocity. One force evaluation per step.
#[derive(Default)]
pub struct Boris {
    acc: Vec<Vec3>,
    b: Vec<Vec3>,
}

impl Integrator for Boris {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        let h = 0.5 * dt;
        for (x, v) in s.pos.iter_mut().zip(&s.vel) {
            *x += *v * h;
        }
        let t_half = s.t + h;
        forces.electric_and_magnetic(
            t_half,
            &s.pos,
            &s.vel,
            &s.mass,
            &s.charge,
            &s.radius,
            &s.pinned,
            &mut self.acc,
            &mut self.b,
        )?;
        for i in 0..s.len() {
            if s.pinned[i] {
                continue;
            }
            let (q, m) = (s.charge[i], s.mass[i]);
            let mut v = s.vel[i] + self.acc[i] * h;
            if q != 0.0 && self.b[i] != Vec3::ZERO {
                if m == 0.0 {
                    return invalid(format!("boris: particle {i} is charged but massless"));
                }
                let t = self.b[i] * (q / m * h);
                let s_vec = t * (2.0 / (1.0 + t.norm_squared()));
                let v_prime = v + v.cross(t);
                v += v_prime.cross(s_vec);
            }
            s.vel[i] = v + self.acc[i] * h;
        }
        for (x, v) in s.pos.iter_mut().zip(&s.vel) {
            *x += *v * h;
        }
        s.t += dt;
        Ok(())
    }
    fn name(&self) -> &str {
        "boris"
    }
    fn order(&self) -> u32 {
        2
    }
    fn symplectic(&self) -> bool {
        false
    }
}
