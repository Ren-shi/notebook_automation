//! Relativistic two-body kinematics for the per-event loop: a port of
//! `physim.nuclear.kinematics.TwoBody` (Python), checked against it in the Python tests.
//!
//! The beam (mass `m1`, kinetic energy `t`) hits a target (`m2`) at rest; the ejectile (`m3`)
//! leaves at CM angle θ* and the recoil (`m4`) at 180° − θ*. Masses and energies in MeV.

/// Lab angle (radians) and kinetic energy (MeV) of one outgoing particle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LabPoint {
    pub theta: f64,
    pub energy: f64,
}

/// One reaction at one beam energy.
#[derive(Clone, Copy, Debug)]
pub struct TwoBody {
    pub m: [f64; 4],
    /// CM velocity (units of c) and Lorentz factor.
    pub beta_cm: f64,
    pub gamma_cm: f64,
    /// CM momentum of the outgoing pair, MeV/c.
    pub p_cm: f64,
    /// Total CM energies of the ejectile and the recoil, MeV.
    pub e_cm: [f64; 2],
    /// √s, MeV.
    pub sqrt_s: f64,
}

impl TwoBody {
    /// `None` below threshold.
    pub fn new(m1: f64, m2: f64, m3: f64, m4: f64, t: f64) -> Option<Self> {
        let p1 = (t * (t + 2.0 * m1)).sqrt();
        let e_total = t + m1 + m2;
        let s = (m1 + m2).powi(2) + 2.0 * m2 * t;
        let sqrt_s = s.sqrt();
        if t <= 0.0 || sqrt_s < m3 + m4 {
            return None;
        }
        // s − (m3 ± m4)², written so nothing large cancels (as in the Python version).
        let q = (m1 + m2) - (m3 + m4);
        let plus = q * (m1 + m2 + m3 + m4) + 2.0 * m2 * t;
        let minus = (m1 + m2 - m3 + m4) * (m1 + m2 + m3 - m4) + 2.0 * m2 * t;
        let p_cm = (plus * minus).max(0.0).sqrt() / (2.0 * sqrt_s);
        Some(TwoBody {
            m: [m1, m2, m3, m4],
            beta_cm: p1 / e_total,
            gamma_cm: e_total / sqrt_s,
            p_cm,
            e_cm: [p_cm.hypot(m3), p_cm.hypot(m4)],
            sqrt_s,
        })
    }

    /// Kinetic energy available in the CM frame before the reaction, MeV.
    pub fn cm_energy(&self) -> f64 {
        self.sqrt_s - self.m[0] - self.m[1]
    }

    /// The ejectile (`recoil = false`) or recoil leaving at CM angle `theta_cm` (radians).
    pub fn at_cm(&self, theta_cm: f64, recoil: bool) -> LabPoint {
        let k = usize::from(recoil);
        let (e, m) = (self.e_cm[k], self.m[2 + k]);
        let (s, c) = theta_cm.sin_cos();
        let g = self.beta_cm * e / self.p_cm;
        let p_par = self.gamma_cm * self.p_cm * (c + g);
        let p_perp = self.p_cm * s;
        LabPoint {
            theta: p_perp.atan2(p_par),
            energy: self.gamma_cm * (e + self.beta_cm * self.p_cm * c) - m,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn energy_and_momentum_are_conserved() {
        // 16O on 208Pb at 64 MeV and 208Pb on 1H at 1 GeV (inverse kinematics), with Q ≠ 0.
        for &(m1, m2, m3, m4, t) in &[
            (14895.08, 193688.0, 14895.08, 193688.0, 64.0),
            (193688.0, 938.272, 193688.0, 938.272 + 2.0, 1000.0),
            (1875.6, 2808.9, 938.272, 3727.4, 10.0),
        ] {
            let r = TwoBody::new(m1, m2, m3, m4, t).unwrap();
            for k in 0..=36 {
                let th = (k as f64 * 5.0).to_radians();
                let a = r.at_cm(th, false);
                let b = r.at_cm(std::f64::consts::PI - th, true);
                let p = |e: f64, m: f64| (e * (e + 2.0 * m)).sqrt();
                let (pa, pb) = (p(a.energy, m3), p(b.energy, m4));
                let p1 = p(t, m1);
                assert!((a.energy + b.energy + m3 + m4 - t - m1 - m2).abs() < 1e-7);
                let pz = pa * a.theta.cos() + pb * b.theta.cos();
                let px = pa * a.theta.sin() - pb * b.theta.sin();
                assert!((pz - p1).abs() < 1e-7 * p1, "{pz} {p1}");
                assert!(px.abs() < 1e-7 * p1, "{px}");
            }
        }
        assert!(TwoBody::new(938.0, 2808.0, 939.0, 2808.0, 0.1).is_none());
    }
}
