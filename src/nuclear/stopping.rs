//! Energy loss through layers from tables that the Python side computes
//! (`physim.nuclear.stopping.Stopping.transport_table`).
//!
//! A table holds, on a grid of kinetic energies E: the CSDA range R(E) (mg/cm²), the stopping
//! power S(E) (MeV per mg/cm²) and the straggling integral W(E) = ∫ (dΩ²/dx) / S³ dE. With
//! these, crossing a path L is
//!
//! E_out = R⁻¹(R(E_in) − L),   σ²_added = S(E_out)² (W(E_in) − W(E_out)),
//!
//! which is the solution of the straggling equation dσ²/dx = −2 S′ σ² + dΩ²/dx: a spread
//! already present is carried by the mapping E_in → E_out itself, and the Bohr straggling added
//! along the path is carried to the exit by the factor S(E_out)² / S(E)².

/// Range, stopping power and straggling integral of one ion in one material.
#[derive(Clone, Debug)]
pub struct Table {
    ln_e: Vec<f64>,
    ln_r: Vec<f64>,
    ln_s: Vec<f64>,
    w: Vec<f64>,
}

/// Index `i` with `xs[i] <= x < xs[i + 1]`, clamped to the table, and the fraction along it.
fn locate(xs: &[f64], x: f64) -> (usize, f64) {
    let n = xs.len();
    if x <= xs[0] {
        return (0, 0.0);
    }
    if x >= xs[n - 1] {
        return (n - 2, 1.0);
    }
    let i = xs.partition_point(|&v| v <= x) - 1;
    let i = i.min(n - 2);
    (i, (x - xs[i]) / (xs[i + 1] - xs[i]))
}

fn lerp(ys: &[f64], (i, f): (usize, f64)) -> f64 {
    ys[i] + f * (ys[i + 1] - ys[i])
}

impl Table {
    /// From energies (MeV, increasing), ranges (mg/cm², increasing), stopping powers
    /// (MeV/(mg/cm²)) and the straggling integral W (see the module docs), all on one grid.
    pub fn new(energy: &[f64], range: &[f64], stopping: &[f64], w: &[f64]) -> Result<Self, String> {
        let n = energy.len();
        if n < 2 || range.len() != n || stopping.len() != n || w.len() != n {
            return Err("transport table arrays must have the same length, at least 2".into());
        }
        let increasing = |v: &[f64]| v.windows(2).all(|p| p[1] > p[0]) && v[0] > 0.0;
        if !increasing(energy) || !increasing(range) {
            return Err(
                "transport table energies and ranges must be positive and increasing".into(),
            );
        }
        if stopping
            .iter()
            .any(|&s| s.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater))
        {
            return Err("transport table stopping powers must be positive".into());
        }
        Ok(Table {
            ln_e: energy.iter().map(|e| e.ln()).collect(),
            ln_r: range.iter().map(|r| r.ln()).collect(),
            ln_s: stopping.iter().map(|s| s.ln()).collect(),
            w: w.to_vec(),
        })
    }

    /// CSDA range from `energy` (MeV), mg/cm²; below the table, extrapolated in log-log from the
    /// first interval.
    pub fn range(&self, energy: f64) -> f64 {
        if energy <= 0.0 {
            return 0.0;
        }
        let x = energy.ln();
        if x < self.ln_e[0] {
            let slope = (self.ln_r[1] - self.ln_r[0]) / (self.ln_e[1] - self.ln_e[0]);
            return (self.ln_r[0] + slope * (x - self.ln_e[0])).exp();
        }
        lerp(&self.ln_r, locate(&self.ln_e, x)).exp()
    }

    /// Energy (MeV) with range `range` (mg/cm²).
    pub fn energy_for_range(&self, range: f64) -> f64 {
        if range <= 0.0 {
            return 0.0;
        }
        let x = range.ln();
        if x < self.ln_r[0] {
            let slope = (self.ln_e[1] - self.ln_e[0]) / (self.ln_r[1] - self.ln_r[0]);
            return (self.ln_e[0] + slope * (x - self.ln_r[0])).exp();
        }
        lerp(&self.ln_e, locate(&self.ln_r, x)).exp()
    }

    /// Stopping power at `energy`, MeV/(mg/cm²).
    pub fn stopping(&self, energy: f64) -> f64 {
        lerp(&self.ln_s, locate(&self.ln_e, energy.max(1e-300).ln())).exp()
    }

    fn w_at(&self, energy: f64) -> f64 {
        lerp(&self.w, locate(&self.ln_e, energy.max(1e-300).ln()))
    }

    /// Mean energy after a path `path` (mg/cm²); 0 if the ion stops.
    pub fn energy_after(&self, energy: f64, path: f64) -> f64 {
        let left = self.range(energy) - path;
        if left > 0.0 {
            self.energy_for_range(left)
        } else {
            0.0
        }
    }

    /// Energy after `path` and the standard deviation of the straggling added on the way, MeV.
    pub fn cross(&self, energy: f64, path: f64) -> (f64, f64) {
        let out = self.energy_after(energy, path);
        if out <= 0.0 {
            return (0.0, 0.0);
        }
        let var = self.stopping(out).powi(2) * (self.w_at(energy) - self.w_at(out));
        (out, var.max(0.0).sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A power law, S = c E^(−½), has R = (2/3) E^(3/2) / c, which log-log interpolation
    /// reproduces exactly.
    #[test]
    fn power_law_table_is_exact() {
        let c = 0.5;
        let e: Vec<f64> = (0..400).map(|k| 1e-3 * 1.03f64.powi(k)).collect();
        let r: Vec<f64> = e.iter().map(|e| 2.0 / 3.0 * e.powf(1.5) / c).collect();
        let s: Vec<f64> = e.iter().map(|e| c / e.sqrt()).collect();
        let w = vec![0.0; e.len()];
        let t = Table::new(&e, &r, &s, &w).unwrap();
        for &e0 in &[0.01, 0.5, 3.0, 40.0] {
            let path = 0.3 * t.range(e0);
            let exact = (e0.powf(1.5) - 1.5 * c * path).powf(2.0 / 3.0);
            assert!((t.energy_after(e0, path) / exact - 1.0).abs() < 1e-9);
        }
        assert_eq!(t.energy_after(1.0, 10.0), 0.0);
        assert!(Table::new(&e, &r, &s, &w[1..]).is_err());
    }
}
