//! Counter-based random numbers: every value is a pure function of `(seed, counter, index)`,
//! so stochastic runs are reproducible, independent of thread count and of the order in which
//! values are drawn, and restart exactly from a checkpoint that stores the seed and the counter.
//!
//! The key layout used in physim:
//! - `seed`: [`crate::World::seed`] (or a seed passed explicitly, as for the Langevin thermostat).
//! - `counter`: the step number for draws made during a run (the Langevin thermostat), or
//!   [`SETUP`] + a call number for draws made while setting a system up (random initial
//!   conditions), or `crate::nuclear::events::EVENTS` + an event number for the nuclear event
//!   generator, so these never overlap.
//! - `index`: which value within that draw, e.g. `3 × particle + component`. Each particle has
//!   its own run of indices, i.e. its own stream.
//!
//! The generator is SplitMix64's finaliser applied to the key three times; consecutive keys
//! give statistically independent values (moments, uniformity and correlations are tested).

/// Counter offset for draws made while setting a system up (e.g. random initial conditions).
pub const SETUP: u64 = 1 << 62;

/// SplitMix64 finaliser: a bijective 64-bit mixer with good avalanche.
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// Uniform deviate in (0, 1), never exactly 0 or 1.
pub fn uniform(seed: u64, counter: u64, index: u64) -> f64 {
    let h = mix(mix(mix(seed) ^ counter) ^ index);
    ((h >> 11) as f64 + 0.5) / (1u64 << 53) as f64
}

/// Standard normal deviate number `index` of draw `counter` (Box-Muller on uniforms `2 index`
/// and `2 index + 1`).
pub fn normal(seed: u64, counter: u64, index: u64) -> f64 {
    let u1 = uniform(seed, counter, 2 * index);
    let u2 = uniform(seed, counter, 2 * index + 1);
    (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normals_have_the_right_moments() {
        let n = 200_000;
        let x: Vec<f64> = (0..n)
            .map(|k| normal(42, (k / 100) as u64, (k % 100) as u64))
            .collect();
        let mean = x.iter().sum::<f64>() / n as f64;
        let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
        let kurt = x.iter().map(|v| (v - mean).powi(4)).sum::<f64>() / n as f64 / (var * var);
        assert!(mean.abs() < 0.01, "{mean}");
        assert!((var - 1.0).abs() < 0.01, "{var}");
        assert!((kurt - 3.0).abs() < 0.05, "{kurt}");
        assert_ne!(normal(1, 0, 0), normal(2, 0, 0));
        assert_eq!(normal(7, 3, 5), normal(7, 3, 5));
    }

    #[test]
    fn uniforms_are_uniform_and_uncorrelated() {
        // χ² over 100 bins for 10⁶ values along the index, and along the counter.
        let n = 1_000_000u64;
        for along_counter in [false, true] {
            let mut bins = [0u64; 100];
            let mut prev = 0.5;
            let mut lag1 = 0.0;
            for k in 0..n {
                let u = if along_counter {
                    uniform(9, k, 3)
                } else {
                    uniform(9, 3, k)
                };
                assert!(u > 0.0 && u < 1.0);
                bins[(u * 100.0) as usize] += 1;
                lag1 += (u - 0.5) * (prev - 0.5);
                prev = u;
            }
            let expect = n as f64 / 100.0;
            let chi2: f64 = bins
                .iter()
                .map(|&b| (b as f64 - expect).powi(2) / expect)
                .sum();
            // 99 degrees of freedom: mean 99, standard deviation 14; 160 is beyond 4σ.
            assert!(chi2 < 160.0, "chi2 = {chi2}");
            // Lag-1 autocorrelation of uniforms, normalised by its variance 1/12: ~N(0, 1/n).
            let rho = lag1 / n as f64 * 12.0;
            assert!(rho.abs() < 5.0 / (n as f64).sqrt(), "rho = {rho}");
        }
        // Neighbouring seeds give unrelated streams.
        let m = 100_000u64;
        let corr: f64 = (0..m)
            .map(|k| (uniform(1, 0, k) - 0.5) * (uniform(2, 0, k) - 0.5))
            .sum::<f64>()
            / m as f64
            * 12.0;
        assert!(corr.abs() < 5.0 / (m as f64).sqrt(), "{corr}");
    }
}
