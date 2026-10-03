//! Counter-based random numbers: every value is a pure function of (seed, counter, index),
//! so stochastic runs are reproducible, independent of thread count, and restart exactly
//! from a checkpoint that stores the seed and the counter.

/// SplitMix64 finaliser: a bijective 64-bit mixer with good avalanche.
fn mix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// Uniform in (0, 1).
fn uniform(seed: u64, counter: u64, index: u64) -> f64 {
    let h = mix(mix(mix(seed) ^ counter) ^ index);
    ((h >> 11) as f64 + 0.5) / (1u64 << 53) as f64
}

/// Standard normal deviate number `index` of draw `counter` (Box-Muller).
pub(crate) fn normal(seed: u64, counter: u64, index: u64) -> f64 {
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
}
