//! Particle-mesh gravity against Barnes-Hut and direct summation on a Plummer sphere:
//! accuracy (against an exact direct sum for a sample of particles) and time per evaluation.
//!
//! Run with `cargo run --release --example particle_mesh`.

use physim::*;
use std::time::Instant;

fn main() -> Result<()> {
    let n = 100_000;
    let mut s = 12345u64;
    let mut u = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    // Plummer sphere (a = 1), truncated at r = 4.
    let pos: Vec<Vec3> = (0..n)
        .map(|_| {
            let r = 1.0 / ((u() * 0.85).powf(-2.0 / 3.0) - 1.0).sqrt();
            let z = 2.0 * u() - 1.0;
            let phi = 2.0 * std::f64::consts::PI * u();
            let q = (1.0 - z * z).sqrt();
            Vec3::new(q * phi.cos(), q * phi.sin(), z) * r.min(4.0)
        })
        .collect();
    let mass = vec![1.0 / n as f64; n];
    let vel = vec![Vec3::ZERO; n];
    let sample: Vec<usize> = (0..n).step_by(n / 1000).collect();

    println!("N = {n} Plummer sphere; errors against direct summation with the same softening");
    println!(
        "{:<34} {:>10} {:>12} {:>12}",
        "method", "time", "median err", "90% err"
    );
    let report = |name: &str, eps: f64, f: &dyn Force| -> Result<()> {
        let mut acc = vec![Vec3::ZERO; n];
        f.accumulate(0.0, &pos, &vel, &mass, &mut acc)?; // warm-up (kernels, trees)
        let t0 = Instant::now();
        acc.iter_mut().for_each(|a| *a = Vec3::ZERO);
        f.accumulate(0.0, &pos, &vel, &mass, &mut acc)?;
        let dt = t0.elapsed();
        let mut errs: Vec<f64> = sample
            .iter()
            .map(|&i| {
                let mut exact = Vec3::ZERO;
                for j in 0..n {
                    let d = pos[j] - pos[i];
                    let s2 = d.norm_squared() + eps * eps;
                    exact += d * (mass[j] / (s2 * s2.sqrt()));
                }
                (acc[i] - exact).norm() / exact.norm()
            })
            .collect();
        errs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "{name:<34} {:>8.1} ms {:>12.2e} {:>12.2e}",
            dt.as_secs_f64() * 1e3,
            errs[errs.len() / 2],
            errs[errs.len() * 9 / 10]
        );
        Ok(())
    };
    for cells in [32usize, 64] {
        let pm = ParticleMesh::new(1.0, cells, 10.0, Vec3::ZERO, false, None);
        let eps = pm.spacing();
        report(
            &format!("particle mesh {cells}³ (eps = h = {eps:.3})"),
            eps,
            &pm,
        )?;
        let tree = TreeGravity {
            g: 1.0,
            softening: eps,
            theta: 0.5,
            quadrupole: false,
        };
        report(
            &format!("Barnes-Hut theta 0.5 (eps = {eps:.3})"),
            eps,
            &tree,
        )?;
    }
    Ok(())
}
