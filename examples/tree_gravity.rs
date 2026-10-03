//! Barnes-Hut tree gravity against direct summation on a Plummer sphere: force error as a
//! function of the opening angle θ, and the time per force evaluation.
//!
//! cargo run --release --example tree_gravity [N]

use std::time::Instant;

use physim::*;

/// `n` equal-mass particles from a Plummer sphere (scale radius 1, total mass 1).
fn plummer(n: usize, seed: u64) -> (Vec<Vec3>, Vec<f64>) {
    let mut s = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
    let mut rnd = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut pos = Vec::with_capacity(n);
    while pos.len() < n {
        let x: f64 = rnd();
        if x < 1e-9 {
            continue;
        }
        let r = 1.0 / (x.powf(-2.0 / 3.0) - 1.0).sqrt();
        if r > 30.0 {
            continue;
        }
        let z = 2.0 * rnd() - 1.0;
        let phi = std::f64::consts::TAU * rnd();
        let st = (1.0 - z * z).sqrt();
        pos.push(Vec3::new(r * st * phi.cos(), r * st * phi.sin(), r * z));
    }
    (pos, vec![1.0 / n as f64; n])
}

fn accelerations(f: &dyn Force, pos: &[Vec3], mass: &[f64]) -> (Vec<Vec3>, f64) {
    let vel = vec![Vec3::ZERO; pos.len()];
    let mut acc = vec![Vec3::ZERO; pos.len()];
    let start = Instant::now();
    f.accumulate(0.0, pos, &vel, mass, &mut acc).unwrap();
    (acc, start.elapsed().as_secs_f64())
}

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .map_or(100_000, |a| a.parse().expect("N"));
    let (pos, mass) = plummer(n, 3);
    let (direct, t_direct) = accelerations(
        &NewtonianGravity {
            g: 1.0,
            softening: 0.0,
        },
        &pos,
        &mass,
    );
    println!("N = {n}: direct summation {:.0} ms\n", t_direct * 1e3);
    println!(
        "{:>10} {:>6} {:>12} {:>12} {:>10} {:>9}",
        "terms", "θ", "rms error", "99% error", "time", "speedup"
    );
    for quadrupole in [false, true] {
        for theta in [0.2, 0.3, 0.5, 0.7, 1.0] {
            let tree = TreeGravity {
                g: 1.0,
                softening: 0.0,
                theta,
                quadrupole,
            };
            let (acc, t) = accelerations(&tree, &pos, &mass);
            let mut err: Vec<f64> = acc
                .iter()
                .zip(&direct)
                .map(|(a, b)| (*a - *b).norm() / b.norm())
                .collect();
            let rms = (err.iter().map(|e| e * e).sum::<f64>() / n as f64).sqrt();
            err.sort_by(f64::total_cmp);
            println!(
                "{:>10} {theta:>6} {rms:>12.2e} {:>12.2e} {:>8.0}ms {:>8.0}x",
                if quadrupole { "quadrupole" } else { "monopole" },
                err[err.len() * 99 / 100],
                t * 1e3,
                t_direct / t
            );
        }
    }
}
