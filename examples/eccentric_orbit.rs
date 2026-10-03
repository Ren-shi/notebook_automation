//! Force evaluations needed to integrate an e = 0.99 Kepler orbit to a given energy error,
//! adaptive Dormand-Prince versus every fixed-step integrator.
//!
//! cargo run --release --example eccentric_orbit

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use physim::*;

const E: f64 = 0.99;
const ORBITS: f64 = 3.0;

/// Gravity that counts its evaluations.
struct Counted {
    inner: NewtonianGravity,
    count: Arc<AtomicU64>,
}

impl Force for Counted {
    fn accumulate(&self, t: f64, p: &[Vec3], v: &[Vec3], m: &[f64], a: &mut [Vec3]) -> Result<()> {
        self.count.fetch_add(1, Ordering::Relaxed);
        self.inner.accumulate(t, p, v, m, a)
    }
    fn potential(&self, t: f64, p: &[Vec3], m: &[f64]) -> Result<Option<f64>> {
        self.inner.potential(t, p, m)
    }
    fn name(&self) -> String {
        "CountedGravity".into()
    }
}

/// Unit-mass sun pinned at the origin, light planet at apoapsis of an a = 1 orbit (period 2π).
fn orbit(integrator: &str) -> (World, Arc<AtomicU64>) {
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.pin(0, true).unwrap();
    let v = ((1.0 - E) / (1.0 + E)).sqrt();
    w.add_particle(Vec3::new(1.0 + E, 0.0, 0.0), Vec3::new(0.0, v, 0.0), 1e-6)
        .unwrap();
    let count = Arc::new(AtomicU64::new(0));
    w.add_force(Counted {
        inner: NewtonianGravity {
            g: 1.0,
            softening: 0.0,
        },
        count: count.clone(),
    });
    (w, count)
}

fn energy_error(traj: &Trajectory) -> f64 {
    let e = traj.total_energy();
    e.iter()
        .map(|x| ((x - e[0]) / e[0]).abs())
        .fold(0.0, f64::max)
}

/// Max energy error and force evaluations for a fixed-step run with `steps` steps.
fn fixed(integrator: &str, steps: usize) -> (f64, u64) {
    let (mut w, count) = orbit(integrator);
    let t_end = ORBITS * std::f64::consts::TAU;
    match w.run(t_end / steps as f64, steps, 1) {
        Ok(traj) => (energy_error(&traj), count.load(Ordering::Relaxed)),
        Err(_) => (f64::INFINITY, count.load(Ordering::Relaxed)),
    }
}

fn main() {
    let t_end = ORBITS * std::f64::consts::TAU;
    println!("e = {E}, {ORBITS} orbits; max relative energy error over the run\n");
    println!(
        "{:>10} {:>12} {:>10} {:>9}",
        "rtol", "energy err", "evals", "rejected"
    );
    let mut targets = Vec::new();
    for rtol in [1e-6, 1e-8, 1e-10, 1e-12] {
        let (mut w, count) = orbit("verlet");
        let run = w
            .run_adaptive(
                t_end,
                &AdaptiveOptions {
                    rtol,
                    atol: rtol * 1e-3,
                    ..Default::default()
                },
            )
            .unwrap();
        let err = energy_error(&run.trajectory);
        let evals = count.load(Ordering::Relaxed);
        assert_eq!(evals, run.stats.evaluations);
        println!(
            "{rtol:>10.0e} {err:>12.2e} {evals:>10} {:>9}",
            run.stats.rejected
        );
        targets.push((err, evals));
    }

    println!("\nFixed steps needed to match each adaptive error (doubling search):");
    print!("{:>12}", "target");
    let names = ["verlet", "yoshida4", "rk4", "dopri5"];
    for n in names {
        print!(" {n:>14}");
    }
    println!();
    for &(target, adaptive_evals) in &targets {
        print!("{target:>12.2e}");
        for name in names {
            let mut steps = 1000;
            let found = loop {
                let (err, evals) = fixed(name, steps);
                if err <= target {
                    break Some(evals);
                }
                if steps > 1 << 25 {
                    break None;
                }
                steps *= 2;
            };
            match found {
                Some(evals) => print!(
                    " {:>14}",
                    format!("{evals} ({:.0}x)", evals as f64 / adaptive_evals as f64)
                ),
                None => print!(" {:>14}", "> 3e7 steps"),
            }
        }
        println!();
    }
}
