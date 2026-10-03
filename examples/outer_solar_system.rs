//! The outer solar system (Sun with the inner planets' mass, Jupiter to Pluto) for a million
//! years: Wisdom-Holman against other integrators at the same step.
//! Initial conditions from Hairer, Lubich & Wanner, Geometric Numerical Integration, §I.2.4
//! (AU, days, solar masses).
//!
//! cargo run --release --example outer_solar_system [years]

use std::time::Instant;

use physim::*;

const G: f64 = 2.95912208286e-4;

fn outer_solar_system(integrator: &str) -> World {
    let bodies: [(f64, [f64; 3], [f64; 3]); 6] = [
        (1.00000597682, [0.0; 3], [0.0; 3]),
        (
            0.000954786104043,
            [-3.5023653, -3.8169847, -1.5507963],
            [0.00565429, -0.00412490, -0.00190589],
        ),
        (
            0.000285583733151,
            [9.0755314, -3.0458353, -1.6483708],
            [0.00168318, 0.00483525, 0.00192462],
        ),
        (
            0.0000437273164546,
            [8.3101420, -16.2901086, -7.2521278],
            [0.00354178, 0.00137102, 0.00055029],
        ),
        (
            0.0000517759138449,
            [11.4707666, -25.7294829, -10.8169456],
            [0.00288930, 0.00114527, 0.00039677],
        ),
        (
            1.0 / 1.3e8,
            [-15.5387357, -25.2225594, -3.1902382],
            [0.00276725, -0.00170702, -0.00136504],
        ),
    ];
    let mut w = World::with_integrator(integrator).unwrap();
    for (m, x, v) in bodies {
        w.add_particle(Vec3::from(x), Vec3::from(v), m).unwrap();
    }
    w.add_force(NewtonianGravity {
        g: G,
        softening: 0.0,
    });
    w
}

fn main() {
    let years: f64 = std::env::args()
        .nth(1)
        .map_or(1e6, |a| a.parse().expect("years"));
    let dt = 4332.59 / 20.0; // 1/20 of Jupiter's period, in days
    let steps = (years * 365.25 / dt) as usize;
    let record = (steps / 2000).max(1);
    println!("{years:e} years, dt = {dt:.1} days, {steps} steps");
    println!(
        "{:>14} {:>12} {:>12} {:>12} {:>9}",
        "integrator", "max |dE/E|", "first 10%", "last 10%", "time"
    );
    for name in ["wisdom_holman", "verlet", "yoshida4", "blanes_moan4"] {
        let mut w = outer_solar_system(name);
        let start = Instant::now();
        let result = w.run_with_options(
            dt,
            steps,
            &RunOptions {
                record_every: record,
                events: &[],
                energies: true,
            },
        );
        let elapsed = start.elapsed().as_secs_f64();
        let traj = match result {
            Ok(t) => t,
            Err(f) => {
                println!("{name:>14} failed: {}", f.error);
                continue;
            }
        };
        let e = traj.total_energy();
        let err: Vec<f64> = e.iter().map(|x| ((x - e[0]) / e[0]).abs()).collect();
        let k = err.len() / 10;
        let max = |s: &[f64]| s.iter().copied().fold(0.0, f64::max);
        println!(
            "{name:>14} {:>12.2e} {:>12.2e} {:>12.2e} {:>8.1}s",
            max(&err),
            max(&err[..k]),
            max(&err[err.len() - k..]),
            elapsed
        );
    }
}
