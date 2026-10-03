//! Lennard-Jones liquid near the triple point (ρ = 0.8442, T ≈ 0.72, rc = 2.5σ shifted):
//! NVE energy conservation, Langevin and Nosé-Hoover thermostats, and g(r).
//!
//! cargo run --release --example lj_liquid [steps]

use std::time::Instant;

use physim::*;

/// N = 4 k³ atoms on an fcc lattice filling a cubic box at density `rho`, with velocities of
/// temperature `t` and zero total momentum. Returns the world and the box side.
fn lj_liquid(k: usize, rho: f64, t: f64, integrator: &str) -> (World, f64) {
    let n = 4 * k * k * k;
    let side = (n as f64 / rho).cbrt();
    let a = side / k as f64;
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut w = World::with_integrator(integrator).unwrap();
    let mut s = 12345u64;
    let mut gauss = || {
        let mut u = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        };
        let (u1, u2) = (u(), u());
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    };
    for x in 0..k {
        for y in 0..k {
            for z in 0..k {
                for b in basis {
                    let p = Vec3::new(
                        (x as f64 + b[0]) * a,
                        (y as f64 + b[1]) * a,
                        (z as f64 + b[2]) * a,
                    );
                    let v = Vec3::new(gauss(), gauss(), gauss()) * t.sqrt();
                    w.add_particle(p, v, 1.0).unwrap();
                }
            }
        }
    }
    let p = w.state.momentum() / n as f64;
    let vel: Vec<Vec3> = w.state.vel.iter().map(|v| *v - p).collect();
    w.set_velocities(vel).unwrap();
    let scale = (t / w.temperature()).sqrt();
    let vel: Vec<Vec3> = w.state.vel.iter().map(|v| *v * scale).collect();
    w.set_velocities(vel).unwrap();
    w.add_force(PairPotential::new(
        PairKind::LennardJones {
            epsilon: 1.0,
            sigma: 1.0,
        },
        2.5,
        true,
        Some(Vec3::new(side, side, side)),
    ));
    (w, side)
}

/// g(r) histogram accumulated over snapshots (minimum image).
fn rdf(pos: &[Vec3], side: f64, bins: &mut [f64], dr: f64) {
    let n = pos.len();
    for i in 0..n {
        for j in i + 1..n {
            let mut d = pos[j] - pos[i];
            d = Vec3::new(
                d.x - side * (d.x / side).round(),
                d.y - side * (d.y / side).round(),
                d.z - side * (d.z / side).round(),
            );
            let k = (d.norm() / dr) as usize;
            if k < bins.len() {
                bins[k] += 2.0;
            }
        }
    }
}

fn main() {
    let steps: usize = std::env::args()
        .nth(1)
        .map_or(100_000, |a| a.parse().unwrap());
    let (rho, t0, dt) = (0.8442, 0.722, 0.005);

    // NVE: melt the lattice with a short Langevin run at T, then switch to Verlet.
    let (mut w, side) = lj_liquid(5, rho, t0, "verlet");
    let n = w.state.len();
    w.set_integrator(Box::new(integrators::Langevin::new(t0, 1.0, 7).unwrap()));
    w.run(dt, 4000, 4000).unwrap();
    w.set_integrator(integrators::by_name("verlet").unwrap());
    let start = Instant::now();
    let traj = w
        .run_with_options(
            dt,
            steps,
            &RunOptions {
                record_every: steps / 100,
                events: &[],
                energies: true,
            },
        )
        .unwrap();
    let elapsed = start.elapsed().as_secs_f64();
    let e = traj.total_energy();
    let drift = (e[e.len() - 1] - e[0]) / e[0].abs();
    let spread = e
        .iter()
        .map(|x| ((x - e[0]) / e[0]).abs())
        .fold(0.0, f64::max);
    let temps: Vec<f64> = traj
        .kinetic
        .iter()
        .map(|k| 2.0 * k / (3.0 * n as f64))
        .collect();
    println!(
        "NVE, N = {n}, {steps} steps of {dt}: {:.1} µs/step; E/N = {:.4}; relative energy drift {drift:.1e}, max deviation {spread:.1e}; <T> = {:.3}",
        elapsed / steps as f64 * 1e6,
        e[0] / n as f64,
        temps.iter().sum::<f64>() / temps.len() as f64
    );
    println!("P = {:.3}", w.pressure(side * side * side).unwrap());

    // g(r) over the run's last frames.
    let dr = 0.02;
    let mut bins = vec![0.0; (side / 2.0 / dr) as usize];
    let frames = 50;
    for f in 0..frames {
        let frame = traj.n_frames() - 1 - f;
        rdf(&traj.pos[frame * n..(frame + 1) * n], side, &mut bins, dr);
    }
    let g: Vec<f64> = bins
        .iter()
        .enumerate()
        .map(|(k, c)| {
            let (r0, r1) = (k as f64 * dr, (k + 1) as f64 * dr);
            let shell = 4.0 / 3.0 * std::f64::consts::PI * (r1.powi(3) - r0.powi(3));
            c / (frames as f64 * n as f64 * rho * shell)
        })
        .collect();
    let peak = (0..g.len()).max_by(|&a, &b| g[a].total_cmp(&g[b])).unwrap();
    let after = peak + (0.2 / dr) as usize;
    let min = (after..(2.0 / dr) as usize)
        .min_by(|&a, &b| g[a].total_cmp(&g[b]))
        .unwrap();
    let second = (min..(2.6 / dr) as usize)
        .max_by(|&a, &b| g[a].total_cmp(&g[b]))
        .unwrap();
    let r = |k: usize| (k as f64 + 0.5) * dr;
    println!(
        "g(r): first peak {:.2} at r = {:.2}, minimum {:.2} at {:.2}, second peak {:.2} at {:.2}",
        g[peak],
        r(peak),
        g[min],
        r(min),
        g[second],
        r(second)
    );

    // Thermostats: start cold, target T = 1.
    for (name, integrator) in [
        (
            "langevin",
            Box::new(integrators::Langevin::new(1.0, 1.0, 3).unwrap()) as Box<dyn Integrator>,
        ),
        (
            "nose_hoover",
            Box::new(integrators::NoseHoover::new(1.0, 0.5).unwrap()),
        ),
    ] {
        let (mut w, _) = lj_liquid(4, rho, 0.3, "verlet");
        w.set_integrator(integrator);
        w.run(dt, 4000, 4000).unwrap();
        let e0 = w.total_energy().unwrap() + w.thermostat_energy();
        let traj = w.run(dt, 20_000, 20).unwrap();
        let n = w.state.len() as f64;
        let t: Vec<f64> = traj.kinetic.iter().map(|k| 2.0 * k / (3.0 * n)).collect();
        let mean = t.iter().sum::<f64>() / t.len() as f64;
        let var = t.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / t.len() as f64;
        let ext = w.total_energy().unwrap() + w.thermostat_energy();
        println!(
            "{name}: <T> = {mean:.4} (target 1), relative fluctuation {:.3} (canonical kinetic: sqrt(2/3N) = {:.3}); extended energy change {:.1e}",
            var.sqrt() / mean,
            (2.0 / (3.0 * n)).sqrt(),
            (ext - e0) / e0.abs()
        );
    }
}
