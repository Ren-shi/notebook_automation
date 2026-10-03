//! Blocked (multi-threaded) gravity must agree with a plain direct sum and give
//! bit-identical results for any thread count.

use physim::*;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// N = 1500 gives ~1.1 M pairs, well inside the multi-block regime.
fn cluster(n: usize) -> World {
    let mut rng = Lcg(7);
    let mut w = World::with_integrator("verlet").unwrap();
    for _ in 0..n {
        let r = Vec3::new(rng.next(), rng.next(), rng.next());
        let v = Vec3::new(rng.next() - 0.5, rng.next() - 0.5, rng.next() - 0.5) * 0.1;
        w.add_particle(r, v, 0.5 + rng.next()).unwrap();
    }
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.01,
    });
    w
}

fn direct_sum(w: &World) -> (Vec<Vec3>, f64) {
    let (pos, mass) = (&w.state.pos, &w.state.mass);
    let eps2 = 0.01 * 0.01;
    let mut acc = vec![Vec3::ZERO; pos.len()];
    let mut u = 0.0;
    for i in 0..pos.len() {
        for j in 0..pos.len() {
            if i != j {
                let d = pos[j] - pos[i];
                let r2 = d.norm_squared() + eps2;
                acc[i] += d * (mass[j] / (r2 * r2.sqrt()));
                if i < j {
                    u -= mass[i] * mass[j] / r2.sqrt();
                }
            }
        }
    }
    (acc, u)
}

#[test]
fn blocked_gravity_matches_direct_sum() {
    let w = cluster(1500);
    let (expected, u_expected) = direct_sum(&w);
    let acc = w.accelerations().unwrap();
    let scale = expected.iter().map(|a| a.norm()).fold(0.0, f64::max);
    for (a, e) in acc.iter().zip(&expected) {
        assert!((*a - *e).norm() <= 1e-12 * scale);
    }
    let u = w.potential_energy().unwrap();
    assert!((u - u_expected).abs() <= 1e-12 * u_expected.abs());
}

#[test]
fn blocked_gravity_conserves_momentum() {
    let mut w = cluster(1500);
    let p0 = w.state.momentum();
    w.run(1e-4, 20, 20).unwrap();
    let scale: f64 = w
        .state
        .vel
        .iter()
        .zip(&w.state.mass)
        .map(|(v, m)| m * v.norm())
        .sum();
    assert!((w.state.momentum() - p0).norm() < 1e-13 * scale);
}

#[test]
fn results_do_not_depend_on_thread_count() {
    let run = |threads: usize| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let mut w = cluster(1500);
            w.run(1e-4, 5, 5).unwrap();
            let u = w.potential_energy().unwrap();
            (w.state.pos, u)
        })
    };
    let (pos1, u1) = run(1);
    let (pos4, u4) = run(4);
    assert_eq!(pos1, pos4);
    assert_eq!(u1.to_bits(), u4.to_bits());
}

#[test]
fn per_particle_forces_are_unchanged_in_parallel() {
    // 20 000 particles: above the per-particle parallel threshold.
    let mut w = World::with_integrator("rk4").unwrap();
    for i in 0..20_000 {
        let x = i as f64;
        w.add_particle(Vec3::new(x, 0.0, 0.0), Vec3::new(0.0, x * 1e-3, 0.0), 1.0)
            .unwrap();
    }
    w.add_force(UniformField {
        g: Vec3::new(0.0, -9.81, 0.0),
    });
    w.add_force(LinearDrag { gamma: 0.5 });
    let acc = w.accelerations().unwrap();
    for (i, a) in acc.iter().enumerate() {
        let expected = Vec3::new(0.0, -9.81, 0.0) - Vec3::new(0.0, i as f64 * 1e-3, 0.0) * 0.5;
        assert_eq!(*a, expected);
    }
}

#[test]
fn particle_mesh_does_not_depend_on_thread_count() {
    let run = |threads: usize, periodic: bool| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let mut w = cluster(500);
            w.forces.clear();
            let center = Vec3::new(0.5, 0.5, 0.5);
            w.add_force(ParticleMesh::new(1.0, 16, 2.0, center, periodic, None));
            w.run(1e-3, 3, 3).unwrap();
            w.state.pos
        })
    };
    for periodic in [false, true] {
        assert_eq!(run(1, periodic), run(4, periodic));
    }
}
