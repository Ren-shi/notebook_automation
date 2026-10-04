//! Noisy runs are reproducible: bit-identical for any thread count, and across checkpoint
//! restarts, because every random number is a pure function of (seed, step, particle).

use physim::integrators::Langevin;
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

/// A Langevin-thermostatted cluster large enough for the blocked, multi-threaded gravity sum.
fn noisy_world(n: usize, seed: u64) -> World {
    let mut rng = Lcg(11);
    let mut w = World::with_integrator("verlet").unwrap();
    for _ in 0..n {
        let r = Vec3::new(rng.next(), rng.next(), rng.next());
        w.add_particle(r, Vec3::ZERO, 0.5 + rng.next()).unwrap();
    }
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.02,
    });
    w.seed = seed;
    w.set_integrator(Box::new(Langevin::new(0.5, 2.0, w.seed).unwrap()));
    w
}

fn bits(w: &World) -> Vec<u64> {
    w.state
        .pos
        .iter()
        .chain(&w.state.vel)
        .flat_map(|v| [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()])
        .collect()
}

#[test]
fn noisy_run_is_identical_for_any_thread_count() {
    let mut results = Vec::new();
    for threads in [1, 4] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let out = pool.install(|| {
            let mut w = noisy_world(1500, 7);
            for _ in 0..5 {
                w.step(1e-3).unwrap();
            }
            bits(&w)
        });
        results.push(out);
    }
    assert_eq!(results[0], results[1]);
    // The noise really acts: another seed gives another run.
    let mut other = noisy_world(1500, 8);
    for _ in 0..5 {
        other.step(1e-3).unwrap();
    }
    assert_ne!(bits(&other), results[0]);
}

#[test]
fn noisy_run_restarts_exactly_from_a_checkpoint() {
    let mut straight = noisy_world(200, 3);
    for _ in 0..40 {
        straight.step(1e-3).unwrap();
    }

    let mut first = noisy_world(200, 3);
    for _ in 0..17 {
        first.step(1e-3).unwrap();
    }
    let saved = first.checkpoint();
    assert_eq!(saved.seed, 3);
    let mut resumed = World::from_checkpoint(saved, |_, _| unreachable!()).unwrap();
    assert_eq!(resumed.seed, 3);
    for _ in 17..40 {
        resumed.step(1e-3).unwrap();
    }
    assert_eq!(bits(&resumed), bits(&straight));
}

#[test]
fn random_numbers_are_a_function_of_the_key() {
    assert_eq!(rng::normal(5, 2, 9), rng::normal(5, 2, 9));
    assert_ne!(rng::uniform(5, 2, 9), rng::uniform(5, 2, 10));
    assert_ne!(rng::uniform(5, 2, 9), rng::uniform(5, rng::SETUP + 2, 9));
}
