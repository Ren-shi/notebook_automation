//! Engine benchmarks: `cargo bench` (HTML reports in target/criterion/).
//! Baseline numbers are recorded in backlog/02-benchmarks.md.

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use physim::integrators::NAMES;
use physim::*;

/// Deterministic pseudo-random numbers in [0, 1) (64-bit LCG), so runs are comparable.
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

/// N bodies of total mass 1 uniformly in the unit cube, at rest, softened gravity.
fn cluster(n: usize, integrator: &str) -> World {
    let mut rng = Lcg(42);
    let mut w = World::with_integrator(integrator).unwrap();
    for _ in 0..n {
        let r = Vec3::new(rng.next(), rng.next(), rng.next());
        w.add_particle(r, Vec3::ZERO, 1.0 / n as f64).unwrap();
    }
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.01,
    });
    w
}

/// Chain of n unit masses joined by springs, ends anchored: cheap O(N) forces.
fn chain(n: usize, integrator: &str) -> World {
    let mut w = World::with_integrator(integrator).unwrap();
    for i in 0..n {
        let y = if i == n / 2 { 0.1 } else { 0.0 };
        w.add_particle(Vec3::new(i as f64, y, 0.0), Vec3::ZERO, 1.0)
            .unwrap();
    }
    for i in 0..n - 1 {
        w.add_force(Spring {
            i,
            j: i + 1,
            k: 100.0,
            rest_length: 0.8,
        });
    }
    w.add_force(AnchorSpring {
        i: 0,
        anchor: Vec3::new(-1.0, 0.0, 0.0),
        k: 100.0,
        rest_length: 0.8,
    });
    w.add_force(AnchorSpring {
        i: n - 1,
        anchor: Vec3::new(n as f64, 0.0, 0.0),
        k: 100.0,
        rest_length: 0.8,
    });
    w
}

/// Direct-sum gravity: one Verlet step (2 force evaluations). Throughput is in pair interactions.
fn gravity(c: &mut Criterion) {
    let mut g = c.benchmark_group("gravity_verlet_step");
    for n in [100usize, 1_000, 5_000] {
        if n >= 1_000 {
            g.sample_size(10);
        }
        g.throughput(Throughput::Elements((n * (n - 1) / 2) as u64));
        g.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            let mut w = cluster(n, "verlet");
            b.iter(|| w.step(black_box(1e-4)).unwrap());
        });
    }
    g.finish();
}

/// Barnes-Hut tree gravity, one force evaluation per Verlet step (compare `gravity_verlet_step`;
/// `cargo run --release --example tree_gravity` compares against direct summation at N = 1e5).
fn tree_gravity(c: &mut Criterion) {
    let mut g = c.benchmark_group("tree_gravity_verlet_step");
    g.sample_size(10);
    for n in [5_000usize, 100_000] {
        g.throughput(Throughput::Elements(n as u64));
        g.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            let mut w = cluster(n, "verlet");
            w.forces.clear();
            w.add_force(TreeGravity {
                g: 1.0,
                softening: 0.01,
                theta: 0.5,
                quadrupole: false,
            });
            b.iter(|| w.step(black_box(1e-4)).unwrap());
        });
    }
    g.finish();
}

/// Each integrator on cheap forces (integrator overhead) and on gravity (force-dominated).
fn integrators(c: &mut Criterion) {
    let mut g = c.benchmark_group("integrator_step_chain10k");
    g.throughput(Throughput::Elements(10_000));
    for &name in NAMES {
        let mut w = chain(10_000, name);
        g.bench_function(name, |b| b.iter(|| w.step(black_box(1e-3)).unwrap()));
    }
    g.finish();

    let mut g = c.benchmark_group("integrator_step_gravity200");
    for &name in NAMES {
        let mut w = cluster(200, name);
        g.bench_function(name, |b| b.iter(|| w.step(black_box(1e-4)).unwrap()));
    }
    g.finish();
}

/// Cost of recording frames: 1 000 steps of a 1 000-particle chain, recording every step vs
/// only the ends, with and without per-frame energies.
fn recording(c: &mut Criterion) {
    let mut g = c.benchmark_group("run_chain1k_1000steps");
    g.sample_size(20);
    for every in [1usize, 1_000] {
        g.bench_with_input(
            BenchmarkId::new("record_every", every),
            &every,
            |b, &every| {
                let mut w = chain(1_000, "verlet");
                b.iter(|| w.run(1e-3, 1_000, every).unwrap());
            },
        );
    }
    g.bench_function("record_every/1/no_energies", |b| {
        let mut w = chain(1_000, "verlet");
        let options = RunOptions {
            energies: false,
            ..RunOptions::default()
        };
        b.iter(|| w.run_with_options(1e-3, 1_000, &options).unwrap());
    });
    g.finish();
}

/// One Verlet step of a 10 000-particle chain: one `Spring` per bond versus a single
/// `SpringNetwork` holding all bonds.
fn springs(c: &mut Criterion) {
    let n = 10_000;
    let mut g = c.benchmark_group("springs_chain10k_verlet_step");
    g.throughput(Throughput::Elements((n - 1) as u64));
    let mut separate = chain(n, "verlet");
    g.bench_function("separate", |b| {
        b.iter(|| separate.step(black_box(1e-3)).unwrap())
    });
    let mut network = chain(n, "verlet");
    network.forces.clear();
    network.add_force(
        SpringNetwork::new(
            (0..n - 1).collect(),
            (1..n).collect(),
            vec![100.0; n - 1],
            vec![0.8; n - 1],
        )
        .unwrap(),
    );
    g.bench_function("network", |b| {
        b.iter(|| network.step(black_box(1e-3)).unwrap())
    });
    g.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default().warm_up_time(Duration::from_secs(1)).measurement_time(Duration::from_secs(3));
    targets = gravity, tree_gravity, integrators, recording, springs
}
criterion_main!(benches);
