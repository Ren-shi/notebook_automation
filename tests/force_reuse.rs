//! The Verlet-family integrators reuse the end-of-step acceleration. These tests pin down
//! when that is allowed and check that results are unchanged by it.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use physim::*;

/// Harmonic trap that counts its evaluations.
struct Counted {
    calls: Arc<AtomicUsize>,
    velocity_dependent: bool,
}

impl Force for Counted {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _v: &[Vec3],
        _m: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        for (a, r) in acc.iter_mut().zip(pos) {
            *a -= *r;
        }
        Ok(())
    }
    fn velocity_dependent(&self) -> bool {
        self.velocity_dependent
    }
    fn name(&self) -> String {
        "Counted".into()
    }
}

fn counted_world(integrator: &str, velocity_dependent: bool) -> (World, Arc<AtomicUsize>) {
    let calls = Arc::new(AtomicUsize::new(0));
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.5, 0.0), 1.0)
        .unwrap();
    w.add_force(Counted {
        calls: calls.clone(),
        velocity_dependent,
    });
    (w, calls)
}

fn evals(integrator: &str, velocity_dependent: bool, steps: usize) -> usize {
    let (mut w, calls) = counted_world(integrator, velocity_dependent);
    for _ in 0..steps {
        w.step(0.01).unwrap();
    }
    calls.load(Ordering::Relaxed)
}

#[test]
fn verlet_and_yoshida_reuse_accelerations_for_conservative_forces() {
    assert_eq!(evals("verlet", false, 100), 101);
    assert_eq!(evals("yoshida4", false, 100), 301);
    // Velocity-dependent forces are re-evaluated every time.
    assert_eq!(evals("verlet", true, 100), 200);
    assert_eq!(evals("yoshida4", true, 100), 600);
}

#[test]
fn external_changes_invalidate_the_cache() {
    let (mut w, calls) = counted_world("verlet", false);
    w.step(0.01).unwrap();
    let n = calls.load(Ordering::Relaxed);

    w.state.pos[0].x += 0.1; // edited positions
    w.step(0.01).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), n + 2);

    w.state.t = 5.0; // edited time
    w.step(0.01).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), n + 4);

    w.set_masses(vec![2.0]).unwrap(); // edited masses
    w.step(0.01).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), n + 6);

    let g = w.add_force(UniformField {
        g: Vec3::new(0.0, -1.0, 0.0),
    }); // new force
    w.step(0.01).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), n + 8);

    w.set_force_params(g, &[("g".into(), Param::Vector(Vec3::ZERO))])
        .unwrap(); // new parameters
    w.step(0.01).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), n + 10);

    w.step(0.01).unwrap(); // nothing changed: one evaluation
    assert_eq!(calls.load(Ordering::Relaxed), n + 11);
}

/// Reference: velocity Verlet written out with two fresh evaluations per step.
fn reference_verlet(w: &mut World, dt: f64, steps: usize) {
    for _ in 0..steps {
        let a = w.accelerations().unwrap();
        for ((x, v), a) in w.state.pos.iter_mut().zip(w.state.vel.iter_mut()).zip(&a) {
            *v += *a * (0.5 * dt);
            *x += *v * dt;
        }
        w.state.t += dt;
        let a = w.accelerations().unwrap();
        for (v, a) in w.state.vel.iter_mut().zip(&a) {
            *v += *a * (0.5 * dt);
        }
    }
}

#[test]
fn reuse_gives_bit_identical_trajectories() {
    let build = || {
        let mut w = World::with_integrator("verlet").unwrap();
        w.add_particle(Vec3::new(0.5, 0.0, 0.0), Vec3::new(0.0, 0.4, 0.1), 0.5)
            .unwrap();
        w.add_particle(Vec3::new(-0.5, 0.0, 0.0), Vec3::new(0.0, -0.4, 0.0), 0.5)
            .unwrap();
        w.add_particle(Vec3::new(0.0, 1.5, 0.0), Vec3::new(-0.3, 0.0, 0.0), 0.01)
            .unwrap();
        w.add_force(NewtonianGravity {
            g: 1.0,
            softening: 0.0,
        });
        w.add_force(UniformField {
            g: Vec3::new(0.0, 0.0, -0.01),
        });
        w
    };
    let mut cached = build();
    cached.run(1e-3, 5_000, 5_000).unwrap();
    let mut reference = build();
    reference_verlet(&mut reference, 1e-3, 5_000);
    assert_eq!(cached.state.pos, reference.state.pos);
    assert_eq!(cached.state.vel, reference.state.vel);
}
