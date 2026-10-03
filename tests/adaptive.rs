use std::f64::consts::{PI, TAU};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use physim::events::{CoordinateCrossing, RadialVelocity};
use physim::*;

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

/// Unit-mass sun pinned at the origin and a light planet at apoapsis of an a = 1 orbit
/// with eccentricity `e` (period 2π).
fn orbit(integrator: &str, e: f64) -> (World, Arc<AtomicU64>) {
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.pin(0, true).unwrap();
    let v = ((1.0 - e) / (1.0 + e)).sqrt();
    w.add_particle(Vec3::new(1.0 + e, 0.0, 0.0), Vec3::new(0.0, v, 0.0), 1e-6)
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

fn options(rtol: f64) -> AdaptiveOptions<'static> {
    AdaptiveOptions {
        rtol,
        atol: rtol * 1e-3,
        ..Default::default()
    }
}

#[test]
fn eccentric_orbit_needs_far_fewer_evaluations_than_fixed_steps() {
    let t_end = 3.0 * TAU;
    let (mut w, count) = orbit("verlet", 0.99);
    let run = w.run_adaptive(t_end, &options(1e-10)).unwrap();
    let target = energy_error(&run.trajectory);
    let evals = count.load(Ordering::Relaxed);
    assert!(target < 1e-7, "adaptive energy error {target:e}");
    assert_eq!(evals, run.stats.evaluations);
    assert_eq!(w.state.t, t_end);
    // Every step costs six evaluations (FSAL), plus the first stage and the step guess.
    let s = run.stats;
    assert_eq!(evals, 6 * (s.accepted + s.rejected) as u64 + 2);
    assert_eq!(run.trajectory.n_frames(), s.accepted + 1);

    // With 20 times the evaluations, no fixed-step method gets there
    // (`cargo run --release --example eccentric_orbit` measures 500-15000 times).
    for (name, per_step) in [("verlet", 1), ("yoshida4", 3), ("rk4", 4), ("dopri5", 6)] {
        let steps = 20 * evals as usize / per_step;
        let (mut w, count) = orbit(name, 0.99);
        let err = match w.run(t_end / steps as f64, steps, 1) {
            Ok(traj) => energy_error(&traj),
            Err(_) => f64::INFINITY,
        };
        let used = count.load(Ordering::Relaxed);
        assert!(used <= 20 * evals + 1, "{name}: {used} evaluations");
        assert!(
            err > 10.0 * target,
            "{name}: {err:e} vs adaptive {target:e}"
        );
    }
}

#[test]
fn tolerance_controls_the_error() {
    // Harmonic oscillator x = cos t: the global error tracks the tolerance.
    let mut previous = f64::INFINITY;
    for rtol in [1e-6, 1e-9, 1e-12] {
        let mut w = World::with_integrator("rk4").unwrap();
        w.add_particle(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0)
            .unwrap();
        w.add_force(AnchorSpring {
            i: 0,
            anchor: Vec3::ZERO,
            k: 1.0,
            rest_length: 0.0,
        });
        w.run_adaptive(10.0, &options(rtol)).unwrap();
        let err = (w.state.pos[0].x - 10f64.cos()).abs() + (w.state.vel[0].x + 10f64.sin()).abs();
        assert!(err < 100.0 * rtol, "rtol {rtol:e}: error {err:e}");
        assert!(err < previous / 100.0);
        previous = err;
    }
}

#[test]
fn output_at_requested_times_uses_dense_output() {
    // A damped oscillator (velocity-dependent force): x'' = -x - 0.2 x'.
    let damped = || {
        let mut w = World::with_integrator("verlet").unwrap();
        w.add_particle(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0)
            .unwrap();
        w.add_force(AnchorSpring {
            i: 0,
            anchor: Vec3::ZERO,
            k: 1.0,
            rest_length: 0.0,
        });
        w.add_force(LinearDrag { gamma: 0.2 });
        w
    };
    let mut w = damped();
    let (g, wd) = (0.1, (1.0f64 - 0.01).sqrt());
    let exact = |t: f64| (-g * t).exp() * ((wd * t).cos() + g / wd * (wd * t).sin());
    let times: Vec<f64> = (0..=200).map(|k| k as f64 * 0.05).collect();
    let run = w
        .run_adaptive(
            10.0,
            &AdaptiveOptions {
                output: Output::Times(&times),
                ..options(1e-10)
            },
        )
        .unwrap();
    let traj = &run.trajectory;
    assert_eq!(traj.t, times);
    // Output times do not affect the steps taken.
    let steps = damped().run_adaptive(10.0, &options(1e-10)).unwrap();
    assert_eq!(run.stats, steps.stats);
    assert_eq!(steps.trajectory.pos.last(), traj.pos.last());
    for (k, &t) in times.iter().enumerate() {
        let x = traj.pos[k].x;
        assert!((x - exact(t)).abs() < 1e-8, "t = {t}: {x} vs {}", exact(t));
    }
    assert_eq!(w.state.t, 10.0);
}

#[test]
fn events_are_located_on_the_continuous_extension() {
    // Periapsis passages of an e = 0.9 orbit come exactly one period (2π) apart.
    let (mut w, _) = orbit("verlet", 0.9);
    let events = [Event::new(
        RadialVelocity { i: 1, j: Some(0) },
        Direction::Rising,
        false,
    )];
    let run = w
        .run_adaptive(
            5.0 * TAU,
            &AdaptiveOptions {
                events: &events,
                ..options(1e-12)
            },
        )
        .unwrap();
    let hits = &run.trajectory.events;
    assert_eq!(hits.len(), 5);
    for (k, hit) in hits.iter().enumerate() {
        let expected = (k as f64 + 0.5) * TAU;
        assert!((hit.t - expected).abs() < 1e-8, "{}: {}", k, hit.t);
        assert!((hit.pos[1].norm() - 0.1).abs() < 1e-9);
    }

    // A terminal event stops the run at the event, recorded as the last frame.
    let (mut w, _) = orbit("verlet", 0.9);
    let events = [Event::new(
        CoordinateCrossing {
            i: 1,
            axis: 0,
            value: 0.0,
        },
        Direction::Falling,
        true,
    )];
    let run = w
        .run_adaptive(
            TAU,
            &AdaptiveOptions {
                events: &events,
                ..options(1e-12)
            },
        )
        .unwrap();
    assert_eq!(run.trajectory.terminated_by, Some(0));
    assert!(w.state.pos[1].x.abs() < 1e-12);
    assert_eq!(*run.trajectory.t.last().unwrap(), w.state.t);
    // From apoapsis, x = a (e - cos E) reaches zero at cos E = e with E in (π, 2π);
    // Kepler's equation M = E - e sin E then gives the time since apoapsis.
    let e: f64 = 0.9;
    let expected = PI - e.acos() + e * (1.0 - e * e).sqrt();
    assert!(
        (w.state.t - expected).abs() < 1e-9,
        "{} vs {expected}",
        w.state.t
    );
}

#[test]
fn integrates_backwards_to_the_start() {
    let (mut w, _) = orbit("verlet", 0.7);
    let start = w.state.clone();
    w.run_adaptive(4.0, &options(1e-12)).unwrap();
    assert_eq!(w.state.t, 4.0);
    let times = [4.0, 2.0, 0.0];
    let run = w
        .run_adaptive(
            0.0,
            &AdaptiveOptions {
                output: Output::Times(&times),
                ..options(1e-12)
            },
        )
        .unwrap();
    assert_eq!(run.trajectory.t, times);
    assert_eq!(w.state.t, 0.0);
    assert!((w.state.pos[1] - start.pos[1]).norm() < 1e-9);
    assert!((w.state.vel[1] - start.vel[1]).norm() < 1e-9);
}

#[test]
fn streaming_matches_in_memory() {
    let (mut a, _) = orbit("verlet", 0.5);
    let (mut b, _) = orbit("verlet", 0.5);
    let run = a.run_adaptive(TAU, &options(1e-9)).unwrap();
    let mut traj = Trajectory::new(2, true);
    let outcome = b.run_adaptive_into(TAU, &options(1e-9), &mut traj).unwrap();
    assert_eq!(outcome.stats, run.stats);
    assert_eq!(traj.t, run.trajectory.t);
    assert_eq!(traj.pos, run.trajectory.pos);
    assert_eq!(a.state.pos, b.state.pos);
}

#[test]
fn rejects_bad_input_and_reports_failure() {
    let (mut w, _) = orbit("verlet", 0.5);
    let bad = |w: &mut World, t_end: f64, o: AdaptiveOptions<'_>| {
        w.run_adaptive(t_end, &o).err().map(|f| f.error.to_string())
    };
    assert!(bad(&mut w, f64::NAN, options(1e-9))
        .unwrap()
        .contains("t_end"));
    for o in [
        options(-1.0),
        AdaptiveOptions {
            rtol: 0.0,
            atol: 0.0,
            ..Default::default()
        },
        options(1e-17),
        AdaptiveOptions {
            max_step: 0.0,
            ..Default::default()
        },
        AdaptiveOptions {
            first_step: Some(-0.1),
            ..Default::default()
        },
        AdaptiveOptions {
            output: Output::Times(&[0.5, 0.2]),
            ..Default::default()
        },
        AdaptiveOptions {
            output: Output::Times(&[2.0]),
            ..Default::default()
        },
    ] {
        assert!(bad(&mut w, 1.0, o).is_some());
    }
    assert_eq!(w.state.t, 0.0);

    // max_steps: the world is left at the last accepted step, which was recorded.
    let failure = w
        .run_adaptive(
            10.0,
            &AdaptiveOptions {
                max_steps: 5,
                ..Default::default()
            },
        )
        .err()
        .unwrap();
    assert!(failure.error.to_string().contains("max_steps"));
    assert_eq!(failure.trajectory.n_frames(), 6);
    assert_eq!(*failure.trajectory.t.last().unwrap(), w.state.t);

    // A collision course: the step shrinks until it underflows.
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.add_particle(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0)
        .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    let failure = w.run_adaptive(10.0, &options(1e-9)).err().unwrap();
    let msg = failure.error.to_string();
    assert!(
        msg.contains("underflow") || msg.contains("max_steps"),
        "{msg}"
    );
    // Head-on free fall from rest at separation 1 with G(m1 + m2) = 2 collides at t = π/4.
    assert!((w.state.t - PI / 4.0).abs() < 1e-8, "{}", w.state.t);

    // Constraints are not supported.
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0)
        .unwrap();
    w.add_rod(0, Anchor::Point(Vec3::ZERO), None).unwrap();
    assert!(bad(&mut w, 1.0, options(1e-9))
        .unwrap()
        .contains("constraints"));
}
