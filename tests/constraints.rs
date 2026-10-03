use std::f64::consts::PI;

use physim::events::CoordinateCrossing;
use physim::*;

fn p(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

const G: f64 = 9.81;

/// Complete elliptic integral of the first kind K(k), by the arithmetic-geometric mean.
fn ellip_k(k: f64) -> f64 {
    let (mut a, mut b) = (1.0, (1.0 - k * k).sqrt());
    // Converges quadratically; a fixed count avoids looping on a one-ulp difference.
    for _ in 0..40 {
        (a, b) = (0.5 * (a + b), (a * b).sqrt());
    }
    PI / (2.0 * a)
}

/// Exact period of a rigid pendulum released from rest at angle `theta0`.
fn exact_period(length: f64, theta0: f64) -> f64 {
    4.0 * (length / G).sqrt() * ellip_k((theta0 / 2.0).sin())
}

/// A bob on a rod of `length` from the origin, released from rest at `theta0` from the
/// downward vertical, under gravity -y.
fn pendulum(integrator: &str, length: f64, theta0: f64) -> World {
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(
        p(length * theta0.sin(), -length * theta0.cos(), 0.0),
        Vec3::ZERO,
        2.0,
    )
    .unwrap();
    w.add_force(UniformField { g: p(0.0, -G, 0.0) });
    w.add_rod(0, Anchor::Point(Vec3::ZERO), Some(length))
        .unwrap();
    w
}

/// Period from the times of successive returns to the release side (x rising through 0).
fn measured_period(w: &mut World, dt: f64, periods: f64, estimate: f64) -> f64 {
    let crossing = Event {
        function: Box::new(CoordinateCrossing {
            i: 0,
            axis: 0,
            value: 0.0,
        }),
        direction: Direction::Falling,
        terminal: false,
    };
    let steps = (periods * estimate / dt) as usize;
    let traj = w.run_with_events(dt, steps, steps, &[crossing]).unwrap();
    let t: Vec<f64> = traj.events.iter().map(|h| h.t).collect();
    assert!(t.len() >= 2, "only {} crossings", t.len());
    (t[t.len() - 1] - t[0]) / (t.len() - 1) as f64
}

#[test]
fn small_amplitude_period_matches_two_pi_sqrt_l_over_g() {
    let (length, theta0) = (1.5, 1e-3);
    // yoshida4: velocity Verlet's own phase error, (omega dt)^2 / 24 ~ 3e-7 at this dt,
    // would hide the comparison.
    let mut w = pendulum("yoshida4", length, theta0);
    let small = 2.0 * PI * (length / G).sqrt();
    let t = measured_period(&mut w, 1e-3, 20.0, small);
    // The finite-amplitude correction is theta0^2 / 16 ~ 6e-8.
    assert!(
        (t / small - 1.0).abs() < 1e-7,
        "period {t}, expected {small}"
    );
    let exact = exact_period(length, theta0);
    assert!((t / exact - 1.0).abs() < 1e-10, "period {t}, exact {exact}");
}

#[test]
fn large_amplitude_period_matches_elliptic_integral() {
    let length = 0.8;
    for theta0 in [1.0, 2.0, 3.0] {
        let exact = exact_period(length, theta0);
        for (integrator, dt, tol) in [("verlet", 2e-4, 1e-6), ("yoshida4", 1e-3, 1e-8)] {
            let mut w = pendulum(integrator, length, theta0);
            let t = measured_period(&mut w, dt, 10.0, exact);
            let err = (t / exact - 1.0).abs();
            assert!(
                err < tol,
                "{integrator}, theta0 = {theta0}: period {t}, exact {exact}, error {err:.2e}"
            );
        }
    }
}

/// Checks the rods and the energy at every step of a run, keeping the worst energy error
/// of each block of `block` frames.
struct Watch {
    rods: Vec<Rod>,
    e0: f64,
    block: usize,
    frames: usize,
    worst_violation: f64,
    worst_energy: Vec<f64>,
}

impl Recorder for Watch {
    fn frame(&mut self, f: &Frame<'_>) -> Result<()> {
        for r in &self.rods {
            let v = (r.separation(f.pos).norm() / r.length - 1.0).abs();
            self.worst_violation = self.worst_violation.max(v);
        }
        let e = f.kinetic.unwrap() + f.potential.unwrap();
        if self.frames.is_multiple_of(self.block) {
            self.worst_energy.push(0.0);
        }
        let worst = self.worst_energy.last_mut().unwrap();
        *worst = worst.max(((e - self.e0) / self.e0).abs());
        self.frames += 1;
        Ok(())
    }
}

#[test]
fn double_pendulum_holds_lengths_for_a_million_steps_with_bounded_energy() {
    for (integrator, energy_bound) in [("verlet", 5e-4), ("yoshida4", 1e-6)] {
        let mut w = World::with_integrator(integrator).unwrap();
        w.add_particle(p(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
        w.add_particle(p(1.0, 0.7, 0.0), Vec3::ZERO, 0.5).unwrap();
        w.add_force(UniformField { g: p(0.0, -G, 0.0) });
        w.add_rod(0, Anchor::Point(Vec3::ZERO), None).unwrap();
        w.add_rod(1, Anchor::Particle(0), None).unwrap();
        let mut watch = Watch {
            rods: w.constraints.iter().map(|(_, r)| *r).collect(),
            e0: w.total_energy().unwrap(),
            block: 100_000,
            frames: 0,
            worst_violation: 0.0,
            worst_energy: Vec::new(),
        };
        // 10^6 steps of a chaotic double pendulum, every step checked.
        w.run_into(1e-3, 1_000_000, &RunOptions::default(), &mut watch)
            .unwrap();
        let tol = w.constraints.tolerance;
        assert!(
            watch.worst_violation <= tol,
            "{integrator}: worst violation {:.2e}",
            watch.worst_violation
        );
        let e = &watch.worst_energy;
        assert_eq!(e.len(), 11);
        // Bounded, not drifting: the last block is no worse than the first.
        assert!(
            e.iter().all(|&x| x < energy_bound) && e[9] < 2.0 * e[0],
            "{integrator}: worst energy error per 1e5 steps {e:?}"
        );
    }
}

#[test]
fn tension_balances_gravity_at_rest_and_adds_centripetal_force_in_a_swing() {
    let mut w = pendulum("verlet", 1.0, 0.0);
    w.run(1e-3, 100, 100).unwrap();
    let t = w.constraint_tensions()[0];
    assert!((t / (2.0 * G) - 1.0).abs() < 1e-9, "hanging tension {t}");

    // Released horizontally, the bob reaches the bottom with v^2 = 2 g L: T = 3 m g there.
    let mut w = pendulum("yoshida4", 1.0, PI / 2.0);
    let bottom = Event {
        function: Box::new(CoordinateCrossing {
            i: 0,
            axis: 0,
            value: 0.0,
        }),
        direction: Direction::Falling,
        terminal: true,
    };
    w.run_with_events(1e-3, 10_000, 10_000, &[bottom]).unwrap();
    let t = w.constraint_tensions()[0];
    // The reported tension is the multiplier of the step's last velocity projection, accurate
    // to O(h^2) in that (here partial, event-located) substep: ~1e-5 at this dt.
    assert!(
        (t / (3.0 * 2.0 * G) - 1.0).abs() < 5e-5,
        "tension at the bottom {t}"
    );
}

#[test]
fn free_spinning_dumbbell_keeps_length_momentum_and_centripetal_tension() {
    // Masses 1 and 2 on a rod of length 3 spinning at omega about their centre of mass.
    let (m1, m2, length, omega) = (1.0, 2.0, 3.0, 0.7);
    let (r1, r2) = (length * m2 / (m1 + m2), length * m1 / (m1 + m2));
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(p(r1, 0.0, 0.0), p(0.0, omega * r1, 0.0), m1)
        .unwrap();
    w.add_particle(p(-r2, 0.0, 0.0), p(0.0, -omega * r2, 0.0), m2)
        .unwrap();
    w.add_rod(0, Anchor::Particle(1), None).unwrap();
    let l0 = w.state.angular_momentum();
    w.run(1e-3, 20_000, 20_000).unwrap();
    assert!(w.constraints.max_violation(&w.state.pos) <= 1e-10);
    assert!((w.state.angular_momentum() - l0).norm() < 1e-12);
    assert!(w.state.momentum().norm() < 1e-12);
    let expected = m1 * omega * omega * r1;
    let t = w.constraint_tensions()[0];
    assert!(
        (t / expected - 1.0).abs() < 1e-6,
        "tension {t}, expected {expected}"
    );
}

#[test]
fn constrained_integrators_converge_at_their_order() {
    // Angle error at t = 2 against a much finer yoshida4 reference.
    let (length, theta0, t_end) = (1.0, 2.0, 2.0);
    let angle = |integrator: &str, dt: f64| {
        let mut w = pendulum(integrator, length, theta0);
        w.run(dt, (t_end / dt).round() as usize, 1_000_000).unwrap();
        let r = w.state.pos[0];
        r.x.atan2(-r.y)
    };
    let reference = angle("yoshida4", 1e-4);
    for (integrator, order) in [("verlet", 2.0), ("yoshida4", 4.0)] {
        let e1 = (angle(integrator, 0.02) - reference).abs();
        let e2 = (angle(integrator, 0.01) - reference).abs();
        let slope = (e1 / e2).log2();
        assert!(
            (slope - order).abs() < 0.3,
            "{integrator}: slope {slope:.2}"
        );
    }
}

#[test]
fn constraint_management_and_errors() {
    let mut w = World::with_integrator("rk4").unwrap();
    w.add_particle(p(0.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_particle(p(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_particle(p(1.0, 1.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_particle(p(5.0, 1.0, 0.0), Vec3::ZERO, 0.0).unwrap();

    let err = |r: Result<ConstraintId>| r.unwrap_err().to_string();
    assert!(err(w.add_rod(0, Anchor::Particle(1), Some(2.0))).contains("place them first"));
    assert!(err(w.add_rod(0, Anchor::Particle(0), None)).contains("two different"));
    assert!(err(w.add_rod(3, Anchor::Particle(0), None)).contains("massless"));
    assert!(err(w.add_rod(9, Anchor::Particle(0), None)).contains("out of range"));

    let a = w.add_rod(0, Anchor::Particle(1), None).unwrap();
    let b = w.add_rod(2, Anchor::Particle(1), Some(1.0)).unwrap();
    assert!(w
        .step(1e-3)
        .unwrap_err()
        .to_string()
        .contains("verlet or yoshida4"));
    assert_eq!(w.state.t, 0.0, "a refused step changes nothing");
    w.set_integrator(integrators::by_name("verlet").unwrap());
    w.step(1e-3).unwrap();

    assert!(
        w.set_masses(vec![1.0, 0.0, 1.0, 0.0]).is_err(),
        "constrained particle made massless"
    );
    assert_eq!(
        w.state.mass[1], 1.0,
        "masses unchanged after a refused update"
    );
    assert!(w
        .remove_particle(1)
        .unwrap_err()
        .to_string()
        .contains("Rod(0, 1)"));
    w.remove_particle(3).unwrap();
    w.remove_constraint(a).unwrap();
    w.remove_particle(0).unwrap();
    // Rod(2, 1) is now Rod(1, 0).
    assert_eq!(w.constraints.get(b).unwrap().name(), "Rod(1, 0)");
    w.step(1e-3).unwrap();
    assert_eq!(w.constraint_tensions().len(), 1);
}

#[test]
fn pinned_end_acts_as_a_fixed_point() {
    // A rod to a pinned particle behaves like a rod to a fixed point.
    let mut a = pendulum("verlet", 1.0, 0.5);
    let mut b = World::with_integrator("verlet").unwrap();
    b.add_particle(a.state.pos[0], Vec3::ZERO, 2.0).unwrap();
    b.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    b.pin(1, true).unwrap();
    b.add_force(UniformField { g: p(0.0, -G, 0.0) });
    b.add_rod(0, Anchor::Particle(1), None).unwrap();
    a.run(1e-3, 3000, 3000).unwrap();
    b.run(1e-3, 3000, 3000).unwrap();
    assert!((a.state.pos[0] - b.state.pos[0]).norm() < 1e-12);
    assert_eq!(b.state.pos[1], Vec3::ZERO);
}

#[test]
fn checkpoint_restores_constraints_bit_for_bit() {
    let mut straight = pendulum("yoshida4", 1.0, 2.5);
    let mut interrupted = pendulum("yoshida4", 1.0, 2.5);
    straight.run(1e-3, 2000, 2000).unwrap();
    interrupted.run(1e-3, 700, 700).unwrap();
    let mut resumed =
        World::from_checkpoint(interrupted.checkpoint(), |_, _| unreachable!()).unwrap();
    resumed.run(1e-3, 1300, 1300).unwrap();
    assert_eq!(straight.state.pos, resumed.state.pos);
    assert_eq!(straight.state.vel, resumed.state.vel);
    assert_eq!(
        straight.constraint_tensions(),
        resumed.constraint_tensions()
    );
}
