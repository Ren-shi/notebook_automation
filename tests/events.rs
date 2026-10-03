use std::f64::consts::PI;

use physim::events::{CoordinateCrossing, RadialVelocity, Separation};
use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

#[test]
fn terminal_ground_hit_stops_at_the_exact_time() {
    let (y0, vy, g) = (1.0, 4.0, 9.81);
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(v(0.0, y0, 0.0), v(2.0, vy, 0.0), 1.0)
        .unwrap();
    w.add_force(UniformField { g: v(0.0, -g, 0.0) });
    let ground = Event::new(
        CoordinateCrossing {
            i: 0,
            axis: 1,
            value: 0.0,
        },
        Direction::Falling,
        true,
    );
    let traj = w.run_with_events(0.01, 10_000, 10, &[ground]).unwrap();

    // Verlet is exact for constant acceleration, so the event time is limited only by root finding.
    let t_hit = (vy + (vy * vy + 2.0 * g * y0).sqrt()) / g;
    assert!(
        (w.state.t - t_hit).abs() < 1e-12,
        "t = {}, expected {t_hit}",
        w.state.t
    );
    assert!(w.state.pos[0].y.abs() < 1e-10);
    assert_eq!(traj.terminated_by, Some(0));
    assert_eq!(traj.events.len(), 1);
    assert_eq!(traj.events[0].t, w.state.t);
    assert_eq!(*traj.t.last().unwrap(), w.state.t);
}

#[test]
fn kepler_periapses_and_apoapses_land_on_half_periods() {
    // Relative orbit: r = 1, v = 0.8 perpendicular, GM = 1 => starts at apoapsis.
    let mut w = World::with_integrator("yoshida4").unwrap();
    w.add_particle(v(0.5, 0.0, 0.0), v(0.0, 0.4, 0.0), 0.5)
        .unwrap();
    w.add_particle(v(-0.5, 0.0, 0.0), v(0.0, -0.4, 0.0), 0.5)
        .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    let a: f64 = 1.0 / (2.0 - 0.64);
    let e = 1.0 / a - 1.0;
    let period = 2.0 * PI * a.powf(1.5);

    let peri = Event::new(
        RadialVelocity { i: 0, j: Some(1) },
        Direction::Rising,
        false,
    );
    let apo = Event::new(
        RadialVelocity { i: 0, j: Some(1) },
        Direction::Falling,
        false,
    );
    let traj = w
        .run_with_events(1e-3, (10.0 * period / 1e-3) as usize, 1000, &[peri, apo])
        .unwrap();

    let (mut n_peri, mut n_apo) = (0, 0);
    for hit in &traj.events {
        let r = (hit.pos[0] - hit.pos[1]).norm();
        let (k, expected_r) = if hit.event == 0 {
            n_peri += 1;
            (n_peri as f64 - 0.5, a * (1.0 - e))
        } else {
            n_apo += 1;
            (n_apo as f64, a * (1.0 + e))
        };
        assert!(
            (hit.t - k * period).abs() < 1e-9,
            "event {} at {}, expected {}",
            hit.event,
            hit.t,
            k * period
        );
        assert!(
            (r - expected_r).abs() < 1e-10,
            "r = {r}, expected {expected_r}"
        );
    }
    assert_eq!((n_peri, n_apo), (10, 9)); // the run ends just before the 10th apoapsis
    assert!(traj.events.windows(2).all(|p| p[0].t < p[1].t));
}

fn oscillator() -> World {
    let mut w = World::with_integrator("rk4").unwrap();
    w.add_particle(v(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_force(AnchorSpring {
        i: 0,
        anchor: Vec3::ZERO,
        k: 1.0,
        rest_length: 0.0,
    });
    w
}

#[test]
fn both_directions_find_every_zero_crossing() {
    let mut w = oscillator();
    let zero = Event::new(
        CoordinateCrossing {
            i: 0,
            axis: 0,
            value: 0.0,
        },
        Direction::Either,
        false,
    );
    let traj = w.run_with_events(1e-3, 10_000, 10_000, &[zero]).unwrap();
    assert_eq!(traj.events.len(), 3); // pi/2, 3pi/2, 5pi/2 < 10
    for (k, hit) in traj.events.iter().enumerate() {
        assert!(
            (hit.t - (PI / 2.0 + k as f64 * PI)).abs() < 1e-10,
            "{}",
            hit.t
        );
    }
}

#[test]
fn events_in_one_step_are_ordered_and_a_terminal_event_ends_the_run() {
    let mut w = oscillator();
    let events = [
        Event::new(
            CoordinateCrossing {
                i: 0,
                axis: 0,
                value: 0.3,
            },
            Direction::Falling,
            false,
        ),
        Event::new(
            CoordinateCrossing {
                i: 0,
                axis: 0,
                value: 0.4,
            },
            Direction::Falling,
            true,
        ),
        Event::new(
            CoordinateCrossing {
                i: 0,
                axis: 0,
                value: 0.5,
            },
            Direction::Falling,
            false,
        ),
    ];
    // One big step covers all three crossings (x = cos t falls through 0.5, 0.4, 0.3).
    let traj = w.run_with_events(1.4, 1, 1, &events).unwrap();
    let fired: Vec<usize> = traj.events.iter().map(|h| h.event).collect();
    assert_eq!(fired, vec![2, 1]); // 0.5 first, then the terminal 0.4; 0.3 is never reached
    assert_eq!(traj.terminated_by, Some(1));
    assert!((w.state.pos[0].x - 0.4).abs() < 1e-10);
}

#[test]
fn starting_exactly_on_zero_is_not_an_event() {
    let mut w = oscillator();
    let at_start = Event::new(
        CoordinateCrossing {
            i: 0,
            axis: 0,
            value: 1.0,
        },
        Direction::Either,
        true,
    );
    let traj = w.run_with_events(1e-2, 100, 100, &[at_start]).unwrap();
    assert!(traj.events.is_empty());
    assert_eq!(traj.terminated_by, None);
}

#[test]
fn separation_events_detect_approach() {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(v(-5.0, 0.0, 0.0), v(1.0, 0.0, 0.0), 1.0)
        .unwrap();
    w.add_particle(v(5.0, 0.0, 0.0), v(-1.0, 0.0, 0.0), 1.0)
        .unwrap();
    let contact = Event::new(
        Separation {
            i: 0,
            j: 1,
            distance: 1.0,
        },
        Direction::Falling,
        true,
    );
    w.run_with_events(0.1, 1000, 1000, &[contact]).unwrap();
    assert!((w.state.t - 4.5).abs() < 1e-12);
}

#[test]
fn bad_event_functions_fail_cleanly() {
    let mut w = oscillator();
    let bad = Event::new(
        CoordinateCrossing {
            i: 3,
            axis: 0,
            value: 0.0,
        },
        Direction::Either,
        false,
    );
    let failure = w.run_with_events(1e-2, 10, 1, &[bad]).unwrap_err();
    assert!(failure.error.to_string().contains("out of range"));
    assert_eq!(w.state.t, 0.0);
    assert_eq!(failure.trajectory.n_frames(), 1);
}
