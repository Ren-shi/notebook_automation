use std::f64::consts::TAU;

use physim::integrators::{Composition, Op, Scheme, Splitting};
use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

const G_SOLAR: f64 = 2.95912208286e-4;

/// Sun (with the inner planets' mass) and Jupiter to Pluto; AU, days, solar masses.
/// Hairer, Lubich & Wanner, Geometric Numerical Integration, §I.2.4.
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
    for (m, x, u) in bodies {
        w.add_particle(Vec3::from(x), Vec3::from(u), m).unwrap();
    }
    w.add_force(NewtonianGravity {
        g: G_SOLAR,
        softening: 0.0,
    });
    w
}

fn energy_errors(traj: &Trajectory) -> Vec<f64> {
    let e = traj.total_energy();
    e.iter().map(|x| ((x - e[0]) / e[0]).abs()).collect()
}

#[test]
fn wisdom_holman_keeps_the_outer_solar_system_bounded_for_a_million_years() {
    let dt = 4332.59 / 20.0; // 1/20 of Jupiter's period, in days
    let steps = (1e6 * 365.25 / dt) as usize;
    let mut w = outer_solar_system("wisdom_holman");
    let traj = w
        .run_with_options(
            dt,
            steps,
            &RunOptions {
                record_every: steps / 1000,
                events: &[],
                energies: true,
            },
        )
        .unwrap();
    let err = energy_errors(&traj);
    let tenth = err.len() / 10;
    let max = |s: &[f64]| s.iter().copied().fold(0.0, f64::max);
    let (first, last) = (max(&err[..tenth]), max(&err[err.len() - tenth..]));
    assert!(max(&err) < 1e-5, "max energy error {:e}", max(&err));
    // Bounded: no secular growth between the first and last 100 000 years.
    assert!(last < 1.5 * first, "first {first:e}, last {last:e}");

    // Verlet at the same step is ~500 times worse.
    let mut w = outer_solar_system("verlet");
    let traj = w.run(dt, steps / 10, steps / 1000).unwrap();
    assert!(max(&energy_errors(&traj)) > 100.0 * max(&err));
}

/// Two planets on eccentric, inclined orbits around a unit-mass star (G = 1).
fn two_planets(integrator: &str, planet_mass: f64) -> World {
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.add_particle(v(1.0, 0.0, 0.0), v(0.0, 1.1, 0.05), planet_mass)
        .unwrap();
    w.add_particle(v(-1.6, 0.3, 0.1), v(0.1, -0.75, 0.0), planet_mass)
        .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    w
}

fn reference(planet_mass: f64, t: f64) -> State {
    let mut w = two_planets("verlet", planet_mass);
    w.run_adaptive(
        t,
        &AdaptiveOptions {
            rtol: 1e-13,
            atol: 1e-16,
            output: Output::Times(&[]),
            ..Default::default()
        },
    )
    .unwrap();
    w.state
}

fn position_error(a: &State, b: &State) -> f64 {
    a.pos
        .iter()
        .zip(&b.pos)
        .map(|(x, y)| (*x - *y).norm())
        .fold(0.0, f64::max)
}

#[test]
fn wisdom_holman_is_second_order_with_error_proportional_to_the_mass_ratio() {
    let t = 10.0 * TAU;
    let err = |mass: f64, dt: f64| {
        let mut w = two_planets("wisdom_holman", mass);
        let steps = (t / dt).round() as usize;
        w.run(t / steps as f64, steps, steps).unwrap();
        position_error(&w.state, &reference(mass, t))
    };
    let (e1, e2) = (err(1e-3, 0.05), err(1e-3, 0.025));
    let order = (e1 / e2).log2();
    assert!((order - 2.0).abs() < 0.2, "order {order}: {e1:e} {e2:e}");
    // Ten times lighter planets: ten times smaller error at the same step, once the
    // planets are light enough for the O(ε²) terms to be negligible.
    let ratio = err(1e-4, 0.05) / err(1e-5, 0.05);
    assert!((ratio / 10.0 - 1.0).abs() < 0.2, "mass scaling {ratio}");
    // Verlet at the same step is far less accurate.
    let mut w = two_planets("verlet", 1e-3);
    let steps = (t / 0.05).round() as usize;
    w.run(t / steps as f64, steps, steps).unwrap();
    assert!(position_error(&w.state, &reference(1e-3, t)) > 100.0 * e1);
}

#[test]
fn wisdom_holman_handles_perturbations_restarts_and_bad_input() {
    // J2 on a free central body acts as a perturbation in the interaction kick:
    // the node regresses at the secular rate.
    let (a, incl, j2, radius): (f64, f64, f64, f64) = (2.0, 0.5, 1e-3, 1.0);
    let mut w = World::with_integrator("wisdom_holman").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    let speed = (1.0 / a).sqrt();
    w.add_particle(
        v(a, 0.0, 0.0),
        v(0.0, speed * incl.cos(), speed * incl.sin()),
        1e-9,
    )
    .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    w.add_force(J2Oblateness {
        central: 0,
        g: 1.0,
        j2,
        radius,
        axis: v(0.0, 0.0, 1.0),
    });
    let n = a.powf(-1.5);
    let period = TAU / n;
    let steps_per_orbit = 40;
    let orbits = 60;
    let traj = w
        .run(period / steps_per_orbit as f64, orbits * steps_per_orbit, 4)
        .unwrap();
    let node: Vec<f64> = (0..traj.n_frames())
        .map(|k| {
            let rel = traj.pos[2 * k + 1] - traj.pos[2 * k];
            let vel = traj.vel[2 * k + 1] - traj.vel[2 * k];
            let h = rel.cross(vel);
            h.x.atan2(-h.y)
        })
        .collect();
    let (tm, nm) = (
        traj.t.iter().sum::<f64>() / traj.t.len() as f64,
        node.iter().sum::<f64>() / node.len() as f64,
    );
    let slope = traj
        .t
        .iter()
        .zip(&node)
        .map(|(t, o)| (t - tm) * (o - nm))
        .sum::<f64>()
        / traj.t.iter().map(|t| (t - tm).powi(2)).sum::<f64>();
    let expected = -1.5 * n * j2 * (radius / a).powi(2) * incl.cos();
    assert!(
        (slope / expected - 1.0).abs() < 5e-3,
        "dΩ/dt = {slope:e} vs {expected:e}"
    );

    // Restarts are bit-exact.
    let mut straight = two_planets("wisdom_holman", 1e-3);
    let mut interrupted = two_planets("wisdom_holman", 1e-3);
    straight.run(0.05, 300, 300).unwrap();
    interrupted.run(0.05, 120, 120).unwrap();
    let mut resumed = World::from_checkpoint(interrupted.checkpoint(), |_, _| {
        unreachable!("no external forces")
    })
    .unwrap();
    resumed.run(0.05, 180, 180).unwrap();
    assert_eq!(straight.state.pos, resumed.state.pos);
    assert_eq!(straight.state.vel, resumed.state.vel);

    // Errors.
    let mut w = two_planets("wisdom_holman", 1e-3);
    w.forces.clear();
    assert!(w
        .step(0.1)
        .unwrap_err()
        .to_string()
        .contains("NewtonianGravity"));
    let mut w = two_planets("wisdom_holman", 1e-3);
    w.add_force(LinearDrag { gamma: 0.1 });
    assert!(w.step(0.1).unwrap_err().to_string().contains("velocity"));
    let mut w = two_planets("wisdom_holman", 1e-3);
    w.forces.clear();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.1,
    });
    assert!(w.step(0.1).unwrap_err().to_string().contains("softening"));
}

fn oscillator(integrator: Box<dyn Integrator>) -> World {
    let mut w = World::new(integrator);
    w.add_particle(v(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_force(AnchorSpring {
        i: 0,
        anchor: Vec3::ZERO,
        k: 1.0,
        rest_length: 0.0,
    });
    w
}

fn oscillator_error(w: &mut World, dt: f64, t_end: f64) -> f64 {
    let steps = (t_end / dt).round() as usize;
    w.run(dt, steps, steps).unwrap();
    let (x, u) = (w.state.pos[0].x, w.state.vel[0].x);
    (x - t_end.cos()).hypot(u + t_end.sin())
}

#[test]
fn user_defined_schemes() {
    // The triple jump as a user composition reproduces yoshida4 bit for bit.
    let c = 2f64.cbrt();
    let w1 = 1.0 / (2.0 - c);
    let custom = Scheme::Composition {
        name: "triple_jump".into(),
        order: 4,
        weights: vec![w1, -c / (2.0 - c), w1],
    };
    let mut a = oscillator(custom.build().unwrap());
    let mut b = oscillator(integrators::by_name("yoshida4").unwrap());
    a.run(1e-2, 500, 500).unwrap();
    b.run(1e-2, 500, 500).unwrap();
    assert_eq!(a.state.pos, b.state.pos);
    assert_eq!(a.integrator().name(), "triple_jump");

    // Leapfrog written as a splitting is second order; it survives a checkpoint.
    let leapfrog = Scheme::Splitting {
        name: "my_leapfrog".into(),
        order: 2,
        ops: vec![Op::Kick(0.5), Op::Drift(1.0), Op::Kick(0.5)],
    };
    let mut w = oscillator(leapfrog.build().unwrap());
    let e1 = oscillator_error(&mut w, 0.02, 2.0);
    let mut w = oscillator(leapfrog.build().unwrap());
    let e2 = oscillator_error(&mut w, 0.01, 2.0);
    assert!(((e1 / e2).log2() - 2.0).abs() < 0.1);
    let cp = w.checkpoint();
    assert_eq!(cp.integrator_scheme, Some(leapfrog.clone()));
    let mut r = World::from_checkpoint(cp, |_, _| unreachable!()).unwrap();
    w.run(0.01, 10, 10).unwrap();
    r.run(0.01, 10, 10).unwrap();
    assert_eq!(w.state.pos, r.state.pos);
    assert_eq!(r.integrator().name(), "my_leapfrog");

    // Built-in schemes are rebuilt by name.
    assert_eq!(Composition::yoshida6().scheme(), None);
    assert_eq!(Splitting::pefrl().scheme(), None);

    // Validation.
    for bad in [
        Scheme::Composition {
            name: "x".into(),
            order: 2,
            weights: vec![0.5, 0.4],
        },
        Scheme::Splitting {
            name: "x".into(),
            order: 2,
            ops: vec![Op::Kick(1.0), Op::Drift(0.9)],
        },
        Scheme::Composition {
            name: "verlet".into(),
            order: 2,
            weights: vec![1.0],
        },
    ] {
        assert!(bad.build().is_err());
    }
}

#[test]
fn optimised_fourth_order_schemes_beat_yoshida_at_equal_cost() {
    // Same number of force evaluations: yoshida4 (3 per step) at dt, pefrl (4) at 4/3 dt,
    // blanes_moan4 (6) at 2 dt.
    let dt = 0.01;
    let mut y = oscillator(integrators::by_name("yoshida4").unwrap());
    let ey = oscillator_error(&mut y, dt, 4.0);
    let mut p = oscillator(integrators::by_name("pefrl").unwrap());
    let ep = oscillator_error(&mut p, 4.0 / 300.0, 4.0);
    let mut b = oscillator(integrators::by_name("blanes_moan4").unwrap());
    let eb = oscillator_error(&mut b, 2.0 * dt, 4.0);
    assert!(ep < ey / 50.0, "pefrl {ep:e} vs yoshida4 {ey:e}");
    assert!(eb < ey / 150.0, "blanes_moan4 {eb:e} vs yoshida4 {ey:e}");
}

#[test]
fn gauss_legendre_conserves_quadratic_invariants_and_is_reversible() {
    // a = v × B: a velocity-dependent force that conserves |v|², a quadratic invariant.
    let b = v(0.3, -0.2, 1.0);
    let gyro = || ClosureForce::per_particle("Gyro", move |_t, _r, u, _m| u.cross(b));
    for name in ["gauss2", "gauss4", "gauss6", "rk4"] {
        let mut w = World::with_integrator(name).unwrap();
        w.add_particle(Vec3::ZERO, v(1.0, 0.5, -0.3), 1.0).unwrap();
        w.add_force(gyro());
        let speed2 = w.state.vel[0].norm_squared();
        w.run(0.2, 5000, 5000).unwrap();
        let drift = (w.state.vel[0].norm_squared() / speed2 - 1.0).abs();
        if name == "rk4" {
            assert!(drift > 1e-6, "rk4 should drift: {drift:e}");
        } else {
            assert!(drift < 1e-13, "{name}: |v|² drifted by {drift:e}");
            // Symmetric: stepping back returns to the start.
            w.run(-0.2, 5000, 5000).unwrap();
            assert!(
                w.state.pos[0].norm() < 1e-10,
                "{name}: {:?}",
                w.state.pos[0]
            );
        }
    }

    // A step too large for the fixed-point iteration fails cleanly.
    let mut w = World::with_integrator("gauss4").unwrap();
    w.add_particle(v(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_force(AnchorSpring {
        i: 0,
        anchor: Vec3::ZERO,
        k: 1.0,
        rest_length: 0.0,
    });
    assert!(w.step(5.0).unwrap_err().to_string().contains("converge"));
    assert_eq!(w.state.t, 0.0);
}

#[test]
fn new_symplectic_schemes_bound_the_energy_error() {
    for name in [
        "yoshida6",
        "yoshida8",
        "pefrl",
        "blanes_moan4",
        "gauss4",
        "gauss6",
    ] {
        let mut w = World::with_integrator(name).unwrap();
        w.add_particle(v(0.5, 0.0, 0.0), v(0.0, 0.4, 0.0), 0.5)
            .unwrap();
        w.add_particle(v(-0.5, 0.0, 0.0), v(0.0, -0.4, 0.0), 0.5)
            .unwrap();
        w.add_force(NewtonianGravity {
            g: 1.0,
            softening: 0.0,
        });
        let traj = w.run(0.01, 20_000, 20).unwrap();
        let err = energy_errors(&traj);
        let half = err.len() / 2;
        let max = |s: &[f64]| s.iter().copied().fold(0.0, f64::max);
        let (first, second) = (max(&err[..half]), max(&err[half..]));
        assert!(max(&err) < 1e-5, "{name}: {:e}", max(&err));
        assert!(
            second < 1.5 * first + 1e-14,
            "{name}: {first:e} then {second:e}"
        );
    }
}
