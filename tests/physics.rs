use physim::integrators::{self, NAMES};
use physim::*;

/// Unit harmonic oscillator x(0)=1, v(0)=0 => x(t) = cos t.
fn oscillator_error(integrator: &str, dt: f64) -> f64 {
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0)
        .unwrap();
    w.add_force(AnchorSpring {
        i: 0,
        anchor: Vec3::ZERO,
        k: 1.0,
        rest_length: 0.0,
    });
    let t_end = 2.0;
    let steps = (t_end / dt).round() as usize;
    w.run(dt, steps, steps).unwrap();
    let x = w.state.pos[0];
    let v = w.state.vel[0];
    ((x.x - t_end.cos()).powi(2) + (v.x + t_end.sin()).powi(2)).sqrt()
}

#[test]
fn integrators_converge_at_their_stated_order() {
    for &name in NAMES {
        let order = integrators::by_name(name).unwrap().order() as f64;
        let (e1, e2) = (oscillator_error(name, 0.02), oscillator_error(name, 0.01));
        let measured = (e1 / e2).log2();
        assert!(
            (measured - order).abs() < 0.25,
            "{name}: expected order {order}, measured {measured:.3}"
        );
    }
}

/// Eccentric two-body orbit in the centre-of-mass frame with G = 1.
fn kepler_world(integrator: &str) -> World {
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(Vec3::new(0.5, 0.0, 0.0), Vec3::new(0.0, 0.4, 0.0), 0.5)
        .unwrap();
    w.add_particle(Vec3::new(-0.5, 0.0, 0.0), Vec3::new(0.0, -0.4, 0.0), 0.5)
        .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    w
}

fn max_relative_energy_error(traj: &Trajectory) -> f64 {
    let e = traj.total_energy();
    e.iter()
        .map(|x| ((x - e[0]) / e[0]).abs())
        .fold(0.0, f64::max)
}

#[test]
fn symplectic_integrators_bound_energy_error_on_kepler_orbit() {
    // ~50 orbits.
    let (dt, steps) = (1e-3, 200_000);
    let verlet = max_relative_energy_error(&kepler_world("verlet").run(dt, steps, 100).unwrap());
    let yoshida = max_relative_energy_error(&kepler_world("yoshida4").run(dt, steps, 100).unwrap());
    let euler =
        max_relative_energy_error(&kepler_world("explicit_euler").run(dt, steps, 100).unwrap());
    assert!(verlet < 1e-4, "verlet energy error {verlet:e}");
    assert!(yoshida < 1e-8, "yoshida4 energy error {yoshida:e}");
    assert!(
        euler > 100.0 * verlet,
        "explicit Euler should drift: {euler:e} vs {verlet:e}"
    );
}

#[test]
fn leapfrog_conserves_momentum_and_angular_momentum_in_nbody() {
    let mut w = World::with_integrator("verlet").unwrap();
    let bodies = [
        ([1.0, 0.0, 0.1], [0.0, 0.5, 0.0], 1.0),
        ([-0.6, 0.8, 0.0], [-0.3, -0.2, 0.1], 0.7),
        ([-0.3, -0.9, -0.2], [0.4, -0.1, 0.0], 1.3),
        ([0.2, 0.3, 1.1], [0.0, 0.0, -0.2], 0.4),
    ];
    for (r, v, m) in bodies {
        w.add_particle(r.into(), v.into(), m).unwrap();
    }
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.05,
    });
    let (p0, l0) = (w.state.momentum(), w.state.angular_momentum());
    w.run(1e-3, 20_000, 1000).unwrap();
    assert!((w.state.momentum() - p0).norm() < 1e-12);
    assert!((w.state.angular_momentum() - l0).norm() < 1e-11);
}

#[test]
fn verlet_is_exact_for_constant_acceleration() {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::new(3.0, 4.0, 0.0), 2.0)
        .unwrap();
    w.add_force(UniformField {
        g: Vec3::new(0.0, -9.81, 0.0),
    });
    w.run(0.01, 50, 50).unwrap();
    let t = w.state.t;
    let expected = Vec3::new(3.0 * t, 4.0 * t - 0.5 * 9.81 * t * t, 0.0);
    assert!((w.state.pos[0] - expected).norm() < 1e-12);
}

#[test]
fn drag_dissipates_energy() {
    let mut w = World::with_integrator("rk4").unwrap();
    w.add_particle(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0)
        .unwrap();
    w.add_force(AnchorSpring {
        i: 0,
        anchor: Vec3::ZERO,
        k: 1.0,
        rest_length: 0.0,
    });
    w.add_force(LinearDrag { gamma: 0.2 });
    let e = w.run(0.01, 1000, 1).unwrap().total_energy();
    assert!(e.windows(2).all(|p| p[1] <= p[0] + 1e-15));
    // Weakly damped oscillator: E(t) ~ E0 exp(-gamma t), averaged over a period.
    let ratio = e.last().unwrap() / e[0];
    assert!(
        (ratio - (-0.2f64 * 10.0).exp()).abs() < 0.05,
        "energy ratio {ratio}"
    );
}

#[test]
fn invalid_input_is_rejected() {
    assert!(World::with_integrator("nope").is_err());
    let mut w = World::with_integrator("verlet").unwrap();
    assert!(w.add_particle(Vec3::ZERO, Vec3::ZERO, -1.0).is_err());
    assert!(w.add_particle(Vec3::ZERO, Vec3::ZERO, f64::NAN).is_err());
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.add_force(Spring {
        i: 0,
        j: 5,
        k: 1.0,
        rest_length: 1.0,
    });
    assert!(w.step(0.1).is_err());
    assert!(w.run(0.1, 10, 0).is_err());
}
