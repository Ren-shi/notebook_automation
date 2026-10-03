use std::f64::consts::{PI, TAU};

use physim::events::RadialVelocity;
use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Five particles with different masses, none at a symmetric point.
fn cloud() -> State {
    let mut s = State::new();
    let pts = [
        (v(1.0, 0.2, -0.3), 1.0),
        (v(-0.7, 1.1, 0.4), 2.0),
        (v(0.3, -0.9, 1.2), 0.5),
        (v(1.6, 1.4, 0.1), 1.5),
        (v(-1.2, -0.4, -0.8), 0.8),
    ];
    for (p, m) in pts {
        s.add_particle(p, Vec3::ZERO, m);
    }
    s
}

/// m_p a_p = -∂U/∂r_p for every particle, by central differences.
fn assert_gradient(force: &dyn Force, t: f64) {
    let s = cloud();
    let n = s.len();
    let mut acc = vec![Vec3::ZERO; n];
    force
        .accumulate(t, &s.pos, &s.vel, &s.mass, &mut acc)
        .unwrap();
    let u = |pos: &[Vec3]| force.potential(t, pos, &s.mass).unwrap().unwrap();
    let h = 1e-6;
    let units = [v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)];
    for (p, a) in acc.iter().enumerate() {
        for (axis, e) in units.iter().enumerate() {
            let mut plus = s.pos.clone();
            let mut minus = s.pos.clone();
            plus[p] += *e * h;
            minus[p] -= *e * h;
            let grad = (u(&plus) - u(&minus)) / (2.0 * h);
            let f = s.mass[p] * a.dot(*e);
            let scale = grad.abs().max(1e-3);
            assert!(
                (f + grad).abs() < 1e-6 * scale,
                "{}: particle {p} axis {axis}: m a = {f}, -dU/dx = {}",
                force.name(),
                -grad
            );
        }
    }
}

#[test]
fn conservative_forces_are_minus_the_gradient_of_their_potential() {
    let c = v(0.1, -0.2, 0.05);
    let forces: Vec<Box<dyn Force>> = vec![
        Box::new(PowerLaw {
            center: c,
            k: -1.3,
            n: -1.0,
        }),
        Box::new(PowerLaw {
            center: c,
            k: 0.7,
            n: 0.0,
        }),
        Box::new(PowerLaw {
            center: c,
            k: 0.4,
            n: 2.5,
        }),
        Box::new(Yukawa {
            center: c,
            k: 1.2,
            length: 0.8,
        }),
        Box::new(PlummerPotential {
            center: c,
            gm: 2.0,
            a: 0.5,
        }),
        Box::new(HernquistPotential {
            center: c,
            gm: 2.0,
            a: 0.5,
        }),
        Box::new(HarmonicTrap {
            center: c,
            omega: v(1.0, 2.0, 0.5),
        }),
        Box::new(DampedSpring {
            i: 0,
            j: 2,
            k: 3.0,
            rest_length: 0.5,
            c: 0.7,
        }),
        Box::new(ModulatedSpring {
            i: 1,
            to: Anchor::Particle(3),
            k: 2.0,
            depth: 0.3,
            omega: 1.7,
            phase: 0.2,
            rest_length: 0.4,
        }),
        Box::new(ModulatedSpring {
            i: 4,
            to: Anchor::Point(c),
            k: 2.0,
            depth: 0.3,
            omega: 1.7,
            phase: 0.2,
            rest_length: 0.0,
        }),
        Box::new(
            SpringNetwork::new(
                vec![0, 1, 2, 3],
                vec![1, 2, 3, 4],
                vec![1.0, 2.0, 3.0, 4.0],
                vec![0.5, 0.0, 1.0, 2.0],
            )
            .unwrap(),
        ),
        Box::new(J2Oblateness {
            central: 1,
            g: 1.0,
            j2: 0.01,
            radius: 0.5,
            axis: v(0.2, 0.1, 1.0),
        }),
    ];
    for f in &forces {
        assert_gradient(f.as_ref(), 0.9);
    }
}

#[test]
fn spring_network_matches_separate_springs_bit_for_bit() {
    // A 10x10 square lattice with diagonals, perturbed.
    let side = 10;
    let id = |x: usize, y: usize| y * side + x;
    let mut bonds = Vec::new();
    for y in 0..side {
        for x in 0..side {
            if x + 1 < side {
                bonds.push((id(x, y), id(x + 1, y), 1.0));
            }
            if y + 1 < side {
                bonds.push((id(x, y), id(x, y + 1), 1.0));
            }
            if x + 1 < side && y + 1 < side {
                bonds.push((id(x, y), id(x + 1, y + 1), 2f64.sqrt()));
            }
        }
    }
    let build = |network: bool| {
        let mut w = World::with_integrator("verlet").unwrap();
        for y in 0..side {
            for x in 0..side {
                let jitter = 0.01 * ((x * 7 + y * 13) % 5) as f64;
                w.add_particle(v(x as f64 + jitter, y as f64, 0.0), Vec3::ZERO, 1.0)
                    .unwrap();
            }
        }
        if network {
            let (i, j, l): (Vec<_>, Vec<_>, Vec<_>) = bonds.iter().fold(
                (vec![], vec![], vec![]),
                |(mut i, mut j, mut l), &(a, b, r)| {
                    i.push(a);
                    j.push(b);
                    l.push(r);
                    (i, j, l)
                },
            );
            let k = vec![50.0; i.len()];
            w.add_force(SpringNetwork::new(i, j, k, l).unwrap());
        } else {
            for &(i, j, rest_length) in &bonds {
                w.add_force(Spring {
                    i,
                    j,
                    k: 50.0,
                    rest_length,
                });
            }
        }
        w
    };
    let (mut a, mut b) = (build(true), build(false));
    let ta = a.run(1e-3, 500, 100).unwrap();
    let tb = b.run(1e-3, 500, 100).unwrap();
    assert_eq!(ta.pos, tb.pos);
    assert_eq!(ta.potential, tb.potential);
}

/// Projects `x(t)` sampled over whole periods onto cos and sin of `omega t`.
fn fourier(t: &[f64], x: &[f64], omega: f64) -> (f64, f64) {
    let n = t.len() - 1; // the last sample repeats the first phase
    let (mut a, mut b) = (0.0, 0.0);
    for k in 0..n {
        a += x[k] * (omega * t[k]).cos();
        b += x[k] * (omega * t[k]).sin();
    }
    (2.0 * a / n as f64, 2.0 * b / n as f64)
}

#[test]
fn driven_damped_oscillator_follows_the_resonance_curve() {
    // m x'' = -k (x - L) - c x' + F0 cos(ω t), built from a DampedSpring to a pinned particle.
    let (m, k, c, f0, l): (f64, f64, f64, f64, f64) = (1.0, 1.0, 0.2, 0.1, 2.0);
    let (w0, gamma) = ((k / m).sqrt(), c / m);
    for omega in [0.3, 0.8, 0.95, 1.0, 1.05, 1.5, 3.0] {
        let mut w = World::with_integrator("verlet").unwrap();
        w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
        w.pin(0, true).unwrap();
        w.add_particle(v(l, 0.0, 0.0), Vec3::ZERO, m).unwrap();
        w.add_force(DampedSpring {
            i: 0,
            j: 1,
            k,
            rest_length: l,
            c,
        });
        w.add_force(PeriodicForce {
            i: 1,
            amplitude: v(f0, 0.0, 0.0),
            omega,
            phase: 0.0,
        });
        // Transients decay as exp(-γ t / 2): wait 200 time units, then sample 10 periods.
        let period = TAU / omega;
        let t0 = (200.0 / period).ceil() * period;
        let samples = 2000;
        let times: Vec<f64> = (0..=samples)
            .map(|s| t0 + 10.0 * period * s as f64 / samples as f64)
            .collect();
        let run = w
            .run_adaptive(
                *times.last().unwrap(),
                &AdaptiveOptions {
                    rtol: 1e-11,
                    atol: 1e-14,
                    output: Output::Times(&times),
                    ..Default::default()
                },
            )
            .unwrap();
        let x: Vec<f64> = run.trajectory.pos.chunks(2).map(|p| p[1].x - l).collect();
        let (a, b) = fourier(&times, &x, omega);
        let amplitude = a.hypot(b);
        let phase = b.atan2(a); // x = A cos(ω t - δ)
        let denom = ((w0 * w0 - omega * omega).powi(2) + (gamma * omega).powi(2)).sqrt();
        let expected = f0 / m / denom;
        let delta = (gamma * omega).atan2(w0 * w0 - omega * omega);
        assert!(
            (amplitude / expected - 1.0).abs() < 1e-8,
            "ω = {omega}: amplitude {amplitude} vs {expected}"
        );
        assert!(
            (phase - delta).abs() < 1e-6,
            "ω = {omega}: phase {phase} vs {delta}"
        );
    }
}

/// Unit-mass central body pinned at the origin, light planet at apoapsis of an orbit with
/// semi-major axis `a` and eccentricity `e`, inclined by `incl` about the x axis.
fn planet(a: f64, e: f64, incl: f64) -> World {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.pin(0, true).unwrap();
    let speed = ((1.0 - e) / (a * (1.0 + e))).sqrt();
    w.add_particle(
        v(a * (1.0 + e), 0.0, 0.0),
        v(0.0, speed * incl.cos(), speed * incl.sin()),
        1e-9,
    )
    .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    w
}

#[test]
fn post_newtonian_periapsis_advance_matches_general_relativity() {
    let (a, e, c) = (1.0, 0.5, 100.0);
    let mut w = planet(a, e, 0.0);
    w.add_force(PostNewtonian {
        central: 0,
        g: 1.0,
        c,
    });
    let events = [Event::new(
        RadialVelocity { i: 1, j: Some(0) },
        Direction::Rising,
        false,
    )];
    let orbits = 20;
    let run = w
        .run_adaptive(
            orbits as f64 * TAU * a.powf(1.5),
            &AdaptiveOptions {
                rtol: 1e-12,
                atol: 1e-15,
                events: &events,
                output: Output::Times(&[]),
                ..Default::default()
            },
        )
        .unwrap();
    let angles: Vec<f64> = run
        .trajectory
        .events
        .iter()
        .map(|h| h.pos[1].y.atan2(h.pos[1].x))
        .collect();
    assert_eq!(angles.len(), orbits);
    let mut advance: Vec<f64> = angles.windows(2).map(|p| p[1] - p[0]).collect();
    for d in &mut advance {
        *d = (*d + PI).rem_euclid(TAU) - PI;
    }
    let measured = advance.iter().sum::<f64>() / advance.len() as f64;
    let expected = 6.0 * PI / (c * c * a * (1.0 - e * e));
    // 1PN is the leading order: corrections are O(GM / (c² a)) = 1e-4 relative.
    assert!(
        (measured / expected - 1.0).abs() < 1e-3,
        "advance per orbit {measured:e} vs {expected:e}"
    );
}

#[test]
fn j2_nodal_regression_matches_the_secular_rate() {
    let (a, incl, j2, radius) = (2.0, 0.5, 1e-3, 1.0);
    let mut w = planet(a, 0.0, incl);
    w.add_force(J2Oblateness {
        central: 0,
        g: 1.0,
        j2,
        radius,
        axis: v(0.0, 0.0, 1.0),
    });
    let n = a.powf(-1.5);
    let period = TAU / n;
    let t_end = 60.0 * period;
    let times: Vec<f64> = (0..=600).map(|k| t_end * k as f64 / 600.0).collect();
    let run = w
        .run_adaptive(
            t_end,
            &AdaptiveOptions {
                rtol: 1e-11,
                atol: 1e-14,
                output: Output::Times(&times),
                ..Default::default()
            },
        )
        .unwrap();
    let traj = &run.trajectory;
    // Longitude of the ascending node from the angular momentum vector.
    let node: Vec<f64> = (0..times.len())
        .map(|k| {
            let h = traj.pos[2 * k + 1].cross(traj.vel[2 * k + 1]);
            h.x.atan2(-h.y)
        })
        .collect();
    let (tm, nm) = (
        times.iter().sum::<f64>() / times.len() as f64,
        node.iter().sum::<f64>() / node.len() as f64,
    );
    let slope = times
        .iter()
        .zip(&node)
        .map(|(t, o)| (t - tm) * (o - nm))
        .sum::<f64>()
        / times.iter().map(|t| (t - tm).powi(2)).sum::<f64>();
    let expected = -1.5 * n * j2 * (radius / a).powi(2) * incl.cos();
    assert!(
        (slope / expected - 1.0).abs() < 3e-3,
        "dΩ/dt = {slope:e} vs {expected:e}"
    );
}

#[test]
fn j2_conserves_momentum_with_a_free_central_body() {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.add_particle(v(2.0, 0.0, 0.3), v(0.0, 0.7, 0.1), 0.01)
        .unwrap();
    w.add_particle(v(-1.5, 0.5, -0.2), v(0.1, -0.8, 0.0), 0.02)
        .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    w.add_force(J2Oblateness {
        central: 0,
        g: 1.0,
        j2: 0.05,
        radius: 0.5,
        axis: v(0.0, 0.3, 1.0),
    });
    let p0 = w.state.momentum();
    let traj = w.run(1e-3, 5000, 5000).unwrap();
    assert!((w.state.momentum() - p0).norm() < 1e-14);
    let e = traj.total_energy();
    assert!(((e[1] - e[0]) / e[0]).abs() < 1e-6);
}

#[test]
fn modulated_spring_grows_at_the_mathieu_rate() {
    // x'' + ω0² (1 + ε cos 2ω0 t) x = 0 is unstable at the principal parametric resonance,
    // growing as exp(μ t) with μ = ε ω0 / 4 to first order in ε.
    let (w0, eps) = (1.0, 0.04);
    let monodromy = |x0: f64, v0: f64| {
        let mut w = World::with_integrator("yoshida4").unwrap();
        w.add_particle(v(x0, 0.0, 0.0), v(v0, 0.0, 0.0), 1.0)
            .unwrap();
        w.add_force(ModulatedSpring {
            i: 0,
            to: Anchor::Point(Vec3::ZERO),
            k: w0 * w0,
            depth: eps,
            omega: 2.0 * w0,
            phase: 0.0,
            rest_length: 0.0,
        });
        let period = PI / w0;
        let steps = 2000;
        w.run(period / steps as f64, steps, steps).unwrap();
        (w.state.pos[0].x, w.state.vel[0].x)
    };
    let (a, c) = monodromy(1.0, 0.0);
    let (b, d) = monodromy(0.0, 1.0);
    // Largest Floquet multiplier of [[a, b], [c, d]] (determinant 1).
    let tr = a + d;
    let lambda = 0.5 * (tr.abs() + (tr * tr - 4.0).max(0.0).sqrt());
    let mu = lambda.ln() / (PI / w0);
    let expected = eps * w0 / 4.0;
    assert!((a * d - b * c - 1.0).abs() < 1e-10);
    assert!((mu / expected - 1.0).abs() < 1e-3, "μ = {mu} vs {expected}");
}

#[test]
fn central_potentials_have_the_expected_circular_orbits() {
    // Circular speed v_c² = r Φ'(r): the orbit radius must stay constant.
    let r = 1.3;
    let cases: Vec<(Box<dyn Force>, f64)> = vec![
        (
            Box::new(PlummerPotential {
                center: Vec3::ZERO,
                gm: 2.0,
                a: 0.6,
            }),
            2.0 * r * r / (r * r + 0.36f64).powf(1.5),
        ),
        (
            Box::new(HernquistPotential {
                center: Vec3::ZERO,
                gm: 2.0,
                a: 0.6,
            }),
            2.0 * r / (r + 0.6f64).powi(2),
        ),
        (
            Box::new(PowerLaw {
                center: Vec3::ZERO,
                k: 0.8,
                n: 0.0,
            }),
            0.8,
        ),
        (
            Box::new(Yukawa {
                center: Vec3::ZERO,
                k: 1.0,
                length: 2.0,
            }),
            (-r / 2.0f64).exp() * (1.0 / r + 0.5),
        ),
    ];
    for (force, vc2) in cases {
        let name = force.name();
        let mut w = World::with_integrator("yoshida4").unwrap();
        w.add_particle(v(r, 0.0, 0.0), v(0.0, vc2.sqrt(), 0.0), 1.0)
            .unwrap();
        w.forces.add(force);
        let traj = w.run(1e-3, 20_000, 100).unwrap();
        let worst = traj
            .pos
            .iter()
            .map(|p| (p.norm() - r).abs())
            .fold(0.0, f64::max);
        assert!(worst < 1e-9, "{name}: radius drifts by {worst:e}");
    }

    // An anisotropic trap separates into independent oscillators.
    let mut w = World::with_integrator("yoshida4").unwrap();
    w.add_particle(v(1.0, 1.0, 1.0), Vec3::ZERO, 3.0).unwrap();
    w.add_force(HarmonicTrap {
        center: Vec3::ZERO,
        omega: v(1.0, 2.0, 3.0),
    });
    w.run(1e-3, 3000, 3000).unwrap();
    let p = w.state.pos[0];
    assert!((p - v(3f64.cos(), 6f64.cos(), 9f64.cos())).norm() < 1e-10);

    // n = -1, k = -GM is a Kepler centre.
    let mut a = planet(1.0, 0.6, 0.0);
    let mut b = planet(1.0, 0.6, 0.0);
    b.forces.clear();
    b.add_force(PowerLaw {
        center: Vec3::ZERO,
        k: -1.0,
        n: -1.0,
    });
    a.run(1e-3, 6283, 6283).unwrap();
    b.run(1e-3, 6283, 6283).unwrap();
    assert!((a.state.pos[1] - b.state.pos[1]).norm() < 1e-10);
}

#[test]
fn closure_forces() {
    let mut a = World::with_integrator("verlet").unwrap();
    let mut b = World::with_integrator("verlet").unwrap();
    for w in [&mut a, &mut b] {
        w.add_particle(v(0.0, 1.0, 0.0), v(1.0, 0.0, 0.0), 2.0)
            .unwrap();
    }
    a.add_force(UniformField {
        g: v(0.0, -9.81, 0.0),
    });
    b.add_force(
        ClosureForce::per_particle("Down", |_t, _r, _v, _m| v(0.0, -9.81, 0.0))
            .with_potential(|_t, pos, mass| {
                Ok(pos.iter().zip(mass).map(|(r, m)| 9.81 * m * r.y).sum())
            })
            .velocity_independent(),
    );
    let ta = a.run(1e-3, 1000, 100).unwrap();
    let tb = b.run(1e-3, 1000, 100).unwrap();
    assert_eq!(ta.pos, tb.pos);
    for (x, y) in ta.total_energy().iter().zip(tb.total_energy()) {
        assert!((x - y).abs() < 1e-12);
    }
    assert!(!b.forces.velocity_dependent());

    // The whole-system form, with an error path.
    let mut w = World::with_integrator("rk4").unwrap();
    w.add_particle(Vec3::ZERO, v(1.0, 0.0, 0.0), 1.0).unwrap();
    w.add_force(ClosureForce::new("Brake", |_t, _p, vel, _m, acc| {
        for (a, v) in acc.iter_mut().zip(vel) {
            *a -= *v;
        }
        Ok(())
    }));
    w.run(1e-3, 1000, 1000).unwrap();
    assert!((w.state.vel[0].x - (-1.0f64).exp()).abs() < 1e-12);
    assert!(w.forces.velocity_dependent());
}

#[test]
fn invalid_parameters_are_rejected() {
    assert!(SpringNetwork::new(vec![0], vec![0], vec![1.0], vec![1.0]).is_err());
    assert!(SpringNetwork::new(vec![0, 1], vec![1], vec![1.0], vec![1.0]).is_err());
    assert!(SpringNetwork::new(vec![0], vec![1], vec![1.0], vec![-1.0]).is_err());
    let mut y = Yukawa {
        center: Vec3::ZERO,
        k: 1.0,
        length: 1.0,
    };
    assert!(y.set_param("length", Param::Scalar(0.0)).is_err());
    let mut j = J2Oblateness {
        central: 0,
        g: 1.0,
        j2: 1e-3,
        radius: 1.0,
        axis: v(0.0, 0.0, 1.0),
    };
    assert!(j.set_param("axis", Param::Vector(Vec3::ZERO)).is_err());
    assert_eq!(j.axis, v(0.0, 0.0, 1.0));
    let mut d = DampedSpring {
        i: 0,
        j: 1,
        k: 1.0,
        rest_length: 1.0,
        c: 0.1,
    };
    assert!(d.set_param("c", Param::Scalar(-0.1)).is_err());
    let mut p = PostNewtonian {
        central: 0,
        g: 1.0,
        c: 1.0,
    };
    assert!(p.set_param("c", Param::Scalar(0.0)).is_err());

    // Particle bookkeeping: forces that name particles block their removal and renumber.
    let mut w = World::with_integrator("verlet").unwrap();
    for k in 0..4 {
        w.add_particle(v(k as f64, 0.0, 0.0), Vec3::ZERO, 1.0)
            .unwrap();
    }
    let net = w
        .add_force(SpringNetwork::new(vec![1, 2], vec![2, 3], vec![1.0; 2], vec![1.0; 2]).unwrap());
    w.remove_particle(0).unwrap();
    assert!(w.remove_particle(1).is_err());
    match w.forces.get(net).unwrap().builtin() {
        Some(BuiltinForce::SpringNetwork(n)) => {
            assert_eq!(n.i(), &[0, 1]);
            assert_eq!(n.j(), &[1, 2]);
        }
        _ => panic!("expected a SpringNetwork"),
    }
}
