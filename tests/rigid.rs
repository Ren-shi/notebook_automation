use physim::rigid::inertia;
use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Complete elliptic integral of the first kind K(m), m = k², by the AGM.
fn elliptic_k(m: f64) -> f64 {
    let (mut a, mut b) = (1.0, (1.0 - m).sqrt());
    while (a - b).abs() > 1e-15 * a {
        (a, b) = (0.5 * (a + b), (a * b).sqrt());
    }
    std::f64::consts::PI / (2.0 * a)
}

/// Period of the body-frame angular velocity of a torque-free top (Landau & Lifshitz §37).
fn euler_period(i: [f64; 3], omega: [f64; 3]) -> f64 {
    let l = [i[0] * omega[0], i[1] * omega[1], i[2] * omega[2]];
    let m2: f64 = l.iter().map(|x| x * x).sum();
    let e2: f64 = (0..3).map(|k| l[k] * l[k] / i[k]).sum();
    let [i1, i2, i3] = i;
    let (k2, scale) = if m2 > e2 * i2 {
        (
            (i2 - i1) * (e2 * i3 - m2) / ((i3 - i2) * (m2 - e2 * i1)),
            (i1 * i2 * i3 / ((i3 - i2) * (m2 - e2 * i1))).sqrt(),
        )
    } else {
        (
            (i3 - i2) * (m2 - e2 * i1) / ((i2 - i1) * (e2 * i3 - m2)),
            (i1 * i2 * i3 / ((i2 - i1) * (e2 * i3 - m2))).sqrt(),
        )
    };
    4.0 * elliptic_k(k2) * scale
}

fn free_top(omega: Vec3, order: u8) -> RigidSystem {
    let mut s = RigidSystem::new();
    s.set_order(order).unwrap();
    let mut b = RigidBody::new(1.0, v(1.0, 2.0, 3.0));
    b.orientation = Quat::from_axis_angle(v(1.0, 1.0, 0.3), 0.7);
    b.set_angular_velocity(b.to_space(omega));
    s.add_body(b).unwrap();
    s
}

#[test]
fn torque_free_top_flips_about_its_intermediate_axis() {
    let omega = [0.01, 1.0, 0.02];
    let mut s = free_top(Vec3::from(omega), 2);
    let l0 = s.bodies()[0].ang_mom;
    let e0 = s.total_energy().unwrap();
    let dt = 1e-3;
    let traj = s.run(dt, 200_000, 1).unwrap();
    // Body-frame angular velocity about the intermediate axis changes sign periodically.
    let w2: Vec<f64> = (0..traj.n_frames())
        .map(|f| {
            let q = traj.orientation[f];
            q.rotate_inverse(traj.angular_velocity[f]).y
        })
        .collect();
    let mut ups = Vec::new();
    for f in 1..w2.len() {
        if w2[f - 1] < 0.0 && w2[f] >= 0.0 {
            ups.push(traj.t[f - 1] + dt * w2[f - 1] / (w2[f - 1] - w2[f]));
        }
    }
    assert!(ups.len() >= 3, "only {} flips", ups.len());
    let expected = euler_period([1.0, 2.0, 3.0], omega);
    for p in ups.windows(2) {
        let period = p[1] - p[0];
        assert!(
            ((period - expected) / expected).abs() < 1e-4,
            "flip period {period} vs {expected}"
        );
    }
    // |L| (indeed L itself) is conserved exactly; energy to the splitting error.
    for (f, e) in traj.total_energy().iter().enumerate() {
        assert!((traj.ang_mom[f] - l0).norm() < 1e-14);
        assert!(((e - e0) / e0).abs() < 1e-6, "energy {e} vs {e0}");
    }
}

#[test]
fn rotation_step_converges_at_its_order() {
    let omega = v(0.3, 1.0, 0.4);
    let reference = {
        let mut s = free_top(omega, 4);
        s.run(1e-4, 50_000, 50_000).unwrap();
        s.bodies()[0].orientation
    };
    for (order, expect) in [(2u8, 2.0), (4, 4.0)] {
        let err = |dt: f64| {
            let mut s = free_top(omega, order);
            s.run(dt, (5.0 / dt).round() as usize, 1_000_000).unwrap();
            let q = s.bodies()[0].orientation;
            // q and -q are the same rotation.
            let d = q.w * reference.w + q.x * reference.x + q.y * reference.y + q.z * reference.z;
            let sgn = d.signum();
            let diff = [
                q.w - sgn * reference.w,
                q.x - sgn * reference.x,
                q.y - sgn * reference.y,
                q.z - sgn * reference.z,
            ];
            diff.iter().map(|c| c * c).sum::<f64>().sqrt()
        };
        let (a, b) = (err(0.04), err(0.02));
        let rate = (a / b).log2();
        assert!(
            (rate - expect).abs() < 0.3,
            "order {order}: observed {rate} ({a:e} -> {b:e})"
        );
    }
}

#[test]
fn fast_heavy_top_precesses_at_mgl_over_spin() {
    let (m, g, l, i3, spin) = (1.0, 1.0, 1.0, 0.5, 50.0);
    let mut s = RigidSystem::new();
    let mut top = RigidBody::pivoted(m, v(0.5, 0.5, i3), Vec3::ZERO, v(0.0, 0.0, l)).unwrap();
    assert!((top.inertia.x - 1.5).abs() < 1e-15 && top.inertia.z == i3);
    let tilt = 0.5;
    top.orientation = Quat::from_axis_angle(v(1.0, 0.0, 0.0), tilt);
    top.set_angular_velocity(top.to_space(v(0.0, 0.0, spin)));
    s.add_body(top).unwrap();
    s.add_force(Box::new(BodyGravity { g: v(0.0, 0.0, -g) }))
        .unwrap();
    let e0 = s.total_energy().unwrap();
    let lz0 = s.bodies()[0].ang_mom.z;
    let traj = s.run(1e-3, 50_000, 10).unwrap();
    let axis: Vec<Vec3> = traj
        .orientation
        .iter()
        .map(|q| q.rotate(v(0.0, 0.0, 1.0)))
        .collect();
    let mut phi: Vec<f64> = axis.iter().map(|a| a.y.atan2(a.x)).collect();
    for k in 1..phi.len() {
        while phi[k] - phi[k - 1] > std::f64::consts::PI {
            phi[k] -= 2.0 * std::f64::consts::PI;
        }
        while phi[k] - phi[k - 1] < -std::f64::consts::PI {
            phi[k] += 2.0 * std::f64::consts::PI;
        }
    }
    let rate = (phi[phi.len() - 1] - phi[0]) / (traj.t[traj.t.len() - 1] - traj.t[0]);
    let fast_top = m * g * l / (i3 * spin);
    assert!(
        ((rate - fast_top) / fast_top).abs() < 0.01,
        "precession {rate} vs {fast_top}"
    );
    // The tilt only nutates slightly, and energy and L_z (gravity has no z torque) hold.
    for a in &axis {
        assert!((a.z.acos() - tilt).abs() < 0.01);
    }
    let e1 = s.total_energy().unwrap();
    assert!(((e1 - e0) / e0).abs() < 1e-6, "energy {e0} -> {e1}");
    assert!((s.bodies()[0].ang_mom.z - lz0).abs() < 1e-10 * lz0.abs());
}

#[test]
fn torques_are_minus_the_angular_gradient() {
    let mut s = RigidSystem::new();
    let mut a = RigidBody::new(2.0, inertia::solid_box(2.0, 1.0, 0.6, 0.3));
    a.orientation = Quat::from_axis_angle(v(0.2, 1.0, -0.4), 1.1);
    let mut b =
        RigidBody::pivoted(1.5, v(0.2, 0.3, 0.4), v(2.0, 0.5, 0.0), v(0.0, 0.7, 0.0)).unwrap();
    b.orientation = Quat::from_axis_angle(v(1.0, -0.3, 0.2), 0.8);
    s.add_body(a).unwrap();
    s.add_body(b).unwrap();
    s.add_force(Box::new(BodyGravity {
        g: v(0.1, -9.8, 0.3),
    }))
    .unwrap();
    s.add_force(Box::new(BodySpring {
        a: Attachment::Body {
            body: 0,
            point: v(0.5, 0.3, -0.1),
        },
        b: Attachment::Body {
            body: 1,
            point: v(0.1, 0.0, 0.4),
        },
        k: 7.0,
        rest: 0.4,
    }))
    .unwrap();
    s.add_force(Box::new(BodySpring {
        a: Attachment::Fixed(v(-1.0, 1.0, 1.0)),
        b: Attachment::Body {
            body: 0,
            point: v(-0.5, 0.0, 0.2),
        },
        k: 3.0,
        rest: 1.0,
    }))
    .unwrap();
    let w = s.wrenches().unwrap();
    let h = 1e-6;
    for (i, wi) in w.iter().enumerate() {
        for axis in 0..3 {
            let mut e = [0.0; 3];
            e[axis] = 1.0;
            let n = Vec3::from(e);
            let energy_at = |rot: f64, shift: f64| {
                let mut t = RigidSystem::new();
                for (j, body) in s.bodies().iter().enumerate() {
                    let mut body = body.clone();
                    if j == i {
                        body.orientation = Quat::from_axis_angle(n, rot) * body.orientation;
                        body.pos += n * shift;
                    }
                    t.add_body(body).unwrap();
                }
                let mut u = 0.0;
                for f in s.forces() {
                    u += f.potential(0.0, t.bodies()).unwrap().unwrap();
                }
                u
            };
            let torque = -(energy_at(h, 0.0) - energy_at(-h, 0.0)) / (2.0 * h);
            assert!((wi.torque[axis] - torque).abs() < 1e-6, "{i} {axis}");
            if !s.bodies()[i].pivot {
                let force = -(energy_at(0.0, h) - energy_at(0.0, -h)) / (2.0 * h);
                assert!((wi.force[axis] - force).abs() < 1e-6, "{i} {axis}");
            }
        }
    }
}

#[test]
fn coupled_free_bodies_conserve_momentum_and_angular_momentum() {
    let mut s = RigidSystem::new();
    s.set_order(4).unwrap();
    let mut a = RigidBody::new(1.0, inertia::solid_ellipsoid(1.0, 0.5, 0.3, 0.2));
    a.set_angular_velocity(v(0.5, 2.0, -1.0));
    a.vel = v(0.1, 0.0, 0.0);
    let mut b = RigidBody::new(3.0, inertia::solid_cylinder(3.0, 0.2, 1.0));
    b.pos = v(1.5, 0.2, 0.0);
    b.orientation = Quat::from_axis_angle(v(0.0, 1.0, 0.0), 0.4);
    b.set_angular_velocity(v(-1.0, 0.0, 3.0));
    s.add_body(a).unwrap();
    s.add_body(b).unwrap();
    s.add_force(Box::new(BodySpring {
        a: Attachment::Body {
            body: 0,
            point: v(0.4, 0.1, 0.0),
        },
        b: Attachment::Body {
            body: 1,
            point: v(0.0, 0.1, -0.4),
        },
        k: 20.0,
        rest: 0.6,
    }))
    .unwrap();
    let (p0, l0, e0) = (
        s.momentum(),
        s.angular_momentum(),
        s.total_energy().unwrap(),
    );
    let traj = s.run(2e-3, 10_000, 100).unwrap();
    assert!((s.momentum() - p0).norm() < 1e-12);
    assert!((s.angular_momentum() - l0).norm() < 1e-11);
    for e in traj.total_energy() {
        assert!(((e - e0) / e0).abs() < 1e-5, "energy {e} vs {e0}");
    }

    // A free body under gravity falls on a parabola while spinning freely.
    let mut s = RigidSystem::new();
    let mut c = RigidBody::new(2.0, inertia::solid_sphere(2.0, 0.1));
    c.vel = v(1.0, 0.0, 3.0);
    c.set_angular_velocity(v(0.0, 4.0, 0.0));
    s.add_body(c).unwrap();
    s.add_force(Box::new(BodyGravity {
        g: v(0.0, 0.0, -9.8),
    }))
    .unwrap();
    s.run(0.01, 50, 50).unwrap();
    let c = &s.bodies()[0];
    assert!((c.pos - v(0.5, 0.0, 3.0 * 0.5 - 4.9 * 0.25)).norm() < 1e-12);
    assert!((c.angular_velocity() - v(0.0, 4.0, 0.0)).norm() < 1e-12);
}

#[test]
fn bad_bodies_are_rejected() {
    let mut s = RigidSystem::new();
    assert!(s.add_body(RigidBody::new(1.0, v(1.0, 1.0, 3.0))).is_err());
    assert!(s.add_body(RigidBody::new(0.0, v(1.0, 1.0, 1.0))).is_err());
    let mut b = RigidBody::new(1.0, v(1.0, 1.0, 1.0));
    b.com = v(0.0, 0.0, 1.0);
    assert!(s.add_body(b).is_err());
    assert!(RigidBody::pivoted(1.0, v(1.0, 1.0, 1.0), Vec3::ZERO, v(1.0, 1.0, 0.0)).is_err());
    assert!(s.set_order(3).is_err());
    assert!(s
        .add_force(Box::new(BodySpring {
            a: Attachment::Body {
                body: 0,
                point: Vec3::ZERO
            },
            b: Attachment::Fixed(Vec3::ZERO),
            k: 1.0,
            rest: 0.0,
        }))
        .is_err());
    // Quaternions compose and rotate like rotation matrices.
    let q = Quat::from_axis_angle(v(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
    assert!((q.rotate(v(1.0, 0.0, 0.0)) - v(0.0, 1.0, 0.0)).norm() < 1e-15);
    let r = Quat::from_axis_angle(v(1.0, 0.0, 0.0), 0.3);
    let x = v(0.2, -0.5, 0.9);
    assert!(((q * r).rotate(x) - q.rotate(r.rotate(x))).norm() < 1e-15);
    let m = q.to_matrix();
    assert!((m[1][0] - 1.0).abs() < 1e-15);
}
