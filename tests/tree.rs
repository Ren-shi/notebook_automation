use physim::*;

fn plummer(n: usize, seed: u64) -> (Vec<Vec3>, Vec<f64>) {
    let mut s = seed.wrapping_mul(0x9E3779B97F4A7C15) | 1;
    let mut rnd = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut pos = Vec::with_capacity(n);
    while pos.len() < n {
        let x: f64 = rnd();
        if x < 1e-9 {
            continue;
        }
        let r = 1.0 / (x.powf(-2.0 / 3.0) - 1.0).sqrt();
        if r > 30.0 {
            continue;
        }
        let z = 2.0 * rnd() - 1.0;
        let phi = std::f64::consts::TAU * rnd();
        let st = (1.0 - z * z).sqrt();
        pos.push(Vec3::new(r * st * phi.cos(), r * st * phi.sin(), r * z));
    }
    (pos, vec![1.0 / n as f64; n])
}

fn accel(f: &dyn Force, pos: &[Vec3], mass: &[f64]) -> Vec<Vec3> {
    let vel = vec![Vec3::ZERO; pos.len()];
    let mut acc = vec![Vec3::ZERO; pos.len()];
    f.accumulate(0.0, pos, &vel, mass, &mut acc).unwrap();
    acc
}

fn rms_error(a: &[Vec3], b: &[Vec3]) -> f64 {
    (a.iter()
        .zip(b)
        .map(|(x, y)| ((*x - *y).norm() / y.norm()).powi(2))
        .sum::<f64>()
        / a.len() as f64)
        .sqrt()
}

fn tree(theta: f64, quadrupole: bool, softening: f64) -> TreeGravity {
    TreeGravity {
        g: 1.0,
        softening,
        theta,
        quadrupole,
    }
}

#[test]
fn force_error_scales_with_the_opening_angle() {
    let (pos, mass) = plummer(4000, 1);
    let direct = accel(
        &NewtonianGravity {
            g: 1.0,
            softening: 0.0,
        },
        &pos,
        &mass,
    );
    let err = |theta, quad| rms_error(&accel(&tree(theta, quad, 0.0), &pos, &mass), &direct);
    let (m3, m6) = (err(0.3, false), err(0.6, false));
    let (q3, q6) = (err(0.3, true), err(0.6, true));
    assert!(
        err(0.5, false) < 2e-3,
        "monopole θ = 0.5: {:e}",
        err(0.5, false)
    );
    // Monopole error grows about as θ^3, quadrupole about as θ^4 (one power more).
    let (pm, pq) = ((m6 / m3).log2(), (q6 / q3).log2());
    assert!(pm > 2.0 && pm < 3.8, "monopole exponent {pm}");
    assert!(pq > pm + 0.5, "quadrupole exponent {pq} vs monopole {pm}");
    assert!(q3 < m3 / 3.0);
}

#[test]
fn matches_direct_summation_details() {
    // Softening, massless test particles (they feel but do not source gravity), G scaling.
    let (mut pos, mut mass) = plummer(500, 2);
    pos.push(Vec3::new(0.3, -0.2, 0.1));
    mass.push(0.0);
    let direct = accel(
        &NewtonianGravity {
            g: 2.5,
            softening: 0.05,
        },
        &pos,
        &mass,
    );
    let mut t = tree(0.4, true, 0.05);
    t.g = 2.5;
    let a = accel(&t, &pos, &mass);
    assert!(rms_error(&a, &direct) < 1e-3);
    let mut light = mass.clone();
    light[500] = 0.0;
    assert!(a[500].norm() > 0.0);

    // Potential energy agrees with direct summation to the force accuracy.
    let u_direct = NewtonianGravity {
        g: 2.5,
        softening: 0.05,
    }
    .potential(0.0, &pos, &mass)
    .unwrap()
    .unwrap();
    let u_tree = t.potential(0.0, &pos, &mass).unwrap().unwrap();
    assert!(
        ((u_tree - u_direct) / u_direct).abs() < 1e-3,
        "{u_tree} vs {u_direct}"
    );

    // Coincident particles and tiny systems are handled.
    let pos = vec![Vec3::new(1.0, 1.0, 1.0); 40];
    let mass = vec![1.0; 40];
    let a = accel(&tree(0.5, false, 0.1), &pos, &mass);
    assert!(a.iter().all(|x| x.norm() == 0.0));
    let a = accel(&tree(0.5, false, 0.0), &[Vec3::ZERO], &[1.0]);
    assert_eq!(a[0], Vec3::ZERO);
}

#[test]
fn runs_a_cluster_with_good_energy_conservation_and_exact_restarts() {
    let (pos, mass) = plummer(1500, 3);
    let build = |force: Box<dyn Force>| {
        let mut w = World::with_integrator("verlet").unwrap();
        for (p, m) in pos.iter().zip(&mass) {
            // Roughly virialised: circular-ish speeds around the centre.
            let r = p.norm().max(1e-3);
            let speed = (r * r / (r * r + 1.0).powf(1.5)).sqrt();
            let dir = Vec3::new(-p.y, p.x, 0.0);
            let dir = if dir.norm() > 0.0 {
                dir / dir.norm()
            } else {
                dir
            };
            w.add_particle(*p, dir * speed, *m).unwrap();
        }
        w.forces.add(force);
        w
    };
    let mut t = build(Box::new(tree(0.5, false, 0.05)));
    let p0 = t.state.momentum();
    let traj = t.run(1e-3, 300, 100).unwrap();
    let e = traj.total_energy();
    let drift = ((e[e.len() - 1] - e[0]) / e[0]).abs();
    assert!(drift < 1e-3, "energy drift {drift:e}");
    // Momentum is conserved only to the force accuracy.
    let dp = t.state.momentum() - p0;
    assert!(dp.norm() < 1e-4, "momentum change {dp:?}");

    let mut a = build(Box::new(tree(0.5, true, 0.05)));
    let mut b = build(Box::new(tree(0.5, true, 0.05)));
    a.run(1e-3, 50, 50).unwrap();
    b.run(1e-3, 20, 20).unwrap();
    let mut r = World::from_checkpoint(b.checkpoint(), |_, _| unreachable!()).unwrap();
    r.run(1e-3, 30, 30).unwrap();
    assert_eq!(a.state.pos, r.state.pos);
}

#[test]
fn rejects_bad_theta() {
    let mut t = tree(0.5, false, 0.0);
    assert!(t.set_param("theta", Param::Scalar(1.5)).is_err());
    assert!(t.set_param("theta", Param::Scalar(0.0)).is_err());
    t.theta = 2.0;
    let pos = vec![Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)];
    let vel = vec![Vec3::ZERO; 2];
    let mut acc = vec![Vec3::ZERO; 2];
    assert!(t
        .accumulate(0.0, &pos, &vel, &[1.0, 1.0], &mut acc)
        .is_err());
}
