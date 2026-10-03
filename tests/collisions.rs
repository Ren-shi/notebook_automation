use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

fn hard(restitution: f64, walls: Vec<Wall>) -> Collisions {
    Collisions {
        restitution,
        walls,
        ..Default::default()
    }
}

fn kinetic(w: &World) -> f64 {
    w.state.kinetic_energy()
}

#[test]
fn elastic_and_inelastic_pair_collisions() {
    // Unequal masses head-on: the textbook 1D result, independent of the step size.
    let (m1, m2, u1, u2) = (1.0, 3.0, 2.0, -1.0);
    for dt in [0.37, 0.01] {
        let mut w = World::with_integrator("verlet").unwrap();
        w.add_particle(v(0.0, 0.0, 0.0), v(u1, 0.0, 0.0), m1)
            .unwrap();
        w.add_particle(v(3.0, 0.0, 0.0), v(u2, 0.0, 0.0), m2)
            .unwrap();
        w.set_radii(vec![0.5, 0.5]).unwrap();
        w.set_collisions(Some(hard(1.0, vec![]))).unwrap();
        let (p0, e0) = (w.state.momentum(), kinetic(&w));
        w.run(dt, (3.0 / dt) as usize, 1000).unwrap();
        let v1 = ((m1 - m2) * u1 + 2.0 * m2 * u2) / (m1 + m2);
        let v2 = ((m2 - m1) * u2 + 2.0 * m1 * u1) / (m1 + m2);
        assert!(
            (w.state.vel[0].x - v1).abs() < 1e-12,
            "dt {dt}: {}",
            w.state.vel[0].x
        );
        assert!((w.state.vel[1].x - v2).abs() < 1e-12);
        assert!((w.state.momentum() - p0).norm() < 1e-12);
        assert!((kinetic(&w) - e0).abs() < 1e-12);
        assert_eq!(w.collision_count(), 1);
        // Contact happened exactly at touching distance: at t = 2/3 the gap closed.
        let gap = (w.state.pos[1] - w.state.pos[0]).norm();
        assert!(gap >= 1.0 - 1e-12);
    }

    // Oblique, equal masses, elastic, one at rest: outgoing velocities are perpendicular.
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(v(-3.0, 0.4, 0.0), v(1.0, 0.0, 0.0), 1.0)
        .unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.set_radii(vec![0.5, 0.5]).unwrap();
    w.set_collisions(Some(hard(1.0, vec![]))).unwrap();
    w.run(0.05, 120, 120).unwrap();
    assert!(w.state.vel[0].dot(w.state.vel[1]).abs() < 1e-12);
    assert!((kinetic(&w) - 0.5).abs() < 1e-12);

    // Restitution e: the kinetic energy lost is (1 - e²) μ u_n² / 2.
    let e = 0.6;
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, v(1.5, 0.0, 0.0), 2.0).unwrap();
    w.add_particle(v(2.0, 0.0, 0.0), v(-0.5, 0.0, 0.0), 1.0)
        .unwrap();
    w.set_radii(vec![0.4, 0.6]).unwrap();
    w.set_collisions(Some(hard(e, vec![]))).unwrap();
    let e0 = kinetic(&w);
    w.run(0.1, 20, 20).unwrap();
    let mu = 2.0 * 1.0 / 3.0;
    let loss = (1.0 - e * e) * mu * 2.0 * 2.0 / 2.0;
    assert!((e0 - kinetic(&w) - loss).abs() < 1e-12);
    assert!((w.state.vel[1].x - w.state.vel[0].x - e * 2.0).abs() < 1e-12);
}

#[test]
fn bouncing_ball_loses_height_by_e_squared() {
    // A ball dropped on a floor under gravity rebounds to e² of its height each bounce.
    let (e, g, h0, r) = (0.8, 9.81, 2.0, 0.1);
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(v(0.0, h0 + r, 0.0), Vec3::ZERO, 1.0)
        .unwrap();
    w.set_radius(0, r).unwrap();
    w.add_force(UniformField { g: v(0.0, -g, 0.0) });
    w.set_collisions(Some(hard(
        e,
        vec![Wall {
            normal: v(0.0, 1.0, 0.0),
            offset: 0.0,
        }],
    )))
    .unwrap();
    let traj = w.run(1e-3, 4000, 1).unwrap();
    // Apex heights: local maxima of y.
    let y: Vec<f64> = traj.pos.iter().map(|p| p.y - r).collect();
    let apexes: Vec<f64> = (1..y.len() - 1)
        .filter(|&k| y[k] > y[k - 1] && y[k] >= y[k + 1] && y[k] > 1e-3)
        .map(|k| y[k])
        .collect();
    assert!(apexes.len() >= 3, "{apexes:?}");
    for pair in apexes.windows(2) {
        let ratio = pair[1] / pair[0];
        // Apex sampled at the step: error O(g dt²).
        assert!((ratio - e * e).abs() < 1e-5, "ratio {ratio}");
    }
    assert!(y.iter().all(|&h| h > -1e-9), "ball went through the floor");
}

#[test]
fn hard_sphere_gas_relaxes_to_maxwell_boltzmann() {
    // 500 spheres in a box, all with the same speed in random directions. Elastic
    // collisions conserve energy exactly and thermalise the speeds: for a 3D Maxwellian
    // ⟨v⁴⟩/⟨v²⟩² = 5/3 (it starts at 1).
    let n = 500;
    let side = 10.0;
    let mut s = 11u64;
    let mut rnd = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut w = World::with_integrator("verlet").unwrap();
    let per_side = 8;
    for k in 0..n {
        let (a, b, c) = (
            k % per_side,
            (k / per_side) % per_side,
            k / (per_side * per_side),
        );
        let p =
            v(a as f64 + 0.6, b as f64 + 0.6, c as f64 + 0.6) * (side / (per_side as f64 + 0.2));
        let z = 2.0 * rnd() - 1.0;
        let phi = std::f64::consts::TAU * rnd();
        let st = (1.0 - z * z).sqrt();
        w.add_particle(p, v(st * phi.cos(), st * phi.sin(), z), 1.0)
            .unwrap();
    }
    w.set_radii(vec![0.25; n]).unwrap();
    w.set_collisions(Some(hard(
        1.0,
        Wall::box_walls(Vec3::ZERO, v(side, side, side)),
    )))
    .unwrap();
    let e0 = kinetic(&w);
    let moment_ratio = |w: &World| {
        let v2: Vec<f64> = w.state.vel.iter().map(|u| u.norm_squared()).collect();
        let m2 = v2.iter().sum::<f64>() / n as f64;
        let m4 = v2.iter().map(|x| x * x).sum::<f64>() / n as f64;
        m4 / (m2 * m2)
    };
    assert!((moment_ratio(&w) - 1.0).abs() < 1e-12);
    // Mean free time ~ 1 / (n_density · π d² · √2 v) ≈ 1.1; run ~30 of them.
    let mut ratios = Vec::new();
    for _ in 0..40 {
        w.run(0.05, 20, 20).unwrap();
        ratios.push(moment_ratio(&w));
    }
    assert!(((kinetic(&w) - e0) / e0).abs() < 1e-10, "energy drift");
    let late: f64 = ratios[20..].iter().sum::<f64>() / 20.0;
    // Finite-N sampling noise in ⟨v⁴⟩/⟨v²⟩² is ~0.1 for 500 particles.
    assert!((late - 5.0 / 3.0).abs() < 0.15, "moment ratio {late}");
    assert!(w.collision_count() > 5 * n as u64);
    // Every particle is still inside the box.
    for p in &w.state.pos {
        for axis in 0..3 {
            assert!(
                p[axis] >= 0.25 - 1e-9 && p[axis] <= side - 0.25 + 1e-9,
                "{p:?}"
            );
        }
    }
}

#[test]
fn soft_contact_matches_contact_theory() {
    let collide = |law: ContactLaw, k: f64, damping: f64| {
        let mut w = World::with_integrator("verlet").unwrap();
        w.add_particle(v(-1.0, 0.0, 0.0), v(1.0, 0.0, 0.0), 1.0)
            .unwrap();
        w.add_particle(v(1.0, 0.0, 0.0), v(-1.0, 0.0, 0.0), 1.0)
            .unwrap();
        w.set_radii(vec![0.5, 0.5]).unwrap();
        w.add_force(SoftContact { k, damping, law });
        let dt = 1e-5;
        let traj = w.run(dt, 200_000, 1).unwrap();
        let in_contact = traj
            .pos
            .chunks(2)
            .filter(|p| (p[1] - p[0]).norm() < 1.0)
            .count();
        (in_contact as f64 * dt, w.state.vel[1].x, traj)
    };
    let mu = 0.5; // reduced mass
    let speed = 2.0; // approach speed
                     // Linear spring: contact lasts half a period, π √(μ/k); elastic.
    let k = 1e4;
    let (duration, v_out, traj) = collide(ContactLaw::Linear, k, 0.0);
    let expected = std::f64::consts::PI * (mu / k).sqrt();
    assert!(
        (duration - expected).abs() < 2e-5,
        "{duration} vs {expected}"
    );
    assert!((v_out - 1.0).abs() < 1e-6);
    let e = traj.total_energy();
    assert!(((e[e.len() - 1] - e[0]) / e[0]).abs() < 1e-6);
    // Hertz: duration 2.94326 δ_max / v with δ_max = (5 μ v² / (4 k))^(2/5).
    let k = 1e5;
    let (duration, _, _) = collide(ContactLaw::Hertz, k, 0.0);
    let delta_max = (5.0 * mu * speed * speed / (4.0 * k)).powf(0.4);
    let expected = 2.94326 * delta_max / speed;
    assert!(
        (duration / expected - 1.0).abs() < 1e-3,
        "{duration} vs {expected}"
    );
    // Linear spring with a dashpot: restitution exp(-π ζ / √(1 - ζ²)), ζ = c / (2 √(μ k)).
    let (k, c) = (1e4, 10.0);
    let (_, v_out, _) = collide(ContactLaw::Linear, k, c);
    let zeta = c / (2.0 * (mu * k).sqrt());
    let restitution = (-std::f64::consts::PI * zeta / (1.0 - zeta * zeta).sqrt()).exp();
    // The never-attractive clamp ends contact slightly early, so e is a little larger.
    assert!(
        (v_out - restitution).abs() < 0.02,
        "{v_out} vs {restitution}"
    );
}

#[test]
fn collisions_survive_checkpoints_and_reject_bad_settings() {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, v(1.0, 0.3, 0.0), 1.0).unwrap();
    w.add_particle(v(2.0, 0.2, 0.0), v(-1.0, 0.0, 0.0), 2.0)
        .unwrap();
    w.set_radii(vec![0.4, 0.3]).unwrap();
    w.set_collisions(Some(hard(
        0.9,
        Wall::box_walls(v(-3.0, -3.0, -3.0), v(3.0, 3.0, 3.0)),
    )))
    .unwrap();
    let mut straight = World::from_checkpoint(w.checkpoint(), |_, _| unreachable!()).unwrap();
    straight.run(0.01, 1000, 1000).unwrap();
    w.run(0.01, 400, 400).unwrap();
    let mut resumed = World::from_checkpoint(w.checkpoint(), |_, _| unreachable!()).unwrap();
    resumed.run(0.01, 600, 600).unwrap();
    assert_eq!(straight.state.pos, resumed.state.pos);
    assert_eq!(resumed.state.radius, vec![0.4, 0.3]);
    assert!(resumed.collision_count() > 0);

    assert!(w.set_collisions(Some(hard(1.5, vec![]))).is_err());
    assert!(w
        .set_collisions(Some(hard(
            1.0,
            vec![Wall {
                normal: v(0.0, 2.0, 0.0),
                offset: 0.0
            }]
        )))
        .is_err());
    assert!(w.set_radii(vec![-1.0, 0.0]).is_err());
    assert!(w.run_adaptive(1.0, &AdaptiveOptions::default()).is_err());
}
