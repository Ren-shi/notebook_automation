use std::f64::consts::{PI, TAU};

use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Unit charge and mass with speed 1 perpendicular to B = (0, 0, 1): ω = 1, radius 1.
fn gyrating(integrator: &str) -> World {
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(v(0.0, 1.0, 0.0), v(1.0, 0.0, 0.0), 1.0)
        .unwrap();
    w.set_charge(0, 1.0).unwrap();
    w.add_force(MagneticField {
        b: v(0.0, 0.0, 1.0),
    });
    w
}

#[test]
fn boris_gyration_matches_qb_over_m_and_keeps_energy_to_round_off() {
    let h = 0.01;
    let mut w = gyrating("boris");
    let steps = 1_000_000;
    let traj = w.run(h, steps, 1000).unwrap();
    // |v| (hence kinetic energy) is conserved to round-off for any step.
    let worst = traj
        .vel
        .iter()
        .map(|u| (u.norm() - 1.0).abs())
        .fold(0.0, f64::max);
    assert!(worst < 1e-13, "speed drift {worst:e}");
    // The velocity turns by exactly 2 atan(ω h / 2) per step (ω = qB/m = 1, clockwise for
    // a positive charge in +z), so the frequency error is (ω h)²/12.
    let angle = steps as f64 * 2.0 * (0.5 * h).atan();
    let expected = v(angle.cos(), -angle.sin(), 0.0);
    assert!(
        (w.state.vel[0] - expected).norm() < 1e-9,
        "{:?}",
        w.state.vel[0]
    );
    let omega_num = 2.0 * (0.5 * h).atan() / h;
    assert!((omega_num - 1.0 + h * h / 12.0).abs() < 1e-8);
    // Radius m v / (q B) = 1 about the guiding centre at the origin.
    let radius = traj.pos.iter().map(|p| p.norm()).fold(0.0, f64::max);
    assert!((radius - 1.0).abs() < 1e-10, "radius {radius}");

    // RK4 at the same step drifts in energy.
    let mut w = gyrating("rk4");
    w.run(h, 100_000, 100_000).unwrap();
    assert!((w.state.vel[0].norm() - 1.0).abs() > 1e-10);
}

#[test]
fn e_cross_b_drift() {
    // Choose the step so the numerical gyration closes in exactly 50 steps; then the mean
    // velocity over whole gyrations is the guiding-centre drift.
    let n = 50;
    let theta = TAU / n as f64;
    let h = 2.0 * (0.5 * theta).tan();
    let (e, b) = (v(0.0, 0.1, 0.0), v(0.0, 0.0, 2.0));
    let mut w = World::with_integrator("boris").unwrap();
    w.add_particle(Vec3::ZERO, v(0.3, 0.0, 0.0), 2.0).unwrap();
    w.set_charge(0, 4.0).unwrap(); // q B / m = 4: rotation 2 atan(4 h / 2) per step
    w.add_force(ElectricField { e });
    w.add_force(MagneticField { b });
    // With qB/m = 4, use a quarter of h so the per-step angle is theta.
    let h = h / 4.0;
    let periods = 40;
    let start = w.state.pos[0];
    w.run(h, n * periods, n * periods).unwrap();
    let mean_velocity = (w.state.pos[0] - start) / (h * (n * periods) as f64);
    let drift = e.cross(b) / b.norm_squared();
    assert!(
        (mean_velocity - drift).norm() < 1e-12,
        "{mean_velocity:?} vs {drift:?}"
    );
}

#[test]
fn rutherford_scattering_and_coulomb_balance() {
    // A light charge q passing a pinned charge Q: tan(θ/2) = k q Q / (m v∞² b).
    let (k, q, big_q, m): (f64, f64, f64, f64) = (1.0, 1.0, 2.0, 1.0);
    let (b, speed, x0): (f64, f64, f64) = (0.5, 1.5, -1e5);
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.pin(0, true).unwrap();
    // Start far away with the speed that has the right energy at infinity.
    let r0 = (x0 * x0 + b * b).sqrt();
    let v0 = (speed * speed - 2.0 * k * q * big_q / (m * r0)).sqrt();
    w.add_particle(v(x0, b, 0.0), v(v0, 0.0, 0.0), m).unwrap();
    w.set_charges(vec![big_q, q]).unwrap();
    w.add_force(Coulomb { k, softening: 0.0 });
    let t_end = 2.0 * (-x0) / speed;
    let e0 = w.total_energy().unwrap();
    w.run_adaptive(
        t_end,
        &AdaptiveOptions {
            rtol: 1e-12,
            atol: 1e-15,
            output: Output::Times(&[]),
            ..Default::default()
        },
    )
    .unwrap();
    let u = w.state.vel[1];
    let deflection = u.y.atan2(u.x);
    let expected = 2.0 * (k * q * big_q / (m * speed * speed * b)).atan();
    // Corrections from starting and ending at a finite distance are O(kqQ/(m v² r)) ~ 1e-5.
    assert!(
        (deflection - expected).abs() < 5e-5,
        "deflection {deflection} vs {expected}"
    );
    assert!(((w.total_energy().unwrap() - e0) / e0).abs() < 1e-10);

    // Gravity and Coulomb cancel exactly when G m² = k q².
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(v(-0.5, 0.0, 0.0), Vec3::ZERO, 2.0).unwrap();
    w.add_particle(v(0.5, 0.0, 0.0), Vec3::ZERO, 2.0).unwrap();
    w.set_charges(vec![3.0, 3.0]).unwrap();
    w.add_force(NewtonianGravity {
        g: 9.0 / 4.0,
        softening: 0.0,
    });
    w.add_force(Coulomb {
        k: 1.0,
        softening: 0.0,
    });
    let a = w.accelerations().unwrap();
    assert!(a[0].norm() < 1e-15 && a[1].norm() < 1e-15, "{a:?}");
    assert!(w.potential_energy().unwrap().abs() < 1e-15);
}

#[test]
fn magnetic_mirror_conserves_the_magnetic_moment() {
    // B = B0 (-x z / L², -y z / L², 1 + z² / L²) (divergence-free to first order near the
    // axis): a particle with enough perpendicular speed bounces between the mirrors while
    // μ = m v⊥² / (2 |B|) stays nearly constant (adiabatic invariant).
    let (b0, l) = (5.0, 4.0);
    let field = move |_t: f64, pos: &[Vec3], out: &mut [Vec3]| {
        for (o, p) in out.iter_mut().zip(pos) {
            *o += v(
                -p.x * p.z / (l * l),
                -p.y * p.z / (l * l),
                1.0 + p.z * p.z / (l * l),
            ) * b0;
        }
        Ok(())
    };
    let bfield = |p: Vec3| {
        v(
            -p.x * p.z / (l * l),
            -p.y * p.z / (l * l),
            1.0 + p.z * p.z / (l * l),
        ) * b0
    };
    let mut w = World::with_integrator("boris").unwrap();
    w.add_particle(v(0.0, 0.2, 0.0), v(1.0, 0.0, 0.5), 1.0)
        .unwrap();
    w.set_charge(0, 1.0).unwrap();
    w.add_force(FieldFunctions::new("Mirror").magnetic(field));
    let traj = w.run(1e-3, 200_000, 50).unwrap();
    let mu: Vec<f64> = traj
        .pos
        .iter()
        .zip(&traj.vel)
        .map(|(p, u)| {
            let b = bfield(*p);
            let bhat = b / b.norm();
            let vperp2 = u.norm_squared() - u.dot(bhat).powi(2);
            0.5 * vperp2 / b.norm()
        })
        .collect();
    let (lo, hi) = mu
        .iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
    let z_max = traj.pos.iter().map(|p| p.z).fold(f64::MIN, f64::max);
    let z_min = traj.pos.iter().map(|p| p.z).fold(f64::MAX, f64::min);
    assert!(
        z_max > 1.0 && z_min < -1.0,
        "should bounce: z in [{z_min}, {z_max}]"
    );
    assert!(z_max < 3.0, "should be reflected: z_max {z_max}");
    assert!((hi - lo) / lo < 0.05, "μ varies by {:.3}", (hi - lo) / lo);
    // Magnetic forces do no work.
    let speed = |u: &Vec3| u.norm();
    let s0 = speed(&traj.vel[0]);
    assert!(traj.vel.iter().all(|u| (speed(u) - s0).abs() < 1e-12));
}

/// The same force without its analytic Jacobian.
struct Opaque(Box<dyn Force>);

impl Force for Opaque {
    fn accumulate(&self, t: f64, p: &[Vec3], u: &[Vec3], m: &[f64], a: &mut [Vec3]) -> Result<()> {
        self.0.accumulate(t, p, u, m, a)
    }
    fn accumulate_charged(
        &self,
        t: f64,
        p: &[Vec3],
        u: &[Vec3],
        m: &[f64],
        q: &[f64],
        a: &mut [Vec3],
    ) -> Result<()> {
        self.0.accumulate_charged(t, p, u, m, q, a)
    }
    fn velocity_dependent(&self) -> bool {
        self.0.velocity_dependent()
    }
    fn name(&self) -> String {
        format!("Opaque({})", self.0.name())
    }
}

#[test]
fn em_forces_have_correct_jacobians_and_survive_checkpoints() {
    let pos = [v(1.0, 0.2, -0.3), v(-0.7, 1.1, 0.4), v(0.3, -0.9, 1.2)];
    let vel = [v(0.1, 0.0, 0.2), v(-0.3, 0.5, 0.0), v(0.0, 0.0, -0.4)];
    let mass = [1.0, 2.0, 0.5];
    let charge = [1.0, -2.0, 0.5];
    let dpos = [v(0.3, -0.1, 0.2), v(0.0, 0.4, -0.2), v(-0.5, 0.1, 0.0)];
    let dvel = [v(-0.2, 0.3, 0.0), v(0.1, 0.0, 0.5), v(0.0, -0.3, 0.2)];
    let make: Vec<Box<dyn Fn() -> Box<dyn Force>>> = vec![
        Box::new(|| {
            Box::new(Coulomb {
                k: 1.3,
                softening: 0.1,
            })
        }),
        Box::new(|| {
            Box::new(MagneticField {
                b: v(0.2, -0.4, 1.0),
            })
        }),
        Box::new(|| {
            Box::new(ElectricField {
                e: v(0.2, -0.4, 1.0),
            })
        }),
    ];
    for f in &make {
        let (mut a, mut b) = (Vec::new(), Vec::new());
        let mut analytic = ForceSet::new();
        analytic.add(f());
        let mut numeric = ForceSet::new();
        numeric.add(Box::new(Opaque(f())));
        let pinned = [false; 3];
        analytic
            .jacobian_vector(
                0.0, &pos, &vel, &mass, &charge, &[0.0; 3], &pinned, &dpos, &dvel, &mut a,
            )
            .unwrap();
        numeric
            .jacobian_vector(
                0.0, &pos, &vel, &mass, &charge, &[0.0; 3], &pinned, &dpos, &dvel, &mut b,
            )
            .unwrap();
        for (x, y) in a.iter().zip(&b) {
            assert!((*x - *y).norm() < 1e-7, "{}: {x:?} vs {y:?}", f().name());
        }
    }

    // Charges and EM forces are saved; restarts are bit-exact.
    let mut w = gyrating("boris");
    w.add_particle(v(2.0, 0.0, 0.0), v(0.0, 0.5, 0.0), 3.0)
        .unwrap();
    w.set_charge(1, -0.5).unwrap();
    w.add_force(ElectricField {
        e: v(0.0, 0.05, 0.0),
    });
    w.add_force(Coulomb {
        k: 0.1,
        softening: 0.01,
    });
    let mut straight = World::from_checkpoint(w.checkpoint(), |_, _| unreachable!()).unwrap();
    assert_eq!(straight.state.charge, vec![1.0, -0.5]);
    straight.run(0.01, 500, 500).unwrap();
    w.run(0.01, 200, 200).unwrap();
    let mut resumed = World::from_checkpoint(w.checkpoint(), |_, _| unreachable!()).unwrap();
    resumed.run(0.01, 300, 300).unwrap();
    assert_eq!(straight.state.pos, resumed.state.pos);
}

#[test]
fn errors() {
    // Boris rejects velocity-dependent forces other than magnetic fields.
    let mut w = gyrating("boris");
    w.add_force(LinearDrag { gamma: 0.1 });
    assert!(w.step(0.01).unwrap_err().to_string().contains("LinearDrag"));
    // A charged massless particle cannot be accelerated by a field.
    let mut w = gyrating("verlet");
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 0.0).unwrap();
    w.set_charge(1, 1.0).unwrap();
    assert!(w.step(0.01).unwrap_err().to_string().contains("massless"));
    // Charges must be finite and one per particle.
    assert!(w.set_charges(vec![1.0]).is_err());
    assert!(w.set_charge(0, f64::NAN).is_err());
    // Uncharged particles ignore fields entirely; Boris is then leapfrog.
    let mut w = World::with_integrator("boris").unwrap();
    w.add_particle(v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), 1.0)
        .unwrap();
    w.add_force(MagneticField {
        b: v(0.0, 0.0, 1.0),
    });
    w.run(0.1, 10, 10).unwrap();
    assert!((w.state.pos[0] - v(1.0, 1.0, 0.0)).norm() < 1e-14);
    let _ = PI;
}
