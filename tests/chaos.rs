use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// Hénon-Heiles particle at x = 0 with energy `e` (px from the energy).
fn henon_heiles(y: f64, py: f64, e: f64, integrator: &str) -> World {
    let potential = 0.5 * y * y - y * y * y / 3.0;
    let px = (2.0 * (e - potential) - py * py).sqrt();
    let mut w = World::with_integrator(integrator).unwrap();
    w.add_particle(v(0.0, y, 0.0), v(px, py, 0.0), 1.0).unwrap();
    w.add_force(HenonHeiles {
        center: Vec3::ZERO,
        lambda: 1.0,
    });
    w
}

fn options(n: usize, record_every: usize) -> LyapunovOptions {
    LyapunovOptions {
        n_exponents: n,
        record_every,
        ..Default::default()
    }
}

/// The same force without its analytic Jacobian, so the finite-difference path is used.
struct Opaque(Box<dyn Force>);

impl Force for Opaque {
    fn accumulate(&self, t: f64, p: &[Vec3], u: &[Vec3], m: &[f64], a: &mut [Vec3]) -> Result<()> {
        self.0.accumulate(t, p, u, m, a)
    }
    fn velocity_dependent(&self) -> bool {
        self.0.velocity_dependent()
    }
    fn name(&self) -> String {
        format!("Opaque({})", self.0.name())
    }
}

#[test]
fn analytic_jacobians_match_finite_differences() {
    let c = v(0.1, -0.2, 0.05);
    let make: Vec<Box<dyn Fn() -> Box<dyn Force>>> = vec![
        Box::new(|| {
            Box::new(NewtonianGravity {
                g: 1.3,
                softening: 0.1,
            })
        }),
        Box::new(|| {
            Box::new(Spring {
                i: 0,
                j: 2,
                k: 3.0,
                rest_length: 0.5,
            })
        }),
        Box::new(move || {
            Box::new(AnchorSpring {
                i: 1,
                anchor: c,
                k: 2.0,
                rest_length: 0.7,
            })
        }),
        Box::new(move || {
            Box::new(AnchorSpring {
                i: 3,
                anchor: c,
                k: 2.0,
                rest_length: 0.0,
            })
        }),
        Box::new(|| Box::new(LinearDrag { gamma: 0.3 })),
        Box::new(|| {
            Box::new(UniformField {
                g: v(0.0, -9.8, 0.0),
            })
        }),
        Box::new(move || {
            Box::new(PowerLaw {
                center: c,
                k: -1.2,
                n: -1.0,
            })
        }),
        Box::new(move || {
            Box::new(PowerLaw {
                center: c,
                k: 0.6,
                n: 0.0,
            })
        }),
        Box::new(move || {
            Box::new(PowerLaw {
                center: c,
                k: 0.3,
                n: 2.7,
            })
        }),
        Box::new(move || {
            Box::new(Yukawa {
                center: c,
                k: 1.1,
                length: 0.7,
            })
        }),
        Box::new(move || {
            Box::new(PlummerPotential {
                center: c,
                gm: 2.0,
                a: 0.4,
            })
        }),
        Box::new(move || {
            Box::new(HernquistPotential {
                center: c,
                gm: 2.0,
                a: 0.4,
            })
        }),
        Box::new(move || {
            Box::new(HarmonicTrap {
                center: c,
                omega: v(1.0, 2.0, 0.5),
            })
        }),
        Box::new(move || {
            Box::new(HenonHeiles {
                center: c,
                lambda: 0.8,
            })
        }),
        Box::new(|| {
            Box::new(PeriodicForce {
                i: 2,
                amplitude: v(1.0, 0.0, 0.0),
                omega: 2.0,
                phase: 0.1,
            })
        }),
    ];
    let pos = [
        v(1.0, 0.2, -0.3),
        v(-0.7, 1.1, 0.4),
        v(0.3, -0.9, 1.2),
        v(1.6, 1.4, 0.1),
    ];
    let vel = [
        v(0.1, 0.0, 0.2),
        v(-0.3, 0.5, 0.0),
        v(0.0, 0.0, -0.4),
        v(0.2, 0.1, 0.3),
    ];
    let mass = [1.0, 2.0, 0.5, 1.5];
    let pinned = [false; 4];
    let dpos = [
        v(0.3, -0.1, 0.2),
        v(0.0, 0.4, -0.2),
        v(-0.5, 0.1, 0.0),
        v(0.1, 0.1, 0.1),
    ];
    let dvel = [
        v(-0.2, 0.3, 0.0),
        v(0.1, 0.0, 0.5),
        v(0.0, -0.3, 0.2),
        v(0.4, 0.0, -0.1),
    ];
    for f in &make {
        let mut analytic = ForceSet::new();
        analytic.add(f());
        let mut numeric = ForceSet::new();
        numeric.add(Box::new(Opaque(f())));
        let (mut a, mut b) = (Vec::new(), Vec::new());
        analytic
            .jacobian_vector(
                0.7, &pos, &vel, &mass, &[0.0; 4], &[0.0; 4], &pinned, &dpos, &dvel, &mut a,
            )
            .unwrap();
        numeric
            .jacobian_vector(
                0.7, &pos, &vel, &mass, &[0.0; 4], &[0.0; 4], &pinned, &dpos, &dvel, &mut b,
            )
            .unwrap();
        let scale = b.iter().map(|x| x.norm()).fold(1e-3, f64::max);
        for (x, y) in a.iter().zip(&b) {
            assert!(
                (*x - *y).norm() < 1e-7 * scale,
                "{}: analytic {x:?} vs finite differences {y:?}",
                f().name()
            );
        }
    }
}

#[test]
fn henon_heiles_chaos_and_regularity() {
    // Near the escape energy E = 1/6 most orbits are chaotic; at E = 1/12 they are regular.
    let mut chaotic = henon_heiles(0.0, 0.0, 1.0 / 6.0, "yoshida4");
    let c = chaotic
        .lyapunov(0.01, 200_000, &options(1, 10_000))
        .unwrap();
    let mut regular = henon_heiles(0.0, 0.0, 1.0 / 12.0, "yoshida4");
    let r = regular
        .lyapunov(0.01, 200_000, &options(1, 10_000))
        .unwrap();
    assert!(c.exponents[0] > 0.05, "chaotic λ = {}", c.exponents[0]);
    assert!(r.exponents[0] < 0.01, "regular λ = {}", r.exponents[0]);
    // MEGNO: ⟨Y⟩ → 2 for quasi-periodic motion, grows linearly for chaos.
    let (cy, ry) = (c.mean_megno.last().unwrap(), r.mean_megno.last().unwrap());
    assert!((ry - 2.0).abs() < 0.1, "regular ⟨Y⟩ = {ry}");
    assert!(*cy > 20.0, "chaotic ⟨Y⟩ = {cy}");
    // The regular estimate decays like ln t / t; the chaotic one settles.
    assert!(r.running[19][0] < 0.5 * r.running[4][0]);
    assert_eq!(c.t.len(), 20);
}

#[test]
fn kepler_motion_is_regular() {
    let mut w = World::with_integrator("yoshida4").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    w.pin(0, true).unwrap();
    w.add_particle(v(1.5, 0.0, 0.0), v(0.0, (0.5f64 / 1.5).sqrt(), 0.0), 1e-6)
        .unwrap();
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    let run = w.lyapunov(0.01, 100_000, &options(1, 10_000)).unwrap();
    assert!(run.exponents[0] < 0.01, "λ = {}", run.exponents[0]);
    assert!((run.mean_megno.last().unwrap() - 2.0).abs() < 0.05);
}

#[test]
fn full_spectrum_is_symplectic() {
    // Hamiltonian flow: exponents come in ± pairs and sum to zero.
    let mut w = henon_heiles(-0.1, 0.0, 1.0 / 6.0, "yoshida4");
    let run = w.lyapunov(0.01, 100_000, &options(6, 10_000)).unwrap();
    let l = &run.exponents;
    assert_eq!(l.len(), 6);
    assert!(
        l.iter().sum::<f64>().abs() < 1e-12,
        "sum {}",
        l.iter().sum::<f64>()
    );
    assert!((l[0] + l[5]).abs() < 0.02 * l[0], "{l:?}");
    assert!(l[1..5].iter().all(|x| x.abs() < 0.01), "{l:?}");
}

#[test]
fn integrators_agree_on_the_classification_and_the_trajectory_is_unchanged() {
    // Chaotic trajectories from different integrators diverge, so finite-time exponents differ;
    // what must agree is the classification: clearly positive in the chaotic sea, ~0 on a torus.
    for name in ["verlet", "yoshida4", "rk4", "gauss4", "pefrl", "dopri5"] {
        let mut w = henon_heiles(0.0, 0.0, 1.0 / 6.0, name);
        let chaotic = w.lyapunov(0.005, 200_000, &options(1, 200_000)).unwrap();
        let mut w = henon_heiles(0.0, 0.0, 1.0 / 12.0, name);
        let regular = w.lyapunov(0.005, 200_000, &options(1, 200_000)).unwrap();
        assert!(
            chaotic.exponents[0] > 0.02,
            "{name}: chaotic {}",
            chaotic.exponents[0]
        );
        assert!(
            regular.exponents[0] < 0.01,
            "{name}: regular {}",
            regular.exponents[0]
        );
    }

    // The real particles follow exactly the trajectory of a plain run.
    let mut a = henon_heiles(0.1, 0.1, 1.0 / 8.0, "yoshida4");
    let mut b = henon_heiles(0.1, 0.1, 1.0 / 8.0, "yoshida4");
    a.lyapunov(0.01, 1000, &options(2, 100)).unwrap();
    b.run(0.01, 1000, 1000).unwrap();
    assert_eq!(a.state.pos, b.state.pos);
    assert_eq!(a.state.vel, b.state.vel);
    assert_eq!(a.state.t, b.state.t);
    assert_eq!(a.forces.len(), 1); // forces restored
}

#[test]
fn rejects_bad_input_and_restores_the_world() {
    let mut w = henon_heiles(0.0, 0.0, 0.1, "yoshida4");
    assert!(w.lyapunov(0.01, 10, &options(7, 1)).is_err()); // more than 6 N
    assert!(w.lyapunov(0.01, 10, &options(0, 1)).is_err());
    assert!(w.lyapunov(0.0, 10, &options(1, 1)).is_err());
    w.set_integrator(integrators::by_name("wisdom_holman").unwrap());
    assert!(w.lyapunov(0.01, 10, &options(1, 1)).is_err());
    assert_eq!(w.forces.len(), 1);

    // A failing force: the world stops at the last good step with its forces back.
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(v(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_force(ClosureForce::new("Breaks", |t, _p, _v, _m, _a| {
        if t > 0.05 {
            Err(SimError::Invalid("broken".into()))
        } else {
            Ok(())
        }
    }));
    let err = w.lyapunov(0.01, 100, &options(1, 10)).unwrap_err();
    assert!(err.to_string().contains("broken"));
    assert!(w.state.t > 0.0 && w.state.t <= 0.06, "{}", w.state.t);
    assert_eq!(w.forces.len(), 1);
}
