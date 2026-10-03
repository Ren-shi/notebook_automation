use physim::integrators::NAMES;
use physim::*;

fn p(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

#[test]
fn removing_a_particle_renumbers_index_based_forces() {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(p(-1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_particle(p(0.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_particle(p(1.5, 0.0, 0.0), Vec3::ZERO, 2.0).unwrap();
    w.add_force(Spring {
        i: 1,
        j: 2,
        k: 3.0,
        rest_length: 1.0,
    });
    w.add_force(AnchorSpring {
        i: 2,
        anchor: p(3.0, 0.0, 0.0),
        k: 1.0,
        rest_length: 0.0,
    });

    let err = w.remove_particle(1).unwrap_err().to_string();
    assert!(err.contains("Spring(1, 2)"), "{err}");
    w.remove_particle(0).unwrap();

    let mut expected = World::with_integrator("verlet").unwrap();
    expected
        .add_particle(p(0.0, 0.0, 0.0), Vec3::ZERO, 1.0)
        .unwrap();
    expected
        .add_particle(p(1.5, 0.0, 0.0), Vec3::ZERO, 2.0)
        .unwrap();
    expected.add_force(Spring {
        i: 0,
        j: 1,
        k: 3.0,
        rest_length: 1.0,
    });
    expected.add_force(AnchorSpring {
        i: 1,
        anchor: p(3.0, 0.0, 0.0),
        k: 1.0,
        rest_length: 0.0,
    });
    assert_eq!(
        w.accelerations().unwrap(),
        expected.accelerations().unwrap()
    );
    assert!(w.remove_particle(5).is_err());
}

#[test]
fn force_ids_support_remove_replace_and_parameter_changes() {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(p(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    let grav = w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    let field = w.add_force(UniformField {
        g: p(0.0, -1.0, 0.0),
    });
    let a1 = w.accelerations().unwrap()[0];

    w.set_force_params(grav, &[("G".into(), Param::Scalar(2.0))])
        .unwrap();
    let a2 = w.accelerations().unwrap()[0];
    assert_eq!(a2.x, 2.0 * a1.x);

    // All-or-nothing: the invalid softening rolls back the G change.
    let bad = [
        ("G".into(), Param::Scalar(5.0)),
        ("softening".into(), Param::Scalar(-1.0)),
    ];
    assert!(w.set_force_params(grav, &bad).is_err());
    assert_eq!(
        w.forces.get(grav).unwrap().params()[0],
        ("G", Param::Scalar(2.0))
    );
    assert!(w
        .set_force_params(grav, &[("nope".into(), Param::Scalar(1.0))])
        .is_err());
    assert!(w
        .set_force_params(field, &[("g".into(), Param::Scalar(1.0))])
        .is_err());

    w.forces
        .replace(field, Box::new(UniformField { g: Vec3::ZERO }))
        .unwrap();
    w.remove_force(grav).unwrap();
    assert_eq!(w.accelerations().unwrap()[0], Vec3::ZERO);
    assert!(w.remove_force(grav).is_err());
    let ids: Vec<_> = w.forces.iter().map(|(id, _)| id).collect();
    assert_eq!(ids, vec![field]);
}

#[test]
fn massless_tracers_feel_but_do_not_source_gravity() {
    let mut w = World::with_integrator("yoshida4").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 1.0).unwrap();
    for k in 0..5 {
        let r = 1.0 + k as f64;
        w.add_particle(p(r, 0.0, 0.0), p(0.0, (1.0 / r).sqrt(), 0.0), 0.0)
            .unwrap();
    }
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.0,
    });
    w.run(1e-3, 10_000, 10_000).unwrap();
    assert_eq!(w.state.pos[0], Vec3::ZERO);
    assert_eq!(w.state.vel[0], Vec3::ZERO);
    for k in 0..5 {
        let r = w.state.pos[k + 1].norm();
        assert!((r - (1.0 + k as f64)).abs() < 1e-9, "tracer {k} radius {r}");
    }
}

#[test]
fn forces_defined_by_force_reject_massless_particles() {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(Vec3::ZERO, Vec3::ZERO, 0.0).unwrap();
    w.add_force(AnchorSpring {
        i: 0,
        anchor: p(1.0, 0.0, 0.0),
        k: 1.0,
        rest_length: 0.0,
    });
    assert!(w.step(0.1).is_err());
}

#[test]
fn pinned_particles_stay_put_with_every_integrator() {
    for &name in NAMES {
        let mut w = World::with_integrator(name).unwrap();
        w.add_particle(p(0.0, 2.0, 0.0), Vec3::ZERO, 1.0).unwrap();
        w.add_particle(p(1.0, 0.0, 0.0), p(0.0, 0.3, 0.0), 1.0)
            .unwrap();
        w.add_force(Spring {
            i: 0,
            j: 1,
            k: 5.0,
            rest_length: 1.0,
        });
        w.add_force(UniformField {
            g: p(0.0, -1.0, 0.0),
        });
        w.pin(0, true).unwrap();
        assert!(w
            .set_velocities(vec![p(1.0, 0.0, 0.0), Vec3::ZERO])
            .is_err());
        w.run(1e-3, 2000, 2000).unwrap();
        assert_eq!(w.state.pos[0], p(0.0, 2.0, 0.0), "{name}");
        assert!(
            (w.state.pos[1] - p(1.0, 0.0, 0.0)).norm() > 1e-3,
            "{name}: free particle should move"
        );

        w.pin(0, false).unwrap();
        w.run(1e-3, 100, 100).unwrap();
        assert_ne!(
            w.state.pos[0],
            p(0.0, 2.0, 0.0),
            "{name}: released particle should move"
        );
    }
}

/// Fails (error or NaN) once `t` passes `t_fail`.
struct FailsAfter {
    t_fail: f64,
    nan: bool,
}

impl Force for FailsAfter {
    fn accumulate(
        &self,
        t: f64,
        _: &[Vec3],
        _: &[Vec3],
        _: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        if t > self.t_fail {
            if self.nan {
                acc[0] = p(f64::NAN, 0.0, 0.0);
            } else {
                return Err(SimError::Invalid("boom".into()));
            }
        }
        Ok(())
    }
    fn name(&self) -> String {
        "FailsAfter".into()
    }
}

#[test]
fn failed_run_keeps_last_good_state_and_partial_trajectory() {
    for nan in [false, true] {
        let mut w = World::with_integrator("verlet").unwrap();
        w.add_particle(Vec3::ZERO, p(1.0, 0.0, 0.0), 1.0).unwrap();
        w.add_force(FailsAfter { t_fail: 0.55, nan });
        let failure = w.run(0.1, 100, 2).unwrap_err();
        // The step from t = 0.5 evaluates forces at t = 0.6 and fails.
        assert!((w.state.t - 0.5).abs() < 1e-12, "t = {}", w.state.t);
        assert!((w.state.pos[0].x - 0.5).abs() < 1e-12);
        let traj = &failure.trajectory;
        assert_eq!(traj.n_frames(), 4); // t = 0, 0.2, 0.4, then the last good state 0.5
        assert_eq!(traj.t.last(), Some(&w.state.t));
        assert_eq!(traj.kinetic.len(), traj.potential.len());
    }
}
