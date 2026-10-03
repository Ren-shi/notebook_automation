use physim::integrators::NAMES;
use physim::*;

fn p(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

/// A force the checkpoint cannot describe, standing in for a Python `CustomForce`.
struct Radial {
    k: f64,
}

impl Force for Radial {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _v: &[Vec3],
        _m: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        for (a, r) in acc.iter_mut().zip(pos) {
            *a -= *r * self.k;
        }
        Ok(())
    }
    fn name(&self) -> String {
        "radial".into()
    }
}

/// Every kind of built-in force, a pinned particle, a removed force (so ids have a gap)
/// and a force that has to be supplied again on restore.
fn busy_world(integrator: &str, n: usize) -> World {
    let mut w = World::with_integrator(integrator).unwrap();
    for i in 0..n {
        let a = i as f64 * 2.399963;
        let r = 1.0 + 0.01 * i as f64;
        w.add_particle(
            p(r * a.cos(), r * a.sin(), 0.05 * (i % 7) as f64),
            p(-0.3 * a.sin(), 0.3 * a.cos(), 0.0),
            1.0 / n as f64,
        )
        .unwrap();
    }
    w.pin(0, true).unwrap();
    let doomed = w.add_force(UniformField {
        g: p(0.0, 0.0, -1.0),
    });
    w.add_force(NewtonianGravity {
        g: 1.0,
        softening: 0.05,
    });
    w.add_force(Spring {
        i: 1,
        j: 2,
        k: 3.0,
        rest_length: 0.5,
    });
    w.add_force(AnchorSpring {
        i: 4,
        anchor: p(0.0, 0.0, 1.0),
        k: 0.5,
        rest_length: 0.2,
    });
    w.add_force(LinearDrag { gamma: 0.01 });
    w.add_force(QuadraticDrag { c: 0.001 });
    w.add_force(Radial { k: 0.1 });
    w.remove_force(doomed).unwrap();
    w.add_force(UniformField {
        g: p(0.0, -0.1, 0.0),
    });
    w
}

fn supply(_id: ForceId, name: &str) -> Result<Box<dyn Force>> {
    assert_eq!(name, "radial");
    Ok(Box::new(Radial { k: 0.1 }))
}

fn assert_same_bits(a: &World, b: &World) {
    assert_eq!(a.state.t.to_bits(), b.state.t.to_bits());
    let bits = |v: &[Vec3]| -> Vec<u64> {
        v.iter()
            .flat_map(|x| x.to_array())
            .map(f64::to_bits)
            .collect()
    };
    assert_eq!(bits(&a.state.pos), bits(&b.state.pos));
    assert_eq!(bits(&a.state.vel), bits(&b.state.vel));
}

#[test]
fn restart_from_checkpoint_matches_uninterrupted_run_bit_for_bit() {
    for &name in NAMES {
        if name == "wisdom_holman" {
            continue; // needs pure gravity and no pins; see tests/integrators.rs
        }
        let mut straight = busy_world(name, 12);
        straight.step(1e-3).unwrap(); // warm the acceleration cache, as a real run would be
        let mut interrupted = busy_world(name, 12);
        interrupted.step(1e-3).unwrap();

        straight.run(1e-3, 600, 100).unwrap();
        interrupted.run(1e-3, 250, 100).unwrap();
        let mut resumed = World::from_checkpoint(interrupted.checkpoint(), supply).unwrap();
        resumed.run(1e-3, 350, 100).unwrap();

        assert_same_bits(&straight, &resumed);
        assert_eq!(resumed.integrator().name(), name);
    }
}

#[test]
fn restart_is_exact_on_the_parallel_gravity_path() {
    // Enough bodies for the force loops to split into parallel blocks.
    let mut straight = busy_world("yoshida4", 500);
    let mut interrupted = busy_world("yoshida4", 500);
    straight.step(1e-3).unwrap();
    straight.step(1e-3).unwrap();
    interrupted.step(1e-3).unwrap();
    let mut resumed = World::from_checkpoint(interrupted.checkpoint(), supply).unwrap();
    resumed.step(1e-3).unwrap();
    assert_same_bits(&straight, &resumed);
}

#[test]
fn checkpoint_keeps_force_ids_parameters_and_pins() {
    let mut w = busy_world("verlet", 6);
    w.set_force_params(1, &[("G".into(), Param::Scalar(2.5))])
        .unwrap();
    let ids: Vec<ForceId> = w.forces.iter().map(|(id, _)| id).collect();
    let mut r = World::from_checkpoint(w.checkpoint(), supply).unwrap();
    assert_eq!(r.forces.iter().map(|(id, _)| id).collect::<Vec<_>>(), ids);
    assert_eq!(ids[0], 1, "removed force 0 must not come back");
    assert_eq!(
        r.forces.get(1).unwrap().params()[0],
        ("G", Param::Scalar(2.5))
    );
    assert_eq!(r.state.pinned, w.state.pinned);
    assert_eq!(r.state.mass, w.state.mass);
    // New forces continue the id sequence instead of reusing an old id.
    assert_eq!(
        r.add_force(LinearDrag { gamma: 1.0 }),
        w.add_force(LinearDrag { gamma: 1.0 })
    );
}

#[test]
fn restoring_needs_external_forces_and_rejects_corrupt_data() {
    let w = busy_world("rk4", 6);
    let err = World::from_checkpoint(w.checkpoint(), |_, name| {
        Err(SimError::Invalid(format!("no force {name}")))
    })
    .err()
    .unwrap();
    assert!(err.to_string().contains("radial"), "{err}");

    let mut bad = w.checkpoint();
    bad.state.vel[0] = p(1.0, 0.0, 0.0); // particle 0 is pinned
    assert!(World::from_checkpoint(bad, supply).is_err());

    let mut bad = w.checkpoint();
    bad.state.mass.pop();
    assert!(World::from_checkpoint(bad, supply).is_err());

    let mut bad = w.checkpoint();
    let dup = bad.forces[0].clone();
    bad.forces.push(dup);
    assert!(World::from_checkpoint(bad, supply).is_err());

    let mut bad = w.checkpoint();
    bad.integrator = "leapfrogger".into();
    assert!(World::from_checkpoint(bad, supply).is_err());

    let mut bad = w.checkpoint();
    bad.state.mass[2] = -1.0;
    assert!(World::from_checkpoint(bad, supply).is_err());
}

#[derive(Default)]
struct Counter {
    frames: usize,
    events: usize,
    last_t: f64,
    saw_energy: bool,
}

impl Recorder for Counter {
    fn frame(&mut self, f: &Frame<'_>) -> Result<()> {
        self.frames += 1;
        self.last_t = f.t;
        self.saw_energy |= f.kinetic.is_some() || f.potential.is_some();
        Ok(())
    }
    fn event(&mut self, _hit: EventHit) -> Result<()> {
        self.events += 1;
        Ok(())
    }
}

fn oscillator() -> World {
    let mut w = World::with_integrator("verlet").unwrap();
    w.add_particle(p(1.0, 0.0, 0.0), Vec3::ZERO, 1.0).unwrap();
    w.add_force(AnchorSpring {
        i: 0,
        anchor: Vec3::ZERO,
        k: 1.0,
        rest_length: 0.0,
    });
    w
}

#[test]
fn streaming_run_hands_every_frame_and_event_to_the_recorder() {
    let events = [Event {
        function: Box::new(events::CoordinateCrossing {
            i: 0,
            axis: 0,
            value: 0.0,
        }),
        direction: Direction::Either,
        terminal: false,
    }];
    let options = RunOptions {
        record_every: 1,
        events: &events,
        energies: false,
    };
    let steps = 1_000_000;
    let mut counter = Counter::default();
    let mut w = oscillator();
    let terminated = w.run_into(1e-3, steps, &options, &mut counter).unwrap();
    assert_eq!(terminated, None);
    assert_eq!(counter.frames, steps + 1);
    assert_eq!(counter.last_t, w.state.t);
    assert!(!counter.saw_energy);
    // x = cos t crosses zero at (k + 1/2) pi.
    assert_eq!(
        counter.events,
        (1000.0 / std::f64::consts::PI + 0.5) as usize
    );

    // Same physics with energies on, and the in-memory trajectory agrees with the stream.
    let mut v = oscillator();
    let traj = v
        .run_with_options(
            1e-3,
            1000,
            &RunOptions {
                energies: true,
                ..options
            },
        )
        .unwrap();
    let mut u = oscillator();
    let lean = u.run_with_options(1e-3, 1000, &options).unwrap();
    assert_eq!(traj.pos, lean.pos);
    assert_eq!(traj.kinetic.len(), traj.n_frames());
    assert!(lean.kinetic.is_empty() && lean.potential.is_empty() && !lean.energies);
}

#[test]
fn a_failing_recorder_stops_the_run() {
    struct Full(usize);
    impl Recorder for Full {
        fn frame(&mut self, _f: &Frame<'_>) -> Result<()> {
            self.0 += 1;
            if self.0 > 10 {
                return Err(SimError::Invalid("disk full".into()));
            }
            Ok(())
        }
    }
    let mut w = oscillator();
    let err = w
        .run_into(1e-2, 100, &RunOptions::default(), &mut Full(0))
        .unwrap_err();
    assert!(err.to_string().contains("disk full"));
    assert!(
        (w.state.t - 0.1).abs() < 1e-12,
        "stopped after the 11th frame, t = {}",
        w.state.t
    );
}
