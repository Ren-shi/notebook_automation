use physim::integrators::{Langevin, NoseHoover};
use physim::*;

fn v(x: f64, y: f64, z: f64) -> Vec3 {
    Vec3::new(x, y, z)
}

fn lj(cutoff: f64, period: Option<Vec3>) -> PairPotential {
    PairPotential::new(
        PairKind::LennardJones {
            epsilon: 1.0,
            sigma: 1.0,
        },
        cutoff,
        true,
        period,
    )
}

/// N = 4 k³ atoms on an fcc lattice at density `rho` with temperature `t`.
fn liquid(k: usize, rho: f64, t: f64) -> (World, f64) {
    let n = 4 * k * k * k;
    let side = (n as f64 / rho).cbrt();
    let a = side / k as f64;
    let basis = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ];
    let mut w = World::with_integrator("verlet").unwrap();
    let mut s = 99u64;
    let mut u = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64 - 0.5
    };
    for x in 0..k {
        for y in 0..k {
            for z in 0..k {
                for b in basis {
                    let p = v(
                        (x as f64 + b[0]) * a,
                        (y as f64 + b[1]) * a,
                        (z as f64 + b[2]) * a,
                    );
                    w.add_particle(p, v(u(), u(), u()), 1.0).unwrap();
                }
            }
        }
    }
    let p = w.state.momentum() / n as f64;
    let scale = (t / {
        let vel: Vec<Vec3> = w.state.vel.iter().map(|q| *q - p).collect();
        w.set_velocities(vel).unwrap();
        w.temperature()
    })
    .sqrt();
    let vel: Vec<Vec3> = w.state.vel.iter().map(|q| *q * scale).collect();
    w.set_velocities(vel).unwrap();
    w.add_force(lj(2.5, Some(v(side, side, side))));
    (w, side)
}

#[test]
fn lennard_jones_liquid_conserves_energy() {
    // 108 atoms at the triple point, melted with Langevin, then NVE with Verlet.
    let (mut w, _) = liquid(3, 0.8442, 0.722);
    w.set_integrator(Box::new(Langevin::new(0.722, 1.0, 1).unwrap()));
    w.run(0.005, 2000, 2000).unwrap();
    w.set_integrator(integrators::by_name("verlet").unwrap());
    let p_start = w.state.momentum();
    let traj = w.run(0.005, 20_000, 100).unwrap();
    let e = traj.total_energy();
    let worst = e
        .iter()
        .map(|x| ((x - e[0]) / e[0]).abs())
        .fold(0.0, f64::max);
    assert!(worst < 1e-3, "energy deviation {worst:e}");
    let dp = w.state.momentum() - p_start;
    assert!(dp.norm() < 1e-10, "momentum change {dp:?}");
}

#[test]
fn pair_forces_are_minus_the_gradient_and_respect_the_box() {
    let period = Some(v(5.0, 6.0, 7.0));
    let pots = vec![
        lj(2.4, period),
        PairPotential::new(
            PairKind::Morse {
                depth: 1.3,
                a: 1.7,
                r0: 1.1,
            },
            2.4,
            false,
            period,
        ),
    ];
    // Particles straddling the periodic boundary.
    let pos = vec![
        v(0.1, 0.2, 0.3),
        v(4.9, 0.1, 6.8),
        v(1.0, 5.8, 0.5),
        v(0.4, 1.2, 1.1),
        v(3.0, 3.0, 3.0),
    ];
    let mass = vec![1.0, 2.0, 1.5, 0.7, 1.0];
    let vel = vec![Vec3::ZERO; 5];
    for f in &pots {
        let mut acc = vec![Vec3::ZERO; 5];
        f.accumulate(0.0, &pos, &vel, &mass, &mut acc).unwrap();
        let u = |p: &[Vec3]| f.potential(0.0, p, &mass).unwrap().unwrap();
        let h = 1e-6;
        for i in 0..5 {
            for axis in 0..3 {
                let mut e = [0.0; 3];
                e[axis] = h;
                let (mut a, mut b) = (pos.clone(), pos.clone());
                a[i] += Vec3::from(e);
                b[i] -= Vec3::from(e);
                let grad = (u(&a) - u(&b)) / (2.0 * h);
                let f_axis = mass[i] * [acc[i].x, acc[i].y, acc[i].z][axis];
                assert!(
                    (f_axis + grad).abs() < 1e-6 * grad.abs().max(1.0),
                    "{}: {i} {axis}: {f_axis} vs {}",
                    f.name(),
                    -grad
                );
            }
        }
        // Translating everything by a box vector changes nothing.
        let shifted: Vec<Vec3> = pos.iter().map(|p| *p + v(5.0, -6.0, 14.0)).collect();
        let mut acc2 = vec![Vec3::ZERO; 5];
        f.accumulate(0.0, &shifted, &vel, &mass, &mut acc2).unwrap();
        for (a, b) in acc.iter().zip(&acc2) {
            assert!((*a - *b).norm() < 1e-9);
        }
    }

    // A table built from Lennard-Jones reproduces it.
    let (r_min, dr, rc): (f64, f64, f64) = (0.5, 0.001, 2.5);
    let m = ((rc - r_min) / dr).round() as usize + 1;
    let r: Vec<f64> = (0..m).map(|k| r_min + k as f64 * dr).collect();
    let vtab: Vec<f64> = r.iter().map(|x| 4.0 * (x.powi(-12) - x.powi(-6))).collect();
    let dvtab: Vec<f64> = r
        .iter()
        .map(|x| 4.0 * (-12.0 * x.powi(-13) + 6.0 * x.powi(-7)))
        .collect();
    let table = PairPotential::new(
        PairKind::Table {
            r_min,
            dr,
            v: vtab,
            dv: dvtab,
        },
        rc,
        true,
        period,
    );
    let (mut a, mut b) = (vec![Vec3::ZERO; 5], vec![Vec3::ZERO; 5]);
    table.accumulate(0.0, &pos, &vel, &mass, &mut a).unwrap();
    lj(rc, period)
        .accumulate(0.0, &pos, &vel, &mass, &mut b)
        .unwrap();
    for (x, y) in a.iter().zip(&b) {
        assert!(
            (*x - *y).norm() < 1e-6 * y.norm().max(1.0),
            "{x:?} vs {y:?}"
        );
    }
}

#[test]
fn virial_pressure_matches_the_volume_derivative() {
    let (w, side) = liquid(3, 0.8, 1.0);
    let volume = side.powi(3);
    let p = w.pressure(volume).unwrap();
    // Scale positions and box by (1 ± ε): −dU/dV from the configurational energy.
    let energy = |s: f64| {
        let f = lj(2.5, Some(v(side * s, side * s, side * s)));
        let pos: Vec<Vec3> = w.state.pos.iter().map(|q| *q * s).collect();
        f.potential(0.0, &pos, &w.state.mass).unwrap().unwrap()
    };
    let eps = 1e-6;
    let du_dv = (energy(1.0 + eps) - energy(1.0 - eps))
        / (volume * ((1.0 + eps).powi(3) - (1.0 - eps).powi(3)));
    let kinetic = 2.0 * w.kinetic_energy() / (3.0 * volume);
    // The energy shift is constant, so this is exact apart from pairs crossing the cutoff.
    assert!(
        (p - (kinetic - du_dv)).abs() < 1e-5 * p.abs().max(1.0),
        "{p} vs {}",
        kinetic - du_dv
    );
}

#[test]
fn neighbour_lists_never_change_results() {
    // A run that rebuilds its list many times gives exactly the forces of a fresh evaluation.
    let (mut w, side) = liquid(3, 0.8442, 1.5);
    w.run(0.005, 300, 300).unwrap();
    let a = w.accelerations().unwrap();
    let fresh = lj(2.5, Some(v(side, side, side)));
    let mut b = vec![Vec3::ZERO; w.state.len()];
    fresh
        .accumulate(0.0, &w.state.pos, &w.state.vel, &w.state.mass, &mut b)
        .unwrap();
    assert_eq!(a, b);
    // And restarts from checkpoints are bit-exact.
    let (mut x, _) = liquid(3, 0.8442, 1.5);
    let (mut y, _) = liquid(3, 0.8442, 1.5);
    x.run(0.005, 200, 200).unwrap();
    y.run(0.005, 80, 80).unwrap();
    let mut z = World::from_checkpoint(y.checkpoint(), |_, _| unreachable!()).unwrap();
    z.run(0.005, 120, 120).unwrap();
    assert_eq!(x.state.pos, z.state.pos);
}

#[test]
fn thermostats_reach_their_temperature() {
    for (name, integrator) in [
        (
            "langevin",
            Box::new(Langevin::new(1.2, 1.0, 5).unwrap()) as Box<dyn Integrator>,
        ),
        ("nose_hoover", Box::new(NoseHoover::new(1.2, 0.5).unwrap())),
    ] {
        let (mut w, _) = liquid(3, 0.8442, 0.3);
        w.set_integrator(integrator);
        w.run(0.005, 3000, 3000).unwrap();
        let e0 = w.total_energy().unwrap() + w.thermostat_energy();
        let traj = w.run(0.005, 20_000, 10).unwrap();
        let n = w.state.len() as f64;
        let mean: f64 = traj
            .kinetic
            .iter()
            .map(|k| 2.0 * k / (3.0 * n))
            .sum::<f64>()
            / traj.kinetic.len() as f64;
        assert!((mean / 1.2 - 1.0).abs() < 0.03, "{name}: <T> = {mean}");
        if name == "nose_hoover" {
            let e1 = w.total_energy().unwrap() + w.thermostat_energy();
            assert!(
                ((e1 - e0) / e0).abs() < 1e-3,
                "extended energy {e0} -> {e1}"
            );
        }
    }

    // Stochastic runs restart exactly from checkpoints.
    let (mut x, _) = liquid(3, 0.8, 1.0);
    let (mut y, _) = liquid(3, 0.8, 1.0);
    for w in [&mut x, &mut y] {
        w.set_integrator(Box::new(Langevin::new(1.0, 2.0, 9).unwrap()));
    }
    x.run(0.005, 100, 100).unwrap();
    y.run(0.005, 40, 40).unwrap();
    let mut z = World::from_checkpoint(y.checkpoint(), |_, _| unreachable!()).unwrap();
    assert_eq!(z.integrator().name(), "langevin");
    z.run(0.005, 60, 60).unwrap();
    assert_eq!(x.state.vel, z.state.vel);

    assert!(Langevin::new(-1.0, 1.0, 0).is_err());
    assert!(NoseHoover::new(1.0, 0.0).is_err());
    assert!(lj(3.0, Some(v(5.0, 5.0, 5.0)))
        .potential(0.0, &[Vec3::ZERO; 2], &[1.0; 2])
        .is_err());
}
