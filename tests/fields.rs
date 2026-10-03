use physim::*;
use std::f64::consts::PI;

fn order(e_coarse: f64, e_fine: f64) -> f64 {
    (e_coarse / e_fine).log2()
}

fn max_err(a: &[f64], b: impl Fn(usize) -> f64) -> f64 {
    a.iter()
        .enumerate()
        .map(|(i, v)| (v - b(i)).abs())
        .fold(0.0, f64::max)
}

#[test]
fn wave_equation_is_second_order_and_conserves_energy() {
    // 1D string fixed at both ends: u = sin(πx) cos(πct).
    let c = 1.3;
    let run = |n: usize| {
        let h = 1.0 / (n - 1) as f64;
        let grid = Grid::new(&[n], h, Boundary::Dirichlet).unwrap();
        let mut w = Wave::new(grid, c).unwrap();
        let u0: Vec<f64> = (0..n).map(|i| (PI * i as f64 * h).sin()).collect();
        w.set_state(u0.clone(), vec![0.0; n]).unwrap();
        let t_end = 1.0;
        let steps = (t_end / (0.5 * h / c)).round() as usize;
        let e0 = w.energy();
        w.step(t_end / steps as f64, steps).unwrap();
        assert!((w.energy() / e0 - 1.0).abs() < 1e-3);
        max_err(&w.u, |i| u0[i] * (PI * c * t_end).cos())
    };
    let (a, b) = (run(33), run(65));
    assert!((order(a, b) - 2.0).abs() < 0.1, "1D order {}", order(a, b));

    // 2D periodic membrane: u = sin(2πx) sin(2πy) cos(2√2 π c t).
    let run2 = |n: usize| {
        let h = 1.0 / n as f64;
        let grid = Grid::new(&[n, n], h, Boundary::Periodic).unwrap();
        let mut w = Wave::new(grid, c).unwrap();
        let u0: Vec<f64> = (0..n * n)
            .map(|idx| {
                let x = grid.coordinates(idx);
                (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).sin()
            })
            .collect();
        w.set_state(u0.clone(), vec![0.0; n * n]).unwrap();
        let t_end = 0.5;
        let steps = (t_end / (0.5 * w.max_stable_dt())).ceil() as usize;
        w.step(t_end / steps as f64, steps).unwrap();
        let omega = 2.0 * 2f64.sqrt() * PI * c;
        max_err(&w.u, |i| u0[i] * (omega * t_end).cos())
    };
    let (a, b) = (run2(16), run2(32));
    assert!((order(a, b) - 2.0).abs() < 0.15, "2D order {}", order(a, b));

    // Steps above the CFL limit are refused.
    let grid = Grid::new(&[10, 10], 0.1, Boundary::Periodic).unwrap();
    let mut w = Wave::new(grid, 1.0).unwrap();
    assert!(w.step(0.1, 1).is_err());
}

#[test]
fn heat_equation_explicit_and_crank_nicolson_converge() {
    let d = 0.7;
    // Explicit, 1D Dirichlet: u = sin(πx) e^{-π² D t}, dt ∝ h² so the error is O(h²).
    let explicit = |n: usize| {
        let h = 1.0 / (n - 1) as f64;
        let grid = Grid::new(&[n], h, Boundary::Dirichlet).unwrap();
        let mut s = Heat::new(grid, d, HeatMethod::Explicit).unwrap();
        let u0: Vec<f64> = (0..n).map(|i| (PI * i as f64 * h).sin()).collect();
        s.set_u(u0.clone()).unwrap();
        let t_end = 0.1;
        let steps = (t_end / (0.4 * h * h / d)).ceil() as usize;
        s.step(t_end / steps as f64, steps).unwrap();
        max_err(&s.u, |i| u0[i] * (-PI * PI * d * t_end).exp())
    };
    let (a, b) = (explicit(17), explicit(33));
    assert!(
        (order(a, b) - 2.0).abs() < 0.1,
        "explicit order {}",
        order(a, b)
    );

    // Crank-Nicolson, 2D Neumann (cell-centred). cos(πx) cos(πy) is an exact eigenvector of
    // the discrete Laplacian with eigenvalue μ = 2 (2 cos(πh) - 2) / h², so the semi-discrete
    // solution e^{μ D t} isolates the time error, and the continuous e^{-2π² D t} the space error.
    let cn = |n: usize, dt: f64, continuous: bool| {
        let h = 1.0 / n as f64;
        let grid = Grid::new(&[n, n], h, Boundary::Neumann).unwrap();
        let mut s = Heat::new(grid, d, HeatMethod::CrankNicolson).unwrap();
        let mode: Vec<f64> = (0..n * n)
            .map(|idx| {
                let x = grid.coordinates(idx);
                (PI * x[0]).cos() * (PI * x[1]).cos()
            })
            .collect();
        s.set_u(mode.iter().map(|m| m + 1.0).collect()).unwrap();
        let total = s.total();
        let t_end = 0.1;
        let steps = (t_end / dt).round() as usize;
        s.step(t_end / steps as f64, steps).unwrap();
        assert!((s.total() / total - 1.0).abs() < 1e-10, "heat is conserved");
        let rate = if continuous {
            -2.0 * PI * PI
        } else {
            2.0 * (2.0 * (PI * h).cos() - 2.0) / (h * h)
        };
        max_err(&s.u, |i| mode[i] * (rate * d * t_end).exp() + 1.0)
    };
    // Steps of 0.01 and 0.005 are 30-60x the explicit limit on this grid.
    let (a, b) = (cn(32, 0.01, false), cn(32, 0.005, false));
    assert!(
        (order(a, b) - 2.0).abs() < 0.1,
        "CN time order {}",
        order(a, b)
    );
    let (a, b) = (cn(16, 1e-4, true), cn(32, 1e-4, true));
    assert!(
        (order(a, b) - 2.0).abs() < 0.1,
        "CN space order {}",
        order(a, b)
    );

    // The explicit scheme refuses unstable steps.
    let grid = Grid::new(&[10], 0.1, Boundary::Periodic).unwrap();
    let mut s = Heat::new(grid, 1.0, HeatMethod::Explicit).unwrap();
    assert!(s.step(0.01, 1).is_err());
}

#[test]
fn poisson_solvers_are_second_order() {
    // Periodic 3D by FFT: φ = sin 2πx sin 2πy sin 2πz, f = -12π² φ.
    let periodic = |n: usize| {
        let h = 1.0 / n as f64;
        let grid = Grid::new(&[n, n, n], h, Boundary::Periodic).unwrap();
        let exact: Vec<f64> = (0..grid.len())
            .map(|idx| {
                let x = grid.coordinates(idx);
                (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).sin() * (2.0 * PI * x[2]).sin()
            })
            .collect();
        let f: Vec<f64> = exact.iter().map(|p| -12.0 * PI * PI * p).collect();
        let mut phi = vec![0.0; grid.len()];
        poisson(&grid, &f, &mut phi).unwrap();
        // The FFT solves the discrete equations exactly.
        let mut lap = vec![0.0; grid.len()];
        grid.laplacian(&phi, &mut lap);
        assert!(max_err(&lap, |i| f[i]) < 1e-9);
        max_err(&phi, |i| exact[i])
    };
    let (a, b) = (periodic(16), periodic(32));
    assert!(
        (order(a, b) - 2.0).abs() < 0.1,
        "periodic order {}",
        order(a, b)
    );

    // Dirichlet 2D by conjugate gradients, with non-zero boundary values:
    // φ = sin(πx) sinh(πy) / sinh(π) + x y, ∇²φ = 0.
    let dirichlet = |n: usize| {
        let h = 1.0 / (n - 1) as f64;
        let grid = Grid::new(&[n, n], h, Boundary::Dirichlet).unwrap();
        let exact: Vec<f64> = (0..grid.len())
            .map(|idx| {
                let x = grid.coordinates(idx);
                (PI * x[0]).sin() * (PI * x[1]).sinh() / PI.sinh() + x[0] * x[1]
            })
            .collect();
        let mut phi: Vec<f64> = (0..grid.len())
            .map(|i| if grid.is_fixed(i) { exact[i] } else { 0.0 })
            .collect();
        poisson(&grid, &vec![0.0; grid.len()], &mut phi).unwrap();
        max_err(&phi, |i| exact[i])
    };
    let (a, b) = (dirichlet(17), dirichlet(33));
    assert!(
        (order(a, b) - 2.0).abs() < 0.15,
        "Dirichlet order {}",
        order(a, b)
    );

    // Neumann 1D: φ = cos(πx), f = -π² cos(πx).
    let neumann = |n: usize| {
        let h = 1.0 / n as f64;
        let grid = Grid::new(&[n], h, Boundary::Neumann).unwrap();
        let exact: Vec<f64> = (0..n)
            .map(|i| (PI * grid.coordinates(i)[0]).cos())
            .collect();
        let f: Vec<f64> = exact.iter().map(|p| -PI * PI * p).collect();
        let mut phi = vec![0.0; n];
        poisson(&grid, &f, &mut phi).unwrap();
        let mean = exact.iter().sum::<f64>() / n as f64;
        max_err(&phi, |i| exact[i] - mean)
    };
    let (a, b) = (neumann(32), neumann(64));
    assert!(
        (order(a, b) - 2.0).abs() < 0.15,
        "Neumann order {}",
        order(a, b)
    );
}

fn plummer(n: usize, seed: u64) -> (Vec<Vec3>, Vec<f64>) {
    let mut s = seed;
    let mut u = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        ((s >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    let pos = (0..n)
        .map(|_| {
            let r = 1.0 / ((u() * 0.99).powf(-2.0 / 3.0) - 1.0).sqrt();
            let z = 2.0 * u() - 1.0;
            let phi = 2.0 * PI * u();
            let s = (1.0 - z * z).sqrt();
            Vec3::new(s * phi.cos(), s * phi.sin(), z) * r.min(4.0)
        })
        .collect();
    (pos, vec![1.0 / n as f64; n])
}

#[test]
fn isolated_particle_mesh_matches_direct_summation() {
    let (pos, mass) = plummer(4000, 7);
    let vel = vec![Vec3::ZERO; pos.len()];
    let pm = ParticleMesh::new(1.0, 64, 10.0, Vec3::ZERO, false, None);
    let h = pm.spacing();
    let direct = NewtonianGravity {
        g: 1.0,
        softening: h,
    };
    let mut a_pm = vec![Vec3::ZERO; pos.len()];
    let mut a_dir = vec![Vec3::ZERO; pos.len()];
    pm.accumulate(0.0, &pos, &vel, &mass, &mut a_pm).unwrap();
    direct
        .accumulate(0.0, &pos, &vel, &mass, &mut a_dir)
        .unwrap();
    let mut errs: Vec<f64> = a_pm
        .iter()
        .zip(&a_dir)
        .map(|(p, d)| (*p - *d).norm() / d.norm())
        .collect();
    errs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = errs[errs.len() / 2];
    let p90 = errs[errs.len() * 9 / 10];
    assert!(
        median < 0.02 && p90 < 0.05,
        "median {median:e}, p90 {p90:e}"
    );
    // Momentum is conserved (no self-force, antisymmetric kernel).
    let total = a_pm
        .iter()
        .zip(&mass)
        .fold(Vec3::ZERO, |s, (a, m)| s + *a * *m);
    assert!(total.norm() < 1e-12, "{total:?}");
    // The potential energy agrees too.
    let u_pm = pm.potential(0.0, &pos, &mass).unwrap().unwrap();
    let u_dir = direct.potential(0.0, &pos, &mass).unwrap().unwrap();
    assert!((u_pm / u_dir - 1.0).abs() < 0.01, "{u_pm} vs {u_dir}");

    // Two particles: cloud-in-cell smoothing weakens the force at short range; it converges to
    // the softened 1/r² law within a few cells (measured: -3.8% at 4 cells, -1.0% at 8).
    for (r, tol) in [(4.0, 0.05), (8.0, 0.015), (16.0, 0.005)] {
        let p = vec![
            Vec3::new(-0.5 * r * h, 0.1, 0.0),
            Vec3::new(0.5 * r * h, 0.1, 0.0),
        ];
        let mut a = vec![Vec3::ZERO; 2];
        pm.accumulate(0.0, &p, &[Vec3::ZERO; 2], &[1.0, 1.0], &mut a)
            .unwrap();
        let d = r * h;
        let expected = d / (d * d + h * h).powf(1.5);
        assert!(
            (a[0].x / expected - 1.0).abs() < tol,
            "r = {r} cells: {} vs {expected}",
            a[0].x
        );
    }

    // Particles outside the box are refused.
    let mut a = vec![Vec3::ZERO; 1];
    assert!(pm
        .accumulate(
            0.0,
            &[Vec3::new(6.0, 0.0, 0.0)],
            &[Vec3::ZERO],
            &[1.0],
            &mut a
        )
        .is_err());
}

#[test]
fn periodic_particle_mesh_matches_the_linear_plane_wave() {
    // A lattice displaced by A sin(kq) along x has density contrast -kA cos(kx) to first
    // order, whose field is g_x = 4πGρ̄ A sin(kx). The grid matches the lattice spacing (a
    // finer grid than the particle spacing aliases the lattice into the force).
    let (n, l) = (32usize, 1.0);
    let k = 2.0 * PI / l;
    let amp = 1e-4;
    let mut pos = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for m in 0..n {
                let q = Vec3::new(i as f64 + 0.5, j as f64 + 0.5, m as f64 + 0.5) * (l / n as f64);
                pos.push(q + Vec3::new(amp * (k * q.x).sin(), 0.0, 0.0));
            }
        }
    }
    let total_mass = 1.0;
    let mass = vec![total_mass / pos.len() as f64; pos.len()];
    let rho = total_mass / l.powi(3);
    let pm = ParticleMesh::new(1.0, n, l, Vec3::new(0.5, 0.5, 0.5) * l, true, None);
    let mut a = vec![Vec3::ZERO; pos.len()];
    pm.accumulate(0.0, &pos, &vec![Vec3::ZERO; pos.len()], &mass, &mut a)
        .unwrap();
    let peak = 4.0 * PI * rho * amp;
    let worst = pos
        .iter()
        .zip(&a)
        .map(|(x, a)| (a.x - peak * (k * x.x).sin()).abs() / peak)
        .fold(0.0, f64::max);
    assert!(worst < 0.03, "relative error {worst:e}");
    // Translating everything by a box length changes nothing.
    let shifted: Vec<Vec3> = pos
        .iter()
        .map(|x| *x + Vec3::new(l, -2.0 * l, 0.0))
        .collect();
    let mut b = vec![Vec3::ZERO; pos.len()];
    pm.accumulate(0.0, &shifted, &vec![Vec3::ZERO; pos.len()], &mass, &mut b)
        .unwrap();
    assert!(a
        .iter()
        .zip(&b)
        .all(|(x, y)| (*x - *y).norm() < 1e-9 * peak));
}
