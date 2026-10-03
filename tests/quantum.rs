use physim::*;
use std::f64::consts::PI;

fn order(e_coarse: f64, e_fine: f64) -> f64 {
    (e_coarse / e_fine).log2()
}

/// A 1D solver on `n` points spanning `[-length/2, length/2)` (periodic) or
/// `[-length/2, length/2]` (Dirichlet, edges included).
fn line(n: usize, length: f64, method: QuantumMethod) -> Schrodinger {
    let (boundary, h) = match method {
        QuantumMethod::SplitStep => (Boundary::Periodic, length / n as f64),
        QuantumMethod::CrankNicolson => (Boundary::Dirichlet, length / (n - 1) as f64),
    };
    let grid = Grid::new(&[n], h, boundary).unwrap();
    let mut s = Schrodinger::new(grid, 1.0, 1.0, method).unwrap();
    s.origin = [-0.5 * length, 0.0, 0.0];
    s
}

fn harmonic(s: &Schrodinger, omega: f64) -> Vec<f64> {
    (0..s.grid.len())
        .map(|i| {
            let x = s.coordinates(i);
            0.5 * s.mass * omega * omega * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
        })
        .collect()
}

fn max_diff(a: &[(f64, f64)], b: &[(f64, f64)]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(p, q)| ((p.0 - q.0).powi(2) + (p.1 - q.1).powi(2)).sqrt())
        .fold(0.0, f64::max)
}

#[test]
fn norm_and_energy_are_conserved() {
    for method in [QuantumMethod::SplitStep, QuantumMethod::CrankNicolson] {
        let mut s = line(256, 20.0, method);
        let v = harmonic(&s, 1.0);
        s.set_potential(v).unwrap();
        s.set_gaussian([2.0, 0.0, 0.0], [0.5, 1.0, 1.0], [1.0, 0.0, 0.0])
            .unwrap();
        let e0 = s.energy().unwrap();
        s.step(1e-3, 100_000).unwrap();
        assert!(
            (s.norm() - 1.0).abs() < 1e-11,
            "{method:?} norm {}",
            s.norm()
        );
        let drift = (s.energy().unwrap() / e0 - 1.0).abs();
        // Crank-Nicolson conserves its discrete energy exactly; split-step to O(dt²).
        let limit = if method == QuantumMethod::SplitStep {
            1e-6
        } else {
            1e-10
        };
        assert!(drift < limit, "{method:?} energy drift {drift}");
        assert!((s.t - 100.0).abs() < 1e-9);
    }

    // 2D Crank-Nicolson (BiCGSTAB) is unitary up to the solver tolerance.
    let grid = Grid::new(&[48, 48], 0.25, Boundary::Dirichlet).unwrap();
    let mut s = Schrodinger::new(grid, 1.0, 1.0, QuantumMethod::CrankNicolson).unwrap();
    s.origin = [-5.875, -5.875, 0.0];
    let v = harmonic(&s, 1.0);
    s.set_potential(v).unwrap();
    s.set_gaussian([1.0, -1.0, 0.0], [0.7, 0.7, 1.0], [0.5, 1.0, 0.0])
        .unwrap();
    let e0 = s.energy().unwrap();
    s.step(0.02, 200).unwrap();
    assert!(s.last_iterations > 0);
    assert!((s.norm() - 1.0).abs() < 1e-9, "2D norm {}", s.norm());
    assert!((s.energy().unwrap() / e0 - 1.0).abs() < 1e-9);
}

#[test]
fn free_packet_spreads_at_the_analytic_rate() {
    let (sigma0, p0, t_end) = (1.0, 2.0, 5.0);
    let width = |t: f64| sigma0 * (1.0 + (t / (2.0 * sigma0 * sigma0)).powi(2)).sqrt();
    // Errors in centre, width, ⟨p⟩ and ⟨T⟩ against the exact free evolution.
    let errors = |method: QuantumMethod, n: usize| {
        let mut s = line(n, 80.0, method);
        s.set_gaussian([-10.0, 0.0, 0.0], [sigma0, 1.0, 1.0], [p0, 0.0, 0.0])
            .unwrap();
        s.step(0.005, (t_end / 0.005) as usize).unwrap();
        let x = s.position().unwrap()[0];
        let sd = s.position_variance().unwrap()[0].sqrt();
        let kinetic = (p0 * p0 + 0.25 / (sigma0 * sigma0)) / 2.0;
        [
            (x - (-10.0 + p0 * t_end)).abs(),
            (sd / width(t_end) - 1.0).abs(),
            (s.momentum().unwrap()[0] - p0).abs(),
            (s.kinetic_energy().unwrap() - kinetic).abs(),
        ]
    };
    // Split-step: the free evolution is exact in Fourier space.
    for e in errors(QuantumMethod::SplitStep, 1024) {
        assert!(e < 1e-9, "split-step error {e}");
    }
    // Crank-Nicolson: dispersion error of the finite differences, second order in h.
    let (coarse, fine) = (
        errors(QuantumMethod::CrankNicolson, 2001),
        errors(QuantumMethod::CrankNicolson, 4001),
    );
    for (a, b) in coarse.iter().zip(&fine) {
        assert!(
            order(*a, *b) > 1.8,
            "Crank-Nicolson order {}",
            order(*a, *b)
        );
        assert!(*b < 4e-3, "{b}");
    }
}

#[test]
fn coherent_state_follows_the_classical_orbit() {
    let omega = 1.0;
    let mut s = line(256, 24.0, QuantumMethod::SplitStep);
    let v = harmonic(&s, omega);
    s.set_potential(v).unwrap();
    s.set_order(4).unwrap();
    let sigma = (1.0 / (2.0 * omega)).sqrt();
    s.set_gaussian([3.0, 0.0, 0.0], [sigma, 1.0, 1.0], [0.0; 3])
        .unwrap();
    for _ in 0..20 {
        s.step(0.01, 50).unwrap();
        let x = s.position().unwrap()[0];
        let p = s.momentum().unwrap()[0];
        // Ehrenfest is exact here; what remains is the O(dt⁴) splitting error.
        assert!(
            (x - 3.0 * (omega * s.t).cos()).abs() < 1e-7,
            "x {x} at {}",
            s.t
        );
        assert!((p + 3.0 * (omega * s.t).sin()).abs() < 1e-7);
        // The shape does not change.
        assert!((s.position_variance().unwrap()[0] - sigma * sigma).abs() < 1e-7);
    }
}

#[test]
fn time_stepping_converges_at_the_stated_orders() {
    // An anharmonic well, so nothing is exact; errors against a much finer step of the same
    // scheme isolate the time error.
    let run = |method: QuantumMethod, order: usize, dt: f64| {
        let mut s = line(128, 16.0, method);
        s.set_order(order).unwrap();
        let v: Vec<f64> = (0..128)
            .map(|i| {
                let x = s.coordinates(i)[0];
                0.25 * x.powi(4) - x * x
            })
            .collect();
        s.set_potential(v).unwrap();
        s.set_gaussian([1.0, 0.0, 0.0], [0.6, 1.0, 1.0], [0.5, 0.0, 0.0])
            .unwrap();
        let t_end = 1.0;
        s.step(dt, (t_end / dt).round() as usize).unwrap();
        s.psi
    };
    for (method, ord, dts) in [
        (QuantumMethod::SplitStep, 2, [0.02, 0.01]),
        (QuantumMethod::SplitStep, 4, [0.05, 0.025]),
        (QuantumMethod::CrankNicolson, 2, [0.005, 0.0025]),
    ] {
        let reference = run(method, ord, dts[1] / 32.0);
        let e1 = max_diff(&run(method, ord, dts[0]), &reference);
        let e2 = max_diff(&run(method, ord, dts[1]), &reference);
        let measured = order(e1, e2);
        assert!(
            (measured - ord as f64).abs() < 0.15,
            "{method:?} order {ord}: measured {measured} ({e1:.2e}, {e2:.2e})"
        );
    }
}

#[test]
fn eigenstates_match_the_analytic_spectra() {
    // Spectral (split-step) harmonic oscillator: exponentially accurate.
    let mut s = line(128, 20.0, QuantumMethod::SplitStep);
    let v = harmonic(&s, 1.0);
    s.set_potential(v).unwrap();
    let (e, states) = s.eigenstates(5, 1e-10, 500).unwrap();
    for (n, &en) in e.iter().enumerate() {
        assert!((en - (n as f64 + 0.5)).abs() < 1e-9, "E{n} = {en}");
    }
    let dv = s.cell_volume();
    for a in 0..5 {
        for b in 0..5 {
            let overlap: f64 = states[a].iter().zip(&states[b]).map(|(x, y)| x * y).sum();
            let expected = if a == b { 1.0 } else { 0.0 };
            assert!((overlap * dv - expected).abs() < 1e-9);
        }
    }
    // Ground state is the Gaussian π^{-1/4} e^{-x²/2}.
    let gs_err = (0..128)
        .map(|i| {
            let x = s.coordinates(i)[0];
            (states[0][i] - PI.powf(-0.25) * (-0.5 * x * x).exp()).abs()
        })
        .fold(0.0, f64::max);
    assert!(gs_err < 1e-8, "{gs_err}");
    // Using an eigenstate as ψ gives its energy, and it is stationary (up to the splitting
    // error, small at fourth order).
    s.set_order(4).unwrap();
    s.set_psi(states[2].iter().map(|&v| (v, 0.0)).collect())
        .unwrap();
    assert!((s.energy().unwrap() - 2.5).abs() < 1e-9);
    let before = s.psi.clone();
    s.step(0.01, 100).unwrap();
    let phase = (-2.5f64).sin_cos(); // e^{-iEt}, t = 1
    let rotated: Vec<(f64, f64)> = before
        .iter()
        .map(|&(re, _)| (re * phase.1, re * phase.0))
        .collect();
    assert!(max_diff(&s.psi, &rotated) < 1e-8);

    // Finite differences (Crank-Nicolson's operator): second order in the grid spacing.
    let fd_error = |n: usize| {
        let mut s = line(n, 20.0, QuantumMethod::CrankNicolson);
        let v = harmonic(&s, 1.0);
        s.set_potential(v).unwrap();
        let (e, _) = s.eigenstates(4, 1e-10, 500).unwrap();
        e.iter()
            .enumerate()
            .map(|(k, en)| (en - (k as f64 + 0.5)).abs())
            .fold(0.0, f64::max)
    };
    let (a, b) = (fd_error(101), fd_error(201));
    assert!((order(a, b) - 2.0).abs() < 0.05, "HO order {}", order(a, b));

    // Infinite square well of width L between the fixed edge nodes: n²π²ħ²/(2mL²).
    let well_error = |n: usize| {
        let mut s = line(n, 1.0, QuantumMethod::CrankNicolson);
        let (e, _) = s.eigenstates(3, 1e-10, 500).unwrap();
        e.iter()
            .enumerate()
            .map(|(k, en)| (en / ((k + 1) as f64 * PI).powi(2) * 2.0 - 1.0).abs())
            .fold(0.0, f64::max)
    };
    let (a, b) = (well_error(51), well_error(101));
    assert!(
        (order(a, b) - 2.0).abs() < 0.05,
        "well order {}",
        order(a, b)
    );
    assert!(b < 1e-3);

    // 2D isotropic oscillator: E = 1, 2, 2, 3, 3, 3 (degenerate levels).
    let grid = Grid::new(&[32, 32], 0.5, Boundary::Periodic).unwrap();
    let mut s = Schrodinger::new(grid, 1.0, 1.0, QuantumMethod::SplitStep).unwrap();
    s.origin = [-8.0, -8.0, 0.0];
    let v = harmonic(&s, 1.0);
    s.set_potential(v).unwrap();
    let (e, _) = s.eigenstates(6, 1e-9, 500).unwrap();
    for (en, exact) in e.iter().zip([1.0, 2.0, 2.0, 3.0, 3.0, 3.0]) {
        assert!((en - exact).abs() < 1e-6, "{e:?}");
    }
}

/// Transmission probability through a rectangular barrier of height `v0` and width `a` at
/// energy `e` (ħ = m = 1).
fn barrier_transmission(e: f64, v0: f64, a: f64) -> f64 {
    if (e - v0).abs() < 1e-12 {
        return 1.0 / (1.0 + a * a * v0 / 2.0);
    }
    let q = (2.0 * (e - v0).abs()).sqrt();
    let s = if e < v0 {
        (q * a).sinh().powi(2)
    } else {
        (q * a).sin().powi(2)
    };
    1.0 / (1.0 + v0 * v0 * s / (4.0 * e * (e - v0).abs()))
}

#[test]
fn barrier_transmission_matches_the_analytic_coefficient() {
    let (v0, width, sigma) = (1.0, 1.0, 8.0);
    for k0 in [1.1, 1.414, 1.8] {
        let n = 8192;
        let length = 600.0;
        let mut s = line(n, length, QuantumMethod::SplitStep);
        // The barrier covers exactly `width` worth of grid cells, centred on zero.
        let h = s.grid.h;
        let v: Vec<f64> = (0..n)
            .map(|i| {
                let x = s.coordinates(i)[0];
                if x.abs() < 0.5 * width - 0.25 * h {
                    v0
                } else {
                    0.0
                }
            })
            .collect();
        let cells = v.iter().filter(|&&x| x > 0.0).count();
        let a = cells as f64 * h;
        s.set_potential(v).unwrap();
        s.set_gaussian([-80.0, 0.0, 0.0], [sigma, 1.0, 1.0], [k0, 0.0, 0.0])
            .unwrap();
        s.step(0.05, (130.0 / k0 / 0.05) as usize).unwrap();
        let transmitted: f64 = (0..n)
            .filter(|&i| s.coordinates(i)[0] > 0.5 * width)
            .map(|i| s.psi[i].0.powi(2) + s.psi[i].1.powi(2))
            .sum::<f64>()
            * h;
        // Average T(k) over the packet's momentum distribution ∝ exp(-2σ²(k - k0)²).
        let (mut num, mut den) = (0.0, 0.0);
        for j in 0..4001 {
            let k = k0 + (j as f64 - 2000.0) * 1e-4 * 2.0;
            let w = (-2.0 * sigma * sigma * (k - k0).powi(2)).exp();
            num += w * barrier_transmission(0.5 * k * k, v0, a);
            den += w;
        }
        let expected = num / den;
        assert!(
            (transmitted - expected).abs() < 2e-3 * expected.max(0.05),
            "k0 = {k0}: transmitted {transmitted}, expected {expected}"
        );
    }
}

#[test]
fn moving_trap_drives_the_centre_like_a_classical_particle() {
    // V = ½ (x - a sin Ωt)²: Ehrenfest's theorem is exact for a harmonic potential.
    let (a, big_omega) = (1.0, 0.6);
    let classical =
        |t: f64| a / (1.0 - big_omega * big_omega) * ((big_omega * t).sin() - big_omega * t.sin());
    for (method, ord, tol) in [
        (QuantumMethod::SplitStep, 2, 1e-4),
        (QuantumMethod::SplitStep, 4, 1e-7),
        (QuantumMethod::CrankNicolson, 2, 5e-3),
    ] {
        let mut s = line(
            if method == QuantumMethod::SplitStep {
                256
            } else {
                801
            },
            24.0,
            method,
        );
        s.set_order(ord).unwrap();
        let xs: Vec<f64> = (0..s.grid.len()).map(|i| s.coordinates(i)[0]).collect();
        s.set_potential_fn(Box::new(move |t, v: &mut [f64]| {
            let c = a * (big_omega * t).sin();
            for (vi, x) in v.iter_mut().zip(&xs) {
                *vi = 0.5 * (x - c).powi(2);
            }
            Ok(())
        }))
        .unwrap();
        assert!(s.time_dependent());
        s.set_gaussian([0.0; 3], [0.5f64.sqrt(), 1.0, 1.0], [0.0; 3])
            .unwrap();
        s.step(0.01, 500).unwrap();
        let x = s.position().unwrap()[0];
        assert!(
            (x - classical(5.0)).abs() < tol,
            "{method:?}/{ord}: {x} vs {}",
            classical(5.0)
        );
        // The potential reported is the one at the current time.
        let c = a * (big_omega * s.t).sin();
        let i = s.grid.len() / 2;
        assert!((s.potential()[i] - 0.5 * (s.coordinates(i)[0] - c).powi(2)).abs() < 1e-12);
    }
}

#[test]
fn absorbing_layer_removes_outgoing_waves() {
    let mut s = line(512, 60.0, QuantumMethod::SplitStep);
    s.set_absorbing_layer(10.0, 2.0).unwrap();
    assert_eq!(s.absorber()[256], 0.0);
    assert!(s.absorber()[0] > 1.9);
    s.set_gaussian([0.0; 3], [1.0, 1.0, 1.0], [3.0, 0.0, 0.0])
        .unwrap();
    s.step(0.01, 2000).unwrap();
    // Without the layer the packet would wrap around and stay; with it almost nothing is left.
    assert!(s.norm() < 1e-3, "{}", s.norm());

    let mut s = line(601, 60.0, QuantumMethod::CrankNicolson);
    s.set_absorbing_layer(10.0, 2.0).unwrap();
    s.set_gaussian([0.0; 3], [1.0, 1.0, 1.0], [3.0, 0.0, 0.0])
        .unwrap();
    s.step(0.01, 2000).unwrap();
    assert!(s.norm() < 1e-3, "{}", s.norm());
}

#[test]
fn neumann_and_periodic_crank_nicolson() {
    // A packet bouncing in a box with reflecting (Neumann) walls keeps its norm, as does a
    // periodic ring, on both the direct (1D) and iterative (2D) solvers.
    for boundary in [Boundary::Neumann, Boundary::Periodic] {
        let grid = Grid::new(&[200], 0.1, boundary).unwrap();
        let mut s = Schrodinger::new(grid, 1.0, 1.0, QuantumMethod::CrankNicolson).unwrap();
        s.set_gaussian([10.0, 0.0, 0.0], [1.0, 1.0, 1.0], [4.0, 0.0, 0.0])
            .unwrap();
        let e0 = s.energy().unwrap();
        s.step(0.01, 1000).unwrap();
        assert!((s.norm() - 1.0).abs() < 1e-11, "{boundary:?}");
        assert!(
            (s.energy().unwrap() / e0 - 1.0).abs() < 1e-10,
            "{boundary:?}"
        );

        let grid = Grid::new(&[24, 24], 0.5, boundary).unwrap();
        let mut s = Schrodinger::new(grid, 1.0, 1.0, QuantumMethod::CrankNicolson).unwrap();
        s.set_gaussian([6.0, 6.0, 0.0], [1.0, 1.0, 1.0], [2.0, -1.0, 0.0])
            .unwrap();
        s.step(0.02, 100).unwrap();
        assert!((s.norm() - 1.0).abs() < 1e-9, "{boundary:?} 2D");
    }
}

#[test]
fn invalid_input_is_refused_and_failed_steps_roll_back() {
    let dirichlet = Grid::new(&[16], 0.1, Boundary::Dirichlet).unwrap();
    assert!(Schrodinger::new(dirichlet, 1.0, 1.0, QuantumMethod::SplitStep).is_err());
    let odd = Grid::new(&[12], 0.1, Boundary::Periodic).unwrap();
    assert!(Schrodinger::new(odd, 1.0, 1.0, QuantumMethod::SplitStep).is_err());
    assert!(Schrodinger::new(odd, 1.0, 1.0, QuantumMethod::CrankNicolson).is_ok());
    assert!(Schrodinger::new(odd, 0.0, 1.0, QuantumMethod::CrankNicolson).is_err());
    let mut cn = Schrodinger::new(dirichlet, 1.0, 1.0, QuantumMethod::CrankNicolson).unwrap();
    assert!(cn.set_order(4).is_err());
    assert!(cn.set_psi(vec![(1.0, 0.0); 3]).is_err());
    assert!(cn.normalize().is_err());
    assert!(cn.set_absorber(vec![-1.0; 16]).is_err());
    assert!(QuantumMethod::parse("euler").is_err());

    // Dirichlet edge values are forced to zero.
    cn.set_psi(vec![(1.0, 1.0); 16]).unwrap();
    assert_eq!(cn.psi[0], (0.0, 0.0));
    assert_eq!(cn.psi[15], (0.0, 0.0));

    // A potential function that fails part-way leaves the last completed step.
    let mut s = line(64, 10.0, QuantumMethod::SplitStep);
    s.set_gaussian([0.0; 3], [1.0, 1.0, 1.0], [1.0, 0.0, 0.0])
        .unwrap();
    s.set_potential_fn(Box::new(|t, v: &mut [f64]| {
        if t > 0.25 {
            return Err(SimError::Invalid("too late".into()));
        }
        v.iter_mut().for_each(|x| *x = t);
        Ok(())
    }))
    .unwrap();
    let err = s.step(0.1, 10).unwrap_err();
    assert!(err.to_string().contains("too late"));
    assert!((s.t - 0.2).abs() < 1e-12, "{}", s.t);
    let at_02 = s.psi.clone();
    assert!(s.step(0.1, 1).is_err());
    assert_eq!(s.psi, at_02);
}
