# 32 · Quantum mechanics: the Schrödinger equation on a grid

**Priority:** P3 · **Size:** L · **Area:** Physics · **Status: Done**

> **Done** (`src/fields/quantum.rs`, `src/python/quantum.rs`, `python/physim/quantum.py`; Python `Schrodinger`,
> `physim.quantum`, `plot.wavefunction`, `plot.animate_wavefunction`; `tests/quantum.rs`,
> `tests/python/test_quantum.py`, `notebooks/08_quantum.ipynb`, `docs/theory/quantum.md`).
> - Split-step Fourier (periodic, power-of-two sizes; Strang or Yoshida order 4) and Crank-Nicolson (any boundary;
>   Thomas/Sherman-Morrison in 1D, Jacobi-preconditioned BiCGSTAB in 2D/3D), static or time-dependent potentials
>   (arrays or Python functions `V(t, x[, y[, z]])`), quadratic absorbing layers, Gaussian packets, `run` recording.
> - Bound states by LOBPCG on the same discretisation as the time stepping (Fourier preconditioner for spectral,
>   exact tridiagonal solve in 1D, inner CG otherwise), instead of the imaginary-time evolution the plan named.
>   LOBPCG was 8-60× faster than the subspace inverse iteration tried first, at the same accuracy.
> - Observables: norm, ⟨x⟩, variance, ⟨p⟩, ⟨T⟩, ⟨V⟩, energy; helpers for |φ(p)|², probability current, region
>   probabilities, analytic barrier transmission (plane wave and packet-averaged) and oscillator/box levels.
>
> Results (all in the tests):
> - Norm conserved to < 1e-11 over 1e5 steps; Crank-Nicolson energy to < 1e-10, split-step < 1e-6 with no drift.
> - Free packet width and centre exact to 1e-9 (split-step); Crank-Nicolson O(h²) dispersion error, measured order
>   > 1.8. Coherent state follows x₀ cos ωt to 1e-7 with order 4.
> - Time orders: Strang 2.0, Yoshida 4.0, Crank-Nicolson 2.0. Spectral oscillator levels to 1e-9 (2D degeneracies
>   correct); finite-difference oscillator and box levels converge at order 2.0.
> - Barrier transmission within 0.15% of the packet-averaged analytic value at three energies (below, near and above
>   the barrier). Driven trap and Ehrenfest checks against the classical solution and a `World` run.
> - Speed (12 logical cores): split-step 0.6 ms/step at 4096 points, 10 ms at 256², 16 ms at 64³; ten 1D eigenstates in
>   0.1-0.3 s.
>
> Left out (as planned): many-body wavefunctions, spin and two-level systems, relativistic equations. Also: the
> split-step kinetic factor needs `dt ≲ 2mh²/πħ` with sharp potentials (documented), and Crank-Nicolson in 2D/3D is
> not multithreaded.

## Why
Everything in the engine is classical. Wave packets, tunnelling, bound states and the quantum-classical
correspondence are standard topics that need a different state: a complex wavefunction ψ(x, t) on a grid rather
than particle positions and velocities. Item 20 already provides most of the machinery (grids in 1-3D, boundary
conditions, a complex FFT, Crank-Nicolson), and the Schrödinger equation is the heat equation in imaginary time, so
this extends `fields` rather than starting from scratch.

## Scope
- Separate solver beside `WaveEquation` and `HeatEquation`, not part of `World`: `Schrodinger` (Rust in
  `src/fields/`, Python binding) for one particle in 1, 2 or 3 dimensions, with ħ and mass as parameters.
- Potential V(x) as an array, or V(x, t) as a vectorised Python function evaluated once per step.
- Time evolution:
  - Split-step Fourier (Strang splitting, second order; optionally Yoshida composition for fourth order) on periodic
    grids, reusing `fft_nd`. Exactly unitary.
  - Crank-Nicolson on Dirichlet (hard-wall) grids. Unitary, second order. Needs a complex solver: Thomas algorithm in
    1D, BiCGSTAB or complex CG on the (complex-symmetric) system in 2-3D, since the existing CG is real-only.
  - An optional complex absorbing potential at the edges to stop outgoing waves reflecting.
- Bound states by imaginary-time evolution with Gram-Schmidt against lower states (ground state plus the first few
  excited states), returning energies and wavefunctions.
- Observables: norm, ⟨x⟩, ⟨p⟩, ⟨x²⟩, ⟨p²⟩, energy, probability current, momentum-space density (one FFT), and
  transmitted/reflected probability across a plane.
- Initial states: Gaussian wave packet (centre, width, mean momentum), eigenstates from the solver above, or any
  complex array.
- Output as complex NumPy arrays recorded every k steps, savable like other field data; `physim.plot` helpers for
  |ψ|² and phase-coloured plots and animations.
- Left out: many-body wavefunctions (cost grows exponentially with particle count), spin and two-level systems
  (could be a later item), relativistic equations (Dirac, Klein-Gordon), and density functional theory.

## Depends on
20.

## Done when
- Norm is conserved to round-off (split-step and Crank-Nicolson) over 10⁵ steps.
- A free Gaussian packet's width follows σ(t) = σ₀ √(1 + (ħt / 2mσ₀²)²) and its centre moves at p₀/m.
- A coherent state in a harmonic oscillator keeps its shape and its ⟨x⟩ follows the classical trajectory.
- Harmonic-oscillator eigenvalues come out as (n + ½)ħω and the infinite square well as n²π²ħ²/(2mL²), each
  converging at the expected order in grid spacing.
- Transmission through a rectangular barrier matches the analytic coefficient for several energies above and below
  the barrier height.
- Measured convergence orders in time: 2 for Strang splitting and Crank-Nicolson (4 if Yoshida composition is added).
- Ehrenfest check: in a smooth anharmonic potential, ⟨x⟩(t) agrees with a classical `World` trajectory for a narrow
  packet at short times and visibly departs once the packet spreads.
- A notebook covering wave-packet spreading, tunnelling, eigenstates, and the quantum vs classical comparison, plus
  a theory note in `docs/theory/`.
