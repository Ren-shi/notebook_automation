# 20 · Fields and continua

**Priority:** P3 · **Size:** L · **Area:** Physics · **Status: Done**

> **Done** (`src/fields/{mod,fft,solvers,pm}.rs`, `src/python/fields.rs`; Python `WaveEquation`, `HeatEquation`,
> `solve_poisson`, `ParticleMesh`; `tests/fields.rs`, `tests/python/test_fields.py`, `notebooks/06_fields.ipynb`,
> `docs/theory/fields.md`, `examples/particle_mesh.rs`).
> - Grids with 1-3 axes and Dirichlet (fixed edge nodes), Neumann (cell-centred, symmetric) or periodic boundaries;
>   second-order Laplacian; a conjugate-gradient solver for the shifted Laplacian.
> - Wave equation by leapfrog (CFL-checked), heat equation explicit (stability-checked) or Crank-Nicolson (CG), Poisson
>   by FFT (periodic, exact for the discrete operator) or CG (Dirichlet with boundary values, Neumann).
> - A radix-2 FFT acting on whole rows along the leading axes, with pruning for zero-padded data and deterministic
>   slab parallelism (bit-identical for any thread count, tested).
> - `ParticleMesh`: cloud-in-cell deposit/interpolation; isolated mode by zero-padded FFT convolution with the softened
>   kernel (cached kernel transforms) or periodic mode by FFT Poisson solve; no self-force, momentum conserved to
>   round-off; checkpointable like any built-in force.
>
> Results:
> - Measured orders against analytic solutions: wave 1D Dirichlet and 2D periodic 2.0; heat explicit (dt ∝ h²) 2.0,
>   Crank-Nicolson 2.0 in time (vs the semi-discrete solution) and 2.0 in space; Poisson 3D periodic, 2D Dirichlet with
>   non-zero boundary values and 1D Neumann 2.0. Crank-Nicolson runs at 30-60× the explicit step limit.
> - Particle mesh vs direct summation (same softening): two-body force -11% at 2 cells, -4% at 4, -1% at 8, -0.1% at
>   32; Plummer sphere (N = 4000, 64³) median error < 2%, potential energy within 1%; total momentum zero to 1e-12.
>   Periodic mode reproduces the linear plane-wave field 4πGρ̄A sin(kx) to 3%.
> - Speed at N = 10⁵ (4 cores): 32 ms (32³) and 169 ms (64³) per evaluation against 270 ms for Barnes-Hut θ = 0.5, at
>   median errors 4.5% and 1.1% against 0.5% and 0.2%: faster but coarser, as expected of particle mesh.
> - Notebook 06: a struck drum, explicit vs Crank-Nicolson diffusion (15 CN steps vs 1689 explicit), a dipole in a
>   grounded box, and the cold collapse of a uniform sphere following the analytic cycloid (half-mass radius 0.854 vs
>   0.836 at t_ff/2 on a 32³ grid).
>
> Left out: multigrid (CG and FFT cover the stated scope), non-uniform spacing, mixed boundary conditions per face,
> and P³M (short-range particle-particle correction to particle mesh).

## Why
Everything so far is particles. Waves, diffusion and fields on a grid were among the original options, and some
ideas need them, either alone or coupled to particles.

## Scope
- Separate module: grids in 1D, 2D and 3D, finite differences, boundary conditions.
- Solvers: wave equation (leapfrog), heat equation (explicit and Crank–Nicolson), Poisson (multigrid or FFT).
- Coupling to particles: particle-mesh gravity or electrostatics (deposit charge, solve Poisson, interpolate the
  force back), which also gives an O(N log N) alternative to Barnes–Hut for smooth distributions.

## Done when
- Each solver converges at its stated order against an analytic solution.
- Particle-mesh gravity reproduces direct-sum forces to the expected grid accuracy.
