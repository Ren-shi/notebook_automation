# 20 · Fields and continua

**Priority:** P3 · **Size:** L · **Area:** Physics

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
