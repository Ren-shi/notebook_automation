# 19 · More built-in forces

**Priority:** P2 · **Size:** M · **Area:** Physics

## Why
Several common forces exist only as Python `CustomForce` code, which is slow. They belong in Rust.

## Scope
- Damped spring (spring + dashpot along the bond).
- Central potentials: power law, Yukawa, Plummer and Hernquist halos, isotropic and anisotropic harmonic traps.
- Time-dependent drives: sinusoidal forcing, parametric modulation of a spring constant.
- Orbital perturbations: 1PN (relativistic) correction, J2 oblateness of the central body.
- A Rust-side "closure force" that takes a function or closure, so new forces in Rust need no boilerplate struct.
- `SpringNetwork`: all bonds of a lattice or polymer in one force, stored as index and parameter arrays.
  Separate `Spring` objects cost ~21 ns per spring per evaluation (item 02 baseline) from dynamic dispatch and
  per-spring checks.
- Expose each one to Python.

## Done when
- Each force has an analytic check: a driven damped oscillator's resonance curve; 1PN perihelion advance of
  6πGM/(c²a(1−e²)) per orbit; J2 nodal regression rate.
