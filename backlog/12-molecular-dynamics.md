# 12 · Molecular dynamics: pair potentials, periodic boxes, thermostats

**Priority:** P2 · **Size:** L · **Area:** Physics

## Why
Statistical-mechanics experiments (phase behaviour, transport, equilibration) need short-range pair
potentials in a periodic box and temperature control.

## Scope
- Generic pair potential with cutoff: Lennard-Jones, Morse, and a user-defined `V(r)` given as a table
  or a vectorised Python function evaluated once per step (avoids per-pair Python calls).
- Periodic boundary conditions with minimum-image convention; cell lists or neighbour lists for O(N).
- Thermostats: Langevin (needs a seeded RNG and a stochastic integrator such as BAOAB), Nosé–Hoover.
- Observables: temperature, pressure (virial), radial distribution function g(r).

## Done when
- NVE Lennard-Jones liquid conserves energy over 10⁵ steps with a shifted potential.
- Langevin and Nosé–Hoover runs reach the target temperature; g(r) matches published LJ data at a
  reference state point.
