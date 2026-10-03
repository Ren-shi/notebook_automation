# 26 · Conserved-quantity monitors

**Priority:** P2 · **Size:** S · **Area:** Analysis

## Why
Energy is the only invariant recorded today. Checking momentum, angular momentum, or an invariant you have
derived yourself is the quickest way to validate a new model.

## Scope
- Built-in observables: total momentum, total angular momentum about a point, centre of mass, virial.
- User observables: a vectorised Python function of (t, pos, vel, mass) evaluated at recorded frames only.
- Recorded in `Trajectory` alongside energies, and saved/streamed by item 06's writers.

## Done when
- A free N-body run reports momentum and angular momentum conserved to round-off with `verlet`.
- A user-defined invariant (e.g. the Laplace–Runge–Lenz vector for Kepler) is recorded and saved with the trajectory.
