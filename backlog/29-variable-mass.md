# 29 · Variable mass

**Priority:** P3 · **Size:** S · **Area:** Physics

## Why
Rockets, accretion and mass-loss problems need masses that change during a run.

## Scope
- Per-particle mass rate ṁ(t, state) and exhaust velocity, with the thrust term handled consistently.
- Integrators that update mass within a step; document which ones stay symplectic (none in general).

## Done when
- A rocket in free space matches the Tsiolkovsky equation to the integrator's order.
- A two-body orbit with slow isotropic mass loss expands as a ∝ 1/M (adiabatic invariant).
