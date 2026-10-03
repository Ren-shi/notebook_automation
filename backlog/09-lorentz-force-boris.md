# 09 · Charged particles: Lorentz force and Boris pusher

**Priority:** P2 · **Size:** M · **Area:** Physics · **Status: Done**

> **Done** (`src/forces/em.rs`, `src/integrators/boris.rs`, charges in `State`). Per-particle charges
> (`World::set_charges`/`set_charge`, Python `add_particle(..., charge=q)` and `w.charges`), saved in checkpoints (older
> checkpoints load as neutral) and run metadata. Forces see charges through new default trait methods
> (`Force::accumulate_charged`, `potential_charged`, `jacobian_vector_charged`), so no existing force changed; the
> acceleration cache keys on charges too. New forces: `ElectricField`, `MagneticField`, `Coulomb` (the gravity pair
> loop generalised to source weights and per-particle response, so gravity's results are bit-identical and Coulomb
> gets the same parallel blocking), and `FieldFunctions` (Rust closures; Python `FieldForce(E=f, B=g)`). The `boris`
> integrator splits forces through `Force::magnetic_field`/`accumulate_electric` and rejects velocity-dependent
> non-magnetic forces.
>
> Results (`tests/em.rs`, `tests/python/test_em.py`, notebook §11):
> - Gyration in uniform B with `boris`, h = 0.01, 10⁶ steps: |v| constant to 4e-14; radius m v/(qB) to 1.4e-13; the
>   velocity turns by exactly 2 atan(ωh/2) per step (ω = qB/m), i.e. a frequency error of (ωh)²/12 as predicted.
>   RK4 at the same step loses energy (2% by t = 10⁴ at h = 0.2, notebook §11).
> - E×B drift: mean velocity over whole numerical gyrations equals E×B/B² to 4e-16.
> - Coulomb: Rutherford deflection tan(θ/2) = kqQ/(m v∞² b) to 7.6e-6 rad (the finite start/end distance);
>   gravity and Coulomb cancel to round-off when G m² = k q²; momentum and energy conserved.
> - Magnetic mirror (Python-free `FieldFunctions` closure): the particle bounces between z = ±2 with μ = m v⊥²/(2B)
>   constant to 1.6e-4.
> - Analytic Jacobians for the new forces match finite differences (so `lyapunov` works with fields); checkpoint
>   restarts with charges and EM forces are bit-exact.
>
> Not done: a relativistic pusher (item 30), self-consistent fields (particle-in-cell, item 20), magnetic dipole
> interactions, and a Boris variant that accepts drag.

## Why
Charged-particle motion in electric and magnetic fields is a standard test bed, and the magnetic force
is velocity-dependent, so the current symplectic integrators handle it poorly.

## Scope
- Per-particle charge in `State`.
- Forces: uniform or user-defined E(r, t) and B(r, t); pairwise Coulomb interaction (shares code with gravity).
- Boris integrator: exactly conserves kinetic energy in a pure magnetic field and is the standard
  scheme for this problem.

## Done when
- Gyration in uniform B: radius and frequency match qB/m; kinetic energy constant to round-off over 10⁶ steps.
- E×B drift velocity matches E×B/B².
