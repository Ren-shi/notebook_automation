# 09 · Charged particles: Lorentz force and Boris pusher

**Priority:** P2 · **Size:** M · **Area:** Physics

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
