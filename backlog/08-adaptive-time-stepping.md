# 08 · Adaptive time stepping

**Priority:** P2 · **Size:** M · **Area:** Integrators

## Why
Highly eccentric orbits, close encounters and stiff transients need small steps only some of the time.
A fixed step wastes work or loses accuracy.

## Scope
- Embedded Runge–Kutta with error control: Dormand–Prince 5(4), with relative and absolute tolerances.
- Dense output (also feeds event location in item 05).
- `run` gains a mode where it integrates to a final time and records at requested output times
  rather than every k-th step.

## Notes
- Naive adaptive steps destroy the long-time behaviour of symplectic methods. For conservative
  problems consider time-symmetric step control or a time transformation (Sundman) instead; document
  the trade-off for users.

## Done when
- A Kepler orbit with eccentricity 0.99 reaches a target energy error with far fewer force
  evaluations than any fixed-step method in the engine.
