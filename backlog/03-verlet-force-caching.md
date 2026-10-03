# 03 · Reuse forces between Verlet steps

**Priority:** P1 · **Size:** S · **Area:** Performance

## Why
`kdk` in `src/integrators.rs` evaluates accelerations at the start and end of every step. For
velocity-independent forces the end-of-step acceleration equals the next step's start-of-step one, so
half the force evaluations are wasted. Force evaluation dominates N-body cost, so this is close to a 2×
speedup for `verlet` and saves 2 of 6 evaluations per `yoshida4` step.

## Scope
- Cache the last acceleration and the state it was computed for.
- Invalidate when the cache cannot be trusted: forces added or removed, particles added, positions,
  velocities or time set from outside, integrator switched, or any force is velocity-dependent.

## Notes
- Simplest robust invalidation: a `generation` counter on `World` bumped by every mutating method, stored
  alongside the cache.
- Keep results bit-identical to the uncached version; the existing tests then cover correctness.

## Done when
- Force evaluations per `verlet` step drop from 2 to 1 for conservative systems (count them in a test).
- All existing tests pass unchanged; benchmark shows the speedup.
