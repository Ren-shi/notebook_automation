# 03 · Reuse forces between Verlet steps

**Priority:** P1 · **Size:** S · **Area:** Performance · **Status: Done**

> **Done.** `AccelCache` in `src/integrators.rs` stores the last acceleration with a snapshot of its inputs
> (time, positions, masses and pinned flags, compared bit-for-bit) and the force set's version counter, which every
> change to the forces bumps. Validity is checked against the actual inputs rather than invalidated by setters, so
> direct edits to `World.state` and rollback after a failed step are handled automatically. Reuse is skipped when
> any force is velocity-dependent. Results are bit-identical to two fresh evaluations per step
> (`tests/force_reuse.rs`). Forces must be pure functions of their inputs; this is now documented on `Force`.
>
> | Benchmark | Before | After |
> |---|---|---|
> | Gravity N = 1 000, Verlet step | 5.9 ms | 3.05 ms (−47%) |
> | Gravity N = 200: Verlet / Yoshida4 | 250 / 740 µs | 126 / 394 µs (Verlet now equals Euler) |
> | Spring chain N = 10 000: Verlet / Yoshida4 | 426 / 1 200 µs | 300 / 805 µs (−30%) |
>
> The spring chain gains less because its forces are so cheap that copying and comparing the position snapshot
> (~55 µs for 10 000 particles) is a visible share of the cost.

## Why
`kdk` in `src/integrators.rs` evaluates accelerations at the start and end of every step. For
velocity-independent forces the end-of-step acceleration equals the next step's start-of-step one, so
half the force evaluations are wasted. Force evaluation dominates N-body cost, so this is close to a 2×
speedup for `verlet` and saves 2 of 6 evaluations per `yoshida4` step.

Measured (item 02 baseline): Verlet costs 1.7–2.0× symplectic Euler per step, and Yoshida 4 is 4.9–6× Euler, in line with 2 and 6 force evaluations.

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
