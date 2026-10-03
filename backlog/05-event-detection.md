# 05 · Event detection and stopping conditions

**Priority:** P1 · **Size:** M · **Area:** Core · **Status: Done**

> **Done** (`src/events.rs`, `World::run_with_events`, Python `ps.Event` and `run(..., events=[...])`). One change
> from the plan: instead of a Hermite interpolant, crossings are located by re-stepping the integrator from the start
> of the step with `theta * dt` (Illinois root finding to ~1e-12 of a step). Event states are therefore the
> integrator's own solution, as accurate as a normal step, with no interpolation error; events are rare, so the extra
> steps cost little. Built-in Rust events: radial velocity (periapsis/apoapsis), coordinate crossing, separation.
>
> Results: the notebook's periapsis advance now agrees with theory to **8.7e-10** relative (was 1.4e-3 from
> scanning frames). Kepler periapses/apoapses land on exact half-period multiples to 1e-9 and at a(1 ± e) to
> 1e-10; a projectile's ground hit under Verlet matches the analytic time to 1e-12 (`tests/events.rs`).

## Why
Many questions are about *when* something happens: periapsis passage, a zero crossing, a collision,
escape. Right now this is found by scanning recorded frames, so accuracy is limited by the sampling
interval (the precession example in the notebook is limited this way).

## Scope
- An event is a function `g(state) -> f64`; an event fires when `g` changes sign during a step.
- Locate the crossing to tolerance by root-finding (bisection or Brent) on a dense-output interpolant
  of the step (cubic Hermite from positions and velocities at both ends is enough).
- Per event: direction (rising, falling, both), and action (record, stop the run).
- Python: `w.run(..., events=[ps.Event(fn, direction=-1, terminal=False)])`, with event times and
  states returned on the trajectory. Built-in Rust events for common cases (radial distance extrema,
  coordinate crossings) avoid Python call cost.

## Done when
- The notebook's periapsis advance is measured from events and agrees with theory to better than 1e-6.
- A terminal event stops `run` at the crossing time within tolerance.
