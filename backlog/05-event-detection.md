# 05 · Event detection and stopping conditions

**Priority:** P1 · **Size:** M · **Area:** Core

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
