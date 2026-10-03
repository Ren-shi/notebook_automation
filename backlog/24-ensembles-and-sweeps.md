# 24 · Ensembles and parameter sweeps

**Priority:** P2 · **Size:** M · **Area:** Analysis

## Why
Many questions need hundreds or thousands of runs: basins of attraction, stability diagrams, ensemble averages,
sensitivity to initial conditions. Looping over `World` objects in Python is serial and repetitive.

## Scope
- `physim.ensemble`: build many copies of a world from a template with varied initial conditions or parameters.
- Run them in parallel in Rust (rayon over worlds, with per-world parallelism turned off).
- Per-run reductions (final state, event times, a user summary) so memory stays small for large sweeps.
- Results as NumPy arrays shaped by the sweep grid.

## Done when
- A 100×100 sweep of a driven pendulum produces a basin-of-attraction map, and runs at least 0.8× the core count
  faster than a serial loop.
- Results are bit-identical between serial and parallel execution.
