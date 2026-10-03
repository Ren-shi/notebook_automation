# 14 · Scenario library and units

**Priority:** P3 · **Size:** M · **Area:** Usability

## Why
Setting up standard systems by hand is repetitive and error-prone, and the engine has no notion of units.

## Scope
- `physim.scenarios`: two-body orbit from orbital elements, solar system from tabulated elements,
  Plummer sphere, figure-eight three-body orbit, pendulum, spring lattice.
- Unit systems: SI, astronomical (AU, year, solar mass, with G chosen so G = 4π²), and dimensionless.
  Keep the engine unit-agnostic; do conversions in Python.

## Done when
- Each scenario has a test checking a known property (the figure-eight returns to its start after one
  period; Plummer sphere satisfies the virial theorem on average).
