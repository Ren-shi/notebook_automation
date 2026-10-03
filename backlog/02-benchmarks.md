# 02 · Benchmarks

**Priority:** P1 · **Size:** S · **Area:** Tooling

## Why
Items 03, 04 and 10 are performance work. Without a baseline there is no way to tell if they helped.

## Scope
`benches/` using `criterion`:
- Direct-sum gravity step for N = 100, 1 000, 5 000
- Spring chain of 10 000 particles (cheap per-particle forces, measures integrator overhead)
- Each integrator on the same system (cost per step and cost per unit accuracy)
- Python overhead: `World.step` from Python vs pure Rust, and `CustomForce` call cost

## Done when
- `cargo bench` produces a report; numbers for the current engine are recorded in this file as the baseline.
