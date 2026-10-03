# 02 · Benchmarks

**Priority:** P1 · **Size:** S · **Area:** Tooling · **Status: Done**

> **Done.** Rust benchmarks in `benches/engine.rs` (`cargo bench`, criterion, HTML reports in
> `target/criterion/`); Python overhead in `benches/python_overhead.py`. Baseline below. Accuracy per unit cost
> is covered by the convergence plot in the example notebook rather than a timing benchmark.

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

## Baseline

Commit `649e586` + benchmarks, release build, single thread, 4-vCPU Intel Xeon @ 2.80 GHz (cloud VM;
expect ±5–10% run-to-run). Median of criterion's estimate.

### Rust (`cargo bench --bench engine`)

| Benchmark | Time per step | Notes |
|---|---|---|
| Gravity, Verlet, N = 100 | 62 µs | 80 M pair interactions/s per step |
| Gravity, Verlet, N = 1 000 | 5.9 ms | 84 M pairs/s |
| Gravity, Verlet, N = 5 000 | 154 ms | 81 M pairs/s; each step is 2 force evaluations |

Integrator cost on the same system:

| Integrator | Force evals / step | Spring chain, N = 10 000 | Gravity, N = 200 |
|---|---|---|---|
| explicit_euler | 1 | 252 µs | 122 µs |
| symplectic_euler | 1 | 244 µs | 122 µs |
| verlet | 2 | 426 µs | 250 µs |
| rk4 | 4 | 1.01 ms | 481 µs |
| yoshida4 | 6 | 1.20 ms | 740 µs |

Recording, 1 000 Verlet steps of a 1 000-particle chain:

| `record_every` | Time |
|---|---|
| 1 000 (ends only) | 46.8 ms |
| 1 (every step) | 107.8 ms |

### Python (`python benches/python_overhead.py`)

| Benchmark | µs / step |
|---|---|
| Gravity N = 10, Rust loop `w.step(dt, n)` | 0.64 |
| Gravity N = 10, Python loop `w.step(dt)` | 0.80 |
| Gravity N = 10, built-in / `CustomForce` (numpy) | 0.66 / 34.6 |
| Gravity N = 100, built-in / `CustomForce` | 59 / 840 |
| Gravity N = 1 000, built-in / `CustomForce` | 5 900 / 72 300 |

### What the numbers say
- **Cost tracks force evaluations.** Step times follow 1 : 2 : 4 : 6 for Euler, Verlet, RK4 and Yoshida. Reusing
  the end-of-step acceleration (item 03) should make Verlet cost the same as Euler, and Yoshida 4 instead of 6
  evaluations.
- **Recording every step more than doubles run time**, mostly because each frame evaluates the potential energy:
  an extra pass over all forces, O(N²) for gravity. Item 06 should make per-frame energies optional.
- **Index-based forces are slow in bulk.** Each `Spring` is its own boxed force, about 21 ns per spring per
  evaluation; a single spring-network force holding all bonds in arrays should be several times faster (item 19).
- **Python loop overhead is negligible** (~0.16 µs per step), so driving runs from Python is fine.
- **`CustomForce` costs ~17 µs per force evaluation in fixed overhead**, and numpy direct-sum gravity is 12–14×
  slower than the Rust force. Fine for prototyping; port to Rust for production runs.
