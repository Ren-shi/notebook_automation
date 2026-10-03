# 04 · Parallel force evaluation

**Priority:** P1 · **Size:** M · **Area:** Performance · **Status: Done**

> **Done**, with a different design from the notes below: instead of parallelising over target particle i (twice
> the pair work), the symmetric i < j loop is kept and split into row blocks of roughly equal pair counts
> (`src/parallel.rs`). Each block writes a buffer covering particles `lo..n`; buffers are summed in block order.
> Block boundaries depend only on N, so results are **bit-identical for any thread count and with the
> `parallel` feature off** (tested with 1 vs 4 threads). Per-particle forces (`UniformField`, `LinearDrag`)
> run in parallel above 8 192 particles. Gravity potential energy is blocked the same way. CI also tests the
> serial build.
>
> Gravity, one Verlet step, 4-vCPU Xeon @ 2.8 GHz:
>
> | N | Serial (before) | 1 thread | 4 threads |
> |---|---|---|---|
> | 100 | 31 µs | 31 µs | 31 µs (single block) |
> | 1 000 | 3.02 ms | 3.30 ms | 1.09 ms (2.8×) |
> | 5 000 | 84.4 ms | 92.4 ms | 23.6 ms (3.6×) |
>
> Single-thread cost is 7–10% above the old serial loop (buffer zeroing and the final sum). `Spring` objects are
> still evaluated one by one; bulk bonds belong to the `SpringNetwork` force in item 19.

## Why
Direct-sum gravity is O(N²) and single-threaded. Multi-core gives a near-linear speedup for N in the
thousands with no change in physics.

## Scope
- `rayon` behind a `parallel` feature (on by default).
- `NewtonianGravity`: parallelise over target particle i, summing over all j. This does twice the pair
  work of the current symmetric i<j loop but has no write conflicts; it wins once cores ≥ 3.
- Per-particle forces (fields, drag, anchors) in parallel over particles.

## Notes
- Summation order changes, so results differ at round-off level; compare with tolerances, not exactly.
- Keep a deterministic serial path (feature off or a runtime switch) for reproducible runs.

## Done when
- Gravity with N = 5 000 scales roughly with core count in the benchmark.
- Momentum conservation test passes with a tolerance appropriate to the non-symmetric sum.
