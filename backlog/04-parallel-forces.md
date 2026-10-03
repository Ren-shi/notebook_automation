# 04 · Parallel force evaluation

**Priority:** P1 · **Size:** M · **Area:** Performance

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
