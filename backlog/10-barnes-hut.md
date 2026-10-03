# 10 · Tree gravity (Barnes–Hut)

**Priority:** P2 · **Size:** L · **Area:** Performance · **Status: Done** (speed target met only at θ = 1; see below)

> **Done** (`src/forces/tree.rs`, Python `TreeGravity(G, softening, theta=0.5, quadrupole=False)`). Octree rebuilt each
> evaluation (leaves of ≤ 8 particles), mass, centre of mass and traceless quadrupole per node. Acceptance: a node of
> size s is used as a multipole when its centre of mass is farther than s/θ from the whole receiving group's bounding
> box and the node does not overlap the group (which rules out self-interaction; θ is limited to (0, 1]). Group walk
> (Barnes 1990): groups of ≤ 32 particles share one interaction list, stored structure-of-arrays and summed in four
> independent lanes; groups are processed in deterministic parallel blocks, so results do not depend on the thread
> count and checkpoint restarts are bit-exact. Same `Force` interface as `NewtonianGravity`, including softening and
> massless test particles; the potential comes from the same expansion.
>
> Results (`tests/tree.rs`, `tests/python/test_tree.py`, `cargo run --release --example tree_gravity`, bench
> `tree_gravity_verlet_step`), Plummer sphere N = 10⁵:
>
> | θ | monopole rms error | time | speedup | quadrupole rms error | time | speedup |
> |---|---|---|---|---|---|---|
> | 0.2 | 8.2e-5 | 2546 ms | 3× | 8.6e-6 | 6312 ms | 1× |
> | 0.3 | 2.5e-4 | 1110 ms | 7× | 4.0e-5 | 2897 ms | 3× |
> | 0.5 | 8.7e-4 | 389 ms | 20× | 2.6e-4 | 1018 ms | 8× |
> | 0.7 | 2.5e-3 | 190 ms | 41× | 1.2e-3 | 473 ms | 17× |
> | 1.0 | 8.2e-3 | 113 ms | 69× | 5.7e-3 | 267 ms | 30× |
>
> (direct summation: 7.9 s per evaluation on the same 4-thread machine.)
> - Error scaling: monopole error grows as θ^2.9, quadrupole as θ^4.0 — one power of θ per multipole order, as expected.
> - The ≥ 50× target at N = 10⁵ is met at θ = 1 (rms 0.8%) but not at the more usual θ = 0.5 (20×). Profiling: ~2 300
>   interactions per particle at θ = 0.5, at about the same cost per interaction as a direct-sum pair (~6 ns per thread),
>   so the speedup is close to its algorithmic bound N / (2 × interactions). The next gain needs explicit SIMD (the
>   inner loop is a pure reduction, unlike the direct sum's scatter), or larger N (the ratio grows ~linearly with N).
>   The quadrupole loop is the slower of the two and only pays off when you need errors below ~3e-4.
> - Getting here: per-particle walks (8 µs/particle) → group walk (2×) → dropping the conservative δ term from the
>   acceptance test in favour of an overlap check (1.2×) → structure-of-arrays lanes (1.25×).
> - In a 1 500-body cluster run, energy drifts by < 1e-3 over 300 steps and momentum changes by < 1e-4.
>
> Not done: SIMD intrinsics, fast multipole method, tree reuse between steps, and an analytic Jacobian (Lyapunov
> exponents fall back to finite differences).

## Why
Direct summation is O(N²). Barnes–Hut is O(N log N) and makes 10⁵–10⁶ bodies practical.

## Scope
- Octree build each step; opening angle θ as a parameter; monopole first, then quadrupole moments.
- Same `Force` interface as `NewtonianGravity`, so it is a drop-in swap.
- Parallel tree walk (item 04).

## Notes
- Tree forces are not exactly pairwise-antisymmetric, so momentum is no longer conserved to round-off;
  document the expected error versus θ.
- Fast multipole method is a possible follow-up if exact momentum conservation matters.

## Done when
- Force error versus direct sum is measured as a function of θ and matches expected scaling.
- N = 10⁵ step time is at least 50× faster than direct sum in the benchmark.
