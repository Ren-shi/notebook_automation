# 10 · Tree gravity (Barnes–Hut)

**Priority:** P2 · **Size:** L · **Area:** Performance

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
