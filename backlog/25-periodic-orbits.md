# 25 · Periodic orbits: finding and continuation

**Priority:** P2 · **Size:** M · **Area:** Analysis

## Why
Periodic orbits organise the dynamics of most systems. Finding them, following them as a parameter changes, and
reading off their stability is a standard way to test an idea.

## Scope
- Newton shooting on the period map, using the monodromy matrix from the variational equations (item 18), with a
  finite-difference fallback.
- Symmetry-reduced shooting (e.g. half-period with a reflection) for speed and robustness.
- Pseudo-arclength continuation in a parameter, detecting folds and period-doubling.
- Floquet multipliers and a stable/unstable classification.

## Depends on
18.

## Done when
- The figure-eight three-body orbit is found from a rough guess to 1e-10 and its Floquet multipliers match the
  published values.
- Continuation of a driven Duffing oscillator's periodic branch reproduces the known fold points.
