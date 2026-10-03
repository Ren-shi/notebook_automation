# 28 · Richer constraints

**Priority:** P3 · **Size:** L · **Area:** Physics

## Why
Item 07 covers fixed-length rods. Many systems need other constraint types.

## Scope
- One-sided constraints: strings that go slack, particles resting on a surface (complementarity / impulse at
  activation).
- Surface and curve constraints: a bead on a wire, a particle on a sphere or a general implicit surface g(r) = 0.
- An O(links) direct solver for chains and trees (follow-up noted in item 07).
- Nonholonomic (rolling) constraints, if the rigid-body work in item 13 makes them useful.

## Depends on
07.

## Done when
- A pendulum on a string that goes slack and re-tightens matches the analytic slack trajectory.
- A bead on a rotating hoop shows the expected bifurcation of the equilibrium at ω² = g/R.
- A 1 000-link chain steps at least 10× faster than the current CG solver.
