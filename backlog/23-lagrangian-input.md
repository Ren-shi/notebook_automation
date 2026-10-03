# 23 · Equations of motion from a Lagrangian or Hamiltonian

**Priority:** P2 · **Size:** L · **Area:** Usability

## Why
Testing a new model usually starts from a Lagrangian or Hamiltonian, not from Cartesian forces. Rewriting it as
forces by hand is slow and error-prone, and many systems (bead on a rotating hoop, spherical pendulum, coupled
oscillators in normal coordinates) are far more natural in generalised coordinates.

## Scope
- `physim.mechanics`: take L(q, q̇, t) or H(q, p, t) as a SymPy expression, derive the equations of motion
  symbolically, and generate fast code for them (SymPy `lambdify` first; optional generated Rust or C later).
- A generalised-coordinate state alongside the particle state, integrated by the existing integrators
  (symplectic ones for separable H; implicit midpoint / Gauss–Legendre from item 17 for non-separable H).
- Conserved quantities from the symbolic form: energy function, and momenta conjugate to cyclic coordinates.

## Done when
- A spherical pendulum built from its Lagrangian matches the same system built from a particle and a rod (item 07).
- A non-separable Hamiltonian (e.g. a charged particle in a magnetic field in canonical momenta) conserves H to
  the integrator's stated order.
