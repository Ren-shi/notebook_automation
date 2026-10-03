# 17 · Higher-order and specialised integrators

**Priority:** P2 · **Size:** M · **Area:** Integrators

## Why
The engine has orders 1, 2 and 4. Precision studies and long planetary integrations need higher order or
schemes built for the problem structure.

## Scope
- Higher-order symplectic compositions: Yoshida 6 and 8, and optimised 4th-order schemes with smaller error
  constants (Forest–Ruth, Omelyan, Blanes–Moan).
- Wisdom–Holman (Kepler drift + interaction kick) for near-Keplerian systems; orders of magnitude cheaper for
  planetary problems.
- Implicit midpoint / Gauss–Legendre Runge–Kutta: symplectic even when forces depend on velocity.
- A generic composition integrator taking a list of coefficients, so new schemes are one line to define.

## Done when
- The convergence-order test covers each new scheme at its stated order.
- Wisdom–Holman integrates the outer solar system for 10⁶ years with bounded energy error at a step of
  ~1/20 of Jupiter's period.
