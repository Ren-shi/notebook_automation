# 17 · Higher-order and specialised integrators

**Priority:** P2 · **Size:** M · **Area:** Integrators · **Status: Done**

> **Done** (`src/integrators/`: `splitting.rs`, `gauss.rs`, `wisdom_holman.rs`). New integrators: `yoshida6`,
> `yoshida8` (compositions of Verlet substeps, so they support rods), `pefrl` (Omelyan et al.), `blanes_moan4`
> (Blanes & Moan 6-stage splitting), `gauss2`/`gauss4`/`gauss6` (Gauss-Legendre implicit RK, fixed-point iteration to
> round-off), and `wisdom_holman` (democratic heliocentric coordinates with a universal-variable Kepler solver; other
> position-dependent forces such as `J2Oblateness` act as perturbations in the interaction kick). Generic schemes:
> `Composition`/`Splitting` in Rust, `w.use_composition(weights, order)` / `w.use_splitting(ops, order)` in Python,
> saved in checkpoints (new optional `integrator_scheme` field). `Integrator::name` now returns `&str` so custom
> schemes can carry their own names. Forest-Ruth is the same triple jump as `yoshida4` (an alias).
>
> Results (`tests/integrators.rs`, `tests/physics.rs`, `tests/python/test_integrators.py`, notebook §3 and §9,
> `cargo run --release --example outer_solar_system`):
> - Convergence orders measured on the harmonic oscillator: 6.0, 8.0 (7.7 over the notebook's wider step range),
>   4.0, 4.0, 2.0, 4.0, 6.0. The order test now uses steps suited to each order.
> - At equal force evaluations (harmonic oscillator), `pefrl` is 62× and `blanes_moan4` 189× more accurate than
>   `yoshida4`.
> - Outer solar system at 1/20 of Jupiter's period for 10⁶ years: `wisdom_holman` energy error 3.7e-6 in the first
>   10% and 4.3e-6 in the last 10% (bounded); `verlet` 2.0e-3, `yoshida4` 8.8e-4, `blanes_moan4` 4.2e-6 (at 6
>   evaluations per step). 5 s in release for 1.7 million steps.
> - `wisdom_holman` is second order with error proportional to the planet mass (10.0× smaller per decade of mass
>   below 1e-4); with `J2Oblateness` it reproduces the J2 nodal regression rate to 5e-3; restarts are bit-exact.
> - Gauss-Legendre conserves |v|² under a magnetic-type force `a = v × B` to 1e-13 over 5 000 steps (RK4 drifts) and
>   retraces its path when run backwards.
> - The Kepler solver had two bugs caught by its unit tests: a poor first guess on long eccentric steps and an
>   endless cycle at round-off; both fixed (mean-motion guess, residual-based convergence test). Rewriting the
>   Stumpff series with Horner and precomputed factorials halved its cost.
>
> Not done: a time-transformed (Sundman) or symmetric-adaptive variant for long eccentric integrations (see item 08),
> WHFast-style symplectic correctors, and Jacobi-coordinate Wisdom-Holman. The Gauss methods use plain fixed-point
> iteration, so stiff problems need a Newton solver (not implemented).

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
