# 08 · Adaptive time stepping

**Priority:** P2 · **Size:** M · **Area:** Integrators · **Status: Done**

> **Done** (`src/adaptive.rs`, `World::run_adaptive`/`run_adaptive_into`, Python `w.run_adaptive(t_end, ...)`, and a
> fixed-step `dopri5` integrator). Dormand-Prince 5(4) with Hairer's error norm, PI step-size controller and
> automatic first step; first-same-as-last, so 6 force evaluations per step. Output at every accepted step or at
> requested times from the fourth-order continuous extension (no extra force evaluations, and the steps taken do not
> depend on the output times). Events are located on the continuous extension. Integrating backwards works.
> Streaming to a sink works as in `run`, and `traj.metadata["adaptive"]` reports accepted/rejected steps and force
> evaluations.
>
> Results (`tests/adaptive.rs`, `tests/python/test_adaptive.py`, notebook §7,
> `cargo run --release --example eccentric_orbit`):
> - e = 0.99 Kepler orbit, 3 orbits: force evaluations to reach the same maximum energy error.
>
>   | energy error | adaptive | verlet | yoshida4 | rk4 | dopri5 (fixed) |
>   |---|---|---|---|---|---|
>   | 1.1e-3 | 2 210 | 1 853× | 348× | 463× | 348× |
>   | 4.1e-6 | 4 160 | 15 754× | 738× | 492× | 369× |
>   | 1.2e-8 | 9 374 | > 3e7 steps | 1 311× | 1 748× | 655× |
>   | 1.8e-10 | 23 504 | > 3e7 steps | 2 091× | 1 394× | 523× |
> - The step varies by a factor of ~5 700 over one e = 0.99 orbit (1.3e-5 at periapsis, 0.077 at apoapsis).
> - Global error tracks the tolerance (harmonic oscillator, rtol 1e-6/1e-9/1e-12). Dense output matches the analytic
>   damped oscillator to 1e-8 at rtol 1e-10. Periapsis events land on the exact times to 1e-8, and a terminal
>   coordinate event matches Kepler's equation to 1e-9. A head-on free fall stops with a step-size-underflow error at
>   t = π/4 (the collision time) to 3e-11.
> - The symplectic trade-off, measured: at rtol 1e-10 the adaptive energy error grows linearly, 1.6e-10 per orbit at
>   e = 0.5 (1.6e-6 after 10 000 orbits), while `yoshida4` with the same number of evaluations stays at 1e-12. At
>   e = 0.9, `yoshida4` at that cost has errors of 1e-2 during each periapsis passage that then recover to 1e-9. The
>   README says when to use which.
>
> Not done: time-symmetric step control and a Sundman time transformation (adaptive steps that keep long-time
> energy behaviour). Worth revisiting with Wisdom-Holman in item 17 if long eccentric integrations are needed.
> Constraints are not supported in adaptive mode (rejected with a clear error).

## Why
Highly eccentric orbits, close encounters and stiff transients need small steps only some of the time.
A fixed step wastes work or loses accuracy.

## Scope
- Embedded Runge–Kutta with error control: Dormand–Prince 5(4), with relative and absolute tolerances.
- Dense output (also feeds event location in item 05).
- `run` gains a mode where it integrates to a final time and records at requested output times
  rather than every k-th step.

## Notes
- Naive adaptive steps destroy the long-time behaviour of symplectic methods. For conservative
  problems consider time-symmetric step control or a time transformation (Sundman) instead; document
  the trade-off for users.

## Done when
- A Kepler orbit with eccentricity 0.99 reaches a target energy error with far fewer force
  evaluations than any fixed-step method in the engine.
