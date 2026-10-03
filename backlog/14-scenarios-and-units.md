# 14 · Scenario library and units

**Priority:** P3 · **Size:** M · **Area:** Usability · **Status: Done**

> **Done** (`python/physim/scenarios.py`, `python/physim/units.py`; pure Python, the engine stays unit-agnostic).
> - Scenarios return configured `World`s: `two_body` (orbital elements, barycentric), `solar_system` (Sun and any of
>   Mercury-Pluto from JPL's J2000 mean elements, in any unit system, Wisdom-Holman by default), `plummer_sphere`
>   (Aarseth-Hénon-Wielen sampling, seeded), `figure_eight`, `pendulum` (rod to a fixed pivot), `spring_lattice`
>   (1D/2D/3D, optional face diagonals, free or fixed boundary). Helpers `elements_to_state`, `state_to_elements`.
> - Units: `SI`, `ASTRONOMICAL` (AU, Gaussian year, the solar mass implied by Gauss's constant, so G = 4π² exactly),
>   `DIMENSIONLESS`, and `scaled(...)` for custom systems; conversion by named quantity or (L, M, T) powers.
>
> Results (`tests/python/test_scenarios.py`, notebook §15):
> - Figure-eight returns to its initial state to 3e-8 after one period (yoshida4, 8000 steps).
> - Plummer sphere (N = 1000): time-averaged virial ratio 2K/|U| = 0.992; energy within sampling noise of -3π/64.
> - Pendulum at 2 rad: first crossing at 3T/4 of the elliptic-integral period to 1e-8.
> - Spring chain: the lowest mode returns after 2π/ω with ω = 2√(k/m) sin(π/2(n-1)) to 1e-6.
> - Two-body: energy -G m1 m2 / 2a to 1e-12, return after one period; solar system: Earth 0.983 AU from the Sun at
>   J2000, Jupiter's period 11.86 yr, outer planets conserve energy to 3e-6 over 1000 years; SI and astronomical runs
>   agree to 1e-9 AU.

## Why
Setting up standard systems by hand is repetitive and error-prone, and the engine has no notion of units.

## Scope
- `physim.scenarios`: two-body orbit from orbital elements, solar system from tabulated elements,
  Plummer sphere, figure-eight three-body orbit, pendulum, spring lattice.
- Unit systems: SI, astronomical (AU, year, solar mass, with G chosen so G = 4π²), and dimensionless.
  Keep the engine unit-agnostic; do conversions in Python.

## Done when
- Each scenario has a test checking a known property (the figure-eight returns to its start after one
  period; Plummer sphere satisfies the virial theorem on average).
