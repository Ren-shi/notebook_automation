# 18 · Chaos and stability analysis

**Priority:** P2 · **Size:** M · **Area:** Analysis · **Status: Done**

> **Done** (`src/chaos.rs`, `Force::jacobian_vector`/`ForceSet::jacobian_vector`, `HenonHeiles` force, Python
> `w.lyapunov(dt, steps, n)` and `ps.poincare_section`). Tangent vectors are carried as blocks of pseudo-particles
> whose acceleration is the Jacobian-vector product, so the world's own integrator advances state and tangents
> together (for Verlet-type schemes this is the exact linearisation of the numerical map, and with explicit
> integrators the real particles follow bit-for-bit the same trajectory as a plain `run`; the implicit Gauss methods'
> convergence test also sees the tangents, so their trajectory can differ at round-off). Analytic Jacobian-vector products for gravity, springs,
> anchor springs, drag, uniform fields, all central potentials, the harmonic trap, Hénon-Heiles and periodic forcing;
> central finite differences (relative error ~1e-10) for every other force, including `CustomForce`. Full spectrum by
> modified Gram-Schmidt (QR); MEGNO from the leading tangent vector's growth each step.
>
> Results (`tests/chaos.rs`, `tests/python/test_chaos.py`, notebook §10):
> - Hénon-Heiles: largest exponent 0.10-0.13 for chaotic orbits at E = 1/6 and 0.06 at E = 1/8; regular orbits
>   0.003 at t = 2000 and falling like ln t / t; MEGNO ⟨Y⟩ = 2.00 ± 0.02 on regular orbits, 50-120 (growing) on chaotic
>   ones. Every integrator tested (verlet, yoshida4, rk4, gauss4, pefrl, dopri5) classifies the orbits the same way.
> - Kepler: λ_max = 0.004 at t = 2000 (→ 0) and ⟨Y⟩ = 1.995.
> - The full six-exponent spectrum comes in ± pairs (0.1251 / −0.1256) summing to 2e-16, as a symplectic flow requires.
> - Analytic Jacobians agree with finite differences to 1e-7 for every implementation.
> - The Hénon-Heiles Poincaré section at E = 1/8 shows the familiar islands and chaotic sea (notebook §10), with
>   crossings located to 1e-12.
>
> Not done: fast Lyapunov indicators and frequency analysis, Lyapunov exponents with adaptive steps, and
> Jacobian-vector products for constraints. Parallel blocking for the gravity Jacobian (currently serial O(N²)).

## Why
Testing dynamical ideas often means asking whether motion is regular or chaotic, and how fast nearby
trajectories separate. None of that is available yet.

## Scope
- Variational equations: integrate the tangent map alongside the state (needs the force Jacobian; provide
  analytic Jacobians for built-in forces and finite differences as a fallback).
- Lyapunov exponents: the largest via renormalised tangent vectors, the full spectrum via QR (Benettin).
- MEGNO as a faster chaos indicator.
- Poincaré sections, using event detection (item 05) to find crossings of the section exactly.

## Depends on
05 (for Poincaré sections).

## Done when
- The largest Lyapunov exponent is zero (to tolerance) for the Kepler problem and positive for the
  Hénon–Heiles system at E = 1/6.
- A Hénon–Heiles Poincaré section reproduces the familiar mix of islands and chaotic sea.
