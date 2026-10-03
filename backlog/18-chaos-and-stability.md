# 18 · Chaos and stability analysis

**Priority:** P2 · **Size:** M · **Area:** Analysis

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
