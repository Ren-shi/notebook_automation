# 35 · Two-body reaction kinematics

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner

## Why
Every output of the planner starts from kinematics: which energy a scattered particle has at a given angle, where
the recoil goes, how lab and centre-of-mass quantities convert. Writing it for a general reaction a(b, c)d with a
Q-value now means inelastic scattering (item 43) and transfer reactions only add a term, not a rewrite.

## Scope
- Fully relativistic two-body kinematics for a(b, c)d with Q = (m_a + m_b − m_c − m_d)c² − E*; elastic scattering is
  the case c = b, d = a, E* = 0.
- Lab energy and angle of the ejectile and the recoil, as functions of lab or CM angle.
- Both solutions where the kinematics are double-valued (beam heavier than the target, or endothermic reactions near
  threshold), and the maximum lab angle.
- Lab ↔ CM angle conversion and the solid-angle Jacobian dΩ_cm/dΩ_lab (needed by item 37).
- Kinematic broadening dE/dθ, used to estimate peak widths from a detector's angular size.
- Threshold energy for reactions with Q < 0.
- Prototype in Python, then port to Rust (the `docs/new-ideas.md` workflow) once validated, because item 39 calls it
  per event.

## Depends on
34.

## Done when
- The non-relativistic limit matches the analytic elastic formulas (kinematic factor K) to 1e-6 at low energy.
- Energies and angles match LISE++ (or CATKIN) for at least five cases, including inverse kinematics and a
  double-valued case, within 1 keV and 0.01°. Reference files stored as described in item 40.
- Jacobian checked against a numerical derivative of the angle conversion.
- Theory note in `docs/theory/`, and an entry in the physics register.
