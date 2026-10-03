# 07 · Holonomic constraints (SHAKE/RATTLE)

**Priority:** P1 · **Size:** M · **Area:** Physics

## Why
Rigid rods, pendulums and fixed bond lengths are currently faked with stiff springs, which forces a
tiny time step and adds spurious high-frequency motion.

## Scope
- Distance constraints between two particles, and between a particle and a fixed point.
- RATTLE (constrained velocity Verlet): solve for Lagrange multipliers iteratively after the drift
  (positions) and after the final kick (velocities), to a set tolerance.
- Report constraint forces so they can be inspected.

## Notes
- RATTLE stays symplectic and second order; composing it into `yoshida4` needs care, so start with Verlet only.
- Remove the constrained degrees of freedom from any temperature or equipartition calculations later (item 12).

## Done when
- A rigid pendulum's small-amplitude period matches 2π√(L/g), and large-amplitude period matches the
  elliptic-integral result.
- Constraint lengths hold to the set tolerance over 10⁶ steps; energy error stays bounded.
