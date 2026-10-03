# 07 · Holonomic constraints (SHAKE/RATTLE)

**Priority:** P1 · **Size:** M · **Area:** Physics · **Status: Done**

> **Done** (`src/constraints.rs`, `World::add_rod`, Python `w.add_rod(i, to, length=None)`, `traj.tension`).
> Rods between particles or to fixed points, integrated with RATTLE in `verlet` and composed into `yoshida4` (fourth
> order, still symplectic; the composition was straightforward, so both are supported). Tensions are reported per
> rod after every step and recorded in trajectories. Constraints are included in checkpoints and run metadata.
>
> Two changes from the plan, both found by testing:
> - **Solver.** Rod-by-rod (Gauss-Seidel) SHAKE took 94 µs/step for a 10-link chain and did not converge at all for
>   30+ links released straight (its convergence factor is ~1 − π²/n² for a collinear chain). All rods are now solved
>   together: conjugate gradients on `G M⁻¹ Gᵀ` (matrix-free, Jacobi preconditioned), with a chord-Newton outer loop
>   for positions (M-SHAKE). 10 links: 6 µs/step; 30: 34 µs; 100: 0.36 ms; 300: 3.3 ms; 1 000: 43 ms, lengths
>   held to ~1e-13 in all cases.
> - **Tension accuracy.** A position residual δ left by the solver turns into a tension error of order δ/h², which
>   showed up as 4.5e-6 under yoshida4 at the default tolerance. One extra correction after convergence squares the
>   residual; a hanging bob's tension is now m g to 6e-12.
>
> Results (`tests/constraints.rs`, `tests/python/test_constraints.py`, notebook §6):
> - Small-amplitude period matches 2π√(L/g) (to the θ₀²/16 correction) and the elliptic result to 1e-10;
>   large-amplitude periods (θ₀ = 1, 2, 3; and 15 angles up to 3.0 in the notebook) match 4√(L/g) K(sin(θ₀/2))
>   to 1.2e-10 with yoshida4 and to 1e-6 with Verlet at dt = 2e-4.
> - A chaotic double pendulum over 10⁶ steps: every step within the 1e-10 length tolerance; worst |ΔE/E| per
>   10⁵-step block stays flat at ~2.4e-4 (Verlet, dt = 1e-3) and ~1.2e-7 (yoshida4): bounded, no drift.
> - Measured convergence orders: 2.0 (Verlet) and 4.0 (yoshida4) on a large-amplitude pendulum.
> - Tensions: hanging m g to 6e-12; at the bottom of a swing from horizontal 3 m g to 2e-5 (O(dt²));
>   a free spinning dumbbell gives the centripetal tension and conserves momentum and angular momentum to 1e-12.
>
> Follow-up: chains cost O(links²) per step. The constraint graph of a chain or tree gives a matrix with a
> zero-fill elimination order (leaf rods first), so a direct solve would be O(links); worth doing if long ropes or
> polymers become a use case (item 12).

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
