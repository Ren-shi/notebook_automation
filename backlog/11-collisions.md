# 11 · Collisions and contact

**Priority:** P2 · **Size:** L · **Area:** Physics · **Status: Done**

> **Done** (`src/collisions.rs`, `src/broadphase.rs`, `src/forces/contact.rs`; radii in `State`). Per-particle radii
> (`World::set_radii`/`set_radius`, Python `add_particle(..., radius=r)` and `w.radii`) reach forces through a new
> default method `Force::accumulate_full` (radii are threaded through `ForceSet` like charges; the acceleration cache keys
> on them). Broad phase: uniform grid with cell size 2 × max reach (O(N)). Two response models:
> - **Hard, event-driven** (`World::set_collisions(Collisions { restitution, walls, between_particles, max_per_step })`,
>   Python `w.set_collisions(...)`): after each step the earliest contact is found from swept spheres and straight-path
>   contact times (exact without forces; refined by Illinois root finding on the integrator's own trajectory when forces
>   act), the world is re-stepped to that instant, the impulse with restitution `e` is applied along the line of
>   centres (or the wall normal; pinned particles are immovable), and the rest of the step continues. Planar walls and
>   `Wall::box_walls`/`ps.box_walls`. Collision counts are reported; settings are saved in checkpoints (restarts
>   bit-exact).
> - **Soft** `SoftContact { k, damping, law: Linear | Hertz }`: an ordinary force on overlapping spheres, never attractive.
>
> Results (`tests/collisions.rs`, `tests/python/test_collisions.py`, notebook §12):
> - Head-on collisions of unequal masses give the textbook velocities to 1e-12 at any step size (dt = 0.37 or 0.01), with
>   momentum and energy conserved to 1e-12; equal-mass oblique elastic collisions leave at 90°; restitution `e` removes
>   exactly (1 − e²) μ u_n²/2.
> - A ball bouncing under gravity: apex-height ratios e² to ~1e-6 (the apex sampling error), and it never enters the floor.
> - Hard-sphere gas, 500 spheres in a box, all starting at the same speed: 9 400 collisions, energy conserved to 3e-15,
>   ⟨v⁴⟩/⟨v²⟩² goes from 1 to 1.669 (Maxwell-Boltzmann: 5/3), every sphere stays inside.
> - Soft contact: linear-spring contact lasts π√(μ/k) (to 2e-4 at dt = 1e-5) and rebounds elastically; Hertz contact
>   duration 2.943 δ_max/v to 5.7e-5; a linear dashpot gives e = 0.8085 against exp(−πζ/√(1−ζ²)) = 0.8004 (the
>   never-attractive clamp ends contact slightly early, a known effect).
>
> Limitations: events that fall inside a step with a collision are located on the collision-free trajectory; adaptive
> runs and Lyapunov exponents refuse hard collisions; a resting pile keeps colliding every step (use soft contact or
> restitution < 1 with a cap). Not done: friction and rotation (needs item 13), sweep-and-prune, soft walls.

## Why
Granular systems, gases of hard spheres and impact problems need particles with size that collide.

## Scope
- Particle radius in `State`.
- Broad phase: uniform grid or sweep-and-prune to avoid O(N²) checks.
- Two response models:
  - Impulse-based with coefficient of restitution (event-driven, uses item 05 to find contact time).
  - Soft contact (Hertzian or linear spring-dashpot) as an ordinary force, simpler and parallel-friendly.
- Walls and boxes as static boundaries.

## Done when
- Elastic collisions conserve energy and momentum to tolerance; restitution e gives the expected energy loss.
- A hard-sphere gas relaxes to a Maxwell–Boltzmann speed distribution.
