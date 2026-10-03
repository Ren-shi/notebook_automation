# 11 · Collisions and contact

**Priority:** P2 · **Size:** L · **Area:** Physics

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
