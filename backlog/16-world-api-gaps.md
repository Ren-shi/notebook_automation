# 16 · World API gaps: particle and force management

**Priority:** P1 · **Size:** S · **Area:** Core · **Status: Done**

> **Done.** Implemented as described below, with one design choice: removing a particle that a force still
> uses is refused (rather than deleting the force), and higher indices are renumbered. Force parameters are
> changed with `set_force_params(id, name=value)`, all-or-nothing. Pinned particles have zero velocity and
> acceleration, so every integrator keeps them fixed exactly. A failed step rolls the world back to the last
> completed step. Tests: `tests/world_api.rs`, `tests/python/test_world_api.py`.

## Why
Some basic operations are missing, which gets in the way of exploratory work:
- Particles can be added but not removed, and masses cannot be changed after creation (`World.masses` is read-only).
- Forces can only be added or cleared all at once; there is no way to remove one force, read it back, or change a
  parameter (e.g. sweep `G` or a spring constant) without rebuilding the world.
- Every particle needs a positive mass, so massless test particles (tracers in a gravitational field) are impossible,
  and there is no way to pin a particle in place.
- If a run fails partway (a Python exception, or the state becomes non-finite), `World::run` returns the error
  and the frames recorded so far are lost.
- Integrator order and whether it is symplectic exist in Rust but are not visible from Python.

## Scope
- `remove_particle(i)` with clear rules for forces that reference particle indices (`Spring`, `AnchorSpring`):
  reject, or renumber them.
- `masses` setter.
- Forces get handles: `add_force` returns an id; `remove_force(id)`, `get_force(id)`, and mutable parameters.
- Massless test particles: feel forces but do not source gravity. Pinned particles: never move.
- On failure, return the partial trajectory with the error (e.g. attach it to the Python exception).
- `World.integrator_info` → name, order, symplectic.

## Done when
- A parameter sweep over `G` reuses one `World`.
- Tracer particles orbit a mass without perturbing it (the central body's momentum stays exactly zero).
- A run that fails at step k gives back frames 0..k.
