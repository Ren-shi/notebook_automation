# Backlog

Work not yet built. One file per item: why it matters, scope, design notes and what "done" means.
Priority: **P1** = next up / unblocks other work, **P2** = important capability, **P3** = nice to have.
Size: **S** ≈ under a day, **M** ≈ a few days, **L** ≈ a week or more.

| # | Item | Area | Priority | Size | Depends on |
|---|---|---|---|---|---|
| 01 | [Continuous integration](01-continuous-integration.md) | Tooling | P1 | S | — |
| 02 | [Benchmarks](02-benchmarks.md) | Tooling | P1 | S | — |
| 03 | [Reuse forces between Verlet steps](03-verlet-force-caching.md) | Performance | P1 | S | 02 |
| 04 | [Parallel force evaluation](04-parallel-forces.md) | Performance | P1 | M | 02 |
| 05 | [Event detection and stopping conditions](05-event-detection.md) | Core | P1 | M | — |
| 06 | [Save, load and stream simulation data](06-save-load-and-streaming.md) | Core | P1 | M | — |
| 07 | [Holonomic constraints (SHAKE/RATTLE)](07-constraints.md) | Physics | P1 | M | — |
| 08 | [Adaptive time stepping](08-adaptive-time-stepping.md) | Integrators | P2 | M | 05 |
| 09 | [Charged particles: Lorentz force and Boris pusher](09-lorentz-force-boris.md) | Physics | P2 | M | — |
| 10 | [Tree gravity (Barnes–Hut)](10-barnes-hut.md) | Performance | P2 | L | 02, 04 |
| 11 | [Collisions and contact](11-collisions.md) | Physics | P2 | L | 05 |
| 12 | [Molecular dynamics: pair potentials, periodic boxes, thermostats](12-molecular-dynamics.md) | Physics | P2 | L | 04 |
| 13 | [Rigid bodies](13-rigid-bodies.md) | Physics | P3 | L | 07 |
| 14 | [Scenario library and units](14-scenarios-and-units.md) | Usability | P3 | M | — |
| 15 | [Visualization and animation helpers](15-visualization.md) | Usability | P3 | S | — |

Known limits of the current engine (worth keeping in mind until the items above land):
- Gravity is direct O(N²) summation, so practical up to a few thousand bodies.
- `verlet`/`yoshida4` are only symplectic for velocity-independent forces; use `rk4` with drag.
- No constraints: a pendulum is a stiff `AnchorSpring`, which forces a small time step.
- `World::run` keeps every recorded frame in memory.
- `CustomForce` costs one Python call per force evaluation.

To add an item: copy any file, give it the next number, and add a row to the table.
