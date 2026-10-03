# Backlog

Work not yet built. One file per item: why it matters, scope, design notes and what "done" means.
Priority: **P1** = next up / unblocks other work, **P2** = important capability, **P3** = nice to have.
Size: **S** ≈ under a day, **M** ≈ a few days, **L** ≈ a week or more.

| # | Item | Area | Priority | Size | Depends on |
|---|---|---|---|---|---|
| 01 | ~~[Continuous integration](01-continuous-integration.md)~~ **Done** | Tooling | P1 | S | — |
| 02 | ~~[Benchmarks](02-benchmarks.md)~~ **Done** | Tooling | P1 | S | — |
| 03 | ~~[Reuse forces between Verlet steps](03-verlet-force-caching.md)~~ **Done** | Performance | P1 | S | 02 |
| 04 | ~~[Parallel force evaluation](04-parallel-forces.md)~~ **Done** | Performance | P1 | M | 02 |
| 05 | ~~[Event detection and stopping conditions](05-event-detection.md)~~ **Done** | Core | P1 | M | — |
| 06 | ~~[Save, load and stream simulation data](06-save-load-and-streaming.md)~~ **Done** | Core | P1 | M | — |
| 07 | ~~[Holonomic constraints (SHAKE/RATTLE)](07-constraints.md)~~ **Done** | Physics | P1 | M | — |
| 08 | ~~[Adaptive time stepping](08-adaptive-time-stepping.md)~~ **Done** | Integrators | P2 | M | 05 |
| 09 | ~~[Charged particles: Lorentz force and Boris pusher](09-lorentz-force-boris.md)~~ **Done** | Physics | P2 | M | — |
| 10 | ~~[Tree gravity (Barnes–Hut)](10-barnes-hut.md)~~ **Done** | Performance | P2 | L | 02, 04 |
| 11 | ~~[Collisions and contact](11-collisions.md)~~ **Done** | Physics | P2 | L | 05 |
| 12 | ~~[Molecular dynamics: pair potentials, periodic boxes, thermostats](12-molecular-dynamics.md)~~ **Done** | Physics | P2 | L | 04 |
| 13 | ~~[Rigid bodies](13-rigid-bodies.md)~~ **Done** | Physics | P3 | L | 07 |
| 14 | [Scenario library and units](14-scenarios-and-units.md) | Usability | P3 | M | — |
| 15 | [Visualization and animation helpers](15-visualization.md) | Usability | P3 | S | — |
| 16 | ~~[World API gaps: particle and force management](16-world-api-gaps.md)~~ **Done** | Core | P1 | S | — |
| 17 | ~~[Higher-order and specialised integrators](17-more-integrators.md)~~ **Done** | Integrators | P2 | M | — |
| 18 | ~~[Chaos and stability analysis](18-chaos-and-stability.md)~~ **Done** | Analysis | P2 | M | 05 |
| 19 | ~~[More built-in forces](19-more-forces.md)~~ **Done** | Physics | P2 | M | — |
| 20 | [Fields and continua](20-fields-and-continua.md) | Physics | P3 | L | — |
| 21 | [Packaging and release](21-packaging-and-release.md) | Tooling | P3 | S | 01 |
| 22 | [Documentation](22-documentation.md) | Usability | P3 | M | — |
| 23 | [Equations of motion from a Lagrangian or Hamiltonian](23-lagrangian-input.md) | Usability | P2 | L | 17 |
| 24 | [Ensembles and parameter sweeps](24-ensembles-and-sweeps.md) | Analysis | P2 | M | — |
| 25 | [Periodic orbits: finding and continuation](25-periodic-orbits.md) | Analysis | P2 | M | 18 |
| 26 | [Conserved-quantity monitors](26-conserved-quantity-monitors.md) | Analysis | P2 | S | — |
| 27 | [Rotating reference frames](27-rotating-frames.md) | Physics | P3 | S | — |
| 28 | [Richer constraints](28-richer-constraints.md) | Physics | P3 | L | 07 |
| 29 | [Variable mass](29-variable-mass.md) | Physics | P3 | S | — |
| 30 | [Special-relativistic particle dynamics](30-special-relativity.md) | Physics | P3 | M | 09 |
| 31 | [Seeded, reproducible randomness](31-seeded-randomness.md) | Core | P2 | S | — |

Known limits of the current engine (worth keeping in mind until the items above land):
- `NewtonianGravity` is direct O(N²); `TreeGravity` (item 10) is O(N log N) but its per-interaction cost is not yet
  SIMD-vectorised, which caps its speedup at ~N / (2 × interactions per particle).
- `verlet`/`yoshida4` are only symplectic for velocity-independent forces; use `rk4` with drag.
- Rods (item 07) need `verlet` or `yoshida4`; long chains use an iterative solver (~43 ms/step at 1 000 links).
- `CustomForce` costs one Python call per force evaluation.
- Every step copies the whole state so a failed step can be rolled back: ~20 ns per particle on the CI-class machines
  used here, which dominates for cheap forces (a 10 000-particle spring chain spends ~200 of ~235 µs per Verlet step
  outside the force).

Finished items stay in the table, struck through and marked **Done**, with a note at the top of their file.

Items 23–31 were added after the original plan and are deferred: finish 08–22 first.

To add an item: copy any file, give it the next number, and add a row to the table.
