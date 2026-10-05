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
| 14 | ~~[Scenario library and units](14-scenarios-and-units.md)~~ **Done** | Usability | P3 | M | — |
| 15 | ~~[Visualization and animation helpers](15-visualization.md)~~ **Done** | Usability | P3 | S | — |
| 16 | ~~[World API gaps: particle and force management](16-world-api-gaps.md)~~ **Done** | Core | P1 | S | — |
| 17 | ~~[Higher-order and specialised integrators](17-more-integrators.md)~~ **Done** | Integrators | P2 | M | — |
| 18 | ~~[Chaos and stability analysis](18-chaos-and-stability.md)~~ **Done** | Analysis | P2 | M | 05 |
| 19 | ~~[More built-in forces](19-more-forces.md)~~ **Done** | Physics | P2 | M | — |
| 20 | ~~[Fields and continua](20-fields-and-continua.md)~~ **Done** | Physics | P3 | L | — |
| 21 | ~~[Packaging and release](21-packaging-and-release.md)~~ **Done** | Tooling | P3 | S | 01 |
| 22 | ~~[Documentation](22-documentation.md)~~ **Done** | Usability | P3 | M | — |
| 23 | [Equations of motion from a Lagrangian or Hamiltonian](23-lagrangian-input.md) | Usability | P2 | L | 17 |
| 24 | [Ensembles and parameter sweeps](24-ensembles-and-sweeps.md) | Analysis | P2 | M | — |
| 25 | [Periodic orbits: finding and continuation](25-periodic-orbits.md) | Analysis | P2 | M | 18 |
| 26 | [Conserved-quantity monitors](26-conserved-quantity-monitors.md) | Analysis | P2 | S | — |
| 27 | [Rotating reference frames](27-rotating-frames.md) | Physics | P3 | S | — |
| 28 | [Richer constraints](28-richer-constraints.md) | Physics | P3 | L | 07 |
| 29 | [Variable mass](29-variable-mass.md) | Physics | P3 | S | — |
| 30 | [Special-relativistic particle dynamics](30-special-relativity.md) | Physics | P3 | M | 09 |
| 31 | ~~[Seeded, reproducible randomness](31-seeded-randomness.md)~~ **Done** | Core | P2 | S | — |
| 32 | ~~[Quantum mechanics: the Schrödinger equation on a grid](32-quantum-mechanics.md)~~ **Done** | Physics | P3 | L | 20 |
| 33 | ~~[Experiment definition: the setup file](33-experiment-definition.md)~~ **Done** | Nuclear planner | P1 | M | — |
| 34 | ~~[Nuclear and material data](34-nuclear-and-material-data.md)~~ **Done** | Nuclear planner | P1 | M | — |
| 35 | ~~[Two-body reaction kinematics](35-reaction-kinematics.md)~~ **Done** | Nuclear planner | P1 | M | 34 |
| 36 | ~~[Stopping power and energy loss](36-stopping-power-and-energy-loss.md)~~ **Done** | Nuclear planner | P1 | L | 34 |
| 37 | ~~[Rutherford scattering, closest approach and Coulomb trajectories](37-rutherford-scattering.md)~~ **Done** | Nuclear planner | P1 | M | 34, 35 |
| 38 | ~~[Detector geometry and response](38-detector-geometry-and-response.md)~~ **Done** | Nuclear planner | P1 | M | 33, 36 |
| 39 | ~~[Count rates, beam time and the Monte Carlo event generator](39-rates-and-event-generator.md)~~ **Done** | Nuclear planner | P1 | L | 31, 35–38 |
| 40 | ~~[Validation suite against SRIM, LISE++ and literature](40-validation-suite.md)~~ **Done** (user-run references pending) | Nuclear planner | P1 | M | 33 |
| 41 | ~~[Experiment planner web app](41-planner-web-app.md)~~ **Done** (student trial pending) | Nuclear planner | P1 | L | 33, 35–39 |
| 42 | ~~[Beam-time report and data exports](42-report-and-exports.md)~~ **Done** (ROOT-side check pending) | Nuclear planner | P1 | M | 33, 39 |
| 43 | ~~[Inelastic scattering and Coulomb excitation (second slice)](43-inelastic-scattering.md)~~ **Done** (GOSIA comparison pending) | Nuclear planner | P2 | L | 33–42 |
| 44 | [One-click installer for the planner app](44-one-click-installer.md) | Nuclear planner | P2 | M | 41 |
| 45 | ~~[Coulomb excitation set up entirely in the app](45-coulex-in-the-app.md)~~ **Done** | Nuclear planner | P1 | S | 41, 43 |
| 46 | ~~[Guided workflow, field help and "how to read this"](46-guided-workflow.md)~~ **Done** (new-user trial pending) | Nuclear planner | P1 | M | 41, 45 |
| 47 | ~~[Backward detectors and the Coulomb-excitation beam time](47-backward-detectors-and-coulex-beam-time.md)~~ **Done** | Nuclear planner | P1 | S | 45, 46 |
| 48 | ~~[Hidden detectors counted in the analytic rates](48-hidden-detectors-in-the-rates.md)~~ **Done** | Nuclear planner | P1 | S | 38, 39 |
| 49 | ~~[Light and dark themes, and paper figures in a journal's style](49-themes-and-paper-figures.md)~~ **Done** | Nuclear planner | P1 | M | 41, 42 |
| 50 | ~~[Level schemes from ENSDF](50-level-schemes-from-ensdf.md)~~ **Done** (NNDC permission pending) | Experiment workbench | P1 | L | 34 |
| 51 | ~~[Detector catalogue: S3, clover and LaBr₃ as solids](51-detector-catalogue.md)~~ **Done** (typical values to confirm) | Experiment workbench | P1 | M | 38, 48 |
| 52 | [The interactive scene](52-interactive-scene.md) | Experiment workbench | P1 | L | 41, 51 |
| 53 | [γ-detector response and calibration sources](53-gamma-detector-response.md) | Experiment workbench | P1 | L | 51 |
| 54 | [Orientation of excited states, particle–γ correlation and decay](54-angular-correlation-and-decay.md) | Experiment workbench | P1 | L | 43, 50 |
| 55 | [Particle–γ events, Doppler correction, coincidences and background](55-coincidence-events-and-background.md) | Experiment workbench | P1 | L | 39, 53, 54 |
| 56 | [Automatic analysis: from peak areas to B(E2) and the shape](56-automatic-analysis.md) | Experiment workbench | P1 | L | 55 |
| 57 | [The run record: every number explains itself](57-run-record.md) | Experiment workbench | P1 | M | 42, 56 |
| 58 | [Multi-step Coulomb excitation, reorientation and the shape of the nucleus](58-multi-step-coulomb-excitation.md) | Experiment workbench | P1 | L | 50, 54, 56 |
| 59 | [Particle and γ-ray tracks in the scene](59-particle-tracks.md) | Experiment workbench | P2 | M | 52, 55 |
| 60 | [A movable target and sources placed anywhere](60-movable-target-and-sources.md) | Experiment workbench | P2 | M | 52, 53 |
| 61 | [The planner in its own desktop window](61-desktop-window.md) | Experiment workbench | P2 | S | 44, 52 |
| 62 | [Add-back and Compton suppression for clovers](62-add-back-and-compton-suppression.md) | Experiment workbench | P3 | M | 53, 55 |

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

Items 23–32 were added after the original plan and are deferred: finish 08–22 first.

## Nuclear experiment planner (items 33–43)

The next direction: a browser-based tool for planning nuclear physics experiments, aimed first at graduate students
(who often do not want to write code) and later at undergraduates learning the physics. The user describes a beam,
a target and a set of detectors; Physim returns kinematics, cross sections, energy loss, count rates, beam time and
simulated spectra, and writes a **beam-time report** (the main deliverable) with CSV and ROOT exports.

Physim brings together calculations that students currently spread across LISE++, SRIM and kinematics calculators.
It does not replace those tools: each capability is checked against them first, then against published data, and the
result is recorded in the [physics register](../docs/physics-register/README.md). Exports to SRIM and LISE++ formats
let users cross-check any result in the trusted tool.

Slices, one at a time:
1. **Elastic scattering** (33–42), the first vertical slice. Suggested order: 33 and 34 → 35 → 36 and 37 → 38 →
   31 and 39 → 40 throughout → 41 and 42 as the physics lands.
2. Inelastic scattering and Coulomb excitation (43).
3. Transfer reactions, then fusion-evaporation (not yet written up).

Items 33–42 take priority over the deferred items 23–30. Item 31 (seeded randomness) is needed by 39 and moves up
with them.

## Experiment workbench (items 50–62)

The direction after the planner: a simulated experiment that the user sees and rearranges, aimed first at an
experimentalist planning beam time. Detectors are drawn to scale and can be dragged within what a chamber allows;
level schemes come from ENSDF; particle–γ events are analysed to a Doppler-corrected spectrum, a B(E2) and the shape
of the nucleus; and every number explains how it was obtained.

The aim, every decision and what is parked are in the [overview](workbench-overview.md). Suggested order: 50 and
51 → 52 → 53 and 54 → 55 → 56 → 57 → 58 → 59–62.

To add an item: copy any file, give it the next number, and add a row to the table.
