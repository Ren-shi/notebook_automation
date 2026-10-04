# Physics register

What physics Physim can simulate today, how far each capability has got, and the evidence that it is right. Written
for anyone auditing the project: students deciding whether to trust a number, supervisors, reviewers.

Physim brings together calculations usually done with LISE++, SRIM and similar tools. It does not replace them.
Every capability is checked against those tools first, then against published measurements, and this register
records where each one stands.

## Status levels

| Column | Meaning |
|---|---|
| **Engine** | Implemented and usable from Python. |
| **App** | Available in the experiment planner web app, without writing code. |
| **vs tools** | Agrees with SRIM, LISE++ or another established code within the stated tolerance. |
| **vs literature** | Agrees with published data or analytic results within the stated tolerance. |

✅ done · 🟡 in progress · ⬜ planned · — not applicable

## Rules

1. **Every ✅ in a validation column points to an automated test** that checks the stated tolerance and runs in CI.
   If the test fails, the entry is no longer ✅. A claim without a test is 🟡 at best.
2. **Tolerances are stated, not implied.** "Agrees with SRIM" is not an entry; "agrees with SRIM within 3% for
   α in Au from 1 to 50 MeV" is.
3. **Limits are listed as clearly as capabilities.** Each page says where its model stops being valid.
4. **Reference data has provenance**: tool and version (or paper and table), exact inputs, date (see backlog item 40).

## Nuclear experiment planner

First slice: elastic scattering (backlog items 33–42).

| Capability | Engine | App | vs tools | vs literature | Page | Backlog |
|---|---|---|---|---|---|---|
| Experiment setup file | ✅ | ⬜ | — | — | [Guide](../nuclear-setup.md) | 33 |
| Atomic masses and material data | ✅ | ⬜ | — | ✅ Q-values within 1 keV; atomic weights vs NIST | [Data](nuclear-data.md) | 34 |
| Two-body reaction kinematics | ⬜ | ⬜ | ⬜ LISE++ | ⬜ analytic limit | — | 35 |
| Stopping power and range | ⬜ | ⬜ | ⬜ SRIM | ⬜ NIST PSTAR/ASTAR | — | 36 |
| Energy and angular straggling | ⬜ | ⬜ | ⬜ TRIM | ⬜ | — | 36 |
| Rutherford cross section | ⬜ | ⬜ | — | ⬜ analytic, Geiger–Marsden | — | 37 |
| Distance of closest approach and validity checks | ⬜ | ⬜ | ⬜ LISE++ | ⬜ | — | 37 |
| Coulomb trajectories | ⬜ | ⬜ | — | ⬜ analytic b(θ) | — | 37 |
| Detector solid angles and response | ⬜ | ⬜ | — | ⬜ analytic | — | 38 |
| Count rates and beam time | ⬜ | ⬜ | ⬜ LISE++ | ⬜ | — | 39 |
| Monte Carlo spectra | ⬜ | ⬜ | ⬜ | ⬜ | — | 39 |
| Beam-time report, CSV and ROOT export | ⬜ | ⬜ | — | — | — | 42 |

Next slices: inelastic scattering and Coulomb excitation (43), then transfer reactions and fusion-evaporation.

## Foundations: the classical and quantum engine

The general-purpose engine underneath the planner. Its validation is against analytic solutions and conservation
laws rather than nuclear-physics tools; the results are listed in the "Done" note at the top of each backlog item,
and the tests are in `tests/` (Rust) and `tests/python/`. None of it is in the planner app yet.

| Capability | Engine | Validated against | Backlog |
|---|---|---|---|
| Integrators (symplectic, high-order, adaptive) | ✅ | Convergence orders, energy conservation | 08, 17 |
| Gravity (direct and Barnes–Hut), springs, built-in forces | ✅ | Tree error scaling vs direct sum; resonance curve, 1PN periapsis advance | 10, 19 |
| Charged particles: Lorentz force, Boris pusher, Coulomb | ✅ | Analytic gyration, E×B drift, Rutherford deflection angle | 09 |
| Constraints, collisions, rigid bodies | ✅ | Pendulum periods, textbook collision velocities, Euler-top flip period | 07, 11, 13 |
| Molecular dynamics | ✅ | Thermostat statistics, conservation | 12 |
| Chaos indicators | ✅ | Kepler λ → 0, symplectic ± pairing of Lyapunov spectra | 18 |
| Fields on grids: waves, diffusion, Poisson | ✅ | Analytic solutions, convergence orders | 20 |
| Single-particle quantum mechanics | ✅ | Analytic levels, transmission, norm conservation | 32 |

## Adding an entry

When a capability reaches **Engine**, give it a page in this folder (`kinematics.md`, `stopping-power.md`, …) with
these sections, and link it from the table:

- **Model**: the equations used, and any data (with version).
- **Assumptions and range of validity**: where the model holds and where it does not.
- **Validation**: for each comparison, the reference (tool and version, or paper), the cases compared, the
  tolerance, a plot or a link to the validation notebook, the test that enforces it, and the date and commit.
- **Known limitations**: what is missing or approximate.
- **References**: papers and data sources to cite.
