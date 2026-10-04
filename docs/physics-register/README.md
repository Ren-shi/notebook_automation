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
5. **The table is checked by a test.** `physim.nuclear.validation` holds every comparison behind this table, and
   `tests/python/test_nuclear_validation.py` fails if a ✅ here is not backed by passing checks, or if a deliberately
   broken formula slips through. `python -c "from physim.nuclear import validation; print(validation.report())"`
   prints the current state, and the {doc}`validation notebooks <../validation>` plot every comparison.

## Nuclear experiment planner

First slice: elastic scattering (backlog items 33–42).

| Capability | Engine | App | vs tools | vs literature | Page | Backlog |
|---|---|---|---|---|---|---|
| Experiment setup file | ✅ | ⬜ | — | — | [Guide](../nuclear-setup.md) | 33 |
| Atomic masses and material data | ✅ | ⬜ | ✅ LISE++ masses, 46 nuclides, within 2 keV or the AME2020 uncertainty | ✅ Q-values within 1 keV; atomic weights vs NIST | [Data](nuclear-data.md) | 34 |
| Two-body reaction kinematics | ✅ | ⬜ | 🟡 LISE++: 1 case passes, 8 requested | ✅ analytic limit, Lorentz boost, threshold | [Kinematics](kinematics.md) | 35 |
| Stopping power and range | ✅ | ⬜ | ✅ LISE++ (1680 points, tolerances by regime); 🟡 SRIM requested | ✅ NIST CSDA ranges within 1% | [Stopping](stopping.md) | 36 |
| Energy and angular straggling | ✅ | ⬜ | ✅ LISE++ ≥ 30 MeV/u (25%); 🟡 TRIM requested | ✅ Bohr and Highland formulas | [Stopping](stopping.md) | 36 |
| Rutherford cross section | ✅ | ⬜ | ✅ LISE++ (1 case, 0.03%) | ✅ analytic, Geiger–Marsden 1913 | [Rutherford](rutherford.md) | 37 |
| Distance of closest approach and validity checks | ✅ | ⬜ | ✅ LISE++ grazing angle (1 case) | 🟡 above-barrier data requested | [Rutherford](rutherford.md) | 37 |
| Coulomb trajectories | ✅ | ⬜ | — exact closed form instead | ✅ analytic b(θ) to 1e-6° | [Rutherford](rutherford.md) | 37 |
| Detector solid angles and response | ✅ | ⬜ | — exact closed forms instead | ✅ closed forms to 1e-10, Monte Carlo | [Detectors](detectors.md) | 38 |
| Count rates and beam time | ✅ | ⬜ | 🟡 LISE++ cross sections at the detector angles requested | ✅ I n (dσ/dΩ) Ω within 1%; Monte Carlo within 4σ | [Rates and events](rates-and-events.md) | 39 |
| Monte Carlo spectra | ✅ | ⬜ | 🟡 TRIM transmitted spectra requested | ✅ peak positions and widths vs analytic (5%); seeded, thread-independent | [Rates and events](rates-and-events.md) | 39 |
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
| Seeded random numbers (thermostats, initial conditions) | ✅ | χ² uniformity, correlations; bit-identical across thread counts and restarts | 31 |

## Adding an entry

When a capability reaches **Engine**, give it a page in this folder (`kinematics.md`, `stopping-power.md`, …) with
these sections, and link it from the table:

- **Model**: the equations used, and any data (with version).
- **Assumptions and range of validity**: where the model holds and where it does not.
- **Validation**: for each comparison, the reference (tool and version, or paper), the cases compared, the
  tolerance, a plot or a link to the validation notebook, the test that enforces it, and the date and commit.
- **Known limitations**: what is missing or approximate.
- **References**: papers and data sources to cite.
