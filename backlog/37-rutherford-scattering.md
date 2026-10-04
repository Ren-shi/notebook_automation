# 37 · Rutherford scattering, closest approach and Coulomb trajectories

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner

## Why
The first slice is elastic scattering, and below the Coulomb barrier that is Rutherford scattering. The planner has
to give the cross section at each detector, the distance of closest approach, and, just as important, a clear
warning when the setup is outside the range where Rutherford's formula holds. It is also the natural bridge to the
existing engine: the Rust core can integrate the actual hyperbolic Coulomb trajectories, which turns the formula
into something a student can see.

## Scope
- Rutherford cross section dσ/dΩ in the CM frame and converted to the lab frame with the Jacobian from item 35.
- Distance of closest approach as a function of scattering angle, with the head-on value.
- Validity checks, shown as warnings in the app and the report:
  - closest approach compared with an interaction radius R = r₀(A₁^⅓ + A₂^⅓) plus a margin (r₀ and the margin
    configurable, defaults documented), so the user sees where nuclear forces start to matter;
  - Coulomb barrier height compared with the beam energy, and the grazing angle;
  - the Sommerfeld parameter η, to flag when the classical picture is unreliable;
  - electron screening at very low energies and small angles.
- Integrated cross section over a detector's angular acceptance (used by item 39).
- Trajectories: a set of projectiles at chosen impact parameters integrated with the existing engine's Coulomb force,
  for the app's trajectory view, with the impact parameter ↔ angle relation b = (d₀/2) cot(θ/2).
- Left out: optical-model elastic scattering above the barrier (a later item, needed once the warnings above fire);
  Mott scattering of identical particles (noted as a limitation until added).

## Depends on
34, 35.

## Done when
- Cross sections match the analytic formula exactly in the CM frame and the lab conversion matches an independent
  numerical calculation.
- Deflection angles from integrated trajectories match the analytic b(θ) relation to 1e-6.
- The relative angular distribution of the Geiger–Marsden α + Au data is reproduced, and one published elastic
  data set above the barrier is shown to deviate where the warning says it should.
- Theory note in `docs/theory/`, and an entry in the physics register.
