# 45 · Coulomb excitation set up entirely in the app

**Priority:** P1 · **Size:** S · **Area:** Nuclear planner

**Status: Done.**

## Why
From the first look at the app (2026-10-05): a Coulomb-excitation setup could only be written in a setup file.
- The reaction settings and γ-ray detectors were not in the setup panel.
- The γ detectors did not appear in the geometry.
- The silicon particle detectors were not labelled as such.
- The particle energies needed for the Doppler correction were spread over two tabs.

## Scope
- **Reaction section** in the setup panel:
  - elastic or Coulomb excitation;
  - the excited nucleus (target or beam), multipolarity, state energy and B(Eλ↑).

  Switching to Coulomb excitation fills in E2 excitation of the target. Switching back keeps the state for later.
- **γ-ray detectors** section (shown for Coulomb excitation): add, edit, duplicate and remove. They are drawn in the
  3D geometry and listed in a table there.
- **Labels:** "Particle detectors" are labelled as silicon, and γ detectors as germanium.
- **Particle-energy table** on the excitation tab, for each particle detector and each particle that reaches it, at
  the smallest, central and largest angle:
  - the elastic energy;
  - the energy after exciting the state, and the difference;
  - β of the excited nucleus.

## Depends on
41, 43.

## Done when
- A Coulomb-excitation plan can be built from *alpha_on_gold* without editing a file.
- The geometry shows the γ detectors.
- The particle-energy table agrees with the kinematic factor by hand.

> **Done** (`planner.py`: `set("reaction", …)`, `add/remove/duplicate_gamma_detector`, γ detectors in `geometry()`,
> `gamma()["particles"]`; `app.py`: the Reaction and γ-ray-detector sections, the geometry table, the
> particle-energy table).
> - **Tests** (`tests/python/test_nuclear_planner.py`):
>   - `test_coulomb_excitation_is_editable` builds a plan from *alpha_on_gold*, toggles elastic and back, and
>     round-trips the file;
>   - `test_gamma_detectors_in_the_geometry` checks that each outline is the crystal's circle facing the target;
>   - `test_particle_energies_for_coulomb_excitation` checks the elastic energies against the kinematic factor
>     within 1%, a positive energy difference, and β from the excited nucleus's own energy.
> - **Also:** the app draws its Matplotlib figures off screen (Agg). The report builds them in a worker thread,
>   where the default on-screen backend warned that it could fail.
