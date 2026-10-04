# 40 · Validation suite against SRIM, LISE++ and literature

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner · **Status: Done, except references the user has to
produce**

> **Done** (`python/physim/nuclear/validation.py`, `python/physim/nuclear/export.py`,
> `tests/python/test_nuclear_validation.py`, `notebooks/validation/` with `scripts/make_validation_notebooks.py`,
> `docs/validation.md`, `tests/reference/nuclear/masses_lise.csv`, `rutherford_lisepp.csv`, `pending/trim/`,
> `pending/rates_lisepp.md`).
> - **`physim.nuclear.validation`:** 18 checks, one tool and one literature check per register capability.
>   - The exceptions are Coulomb trajectories and detector solid angles, which have exact closed forms; no tool
>     computes these for arbitrary setups.
>   - Each check states its reference, provenance and tolerance, and returns pass, fail or pending with the worst
>     deviation and plot data.
>   - `run()` and `report()` give the results; all 18 run in about 1.5 s.
>   - Shared reference data now lives here too (published Q-values, Geiger–Marsden 1913, the LISE++ stopping
>     tolerances); the older tests import it.
> - **New tool reference: nuclear masses** from LISE++'s `isotope_mass` (47 nuclides, `scripts/lise_reference.py
>   --masses-only`).
>   - 46 agree within 2 keV, or within the AME2020 uncertainty where that is larger. The worst is 57% of the allowed
>     deviation (¹⁰⁰Sn 135 keV of its 240 keV; ²³⁸U 1.1 keV).
>   - ⁷⁸Ni is skipped: AME2020 only estimates it (±400 keV), and LISE++ differs by 990 keV.
>   - Rerunning LISE++'s straggling gave different values from the committed file (the library keeps a straggling
>     setting from its GUI), so that file was not regenerated. This is documented in the script.
> - **LISE++'s Rutherford values** (cross sections, d₀, b(60°), grazing angle) moved from a test into
>   `rutherford_lisepp.csv`, with a provenance header.
> - **The register is tested.**
>   - `test_register_ticks_are_backed_by_passing_checks` parses `docs/physics-register/README.md` and fails if a ✅
>     there is not backed by passing checks.
>   - `test_a_broken_formula_fails_its_check` breaks 12 formulas in turn and requires the matching checks to fail:
>     masses, kinematics, stopping power, range, straggling, the Rutherford normalisation and angular shape, the
>     impact parameter, orbit deflection, solid angle, depth averaging of rates, and generator straggling.
> - **Validation notebooks:** six, one per register page, built (executed) with the docs. Each shows a results table
>   and a plot per comparison.
> - **Exports:**
>   - `export.trim_in` and `trim_in_for(experiment)` write TRIM.IN files: multi-layer targets, enriched isotope
>     masses, transmitted ions saved.
>   - `export.lise_settings(experiment)` lists what to enter in LISE++, alongside physim's energies, cross sections
>     and rates.
>   - The four TRIM cases requested in item 36, and the two example setups, are written to `pending/trim/`.
> - **Waiting on the user (the register shows 🟡):**
>   - the LISE++ kinematics cases (`pending/kinematics_lisepp.md`);
>   - LISE++ cross sections at the detector angles (`pending/rates_lisepp.md`, for rates);
>   - TRIM runs, which also check that the exported TRIM.IN files open and run in SRIM (`pending/stopping_srim.md`,
>     for spectra);
>   - the Farwell–Wegner above-barrier data (`pending/rutherford_above_barrier.md`).
>
>   Each has a check that reads its file as soon as it exists.

## Why
Physicists trust LISE++, SRIM and Geant4 because they have been checked against data for decades. Physim aims to be
a starting point that brings these calculations together, not a replacement, so every number it produces must be
traceable to a comparison with one of those tools or with published data. This item is the machinery behind the
physics register (`docs/physics-register/`).

## Scope
- **Reference data** in `tests/reference/nuclear/`, one file per comparison, each with a header recording: the tool
  and version (or the paper and table), the exact inputs and settings, the date, and who generated it. SRIM and LISE++
  are run by hand (they are closed GUI programs), so the inputs must be complete enough for anyone to reproduce them.
- **Automated tests** in Rust and Python that compare Physim with each reference file and assert the tolerance stated
  in the register. These run in CI, so a regression turns a register entry red.
- **Validation notebooks** in `notebooks/validation/`, one per capability, with comparison plots and residuals.
  Built with the docs so the plots stay current.
- **Export to the trusted tools** so users can check Physim's result in one click: SRIM/TRIM input files (`TRIM.IN`)
  and LISE++ settings for the same setup.
- Process: tool comparisons first to get going, then literature and experimental data as the final standard.
  The register records which level each capability has reached.

## Depends on
33 (reference cases are written as setup files). Grows alongside 35–39.

## Done when
- Every capability in the elastic-scattering slice has at least one tool comparison and one literature or data
  comparison, each with an automated test and a validation notebook.
- A deliberately broken formula makes the matching test and register entry fail.
- Exported `TRIM.IN` files open and run in SRIM.
