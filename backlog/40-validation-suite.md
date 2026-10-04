# 40 · Validation suite against SRIM, LISE++ and literature

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner

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
