# 34 · Nuclear and material data

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner · **Status: Done**

> **Done** (`python/physim/nuclear/data.py`, `python/physim/nuclear/data/` with `SOURCES.md`,
> `scripts/convert_nuclear_data.py`, `tests/python/test_nuclear_data.py`, `docs/physics-register/nuclear-data.md`).
> - AME2020 `mass_1.mas20.txt` shipped unchanged (3558 nuclides, `#` estimates flagged); `nuclide()`, atomic and
>   nuclear masses (Lunney–Pearson–Thibault electron binding), `q_value()` with charge and mass-number checks.
> - NIST isotopic compositions (288 natural isotopes), element densities and mean excitation energies (Z ≤ 92),
>   and 48 compounds, converted to CSV by a script so they can be rebuilt from the downloads.
> - `material()`: elements, isotopes (density scaled by atomic mass), formulas with D/T, NIST compounds by short or
>   full name; `areal_density_mg_cm2`, `thickness_um`, `atoms_per_cm2`; `Material.enriched`.
> - Setup files gain `density` (target, backing, detector) and now reject materials and nuclides missing from the
>   data, neutron beams, and length thicknesses without a known density.
> - Licences: AME2020 is CC BY 3.0; NIST data redistributed with attribution, modification notes and disclaimer.
>
> Results (19 tests): 8 published Q-values reproduced within 0.5 keV (tolerance 1 keV); α and proton nuclear
> masses match CODATA 2018 to 0.2 keV and 1 eV; atomic weights of all 84 natural elements up to U match NIST's
> standard atomic weights (inside the interval or 1σ) except Se (1.5σ) and Hg (2.4σ), whose standard weights were
> revised in 2013, after the abundances used; areal densities for Au, C, Si, CH₂ and CD₂ match hand calculations;
> the data files are in the built wheel.
>
> Left out: nuclear level schemes (ENSDF, for item 43); isotopic enrichment in the setup file (Python only).

## Why
Kinematics need accurate masses, energy loss needs material compositions and densities, and a planner is only as
trustworthy as its input data. The data has to be versioned and cited so an auditor can see exactly where every
number came from.

## Scope
- Atomic mass evaluation (AME2020): mass excess for every nuclide, with nuclear masses derived by removing electron
  masses and binding. Lookup by name (`"197Au"`), by (Z, A), or by element with natural abundance.
- Element data: symbol, Z, standard atomic weight, density, natural isotopic composition.
- Materials: elements, compounds by formula (`"CH2"`, `"CD2"`, `"Mylar"`, `"CsI"`), user-defined compounds by
  mass or atom fraction, isotopically enriched targets. Areal density conversions (mg/cm² ↔ µm ↔ atoms/cm²).
- Each data file stored in the repo with its source, version, retrieval date and licence noted in a header, and
  listed in the physics register.
- Left out for now: nuclear level schemes (ENSDF), needed from item 43 onwards.

## Depends on
—

## Done when
- Q-values for a list of reactions computed from the tables match published values to 1 keV.
- Areal density conversions match hand calculations for Au, C, Si, CH₂ and CD₂.
- Every data file has provenance and licence recorded, and the licences permit redistribution.
