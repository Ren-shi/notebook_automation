# 34 · Nuclear and material data

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner

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
