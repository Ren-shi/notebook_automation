# 48 · Hidden detectors counted in the analytic rates

**Priority:** P1 · **Size:** S · **Area:** Nuclear planner

**Status: Done.**

## Why
Found while testing item 47 (2026-10-05). A backward pad added behind the CD drew the warning "CD hides 100% of B1
from the target", yet the rates table still gave B1 7.9 counts per second and a beam time.
- **Monte Carlo:** the event generator stops each particle in the first face it meets, so it was right.
- **Analytic rates:** they integrated every detector's whole face as if it stood alone.

## Done
- **Visible directions:** `Geometry.distances` (vectorised) and `Array.visible`. `Rates` drops the quadrature
  directions that meet another detector first, in the rates and in the peaks.
- **Finer grid:** a partly hidden detector is integrated on about 48 × 48 points over its face (`SHADOW_ORDER`).
  - With the default 12 × 12 the shadow's edge costs up to a few percent; with the finer grid it is within 0.5%.
  - A segmented detector's grid is refined per segment only as far as needed, so `oxygen_on_lead_array` (DSSD4 0.3%
    hidden) still computes in about 2.5 s.
- **Warnings:**
  - a fully hidden detector says "P records nothing: CD stops every particle before it gets there";
  - shadowing is listed as a warning, not a note.
- **Tests:**
  - `test_visible_part_of_a_hidden_detector` checks the exact difference of two cones within 0.5%;
  - `test_hidden_detectors_count_only_what_reaches_them` checks a half-hidden pad against the same pad alone, and a
    fully hidden one that counts nothing;
  - the Monte Carlo comparison gains a `half_hidden` setup, and agrees within 4σ.
