# 38 · Detector geometry and response

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner · **Status: Done**

> **Done** (`python/physim/nuclear/detectors.py`, `python/physim/nuclear/plot.py`,
> `tests/python/test_nuclear_detectors.py`, `docs/theory/detectors.md`, `docs/physics-register/detectors.md`).
> - `Geometry` (face, axes, segments, `hit`/`hits`, `solid_angle` per face or segment, `solid_angle_monte_carlo`,
>   `theta_range`, `phi_range`, `mean_theta`, `outline`), `Array.from_experiment` (`first_hit`, `shadowing`,
>   `warnings`, `response`), `Response` (dead layer, deposited energy with punch-through, `punch_through_energy`,
>   `measured_energy` with resolution and threshold), `exit_path` through a tilted target with a cap near 90°.
> - Solid angles by 48 × 48 Gauss–Legendre quadrature: exact to round-off for any placement, so the closed forms
>   and Monte Carlo are tests rather than code paths. Whole-face θ extremes found on the edge by dense sampling
>   plus golden-section refinement (the first version sampled only rectangle corners and was 0.01° off).
> - Pictures: `setup_3d` (beam, target, faces, strips, labels) and `coverage` (outlines in θ–φ, unwrapped at ±180°),
>   checked by eye for both examples.
>
> Results (15 tests): solid angles match closed forms (disc, centred and off-centre rectangle, annulus and all 384
> of its segments) to 1e-10 and a Monte Carlo count for tilted faces within 4σ; θ ranges to 1e-9°; shadowing of a
> disc in front of another equals Ω_front/Ω_back to 0.1%; punch-through energies equal the stopping ranges; the
> example array's only warning is that CD hides < 1% of DSSD4.
>
> Left out (as planned): ΔE–E telescopes and particle identification, timing, pulse-height defects, inter-strip
> effects; beam spot size is left to the event generator (item 39).

## Why
A setup is only useful if the planner knows what each detector sees: its angular coverage, solid angle, and how it
turns a particle's energy into a measured signal. These numbers go straight into count rates (item 39) and the
beam-time report (item 42).

## Scope
- Detector shapes: rectangular pads (optionally split into strips), annular strip detectors (CD-type), circular
  apertures. Placement by (θ, φ, distance) with the face normal to the target by default, or a full position and
  orientation.
- Solid angle and angular coverage (θ range, φ range) per detector and per strip; analytic where a formula exists,
  numerical integration otherwise.
- 3D geometry for the app: detector outlines, beam axis, target, and the angular coverage on a unit sphere.
- Response: energy lost in the dead layer (item 36), punch-through when the detector is thinner than the particle's
  range (deposited energy then reported, not full energy), Gaussian energy resolution (FWHM, optionally
  energy-dependent), threshold.
- Overlap and shadowing checks between detectors: warn only, no shadowing calculation.
- Left out: ΔE–E telescopes and particle identification (a later item), timing, pulse-height defects in silicon for
  heavy ions (listed as a known limitation).

## Depends on
33, 36.

## Done when
- Solid angles match analytic results for an on-axis disk (2π(1 − cos α)) and an on-axis rectangle, and Monte Carlo
  estimates for off-axis cases, within 1e-4 relative.
- Punch-through energies agree with range tables from item 36.
- The geometry plot of each example setup is checked by eye and kept as an image test.
