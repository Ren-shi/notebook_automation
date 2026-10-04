# 38 · Detector geometry and response

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner

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
