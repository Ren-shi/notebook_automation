# 39 · Count rates, beam time and the Monte Carlo event generator

**Priority:** P1 · **Size:** L · **Area:** Nuclear planner

## Why
These are the numbers a beam-time proposal needs: how many counts each detector collects per hour, how long to run
for a given statistical precision, and whether any detector will be swamped. The event generator then produces
realistic spectra, so students can see their peaks before they get beam, and test their analysis code on simulated
data (exported to ROOT by item 42).

## Scope
- **Analytic rates:** counts/s = (beam particles/s) × (target atoms/cm²) × ∫dσ/dΩ dΩ over each detector or strip.
  Counts in the planned run, and beam time needed for N counts (relative statistical error 1/√N).
- Rate warnings for detectors above a configurable count rate (pile-up, dead time, detector damage at forward
  angles).
- **Event generator** (Rust, multithreaded):
  1. sample the interaction depth in the target;
  2. slow the beam to that depth with straggling (item 36);
  3. sample the scattering angle from the cross section, restricted to the detectors' acceptance and weighted so the
     rates stay absolute;
  4. compute kinematics (item 35);
  5. slow the outgoing particle through the rest of the target, the backing and the dead layer;
  6. apply detector response (item 38).
  Each event records: detector, strip, true and measured energy, angles, interaction depth.
- Energy spectra per detector and per strip, and θ vs E plots.
- Seeded and reproducible: depends on item 31.
- Port the two-body kinematics of item 35 (Python, NumPy) to Rust for the per-event loop, checked against it.
- Left out: background, random coincidences, beam halo.

## Depends on
31, 35, 36, 37, 38.

## Done when
- Monte Carlo rates agree with the analytic rates within statistical error for every example setup.
- Peak positions agree with kinematics plus mean energy loss; peak widths agree with the combination of resolution,
  straggling, kinematic broadening and target thickness, within stated tolerances.
- Identical output for the same seed with 1 thread and with all threads.
- 10⁶ events in a few seconds on a laptop (benchmark added).
