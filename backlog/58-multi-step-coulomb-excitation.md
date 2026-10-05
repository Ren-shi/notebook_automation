# 58 · Multi-step Coulomb excitation, reorientation and the shape of the nucleus

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

## Why
First-order theory fails once excitation probabilities reach a few percent, and it cannot describe the
reorientation effect at all. Reorientation is how the quadrupole moment of the 2⁺ state is measured, and its sign
separates prolate from oblate. That measurement is the main goal of the experiments this tool plans. Today the app
warns "use GOSIA" at that point.

## Scope
- **Coupled equations:** the amplitudes of every magnetic substate of every level are followed along the orbit,
  coupled by all the matrix elements of the level scheme (item 50), in the Rust core behind the interface of
  item 54. Multi-step excitation, reorientation and interference between paths all follow from this.
- **Size:** 20–30 levels.
- **E1 polarisation correction** (virtual excitation of the giant dipole resonance).
- **Orbit:** the symmetrised orbit, as now. Projectile and target are solved separately.
- **Integration over the experiment:** over the energy loss in the target and the angles each ring covers, on a
  grid of energy and angle that is interpolated.
- **Decay:** cascades through several levels, with γ–γ coincidences.
- **Fit:** up to about three matrix elements are adjusted until the calculated yields match the peak areas of
  item 56, with the others fixed. Uncertainties come from the fit.
- **Prolate or oblate:** the same experiment is run with the quadrupole moment at the prolate rotor value, at zero
  and at the oblate rotor value. The app shows the yield in each ring for the three cases, and whether they differ
  by more than the uncertainties.
- **GOSIA input:** the app writes a GOSIA input file for the setup: levels, matrix elements, detector geometry and
  yields.
- **Left out:** GOSIA's full minimiser; deorientation; lifetimes; nuclear interference above the safe energy;
  mutual excitation of both nuclei.

## Design notes
- Written from the published theory (Alder and Winther; the GOSIA manual as a reference for conventions). No GOSIA
  code is copied, shipped or committed: it has no licence that permits redistribution, and this repository is
  public.
- **Reference numbers:** GOSIA is run on the machine, outside the repository. The repository keeps our input
  files, the output numbers and a note of the GOSIA version; the tests compare against those numbers.
- GOSIA needs a Fortran compiler, which this machine lacks. Installing one needs the user's approval.
- Cygnus (GPL v3) is a possible second reference. It is not embedded either.

## Depends on
50, 54, 56.

## Done when
- With weak coupling, the result equals first order to 1e-6.
- The probabilities of all levels sum to one to 1e-8.
- A two-level case with a diagonal matrix element reproduces the published second-order reorientation formula.
- Populations and yields agree with GOSIA within 1% for three reference cases, including a light beam on ¹⁹⁴Pt at
  backward angles.
- GOSIA, given the simulated yields, recovers the matrix elements that were put in, within the uncertainties.
- A 25-level case solves fast enough for the grid of one experiment to take under a minute on the reference
  machine.
