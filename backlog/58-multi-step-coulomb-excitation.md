# 58 · Multi-step Coulomb excitation, reorientation and the shape of the nucleus

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

**Status: Done, except the comparisons with GOSIA, which need a Fortran compiler (the user's approval) and a
GOSIA installation outside the repository.**

> **Done** (`python/physim/nuclear/coupled.py`, `python/physim/nuclear/multistep.py`, `Planner.multistep` and
> `Planner.fit_matrix_elements`, the "Solve with all orders" block of the app, `tests/python/test_nuclear_multistep.py`,
> guide section "Multi-step excitation and reorientation" in `docs/nuclear-setup.md`).
> - **Decision: the coupled equations are in NumPy, not in the Rust core.** The integration applies the coupling
>   as a few real matrix products per step (the matrices are the Wigner–Eckart factors times the matrix
>   elements, with the orbit's phases as vectors), so 25 levels take 0.4 s per angle in NumPy, which meets the
>   speed asked for without a second language. The interface is item 54's: `CoupledChannels` is an `Excitation`.
> - **The orbit** is the symmetrised one of the ground state and a reference level, with ξ of each level from
>   its own CM speed, so the weak-coupling limit equals first order exactly; matrix elements enter with their
>   signs, the quadrupole moment through the diagonal E2 element.
> - **Integration over the experiment:** a grid of CM angles (4°) at 2 beam energies through the target,
>   interpolated to every quadrature direction of every segment; both the scattered beam and the recoil,
>   both kinematic solutions.
> - **Results:**
>   - With a ten-thousandth of the ⁵⁸Ni strength the result equals first order to 2 × 10⁻⁶ at three angles
>     (the real strength is 1.5% below first order at 170°, from the ground state's depletion).
>   - The probabilities of all levels sum to one to 10⁻⁸ (8 and 25 levels).
>   - With a diagonal E2 element the solver's amplitude equals first plus second order (nested quadrature on
>     the same orbit) to 5% of the second-order term; prolate lowers and oblate raises the yield, as the
>     reorientation effect requires.
>   - A 25-level scheme (265 substates) solves in 0.4 s per angle; the grid of one experiment (45 angles,
>     2 energies) in about 20 s, the prolate–zero–oblate comparison three times that. The real ⁵⁸Ni scheme
>     (6 levels) takes 76 s in the app for the yields, first order and the three shapes together.
>   - The fit recovers ⟨2‖M‖0⟩ to 2% and the quadrupole moment to 20% from counts with 1% errors.
>   - In the ⁵⁸Ni example the planned run tells prolate from oblate by 84 standard deviations over all rings.
> - **Not done here, or pending:**
>   - **GOSIA:** the three reference cases, the 1% comparison and the recovery of matrix elements from the
>     simulated yields need GOSIA on the machine, which needs a Fortran compiler. Nothing of GOSIA is in the
>     repository. The input-file writer follows the manual's format and has not been checked against a run.
>   - The E1 polarisation correction is not included (the formula was not at hand with enough confidence to
>     put a wrong one in).
>   - γ–γ coincidences from cascades: the yields of every transition are there; the angular correlation of two
>     γ rays is not.
>   - The analytic rates (`Rates`) still use first order with the elastic Rutherford factor; the multistep
>     yields use the symmetrised one (6% apart in the example). Which to use everywhere is best settled with
>     the GOSIA comparison.

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
