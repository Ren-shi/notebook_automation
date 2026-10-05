# The experiment workbench: overview of items 50–62

Agreed in discussion on 2026-10-05. This file records the aim and every decision; the numbered items hold the
detail. Nothing here is built yet.

## Aim
The planner so far returns results: tables, spectra and a report. The workbench turns it into a simulated
experiment that the user can see and rearrange.
- **Who it is for:** an experimentalist planning beam time. Students come second.
- **First complete case:** Coulomb excitation of a light beam on a heavy target (normal kinematics), for example on
  ¹⁹⁴Pt or ²⁰⁸Pb.
  - An S3 silicon detector at backward angles records the scattered projectile.
  - Germanium clovers and LaBr₃ detectors record the γ rays.
  - The result is a Doppler-corrected γ spectrum whose peak areas give B(E2), the quadrupole moment of the 2⁺
    state, and from those the shape of the nucleus.
- **Graphs support the scene.** Clicking a detector shows its numbers; clicking the experiment shows the totals and
  coincidences.

## Decisions

### Scene and hardware
- **To scale and segmented:** detectors are solids with their real dimensions and segmentation (the S3's 24 rings
  and 32 sectors, a clover's four crystals).
- **Parametric shapes, no mesh import:** each model is generated from numbers in a catalogue, so the picture and
  the physics use the same geometry.
- **Faithfulness:** active volumes plus simple blocking shapes. Nothing is transported inside a detector.
- **Generic chamber:** the user sets a chamber radius and beam pipe. No facility is built in.
- **Constrained movement:** a detector cannot sit in the beam, overlap another, or cross the chamber wall.
- **Kinematic advice:** the app recommends backward angles for a light beam on a heavy target and forward angles for
  the reverse, and warns when a detector sits where the kinematics send nothing.
- **Movable target, and γ calibration sources** (¹⁵²Eu, ⁶⁰Co) that can stand in for the beam.
- **Tracks:** a sample of simulated particles and γ rays is animated through the scene.

### Nuclear data
- **Level schemes come from ENSDF**, bundled with the app so it works offline, for both beam and target.
- **The user's values take precedence.** Every value is marked as from ENSDF, derived, or supplied by the user.
- **Any ground-state spin** is supported in the data model, though even-even nuclei are the first case.

### Physics
- **Safe Coulomb excitation only,** with Cline's 5 fm criterion shown per ring.
- **Particle–γ angular correlation** by default, with isotropic emission as a switch.
- **All decays in flight after the target,** for now.
- **γ detector response is parametrised:** an efficiency curve scaled by geometry, attenuation in material between
  target and crystal, and a response function with a full-energy peak, Compton continuum and escape peaks.
- **Realistic background:** Compton continua, random coincidences from the rates and the time window, room lines,
  and lines the user adds by hand.
- **Electronics:** thresholds and a coincidence time window. Dead time is a simple correction.
- **Multi-step excitation is our own calculation** (semiclassical coupled channels), sized for 20–30 levels, with
  reorientation and the E1 polarisation correction, and a fit of up to about three matrix elements.
- **Order of the refinements GOSIA includes:** internal conversion first, then γ–γ cascades, then deorientation.

### Analysis and output
- **Automatic analysis in the app:** particle gates, Doppler correction and peak fits, all adjustable by the user.
- **Normalisation:** to target excitation (for example ¹⁹⁴Pt) or to the scattered particles (for example ²⁰⁸Pb,
  which is hardly excited). The app suggests which.
- **Uncertainties:** statistical and systematic from the first version.
- **Shape:** the size of β₂ from B(E2), prolate or oblate from the sign of the quadrupole moment, with the
  statement that these readings depend on the rotor model.
- **Run record:** every calculation carries its formula, the numbers substituted and its physical meaning, pitched
  at an MSc student. Viewable in the app and saved as a PDF.
- **Exports:** a plain ROOT file of events, and a GOSIA input file for the setup.

### GOSIA
- GOSIA is free to use with a citation, but it has no licence that permits redistribution.
- **It is never shipped and never committed.** This repository is public, so GOSIA stays on the machine, outside
  the repository.
- **The repository holds only our own work:** the GOSIA input files for the test cases, the numbers they produced,
  and a note of the GOSIA version.
- **Workflow:**
  1. a known case, where our calculation and GOSIA must agree;
  2. a simulated experiment, from matrix elements we chose;
  3. GOSIA analyses the simulated yields and should recover those matrix elements.
- Running GOSIA needs a Fortran compiler, which this machine does not have. Installing one needs the user's
  approval.

### Application
- All calculation runs on the user's machine. The app stays in the browser for now; a desktop window is the goal.
- **The repository stays public** under the MIT licence. Nothing planned here requires making it private, as long
  as GOSIA stays out.
- **Bundling ENSDF awaits written confirmation from the NNDC.** The data are free to use with a citation, but no
  statement that permits redistribution was found. Item 50 has the details and the fallback.

## Order of work
| Step | Items | What the user gets |
|---|---|---|
| 1 | 50, 51 | Level schemes on choosing beam and target; a catalogue of real detectors |
| 2 | 52 | The to-scale scene with selection and constrained dragging |
| 3 | 53, 54 | Real γ detection and correctly oriented γ emission |
| 4 | 55, 56 | Coincidence events, Doppler-corrected spectra, and B(E2) with uncertainties |
| 5 | 57 | The run record |
| 6 | 58 | Multi-step excitation, reorientation and the prolate–oblate comparison |
| 7 | 59–62 | Tracks, movable target, desktop window, add-back |

## Parked
- **Educational mode:** historical experiments (Rutherford, Compton, Chadwick) as limited-movement scenes with
  modern equipment. It reuses the scene of item 52.
- **Lifetimes:** decays inside the target or a backing, with stopped and shifted components.
- **Deorientation** of the angular correlation in vacuum.
- **Radioactive beams:** mixed beams and decay background from stopped beam. The design leaves room for both.
- **GRSISort export layout.**
- **A facility case study** (iThemba LABS, AFRODITE) as a validation of the whole chain against published data.
- **Triaxial and shape-coexisting nuclei.**
