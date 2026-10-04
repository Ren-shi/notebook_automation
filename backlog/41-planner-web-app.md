# 41 · Experiment planner web app

**Priority:** P1 · **Size:** L · **Area:** Nuclear planner

## Why
The main users are graduate students planning experiments, many of whom do not want to write code. They need a
point-and-click tool: describe the setup, see it, get the numbers and the report. The app is a thin layer over the
setup file (item 33) and the physics (items 35–39), so the same results are available from Python.

## Scope
- Runs locally in the browser, started with one command (`physim app`) or a desktop shortcut. No accounts, no
  server to host.
- **Framework decision first:** Panel or NiceGUI (both Python). Build the same small mock-up (beam form + detector
  table + 3D plot) in each and choose on: interactive 3D plotting, form and table editing, responsiveness with a
  running Monte Carlo, and packaging. Record the decision in this file.
- **Setup panel:** beam, target, detector table (add, remove, duplicate, edit), run conditions; load and save setup
  files; example setups to start from.
- **Result tabs:**
  - Geometry: 3D view of the beam, target and detectors.
  - Kinematics: E vs θ for scattered particles and recoils, with the detector coverage shaded.
  - Cross sections and rates: per detector and strip, with beam time needed.
  - Energy loss: through the target, backing and dead layers.
  - Spectra: from the event generator.
  - Trajectories: Coulomb orbits from the Rust engine.
  - Report: preview and export (item 42).
- **Validity warnings** (item 37) and rate warnings (item 39) shown prominently, not hidden in a tab.
- **"Explain" panel** on every output: the formula, the assumptions, where it stops being valid, and a link to the
  physics register entry. This is the start of the learning track for undergraduates.
- **One-parameter sweep:** vary beam energy, a detector angle or target thickness, and overlay the results.
- Left out: multi-user hosting, 3D game-engine visualisation, undergraduate guided labs (planned after this slice).

## Depends on
33, 35–39 (tabs can be built as each item lands).

## Done when
- From a fresh install, following written steps without writing code, a user can build a four-detector setup,
  see every tab, and export a report.
- Tried by at least one graduate student who has not seen the project, with their feedback recorded and acted on.
- Automated smoke tests load each example setup and render every tab without errors.
- Installation for users who do not use Python is documented, and a one-click installer is at least planned
  (tracked as its own item if not done here).
