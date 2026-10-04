# 41 · Experiment planner web app

**Priority:** P1 · **Size:** L · **Area:** Nuclear planner · **Status: Done, except the trial with a graduate
student (the user)**

> **Done** (`python/physim/nuclear/planner.py`, `python/physim/nuclear/app.py`, `python/physim/__main__.py`,
> `tests/python/test_nuclear_planner.py`, `tests/python/test_nuclear_app.py`, `docs/planner-app.md`).
> - **Framework: NiceGUI**, chosen by the user (2026-10-04); the two-framework mock-up comparison was skipped.
>   NiceGUI 3.17 with Plotly figures gives interactive 3D, editable forms, background tasks for the Monte Carlo
>   (`run.io_bound`), and plain pip packaging (`physim-engine[app]`).
> - **`physim app`** (or `python -m physim app`) starts the planner in the browser.
>   - **Header:** start from an example, load or save the setup file.
>   - **Setup panel:** beam, target, backing and run fields; a detector list with add, duplicate and remove, and
>     fields that follow the shape. A value is applied on Enter or when the field loses focus.
>   - **Warnings banner** above the tabs, errors first.
>   - **Tabs:** Geometry (3D, plus a coverage table), Kinematics (coverage shaded), Rates and beam time (per-strip
>     heat map and a one-parameter sweep), Energy loss, Spectra (re-simulate with any number of events and seed),
>     Trajectories, and Report (downloads the report zip from item 42).
>   - **Explain panel** on every tab.
> - Everything the app shows comes from `Planner` (`physim.nuclear.planner`), so the same results are available
>   from Python.
> - **Tests:**
>   - every figure renders for every example;
>   - `test_the_app_serves_every_example` starts the real server and loads each example page, which builds every
>     tab on the server (an error there gives a 500);
>   - `Planner` tests cover editing, sweeps and every tab.
> - **Checked by hand in a browser** (2026-10-04, the oxygen-on-lead example):
>   - Every tab rendered.
>   - Raising the beam to 6 MeV/u dropped DSSD1 from 14,130/s to 6,262/s, as 1/E² predicts (6,280). The above-barrier
>     warnings for lead appeared.
>   - A negative energy showed the error and kept the last valid results.
>   - Add made a new detector.
>   - The Plotly resize message in the console for plots on hidden tabs is harmless.
> - **Installation without Python** is documented step by step in `docs/planner-app.md`; the one-click installer is
>   planned as item 44.
> - **Waiting on the user:** a trial with a graduate student who has not seen the project, with their feedback
>   recorded here.

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
