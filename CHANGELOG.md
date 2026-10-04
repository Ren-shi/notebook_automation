# Changelog

All notable changes are listed here. The project follows [semantic versioning](https://semver.org/):
before 1.0, minor versions (0.x) may change the API; patch versions only fix bugs.

The version lives in `Cargo.toml` and is shared by the Rust crate, the Python package
(`physim.__version__`) and the PyPI distribution `physim-engine`. To release: bump the version in
`Cargo.toml`, move the "Unreleased" notes under it here, merge, and push a tag `vX.Y.Z`; the release
workflow builds the wheels, tests them and publishes to PyPI.

## Unreleased

### Fixed
- A particle detector hidden behind another no longer reports counts in the analytic rates, peaks and beam times.
  - The integral leaves out directions that meet another detector first, as the event generator always did.
  - A partly hidden detector is integrated on a finer grid, so the shadow's edge is resolved within 0.5%.
  - A fully hidden detector's warning names the detector in front, and shadowing is now a warning, not a note.

### Changed
- The planner app draws its Matplotlib figures off screen (Agg). The report builds them in a worker thread, where
  the on-screen backend could fail.
- `physim app` listens on 127.0.0.1 by default instead of every network interface (NiceGUI's default), so it is not
  reachable from the network and triggers no firewall prompt; `--host 0.0.0.0` restores the old behaviour.

### Added
- Coulomb-excitation beam times count what the measurement uses: excitation events seen in a particle detector
  together with their γ ray.
  - γ-ray detectors gain an optional `efficiency` (full-energy peak); without it, their geometric coverage is used
    as an upper limit.
  - `Rates.rate`, `counts_in_run`, `beam_time_for` and `relative_error` take `what="all" | "excitations" |
    "coincidences" | "measured"`.
  - The rates table and the report say what is counted.
  - For `coulex_ni58` the beam time for 2000 coincidences is 4–6 h, where it was 3 s for 2000 particles of any kind.
- The planner app's **Add** offers a ready-made particle detector for each region: forward strips, a side pad, a
  backward pad, or a backward ring around the beam. The geometry reading points out a setup with nothing backward,
  and which detector sees the most excitation events.
- The planner app opens in a **guided mode**:
  - a step-by-step workflow (what to measure, beam, target, particle detectors, γ-ray detectors, rates and beam
    time, spectra, report), each step showing the results it affects, with **Next** waiting until the step's
    problems are fixed;
  - a **?** help on every input (what it is, a typical value, what raising it does);
  - a **How to read this** box on every result, written with the setup's own numbers.

  The previous layout is the **Expert view**, and the browser remembers the choice. The guidance lives in
  `physim.nuclear.guide` (`HELP`, `STEPS`, `reading`).
- Setup errors suggest a unit that suits the field ("'12' has no unit; write it with one, e.g. '12 h'"), instead of
  always "40 mm".
- Coulomb excitation can be set up entirely in the planner app:
  - a **Reaction** section (elastic or Coulomb excitation, with the excited nucleus, multipolarity, state energy
    and B(Eλ↑));
  - a **γ-ray detectors** section (add, edit, duplicate, remove), shown in the 3D geometry and its table;
  - a table of particle energies per particle detector, elastic and after excitation, with β of the excited nucleus;
  - particle detectors are labelled as silicon.

  `Planner` gains `set("reaction", …)`, `set("gamma detector N", …)` and `add_`/`remove_`/
  `duplicate_gamma_detector`.
- A one-click Windows installer for the planner app (`physim-planner-<version>-windows-x64-setup.exe`, attached to
  each GitHub release): the official embeddable Python with physim and the app's packages, a per-user install with a
  Start-menu shortcut, no Python or administrator rights needed. Built and test-installed by the new `Installer`
  workflow.
- `physim app --desktop` (what the shortcut runs): logs to `%LOCALAPPDATA%\physim\planner.log`, reuses a planner
  that is already running, moves to a free port if 8080 is taken, and stops a minute after the last browser tab
  closes (`--idle-exit`).
- Coulomb excitation (second slice of the nuclear planner):
  - `physim.nuclear.coulex.Coulex`: first-order semiclassical excitation probability and cross sections for E1, E2
    and E3, with the safe-distance check;
  - `physim.nuclear.gamma`: Doppler-shifted γ-ray energies and broadening per pair of particle and γ detector;
  - setup files gain `[reaction] type = "coulex"` and `[[gamma_detectors]]`, with the example `coulex_ni58`;
  - excitation channels in rates, peaks and the Rust event generator (weighted by P(θ));
  - the app's "Excitation and γ rays" tab, a report section with `gamma.csv`, and validation checks.
- The experiment planner web app (NiceGUI): `physim app` opens the planner in the browser.
  - A setup panel with a detector table, a warnings banner, and tabs for geometry (3D), kinematics, rates and beam
    time (per-strip heat maps, parameter sweeps), energy loss, spectra, trajectories and the report.
  - An "Explain" panel on every tab.
  - Install with `pip install physim-engine[app]`; step-by-step guide in `planner-app`. A new `physim` command
    (`physim app`, `physim report`).
- Beam-time report and data exports (`physim.nuclear.report`, `python -m physim.nuclear.report setup.toml`):
  - a self-contained HTML report with print styles for PDF;
  - CSV tables with units in the column names;
  - figures as PNG and PDF, and the setup file;
  - the same setup and seed reproduce every number;
  - with `uproot` installed (`pip install physim-engine[root]`), `events.root`: an event TTree, TH1D spectra with
    Monte Carlo errors, and the setup (`physim.nuclear.rootio`).

  `plot.kinematics` draws E against θ with detector coverage.
- `physim.nuclear.planner.Planner`: the planner app's model. It covers setup editing with validation, every result
  tab as plain data, warnings, "Explain" texts and one-parameter sweeps, so the web app is a thin layer over it.
- Validation suite for the nuclear planner:
  - `physim.nuclear.validation`: one tool comparison and one literature comparison per physics-register capability,
    with `run()`, `report()` and plots; a test keeps the register's ✅ marks honest and checks that broken formulas
    are caught.
  - Validation notebooks in the docs; a LISE++ nuclear-mass reference.
  - `physim.nuclear.export`: `trim_in` and `trim_in_for` write SRIM `TRIM.IN` files, and `lise_settings` lists
    LISE++ inputs with physim's own values.
- Count rates, beam time and Monte Carlo spectra for the nuclear planner:
  - `physim.nuclear.rates.Rates`: rates per detector, strip and channel (ejectile and recoil, target and backing),
    counts in the run, beam time for N counts, expected peaks with their widths broken down, and warnings.
  - `physim.nuclear.events.simulate`: a seeded, multithreaded Rust event generator (`physim::nuclear`) giving
    weighted per-particle events and absolute spectra.
  - `plot.spectra` and `plot.theta_energy`, plus a benchmark.
- Seeded, reproducible randomness: `World.seed` (saved in checkpoints; the Langevin thermostat's default seed), the
  counter-based generator made public (`physim::rng`, `physim.random.uniform`/`normal`), and `physim.random` helpers
  for Maxwell–Boltzmann velocities and uniform positions in a box or ball. Noisy runs are bit-identical for any
  thread count and across checkpoint restarts.
- Nuclear experiment planner, first part: `physim.nuclear.Experiment`, a setup file (TOML) describing the beam,
  target, detectors and run conditions, with every value written with its unit and every problem reported by field
  name. Example setups (`Experiment.example`), the guide page `nuclear-setup`, and the physics register
  (`physics-register`) recording what physics is implemented and how it is validated.
- Nuclear and material data (`physim.nuclear.data`): AME2020 atomic masses for 3558 nuclides, Q-values, nuclear
  masses; natural isotopic compositions, element densities and mean excitation energies, and 48 compounds from NIST;
  materials from elements, isotopes, formulas or compound names, with areal-density and atoms/cm² conversions and
  enriched isotopes. Setup files gain an optional `density` and check that every material and nuclide exists.
  Data sources and licences are in `physim/nuclear/data/SOURCES.md`; `scripts/convert_nuclear_data.py` rebuilds the
  NIST tables from the downloads.
- Relativistic two-body reaction kinematics (`physim.nuclear.kinematics.TwoBody`): lab energies and angles of
  ejectile and recoil from CM or lab angles (both solutions where double-valued), maximum angles, thresholds,
  excitation energies, the solid-angle Jacobian and kinematic broadening dE/dθ; `elastic(experiment)` for every
  target isotope. Theory note `theory/kinematics`; LISE++ comparison cases written to
  `tests/reference/nuclear/pending/` for running by hand.
- Stopping powers and energy loss (`physim.nuclear.stopping.Stopping`): NIST PSTAR/ASTAR tables (shipped, 74
  materials) for protons and α particles, interpolation for other elements, Bragg additivity for compounds,
  effective-charge scaling with a Lindhard–Scharff low-velocity limit for heavy ions, ZBL nuclear stopping; ranges,
  energy after a (tilted) layer, Bohr straggling through thick layers, Highland multiple scattering. Checked against
  NIST ranges and 1680 LISE++ values generated by `scripts/lise_reference.py`; theory note `theory/stopping`.
- Rutherford scattering (`physim.nuclear.rutherford.Rutherford`): CM and lab cross sections (ejectile or recoil,
  both kinematic branches), integrated cross sections, closest approach and impact parameters, grazing angle,
  Coulomb barrier, Sommerfeld parameter, screening and Mott warnings, and Coulomb orbits integrated by the engine.
  Checked against the analytic formulas, LISE++ and Geiger and Marsden's 1913 gold data; theory note
  `theory/rutherford`.
- Detector geometry and response (`physim.nuclear.detectors`): rectangular strip, disc and annular faces from the
  setup file; ray hits with strip/ring/sector and incidence angle; solid angles per detector and segment by
  quadrature; θ/φ coverage; shadowing and beam-blocking warnings; dead layer, punch-through, resolution and
  threshold; the exit path through the target. `physim.nuclear.plot.setup_3d` and `coverage` draw a setup.
- Quantum mechanics: `Schrodinger`, the time-dependent Schrödinger equation for one particle on 1–3D grids, by
  split-step Fourier (spectral in space, second or fourth order in time, exactly unitary) or Crank–Nicolson (any
  boundary, unitary, direct 1D and BiCGSTAB 2D/3D solves); static or time-dependent potentials (arrays or Python
  functions), absorbing layers, Gaussian wave packets, observables, recorded runs, and bound states by LOBPCG.
  `physim.quantum` helpers (momentum densities, probability currents, analytic barrier transmission and energy
  levels), `plot.wavefunction` and `plot.animate_wavefunction`, the theory note `theory/quantum`, and example notebook
  `08_quantum` (wave packets, tunnelling, bound states, Ehrenfest's theorem, a 2D double slit).
- Example notebook `07_accelerators`: electrostatic acceleration, a drift-tube linac with phase stability, classical
  and isochronous cyclotrons (protons and alpha particles, relativistic dynamics via `CustomForce`), and Rutherford
  scattering of alpha particles on gold.
- Fields on grids: wave equation (leapfrog), heat equation (explicit, Crank–Nicolson), Poisson (FFT, conjugate
  gradients) in 1–3 dimensions with Dirichlet, Neumann or periodic boundaries; particle-mesh gravity (isolated or
  periodic) as a force.
- Packaging: wheels for Linux (x86_64, aarch64), macOS (x86_64, arm64) and Windows (x64) and a source
  distribution, tested on Python 3.9 and 3.13 without a Rust toolchain; PyPI publishing on version tags;
  MIT `LICENSE`; this changelog.

## 0.1.0 (first release)

### Engine (Rust)
- `World` with particles (mass, charge, radius, pinning), forces and integrators; fixed-step `run` with
  recording, events (located to round-off inside a step), streaming sinks and checkpoints that restart bit-exact.
- Integrators: Euler, symplectic Euler, Verlet, Yoshida 4/6/8, PEFRL, Blanes–Moan, RK4, Dormand–Prince 5(4)
  (also adaptive with dense output), Gauss–Legendre 2/4/6, Wisdom–Holman, Boris, Langevin (BAOAB),
  Nosé–Hoover; user-defined compositions and splittings.
- Forces: gravity (direct and Barnes–Hut tree), springs and spring networks, drag, central potentials
  (power law, Yukawa, Plummer, Hernquist, Hénon–Heiles), post-Newtonian and J2 corrections, electric and
  magnetic fields, Coulomb, Lennard-Jones / Morse / tabulated pair potentials in periodic boxes with neighbour
  lists, soft contact; custom forces in Rust or Python.
- Rigid rod constraints (RATTLE); hard collisions with walls (event-driven); chaos indicators (Lyapunov
  spectrum, MEGNO); rigid bodies (`RigidSystem`: quaternions, torques, symplectic rotation splitting).
- Deterministic parallel force evaluation (results independent of the thread count).

### Python
- Bindings for everything above, trajectory and checkpoint files (`.npz`, HDF5), analysis helpers
  (energy error, Poincaré sections, radial distribution), `physim.scenarios` (two-body, solar system,
  Plummer sphere, figure-eight, pendulum, spring lattices), `physim.units` (SI, astronomical, custom) and
  `physim.plot` (orbits, energy error, phase space, animations, plotly 3D view).
