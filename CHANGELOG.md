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
- The planner app has one mode (backlog 64): the guided steps, the guided/expert switch and the example drop-down are
  gone. The left panel is the setup as a list of decisions (Beam, Target and reaction, Particle detectors, γ-ray
  detectors, Run conditions), each a card showing its state and its consequence ("CD at 30 mm, behind the target ·
  126°–163°, 2.31 sr, 662/s"), with rarely touched fields behind "More" and the checks always at the top. A status
  strip names the experiment, beam, target, detectors, last run and warnings on every tab. The app opens on
  experiments: **New experiment** asks for the beam, target and what to measure; the examples are templates.
  **Import a setup file** and **Export the setup** replace Load and Save. The tabs are the stages: Setup, Plan,
  Run, Data, Analysis, Report, Physics. `physim.nuclear.workbench` has the cards, strip, checks and templates
  without the page; `physim.nuclear.app` is now a package with one module per tab.
- Particle–γ coincidence rates and beam times now include the angular correlation: each particle detector has its
  own γ-ray efficiency (`gamma_efficiency` in each row of `Planner.rates()`). Particle rates are unchanged.
- A γ-ray detector without `efficiency` now takes its efficiency from the response model, where it took its
  geometric coverage (an upper limit) before. Coincidence rates fall and beam times for Coulomb excitation grow
  accordingly: by a factor of about ten in the ⁵⁸Ni example. `gamma_efficiency_geometric` in `Planner.rates()` is
  now `gamma_efficiency_typical`.
- `resolution` of a γ-ray detector is its value at 1332 keV.
- After an edit, the rates of detectors that did not change, and that nothing shadows, are kept instead of being
  computed again.
- The planner app has a new look: a quieter header and setup panel, numbers in a monospace face, and no fixed white
  page. It still loads no fonts or scripts from the network.
- The planner app draws its Matplotlib figures off screen (Agg). The report builds them in a worker thread, where
  the on-screen backend could fail.
- `physim app` listens on 127.0.0.1 by default instead of every network interface (NiceGUI's default), so it is not
  reachable from the network and triggers no firewall prompt; `--host 0.0.0.0` restores the old behaviour.

### Added
- The Data tab (backlog 67, `physim.nuclear.dataviews`): every spectrum of the current run on screen at once, one
  panel per particle and γ-ray detector (per crystal on request), each enlarged with a click and saved as CSV. The
  operations are the user's: the Doppler correction (off, projectile, recoil, with why it is right or wrong for the
  setup), a particle gate, singles or coincidences with the randoms shown, subtracted or left out, add-back and
  suppression on or off, binning and range, and the whole run (scaled) or its real part. Gates are named objects
  (`Gate`: detectors, rings, the elastic or excited group by energy, an energy window). They are saved with the
  experiment (`gates.json`), and the analysis uses them by name (`Planner.analyse(gate=...)`, and
  `Settings.particle_energy`). The 2D views: energy against ring with each group's kinematic line (drag a box to
  fill in a gate), and γ-ray energy against crystal, raw or corrected. Two runs can be overlaid with their setups'
  differences named. Exports: a ROOT file with the gates as cuts, and the run's `.npz`.
- The Run tab (backlog 66): one Run button for a duration (the beam time by default), as a beam run, a source
  run (which source, activity, position) or an alignment check (the assumed offset). Before the run it shows the
  real part and its cost (`Planner.estimate_run`). While the run is taken: a clock and bar in beam time, a note of
  the part to be scaled, and live counters per detector (counts, rate, busy, live fraction). Per crystal it shows
  the singles and the coincidences; the particle × γ matrix fills in (full-energy peak / any energy), with bars
  for the counts wanted, and the last particles counted. Stop, Extend, Watch events (the run's events as tracks
  in the scene), the run list with the current run, and the list of what is not simulated yet (items 71–74,
  `workbench.NOT_SIMULATED`).
- The Plan tab (backlog 65): one page, answers first. It shows the beam time the counts wanted need against the beam time
  planned; per particle detector its rate, busy fraction (from `[run] dead_time`), share of the excitations and safe
  rings; per γ-ray detector its efficiency at the transition and its coincidence rate; and the particle × γ matrix.
  Each number has its ?. Below come the checks as panels that open: kinematics, energy loss, the orbit and safe
  distance, the correlation and Doppler shifts, the γ efficiency curves, and the rates per strip with the sweep.
  `Planner.plan()` returns the same numbers, and the next run's `summary.json` keeps them as its predictions.
- Experiments and runs (backlog 63, `physim.nuclear.runs`): an experiment is a folder (under
  `~/.physim/experiments`) with the current setup and its runs; each run keeps its exact setup (`setup.toml`), its
  counters and the Plan's predictions (`summary.json`) and its data as compact columns (`events.npz`, about 7 bytes
  per particle). A beam run generates unweighted events, one event being one event, for the first part of the beam
  time (10 min by default, up to an hour) and scales the rest; source and alignment runs too. Stop keeps what is
  accumulated; Extend continues the same seed stream. The setup is locked while a run is taken, and a run is
  marked stale when the physics of the setup changes (not a name). `Planner.create_experiment`, `open_experiment`,
  `start_run`, `stop_run`, `extend_run`, `load_run`, `runs`, `run_status`; the spectra, γ rays, analysis,
  alignment and tracks read the current run instead of simulating their own. A 10-minute run of the ⁵⁸Ni example
  takes about 50 s and 45 MB.
  - The event generator takes `skip_misses` (on for runs): a track that meets no detector is not carried out of
    the target, which saves time; the same seed then gives other events than without it.
  - `simulate_gammas(..., rates=)` reuses rates already computed; `gamma_events.crystals_of` and `backgrounds`
    build the crystals and the singles and randoms without simulating.
- Add-back and Compton suppression for clovers (backlog 62): `addback = true` sums the crystals' energies, so a
  Compton deposit returns to the full-energy peak with the share the `addback_factor` asks for (1.5 at 1332 keV,
  growing with energy) and the γ ray goes to the crystal with the larger deposit; `shield = "BGO"` rejects every
  deposit that is not the full energy with 1 − 1/S(E) (`suppression_factor`, 3 at 1332 keV) and stands as a
  blocking volume round the housing. The efficiency, the peak-to-total ratio, the spectrum shapes and the
  simulated γ rays follow; `simulate_gammas(..., plain=True)` and `Response(..., bare=True)` leave both out. In
  the app, two switches per clover, and the spectrum with and without.
- A movable target and sources placed anywhere (backlog 60): `[target] position` moves the target along the beam
  (every angle, distance and Doppler correction follows; the scene draws it there); `ladder` and `selected` for
  a target ladder; `source_run(..., position=...)` for a source off centre. `physim.nuclear.alignment` shows what a
  misplaced target does (the corrected peak with each geometry, the diagnostic plot of centroid against ring)
  and fits the offset back; **Check the alignment** in the app, with the assumed target as a ghost in the scene;
  the analysis can assume an offset and budgets the target's place.
- Tracks in the scene (`physim.nuclear.tracks`, backlog 59): **Show tracks** animates a sample of simulated
  events, beam, scattered beam, recoil and γ ray, to the segment or crystal each hit, which lights up; filters
  for coincidences or one channel; a sample picked by rate or as generated, stated in the scene; play, pause,
  speed; a selected track shows its event's numbers.
- Multi-step Coulomb excitation and reorientation (`physim.nuclear.coupled`, `physim.nuclear.multistep`, backlog
  58): the coupled equations for every substate of every level of a scheme, with the interface of first order;
  the yields of the planned experiment per detector and ring with all orders; the prolate–zero–oblate comparison
  of the 2⁺ state against the run's statistics; a fit of up to three matrix elements to measured counts; and a
  GOSIA input file for the setup. **Solve with all orders** in the app's "Excitation and γ rays" tab.
- The run record (`physim.nuclear.record`, backlog 57): every number of the plan and of the analysis with its
  formula, the formula with this run's numbers substituted, its physical meaning and its assumptions; the data's
  provenance; the method step by step; what is left out. The **?** beside a number in the app opens its
  explanation; **Show the run record** on the Report tab shows the page, which the report's zip holds as
  `record.html`.
- Automatic analysis of the simulated experiment (`physim.nuclear.analysis`, backlog 56): particle gate, Doppler
  correction, peak fit with random coincidences subtracted, yield, normalisation to the elastic particles or to a
  known transition, B(E2↑) with statistical and systematic uncertainties (a budget, each by varying its input),
  the Monte Carlo's own uncertainty, shape readings (β₂, Weisskopf units, rotor quadrupole moment), the
  comparison with the value put in, and the beam time for a wanted precision. **Analyse** on the app's Spectra
  tab; `Planner.analyse`.
- Particle–γ events (`physim.nuclear.gamma_events`, backlog 55): every simulated excitation whose particle reached a
  detector emits its γ ray with the angular correlation from the moving nucleus; the crystals record it with
  their response.
  - Doppler correction from the segment and crystal that fired, for the scattered beam and for the recoil.
  - True and random particle–γ coincidences from the singles rates and `[run] coincidence_window`; random γ–γ
    rates; `[run] dead_time`; room background (`[run] room_background`) and `[run] extra_lines`; a `threshold`
    on γ-ray detectors.
  - **Simulate the γ rays** on the app's Spectra tab; `write_root(..., gammas=...)` adds a `gammas` tree.
  - Units: `ns`, `us`, `ms` for times and `/s`, `Hz`, `kHz` for rates.
- Orientation of Coulomb-excited states and the particle–γ angular correlation (`physim.nuclear.orientation`,
  `physim.nuclear.angular`, backlog 54).
  - First-order excitation amplitudes for every magnetic substate of every level of a level scheme, for any
    ground-state spin, of the target or of the projectile.
  - Statistical tensors, γ-ray angular distributions in the nucleus's frame and in the laboratory, mixed
    transitions, internal conversion, and feeding from higher levels with the orientation carried down.
  - `[reaction] emission = "isotropic"` turns the correlation off.
  - The app's "Excitation and γ rays" tab shows the correlation factor of every particle-detector and crystal
    pair, and the cross sections of all levels and γ rays of a looked-up level scheme.
- γ-ray detectors have a response (`physim.nuclear.response`, backlog 53).
  - The full-energy-peak efficiency depends on energy, crystal size and distance, from a typical response per
    crystal material (germanium and LaBr₃). A measured `efficiency_curve` in the setup replaces it.
  - `absorbers` on a γ-ray detector, the chamber wall and the housing window attenuate the γ rays, with NIST
    attenuation coefficients.
  - The spectrum of a γ ray has a full-energy peak, a Compton continuum and escape peaks; the resolution depends
    on energy.
  - Calibration sources (²²Na, ⁶⁰Co, ⁸⁸Y, ¹³³Ba, ¹³⁷Cs, ¹⁵²Eu, DDEP data) can be run at the target position:
    `response.source_run`, and **Run the source** in the app, where the efficiency from each peak is drawn on the
    detector's efficiency curve.
  - Activities can be written with units (`37 kBq`, `1 uCi`).
- The planner app's Geometry tab is an interactive scene, drawn to scale (`physim.nuclear.scene`, backlog 52).
  - Detectors are solids with their rings, strips, crystals, circuit boards and housings; the target, the beam and
    the chamber are drawn too.
  - Clicking a detector, a ring or a crystal shows its angles, solid angle and rates in a side panel.
  - A selected detector can be dragged, in angle or in distance; one around the beam slides along it. The panel's
    numbers follow the drag.
  - A drag cannot end in the beam, inside another solid or through the chamber wall.
  - Rings closer than Cline's safe distance are marked, and the panel says where the kinematics send the particles.
- `Planner.move`, `Planner.place`, `Planner.live`, `Planner.selection` and `Planner.safety` do the same from Python.
- A detector catalogue (`physim.nuclear.catalogue`): `model = "S3"`, `"clover"` or `"LaBr3_2x2"` in a setup file fills
  in the real dimensions and segmentation, each with its source and with typical values marked.
  - A γ-ray detector can have real crystals: the four crystals of a clover each get their own Doppler correction.
  - Circuit boards and γ-detector housings stop particles, in the rates and in the simulated events.
  - An optional `[chamber]` section keeps particle detectors inside the chamber and γ-ray detectors outside it.
  - The app's Add menus offer the S3, the clover and the LaBr₃ detector.
- Level schemes (`physim.nuclear.levels`): the levels, γ-ray transitions and reduced matrix elements of the beam and
  target nuclei, read from a local copy of ENSDF (`physim.nuclear.ensdf`) or typed in the setup file as
  `[levels.beam]` and `[levels.target]`.
  - Every value records its source: ENSDF, derived, assumed or the user's. The user's values take precedence.
  - E1, E2 and E3 matrix elements are derived from B(Eλ) in Weisskopf units, or from half-lives and branching.
  - ENSDF is not shipped: `scripts/fetch_ensdf.py` downloads a copy to `~/.physim/ensdf`.
  - In the app, the *Excitation and γ rays* tab looks up a scheme, draws it, and plans excitation of a chosen
    state.
- The planner app has a **light and a dark theme** (**Light** / **Dark** in the header, or `?theme=dark`). It
  follows the system's setting the first time and remembers the choice. Figures are drawn in the theme's colours,
  with boxed axes and a colour-blind-safe palette; `app.themed(fig, theme)` does the same for a notebook.
- **Paper figures in a journal's style.** Every figure in the app has a **Paper figure** button: choose Physical
  Review, Nature, Science, Elsevier or Springer and a column width, and download a PDF, SVG or 600 dpi PNG at
  exactly that width, plus the plotted data as CSV.
  - The new module `physim.nuclear.paper` does the same from Python: `paper.figure`, `paper.export`,
    `paper.data_csv`, and the styles in `paper.JOURNALS`.
  - The figures are white with black, boxed axes; every curve has its own line style as well as its own colour,
    and text stays editable in the PDF and SVG files.
  - The beam-time report's figures use the same styles: `python -m physim.nuclear.report --journal nature --width
    double`, `report.build(exp, journal=..., width=...)`, or the journal choice in the app's **Report** tab. They
    are Physical Review, one column, by default, and the PNG files are now 600 dpi (they were 200 dpi).
- The planner app shows a strip of key results above the tabs: highest rate, longest beam time, closest approach
  and the number of warnings.
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
