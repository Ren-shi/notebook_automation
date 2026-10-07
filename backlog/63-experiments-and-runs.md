# 63 · Experiments and runs: the model behind the workbench

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench (redesign)

## Why
The app grew one block at a time: nine separate "run" buttons (Simulate, Simulate the γ rays, Analyse, Check the
alignment, Solve with all orders, Show tracks, Run the source, Simulate its spectrum, Run sweep), each with its
own event count, seed and cache, each drawing nothing until pressed. In a lab there is one run and everything else
is a view of its data. Nothing today says which setup a result belongs to, and changing a field re-simulates
whatever block happens to be open. This item puts the lab's model under the app; items 64–70 are its tabs.

## Scope
- **An experiment** is the unit of work: a name (default "⁵⁸Ni on ²⁰⁸Pb", editable), a setup and a list of runs.
  On disk, one folder per experiment:

  ```
  58Ni-on-208Pb/
    experiment.toml             the current setup (the setup file as today)
    runs/
      01-152Eu-source/   setup.toml  events.npz  summary.json
      02-beam-10min/     ...
  ```

  `setup.toml` is the exact setup the run was taken with; `summary.json` holds the counters (kind of run, duration,
  real beam time simulated, rates, live fraction, counts per detector and per crystal, the coincidence matrix, and
  the Plan's predictions at that time, for item 68). The experiments live under `~/.physim/experiments` (or a
  folder the user picks); the app lists them on opening.
- **A run** is one of: a **beam run** of a duration; a **source run** (beam off, a calibration source at a
  position, item 53); an **alignment run** (a beam run with a deliberately wrong assumed geometry, item 60).
  Every run records the setup snapshot, the seed, what was simulated for real and what was scaled.
- **Real statistics for the beam time.** A beam run generates *unweighted* events — one event is one event — for
  the first part of the beam time, up to a budget ("simulate the first 10 min for real", default 10 min, choosable
  up to an hour), and scales the rest: histograms beyond the real part are the real part scaled to the full
  duration and redrawn with the full run's Poisson fluctuations. The run always says which part is real. (At the
  ⁵⁸Ni example's 4 800 particles/s, 10 min is 2.9 million events, about 30 s of CPU; an hour is 17 million and
  half a gigabyte; a day as real events is not offered.)
  - The generator's importance sampling and the γ chain's `gammas_per_event` are turned off within the real part
    (weights all 1); the analysis, alignment and uncertainty code keep working on unit weights.
  - **Stop** keeps what is accumulated; **Extend** continues the same seed stream and appends.
- **Compact columns.** Per particle event: detector, segment i, segment j, channel (16-bit integers), measured
  energy, θ_cm, depth (32-bit floats) and the event number; per γ ray: particle row, crystal, deposited and
  measured energy, lab θ and φ. Everything else (Doppler-corrected energies, gates, histograms, the coincidence
  matrix) is derived on demand. About 30 bytes per particle event and 25 per γ ray, stored as a compressed
  `.npz` (named arrays, readable from any Python). ROOT and CSV export stay as exports (item 67).
- **Locking.** While a run is in progress the setup is read-only (the panel greyed, Run becomes Stop). After it,
  editing is free again; earlier runs keep their snapshots, and every view of a run can say "taken with the
  target at 0.5 mg/cm²; the current setup has 1.0". A change to a label (a name) does not invalidate a run; a
  change to the physics does, and the Data tab says so instead of re-simulating.
- **The Planner** keeps its API (`Planner.rates`, `efficiency`, ...) for the Plan tab, and gains
  `Planner.experiment_folder`, `runs()`, `start_run(kind, duration, budget)`, `stop_run()`, `extend_run()`,
  `load_run(n)`; the per-block caches (`gamma_events`, `gamma_spectra`, `analysis`, `alignment`, `tracks`) read
  the current run instead of simulating on their own.
- **Memory:** one run's events loaded at a time; the summaries of all runs are small and always loaded.

## Depends on
55, 56 (the event chains it stores), 60 (alignment runs), 53 (source runs).

## Done when
- A beam run of 10 minutes on the ⁵⁸Ni example gives unweighted events whose counts per detector agree with the
  rates of item 42 within statistics, takes under a minute, and occupies under 50 MB on disk.
- A run folder reloads in a fresh session with the same spectra, counters and setup snapshot.
- Changing a physics field after a run marks the run stale in every view; changing a name does not; nothing
  re-simulates on its own.
- The scaled part of a long run is labelled wherever it is drawn, and a 24 h run's peak counts equal the rates'
  prediction within the run's Poisson fluctuations.
- Item 55's and 56's tests pass on runs with unit weights.
