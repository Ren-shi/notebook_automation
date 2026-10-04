# 39 · Count rates, beam time and the Monte Carlo event generator

**Priority:** P1 · **Size:** L · **Area:** Nuclear planner · **Status: Done**

> **Done** (`python/physim/nuclear/rates.py`, `python/physim/nuclear/events.py`, `python/physim/nuclear/plot.py`,
> `src/nuclear/`, `src/python/nuclear.rs`, `tests/python/test_nuclear_events.py`, `tests/nuclear_events.rs`,
> `benches/engine.rs`, `docs/theory/events.md`, `docs/physics-register/rates-and-events.md`).
> - `Rates(experiment)`:
>   - `rate`, `per_detector`, `per_segment`, `by_channel`, `counts_in_run`, `beam_time_for`, `relative_error`;
>   - `peaks` (mean measured energy and width, with its parts: target thickness, kinematic broadening, beam
>     energy spread, straggling, resolution);
>   - `warnings` (above 5000/s by default, too few counts, beam stopping in the target, Rutherford's limits, layout).
>
>   The target and the backing, every nuclide in each, and both ejectile and recoil are all included. Rates are
>   averaged over target depth.
> - Rust `physim::nuclear`:
>   - `TwoBody`, a port of item 35, equal to Python to 1e-11;
>   - `Table`, energy loss with straggling from range, S and W = ∫(dΩ²/dx)/S³ dE tables, which solve the item 36
>     straggling equation in closed form;
>   - `Face` (detector hits);
>   - `Generator`: multithreaded in fixed blocks; random numbers keyed by (seed, 2⁶¹ + event, draw).
> - Python: `physim.nuclear.events.simulate` returns weighted `Events` (per-particle columns, `rate`, `counts`,
>   `spectrum`). `Stopping` gained `straggling_rate` and `transport_table`. `Geometry.directions` provides the
>   quadrature. `plot.spectra` and `plot.theta_energy` draw the results.
> - **Sampling.** CM angles are drawn over the detectors' acceptance, half from Rutherford's 1/u² and half flat in
>   ln u, with weights f/p. Backward strips get events, and spectra stay absolute.
> - **Finite rates.** Lab angles below 0.5° and particles leaving the reaction below the lowest threshold (at least
>   10 keV) are cut. This removes the 90° Rutherford recoil divergence.
> - **Done-when results:**
>   - Monte Carlo against analytic rates within 4σ for every detector, per channel and in total, for both example
>     setups and a tilted one.
>   - Peak positions agree with kinematics plus mean energy loss (and a hand calculation within 1 keV).
>   - Peak widths agree with the quadrature sum within 5%.
>   - Output is bit-identical on 1 and 4 threads, and for a run split in pieces.
>   - 10⁶ events take 0.13 s in the benchmark, and 0.2–0.4 s from Python, tables included.
> - **Left for later:**
>   - comparing count rates with LISE++ (item 40);
>   - multiple scattering and beam divergence in the generator;
>   - reaction types beyond elastic (item 43).

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
