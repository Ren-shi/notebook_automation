# Count rates, beam time and Monte Carlo spectra

Backlog item 39 · modules `physim.nuclear.rates`, `physim.nuclear.events`, `physim.nuclear.plot`; Rust
`physim::nuclear` · tests `tests/python/test_nuclear_events.py`, `tests/nuclear_events.rs` · theory
{doc}`../theory/events`

## Model

- **Rates:** per detector and per segment, for every nuclide in the target and backing, ejectile and recoil.
  Rutherford's lab cross section is integrated over the face and averaged over target depth.
- **Beam time** for a number of counts (relative error 1/√N); counts in the planned run.
- **Warnings:** counting rate above 5000/s (configurable), too few counts, beam stopping in the target, Rutherford's
  limits.
- **Expected peaks:** mean measured energy, with widths from target thickness, kinematic broadening, beam energy
  spread, straggling and resolution.
- **Event generator (Rust, multithreaded, seeded):**
  - interaction depth, beam energy spread and spot;
  - straggling in and out;
  - CM angle drawn over the detectors' acceptance, with weights;
  - relativistic two-body kinematics (port of item 35);
  - first detector crossed;
  - dead layer, punch-through, resolution, threshold.
- **Spectra:** per detector, strip or channel, and θ against E pictures.

## Assumptions and range of validity

- **Elastic scattering only**, with the Rutherford cross section. The item 37 warnings say where that breaks down.
- **Two cuts keep Rutherford finite:**
  - lab angles below 0.5° are left out;
  - so are particles leaving the reaction below the lowest detector threshold (at least 10 keV).
- **Straight-line tracks:** no multiple scattering and no beam divergence.
- **Gaussian straggling**, carried with the item 36 straggling equation.
- **Analytic counted rates** judge the threshold from mean energies over the whole detector. The Monte Carlo decides
  per particle.

## Validation

| Check | Reference | Tolerance | Result | Test |
|---|---|---|---|---|
| Rust kinematics, 5 reactions (elastic, inverse, Q > 0, inelastic), 181 angles, both particles | Python `TwoBody` (item 35) | 1e-9° and 1e-11 relative | pass | `test_rust_kinematics_matches_python` |
| Energy after a layer, 5 ion/material pairs | `Stopping.energy_after` | 1e-9 | pass | `test_transport_tables_match_the_stopping_module` |
| Straggling from the tables | `Stopping.straggling` (step-by-step integration) | 2% | pass | same |
| Rate into a small disc at 30°, 60° and 135° | I n (dσ/dΩ) Ω | 1% | pass | `test_rate_of_a_small_detector_is_flux_times_cross_section_times_solid_angle` |
| Segment rates sum to the detector rate; every sector of a ring equal | Construction, symmetry | 1e-12, 1e-6 | pass | `test_rates_add_up_and_set_the_beam_time` |
| Monte Carlo against analytic rates: every detector, total and per channel and particle, three setups (one with a tilted target and a backing) | `Rates` | 4σ (2 × 10⁶ events; analytic counted rates +2%) | pass | `test_monte_carlo_rates_agree_with_analytic_rates` |
| Peak positions (α + Au, 20°–60°) | Kinematics plus mean energy loss (`Rates.peaks`); hand calculation at 30° | 4 standard errors or 1 keV | pass | `test_peak_positions_and_widths` |
| Peak widths | Resolution, straggling, kinematic broadening and target thickness in quadrature | 5% | pass | same |
| Same seed, same events; 1 thread and 4 threads; run split in pieces | Bit-identical | exact | pass | `test_same_seed_same_events_and_pieces_add_up`, `tests/nuclear_events.rs` |
| 10⁶ events in a few seconds | Benchmark `nuclear_events/alpha_on_gold_1e6` | — | 0.13 s (12 threads) | `benches/engine.rs` |
| Count rates against LISE++ | LISE++ | — | not yet: item 40 | — |

Validated 2026-10-04.

## Known limitations

- Rates are Rutherford rates even where the warnings say nuclear forces act, for example a carbon backing hit by a
  64 MeV ¹⁶O beam.
- Multiple scattering, beam divergence and non-Gaussian (end-of-range) straggling are not sampled.
- Identical beam and target nuclei are not given Mott interference.
- No background, random coincidences or beam halo.
