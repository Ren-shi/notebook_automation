# Changelog

All notable changes are listed here. The project follows [semantic versioning](https://semver.org/):
before 1.0, minor versions (0.x) may change the API; patch versions only fix bugs.

The version lives in `Cargo.toml` and is shared by the Rust crate, the Python package
(`physim.__version__`) and the PyPI distribution `physim-engine`. To release: bump the version in
`Cargo.toml`, move the "Unreleased" notes under it here, merge, and push a tag `vX.Y.Z`; the release
workflow builds the wheels, tests them and publishes to PyPI.

## Unreleased

### Added
- A reproduced Coulomb-excitation experiment in the validation (Allmond 2014 verification, fix wave 3, S6): Allmond et al., PRC 90, 034309 (2014), ⁶²Ni at 111.6 MeV on carbon, as the example `coulex_ni62_c_allmond` and the check `reproduction_allmond2014` (`physim.nuclear.coulex_reproduction`, `scripts/reproduce_allmond2014.py`; the paper's numbers, the GOSIA references and the bands in `tests/reference/nuclear/coulex/allmond2014.toml`). The thickness gives the measured ⁵⁸Ni/⁶⁰Ni losses; σ_C/σ_R per ring within 0.2% of GOSIA with nuclear masses and DIP off and 1.3% of GOSIA as run (outside its 1%: a known, explained miss, see below); closure; the paper's per-ring B(E2), rebuilt through GOSIA's ratio(M) (the paper gives no yields), recovered 1.2–1.3% low, inside its 4.1%; the Q_s and feeding sensitivities within 0.25 percentage points of GOSIA's.
- Outside programs run out of sight (`physim.external`): GOSIA, SRIM's SR Module and TRIM, PACE4 and GEMINI++
  start on a hidden Windows desktop with no console window, so nothing shows or takes the focus, and in a job
  object that stops them and everything they started at the timeout or when physim exits. One runner holds this;
  a thin module per program feeds it, starts it and collects its output files. The PACE4, GEMINI++ and SRIM
  reference scripts use it, and the GEMINI++ UI Automation driver moved into the package (`gemini_drive.ps1`).
- Fusion-evaporation (backlog 75): `[reaction] type = "fusion-evaporation"` plans the residues of a compound
  nucleus, from the exit channels' cross sections to the residues' Doppler-corrected γ rays.
  - Stage 1 (`physim.nuclear.residues`): the channel table (`[[reaction.channels]]`), the compound nucleus's exact
    E* and speed through the target, events per channel with the particles evaporated one by one and the residue's
    recoil from each, and the residue's cascade entering its known levels by spin, aligned perpendicular to the
    beam. The Plan tab lists the residues, their rates and strongest lines; the Data tab their γ-ray spectra. Example:
    `fusion_dy158` (¹²C on ¹⁵⁰Nd at 60 MeV).
  - The fusion cross section (`physim.nuclear.fusion`): the Bass barrier, Wong's transmission summed over partial
    waves, σ_ℓ and the average over the target. The Plan tab shows the barrier and the excitation function.
  - The evaporation Monte Carlo (`physim.nuclear.evaporation`): one Hauser–Feshbach chain per event (n, p, α,
    giant-dipole γ rays, fission), with Gilbert–Cameron level densities. PACE4 4.36.0 runs in silent mode and was
    compared for five systems: light beams agree within 15%; heavy beams at high spin do not yet (see
    `docs/theory/evaporation.md`).
  - "Fill from physim" in the reaction's card (and the MCP tool `fill_channels`) replaces the channel table with
    the fusion and evaporation results (`physim.nuclear.exit_channels`); the setup keeps its shape.
  - Neutron detectors (`physim.nuclear.neutrons`) and an optional `[neutron_wall]`: the Plan tab's channel
    selection shows which detected class (neutrons, protons, α particles seen) each channel's events fall in, and
    what a gate on each class contains (MCP tool `channel_selection`). The selection applies the particle
    detectors' thresholds (each detector's `threshold`, else `[neutron_wall] particle_threshold`, default 1 MeV,
    after the dead layer and with punch-through) and γ rays firing neutron cells (`gamma_multiplicity`, default
    J/2 + 4 from the entry spin, `gamma_efficiency`, `gamma_as_neutron`), and returns the analytic matrix of the
    same events beside the Monte Carlo.
  - The neutron data file's carbon cross section is ENDF/B-VIII.0 (natural carbon, bin-averaged;
    `scripts/endf_neutron_reference.py`), checked against the evaluation in
    `tests/reference/nuclear/neutron_carbon_endf.csv`. The earlier typical values were up to 31% off (3.4 MeV), which
    moved a 150 mm cell's efficiency by up to 12%.
  - The report of a fusion-evaporation run has a "Fusion and evaporation" section: σ_fus, the residues and the
    channel selection, each with the validation checks of its model beside it (a known failure said plainly, with
    the suite's numbers for the compared system nearest the run's beam and compound nucleus) and the documented
    limits that apply, flagged where the run falls in their region (near or below the barrier, above 10 MeV/u,
    Z₁Z₂ near or above the quasi-fission range, fusion spins above 40 ħ, p and α channels in the selection).
    "Models and their validation" lists every fusion-evaporation check with its checker, version, tolerance and
    outcome, picked by capability, so checks added to the suite appear without a change to the report. The report
    now runs only the checks of the models it shows.
  - A published experiment reproduced end to end through the planner: Usman et al., NIM B 469, 42 (2020),
    ⁴He + ¹⁶⁵Ho stacked-foil activation (`physim.nuclear.reproduction`, `scripts/reproduce_usman2020.py`, new
    example `fusion_tm165`). Each of the 26 Ho foils is the example at the foil's energy, with its channels from
    "Fill from physim", its residues simulated and the production cross sections inferred from them. The new
    literature check `reproduction_usman2020` passes: the peaks of 2n, 3n and 4n within the energy spread + 1 MeV,
    the strongest channel right at every foil the data resolve, the strongest channels within 30% of the band of
    three measurements. The tails above the peaks (pre-equilibrium emission) are reported as the documented limit.
    PACE4 at the same energies is in `pace4_excitation_4He_165Ho.csv`, and a validation notebook draws all four
    excitation functions.
- Heavy-ion stopping from MSTAR 3.12 (Paul and Schinner; the method of ICRU Report 73) for ions from Li to Ar in
  hydrogen, in the ICRU 49 compounds and in their deuterated versions such as CD₂ (`physim.nuclear.mstar`).
  Against 693 measured points in H₂/D₂, polymers and Mylar, the RMS deviation for C–Ar falls from 7–10% to
  4–9%. A new literature check, `stopping_hydrogen_vs_data`, compares heavy-ion stopping with measured data
  from the IAEA Stopping Power Database; it still fails for Kr and Xe in hydrogen-rich targets.
- A planner warning when an ion heavier than Ar, below 10 MeV/u, crosses a layer or detector containing hydrogen
  or helium: there physim's stopping is up to 20% below measured data (`stopping.known_shortfall`).
- The angular correlation between the γ rays of one cascade (backlog 71). Each excitation carries the density
  matrix of its level down the cascade (`physim.nuclear.correlation`): a γ ray is drawn along a direction with
  the probability its level's state gives, and leaves the next level in the state that direction selects, from
  which the next γ ray is drawn. Before, every γ ray of a cascade was drawn on its own from its level's average
  orientation, so the γ–γ coincidences between two crystals missed the correlation.
  - The emission is exact for pure and mixed multipoles (helicity amplitudes, `Emission`), and it matches the
    angular distributions of `Populated.w` and the unobserved feeding of `fed_density` to 1e-15: summed over
    directions nothing changes, so the singles are those of before (63 line × crystal bins of the ⁶²Ni on ¹⁹⁴Pt
    template agree with the previous version, mean pull −0.07, rms 0.92). An unoriented 4⁺ → 2⁺ → 0⁺ cascade
    gives 1 + 0.1020 P₂ + 0.0091 P₄.
  - A cascade's first level starts in the density matrix the collision leaves it in at the event's angle
    (`Populated.density`, `Excitation.density_table`), not the average over the levels feeding it: a 2⁺₁ reached
    directly and one reached through the 2⁺₂ are oriented differently, and each event now has its own.
  - A converted transition hands its orientation on as an unobserved one. `emission = "isotropic"` keeps drawing
    every γ ray on its own, without the correlation, and a setup with one excited state is unchanged bit for bit.
  - Only the last γ ray of a cascade is aimed at the crystals; those before it are drawn over the whole sphere,
    since the next one depends on their direction. The γ rays of a setup with cascades take about twice as long
    (23 s for 300 000 reactions of the ⁶²Ni on ¹⁹⁴Pt template, against 12 s).
- Lifetime measurements, plunger (RDDS) and Doppler-shift attenuation (DSAM), backlog 73: planned, simulated and
  analysed for Coulex and fusion-evaporation, from the setup to τ and B(Eλ).
  - A `[lifetime]` section in the setup for Coulex and fusion-evaporation: a plunger (stopper, degrader, or a
    differential plunger with a second foil, `set_degrader2`) at a distance with an offset, or a DSAM backing,
    with problems for impossible setups (a stopper too thin to stop the recoil, distances not positive). Units
    `ps`, `fs`, `/ns`, `/ps`. Examples after the papers: `coulex_ru98_plunger`, `fusion_cr50_plunger`,
    `fusion_sm143_dsam`.
  - Lifetime physics (`physim.nuclear.lifetimes`): B(Eλ) ↔ τ, decay curves with Bateman feeding, the plunger's
    shifted and unshifted components, DSAM attenuation and line shapes on a slowing-down kernel, and the
    deorientation in vacuum (`lifetime.deorientation`: α₂, α₄ with either the rates or g, C₂, C₄). Theory note
    `docs/theory/lifetimes.md`.
  - Recoil transport (`physim.nuclear.transport`): the recoil's path through target, backing, vacuum gap and foils
    (`rates.recoil_stack`), with electronic stopping and a ZBL scattering Monte Carlo checked against TRIM
    (range, straggling, lateral spread and angles within a few %). The slowing of item 73 moved here unchanged.
  - Coulomb-excitation γ rays through a plunger or DSAM backing, with lifetimes, time dilation and deorientation in
    vacuum (new γ columns `stopped`, `segment`, `time_ps`, `vacuum_ps`, `distance_um`, `stopper_ps`); plunger runs
    at a distance (`runs.start(distance_um=...)`), grouped as distance series.
  - Fusion-evaporation residues carry level lifetimes and side-feeding times (LINESHAPE rotational band,
    `side_feeding_q` and `side_feeding_levels` in `[reaction]`). With a `[lifetime]` section they slow through the
    target, backing and plunger before each γ ray, and there are predicted plunger components and spectra per
    distance and DSAM line shapes.
  - Simulated fusion-evaporation runs (`physim.nuclear.fusion_runs`, or `runs.start(kind="fusion")`): residue
    events with their γ rays, triggered on γ-ray fold and the ancillaries' gates; stored, started, stopped,
    extended and loaded like Coulex runs; γγ matrices and gated spectra; a plunger distance per run.
  - Lifetime analysis (`physim.nuclear.lifetime_analysis`): plunger peak areas, decay-curve fits with Bateman
    feeding, the differential decay-curve method (singles with feeding correction, gated, and three-foil for
    differential plungers, with the degrader and stopper points gated as Beckers et al.'s Eq. 13 requires), DSAM
    line-shape fits, B(Eλ) from τ, and the deorientation correction and fit. Each result has an uncertainty budget
    (statistics, distance, velocity, feeding, stopping).
  - The planner (`physim.nuclear.lifetime_plan`): `Planner.lifetimes` (per level the method that suits it, τ and its
    source, plunger distances and sensitive range, the peak separation per γ-ray detector, F(τ), the beam time per
    distance for a wanted precision, and warnings, including a Coulex plunger with no deorientation),
    `Planner.lifetime_spectra` (predicted spectra per distance and angle) and `Planner.analyse_lifetimes` (τ from
    a distance series by DDCM, three-foil DDCM or decay curve, or from DSAM line shapes, then B(Eλ↓) beside ENSDF,
    the Coulex B(E2) and the shell model); `dataviews.distance_spectra`. The beam-time report has a lifetime
    section with its limits and the lifetime checks; MCP tools `plan_lifetime_measurement`, `set_plunger`,
    `lifetime_spectra` and `analyse_lifetimes`; docs page "Lifetime measurements" (`docs/lifetimes.md`).
  - The planner app: a plunger (RDDS) or DSAM card in Setup for Coulex and fusion (distances, offset, foil kind,
    material and thickness, a differential plunger's second foil, deorientation as rates or g with C₂, C₄), and the
    fusion Plan's "Add a plunger…" and "Line shapes (DSAM)…" buttons. "Lifetimes" is a live kind on the start page
    with three templates. Plan has a lifetimes panel and predicted spectra per distance; Run takes a run at a chosen
    distance, or one run per distance not yet taken, and takes fusion runs; Data groups the distance series
    (spectra per distance and ring) and shows fusion runs (singles, γγ matrix, gated spectra); Analysis shows the
    decay curve and DDCM τ(x), or the line-shape fit, with τ, its budget, B(Eλ) and the comparison rows.
  - Validation: two capabilities, "Lifetimes: recoil distance (RDDS)" and "Lifetimes: Doppler-shift attenuation
    (DSAM)", with analytic, Bateman, deorientation (Radeck et al.) and TRIM checks that pass, and two known
    failures said plainly: against APCAD the DSAM line's shift is 17% off (widths agree), against AhKin its width
    is 25% off (shift inside AhKin's band). The DDCM agrees with NaPyTau on simulated spectra, and τ comes back
    from simulated runs within its uncertainty (closure, RDDS and DSAM). Published experiments reproduced through
    `physim.nuclear.lifetime_reproduction` (`scripts/reproduce_lifetimes.py`): Radeck et al.'s ⁹⁸Ru 2⁺₁ (8.35 ±
    0.27 ps against 8.36(29)) and Beckers et al.'s ⁵⁰Cr 2⁺₁ and 4⁺₁ pass; Rajbanshi et al.'s ¹⁴³Sm DSAM is a known
    failure (fitted with SRIM's stopping, as the paper did, τ comes out about 20% long: SRIM's range of Sm in Au is
    19% longer than physim's). Register page `docs/physics-register/lifetimes.md`, notebook
    `notebooks/validation/lifetimes.ipynb`.

### Documentation
- The README shows the application: what the experiment planner does, step by step, with screenshots
  (`docs/images/`), how to start it (Windows installer or `physim app`), the example setups and `physim report`.
  The engine follows in a shorter summary; its full reference moved from the README into `docs/guide.md`, which
  the documentation site shows unchanged.
- The README's screenshots show the app's new look (the step rail, the dark theme), with a light-theme Setup step
  in `docs/images/planner-setup-light.png`, and the README lists the clover template.

### Fixed
- The Allmond 2014 reproduction says what it misses (Allmond 2014 verification, second fix plan, 9 and 13). The band
  against GOSIA as run is the 1 % the verification's plan had fixed; the 2 % in the reference file had been written
  after about 1.3 % was found, as had the 0.5 % against GOSIA with nuclear masses and DIP off (0.2 %), and the file says
  so now. No band is raised. physim is +1.32 % (ring 2) and +1.16 % (ring 3) above GOSIA as run: the check
  `reproduction_allmond2014` fails, as a known and explained miss (GOSIA's integer masses, E1 polarisation,
  non-relativistic kinematics), worst 1.32 of its band (0.46 without the miss), and the tests assert its size
  (±0.05 %). All orders are what is judged; first order is a second line. The closure is the mean of four seeds fixed
  beforehand (1001–1004) instead of seed 11, to all orders, with its Monte Carlo and scatter errors, within a quarter of
  the paper's statistical error + 3 Monte Carlo errors (it was 1 % + 3 of one seed): −0.25/−0.57/−0.48 % for
  2.6/1.6/1.4 %. New: the closure sample as a stored run, and a run taken in the app, give B(E2) back ring by ring
  (`coulex_reproduction.stored_run`). `Result.monte_carlo` against the scatter of 100 seeds: 1.08(9), 1.12(7), 1.06(7)
  for ring 2, ring 3 and both, no significant difference, so it is unchanged.
- The report's section on the Doppler correction says at what depth of the target it takes the reaction
  (`gamma_events.reaction_depths`) instead of "the middle of the target", and that a segment stands for the mean
  direction of its particles.
- The analysis of a stored run averages the expected excitation over each segment (Allmond 2014 verification, second
  fix plan, 2). A run knows only the segment a particle hit, and ⟨P⟩/B was read at the segment's mean angle: 0.6 and
  1.2 % too large with one ring per detector (⁶²Ni on carbon, 14–28° and 28–44°). Now the cross sections of the
  excitation and of every particle are each integrated over the segment, for each kind of particle, without the
  directions in another detector's shadow and the particles below the threshold; a segment across the largest lab
  angle of a double-valued particle is integrated in the CM system. Within 0.02 % of the cross sections integrated
  over the ring, with seven rings and with one. Weighted simulations keep their particles' angles, bit for bit.
- The analysis's uncertainty budget (fix 12): a "target thickness" term, the change of ⟨P⟩/B with the thickness changed
  by its uncertainty (`Settings.thickness_unc`, a field of the Analysis step; −2.05/+2.10 % for ±1.3 % of the Allmond
  target), and "matrix elements assumed", which was always 0: the change of the γ rays in the peak when the level
  scheme's other matrix elements move by their uncertainties, from the coupled equations (Q_s ± 0.1 eb on ⁶²Ni:
  +2.87/−2.81 % and +2.02/−1.99 %). An input without an uncertainty has no term, and the new `Result.budget_notes`,
  the normalisation step, the record and the report say so. `coulex_ni62_pt194_clovers` (Q_s = 0.05 ± 0.12 eb from
  ENSDF and four transitions) goes from 10 to 26 % systematic. Results kept by earlier versions load as they are.
- A stored run's particle angles and the Doppler correction take the right angle of each segment (Allmond 2014
  verification, second fix plan, 1). `Geometry.segment_centre(weighted=True)` averaged x and y over a sector, which
  put the point 1.1–1.6 % nearer to the axis (a ring without sectors at 0° or 180°), and weighted recoils by the
  scattered beam's 1/sin⁴(θ/2). The new `Geometry.segment_angles` averages the polar angle at the segment's middle
  azimuth, weighted by the particle's own elastic law: 1/cos³θ for recoils (exact for any masses) and Rutherford's
  law in the laboratory for the scattered beam, with the setup's masses, since 1/sin⁴(θ/2) holds for a heavy target
  only (it was 0.4° off for ⁹⁸Ru scattered on ²⁴Mg into one ring). The stored angles of `coulex_ni62_c_allmond` were
  0.39° and 0.37° too small (0.001° now); ring by ring the CD of `coulex_ni58` was 0.08° off. The kinematic lines of
  the Data tab use the same angles, for each line's particle. The Doppler correction counts on a segment's azimuth
  only as far as the segment fixes it.
- The Doppler correction takes the speed the emitters have in a thick target (Allmond 2014 verification, second fix
  plan, 14). It took two-body kinematics at the mid-target energy and let the emitter out from the middle; it takes
  both at one depth now, the excitation-weighted mean depth of the forward model per particle detector and kind of
  particle (`gamma_events.reaction_depths`). In `coulex_ni62_c_allmond` the speed was 8.0 % and 4.7 % too high (inner
  and outer ring; 0.2 % and 0.1 % now), the corrected peak 1.35 and 0.77 keV too high (0.08 and 0.05 keV now) and
  7.8 % and 3.8 % wider (rms). Thin targets, targets the generator keeps in one slice, lifetime measurements and the
  correction for the nucleus the `[reaction]` does not name are unchanged, and so are the events.
  - The tails of the fitted peak end: the wider component is a bell that ends (`analysis._biweight`), where a wider
    Gaussian reached on into the background below the peak and took 0.2 % of the area from there once the speed was
    right. Many counts, against the strict full-energy counts of 12 samples (seeds 401–412): −0.02 ± 0.10, +0.10 ±
    0.07 and +0.04 ± 0.05 % (inner ring, outer ring, both); toys of 60, 140, 1400 and 5000 counts within 0.35 %,
    none of 4000 ran away or failed.
- The live counters of a run, its beam time and its sample of particles go together: `RunTaker.step` shows a chunk
  only when all of it is stored. The beam time was shown first, and for 20–30 ms `run_progress()` had progress while
  `run_sample()` had no particles, or columns of different lengths. The runs themselves are unchanged.
- The peak fit of the analysis at few and at many counts (Allmond 2014 verification, round 3, fixes 3 and 4).
  `analysis.fit_peak` fits by Poisson likelihood with a line that is nowhere negative (least squares with the
  variances 1/max(counts, 1) before; a weighted sample keeps its own variances), and `fit_tails` gives the area of the
  fitted peak, the counts of the whole window less the background, where it divided the counts within ±3σ by the
  share of two Gaussians there.
  - A run of too few counts to fix the tails takes their share, offset and width, relative to the core, from a peak of
    15 000 counts simulated with its own setup, gate, rings, Doppler correction and window, and a seed of its own
    (`Analysis.simulated_peak`, kept with the analysis; about 10 s for nine clovers, once per choice of rings). The
    counts fix the area, centre, core width and background. The step, a note and the record say so, and the budget has
    the line "peak tails (simulated)". A closure test cannot test this assumption.
  - A fit whose width or centre ends at its bound (twice the simulated peak's, where there is one), or that finds no
    counts above the background, is `PeakFit.failed`: reported with the reason, without an area and a B(E2).
  - ⁶²Ni on carbon, toys of 60, 140, 1400 and 5000 counts in the three ring choices: the area is within 0.4 % (it was
    1 to 4 % low below 2000 counts, and 7 of 6000 toys of 60 counts ran away); many counts: +0.10, +0.09, +0.09 %
    against the true counts (−0.11, −0.13, −0.14 %). A Gaussian peak of 1400 counts on 20 per bin: 0.0 % (+6 %).
  - The target's misplacement in the budget refits with the same kind of tails, not with whatever the counts of each
    shifted spectrum allow (3.5 % from the scatter of 68 counts, now 0.9 %).
- The beam-energy systematic raises the energy at the entry of the target and slows the beam down from there; it
  scaled the energy at every depth (fix 5): 5.2 % for 0.5 % on ⁶²Ni through 1 mg/cm² of carbon, as a setup with that
  energy gives, where it said 4.3 %.
- The "all orders" step names the feeding from the levels above only when the scheme has some (fix 8).
- Detectors that only touch are no longer reported as overlapping (Allmond 2014 verification): two rings in one
  plane whose radii meet, pads that share an edge, or detectors stacked face to face. Solids overlap when they reach
  into each other by more than 0.05 mm (`scene.TOUCH_MM`, half the step a dragged detector's place is rounded to;
  less for pieces thinner than 0.2 mm). The check is exact along every edge, where it took nine points of each, so
  thin bars that cross are found too. The check for detectors in the beam is unchanged.
- The γ rays of a thick target come from the emitter of their own depth (Allmond 2014 verification):
  `gamma_events.simulate_gammas` takes the undetected emitter's two-body kinematics, and for a single state its
  orientation, at each event's own beam energy where the generator cuts the target into slices. It took those of
  the middle of the target for every event (emitter β 0.026 instead of 0.029 in the inner ring of
  `coulex_ni62_c_allmond`). Targets of one slice and lifetime measurements are unchanged.
- The peak area of the analysis keeps the tails of a Doppler-corrected peak. The depth of the target broadens the
  peak, 2 % of it lies beyond ±3σ of the Gaussian that fits its core, and the straight line beside it takes the
  tails for background: the area was 2–3 % low. `analysis.fit_tails` fits a core with a wider Gaussian on a
  background that steps down across the peak (maximum likelihood) and divides the counts within ±3σ by the share
  of that peak there. It is used where the tails improve the deviance by more than 25 and the nucleus has no other
  line in the fit window; every other peak keeps its Gaussian. `PeakFit.tails` holds it, and the figures draw it.
  - ⁶²Ni on carbon (rings 14–28°, 28–44°, both; twelve samples of 400 000 events): efficiency × correlation against
    the γ rays detected per excitation 0.0 ± 0.4, 0.0 ± 0.2, +0.1 ± 0.2 %; B(E2) simulated then analysed −0.1 ± 0.4,
    +0.1 ± 0.2, −0.1 ± 0.2 % (−1.6, −3.2, −2.7 % with the area of a Gaussian).
- The particle–γ correlation table follows the beam through a thick target (Allmond 2014 verification):
  `gamma.correlation_table` averages its kinematics and angular correlation over the depth as
  `rates.depth_probability` averages the excitation, each depth at its own beam energy and weighted by the
  Coulomb-excitation cross section there, with the emitting nucleus slowed down on its way out. It took the
  kinematics of the middle of the target, with the emitter at its speed just after the reaction. A thin target (a
  beam energy loss below 0.2 %) keeps its factors to 10⁻¹⁵, as do lifetime measurements and the other nucleus of
  the collision; thicker targets move by up to 0.8 % (`coulex_ni58`, `coulex_ni62_pt194`) and 2.5 %
  (`coulex_ni62_c_allmond`) in single crystals.
  - ⁶²Ni on carbon (rings 14–28°, 28–44°, both; twelve samples of 400 000 events): efficiency × correlation against
    the γ rays detected per excitation −0.3 ± 0.4, −0.2 ± 0.2, −0.1 ± 0.2 % (before −0.9, −1.1, −1.0 %); B(E2)
    simulated then analysed +0.9 ± 0.4, −0.1 ± 0.2, +0.1 ± 0.2 % (before +1.5, +0.8, +0.9 %). The closure test
    (`tests/python/test_nuclear_analysis_closure.py`) now allows the Monte Carlo error alone, not 1 % more, and
    averages five samples of 400 000 events instead of taking one of 2 000 000 (3 GB of memory instead of 10).
- The analysis reads B(E2) against the forward model (Allmond 2014 verification, fixes 1 analysis side, 2 and 7).
  Simulate-then-analyse on 62Ni at 111.6 MeV on 1.04 mg/cm² of carbon read B(E2) 35 % high: ⟨P⟩/B took P at
  the mid-target energy and at the CM angles of the elastic particles, while the forward model excites through
  the depth and gives the excited carbon recoil its own CM angle and lab-to-CM Jacobian (4 % larger here).
  - `Analysis._mean_probability` now takes the excitations expected per particle from the forward model's own
    cross sections (`rates.reaction_at`, as `Rates` integrates them): at each particle's lab angle, the
    Coulomb-excitation cross section over that of every particle reaching the angle, each integrated over the
    depth of the target at the beam energy there and through its own kinematics, then averaged over the particles
    in the rings. A change of the forward cross section (its Rutherford factor, say) reaches the analysis by
    itself. With a level scheme, first order still uses the [reaction] state's first-order P inside the same
    cross section, and the all-orders step scales the scheme's mid-target yields by the forward model's
    excitations per particle over the solution's. The beam-energy systematic raises the beam energy at every
    depth.
  - Each crystal's full-energy efficiency is taken at the Doppler-shifted energy it sees the line at (the measured
    over the corrected energy of the gated γ rays), not at the rest energy: +0.6 % in the efficiency here, where
    most clovers stand backward of 90°.
  - The inelastic particle gate warns (`Result.gate_warning`, also the first note and in the gate step) when the
    elastic and inelastic groups overlap through the depth of the target for most of the excitations in the
    rings, and says to use the gate "all"; it does not switch. `Result.gate` and `Result.gate_overlap` say which
    gate was used and how much the groups overlap. The groups overlap in this sense in most shipped examples too
    (the depth spread is wider than E*'s shift of the energy).
  - Closure on the 62Ni on carbon setup with gate "all" (ring 14–28°, ring 28–44°, both): +1.9/+0.6/+0.9 % with
    4 000 000 weighted events (Monte Carlo 0.4/0.25/0.2 %), +0.9/+0.0/+0.2 % with 2 000 000 and another seed;
    before, +36/+35/+35 % on the wave-1 forward model (+9.7/+8.5/+8.8 % before wave 1). The residual of about
    +0.5 to +1 % is 0.3 % in ⟨P⟩/B against the events' own excited share and 0.3–0.6 % in the angular-correlation
    table, which takes the kinematics of the middle of the target (`gamma.correlation_table`).
    `tests/python/test_nuclear_analysis_closure.py` checks it.
- Coulomb excitation has one Rutherford convention (backlog 58): Alder and Winther's dσ = P (v_f/v_i)(a²/4)/
  sin⁴(θ/2) with the symmetrised orbit, the elastic Rutherford cross section at the initial E_cm times v_i/v_f, as
  in GOSIA. Before, the rates and the event generator used the elastic Rutherford cross section and
  `Coulex.rutherford_cm` the symmetrised a²/4 without v_f/v_i (8.7 % apart for ⁶²Ni on carbon at mid-target).
  - `Coulex.rutherford_cm` and `cross_section_cm` (and `total`, `cross_section_lab`), the rates
    (`rates.rutherford_factor`), the generator's weights (`DepthProbability.slices(..., ratio=True)`), the Doppler
    table, and the all-orders yields of `multistep` and `orientation` all use it.
  - Ratios say which Rutherford they divide by: `Coulex.ratio_to_elastic` and the new `DepthProbability.ratio`
    divide by the elastic one, as an experiment normalises to the elastic singles; `Coulex.elastic_rutherford_cm`
    is that cross section. GOSIA's ratios divide by its own inelastic integrated Rutherford (11–16 % larger on the
    ⁶²Ni rings).
  - Allmond 2014 (⁶²Ni on carbon): the app's integrated Rutherford factor equals GOSIA's to 0.03 %; σ_C/σ_R per
    ring is now +1.3 % and +1.2 % against GOSIA's 1.0080e-3 and 7.699e-4 (was −2.4 % and −2.6 %). The rest is
    GOSIA's integer mass numbers and its E1 polarisation: rerun on the app's input with the nuclear masses and
    DIP off, GOSIA agrees within 0.2 % for both the elastic and the inelastic ratio. `test_nuclear_thick_target`
    tightens to that.
  - The shipped examples' excitation rates rise by v_i/v_f: +3.4 % (`coulex_ni58`), +0.6 % (`coulex_ru98_plunger`),
    +0.3 % (`coulex_ni62_pt194`, `_clovers`); the all-orders yields fall by v_f/v_i from the old symmetrised factor.
- Thick-target Coulomb excitation: the excitation probability P(θ) is now integrated over the depth of the layer,
  each depth at its own beam energy and weighted by the Rutherford cross section there, as the Rutherford factor
  already was. Before, P was taken at the mid-layer energy alone, which for 62Ni at 111.6 MeV on 1.04 mg/cm² of
  carbon (Allmond et al. 2014, a 38 % energy loss, P changing 150-fold) put σ_C/σ_R 21–23 % low.
  - `rates.DepthProbability` (from `rates.depth_probability`) solves P at a few energies through the layer (up to
    8 Chebyshev points, by the energy loss) and interpolates in ln P; a layer the beam loses under 0.2 % of its
    energy in keeps the mid-layer P exactly. The rates take P at each depth point; the event generator cuts such a
    layer into up to 16 slices, each with P averaged over it; the Doppler table weights each depth with its own P.
    `DepthProbability.probability(theta_cm)` is the depth-averaged P at given CM angles, for an analysis.
  - The 62Ni on carbon setup now gives 9.904e-4 (ring 14–28°) and 7.530e-4 (28–44°) in first order, and
    9.835e-4 and 7.500e-4 to all orders: 2.4 % and 2.6 % below GOSIA (1.0080e-3, 7.699e-4). The remainder is the
    Rutherford convention (backlog 58), left to its own fix. `coulex_ni58` (a 7 % energy loss) moves by +0.6 %,
    `coulex_ru98_plunger` by −0.5 %, `coulex_ni62_pt194` by −0.01 %: the depth average itself, which a brute-force
    average over the depth reproduces.
- A single 0⁺ → 2⁺ transition with a diagonal matrix element (Q_s) in the beam's or target's `[levels]` is solved
  with the coupled equations (to all orders, unless `[reaction] order = "first"`). Before, it fell back to the
  first-order one-γ-ray chain and Q_s was ignored; now Q_s = ±0.1 eb moves the 62Ni ratios by +2.9/−2.8 % and
  +2.0/−2.0 %, as GOSIA's +2.7/−2.7 % and +1.9/−1.9 %.
- The GOSIA input writer (`Multistep.gosia_input`) is rewritten and checked against runs of GOSIA 20110524
  (Allmond 2014 verification, fix 5).
  - The sign of the partner Z was backwards for both roles. In GOSIA a negative partner Z means the projectile is
    the investigated nucleus. GOSIA's source says so, and its runs print "PROJECTILE EXCITATION"; a comment in its
    own CMLAB states the rule the other way round. Target excitation (¹⁶O on ⁵⁸Ni) used to stop in GOSIA with
    "maximum scattering angle 14.06°".
  - One experiment per particle detector and per nucleus that reaches it. A detected recoil gets a negative angle
    (the beam's lab angle, from GOSIA's INVKIN), and OP,INTI's mesh is in recoil angles, not OP,INTG's
    two-valued beam angles. Side detectors get IAX = 1 and the azimuth they cover at each angle.
  - The yields are integrated over the target with OP,INTI, using the app's stopping powers for the beam. CONT
    writes `SPL,1.`: the default single polynomial through the stopping table puts a thick target's Rutherford
    17 % low. CONT also writes `DIP,3.8` and the blank line GOSIA needs after `END,`.
  - OP,YIEL is written in GOSIA's real format. The γ detectors come from OP,GDET, in the same file, or from a 4π
    detector (`gamma="4pi"` with `GOSIA_4PI_FORT9` as `fort.9`). The old OP,INTI lines of "measured yields"
    are gone.
  - `physim.external.gosia.integrated` reads the integrated Rutherford and yields back. Ratios built from them
    divide by GOSIA's inelastic Rutherford, 11–16 % above the elastic one for ⁶²Ni on carbon.
  - Checked against the outside-tools GOSIA run for ⁶²Ni on carbon (carbon recoils at 14–28° and 28–44°):
    - γ yield over GOSIA's Rutherford: 9.070e-4 against 9.076e-4 (ring 2), 6.642e-4 against 6.647e-4 (ring 3);
    - over the elastic Rutherford: 1.0077e-3 against 1.0080e-3, 7.696e-4 against 7.699e-4.
- Three misleading outputs found while checking the app against Allmond et al. 2014 (⁶²Ni Coulomb excitation):
  - The lifetime planner no longer names a method when the recoil does not fly forward out of the target (β = 0,
    as for ⁶²Ni backscattered off ¹⁹⁴Pt in `coulex_ni62_pt194`): each level's method is "none", and the warning
    says plainly that no Doppler-shift method (RDDS or DSAM) fits this geometry and why. Before, it named "rdds"
    beside a β of zero.
  - In DSAM mode (the setup's DSAM, or a recoil that stops in the target) the lifetime plan lists no plunger
    distances, sensitive range or counts per distance, even for levels whose τ would suit a plunger.
  - Shell model: an effective charge for a particle type with no valence particles in the chosen space (e_p for
    ⁶²Ni in JUN45, ⁵⁶Ni core) is warned about, since it has no effect: `shellmodel_spaces.charge_warnings`, the
    estimate's `warnings`, and `warnings` in the MCP tools `shell_model_options`, `estimate_shell_model`,
    `run_shell_model` and `shell_model_results`.
- Editing a fusion setup with an ancillary (neutron wall, charged-particle ball, …) could break the app's setup
  panel: `workbench.card_setup` raised for ancillary cards while a change was pending.
- A shell-model job reads as finished, with its finish time, as soon as its results are written; before, a status
  read between the results and the worker's last status write saw it still running (a rare test failure).
- The alignment fit (`alignment.fit_offset`) could come out short of the offset at small statistics: with one
  ring per point and 600 000 events of `coulex_ni58` (about 30 γ rays per point), a 2 mm offset gave −1.33 ±
  0.25 mm over 16 independent samples without refinement (−1.0 ± 0.4 mm with the default reference of 400 000
  events), because each point is made again from the same γ rays at every step of the scan and the noise part
  of χ² changes with the geometry assumed (it moves the minimum of the mean χ² by 0.3 mm). It now merges
  adjacent rings until the median point holds 60 γ rays (`rings_per_point` in the result): −1.87 ± 0.22 mm over
  18 samples, the uncertainty unchanged (1.1 mm), and −0.07 ± 0.27 mm with the target where it is assumed. The
  test of the 2 mm offset, marked as expected to fail, did not fail for a bias alone: at 600 000 events the
  uncertainty is 1.1 mm, so a test within 1.5 of them fails one sample in seven, and the old seed's sample was
  such a one. It now uses 1.8 million events (uncertainty 0.6 mm; −1.98 ± 0.19 mm over 12 samples) and checks
  within 2.5 uncertainties.
- The test of a carbon contaminant's elastic line passed only by its seed: in DSSD-R's forward strips the
  median energy of ¹⁶O scattered on ¹²C sits 0.3 MeV above the line, because the line is drawn at the strip's
  middle while the Rutherford cross section falls across the 4° strip, where the energy falls by 0.65 MeV per
  degree. The energy loss and the generator agree with the line to 20 keV. The test now measures each particle
  against the line at its own angle (`dataviews.kinematic_lines` gives each ring's angle as `"theta"`), with
  medians and median absolute deviations, an allowance about ten times tighter than before; it passes on seeds
  1–10.
- Stopping of slow heavy ions in solids: the effective charge uses ZBL's screening length, the stopping below the
  speed where its fit holds is proportional to the speed (ZBL's rule), and the Lindhard–Scharff floor is scaled by
  0.8. Against 3620 measured points (Ar to U in 12 elemental solids, IAEA database, the new
  `stopping_measured_heavy.csv` and check `stopping_heavy_vs_data`) the median for I–U goes from +19% to +7% below
  0.1 MeV/u and from −16% to −10% at 0.1–1 MeV/u; in light targets (C–Si) below 0.1 MeV/u Ar–Ag ions are now about
  20–30% high (see the register).
- The plan's particle–γ coincidences took the angular correlation of the `[reaction]` state alone, in first
  order, while the simulated events follow the whole level scheme to all orders; and the table was sampled too
  coarsely (4 points each way over a particle detector, 9 over a crystal face: single factors up to 9% off).
  `gamma.correlation_table` now takes the state from the level scheme's solution when the events follow one, as
  the table of the other nucleus does, and uses 8 points each way and 36 per face (within 0.1% of four times as
  many; a quarter of a second, kept with the solution). On the ⁶²Ni on ¹⁹⁴Pt template the state alone sent 14%
  too few γ rays to the crystal at 90° with the CD and 3% too many to those at 45° and 135°: the four crystals
  together happened to cancel, the four clovers near 45° and 135° did not. The plan's nickel coincidences go
  down by 1.7% (CD) and 2.7% (strip detectors) for the four-clover template, and by 0.4% (CD), 0.7% and 2.0%
  (the two strip detectors, mirror images that now agree) for the single crystals. Runs against the plan: nickel
  0.994 ± 0.004 with clovers over an hour of beam (0.973 ± 0.006 before) and 0.999 ± 0.006 with single
  crystals. The plan still counts one γ ray of the state per excitation; in the scheme 0.5% of the excitations
  decay past it. A setup without a level scheme changes only by the finer sampling (⁵⁸Ni template: 0.4% at most
  in a particle detector's coincidences), and so does the platinum table (0.6% at most).
- A run with clovers counted 3% fewer γ rays than the plan and the clover's efficiency say. The faces of a
  clover's four crystals overlap in the model (discs 50 mm across, 45 mm apart), the plan counts each crystal's
  whole face with the data sheet's efficiency, and the simulated events gave a γ ray through the overlap to one
  crystal only. The events now keep the data sheet's efficiency: a crystal's efficiency belongs to the part of
  its face that is its own, 96.9% of it (`gamma_events._own_shares`). The platinum peak of the four-clover
  template is 0.995 ± 0.001 of the plan over an hour of beam (0.97 before). The plan, the response and
  single crystals are unchanged.
- Without a run, B(E2) read up to 23% high on the ⁶²Ni on ¹⁹⁴Pt template (963 to 1090 e²fm⁴ for 884). The
  simulated sample is a few thousand weighted entries scaled to the counts of the beam time, and the peak fit took
  those counts for Poisson counts: it trusted every bin a hundred times too much and ran to its widest peak
  (83 keV for a peak 16 keV wide). The fit now takes the variance of the entries (`analysis.fit_peak(variance=)`).
- The note on a large excitation probability said "first order is not enough … use GOSIA" for a setup whose
  rates and events already follow its level scheme to all orders. It now says that they do, and that the analysis
  still extracts B(E2) with first-order theory; without a level scheme it says to give one, or to use GOSIA.
- γ-ray spectra reach the top of the level scheme. The predicted spectra and the singles stopped at 1.3 × the
  transition energy (1525 keV for ⁶²Ni's 1173 keV), short of the 2302 keV γ ray of its second 2⁺ state. They now
  go to 1.1 × the highest γ-ray energy of the scheme (2532 keV there) at the same width of bin
  (`gamma_events.spectrum_top`); the Data tab and the report use the same top.
- The event generator no longer spends nearly all its events where no detector counts anything. With ⁶²Ni on
  ¹⁹⁴Pt, strip detectors just behind 90° and the target 2 mm upstream, 100,000 events gave 1,400 counted particles
  and 60 excited ones, against 32,000 and 16,000 with the target in place; now about 37,000 and 19,000 wherever
  the target is.
  - The CM angle was drawn over one range, from the smallest to the largest lab angle of any detector. With
    detectors on both sides of 90°, or one reaching just forward of it, that range held the slow recoils of
    distant collisions near 90°, where Rutherford's cross section is hundreds of times larger. Each channel is
    now sampled over the separate ranges its detectors see (`rates.cm_ranges`), without the angles between them.
  - When a channel has more than one range, a short pilot run finds the fraction of each range's particles that
    are counted, and the events are shared out by the rate counted from each range (`events.share_out`). A range
    in which nothing is counted keeps a few percent of its channel's events, so no rate is dropped: the ¹⁹⁴Pt
    recoils on that strip detector's edge are still in `Events.rate(counted=False)`, as in `Rates`, as a few
    rows of large weight.
  - The rates and spectra are the same within their statistical errors, which are smaller: the three examples
    give 1.7 to 2.7 times as many counted particles per event. A seed gives other events than before.
  - `generator_config` lists one entry per range, with `channel_of` naming the channel of each
    (`events.channel_index`); an unweighted run (`runs.BeamRunner`) that is extended after this change makes its
    pilot run again.
- "Check the alignment" with the fit no longer runs out of memory on a run: the reference simulation was asked
  for as many events as the run had reactions (63 million for an hour of ⁶²Ni on ¹⁹⁴Pt, over 14 GB). It is now
  made of simulations of at most 400,000 events each (`alignment.REFERENCE_EVENTS`), one after the other and four
  at most, until it holds as many γ rays as the run, and the check peaks at the 3 GB the run itself takes.
- The fitted target offset is right for a fast emitter (on that run: −2.4 ± 0.3 mm where −2 was put in; before,
  −3.2 ± 0.2 from a jagged scan with χ² six times the number of points).
  - Each centroid's window follows the corrected peak as it stands in the crystal (`alignment.peak_widths`).
    The window was three times the detector's resolution, which suits ⁵⁸Ni recoiling at 3% of c; ⁶²Ni at 9% of c
    leaves a corrected peak of 20 keV FWHM in a crystal of 35 mm radius at 120 mm, and the window cut through
    it. A point needs ten γ rays inside the window, not ten in all.
  - When the fit finds the target more than two uncertainties from where it was assumed, the reference is made
    again with the target where the fit put it (twice at most): left at a geometry 5 mm off, the fit overshot
    by about 1 mm. The uncertainty is scaled by √(χ²/ν) when the centroids scatter by more than their errors.
  - A reference simulation that gives too few γ rays for a pattern is no longer taken for one in silence: the
    fit says it judged against flat lines (`"against"`), and a parabola whose minimum lies outside its own
    points is not used.
  - The overlay's widths and range follow the peak too (it reported a FWHM of 4 keV for that 17 keV peak).
  - The fit takes about 15 s on a planning sample as before, and two minutes on that run.
- The analysis now reaches the reports. Before, it lived only in the app's memory, so the record of the experiment
  said "no analysis" after any edit of the setup, after reopening the experiment, or once a source run was the
  current run.
  - Each analysis is kept in its run's folder (`analysis.json`; `Planner.run_analysis`), with the spectrum it
    fitted; the alignment check is kept with the run (`alignment.json`) and the multi-step solution with the
    experiment (`multistep.json`, used while the setup's physics is unchanged).
  - The record of the experiment shows, for every analysed run: the fitted peak drawn on the Doppler-corrected
    spectrum, the chain step by step, the uncertainty budget, the shape (β₂, Q₀, Weisskopf units), the alignment
    check and, once solved, the multi-step table and the prolate–oblate comparison (`logbook.analysis_html`).
  - Level schemes are drawn in the record and in the beam-time report (`logbook.level_scheme_png`): the setup's
    schemes from ENSDF, or the ground state and the excited state when none has been looked up.
  - The beam-time report's zip also holds `experiment.html` once the experiment has runs.
- The experiment's report now carries every run's spectra: a panel per particle detector and, per γ-ray detector,
  the spectrum as measured and Doppler-corrected for the emitting nucleus (the correction and the gate the
  analysis used), and a source run's crystals; before, no γ-ray spectrum reached any report
  (`logbook.run_spectra_png`, `run_views`).
- The Data tab's γ-ray spectra no longer stop at 1.3 × the transition energy (0.6 to 1.3 × E_γ in the app, which
  cut off every line above it, such as ⁴⁰K at 1461 keV): the histogram now runs from 0 to the highest of
  1.3 × E_γ, the singles spectrum's top and the highest energy recorded (`dataviews.gamma_range`), and "γ view
  from / to" zooms the panels without changing the data or the CSV.
- A particle detector hidden behind another no longer reports counts in the analytic rates, peaks and beam times.
  - The integral leaves out directions that meet another detector first, as the event generator always did.
  - A partly hidden detector is integrated on a finer grid, so the shadow's edge is resolved within 0.5%.
  - A fully hidden detector's warning names the detector in front, and shadowing is now a warning, not a note.
- The app could show one setup's numbers for another: a result still being computed on the worker thread when a
  field was edited was kept for the new setup (the old rates under the new beam energy, until the next edit).
  Each setup's results now have a place of their own, and a computation works on a snapshot of the planner
  (`Planner.snapshot`), so a late result stays with the setup it was computed for. That holds for the buttons
  too (Analyse, Check the alignment, Solve with all orders, the simulations): a multi-step solution that ends
  after a change was applied is kept in `multistep.json` under the setup it was solved for, not the new one.
- A tab left while the numbers were being computed kept its old content when you came back to it. Every tab is
  marked stale by a change, the one on screen too.

### Changed
- The planner app has a new look and layout ("control room").
  - The steps stand in a rail on the left, from the top down, with the experiment's name above them and
    **Physics** and the theme at its foot. There is no header any more. The step on screen has its number on the
    accent colour; a step with its result has a green tick in place of the number.
  - The setup is edited in a window of its own. **Edit the setup** on the Setup step opens it at the right of
    the page, with the scene in view beside it; **Done editing** applies what was changed and closes it, and
    **Cancel** closes it without the changes. The page shows the rail and the step only.
  - Above the scene a strip gives the setup in a few numbers: particles counted, how many with a γ ray, the
    solid angle and the γ-ray efficiency. Errors and warnings of the setup are written out under it; notes
    are a count that opens the window.
  - The Analysis step leads with the result (B(E2), its uncertainty, the value put in, the beam time). The
    spectra and the fitted peaks are panels on the left, with one table of all the peaks found. On the right, the
    chain from the area to B(E2) is four large numbers with what is done between them; the smaller facts are a
    quiet list below, and every step in words is behind a fold. The systematic uncertainties are bars, the largest
    first.
  - Each step uses the width of the window; what was stacked stands side by side where there is room. Setup: the
    tables of angles are beside the scene, under the selection. Plan: the three tables are panels in one row. Run:
    taking a run and the current run's counts on the left, the runs taken on the right. Data: the particle
    spectra share one row and the γ-ray spectra the next, and the two 2D views are side by side. Analysis: the
    comparison beside the choices, the systematics and the shape beside the fitted peaks, the alignment check
    beside the all-orders solution.
  - Colours and controls: a dark page of near-black panels with thin borders and corners of 6 px, numbers in the
    monospace face, an amber accent with dark text on filled buttons, and a light theme of the same design. On
    the dark page the curves of the figures are recoloured to read on it (`figures.DARK_COLORS`).
- One particle per collision. A collision either excites the `[reaction]` nucleus or leaves it as it was, so
  elastic scattering is the Rutherford cross section times 1 − P(θ), in the rates (`rates.reaction_at`) and in
  the events. Elastic and excited particles together are Rutherford scattering; before, the excited particles
  were counted on top of all of it.
  - The particle rates of a Coulomb-excitation setup are lower where the excitation is strong: the CD of the
    ⁶²Ni on ¹⁹⁴Pt template counts 2741 particles per second, not 3292. The excitations, the coincidences and the
    beam time needed are unchanged.
  - The analysis normalises to the particles in the same rings, as a detector counts them. A run taken before
    (no `one_particle` in its summary) is still normalised to its elastic particles, as it was taken.
- Only the γ rays that reach a crystal are drawn (`simulate_gammas(aimed=)`), by a thinning that leaves their
  number and directions as the distribution has them: the γ rays of a setup and seed are other γ rays than
  before, with the same statistics, in under half the time.
- The app is laid out as a workflow. The seven equal tabs are a numbered step bar, **1 Setup → 2 Plan → 3 Run →
  4 Data → 5 Analysis → 6 Report**, that stays at the top of the window (`workbench.STEPS`, `workbench.steps`).
  - Each step has a line of status under its name (the detectors, the beam time needed, the runs, the current
    run and its gates, the B(E2) of the analysed run, what the report holds), a tick when it has its result, and
    the warning colour when it needs attention. A step that needs a run first says so.
  - Each step ends with one button that leads on (`workbench.next_step`): *Next: check the plan*, *Next: take a
    run*, *Next: look at run 2's data*, *Next: analyse run 2*, *Next: the report*.
  - The setup panel stands beside the Setup step only; Plan, Run, Data, Analysis and Report have the whole width.
    The status strip is gone: the experiment, beam and target are in the header, and the detectors and warnings
    under **1 Setup**.
  - Physics is a button at the top right, not a tab: it is a reference, not a step.
  - The experiment's name in the header is a menu (open, new, import a setup file, export the setup) in place of
    four buttons.
  - The Data step without a run says what it needs; the predicted spectra it used to show are a check of the Plan
    (*Predicted spectra: what each detector will record*).
  - The Analysis step opens with a list that jumps to its four parts, and says when it works on a simulated
    sample because there is no run.
- There is one report. The Report tab showed three (the record of the experiment, the beam-time report of the
  setup and the run record), and the spectra of the runs were small panels in one of them. They are now the nine
  sections of a single document (`logbook.experiment_html`): summary; setup; plan; runs and their spectra; gates;
  the Doppler correction; analysis; every number explained; data, method and what is left out.
  - Each beam run is drawn in full: the energy in each particle detector (with the analysis's gate filled in),
    the energy against ring with every gate's window drawn on it (`logbook.run_rings_png`), and every γ-ray
    detector, and all of them summed, as measured and Doppler-corrected with the γ rays of the level scheme
    marked (`run_gamma_png`).
  - The gates have their own section: the corrected γ-ray spectrum under each gate, and a table of the particles
    that pass and the γ rays in coincidence (`run_gates_png`).
  - The Doppler correction has its own section (`doppler_html`): the nucleus corrected for and why, the formula,
    where β and the angle come from event by event, the crystals as the correction sees them, the speed and
    shift for every pair of particle and γ-ray detector, and each run's settings and corrected line widths.
  - The report is saved as a PDF as well as HTML (`logbook.experiment_pdf`, `html_to_pdf`): Edge, Chrome or
    Chromium prints it without a window. Without one of them the app says so; the HTML prints to the same pages.
  - The zip holds `report.html`, `report.pdf`, the setup, the plan's tables under `tables/` and figures under
    `figures/`, and the experiment folder. `experiment.html` and `record.html` are gone from it, and
    `Report.write` writes the one report as `report.html` (`--pdf` adds the PDF). `Report.to_html()` and
    `record_html()` still give the plan and the record alone, from Python.
  - Where the level scheme has levels the simulated events never populate, the report names them: with
    first-order events (`[reaction] order = "first"`) for ⁶²Ni the 0⁺ at 2049 keV and the 4⁺ at 2336 keV.
- The simulated events of a setup with a level scheme populate its levels as the coupled equations give them, to
  all orders, where they excited each level in one step from the ground state. A run of ⁶²Ni (237 MeV) on ¹⁹⁴Pt now
  holds the γ rays of the 4⁺ state (4⁺ → 2⁺, 1163.5 keV, 3% of the 2⁺ → 0⁺ line) and of the second 0⁺ state
  (0⁺ → 2⁺, 875.6 keV, 4%), which first order never reaches; in a 40 min run their counts agree with the multi-step
  yields (`Multistep.yields`) within statistics.
  - The level an excitation starts in is drawn from the populations the collision leaves at its CM angle
    (`gamma.cascade_excitation` is now a `CoupledChannels`, kept per setup). Each γ ray's angular distribution
    comes from the same solution: the 4⁺ state has the alignment its two-step excitation gives it.
  - The excitation probability of the event generator and of `Rates` is the sum of those populations
    (`rates.excitation_probability`), so the two agree as before; in the ⁶²Ni example 3% more is excited than in
    first order. A run's "coincidences in the peak" counts the `[reaction]` state's own γ ray; it counted the
    full-energy γ rays of every transition.
  - `[reaction] order = "first"` keeps first-order events, bit for bit those of before. A setup without a level
    scheme (`coulex_ni58`) is unchanged.
  - The table of 90 angles is solved with all angles together (`CoupledChannels.solve_many`), the forward
    angles where less than 10⁻⁶ is excited staying in first order: 2 s for the six levels of ⁶²Ni, against 16 s
    one angle after the other (first order: 0.2 s). It is made once per setup, beam energy and target.
  - In the corrected spectrum of the ⁶²Ni example the 4⁺ → 2⁺ line does not stand apart: it lies 9.4 keV below
    the 2⁺ line, whose Doppler-corrected width is 15 to 20 keV with crystals of 35 mm radius at 120 mm. The
    0⁺ → 2⁺ line is clear of it.
  - The analysis still extracts B(E2) with first-order theory. On the ⁶²Ni example it returns 5% more on events
    with all orders than on first-order events (3% from the feeding of the 2⁺ state, the rest from the 4⁺ → 2⁺
    line inside the fitted peak). With the ENSDF scheme added to `coulex_ni58`, the 2⁺ yield falls by 5% in the
    CD and 2% in the strip detectors (the reorientation, mostly).
  - Still left out: every excitation has the kinematics of the `[reaction]` state, and elastic scattering is not
    reduced by the excitation probability (`record.LEFT_OUT`, `workbench.NOT_SIMULATED`,
    `docs/theory/coulex.md`).
- The app opens on Coulomb excitation of the beam at backward angles when there is no experiment yet: ⁶²Ni (6⁺,
  237 MeV) on ¹⁹⁴Pt, the first template (`coulex-beam-backward`, example `coulex_ni62_pt194`). Before, it opened
  on ¹⁶O on ⁵⁸Ni, which is still the second template.
  - The annular silicon detector stands 30 mm behind the target (126° to 163°; the beam comes in through its
    hole), and the two strip detectors at 110° (96° to 124°), ending where the annular one begins. No particle
    detector reaches forward of 90°, where the scattered beam and the slow ¹⁹⁴Pt recoils of distant collisions
    arrive. The γ-ray detectors are those of the ⁵⁸Ni template.
  - The reaction excites the projectile's first 2⁺ state (1172.91 keV). The setup carries the level schemes of
    both nuclei from ENSDF: ⁶²Ni with its 4⁺ state at 2336 keV, and ¹⁹⁴Pt with its 2⁺ at 328 keV, for the
    multi-step solution of either.
- The planner app applies the changes of the setup when you say so. Before, every field recomputed everything as
  soon as you left it, so five changes meant five rounds of a second or more each, with the numbers of four
  setups nobody wanted in between.
  - A change now goes into a draft and waits: its field is outlined, its card says *not applied*, and a bar at
    the top of the setup panel lists the changes (*Beam: energy 30 MeV → 31 MeV*) with **Apply** (Ctrl+Enter) and
    **Discard**. Adding, copying, removing and dragging a detector wait in the same way. An edit takes 0.1 s and
    computes nothing; moving from field to field with Tab keeps the cursor, since the panel is updated where it
    stands instead of being rebuilt.
  - What needs no rates follows the draft at once: each card's first line, the problems of a value, the scene,
    the tables of angles, and the quick numbers of the selected detector. Everything else shows the applied
    setup, and a note above the tabs says so.
  - A second after the last change the draft's rates are worked out in the background, so Apply usually shows
    its results in 0.2 to 0.3 s; pressed at once it takes what the rates take (1.5 s on the ⁵⁸Ni example).
  - An edit of the beam current or of the run conditions (beam time, counts wanted, dead time, coincidence
    window) reuses the rate rows instead of integrating again (`Rates`, 5 ms instead of 0.9 s): the run
    conditions do not enter them and the current scales them. Apply then takes 0.4 to 0.5 s.
  - In an experiment the setup is saved on Apply, not after every field; a run is taken of the applied setup.
    Level-scheme actions with a button of their own (look up, remove, a matrix element, **Plan**) and the
    emission pattern apply at once. The switch **Apply each change as I make it** brings the old behaviour back.
  - `Planner.adopt` takes another planner's setup and the results it has for it (Apply: the page's planner adopts
    the draft's; Discard: the reverse); `workbench.changes` lists what differs between two setups,
    `workbench.draft_decisions` gives a draft's cards without computing anything.
- The planner app no longer stops answering while it computes. Before, the page function, every tab switch and
  every rebuild after an edit computed on the server's event loop: the page froze for one to five seconds at a
  time, and for twelve seconds with the connection dropped ("Connection lost") on the Data tab of a long run.
  Measured in a headless browser on the ⁵⁸Ni example and on an experiment whose current run holds an hour of
  events (128 MB), before → after:

  | | before | after |
  |---|---|---|
  | The page opens (first byte) | 2.2 to 5.1 s, server busy throughout | 0.03 to 0.08 s, then a spinner for 1.2 to 1.5 s |
  | A change, Setup tab on screen, to the updated tab | 2.1 to 2.9 s (frozen 1.1 s) | Apply at once 0.4 to 1.5 s, after a pause 0.3 s (longest pause 0.1 s) |
  | A change, Plan tab on screen | 4.2 to 4.3 s (frozen 2.0 s) | Apply at once 1.4 s, after a pause 0.2 s (0.1 s) |
  | A change, Physics tab on screen | 3.9 s (frozen 1.6 s) | Apply at once 2.3 to 2.7 s (0.2 s) |
  | Plan, first shown after a change | 2.3 to 3.5 s (frozen 1.3 to 2.4 s) | 0.4 s: built ahead |
  | Data tab of the long run, first shown | 12.3 s frozen, connection lost | a spinner for 7 to 9 s if chosen at once (the page answers), 0.4 s once loaded ahead |
  | A choice on the Data tab (the Doppler correction) | 3.3 s (frozen 2.8 s) | 0.7 s (0.09 s) |
  | Another run made current, then Data | 1.2 s + 3.7 s (frozen 0.9 + 3.1 s) | 0.2 to 0.3 s + 0.4 s once its events are loaded |

  - Everything a tab shows is computed on a worker thread and kept with the setup's results; the tab is then
    built from what is kept, which takes the event loop tens of milliseconds. After a change the numbers come
    first, then the setup panel, the strip and the tab on screen, then the other tabs in the order one looks at
    them, so a tab is there when it is chosen (a spinner in the tab if it is chosen sooner). The current run's
    events are loaded last. A tab that is up to date is shown as it is.
  - The page is served before anything is computed: a spinner stands in the setup panel and the tab for the
    second the rates take.
  - The Plan tab's checks and the setup panel's cards are built when they are opened. The checks you left open
    stay open when the tab is rebuilt after an edit.
  - The Data tab keeps each view by what it depends on: changing the Doppler correction recomputes the γ-ray
    spectra and leaves the particle spectra and the 2D views (2.4 of the 2.8 s) as they are.
  - Kept instead of computed again: the polar-angle range of a detector face (`Geometry.theta_range`, 25 ms a
    call, asked for forty times to build one page: the warnings, each card and each "How to read this" took 0.1
    to 0.4 s each time they were shown, and the rate integrals 2.7 s instead of 1.4 s), the rate rows of a
    detector (`Rates._select`), the setup's warnings (`Planner.warnings`), the cards' consequence lines, and each
    figure in the page's theme, ready to send (a figure was themed and serialised three times per element).
  - A run's events load in 5.2 s instead of 8.6 s for an hour of beam (on their own; 7 to 9 s beside the page): the particle and the γ-ray columns are
    each inflated only when asked for (the file was inflated whole, twice), and each row's segment centre comes
    from a table over the segments that were hit instead of a sort of all rows.
- The Coulomb-excitation orbit integrals are 20 to 600 times faster at the same accuracy (10⁻⁹): the half-orbit
  symmetry, panels of four radians of phase, and the oscillating tail by its endpoint expansion instead of
  quadrature out to where it dies. A P(θ) table takes 50 ms instead of 2 s (an E1 table, 0.06 s instead of
  30 s), so the planner's rates and the Doppler table after an edit take a fraction of what they did.
- The rates transform the kinematics for all beam energies through the target at once
  (`kinematics.lab_points`, the one-body formulas broadcast) instead of once per energy and branch: a few
  hundred fewer calls per detector, the same numbers.
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
- A template with clovers: `coulex_ni62_pt194_clovers` (in the app, *Coulomb excitation of the beam (62Ni on
  194Pt) with silicon at backward angles and four clovers*). It is the ⁶²Ni on ¹⁹⁴Pt template with four germanium
  clovers in place of its four single crystals: above and below the beam at 45° and 135°, 150 mm from the target,
  with add-back and BGO shields. The γ-ray efficiency at 1173 keV is 1.71% against 0.83%, and the beam time for
  the counts wanted 3.7 min against 9.8 min.
- Both nuclei of the collision are in the simulated events. The `[reaction]` names one; with a level scheme of
  the other in `[levels]` (¹⁹⁴Pt in the ⁶²Ni on ¹⁹⁴Pt template), that nucleus is excited too and its γ rays are
  in the same events and runs.
  - It is excited on the orbits of the scattered particles with the probability its coupled equations give
    (`gamma.other_excitation`), independently of the first nucleus (a collision can excite both), and decays
    by its scheme's cascade. The γ rays carry a `nucleus`
    column (0 the `[reaction]` nucleus, 1 the other), which runs store. The γ rays of the `[reaction]` nucleus
    are unchanged, one for one.
  - Each nucleus's lines are sharp with its own Doppler correction and smeared by the other: ¹⁹⁴Pt at 328.5,
    293.6, 482.8 and 622.0 keV with the recoil correction, ⁶²Ni at 1172.9 keV with the projectile correction.
  - Its γ rays count in the singles, the random coincidences and the dead time (`gamma_events.other_lines`).
  - The plan gives its first peak's coincidences per particle detector (`Rates.other_peak`, 59 per second in
    the CD of the template against 5 for the nickel peak), and predicted against measured has a row for it.
    The counts wanted and the beam time stay those of the `[reaction]` nucleus.
  - The template's ¹⁹⁴Pt scheme holds the seven positive-parity levels up to the 6⁺ state at 1412 keV (it
    stopped at the 4⁺ state at 811 keV, where the excitation piled up).
  - Five minutes of beam on the template take 13 s (7 s without the platinum).
  - Runs taken before hold the `[reaction]` nucleus only and are marked *one nucleus only* in the run list.
    **Take again with both nuclei** (`Planner.retake_run`, `runs.retake`) makes a new run with the old run's
    setup, seed and beam time; the old run stays, and an extension of it stays a run of one nucleus.
- The analysis searches for the peaks (`analysis.find_peaks`, `analysis.identify`). Each Doppler-corrected
  spectrum is searched over its whole range, every peak found is fitted, and one on a transition of the nucleus
  that correction is right for is named. It is a step of the chain ("peak search") and `Result.search` holds
  the spectra and the peaks.
- The analysis integrates the peak of each nucleus and shows it (`Result.lines`): the fit, the area with the
  random coincidences subtracted, the efficiency, the angular correlation, the γ rays emitted and the yield
  per particle. For the other nucleus the yield is compared with what its level scheme gives to all
  orders (0.993 measured over expected on a test run of the template).
- B(E2) can be normalised to the other nucleus's measured peak (`Settings(normalisation="target")` without a
  `reference`): the yield over the number of particles that peak stands for, so the beam current, the
  target thickness, the dead time and the particle count cancel. `Settings.reference` keeps the analytic
  reference for a setup without that scheme.
- B(E2) to all orders (`Settings.order`, `Result.b_first_order`, the "all orders" step). When the events follow
  a level scheme to all orders, the analysis solves the scheme's coupled equations for the B(E2) that gives the
  measured γ rays per particle in the peak's window, with the feeding of the state and every line inside the
  window counted (the 4⁺ → 2⁺ line of ⁶²Ni is 3% of it). First order read 4 to 8% high on the template. The
  excitation saturates (the yield grows as B(E2) to the power 0.67 there), so the statistical uncertainty and
  the budget are multiplied by 1.5. Simulated samples of the template give 859 and 910 e²fm⁴ for 884 put in
  (two seeds), a five-minute test run 947 ± 3.5% (statistical).
- The Analysis step shows the analysis in three parts: the peaks found in each corrected spectrum (figure and
  table), each integrated peak with its fit and its numbers, and the chain to B(E2). The analysis kept with a
  run is shown when the step opens. The report has the same figures and tables.
- The correlation table for the other nucleus (`gamma.correlation_table(nucleus="other")`), the list of γ-ray
  lines of both nuclei (`gamma.transition_lines`) and names of levels with their order (`gamma.level_names`:
  2+₁, 2+₂).
- The scene on the whole screen: a ⛶ button on the Setup tab's scene toolbar hides the setup panel, the status
  strip, the tabs and the panels beside and below the scene, and the scene fills the window (Esc or the button
  brings the page back). `SceneView(..., height=)` sets the scene's height.
- Contaminant reactions, beam spot and halo (backlog 74). `[target] contaminants` (`[["C", "10 ug/cm2"], ...]`) adds
  thin layers on the target's downstream face, each a channel of elastic scattering with its own kinematic line
  on the Data tab's energy-against-ring view; the backing's lines are drawn too. `[beam] halo_fraction` and
  `halo_radius` spread that share of the beam evenly over the frame's aperture, in the event generator. The
  generator now reports where each reaction happened (`x`, `y`). The γ rays leave from that point, so the spot and
  the halo blur both the particle and the γ side of the Doppler correction, and the scattering angles sampled reach
  as far as the halo does. Not included: fusion-evaporation and transfer, and particles that hit the frame itself.
- Lifetimes (backlog 73): a level with a half-life in the setup's level scheme decays after a time drawn from its
  lifetime, with the levels above it adding theirs in a cascade. Meanwhile the nucleus slows down along its path
  through the target and backing (the time to stop is ∫ dR / (ρ v) from the stopping tables and the layer
  densities), so its γ ray is emitted at the velocity it has then. A nucleus that has stopped emits at rest. The
  spectra show the stopped and the shifted components and the slowing-down continuum between them, and γ rays
  carry `stopped`. A short lifetime gives item 55's result. Not included: angular straggling, time dilation, and
  fitting the lifetime from the line shape.
- Summing and pile-up (backlog 72). A calibration source whose lines come in cascade (⁶⁰Co, ⁸⁸Y;
  `response.CASCADES`) sums in a crystal: each line loses the share of decays whose partner also interacts there,
  and a sum peak appears with the product of the peak efficiencies. The efficiency points are corrected for
  summing and say so (`efficiency`, `uncorrected`, `summing_correction`). Two γ rays of one cascade (item 71) that
  interact in one crystal give one signal at their sum. `[run] shaping_time` adds pile-up (`response.pile_up`):
  the share 1 − exp(−2τR) of the signals leaves its bins and reappears above, as a shoulder. It applies to source
  runs, singles and coincidence spectra (a coincident γ ray piles with the crystal's singles).
- Cascades (backlog 71): with a level scheme for the excited nucleus whose decays follow more than one
  transition, each excitation starts in a level drawn from the direct excitation probabilities at its scattering
  angle. It then decays through the scheme's branches, with internal conversion taking its share; each γ ray is
  emitted with its transition's orientation and followed to the crystals. True γ–γ coincidences appear
  (`GammaEvents.gamma_gamma`, `coincidences()["gamma_gamma_true"]`), and a gate on one line selects those in
  coincidence with it (`dataviews.gamma_spectrum(..., gamma_gate=)`). The generator picks excited events by the
  total excitation probability of the scheme. γ rays carry `cascade`, `initial` and `final`, and runs store them.
  One excited state keeps item 55's chain unchanged, event for event. Not included: the angular correlation between
  successive γ rays of one cascade, and the kinematics of each level (the reaction's state sets them).
- The Physics tab (backlog 70): the physics behind the numbers as a reference, in the order of the physics register.
  It covers nuclear data and level schemes, two-body kinematics, energy loss, the Rutherford orbit and safe
  distance, the detectors and the γ-ray response, rates and events, and Coulomb excitation with the correlation and
  Doppler shifts. Each section links to its register page and names the tests that validate it. The ? of a
  number ("See the physics") and the consequence line of a setup card open the section that shows its physics
  (`workbench.PHYSICS`, `physics_for`). Nothing on the tab needs a run.
- The Report tab records the experiment (backlog 69, `physim.nuclear.logbook`). It covers the setup and how it
  changed between runs, the Plan, every run with its summary and whether it was taken with the current setup,
  predicted against measured, the gates, the analysis of the current run, and links to the physics register. It
  can be shown in the app, downloaded as a zip (the report, the run record, the setup, the beam-time report of the
  current setup, and every run's folder), or downloaded as the experiment folder. The run record now also explains
  the run's numbers: the measured rate, the counts over the whole run, busy and live fractions, the coincidences,
  each gate's share and the pull of measured against predicted (`logbook.run_explanations`). Each recomputes from
  the numbers it shows.
- The Analysis tab reads the run (backlog 68). It opens with predicted against measured: the Plan's numbers kept in
  the run's summary beside what the run measured, with whether they agree within three standard deviations. That
  covers the rates and full-energy-peak coincidences per detector, the beam time for the counts wanted, and each
  γ-ray detector's efficiency at the transition, measured by the experiment's latest source run
  (`Planner.compare`). The analysis takes the current run and a gate from the Data tab by name, with no
  simulation of its own; on a run it says its statistics come from the run's real events, and what the whole
  scaled run would give (`Result.statistics_from`, `whole_run_statistical`). The alignment check defaults to an
  alignment run's offset and shows the run's γ rays against crystal in the assumed geometry. The multistep fit
  reads the run's yields (`Planner.run_yields`, `fit_matrix_elements()` with no arguments).
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
- The planner in its own window (backlog 61): `physim app --window`, and the installed shortcut (`--desktop`), open
  the same app in a desktop window on the system's web view (pywebview; the `window` extra, 5 MB). `--browser`
  keeps the browser; without pywebview the browser opens and the command says so.
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
