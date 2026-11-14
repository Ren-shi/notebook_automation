# Changelog

All notable changes are listed here. The project follows [semantic versioning](https://semver.org/):
before 1.0, minor versions (0.x) may change the API; patch versions only fix bugs.

The version lives in `Cargo.toml` and is shared by the Rust crate, the Python package
(`physim.__version__`) and the PyPI distribution `physim-engine`. To release: bump the version in
`Cargo.toml`, move the "Unreleased" notes under it here, merge, and push a tag `vX.Y.Z`; the release
workflow builds the wheels, tests them and publishes to PyPI.

## Unreleased

### Added
- Deorientation in vacuum and the E1 polarisation correction are both on by default (the two entries below). Together
  they change every Coulomb-excitation plan with a level scheme: on the ⁶²Ni on ¹⁹⁴Pt template the ⁶²Ni
  coincidences fall by 6.5% (CD) and 7.7% (strip detectors) and the beam time for the counts wanted rises by 8.3%
  (9.8 to 10.6 min); with clovers by 10.1% and 11.2%, and 12.6% (3.7 to 4.2 min). `coulex_ni58` has no
  half-lives, so only the E1 polarisation acts there: 0.6 to 0.8% fewer coincidences. `deorientation = "none"`
  and `e1_polarisation = 0` give the numbers of before.
- Deorientation in vacuum, on by default (`[reaction] deorientation`: `"vacuum"`, or `"none"` to switch it
  off). An excited nucleus that decays
  after leaving the target is deoriented by the hyperfine field of its highly charged ion: its statistical tensors
  become G_k t_kq, so its γ rays come out more evenly. The model (`physim.nuclear.deorientation`) is GOSIA's: the
  two-state model of Bosch and Spehl (Z. Phys. A 280 (1977) 329) and Brenn et al. (Z. Phys. A 281 (1977) 219)
  with GOSIA's default parameters (manual of 10 May 2012, section 4.4 and table 4.1). Where the manual leaves
  something open, the arithmetic of GOSIA's subroutines GKK and XSTATIC and its table of atomic spins is followed.
  - The time-integral coefficients are GOSIA's closed form, checked against a line-by-line transcription of GKK
    to 10⁻⁶. The time-differential G_k(t) that the events need is the closed-form convolution of the two-state
    model with GOSIA's interpolating function. It is 1 with no slope at t = 0, tends to the static hard core at
    long times, and its average over the decay is GOSIA's coefficient to 2·10⁻⁵.
  - In the events each level is attenuated through its own time in vacuum: from the later of its formation and
    the nucleus leaving the stack, to its decay, at the exit speed. In a cascade the density matrix is attenuated
    before each γ ray is drawn and before it feeds the next level. A converted transition, which has no time
    drawn, takes the average over its lifetime. Levels that decay in the target, or after the nucleus stopped,
    and levels without a half-life are not deoriented.
  - The plan (`correlation_table`, the rates, the analysis) takes each level's G_k averaged over its decay, times
    the chance that it decays outside the target, for both nuclei of the collision. The feeding from above is
    deoriented as well (`Faded`). On the ⁶²Ni on ¹⁹⁴Pt template the events agree with it crystal by crystal for
    both nuclei: 400 000 reactions, the change of each line in each particle × γ-ray pair against the plan's,
    mean pulls +0.03 and +0.07.
  - The template's ¹⁹⁴Pt 2⁺ level (60 ps; the recoils leave the target at β ≈ 0.04) keeps G₂ = 0.20,
    G₄ = 0.08; the ⁶²Ni 2⁺ level (2 ps) G₂ = 0.92, G₄ = 0.78. The plan's ⁶²Ni coincidences drop by 1.5% (CD) and
    2.9% (strip detectors), the ¹⁹⁴Pt ones by 8% and 10%, and the beam time for the counts wanted rises by 3.0%.
    With clovers the drops are 5–7% and 20–23%, and the beam time rises by 7.2%. `coulex_ni58` has no
    half-lives and does not change. Deoriented ⁶²Ni events analysed without deorientation read B(E2) 2.7% low.
    Five-minute runs of both templates agree with the plan on every rate and coincidence row of
    `Planner.compare()`, with deorientation as without it.
  - A setup switch with its guide entry, and a toggle under the angular-correlation table of the app. γ rays
    with deorientation take 25 to 35% longer to make.
- The E1 polarisation correction to Coulomb excitation, as GOSIA has it, on by default
  (`[reaction] e1_polarisation`). The virtual excitation of the giant dipole resonance weakens every E2 coupling,
  transitional and diagonal, by the factor (1 − z a/r) along the orbit, z = d E_cm A / Z² of the excited nucleus
  (GOSIA manual eq. 3.27–3.28, after Alder and Winther, appendix J; for projectile excitation as GOSIA's code takes
  it). The setting is GOSIA's parameter d. Left out it is 0.0038, what GOSIA uses since April 2023 when its input
  does not set it (its manual has 0.005); `e1_polarisation = 0` switches it off. Every setup's numbers change.
  - It is in the coupled equations (`CoupledChannels`), first order (`Excitation`, `Coulex`) and so in the events,
    the plan's rates, the angular correlation, `Multistep` and both B(E2) readings of the analysis; the cache of
    the level schemes' solutions keys on it, and the GOSIA input file writes it as `DIP,`.
  - At 180° in the sudden limit it lowers the probability by (1 − 0.4 z)², as worked by hand, to 1e-6.
  - With the default against none: the ⁶²Ni 2⁺ state of `coulex_ni62_pt194` is excited 3.8% less (z = 0.053), the
    levels it reaches in two steps 6 to 9% less, and the plan's excitations fall 5.1% in the CD and 4.9% in the
    strip detectors; the ¹⁹⁴Pt 2⁺ state is excited 0.6% less (z = 0.022) and the ⁵⁸Ni 2⁺ state of `coulex_ni58`
    0.6% less (z = 0.0064). A ⁶²Ni run made with the correction and analysed without it reads B(E2) 7.2% low. The
    reorientation of the ⁶²Ni 2⁺ state shrinks by 1.2% of itself, but the correction lowers the ratio of its
    populations at 160° and 60° by 1.4%, as a diagonal matrix element 28% smaller would.
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

### Documentation
- The README shows the application: what the experiment planner does, step by step, with screenshots
  (`docs/images/`), how to start it (Windows installer or `physim app`), the example setups and `physim report`.
  The engine follows in a shorter summary; its full reference moved from the README into `docs/guide.md`, which
  the documentation site shows unchanged.
- The README's screenshots show the app's new look (the step rail, the dark theme), with a light-theme Setup step
  in `docs/images/planner-setup-light.png`, and the README lists the clover template.

### Fixed
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
  the counts wanted 3.7 min against 9.8 min (4.2 and 10.6 min with deorientation and the E1 polarisation, on by
  default since).
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
