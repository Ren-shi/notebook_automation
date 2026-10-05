# 55 · Particle–γ events, Doppler correction, coincidences and background

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

**Status: Done, in NumPy on the generator's events (see the decision below).**

> **Done** (`python/physim/nuclear/gamma_events.py`, `[run]` fields `coincidence_window`, `dead_time`,
> `room_background` and `extra_lines`, `threshold` on γ-ray detectors, `Geometry.segment_centre`,
> `Excitation.coefficient_table`, `Planner.gamma_events` and `Planner.gamma_spectra`, the γ block of the app's
> Spectra tab, the `gammas` tree in `rootio.py`, `tests/python/test_nuclear_gamma_events.py`, guide section
> "Particle–γ events" in `docs/nuclear-setup.md`).
> - **Decision: the γ-ray chain runs in NumPy on the particle events the Rust generator makes,** not inside the
>   generator. Each step is a vectorised array operation over the excited events (emission by rejection
>   sampling from the correlation, the boost, the crystal hit, the response, the correction), so there was no
>   per-event loop to move to Rust; the generator stays as it was. The γ rays of 400 000 reactions take about
>   3 s of the 5 s total.
> - **Each excited event emits ten times over,** each γ ray carrying a tenth of the event's rate: coincidences
>   are rare, and this gives the spectra statistics without more reactions. An event whose two particles both
>   reached a detector is a coincidence with each.
> - **Results:**
>   - The recoil-corrected peak of the ⁵⁸Ni example sits at 1454 keV within its uncertainty in every crystal;
>     for a ¹⁶O beam excited at 70 MeV the projectile-corrected peak sits at 6917 keV within its uncertainty.
>     Correcting for the wrong nucleus smears the peak.
>   - With isotropic emission, the spread of the raw γ-ray energy agrees with item 43's table to 25% (and the
>     mean to 1.5 keV) for every pair with enough events. With the correlation the spreads differ, because the
>     table spreads the γ rays evenly over a crystal: the crystal at 90° in coincidence with the backward ring
>     sees a 39% wider peak, as the sin² cos² pattern favours its edges.
>   - Coincidence counts agree with the analytic rates (excitation rate × total efficiency × correlation
>     factor) within 3σ for all twelve particle–γ pairs, with 5% allowed for the factor's coarse quadrature.
>   - Random coincidences equal 2τ × (particle singles) × (γ singles) exactly, and the Poisson draw agrees within
>     statistics.
>   - The same setup and seed give identical γ rays.
> - **Not done here, or to confirm:**
>   - One γ ray per excitation (the single `[reaction]` state): no cascades, so no true γ–γ coincidences.
>   - The room-background line strengths are typical shares, not a measurement of any hall; the user sets the
>     total rate.
>   - Selecting the experiment or a detector in the scene does not yet show the γ matrices and spectra; they
>     are on the Spectra tab.
>   - Summing, pile-up, lifetimes and contaminant reactions, as the scope says.
>   - Events per second on this machine: about 90 000 reactions/s in the generator and 2 000 γ rays/s of
>     the chain (10 emissions each), for item 59.

## Why
The γ results today are analytic: a shifted energy and a width for each pair of particle and γ detector. An
experimentalist works with events: a particle in one ring and sector, a γ ray in one crystal, corrected for the
Doppler shift using what the detectors know. The spectrum must also carry the background a real one has.

## Scope
- **Event chain,** in the Rust generator, for each scattering:
  1. the excitation, from the populations of item 54;
  2. the particle's path to a ring and sector;
  3. the decay, with γ rays emitted according to the angular correlation from the moving nucleus;
  4. each γ ray's path through any material to a crystal, and its recorded energy from the response of item 53.
- **Doppler correction** from the centre of the ring and sector that fired and the centre of the crystal that
  fired, with the velocity of the emitting nucleus reconstructed from two-body kinematics.
  - Corrected for the projectile and, separately, for the target recoil.
  - The remaining peak width follows from the sizes of the segments; it is not put in by hand.
- **Coincidences:** particle–γ and γ–γ, within a time window the user sets. Thresholds apply to every channel.
- **Background:**
  - the Compton continuum of every line, from the response function;
  - random coincidences, from the singles rates and the time window;
  - room background lines (⁴⁰K, the thorium and uranium series) at a rate the user can change;
  - extra lines or a continuum added by hand, for contaminants we cannot predict.
- **Beam current:** entered as electrical current with a charge state, or as particle current.
- **Dead time:** a simple correction from the total rate.
- **In the app:** selecting the experiment shows the particle–γ and γ–γ counts and matrices; selecting a detector
  shows its raw and Doppler-corrected spectra.
- **Export:** the ROOT file of item 42 gains the γ branches and the coincidence events.
- **Left out:** reactions on target contaminants computed from first principles; pile-up; summing; lifetimes.

## Design notes
- Coincidences are rare. Generate only scatterings that reach a particle detector, weight the events, and report
  how many events of a real run the sample stands for.
- Record the events per second the generator reaches on the reference machine; it sets how large a run item 59 can
  animate.

## Depends on
39, 53, 54.

## Done when
- The Doppler-corrected peak sits at the transition energy within its fitted uncertainty, for the projectile and
  for the target.
- The corrected peak's width agrees with the analytic broadening of item 43 for the same geometry.
- The random-coincidence rate equals the product of the singles rates and the window, within statistics.
- Event counts agree with the analytic rates within 3σ for every channel.
- The same setup and seed give identical events.
