# 55 · Particle–γ events, Doppler correction, coincidences and background

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

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
