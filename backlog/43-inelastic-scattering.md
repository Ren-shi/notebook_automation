# 43 · Inelastic scattering and Coulomb excitation (second slice)

**Priority:** P2 · **Size:** L · **Area:** Nuclear planner · **Status: Done, except the GOSIA and published-data
comparisons (the user)**

> **Done** (`python/physim/nuclear/coulex.py`, `python/physim/nuclear/gamma.py`, the excitation channels in
> `rates.py`, `events.py` and `src/nuclear/events.rs`, the γ tab in `planner.py` and `app.py`, the report section,
> `examples/coulex_ni58.toml`, `docs/theory/coulex.md`, `tests/python/test_nuclear_coulex.py` and additions to the
> events, planner, app, report and validation tests, `tests/nuclear_events.rs`).
> - **Physics:** first-order semiclassical Coulomb excitation for E1, E2 and E3 from a 0⁺ ground state, along the
>   symmetrised orbit, with Cline's safe distance and warnings for η and large P. Doppler-shifted γ energies and
>   broadening per pair of particle and γ detector.
> - **Setup file:** `[reaction] type = "coulex"` and `[[gamma_detectors]]`.
> - **Integration:**
>   - excitation channels in `Rates` (Q = −E*, Rutherford × P(θ)), with Coulomb-excitation peaks always listed;
>   - the Rust generator (`Channel.excitation`, `excite_recoil`, `p_table` weighting events by P(θ*), sampled as
>     often as elastic);
>   - inelastic lines in the kinematics;
>   - the app's "Excitation and γ rays" tab;
>   - the report's Coulomb-excitation section and `gamma.csv`;
>   - two register capabilities with validation checks and a notebook.
> - **Results:**
>   - The closed forms at θ = 180°, ξ = 0 hold to 4e-10.
>   - Integration along engine orbits agrees to 5e-8.
>   - The adiabatic cutoff falls monotonically with ξ, below 1% by ξ ≈ 6.
>   - Scaling with B and with the exciting charge is exact.
>   - The Doppler formula equals the four-vector boost to 1e-12.
>   - Monte Carlo agrees with the analytic rates within 1.5σ for every channel of `coulex_ni58`, at 0.7% statistical
>     error.
>   - The safe-distance warning appears exactly at Cline's angle.
>   - ¹⁶O on ⁵⁸Ni at 30 MeV: ξ = 0.88, total 5.9 mb, P(170°) = 0.28%, and the CD records 1.3 excitations per second
>     at 0.1 pnA. Doppler shifts are ±22.5 keV at 45° and 135° for CD events, with 13–18 keV FWHM.
> - **Waiting on the user:** GOSIA and a published measurement (`tests/reference/nuclear/pending/coulex_gosia.md`),
>   which `validation.coulex_vs_gosia` reads as soon as the file exists.
> - **Left out:**
>   - multi-step excitation and reorientation;
>   - nuclear-Coulomb interference;
>   - γ angular distributions, lifetimes and feeding;
>   - an ENSDF lookup of level data, which needs a download approval.

> Refined 2026-10-04, after the elastic slice (33–42). What changed from the outline:
> - Level data come from the setup file (the user types the state's energy and B(E2)); an ENSDF download is left
>   for later because it needs approval.
> - GOSIA and published Coulomb-excitation data are requested from the user, like the LISE++ and SRIM runs.
> - Everything else reuses the elastic slice: kinematics with Q = −E*, rates, the event generator, the app and the
>   report.

## Why
The natural next step after elastic scattering: the beam particle (or the target nucleus) is left in an excited
state, so part of the kinetic energy goes into excitation. Below the Coulomb barrier this is Coulomb excitation, a
standard way of measuring collective properties of nuclei (B(E2) values). It reuses almost all of the elastic slice.

## Scope
- **Setup file:**
  - `[reaction] type = "coulex"`, with one excited state: `excite = "target"` or `"projectile"`, `energy`,
    multipolarity `E1`, `E2` or `E3` (from a 0⁺ ground state, so even-even nuclei), and `b_up`, the reduced
    transition probability B(Eλ↑) in `e2b2` or `e2fm4` (and `e2b`/`e2fm2` for E1, `e2b3`/`e2fm6` for E3).
  - Optional `[[gamma_detectors]]` (θ, φ, distance, radius, resolution).
- **Coulomb excitation** (`physim.nuclear.coulex`): first-order semiclassical theory (Alder and Winther).
  - The excitation probability P(θ) is integrated along the symmetrised hyperbolic orbit for any λ and ξ.
  - dσ/dΩ = P(θ) dσ_R/dΩ, in the CM and the lab; total and per-detector cross sections; the adiabaticity ξ.
- **Safe distance:** Cline's criterion, closest approach ≥ 1.25 (A₁^⅓ + A₂^⅓) + 5 fm, gives a validity warning per
  detector, using the CM angles it covers.
- **Rates, peaks, events, app, report:**
  - inelastic channels next to the elastic ones;
  - kinematics with Q = −E*;
  - excitation probability as an event weight in the Rust generator;
  - inelastic lines in E vs θ and in the spectra.
- **γ rays** (`physim.nuclear.gamma`) from the excited nucleus: Doppler-shifted energy and Doppler broadening for
  each pair of particle and γ detector. The broadening covers the particle detector's angular spread, the γ
  detector's opening angle and the velocity spread through the target.
- **Validation:**
  - closed form at θ = 180°, ξ = 0;
  - independent integration along orbits from physim's engine;
  - adiabatic suppression with ξ;
  - GOSIA (input requested from the user, `pending/coulex_gosia.md`);
  - a published Coulomb-excitation measurement (requested from the user).
- **Left out:**
  - multi-step excitation and reorientation (second order: GOSIA's domain);
  - nuclear-field interference above the barrier;
  - γ angular distributions and internal conversion;
  - feeding from higher states.

## Depends on
33–42.

## Done when
- **Accuracy:**
  - The excitation probability matches the closed form at θ = 180°, ξ = 0 to 1e-6.
  - It matches independent integration along engine orbits at other angles to 1e-4.
  - It falls with ξ as the adiabatic limit requires.
- **Scaling:** cross sections scale exactly with B(E2) and as Z₁² (target excitation).
- **Example:** an example Coulomb-excitation setup gives inelastic lines, rates, spectra and γ-ray energies
  (Doppler shift matching the relativistic formula to 1e-9) in Python, the app and the report.
- **Warnings:** the safe-distance warning appears exactly when Cline's criterion fails.
- **Register:** GOSIA and literature comparisons are 🟡 until the user supplies them.

## After this
Transfer reactions, then fusion-evaporation, each as its own slice.
