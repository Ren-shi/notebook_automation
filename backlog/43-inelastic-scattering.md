# 43 · Inelastic scattering and Coulomb excitation (second slice)

**Priority:** P2 · **Size:** L · **Area:** Nuclear planner · **Status: in progress**

> **Progress (2026-10-04).** The physics is done; integration into rates, events, app and report is next.
> - **Done:**
>   - `physim.nuclear.coulex` (Coulomb excitation for E1, E2 and E3, safe distance, warnings);
>   - `physim.nuclear.gamma` (Doppler table);
>   - setup file `[reaction] type = "coulex"` and `[[gamma_detectors]]`, with the example `coulex_ni58`;
>   - `docs/theory/coulex.md`;
>   - `tests/python/test_nuclear_coulex.py` (18 tests).
>
>   Results so far: the closed forms at 180°, ξ = 0 hold to 4e-10; engine-orbit integration agrees to 5e-8; the
>   Doppler formula equals the four-vector boost to 1e-12. ¹⁶O on ⁵⁸Ni at 30 MeV: ξ = 0.84, P(170°) = 0.28%,
>   σ = 7.7 mb, Doppler shift +22 keV at 45° for CD events.
> - **Next:** inelastic channels in `Rates`, peaks and the Rust event generator (Q = −E*, P(θ) as the event weight);
>   the inelastic lines and Doppler table in the app and report; validation checks in `physim.nuclear.validation`.
> - **Waiting on the user:** GOSIA and a published measurement (`tests/reference/nuclear/pending/coulex_gosia.md`).

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
