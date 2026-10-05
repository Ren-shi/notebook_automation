# 54 · Orientation of excited states, particle–γ correlation and decay

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

**Status: Done, in Python (see the decision below).**

> **Done** (`python/physim/nuclear/orientation.py`, `python/physim/nuclear/angular.py`, `correlation_table` and
> `excitation_of` in `gamma.py`, `emission` in the `[reaction]` section, per-detector γ-ray efficiency in
> `rates.py`, `Planner.correlation` and `Planner.populations`, two blocks in the app's "Excitation and γ rays"
> tab, `tests/python/test_nuclear_orientation.py`, guide section "Orientation and the particle–γ correlation"
> in `docs/nuclear-setup.md`).
> - **Decision: the first-order amplitudes are in Python, not in the Rust core.** They reuse the orbit
>   integrals that `coulex.py` already computes (now kept, so the two share them), which makes the sum over
>   substates equal the existing probability to rounding. The interface is the one the design note asks for:
>   `Excitation(beam, target, energy, scheme)` gives populations and tensors at a scattering angle. Item 58
>   puts the coupled-channel solution behind it, and that integration is what needs Rust.
> - **Amplitudes** for every substate, any ground-state spin, E1, E2 and E3 together, target or projectile.
> - **Decay:** branching ratios, internal conversion, mixing ratios, and feeding with the parent's density
>   matrix handed down.
> - **Laboratory:** the orbit's frame is placed from the particle's angle, and the γ rays are transformed from
>   the moving nucleus's frame with the change of solid angle.
> - **In the rates:** each particle detector has its own γ-ray efficiency, the sum over γ-ray detectors of
>   efficiency × correlation factor. `emission = "isotropic"` leaves only the forward throw.
> - **Results:**
>   - Summed over substates, the probability equals the first-order one to rounding, 3 × 10⁻¹⁶ for the example (E1, E2 and E3; ground
>     states of spin 0, 1/2, 3/2 and 2; target and projectile).
>   - For 0⁺ → 2⁺ → 0⁺ the correlation equals the quadrupole radiation pattern of the amplitudes,
>     |Σ a_M X_2M|², to 10⁻¹⁰ at four scattering angles. In the sudden limit it equals the closed form built
>     from the analytic orbit integrals, to 5 × 10⁻⁶. At 180° it is (15/8π) sin²α cos²α about the beam, that
>     is 1 + (5/7)P₂ − (12/7)P₄, to 10⁻⁹.
>   - The same check with random complex amplitudes fixes the phase convention of the tensor formula.
>   - Integrated over all directions the correlation returns the γ yield to 10⁻¹⁰, and to 10⁻⁶ in the
>     laboratory with a moving emitter.
>   - With conversion the γ yield falls by exactly 1/(1 + α); feeding matches U_k for every rank and component.
>   - The ⁵⁸Ni example's particle rates are unchanged. Its coincidence efficiency is 1.15 times the isotropic
>     value for the backward ring and 0.83 for the forward strip detectors. Ge90 sees 0.42 of the isotropic
>     number in coincidence with the ring; Ge45 and Ge135 see 1.5.
>   - With the real ENSDF scheme of ⁵⁸Ni at 30 MeV: 5.6 mb for the first 2⁺ state, 1.5 × 10⁻⁴ mb for the second.
> - **Not done here, or to confirm:**
>   - The sign of the interference term of a mixed transition follows Krane and Steffen's formula with ENSDF's
>     δ. I have not checked it against a published mixed-transition distribution.
>   - The correlation factor is computed for the single `[reaction]` state (0⁺ ground state assumed) at the
>     mid-target energy, without shadowing, and is applied per particle detector, not per ring.
>   - The table of all levels is first order: a 4⁺ state appears only through feeding. Item 58 changes that.
>   - Simulated events do not carry γ rays yet; item 55 samples them from these distributions.

## Why
The Coulomb-excitation code returns a probability only (item 43). It does not keep how the excited nucleus is
oriented, so γ rays can only be emitted isotropically. The real emission pattern depends on the scattering angle,
and it changes the yield in each crystal. This item also prepares the amplitudes that multi-step excitation
(item 58) needs.

## Scope
- **Amplitudes:** first-order excitation amplitudes for each magnetic substate of each level in the level scheme of
  item 50, where today one probability is computed for one state.
- **Any ground-state spin,** averaged over the initial substates.
- **Orientation:** the statistical tensors of each excited level, as a function of scattering angle.
- **γ emission:** the particle–γ angular correlation from those tensors, in the frame of the emitting nucleus, then
  transformed to the laboratory (the transformation also changes the solid angle).
  - Isotropic emission remains as a switch.
  - Mixed transitions use the mixing ratio.
- **Decay:** each level decays by its branching ratios. Internal conversion removes its share of the γ rays.
  Feeding from higher levels is followed, with the orientation carried down.
- **Excitation of both nuclei:** the projectile and the target are each excited in the field of the other,
  independently.
- **Left out:** multi-step excitation and reorientation (item 58); deorientation in vacuum; lifetimes (all decays
  happen in flight after the target); γ–γ angular correlations.

## Design notes
- The amplitudes live in the Rust core behind one interface that takes a level scheme and an orbit and returns the
  populations and tensors. Item 58 replaces the first-order solution behind the same interface.

## Depends on
43, 50.

## Done when
- Summed over substates, the probability equals the existing first-order result to 1e-8.
- For a 0⁺ → 2⁺ → 0⁺ sequence the angular correlation matches the published closed form at several scattering
  angles.
- Integrated over all emission directions, the correlation leaves the total γ yield unchanged.
- With conversion included, the γ yield of a transition falls by 1/(1 + α).
- The existing `coulex_ni58` example gives unchanged particle rates.
