# 54 · Orientation of excited states, particle–γ correlation and decay

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

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
