# 62 · Add-back and Compton suppression for clovers

**Priority:** P3 · **Size:** M · **Area:** Experiment workbench

## Why
Clovers are normally used with add-back, which sums the energies of neighbouring crystals and raises the
full-energy efficiency at high energy. Many also sit inside BGO shields that reject Compton-scattered events. Both
change the peak areas and the background. The first version works without them.

## Scope
- **Add-back:** a parametrised probability, depending on energy, that a γ ray scatters from one crystal into a
  neighbour. With add-back on, those events return to the full-energy peak.
  - The Doppler correction then uses the crystal with the larger energy deposit.
- **Compton suppression:** a BGO shield as a blocking volume around the clover (item 51), with a suppression
  factor that depends on energy and lowers the continuum.
- **In the app:** both are switches for each clover, with the spectrum shown with and without.
- **Left out:** transport of the scattered photon; the shield's own spectrum.

## Depends on
53, 55.

## Done when
- The add-back factor (the ratio of peak areas with and without) at 1.33 MeV equals the value put in, and lies
  within the range published for clovers.
- Suppression lowers the continuum by the factor put in and leaves the peak area unchanged.
- With both switched off, the results equal those of item 55.
