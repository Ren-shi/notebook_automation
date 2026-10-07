# 74 · Contaminant reactions, beam spot and halo

**Priority:** P3 · **Size:** M · **Area:** Experiment workbench (physics)

## Why
Real runs have more than the reaction planned: scattering on the backing and on contaminants (carbon, oxygen on
the target), beam halo on the frame or the chamber, a beam spot of finite size. They put lines and continuum in
the spectra and particles in detectors that should be clean, and an experimentalist needs to see them before beam.

## Scope
- Elastic scattering on the backing and on named contaminant layers (`[target] contaminants`), each a channel
  with its own kinematic line in the energy-against-ring view (67).
- The beam spot (`[beam] spot_size`, already a field) and a halo fraction on a frame radius: particles from off
  the axis, which blur the kinematic lines and the Doppler correction.
- Left out: fusion-evaporation and transfer channels (a cross-section table could be added later).

## Depends on
63, 67.

## Done when
- A carbon contaminant on the target shows its elastic line at the right energy in every ring; a 3 mm spot
  broadens the corrected γ peak by the amount the geometry gives, within statistics.
