# 53 · γ-detector response and calibration sources

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

## Why
A γ detector today detects with one fixed efficiency number and records only the full energy. Peak areas, and the
B(E2) values taken from them, need a realistic efficiency and a realistic spectrum. A transport code is out of
scope, so the response is parametrised.

## Scope
- **Efficiency:** a full-energy-peak efficiency curve against γ energy for each crystal type, scaled by the solid
  angle the crystal covers from the emission point.
  - Typical published curves for a clover crystal and for LaBr₃ are the defaults.
  - A measured curve supplied by the user replaces the default.
- **Attenuation:** material between the target and the crystal (chamber wall, absorbers, housing) reduces the
  efficiency according to its thickness and attenuation coefficient.
- **Response function:** each detected γ ray gives a full-energy peak, a Compton continuum with its edge, and
  escape peaks above 1.022 MeV, in proportions set by a peak-to-total ratio that depends on energy.
- **Resolution:** depends on energy, with separate defaults for germanium and LaBr₃.
- **Calibration sources:** ¹⁵²Eu and ⁶⁰Co (and others from the decay data) can be placed at the target position in
  place of the beam, with an activity and a counting time. The simulated spectrum lets the user measure the
  efficiency of their own arrangement, as in a real experiment.
- **In the app:** selecting a γ detector shows its efficiency curve, with the points a simulated source run gives.
- **Left out:** transport of photons inside the crystal; add-back and Compton suppression (item 62); summing of
  coincident γ rays in one crystal; timing response.

## Design notes
- The defaults are typical values, not a calibration of any real detector. The report must say so wherever a peak
  area depends on them.
- Cite the source of each default curve and of the attenuation coefficients.

## Depends on
51.

## Done when
- A simulated ¹⁵²Eu run gives back the efficiency curve that was put in, within statistics, at every strong line.
- The peak-to-total ratio of the simulated ⁶⁰Co spectrum matches the value put in.
- Doubling an absorber's thickness changes the efficiency as the attenuation law requires.
- The efficiency at 1.33 MeV of the default clover, at 25 cm, lies within the range published for such detectors.
