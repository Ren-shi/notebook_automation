# 53 · γ-detector response and calibration sources

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

**Status: Done, except the values marked typical, which a measured efficiency curve replaces (the user).**

> **Done** (`python/physim/nuclear/response.py`, `data/nist_attenuation.csv`, `data/calibration_sources.csv`,
> `material`, `absorbers` and `efficiency_curve` on a γ-ray detector in `experiment.py`, activity units in
> `quantity.py`, `Planner.efficiency` and `Planner.source_run`, the γ-ray detector's side panel in `app.py`,
> `tests/python/test_nuclear_response.py`, guide section "γ-ray response" in `docs/nuclear-setup.md`).
> - **Efficiency:** solid angle × attenuation × k (1 − exp(−μL)) × peak-to-total, per crystal.
>   - For germanium, k = 0.687 is fixed so that a 50 × 70 mm crystal has the 21.5% relative efficiency of
>     Mirion's clover sheet (against 1.2 × 10⁻³ for the NaI standard at 25 cm).
>   - A measured `efficiency_curve` in the setup replaces the model; a single `efficiency` still works.
> - **Attenuation:** NIST mass attenuation coefficients (Hubbell and Seltzer) for 22 elements, 5 keV to 20 MeV,
>   read from NIST's tables on 2026-10-05. The chamber wall, the setup's `absorbers` and the housing window are
>   applied, at normal incidence.
> - **Response function:** full-energy peak, Klein–Nishina continuum for one scattering, a flat part between the
>   Compton edge and the peak, and single- and double-escape peaks above 1.022 MeV.
> - **Resolution:** FWHM² = noise² + (FWHM(1332)² − noise²) E/1332 keV. This gives the clover's 1.05 keV at
>   122 keV from its 2.1 keV at 1332 keV, and LaBr₃'s 2.9%, 2.1% and 1.6% at 662, 1332 and 2615 keV.
> - **Calibration sources:** ²²Na, ⁶⁰Co, ⁸⁸Y, ¹³³Ba, ¹³⁷Cs and ¹⁵²Eu, with the DDEP evaluated lines read from
>   LNE-LNHB on 2026-10-05. The expected spectrum is computed and the counts drawn bin by bin from a Poisson
>   distribution, so a run of any length is exact and takes about a second.
> - **In the app:** a selected γ-ray detector shows its efficiency against energy and what the γ rays pass
>   through; **Run the source** adds the spectrum and the efficiency from each strong peak. Checked in the
>   running app with ¹⁵²Eu.
> - **Results:**
>   - A ¹⁵²Eu run (37 kBq, 1 h) gives back the efficiency put in at every strong line that stands alone, within
>     4σ; without counting noise the peak areas agree to 0.2%. The 964 keV line reads 1% high, from its real
>     neighbour at 963.4 keV.
>   - The ⁶⁰Co peak-to-total in the clover is 0.1883 against 0.1882 put in.
>   - Doubling a lead absorber squares its transmission, to rounding.
>   - The default clover at 25 cm has 0.099% at 1332 keV (four crystals, no add-back). The range published
>     for it is the data sheet's 21 to 22% per crystal, that is 0.101 to 0.106% for four, before the 2% the end
>     cap absorbs. I did not compare with a measured clover from the literature.
> - **What changed for existing setups:** a γ-ray detector without `efficiency` used its geometric coverage, an
>   upper limit. It now uses the model, so coincidence rates fall: in the ⁵⁸Ni example the γ-ray efficiency
>   goes from 7.3% to 0.69%. The example's `counts_wanted` went from 2000 to 500 so that its 24 h run still
>   reaches it. The report says when the efficiency is the typical one.
> - **Typical values, not from a document** (listed in `response.CRYSTALS[...].typical`):
>   - the peak-to-total ratio at 1332 keV and its slope (0.18 and 0.68 for germanium, 0.20 and 0.75 for LaBr₃);
>   - LaBr₃'s factor k (1.0), so its absolute efficiency is the least certain number here;
>   - the window thicknesses (1.5 mm and 0.5 mm of aluminium);
>   - the shares of the escape peaks and of the flat part of the continuum.
> - **Not done here:**
>   - LaBr₃'s own background from ¹³⁸La, and the room background (item 55).
>   - The simple peak area misjudges a peak on a Compton edge: the ⁶⁰Co peak-to-total in LaBr₃ reads 3% low.
>     Item 56 fits peaks properly.
>   - Efficiency curves and absorbers are edited as text in the app (absorbers) or in the setup file (curves);
>     there is no curve editor.
>   - The simulated beam spectra do not use the response function yet; that is item 55.

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
