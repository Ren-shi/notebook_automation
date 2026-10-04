# 47 · Backward detectors and the Coulomb-excitation beam time

**Priority:** P1 · **Size:** S · **Area:** Nuclear planner

**Status: Done.**

## Why
From trying the guided mode (2026-10-05): there was no obvious way to put a silicon detector at backward angles.
A detector could be placed there by typing θ, but **Add** always made one at 45°, and nothing explained why one
would go backward.

Behind that sat a real gap. The Coulomb-excitation beam time was worked out from every particle a detector counts,
almost all of them elastic. So `coulex_ni58` asked for 3 seconds of beam, where a real run needs hours, because few
events excite the state and fewer still are seen with their γ ray.

## Scope
- **Counting:** for Coulomb excitation, `counts_wanted`, the counts in the run and the beam time refer to
  particle–γ coincidences: excitation events in a particle detector, times the γ detectors' full-energy-peak
  efficiency. Without γ detectors they refer to excitation events.
- **Efficiency:** γ detectors gain `efficiency`. Without it, the geometric coverage (1 − cos α)/2 is used and
  flagged as an upper limit.
- **Rates table:** for Coulomb excitation it adds excitation and coincidence rates per detector. The report says
  what is counted.
- **Add menu:** ready-made forward strip, side pad, backward pad and backward ring (CD), each with a line on why
  you would put a detector there.
- **Readings:**
  - the geometry reading flags a setup with nothing backward, and names the detector with the largest share of
    excitation events;
  - the rates reading says what is counted and whether the efficiency is only geometric.

## Done when
- A backward detector can be added from the menu.
- The Coulomb-excitation beam time follows from the coincidence rate (tested).

> **Done** (`Rates.measured`, `Rates.gamma_efficiency`, `what=` on the rate methods; `GammaDetector.efficiency`,
> `geometric_efficiency`, `peak_efficiency`; `guide.PLACEMENTS`, `guide.placement`; the Add menu, rates columns and
> readings in the app; the report's counting note).
> - **Tests** (`tests/python/test_nuclear_guide.py`):
>   - `test_placements`: each preset is valid, forward and backward land where they should, and the CD covers
>     about 125–165°;
>   - `test_coulomb_excitation_beam_time_is_set_by_coincidences`: rows equal excitation rate × efficiency, the beam
>     time is `counts_wanted` / coincidence rate and over an hour, a measured efficiency replaces the geometric one,
>     and elastic setups are unchanged;
>   - `test_geometry_reading_points_backward`.
> - **For `coulex_ni58`:**
>   - the four germanium crystals cover 7.3% geometrically;
>   - 2000 coincidences need 4.3 h with the forward strips or 6.0 h with the CD;
>   - with 1.5% measured efficiency per crystal it is 5.2 h;
>   - in the CD, 1 in 518 particles comes from an excitation.
