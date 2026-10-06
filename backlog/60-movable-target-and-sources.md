# 60 · A movable target and sources placed anywhere

**Priority:** P2 · **Size:** M · **Area:** Experiment workbench

**Status: Done.**

> **Done** (`position`, `ladder` and `selected` on the target and target-relative positions in `experiment.py`,
> `python/physim/nuclear/alignment.py`, `recorrect` in `gamma_events.py`, `Response(..., source=)` in
> `response.py`, `target_offset_mm` and the position entry of the budget in `analysis.py`, `show_ghost` in
> `scene_view.py`, `Planner.alignment`, the "Check the alignment" block in the app,
> `tests/python/test_nuclear_alignment.py`, guide section "A misplaced target" in `docs/nuclear-setup.md`).
> - **Target position:** detectors stay placed from the chamber's centre, as the file writes them; the setup
>   gives each an origin, and every position, direction and distance the physics uses is from the target. A
>   setup with the target at the origin is unchanged (the rates are identical).
> - **Misalignment:** the Doppler correction can be made again with any geometry; the overlay gives the shift
>   and the broadening; the diagnostic plot gives the centroid per crystal and ring with the run's statistics;
>   the fit scans the offset and takes the minimum of χ² with Δχ² = 1.
> - **The diagnostic is judged against a simulation with the assumed geometry,** not against flat lines: with
>   the right geometry the corrected centroids are still ring-dependent by about a keV (the 16° crystals and
>   the energy loss), which against flat lines read as a 2.4 mm offset. The reference simulation costs one
>   more simulation per fit (14 s in all for 600 000 reactions).
> - **Results:**
>   - Moving the target by 5 mm changes each CD ring's angles to the hand values, 180° − atan(r/35 mm), to 0.01°.
>   - With no offset the diagnostic is flat within statistics against the reference (χ² below 1.5 per point)
>     and the fit gives −1.0 ± 1.1 mm; with the target assumed 2 mm off, the fit gives it back: −2.9 ± 1.3 mm
>     (600 000 reactions; the ⁵⁸Ni example's slow recoil makes 2 mm worth about a keV).
>   - A source moved 60 mm towards a crystal changes its rate by exactly the solid-angle ratio, for every
>     crystal.
>   - The ladder picks the target in the beam and round-trips through the setup file.
> - **Not done here, or to confirm:**
>   - The sensitivity of the offset fit is that of the simulated sample, not of the run (a run of 24 h has
>     a hundred times the counts): the app says so through the uncertainty; more events sharpen it.
>   - A detector's own misplacement (rather than the target's) is not fitted; the analysis takes the target's
>     place as the uncertain one.
>   - The source off centre keeps the crystals facing the target (the solid angle is exact, the attenuation at
>     normal incidence).
>   - The beam stays on the axis (as the scope says).

## Why
The target is fixed at the origin and the beam runs along +z. Real setups have a target ladder, a target that sits
off centre, and sources placed where the target would be or elsewhere. The user should be able to test those
without rebuilding the setup.

## Scope
- **Target position:** the target can move along the beam axis and be tilted. Angles, solid angles and Doppler
  corrections follow its position.
- **Misalignment:** the analysis may assume a position for the target (or for a detector) that differs from the
  true one. The Doppler correction then uses wrong angles, so each corrected peak shifts and broadens, by an amount
  that differs from ring to ring and from crystal to crystal. Shown in three ways:
  - **In the scene:** the assumed position as a faint outline next to the true one, with a control for the
    offset.
  - **Spectrum overlay:** the corrected peak with the true geometry and with the assumed geometry, with the shift
    of the centroid and the change in width as numbers.
  - **Diagnostic plot:** the corrected peak's centroid against ring number, one line per crystal. The lines are
    flat when the geometry is right and sloped when it is not. This is the plot used on real data to find a
    misplaced target.
- **Effect on the result:** the change in B(E2) that the offset causes, which enters the uncertainty budget of
  item 56 as the position contribution.
- **Finding the offset:** the app fits the offset that flattens the diagnostic plot, as one would with real data.
- **Target ladder:** several targets in one setup, one of them in the beam.
- **Sources:** a calibration source (item 53) can be placed at any point in the chamber.
- **Left out:** α sources for the silicon detectors; a beam that is off axis or at an angle.

## Depends on
52, 53.

## Done when
- Moving the target by 5 mm changes each ring's angle to the value a hand calculation gives.
- With no offset the diagnostic plot is flat within statistics; with a 2 mm offset along the beam the fit recovers
  2 mm within its uncertainty.
- A source moved off centre changes each crystal's count rate as its solid angle changes.
- A setup with the target at the origin gives the results it gave before this item.
