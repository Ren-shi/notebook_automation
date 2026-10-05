# 60 · A movable target and sources placed anywhere

**Priority:** P2 · **Size:** M · **Area:** Experiment workbench

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
