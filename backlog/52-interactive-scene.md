# 52 · The interactive scene

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

## Why
The 3D view is a Plotly figure of outlines (item 41): the target is a square of fixed size, detectors are flat
outlines, and nothing can be selected or moved. The user should see the experiment as it would stand in the chamber
and rearrange it by hand.

## Scope
- **A real scene** replaces the Plotly diagram: the solids of item 51 drawn to scale, with rings, sectors and
  crystals visible, the target at its real size, the beam and the chamber.
- **Selection:** clicking an object selects it, shows its dimensions in the scene, and fills the side panel.
  - A detector shows its angular coverage, solid angle, rates and spectrum. A ring or a crystal can be selected on
    its own.
  - A γ detector shows its efficiency and count rates.
  - Clicking empty space selects the experiment and shows the totals.
- **Constrained dragging:**
  - an S3 slides along the beam axis;
  - a γ detector moves over a sphere around the target and in distance, and keeps facing the target;
  - other detectors move in angle and distance;
  - exact values can always be typed.
- **What a move may not do:** enter the beam, overlap another solid, or cross the chamber wall. The drag stops at
  the limit and says why.
- **Kinematic advice:**
  - backward angles are recommended for a light beam on a heavy target, forward angles for the reverse;
  - a detector placed where the kinematics send nothing is flagged;
  - rings where the scattering is not safe (Cline's criterion) are marked.
- **Live feedback:** angular coverage, solid angle, shadowing and analytic rates update during a drag. The Monte
  Carlo runs again on release or on request.
- **Themes:** light and dark, as item 49.
- **Left out:** particle tracks (item 59); moving the target (item 60); the educational scenes.

## Design notes
- Plotly cannot select or drag objects. NiceGUI's built-in three.js scene is the first candidate; decide between it
  and a custom component, and record the decision here.
- Today every edit clears the planner's whole cache. Dragging needs geometry-only updates that keep what did not
  change.

## Depends on
41, 51.

## Done when
- The example setups open as to-scale scenes; a measured distance on screen matches the setup file.
- Dragging the S3 from 20 mm to 40 mm from the target updates its ring angles and rates to the same values as
  typing the distance.
- No drag can produce an overlap or put a detector in the beam.
- The cheap quantities update within 0.2 s of a move on the reference machine.
- Every panel that existed before is still reachable.
