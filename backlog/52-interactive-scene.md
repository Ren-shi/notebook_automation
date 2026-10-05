# 52 · The interactive scene

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

**Status: Done.**

> **Done** (`python/physim/nuclear/scene.py` for the geometry and the rules, `python/physim/nuclear/scene_view.py`
> for the drawing, the side panel in `app.py`, `Planner.move`, `place`, `live`, `selection` and `safety`, the quick
> and reusing paths of `Rates`, `tests/python/test_nuclear_scene.py`, guide section "The scene" in
> `docs/planner-app.md`).
> - **Decision: NiceGUI's built-in three.js scene, not a custom component.** It draws solids, reports clicks with
>   the object under the mouse, and drags objects with constraints evaluated in the browser. Two things it lacks
>   were added around it rather than by replacing it:
>   - it reports only the start and end of a drag to Python, so the view also listens to its "drag" events
>     (throttled) for the live numbers;
>   - it reports a click after the view has been turned, so a few lines of JavaScript drop clicks whose press and
>     release are more than 4 pixels apart.
> - **Select, then drag.** Only the selected detector can be dragged. The scene has one set of drag constraints,
>   and this way they are always those of the object being moved.
> - **Limits are enforced twice.** The browser holds the drag with the constraints (the beam, the target, the
>   chamber wall, and in distance the other solids). Python applies the same rule again on release and checks for
>   overlaps, which in angle depend on both angles: a detector dropped on another one returns to the last allowed
>   place, with the reason. A test evaluates the JavaScript constraints and compares them with the Python rule.
> - **Overlaps** are found from the edges of each solid against the volume of the other, so two thin detectors
>   that cross are caught whatever their thickness.
> - **Geometry-only updates:** a move keeps the rate rows of every detector that did not change and that nothing
>   shadows. During a drag, only the moved detector is computed, on a coarse grid.
> - **Results:**
>   - The S3 dragged from 20 to 40 mm gives the same setup, ring angles and rates as typing "40 mm".
>   - The live values take 0.10 to 0.14 s per move on the development machine for the examples (a CD of 384
>     segments, a strip detector of 256). The quick rate agrees with the full one to better than 0.1%.
>   - 120 random drags per mode on the Coulomb-excitation example never end in the beam or in an overlap.
>   - After a move in the Coulomb-excitation example the rates take 0.5 s, against 1.0 s for all detectors. Where
>     the moved detector becomes shadowed it is integrated on the fine grid, which takes longer (2.1 s in the lead
>     example).
>   - Checked in the running app, in both themes: selecting a detector and a ring, sliding the CD along the
>     beam, moving a germanium detector over its sphere, and a drop on a neighbour being refused.
> - **Typical sizes used for drawing only:** the target foil's diameter (10 mm), the board's thickness (1.6 mm),
>   and the length of a γ-ray detector without crystal dimensions.
> - **Not done here:**
>   - After a drop, the other tabs still take a few seconds (the Monte Carlo and every figure are redone, as
>     after any edit). The scene and the side panel update first.
>   - Typed values are not held back by the limits; the warnings report a detector in the beam and solids that
>     overlap.
>   - The side panel's spectrum is per detector, not per ring.
>   - The Plotly outline (`figure_geometry`) is kept for scripts; the paper figure of the layout is unchanged.

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
