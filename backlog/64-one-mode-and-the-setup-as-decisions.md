# 64 · One mode: the setup as a list of decisions, the status strip, new experiments

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench (redesign)

## Why
The left panel is the TOML file: sections in file order, every field an input. Guided mode is the same form
again with Next buttons; expert mode is the form without help. Neither tells an experimentalist what to do next or
what a choice costs. The example drop-down at the top ("coulex_ni58", "alpha_on_gold") gets in the way of setting
up one's own experiment.

## Scope
- **One mode.** The stepper and the guided/expert switch go. Guidance comes from the structure: the panel is a
  list of decisions in the order they are made — Beam → Target → Particle detectors → γ-ray detectors → Run
  conditions — each a collapsible card showing its current state in one line and its consequence in a second
  ("Target: ⁵⁸Ni 0.5 mg/cm² · beam loses 1.2 MeV, recoils leave at up to 9°"; "CD at 30 mm behind · 125–165°,
  2.3 sr, 4 800/s"). Opening a card shows its fields; rarely touched fields (dead layer, window gap, crystal
  pitch, housing, absorbers) sit behind a "more" fold. The ? help icons and the "how to read this" notes stay.
- **The checks are always visible** at the top of the panel (warnings, overlaps, safe-angle notes, "setup changed
  since run 3"), not on a tab.
- **The status strip** across the top of the main area, from every tab: experiment name · beam and energy ·
  target · detectors (3 particle + 4 γ) · last run (kind, duration, when) · warning count. Clicking a part opens
  the matching card.
- **Experiments, not examples.** Opening the app shows the experiment list (item 63) and **New experiment**, which
  asks beam, target and what to measure (elastic / Coulomb excitation of the beam or the target) and makes the
  name from them; detectors come after, from the Add menus or a **template** (the present examples, renamed:
  "Coulomb excitation of a target nucleus with a CD and clovers", "Elastic scattering with a DSSD", "Blank").
  The example drop-down disappears; the example files remain as templates and for the tests.
- **Load/Save** become "Import a setup file" and "Export the setup" on the experiment, since the experiment folder
  is the saved state.
- **Tabs** left to right: Setup · Plan · Run · Data · Analysis · Report · Physics (items 65–70). The scene is
  central on Setup; elsewhere it is a small status scene (which detector is which, click to pick) in the side.
- `app.py` is split by tab (one module per tab plus the panel and the strip) so each stays readable.

## Depends on
63 (experiments and runs), 52 (the scene).

## Done when
- A new user can make an experiment from "New experiment" with beam, target and detectors from the Add menus
  without seeing a form of every field, and the status strip names it correctly.
- Every decision card shows its consequence line, and the line changes when the field does (tested on the
  target's thickness and a detector's distance).
- The guided and expert modes, the stepper and the example drop-down are gone; the guide's help entries still
  cover every field (the existing test).
- The app's tabs are the seven stages in that order.
