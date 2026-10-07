# 67 · The Data tab: every spectrum visible, gates and conditions as objects

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench (redesign)

## Why
The γ spectra hide behind a drop-down; a particle spectrum needs a click on the 3D scene and then a button; the
Doppler correction, the gate, add-back and the randoms are fixed by whichever block one is in. Nothing says that
spectra exist. After a run, the data should be in front of you, and the operations should be yours to apply.

## Scope
- **A grid of spectra** for the current run: one small panel per particle detector and one per γ-ray detector
  (per crystal on request), all drawn at once; click a panel to enlarge it and work on it.
- **Operations, applied to a spectrum, not pre-baked:**
  - **Doppler correction:** off / for the projectile / for the recoil, with the kinematic reason beside the
    choice (item 55's `doppler_correct`);
  - **particle gate:** a detector, a ring or strip range, the elastic or excited group (`Settings` of item 56
    becomes a gate the user sees);
  - **coincidence window and randoms:** singles or coincidence, randoms subtracted or shown;
  - **add-back and suppression** on or off for a clover (item 62's `plain`);
  - **binning and range.**
- **Gates and conditions are objects:** named, listed on the tab, editable, saved with the experiment, and reused
  by the Analysis (item 68). A gate drawn on a 2D view (below) becomes one of them.
- **2D views:** energy against ring (or strip) per particle detector with the kinematic lines of each channel
  drawn over it — how one sees the elastic and excited groups; γ energy against crystal, raw and corrected — how
  one sees a misalignment (item 60's diagnostic becomes a view of this).
- **Compare runs:** overlay two runs (with and without add-back; before and after moving a detector), each
  labelled with its setup snapshot.
- **Export** from here: ROOT (`write_root`, with the gates as cuts where ROOT allows), CSV of a spectrum, the
  `.npz` of the run.
- The small status scene in the side: click a detector to select its panels.

## Depends on
63, 66.

## Done when
- After a run every detector's spectrum is on screen without a further click; the γ spectra of the example show
  the Doppler-corrected peak at the transition energy when "for the recoil" is chosen (item 55's test on the
  tab's own spectra).
- A gate made on the energy-against-ring view selects the excited group and the Analysis uses it by name.
- Two runs overlay with their snapshots named.
- ROOT and CSV exports of a run open with the same counts as the panels.
