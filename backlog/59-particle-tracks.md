# 59 · Particle and γ-ray tracks in the scene

**Priority:** P2 · **Size:** M · **Area:** Experiment workbench

**Status: Done, except the smoothness check on the reference machine (the user).**

> **Done** (`python/physim/nuclear/tracks.py`, `show_tracks`, `control_tracks` and `select_track` in
> `scene_view.py` with the animation in the browser, `Planner.tracks`, the tracks toolbar and the event panel
> in `app.py`, `tests/python/test_nuclear_tracks.py`, "Tracks" in `docs/planner-app.md`).
> - **The sample:** events picked by the rate they stand for (a run's look) or as generated (rare large-angle
>   scatterings show); all events, particle–γ coincidences, or one channel; the scene states the choice and
>   the range of rates the tracks stand for.
> - **Paths:** the beam to the target; each detected particle to where its track crosses the detector face,
>   which lies in the segment the event names; the undetected partner to the edge; γ rays to the crystal hit.
> - **Animation in the browser:** the balls move along their paths (the beam first, then particles and γ rays
>   at a common speed), what is hit turns red, and the cycle repeats; play, pause and speed are fields the
>   page sets. No per-frame traffic to Python.
> - **Results:**
>   - Every path ends in the segment or crystal the event record names (46 of 46 hits in a sample of 30;
>     the test checks a sample of 40 with the face's own hit test).
>   - With the coincidence filter every track has a particle hit and a γ hit.
>   - Seen in the running app: 30 tracks, hit strips lit, the animation clock running.
> - **Not done here, or to confirm:**
>   - Smoothness at the default (30 tracks) could not be measured here: the app's browser pane throttles the
>     frame rate when it is not the active window. The default is set below what a laptop draws smoothly in
>     the author's judgement (about 60 tracks), to be confirmed on the reference machine.
>   - Paths are straight (no bending in fields) and stop at the face; nothing inside a detector.
>   - A new sample is needed after any edit (the events are made again).

## Why
Seeing particles leave the target, bend and strike the detectors makes the scene read as an experiment, and it
shows at a glance where the events go.

## Scope
- **A sample of the Monte Carlo events** of item 55 is animated: the beam particle to the target, the scattered
  projectile and the recoil on their paths, and γ rays as straight lines from the moving nucleus.
- **Hits:** the ring, sector or crystal that records a hit lights up. Selecting a track shows that event's
  numbers.
- **Filters:** all events, coincidences only, or one excitation channel.
- **Controls:** play, pause, speed, and the number of tracks shown.
- **Honest sampling:** the scene states how the sample was chosen. Large-angle scatterings are rare, so a sample
  that shows them is weighted, and the scene says so.
- **Left out:** tracks inside a detector; every event of a run.

## Design notes
- Find how many tracks the scene draws smoothly on the reference machine, and set the default below that.

## Depends on
52, 55.

## Done when
- Each track ends in the segment that the event record names.
- The animation stays smooth at the default number of tracks on the reference machine.
- With the coincidence filter on, every track shown has a particle hit and a γ hit.
