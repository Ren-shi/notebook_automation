# 59 · Particle and γ-ray tracks in the scene

**Priority:** P2 · **Size:** M · **Area:** Experiment workbench

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
