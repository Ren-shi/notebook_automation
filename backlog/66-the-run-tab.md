# 66 · The Run tab: one button, a duration, live counters

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench (redesign)

## Why
There is no simulation tab. The simulation is nine buttons on other tabs, each with an event count that means
nothing to an experimentalist. The lab has one action: start the run for a duration, watch the counters, stop or
let it finish.

## Scope
- **Run for ⟨duration⟩**, the beam time of the setup by default, editable ("10 min", "2 h"). Beside it, the
  kind of run: beam; source (which source, activity, position); alignment check (the assumed offset). The
  budget of real events ("simulate the first 10 min for real", item 63) is shown, with the CPU time it implies.
- **Progress in experiment time:** a clock and a bar in beam time, mapped on the events generated; a note when
  the real part ends and the scaled part begins.
- **Live counters per detector**, ticking as the run goes: counts, rate, dead time and live fraction; per γ
  detector and crystal: singles, coincidences; the particle × γ matrix filling in; the counts-wanted bar. These
  are the Plan's numbers now measured.
- **Stop** keeps what is accumulated; **Extend** adds beam time to the same run. The setup panel is read-only
  while a run is in progress (item 63).
- **Watch events:** the tracks of item 59 as a view of the run (a sample of the events as they come, in the small
  scene), with the event panel; clicking a detector in the scene shows what is happening to it (its rate, its
  last events, its spectrum so far).
- **Not simulated** list on the tab, pointing at the backlog items (71–74): cascades, summing and pile-up,
  lifetimes, contaminant reactions, beam spot and halo. The realism is claimed only where it holds.
- The run list (item 63) sits on this tab: every run with its kind, duration, date and a one-line summary;
  select one to make it current for Data and Analysis.

## Depends on
63, 64.

## Done when
- One click runs the setup for the duration and the counters end at the summary's numbers; a 10-minute run on
  the ⁵⁸Ni example finishes in under a minute with progress shown throughout.
- Stop and Extend work and the run list shows the result; a source run and an alignment run appear in the list.
- The setup cannot be edited during a run (tested through the Planner's lock).
- The "not simulated" list names every item in 71–74.
