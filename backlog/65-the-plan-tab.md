# 65 · The Plan tab: what you check before asking for beam

**Priority:** P1 · **Size:** M · **Area:** Experiment workbench (redesign)

## Why
Rates, safe rings, beam time, γ efficiency and the correlation are spread over four tabs (Rates and beam time,
Energy loss, Kinematics, γ) organised by physics. Before beam, an experimentalist checks one page: will this
measure what I want, how long will it take, what is unsafe.

## Scope
- **One page, answers first:** the counts wanted and the beam time it needs; per detector the rate, dead time,
  the share of excitations, the safe rings (greyed where unsafe); per γ detector the efficiency at the transition
  and the coincidence rate; the particle × γ coincidence rates. Each number carries its ? (item 57).
- **Below the answers, the checks:** the kinematics plot and the energy-loss table, the orbit and the safe
  distance, the particle–γ correlation table, the γ efficiency curves and the sweep — the present tabs' content,
  as collapsible panels on this one page, in the order one checks them.
- **The predictions are kept** with the experiment (`summary.json` of the next run records the Plan's numbers at
  the time), so item 68 can show predicted against measured.
- The beam-time sweep stays here; the GOSIA input file and the multistep solver move to Analysis (item 68).

## Depends on
64.

## Done when
- The Plan page shows beam time, rates, safe rings, efficiency and the coincidence matrix without any button
  pressed, for every example as a template.
- The four physics tabs' panels are reachable from Plan (and from Physics, item 70) and nothing of them is lost
  (the existing figure tests pass).
- The plan's numbers are written into the next run's summary.
