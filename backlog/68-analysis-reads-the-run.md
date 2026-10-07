# 68 · The Analysis tab reads the run and the gates; predicted against measured

**Priority:** P1 · **Size:** M · **Area:** Experiment workbench (redesign)

## Why
The analysis (item 56), the alignment check (60) and the multistep solver (58) each simulate their own sample
and keep their own settings. They should read the run the user took and the gates the user drew.

## Scope
- `Analysis` takes a run and a named gate; its chain (gate, correction, fit, randoms, yield, normalisation, B(E2),
  budget, shape) and its history stay, with the ? explanations (57).
- **Predicted against measured:** the Plan's numbers from the run's summary beside the run's: rates, beam time
  for the counts wanted, efficiency, coincidences. The first thing one looks at after a run.
- **Alignment check** as an analysis of an alignment run (66), or of a beam run against an assumed offset, with
  the diagnostic as the 2D view of item 67.
- **Multistep:** the coupled-channel yields and the matrix-element fit read the run's yields; the GOSIA input file
  download moves here from Plan.
- The Monte-Carlo uncertainty of item 56 draws from the run's real events (unit weights) and from the scaled
  part's Poisson statistics, and says which.

## Depends on
63, 67.

## Done when
- The analysis of a 10-minute run of the example gives B(E2) back within its uncertainty, using a gate made on
  the Data tab, with no simulation of its own.
- Predicted and measured rates agree within statistics for every example template, shown side by side.
- The alignment and multistep blocks read the run (their tests pass on run data).
