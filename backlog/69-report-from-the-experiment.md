# 69 · The Report tab: the record of the experiment and its runs

**Priority:** P2 · **Size:** M · **Area:** Experiment workbench (redesign)

## Why
The beam-time report (41) and the run record (57) describe one setup and one simulated sample. With experiments
and runs they should describe the experiment: the setup, each run with its snapshot, the analysis of the run the
user chose.

## Scope
- The report takes the experiment: the setup (with changes between runs noted), the Plan, the run list with each
  run's summary, the Data views and gates used, the Analysis result, and the physics register links.
- The run record (57) covers every number on every tab, including the counters and the gates.
- Exports: the zip as today (report, record, setup, data files), and the experiment folder itself.

## Depends on
63, 68.

## Done when
- The report of an experiment with a source run and a beam run lists both with their snapshots and the analysis
  of the beam run; the record explains the counters and the gate (the reader's check of item 57 applies).
