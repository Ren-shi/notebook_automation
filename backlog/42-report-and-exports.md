# 42 · Beam-time report and data exports

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner · **Status: Done, except opening the ROOT file in ROOT
(the user)**

> **Done** (`python/physim/nuclear/report.py`, `python/physim/nuclear/rootio.py`, `tests/python/test_nuclear_report.py`,
> guide section "The beam-time report and data files" in `docs/nuclear-setup.md`).
> - **Report:** `build(experiment, seed)` collects every number into a `Report`; `write(folder)` saves it.
>   - `report.html` is self-contained, with print styles so the browser's "Save as PDF" gives the PDF from the same
>     source. It covers the setup, all warnings, the detector table (θ, φ, Ω, dσ/dΩ, rate, Monte Carlo rate,
>     counts in the run, beam time), kinematics, expected peaks, energy loss, spectra, each model's register status,
>     data sources, references, and the version and commit.
>   - Five CSV tables with units in the column names.
>   - Figures as PNG (200 dpi) and PDF, and `setup.toml`.
>   - A command line: `python -m physim.nuclear.report setup.toml -o folder --seed N`.
> - **ROOT** (`rootio.write_root` and `read_root`, written with uproot 5.7, approved by the user 2026-10-04; written
>   by `write` when uproot is installed). `events.root` holds:
>   - a TTree `events` with every column;
>   - a TH1D `spectrum_<detector>` per detector, in counts for the run, with Monte Carlo errors √Σw² per bin;
>   - the setup file and run info as TObjStrings.
> - **Results:**
>   - The same setup file and seed give identical numbers (tested).
>   - The ROOT file reads back in uproot unchanged: every branch bit for bit, the spectra's totals and errors, and
>     the setup. A Gaussian on the DSSD3 spectrum recovers the analytic peak.
>   - Both example reports are checked against hand calculations. For α + Au at 45°: rate within 1% of
>     I n (dσ/dΩ) Ω; dσ/dΩ within 0.2%; foil energy loss within 1% of S × t; peak within 2 keV. For ¹⁶O + Pb: DSSD3
>     rate within 3% of I n σ Ω; target loss within 1%; kinematic table equal to `TwoBody`.
>   - Building a report takes about 8 s.
> - **Waiting on the user:** opening the file in ROOT itself, with `TTree::Draw` and `Fit("gaus")`. ROOT is not
>   installed here; the steps are in `tests/reference/nuclear/pending/root_check.md`.
> - **Not done:** a PDF written directly; the browser's print to PDF covers it without another dependency.

## Why
The report is the main deliverable of the planner: what a student attaches to a beam-time proposal or brings to a
meeting with their supervisor. The data exports let them carry the results into their own analysis, and ROOT is the
standard format for that in nuclear physics.

## Scope
- **Report** (HTML, and PDF from the same source):
  - setup summary: beam, target, detectors, run conditions;
  - 3D geometry figure and kinematics plot;
  - per-detector table: θ and φ coverage, solid angle, dσ/dΩ, count rate, counts in the planned run, beam time for
    N counts;
  - energy loss summary, simulated spectra;
  - every validity and rate warning;
  - assumptions, and the physics register status of each model used;
  - Physim version and commit, data versions (item 34), and references to cite.
- **CSV:** per-detector results and kinematic tables, with units in the column headers, documented in the guide.
- **ROOT** (written with `uproot`, no ROOT installation needed): a TTree of Monte Carlo events (item 39) and
  histograms of the spectra. The setup file is stored with it, in the ROOT file if `uproot` can write it, otherwise
  as a sidecar file.
- **Setup file** (item 33) and figures (PNG and PDF, publication quality).
- Left out: custom report templates, LaTeX output.

## Depends on
33, 39 (tables and figures from 35–38).

## Done when
- Regenerating a report from its setup file and seed gives identical numbers.
- The ROOT file opens in ROOT (TTree::Draw and a histogram fit work) and reads back in `uproot` unchanged.
- The report for each example setup is reviewed against a hand calculation of its key numbers.
