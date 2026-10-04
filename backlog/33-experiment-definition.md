# 33 · Experiment definition: the setup file

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner · **Status: Done**

> **Done** (`python/physim/nuclear/` — `experiment.py`, `quantity.py`, `names.py`, `_toml.py`;
> `tests/python/test_nuclear_setup.py`, `examples/nuclear/`, `docs/nuclear-setup.md`).
> - `physim.nuclear.Experiment` with `Beam`, `Target` (+ backing `Layer`), `Detector`, `Run`: `load`, `from_toml`,
>   `from_dict`, `to_toml`, `save` (validates first), `validate` (re-checks edits made in Python).
> - Every physical value is text with a unit (`Quantity`), kept exactly as written; bare numbers are rejected.
>   Particle (`pnA`) and electrical (`enA`) currents are separate units, and plain `nA` is refused as ambiguous.
> - Detectors: rectangle (double-sided strips), annular (rings and sectors), circle; placed by θ/φ/distance or by
>   position, facing the target unless `facing` is given, with optional rotation about the normal.
> - All problems are collected and reported together, each naming the field (`detector 3 (D3): distance must be
>   positive, got '-40 mm'`), with spelling suggestions for unknown fields.
> - TOML is read with `tomllib` (Python 3.11+) or `tomli` (added as a dependency for Python < 3.11) and written by a
>   small built-in writer, so no new package is needed on 3.11+.
>
> Results: 43 tests — both example setups round-trip byte for byte; 30 invalid-input cases each produce a message
> naming the field; derived values (total energy from MeV/u, particle current from electrical current and charge
> state, detector positions) checked by hand.
>
> Left out: beam emittance (only spot size and energy spread); material names are checked for syntax only, with
> densities and compositions coming in item 34. `reaction.type` accepts only `"elastic"` until item 43.

## Why
The nuclear experiment planner (items 33–43) is built around one description of an experiment: the beam, the target,
the detectors and the run conditions. The web app (item 41), the report (item 42), the notebooks and scripts all read
and write the same description, so a setup built in the app can be saved, shared with a supervisor, re-run from
Python and attached to a beam-time proposal. Everything else in the slice depends on this contract, so it comes first.

## Scope
- A plain-text setup file (TOML) with a schema version, and a Python object `physim.nuclear.Experiment` that loads,
  validates, edits and saves it.
- **Beam:** nuclide (Z, A, or a name such as `"4He"`), kinetic energy (MeV total or MeV/u), charge state, current
  (particle-nA or electrical nA together with the charge state), optional energy spread and emittance/spot size.
- **Target:** material (element, compound or isotopically enriched), thickness (mg/cm² or µm), tilt angle, optional
  backing layer.
- **Detectors** (any number): type (rectangular pad, annular strip detector, circular), position (θ, φ, distance,
  or a full position + orientation), active area, number of strips, thickness, dead layer, energy resolution and
  threshold.
- **Run conditions:** beam time (hours), target of N counts per detector (for the "beam time needed" output).
- Units always explicit in the file. Bad input gives a readable message naming the field and the problem
  ("detector 3: distance must be positive, got -40 mm"), never a stack trace.
- Example setups in `examples/nuclear/` (a Geiger–Marsden-style α + Au setup, and a modern four-detector array).
- Left out: reaction choice beyond elastic scattering (a `reaction` field is reserved for item 43 onwards), beam-line
  optics, multi-target setups.

## Depends on
—

## Done when
- Load → save → load round-trips every example setup exactly.
- Each invalid-input case in a test list is rejected with a message that names the field.
- The schema is documented field by field, with units, in the user guide.
