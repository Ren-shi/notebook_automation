# Reference data for the nuclear planner

Results from established codes (LISE++, SRIM/TRIM, ...) that physim is checked against. These programs are closed
and run by hand, so each file records exactly how it was produced, so that anyone can reproduce it.

- `pending/` holds **requests**: the exact cases still to be run, with the settings to use and physim's own values
  for comparison. When a request has been run, put the results in this folder in the format below and delete the
  request.
- `trim/` holds the TRIM.IN files (written by `physim.nuclear.export`) that produced `straggling_trim.csv`; a test
  keeps them identical to what the exporter writes today.
- `stopping_measured_hydrogen.csv` is measured data, not a code's output: heavy-ion stopping in hydrogen-rich
  targets from the IAEA Stopping Power Database (`scripts/iaea_stopping_reference.py`), with SRIM's value per point
  for context. Its `tool:` line names the database and version; cite the original papers by `ref_id` and `doi`.
- `pace4_*.csv` are PACE4 4.36.0 runs (in LISE++), all written by `scripts/pace4_reference.py`, which reruns
  PACE4 in silent mode one deck at a time in a temporary folder. For each of the five systems of
  `docs/theory/evaporation.md`: `pace4_<beam>_<target>_<E>MeV.csv` is the evaporation at one beam energy (20 000
  cascades; E*, σ_fus, σ_J, mean CM energies and multiplicities in the header, the residues as rows), and
  `pace4_fusion_<beam>_<target>.csv` the excitation function (σ_fus, ℓ_grazing and σ_J at 8–10 beam energies from
  the barrier to about 2 V_B, σ_J semicolon-separated in one column). PACE4 is a Monte Carlo code: a rerun changes
  the residue counts within statistics, not σ_fus or σ_J.
- `fusion_*_literature.csv` and `evaporation_*_literature.csv` are measured fusion and residue cross sections
  transcribed from the IAEA EXFOR database (tabulated values only, never digitized plots). Their `tool:` line names
  the database version and entries, `reference:` the papers and tables, and each row its EXFOR subentry.
- Files here are read by the checks in `physim.nuclear.validation` and the tests (`tests/python/test_nuclear_*.py`), which assert the tolerances stated in the
  physics register. Until a file exists, the matching test is skipped and the register shows 🟡.

## File format

CSV with a header block of `#` lines first:

```text
# tool: LISE++ <version>, <module used, e.g. "Kinematics calculator">
# produced by: <name>, <date>
# settings: <anything not in the columns: mass table, relativistic on/off, ...>
beam,target,ejectile,excitation_mev,beam_energy_mev,particle,theta_lab_deg,branch,energy_mev,theta_cm_deg
4He,197Au,,0,5.5,ejectile,30,1,5.4701,30.583
```

Leave `ejectile` empty for elastic or inelastic scattering. `branch` is 1 for the higher-energy solution and 2 for
the second solution where the kinematics are double-valued. Energies are kinetic, MeV; angles in degrees.
