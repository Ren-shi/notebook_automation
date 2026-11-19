# Reference data for the nuclear planner

Results from established codes (LISE++, SRIM/TRIM, ...) that physim is checked against. These programs are closed
and run by hand, so each file records exactly how it was produced, so that anyone can reproduce it.

- `pending/` holds **requests**: the exact cases still to be run, with the settings to use and physim's own values
  for comparison. When a request has been run, put the results in this folder in the format below and delete the
  request.
- `trim/` holds the TRIM.IN files (written by `physim.nuclear.export`) that produced `straggling_trim.csv`; a test
  keeps them identical to what the exporter writes today.
- `stopping_measured_heavy.csv` is measured data too: electronic stopping of 21 ions from Ar to U in 12 elemental
  solids (C to Pb) at 5 keV/u–10 MeV/u, 3620 points from the IAEA Stopping Power Database, with SRIM 2013's value
  per point for context (`scripts/iaea_heavy_ion_reference.py`). Read by the check `stopping_heavy_vs_data` and
  `test_slow_heavy_ions_against_measured_stopping`. Cite the original papers by `ref_id` and `doi`.
- `stopping_measured_hydrogen.csv` is measured data, not a code's output: heavy-ion stopping in hydrogen-rich
  targets from the IAEA Stopping Power Database (`scripts/iaea_stopping_reference.py`), with SRIM's value per point
  for context. Its `tool:` line names the database and version; cite the original papers by `ref_id` and `doi`.
- `pace4_*.csv` are PACE4 4.36.0 runs (in LISE++), all written by `scripts/pace4_reference.py`, which reruns
  PACE4 in silent mode one deck at a time in a temporary folder. For each of the five systems of
  `docs/theory/evaporation.md`: `pace4_<beam>_<target>_<E>MeV.csv` is the evaporation at one beam energy (20 000
  cascades; E*, σ_fus, σ_J, mean CM energies and multiplicities in the header, the residues as rows), and
  `pace4_fusion_<beam>_<target>.csv` the excitation function (σ_fus, ℓ_grazing and σ_J at 8–10 beam energies from
  the barrier to about 2 V_B, σ_J semicolon-separated in one column). `pace4_excitation_4He_165Ho.csv` (`--only
  excitation`) holds the residues of ⁴He + ¹⁶⁵Ho at each Ho foil energy of the reproduced experiment and at
  51–59 MeV (20 000 cascades each, one row per energy and residue), drawn beside physim in its figure. PACE4 is a
  Monte Carlo code: a rerun changes the residue counts within statistics, not σ_fus or σ_J.
- `fusion_*_literature.csv` and `evaporation_*_literature.csv` are measured fusion and residue cross sections
  transcribed from the IAEA EXFOR database (tabulated values only, never digitized plots). Their `tool:` line names
  the database version and entries, `reference:` the papers and tables, and each row its EXFOR subentry. Read by
  the checks `fusion_vs_data` (every `fusion_*_literature.csv`; give each row `e_cm_mev` or `e_lab_mev`) and
  `residues_vs_data` and `reproduction_usman2020` (`evaporation_4He_165Ho_literature.csv`).
- `shellmodel/` holds KSHELL inputs and outputs (²⁸Si/USDB, ⁴⁸Cr/KB3), partition files with their dimensions,
  and the interaction files they use, for testing the shell-model parser and size estimates without running
  KSHELL; its own README records how they were made, the benchmarks they pass and their licence (GPLv3).
  `shellmodel/benchmarks.csv` lists published calculations with the same Hamiltonians (ANTOINE for ²⁸Si/USDB,
  the exact ⁴⁸Cr/KB3 energies, M-scheme dimensions), each with its paper and table, for the check
  `shellmodel_vs_published`, which reads the fixture logs with `physim.nuclear.shellmodel.read_folder`.
  `shellmodel/si28_ensdf.csv` is ENSDF's ²⁸Si as physim reads it (`LevelScheme.from_ensdf`), kept so that
  `shellmodel_dimensions` can report the shell model against experiment (rms of the excitation energies, B(E2))
  as context without the local ENSDF copy.
- `lifetimes/` holds the references for the lifetime checks (plunger RDDS and DSAM): closed-form decay curves and
  attenuation, AhKin's and APCAD's DSAM examples, NaPyTau's DDCM answers, and the published lifetimes of the
  experiments to reproduce, all written by `scripts/lifetime_reference.py`; its own README records where each
  outside program came from, its licence and how it was run.
- `reproduction_usman2020.csv` is physim's own result, not a reference from outside: the reproduction of Usman
  et al., NIM B 469, 42 (2020) (⁴He + ¹⁶⁵Ho stacked-foil activation) through the planner, one row per Ho foil and
  channel (channel table filled by the app, residues made in the irradiation, σ inferred from them, the measured
  value beside it), written by `scripts/reproduce_usman2020.py` (`physim.nuclear.reproduction`). The check
  `reproduction_usman2020` reruns the reproduction rather than reading this file; a test keeps the file equal to
  a rerun (same seed), and the validation notebook *Reproduced experiment* draws it.
- `neutron_carbon_endf.csv` is an evaluated library, not a code run: neutron cross sections of natural carbon (total, elastic, non-elastic), ¹²C(n,α) and H-1 from ENDF/B-VIII.0 at 0.5–20 MeV, pointwise and averaged over ±2.5 % and ±10 % of the energy, written by `scripts/endf_neutron_reference.py <folder with the ENDF-6 files>`. The same script writes physim's carbon table (`--toml`). Read by the check `carbon_vs_endf` and by `test_nuclear_neutrons.py`.
- `gemini_<beam>_<target>_<E>MeV.csv` are GEMINI++ runs (shipped with LISE++) for the same fusion-evaporation
  systems as the PACE4 files: residue cross sections as rows, E*, σ_fus, multiplicities and mean
  CM energies of n, p and α in the header. GEMINI++ has no silent mode; `scripts/gemini_reference.py` drives its
  window through Windows UI Automation (`physim/external/gemini_drive.ps1`) on a hidden desktop. GEMINI++ files exist for
  ⁴He + ¹⁶⁵Ho, ¹²C + ¹⁵⁰Nd and ¹⁶O + ¹⁴⁴Sm (1500 events, see its `note:` line); the LISE++ build hangs on
  ³²S + ⁹²Mo and ⁴⁸Ca + ¹²⁴Sn (`docs/theory/evaporation.md`). Read by the checks `evaporation_channels_vs_gemini`,
  `evaporation_spectra_vs_gemini` (physim from the file's E* and its Fermi spin distribution, `spin:` line) and
  `evaporation_neutron_energy_band` (with the PACE4 file of the same system).
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
