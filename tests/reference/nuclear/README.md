# Reference data for the nuclear planner

Results from established codes (LISE++, SRIM/TRIM, ...) that physim is checked against. These programs are closed
and run by hand, so each file records exactly how it was produced, so that anyone can reproduce it.

- `pending/` holds **requests**: the exact cases still to be run, with the settings to use and physim's own values
  for comparison. When a request has been run, put the results in this folder in the format below and delete the
  request.
- Files here are read by the tests (`tests/python/test_nuclear_*.py`), which assert the tolerances stated in the
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
