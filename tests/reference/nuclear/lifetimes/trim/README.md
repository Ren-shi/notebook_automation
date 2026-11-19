# TRIM references for the recoil transport (lifetimes, wave 1C)

Written by `scripts/transport_reference.py "D:/Projects/srim/SRIM2013"` (SRIM 2013, run unattended through
`physim.external.srim`) on 19 November 2026; read by `tests/python/test_nuclear_transport.py`.

Two recoils of the lifetime experiments in `../published/setups.toml`:

- ⁹⁸Ru at β = 4.5% (92.46 MeV; its speed after the 0.7 mg/cm² Mg target, Radeck et al. 2012) into Cu, the stopper;
- ¹⁴³Sm at β = 1.6% (17.04 MeV, Rajbanshi et al. 2018) into Au, the DSAM backing.

Files:

- `TRIM_<case>_<run>.IN`: the TRIM inputs, written by `physim.nuclear.export.trim_in` (5000 ions, quick
  Kinchin–Pease, random seed 1) with the disk-file line changed. `range` runs (a layer thicker than the range)
  save the stopped ions' positions (`RANGE_3D.txt`); `quarter` and `half` runs (layers of about a quarter and half
  the range) save the transmitted ions (`TRANSMIT.txt`).
- `transport_trim.csv`: the summaries: projected range, its straggling and the lateral spread (µm), and for the
  transmitted ions the fraction, mean energy and its spread, rms and median angle to the normal.
- `srim_stopping.csv`: SRIM's electronic and nuclear stopping for both recoils (SR Module, 10 keV to 1.5 × the
  starting energy), the stopping TRIM slowed them with.

The transport test runs physim's Monte Carlo with SRIM's stopping (from `srim_stopping.csv`) in place of physim's, so
that it and TRIM differ only in how they scatter. With physim's own stopping the ranges are shorter (by about 10% for
⁹⁸Ru in Cu and 30% for ¹⁴³Sm in Au): physim's electronic stopping of slow heavy ions (below about 0.1 MeV/u it
follows Lindhard–Scharff) is 1.4–4 times SRIM's there, which `test_physim_stopping_against_srim_for_the_recoils`
pins.
