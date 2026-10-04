# Request: Rutherford cross sections at the detector angles of the example setups (LISE++)

For backlog items 39 and 40 (see `../README.md` for the file format). Count rates are flux × target atoms ×
∫ dσ/dΩ dΩ. The flux and atom counts are arithmetic, and the solid angles are checked against closed forms. The cross
section at each detector is the part to check against a tool.

In LISE++'s kinematics calculator, reaction "Scattering", read the **lab** Rutherford cross section of the scattered
beam at each lab angle below. Make the target very thin (0.001 µm) or put the reaction point at the target entrance,
so LISE++ uses the full beam energy.

| Beam | Target | Beam energy (MeV) | Lab angles (deg) | physim dσ/dΩ lab (mb/sr) |
|---|---|---|---|---|
| ⁴He | ¹⁹⁷Au | 5.5 | 20, 30, 45, 60, 90, 135 | 1.17619e6, 2.38327e5, 49864.7, 17110.3, 4276.9, 1467.0 |
| ¹⁶O | ²⁰⁸Pb | 64 | 35, 60, 105, 125, 148 | 16645.6, 2176.2, 342.0, 218.3, 157.8 |

The values above are at the angles shown, from `Rutherford(beam, target, E).cross_section_lab(θ)`.
`export.lise_settings` prints the same quantities at each detector's mean angle.

Save as `tests/reference/nuclear/rates_lisepp.csv` with columns
`beam,target,beam_energy_mev,theta_lab_deg,cross_section_lab_mb_sr`. The validation check
`physim.nuclear.validation.rates_vs_lise` reads it: 0.5% per angle.
