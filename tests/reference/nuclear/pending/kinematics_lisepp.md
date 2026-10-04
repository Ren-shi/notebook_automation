# Request: two-body kinematics from LISE++

For backlog item 35 (see `../README.md` for the file format). Please run these in the LISE++ kinematics calculator
with **relativistic** kinematics and the **AME2020** mass table if your version offers it (note which table was
used). Save the results as `tests/reference/nuclear/kinematics_lisepp.csv` with columns
`beam,target,ejectile,excitation_mev,beam_energy_mev,particle,theta_lab_deg,branch,energy_mev,theta_cm_deg`.

The test compares energies to **1 keV** and CM angles to **0.01°**. Physim's values are listed so a mismatch is easy
to spot. If LISE++ disagrees by more than that, keep its numbers anyway; the difference is what we need to see.

| # | Reaction | Beam energy | Particle | Lab angle | physim: energy (MeV) | physim: θ_cm |
|---|---|---|---|---|---|---|
| 1 | ¹⁹⁷Au(α, α) elastic | 5.5 MeV | ejectile | 30°, 90°, 150° | 5.4701, 5.2808, 5.0980 | 30.583°, 91.166°, 150.583° |
| 2 | ²⁰⁸Pb(¹⁶O, ¹⁶O) elastic | 64 MeV | ejectile | 35°, 125° | 62.2398, 50.1922 | 37.538°, 128.624° |
| 3 | ²⁰⁸Pb(¹⁶O, ¹⁶O) elastic | 64 MeV | recoil (²⁰⁸Pb) | 30° | 12.7517 | 60.001° |
| 4 | ²⁰⁸Pb(¹⁶O, ¹⁶O)²⁰⁸Pb*, E* = 2.614 MeV | 64 MeV | ejectile | 90° | 52.4160 | 94.527° |
| 5 | ¹H(²⁰⁸Pb, ²⁰⁸Pb) inverse kinematics | 1000 MeV | ejectile (²⁰⁸Pb), **both solutions** | 0.1°, 0.2° | 999.3479 / 981.4008, 997.0260 / 983.6861 | 21.218° / 158.983°, 46.303° / 134.099° |
| 6 | ¹H(²⁰⁸Pb, ²⁰⁸Pb) inverse kinematics | 1000 MeV | recoil (proton) | 30° | 14.3925 | 60.253° |
| 7 | ¹²C(d, p)¹³C | 10 MeV | ejectile (proton) | 20°, 90° | 12.5065, 10.3617 | 21.967°, 95.753° |
| 8 | ³H(p, n)³He, just above threshold | 1.05 MeV | ejectile (neutron), **both solutions** | 10° | 0.1421 / 0.0166 | 29.766° / 170.236° |

Also useful (not tested yet): the maximum ²⁰⁸Pb angle in case 5 (physim: 0.2776°) and the ³H(p, n) threshold
(physim: 1.0190 MeV).
