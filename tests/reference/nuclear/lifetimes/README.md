# Lifetime fixtures (lifetime measurements, wave 0)

References for the plunger (RDDS) and Doppler-shift attenuation (DSAM) checks, made on 19 November 2026 by
`scripts/lifetime_reference.py` (and, for `published/`, by hand from the papers). Kept byte for byte
(`.gitattributes`: `-text`): AhKin's outputs on Windows end lines with CRLF, its inputs with LF.

## `analytic/` — closed forms, no outside program

Written by `scripts/lifetime_reference.py analytic`, which imports nothing from physim.

- `rdds_one_level.csv`: one level, no feeding, recoil at one speed along the beam. Unshifted share after a flight
  d: exp(−d/(γβcτ)); also without time dilation. τ = 1, 10, 100 ps; β = 2, 5, 10%; d = 1 µm to 5 mm.
- `rdds_two_levels.csv`: a level fed from one level above (Bateman's equations), the unshifted shares of both
  transitions; three lifetime pairs (one with τ₁ ≈ τ₂) and two population splits.
- `dsam_exponential.csv`: stopping proportional to the speed, v(t) = v₀ e^(−t/α): attenuation F = α/(α + τ) and
  the distribution of v/v₀ at emission in 20 bins.

## `ahkin/y83_in_ta/` — AhKin's own DSAM example

AhKin (A. Garzón, W. Rodríguez, F. Cristancho, M. Tao; Comput. Phys. Commun. 2020; Mendeley Data `559kb329p2`,
`AhKin.tar.xz`, MIT licence) built with MinGW-w64 gcc/g++ 16.2 (`stopping` needs `-std=c++17`: it reads into a
`char*` with `>>`, which C++20 removed). Its example: ⁸³Y residues from ³²S + ⁵⁸Ni at 135 MeV slowing in Ta, the
1193 keV line at 52.8°, a ten-level population pattern. `prob_cos` (10⁶ rays), `stopping`, then `AhKin_B -i`
between 0.1 and 0.2 ps give **τ = 0.131308 +0.0339269 −0.0326754 ps**, exactly the value in AhKin's readme.
The folder holds the inputs (`in_*.dat`, the stopping table `Y83_in_Ta.txt`, the measured line shape
`EL1192GFAA53.txt`), the kernel and velocity files, `chisqB.dat`, the best-fitting line shape and `result.json`.
Added for the validation (lifetime wave 2D, 19 November 2026), from the same `AhKin_B -i` run, byte for byte: the
line shapes at the two ends of AhKin's 1σ interval, `ech0.0986328.dat` (τ − 0.0327 ps) and `ech0.165235.dat`
(τ + 0.0339 ps), which `chisqB.dat` lists at χ²_min + 1.

The input formats, from AhKin's `readme.pdf`: `in_casc.dat` gives the number of levels, their mean lives (ps; the
last is the level being fitted, its value unused), the branching matrix B with B[i][j] the share of level j's decays
that go to level i (−1 on the diagonal), and the unnormalised initial populations. `in_AhKin.dat` gives E_γ, the
calibration (keV = cal0 + cal1 × channel), the resolution FWHM = r0 + r1 E + r2 E² and the flat background;
`in_stp.dat` the reaction, the 10% Gaussian spread of the initial energies and the stopping table (keV, keV/µm);
`pcos.dat` the distribution of the cosine between the beam and the γ ray over the detector's face.

## `apcad/` — APCAD's example spectra

APCAD (C. Stahl, J. Leske, M. Lettmann, N. Pietralla; Mendeley Data `4hzzcvx2p7`, `APCADStopSim.tar.gz`,
CC BY 4.0) needs CERN ROOT and Geant4 and was not built. Its `Example/` folder holds two simulated continuous-angle
DSAM spectra (29 angles, 20°–160°, 1 keV bins, 500–1500 keV), converted to CSV (angle, energy, counts):
- `OneShape_example.csv`: a 1000 keV line, **τ = 1 ps**, 7400 events, A₂₀ = 0.4, A₄₀ = −0.8 (the README's second
  "A20" is A40), FWHM 3 keV, background 400 − 0.252x + 5.6·10⁻⁵x².
- `ThreeShapes_example.csv`: the same with contaminants at 1025 keV (τ = 0.5 ps) and 985 keV (τ = 1.6 ps).
The setup that made them is in `Macro_StopSim_example.mac` (StopSim, Geant4): a 400 MeV ¹²⁰Sn beam on 6.64 µm C on
11.21 µm Cu, 10⁶ events, with the cross section against angle in `XSec_example.dat`. Read from StopSim's source
(lifetime wave 2D): `TrackBeamLike true` means the nucleus followed and excited to 1000 keV is the **beam-like
¹²⁰Sn**, not the ¹²C (the transfer reaction gives the excitation to the ejectile); `KillInLayer2` keeps only ions
excited in the first layer (the carbon); a nucleus leaving the target keeps its exit velocity. The thicknesses are
1.5 mg/cm² C at StopSim's 2.26 g/cm³ and 10.0 mg/cm² Cu at 8.92 g/cm³. Each entry of `XSec_example.dat` (beam
energy in MeV, the projectile's CM angle in degrees, cross section) is dσ integrated over its 1° bin, so ∝ dσ/dθ
(physim's own first-order Coulex gives the same shape within 3–6%); StopSim interpolates it linearly in energy.
The stopping is Geant4's `G4IonParametrisedLossModel` (ICRU 73 for ions up to Z = 26, scaled above).
`Detector_example.det` is the detector file, `README.txt` APCAD's own description.

## `napytau/` — NaPyTau's DDCM on the analytic decay curves

NaPyTau (github.com/h-mayr/NaPyTau_Public, commit in the file header, EUPL-1.2; a Python reimplementation of
Cologne's NAPATAU) run from its core functions on `analytic/rdds_one_level.csv` (flight times 0.05–3τ, 1%
errors), polynomials of degree 2, 3, 4 over all distances. Degree 4 gives 10.03 ps (τ = 10 ps, β = 5%),
100.02 ps (τ = 100 ps) and 10.04 ps (τ = 10 ps, β = 10%). Its uncertainty propagation gives meaningless numbers
(10⁴⁰ and more) and is not used. NaPyTau needs numpy and scipy (run with an existing Python 3.12 that had them).

`napytau/ddcm_simulated.csv` (lifetime wave 3C, `scripts/lifetime_reference.py napytau-simulated`): NaPyTau on
physim's own simulated plunger spectra. The ⁹⁸Ru example `coulex_ru98_plunger` without its deorientation section
(one exponential, τ = 8.36 ps) simulated at its 11 distances, 4·10⁶ reactions, γ rays in coincidence with the ²⁴Mg
recoils, the 41.5° ring (`physim.nuclear.lifetime_reproduction.napytau_input`); the unshifted and shifted
intensities are the counts in windows (the stopped peak E₀ ± 2.5 FWHM, the rest of the line out to its kinematic
end, less a flat background), with their errors; β = 4.76% from the shifted line's centroids at 41.5° and 138.5°.
NaPyTau got the intensities over I_u + I_s and the distances plus the 8 µm offset, and gives τ(x) per distance for
one polynomial of degree 2, 3 and 4 over all distances (its uncertainties again not usable). The same areas serve
the one-level closure (`rdds_closure`). This file imports physim (it is physim's simulation that both analyse).

## `published/` — the experiments to reproduce

`lifetimes.csv` (the measured lifetimes, with the other methods' values the papers quote for the same states) and
`setups.toml` (beam, target, foils, distances, velocities, detector angles, analysis), read from the papers:
- ⁹⁸Ru, RDDS after Coulex below the barrier: Radeck et al., PRC 85, 014301 (2012).
- ⁵⁰Cr, RDDS after fusion-evaporation with a three-foil plunger (CoCoDiff): Beckers et al., NIM A 1042, 167418
  (2022), with Coulex, RDDS and DSAM values for the same 2⁺ and 4⁺.
- ¹⁴³Sm, DSAM after fusion-evaporation on a Au backing: Rajbanshi et al., PLB 782, 143 (2018).

## Which validation check reads what

The checks of `physim.nuclear.validation` for the capabilities "Lifetimes: recoil distance (RDDS)" and
"Lifetimes: Doppler-shift attenuation (DSAM)" (register page docs/physics-register/lifetimes.md):

| File | Check |
|---|---|
| `analytic/rdds_one_level.csv` | `rdds_one_level_analytic` |
| `analytic/rdds_two_levels.csv` | `rdds_two_levels_bateman` |
| `analytic/dsam_exponential.csv` | `dsam_exponential_analytic` |
| `published/setups.toml`, `published/lifetimes.csv` (⁹⁸Ru) | `deorientation_ru98` |
| `trim/transport_trim.csv`, `trim/srim_stopping.csv` | `rdds_transport_vs_trim` (⁹⁸Ru in Cu), `dsam_transport_vs_trim` (¹⁴³Sm in Au) |
| `apcad/OneShape_example.csv`, `apcad/XSec_example.dat` | `dsam_lineshape_vs_apcad` |
| `ahkin/y83_in_ta/` | `dsam_lineshape_vs_ahkin` |
| `napytau/ddcm_simulated.csv` | `rdds_ddcm_vs_napytau`, `rdds_closure` (one level) |
| `published/` (⁹⁸Ru, ⁵⁰Cr) | `rdds_ru98_reproduction`, `rdds_cr50_reproduction`, `rdds_closure` (fed from above) |
| `published/` (¹⁴³Sm), `trim/srim_stopping.csv` | `dsam_sm143_reproduction`, `dsam_closure` |

The reproductions (wave 3C, `physim.nuclear.lifetime_reproduction`, rerun with `scripts/reproduce_lifetimes.py`)
simulate the three examples made after the papers and analyse the spectra as the papers did; the side-feeding times
of ¹⁴³Sm come from `published/setups.toml` (a setup cannot hold them).

## Licences

The AhKin files are MIT (Garzón et al.), the APCAD files CC BY 4.0 (Stahl et al.); both are kept here as test data
with that attribution. NaPyTau's output is numbers only. The published values are facts quoted with their sources.
