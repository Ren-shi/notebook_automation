# Request: reproduce a published fusion-evaporation experiment end to end

This is step D10 of the fusion-evaporation validation, the equivalent of `coulex_literature.md` for Coulomb excitation.
When the reproduction has been run and its result saved here in the usual format, delete this request.

## The choice

**A. R. Usman, M. U. Khandaker, H. Haba, N. Otuka and M. Murakami, "Production cross sections of thulium
radioisotopes for alpha-particle induced reactions on holmium", Nucl. Instrum. Methods B 469, 42 (2020), Table 5;
EXFOR entry E2650.** The measured values are already in `../evaporation_4He_165Ho_literature.csv` (`ref = Usman2020`).

Why this one:

- **It is one of the five PACE4 systems** in `docs/theory/evaporation.md` (⁴He + ¹⁶⁵Ho, 50 MeV), so measurement,
  PACE4 and physim can be shown together, and it is a light, low-spin system, where physim's evaporation is
  meant to work (no fission, no quasi-fission, no incomplete fusion of a composite beam).
- **Everything the app needs is published.** Beam, intensity, irradiation time, foil thicknesses, the energy at each
  foil with its spread, the γ line and branching used for each residue, and the error budget are all in the paper
  and the EXFOR entry.
- **Four exit channels across a whole excitation function** (1n–4n, 24–26 energies from 13 to 50 MeV) from one
  experiment with one normalization. That tests where each channel peaks, not just one number.
- **Small errors (8–10% for 3n and 4n), thin foils, and two independent measurements to check it against**
  (Kormazeva 2022, Gadkari 1997, in the same file). Their spread (up to 1.5× at the 4n peak) is the honest size of
  "agreement".
- It is a recent entry proofread by the authors, with an EXFOR compiler (N. Otuka) among them.

Runners-up, for the record: Visakh et al., Phys. Rev. C 104, 054602 (2021) (¹⁶O + ¹⁴²,¹⁵⁰Nd, total residue cross
section with the HIRA separator, normalized to Rutherford monitors at 15.5°: it would exercise physim's monitor
normalization, but gives no channels and needs HIRA's transmission, which physim does not model); Bhattacharjee
et al., Phys. Rev. C 104, 054607 (2021) (¹⁶O + ¹⁶⁵Ho, 3n/4n/5n by activation, mostly below the barrier, where physim
leaves fusion out). No in-beam γ-ray measurement with tabulated channel yields and a fully stated detector setup was
found in this pass; one would test the γ-ray and neutron-detector stages too and is worth a later search.

## The setup, as the app needs it

| Item | Value | Source |
|---|---|---|
| Beam | ⁴He (α), RIKEN AVF cyclotron | paper; EXFOR FACILITY |
| Energy | 50 MeV onto the stack; energy at each Ho foil from 49.6 down to 13.2 MeV (column `e_lab_mev`), spread ±0.7–1.8 MeV half-width (`e_lab_unc_mev`) | Table 5; EXFOR EN, EN-RSL-HW |
| Intensity | 5.45 × 10¹¹ α/s (87 pnA) | EXFOR INC-SOURCE |
| Irradiation | 2.0 h | EXFOR METHOD |
| Target | Metallic Ho foils, 99% pure, 12.29 µm (10.81 mg/cm² at 8.795 g/cm³), 15 × 15 mm | EXFOR SAMPLE |
| Stack | Ho foils interleaved with Cu (9.71 µm, 8.70 mg/cm²) and Ti (10.40 µm, 4.69 mg/cm²) as degraders and monitors; the order and number of foils are in the paper, not in EXFOR | EXFOR METHOD |
| Monitor | natTi(α,x)⁵¹Cr, 26.4 mb at 50 MeV, the value the authors used | EXFOR MONITOR |
| Detection | Off line, HPGe, after cooling: ¹⁶⁵Tm 242.917 keV (35.5%, 30.06 h, 0.7–1.8 d); ¹⁶⁶Tm 778.814 keV (19.1%, 7.70 h, 7.7–20.7 h); ¹⁶⁷Tm 207.801 keV (42%, 9.25 d, 9.5–18 d); ¹⁶⁸Tm 815.989 keV (50.95%, 93.1 d, 19.3–42.6 d) | EXFOR DECAY-DATA |
| Errors | Beam intensity 5%, efficiency 4%, thickness 2%, γ intensity 0.3–19%, statistics 0.1–8%; total per point in `uncertainty_mb` | EXFOR ERR-ANALYS |

How to run it in physim: one experiment per Ho foil (beam ⁴He at the foil's `e_lab_mev`, target ¹⁶⁵Ho 10.81 mg/cm²),
the channel table from physim's own fusion (`Fusion("4He", "165Ho").thick_target`) and evaporation Monte Carlo
(`physim.nuclear.evaporation`), and the residue count per channel from σ × luminosity over 2.0 h at 87 pnA. The
off-line decay counting is outside the app: compare at the level of production cross sections (or, equivalently,
residues made per foil). The stack itself (energy loss through Cu and Ti) is not needed, because the paper gives the
energy at each Ho foil; a check that physim's stopping reproduces those energies from 50 MeV is a useful extra.

## The numbers to compare

All of `evaporation_4He_165Ho_literature.csv` with `ref = Usman2020`. At a glance (σ in mb, measured ± error):

| E_lab (MeV) | 1n ¹⁶⁸Tm | 2n ¹⁶⁷Tm | 3n ¹⁶⁶Tm | 4n ¹⁶⁵Tm | physim σ_fus (Bass/Wong) |
|---|---|---|---|---|---|
| 25.6 | 13.2 ± 0.9 | 837 ± 168 | 0.7 ± 0.1 | — | 1038 |
| 35.1 | 8.85 ± 0.60 | 162 ± 33 | 1070 ± 99 | 10.7 ± 1.3 | 1544 |
| 42.6 | 4.83 ± 0.33 | 87 ± 18 | 561 ± 52 | 708 ± 58 | 1729 |
| 49.6 | 1.91 ± 0.14 | 36 ± 7 | 151 ± 15 | 840 ± 68 | 1680 |

What the comparison should report, with tolerances fixed before running (backlog 75's): each channel's peak energy
within the energy spread plus 1 MeV; the strongest channel at each energy the same; the strongest three channels
within 30% where they exceed 10% of the total; and σ(3n)/σ(4n) at 45–50 MeV, the most sensitive number (below).

What is already known before running (from the fusion module and the PACE4 comparison, no new physics):

- **Fusion.** The measured Tm residues add up to 0.7–0.9 of physim's σ_fus from 25 to 48 MeV, which is plausible
  since (α,pxn), (α,αxn) and 5n are not measured; above 50 MeV the sum falls (0.6 at 49.6 MeV, 0.4 at 59 MeV in
  Kormazeva's data) because 5n (¹⁶⁴Tm, 2 min) takes over and was not measured. So this set cannot check σ_fus on
  its own.
- **Channel shares at 50 MeV.** Measured σ(3n)/σ(4n) is 0.18 (Usman, 49.6 MeV), 0.16 (Kormazeva, 49.7 MeV) and 0.27
  (Gadkari, 50.2 MeV). PACE4 gives 0.057 and physim 0.048 (`docs/theory/evaporation.md`: ¹⁶⁵Tm 0.90/0.85, ¹⁶⁶Tm
  0.051/0.041). Both statistical codes are 3–5× low on 3n at 50 MeV. The usual explanation, and the subject of the
  Gadkari and Mukherjee papers, is pre-equilibrium emission of a fast first neutron, which neither code has. Expect
  physim to match the channels near their peaks and to miss the high-energy tails of 2n and 3n; the tails are the
  documented limit, as they are for PACE4.
