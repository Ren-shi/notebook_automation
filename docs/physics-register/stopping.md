# Stopping power, range, energy loss and straggling

Backlog item 36 · module `physim.nuclear.stopping` · tests `tests/python/test_nuclear_stopping.py` · theory
{doc}`../theory/stopping`

## Model

- **Protons and α particles:** NIST PSTAR/ASTAR tables (ICRU 49/90) for 74 materials. Isotopes are scaled at the
  same velocity.
- **Other elements:** interpolated per electron between tabulated neighbours, linearly in ln I.
- **Other compounds:** Bragg additivity.
- **Heavy ions:** effective-charge scaling of proton stopping (Brandt–Kitagawa with Ziegler's ionisation fit). At
  low velocity, the larger of that and Lindhard–Scharff joined to full-charge stopping.
- **Nuclear stopping:** NIST's tables for p and α in NIST materials, ZBL universal otherwise.
- **Ranges and energy after a layer:** by integrating 1/S. Tilted layers are crossed along t/cos θ.
- **Straggling:** Bohr, with the relativistic factor and Lindhard–Scharff's slow-ion reduction, carried through
  thick layers.
- **Multiple scattering:** Highland's formula.

This choice was made with the user (2026-10-04): own models plus the ICRU 49 tables, rather than formulas alone
(10–20% errors near the Bragg peak) or wrapping an external library.

## Assumptions and range of validity

- Above about 1 MeV/u everything is good to a few percent for p and α, and about 10% for heavy ions.
- Below that, accuracy depends on the ion (see the validation table). Slow heavy ions (Kr, Xe below 1 MeV/u) are
  uncertain by tens of percent, and so are established codes: LISE++'s own three models differ from each other by
  up to 48% for Xe there.
- Untabulated elements near and below the Bragg peak (about 0.1–0.5 MeV/u for protons) can be off by 10–40%:
  their shell structure is not captured.
- Gas targets use NIST's gas-phase values. Channelling, charge-state distributions and target damage are ignored.

## Validation

Tolerances are the measured worst cases rounded up, over 8 elemental targets (C, Al, Si, Ni, Cu, Ag, Ta, Au, two of
them interpolated) and 7 ions.

| Check | Reference | Tolerance | Result | Test |
|---|---|---|---|---|
| p and α ranges, 2 keV to 100 MeV, in C, Si, Cu, Au, Mylar, water, polyethylene | NIST CSDA ranges (tabulated independently of our integration) | 1% | pass (typically < 0.2%) | `test_integrated_range_matches_nist_csda_range` |
| Range–energy consistency, tilted layers, length thicknesses | Exact relations | 1e-6 | pass | `test_energy_after_uses_the_range` |
| Isotope velocity scaling; full stripping at high velocity | Exact scaling; Z₁² S_p | 0.1–1% | pass | `test_isotopes_scale_with_velocity`, `test_heavy_ions_are_fully_stripped_at_high_velocity` |
| Bragg additivity (CD₂) | Hand sum | 0.05% | pass | `test_bragg_additivity_for_a_formula` |
| Interpolation for untabulated elements, p ≥ 0.5 MeV | Leave-one-out against NIST (Al, Cu, Ag, Au, Pb) | 7% | pass (worst 6.8%) | `test_interpolation_for_untabulated_elements` |
| Bohr straggling, thin layer; Highland; X₀(Si) | Formulas; PDG 21.82 g/cm² | 0.2%; 0.1%; 2.5% | pass | `test_bohr_straggling_thin_layer`, `test_highland_multiple_scattering` |
| Stopping power vs LISE++ (ATIMA 1.4), 1680 points | LISE++ library (`stopping_lise.csv`) | p/α: 4% ≥ 10 MeV/u, 7% 1–10, 13% 0.1–1. C–Ar: 4%, 11%, 25%. Kr/Xe: 10%, 12%, 55% | pass | `test_against_lise_atima` |
| Range vs LISE++ | same | p/α: 3.5%, 6%, 21%. C–Ar: 5%, 15%, 39%. Kr/Xe: 9%, 53%, 100% | pass | `test_against_lise_atima` |
| Energy after a fifth of the range vs LISE++ | same | p/α: 0.5%, 1.2%, 5.3%. C–Ar: 0.8%, 4.8%, 11.5%. Kr/Xe: 2.2%, 6%, 19% | pass | `test_against_lise_atima` |
| Straggling (FWHM) ≥ 30 MeV/u vs LISE++ | same | 25% | pass | `test_straggling_against_lise_at_high_energy` |
| Stopping power vs SRIM 2013 (SR Module), 5 ions (⁴He, ¹²C, ¹⁶O, ⁴⁰Ar, ¹³²Xe) × Au, Si × 0.1–10 MeV/u, 50 points | SRIM (`stopping_srim.csv`) | the LISE++ stopping bands above | pass (worst 90% of a band: Xe in Au at 1 MeV/u, −10.7% against 12%) | `stopping_vs_srim` |
| Stopping power vs SRIM 2013, compounds: same ions × Mylar, CD₂, CsI, 75 points | same; compound correction 1 (Bragg additivity) | the same bands (no compound allowance; see below) | **fail**: CD₂ ¹⁶O/⁴⁰Ar 0.1 MeV/u +31%/+39% (25%), ⁴⁰Ar 1 MeV/u −14% (11%), ¹⁶O 10 MeV/u −4.2% (4%), Xe 1 MeV/u −20% (12%); Mylar ⁴⁰Ar 0.1 MeV/u +25.4% (25%). CsI and all ⁴He pass | `stopping_compounds_vs_srim` (a known failure in `test_nuclear_validation.py`) |
| Range vs SRIM | same | not judged | SRIM gives the projected range, physim the path length (CSDA); they differ by up to a factor 2 for slow ions in Au, and by 1–12% at 10 MeV/u | reported in the checks' detail |

The LISE++ values are generated on a machine with LISE++ installed by `scripts/lise_reference.py`; CI reads the
committed file. Below 0.1 MeV/u nothing is tested against LISE++: there, the codes differ from each other by tens
of percent (LISE++'s Hubert/Ziegler/ATIMA by up to 25% for protons and 48% for Xe).

SRIM's values are generated on a machine with SRIM 2013 installed by `scripts/srim_reference.py` (SR Module, no
window). SR Module's output header reads "Calculation using SRIM-2006 / SRIM version ---> SRIM-2012.01": the 2013
installer ships an SR Module built in 2011 whose stopping coefficients (`SCOEF03.dat`) are byte-identical to the
full program's. Every compound was run with compound correction 1 (plain Bragg additivity), as physim's TRIM.IN
export also writes: SR Module takes the factor as an input, and the Core-and-Bond factor SRIM's GUI would apply to
its dictionary entry "Mylar" is computed only by the GUI. physim's Mylar uses the NIST/ICRU proton and α tables,
which include the chemical binding.

**Tolerance for compounds.** The register states no separate tolerance for compounds, and none is proposed: the
model page claims the same accuracy for Bragg-additive compounds as for elements, so the compounds are judged by the
same bands. They fail them. The disagreement is in heavy ions in hydrogen-rich targets (CD₂ most, Mylar at the
low-energy edge), not in α particles (within 8% everywhere) or in CsI. **Where the CD₂ gap comes from (tested 2026-11-14).**

- *Not the Fermi velocity.* physim's effective charge uses v_F = 1 v₀ for every target. SRIM's own values
  (`SCOEF03.dat`, column 9) are H 1.031 and C 1.000, and putting SRIM's per-element v_F into physim moves CD₂ by at
  most 0.3% (Si and Mylar as little; CsI up to 4.3% and Au up to 3.6%, where v_F is far from 1).
- *Hydrogen.* SR Module runs for the same ions in pure C and pure H show that SRIM's CD₂ is exactly its Bragg sum
  of C and H (within 0.05%). physim agrees with SRIM's carbon within about 6% (C, O, Ar; slow Xe up to −27%), but
  not with its hydrogen: +79% to +97% at 0.1 MeV/u, −25% to −35% at 1 MeV/u, −6% to −12% at 10 MeV/u.
- *SRIM treats hydrogen targets differently.* The effective charge SRIM implies, ζ² = S_ion / (Z₁² S_p), matches
  physim's in carbon within a few percent. In hydrogen it departs from SRIM's own carbon value: 0.10 against 0.16
  for O at 0.1 MeV/u, 0.90 against 0.61 at 1 MeV/u, and above 1 (1.09) at 10 MeV/u. So SRIM does not scale its
  proton stopping in hydrogen as it does in other targets.

**Against measurements (2026-11-14).** Measured electronic stopping from the IAEA Stopping Power Database (version
2026-09, read through its API; the database's own "SRIM2013" curve rows dropped) for C, O, Ar, Kr and Xe at
0.05–15 MeV/u, compared point by point with physim and with SR Module (compound correction 1; H₂ and D₂ as gases).
Not yet an automated check:

| Target (measurements) | Points | RMS deviation, physim | RMS deviation, SRIM | What the data say |
|---|---|---|---|---|
| H₂ and D₂ gas, 1–3 MeV/u: Ar and Kr (Bimbot et al. 1989), one C point (1962) | 19 | 23% | 8% | physim too low: Ar −17%, Kr −27% (median). SRIM's special treatment of hydrogen is supported |
| H₂ and D₂ gas, C, Ar, Kr, 3–15 MeV/u | 63 | 11.0% | 9.4% | Kr: physim −16%, SRIM −11% |
| (CH₂)ₙ: polypropylene, polyethylene, C and O ions (6 data sets) | 206 | 10.3% | 9.4% | both about 6–7% low on average. Near 0.1 MeV/u the data lie between the codes: physim +9 to +14%, SRIM −9 to −14% |
| Mylar, C to Xe (15 data sets) | 460 | 8.0% | 4.9% | Kr and Xe: physim −14% to −25%, SRIM about −10% |

One polyethylene set (Diwan et al. 2003, 7 points) lies 20–40% below both codes and the polypropylene data, and is
left out of the (CH₂)ₙ row. There are no measurements for pure hydrogen below 1 MeV/u in this range, and none for
Ar or Xe in CH₂ polymers, so physim's largest disagreements with SRIM in CD₂ (+39% for Ar at 0.1 MeV/u, −20% for Xe
at 1 MeV/u) are not tested directly.

So for CD₂: at 1–3 MeV/u, physim's heavy-ion stopping in hydrogen is too low, by roughly the 15–25% it differs from
SRIM there. At low energy, neither code matches the polymer data to better than about 10–14%. Until physim treats
hydrogen targets better, plan with heavy ions in CD₂ or polyethylene targets with care, and check the result against
SRIM.

Validated 2026-10-04; SRIM comparison 2026-11-14.

## Known limitations

- Heavy ions in hydrogen-rich compounds (CD₂, Mylar) differ from SRIM by up to 20% at 1–3 MeV/u and up to 40%
  at 0.1 MeV/u, beyond the stated bands (see the validation table).
- Slow heavy ions (below 1 MeV/u, especially Kr and heavier): stopping within ~50%, range within a factor of two.
  Plan with care and check against SRIM. Ziegler's per-ion and per-target low-energy coefficient tables would fix
  this; they are not included (a possible later download).
- Straggling below ~10 MeV/u differs from LISE++ by up to a factor 2–4 (slow-ion and charge-exchange straggling
  are not modelled).
- Untabulated elements near the Bragg peak (see above).
- Implemented in Python with NumPy (about 0.1 s to build a table per ion and material, then vectorised); the event
  generator (item 39) may move the hot path to Rust.

## References

- M. J. Berger et al., NIST Standard Reference Database 124 (PSTAR, ASTAR), doi:10.18434/T4NC7P; ICRU Reports 49
  (1993) and 90 (2016).
- J. F. Ziegler, J. P. Biersack and U. Littmark, *The Stopping and Range of Ions in Solids* (Pergamon, 1985).
- W. Brandt and M. Kitagawa, Phys. Rev. B **25**, 5631 (1982).
- J. Lindhard and M. Scharff, Phys. Rev. **124**, 128 (1961); Mat. Fys. Medd. Dan. Vid. Selsk. **27**, no. 15 (1953).
- N. Bohr, Mat. Fys. Medd. Dan. Vid. Selsk. **18**, no. 8 (1948); C. Tschalär, Nucl. Instr. Meth. **61**, 141 (1968).
- V. L. Highland, Nucl. Instr. Meth. **129**, 497 (1975); O. I. Dahl, as given by the Particle Data Group.
