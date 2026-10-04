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
| Compounds, slow heavy ions, TRIM straggling | SRIM/TRIM | — | 🟡 request written (`tests/reference/nuclear/pending/stopping_srim.md`) | — |

The LISE++ values are generated on a machine with LISE++ installed by `scripts/lise_reference.py`; CI reads the
committed file. Below 0.1 MeV/u nothing is tested against LISE++: there, the codes differ from each other by tens
of percent (LISE++'s Hubert/Ziegler/ATIMA by up to 25% for protons and 48% for Xe).

Validated 2026-10-04.

## Known limitations

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
