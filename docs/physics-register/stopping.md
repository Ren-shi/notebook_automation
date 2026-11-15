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
- **Ions from Li to Ar in hydrogen and in ICRU 49 compounds** (since 2026-11-15): MSTAR 3.12 (Paul and Schinner,
  the method of ICRU Report 73), `physim.nuclear.mstar`. It scales ICRU 49 α stopping by fits to measured data:
  S = S_α(4E) · F(Z₁, E) · Z₁²/4, with separate fits for solids, for H₂ and He gas, and per ion for other gases.
  Used for targets made only of hydrogen isotopes (H, D₂), for the ICRU 49 compounds (Mylar, polyethylene,
  water, CsI, air, …), and for a deuterated version of one of them (CD₂ is polyethylene, scaled by electrons per
  gram). Elements (including graphite), other compounds and heavier ions keep the effective-charge model. Above
  50 MeV/u the two are blended in ln E, reaching the effective-charge model at 200 MeV/u (they differ by < 2%
  there).
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
| Stopping power vs SRIM 2013, compounds: same ions × Mylar, CD₂, CsI, 75 points | same; compound correction 1 (Bragg additivity) | the same bands (no compound allowance; see below) | **fail, CD₂ only** (with MSTAR, 2026-11-15): ¹²C/¹⁶O 0.1 MeV/u +26%/+27% (25%), ¹²C/¹⁶O/⁴⁰Ar 1 MeV/u −16%/−15%/−17% (11%), Xe 1 MeV/u −20% (12%). Mylar (worst Xe 1 MeV/u −11.8%), CsI and all ⁴He pass. Before MSTAR: CD₂ ⁴⁰Ar 0.1 MeV/u +39%, Mylar ⁴⁰Ar 0.1 MeV/u +25.4% | `stopping_compounds_vs_srim` (a known failure in `test_nuclear_validation.py`) |
| Range vs SRIM | same | not judged | SRIM gives the projected range, physim the path length (CSDA); they differ by up to a factor 2 for slow ions in Au, and by 1–12% at 10 MeV/u | reported in the checks' detail |
| Electronic stopping vs measurements: C, O, Ar, Kr, Xe in H₂, D₂, (CH₂)ₙ polymers, Mylar, 0.1–15 MeV/u, 744 points | IAEA Stopping Power Database (`stopping_measured_hydrogen.csv`, `scripts/iaea_stopping_reference.py`) | median deviation per target group, ion and energy band inside the LISE++ stopping band above (the existing bands, not set from these data) | **fail** in 3 of 21 cells: Kr in H₂/D₂ 1–10 MeV/u −19.0% (12%, 34 points), Xe in Mylar −17.5% and Kr in Mylar −14.5% at 1–10 MeV/u (12%, 3 points each), all outside MSTAR's range. RMS for C–Ar with MSTAR (before; SRIM): H₂/D₂ 4.1% (7.9%; 7.5%), (CH₂)ₙ 8.5% (10.3%; 9.4%), Mylar 5.9% (6.9%; 4.6%). SRIM by the same rule fails only C in H₂ (+10 to +14%) |
| MSTAR port | MSTAR 3.12's own values for ¹²C in water (solid fit with target correction) and air (gas fit for C), 0.07–10 MeV/u, as extracted from the original Fortran in libdEdx's tests | 0.1% | pass (≤ 0.03%) | `test_heavy_ions_in_icru49_compounds_follow_mstar` |
| MSTAR handoff above 50 MeV/u; CD₂ as polyethylene | Exact: MSTAR below 50 MeV/u, unchanged above 200; electrons per gram 0.8745 | 1e-12; 0.1% | pass | `test_mstar_hands_over_smoothly_at_high_energy`, `test_deuterated_compound_stops_like_its_parent_per_electron` | `stopping_hydrogen_vs_data` (a known failure) |

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
The data and SRIM's value at each point are in `tests/reference/nuclear/stopping_measured_hydrogen.csv`, judged by
`stopping_hydrogen_vs_data` (from 0.1 MeV/u, see the table above); the RMS figures below include 0.05–0.1 MeV/u:

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

**The fix: MSTAR for ions from Li to Ar (2026-11-15).** Splitting physim's old stopping into its carbon and
hydrogen parts showed the carbon part within +7% to −6% of 237 measured points for C and O in carbon, and
the hydrogen part 6–60% low. Three home-made fixes (scaling from α, a factor fitted to the H₂/D₂ data, the
NIST compound proton table) did not hold up; the fitted factor failed a held-out test (3.8% → 6.9% RMS on
the data it was not fitted to). MSTAR, ported from its Fortran (`mspaul.f`, `mstar1.f`, `msid2z.f`) and not
refitted, improves every hydrogen group it covers (RMS against the measured data, ions C–Ar, from
0.1 MeV/u):

| Target | Points | physim before | physim with MSTAR | SRIM 2013 |
|---|---|---|---|---|
| H₂/D₂ gas | 48 | 7.9% | 4.1% | 7.5% |
| (CH₂)ₙ polymers | 206 | 10.3% | 8.5% | 9.4% |
| Mylar | 439 | 6.9% | 5.9% | 4.6% |

Caveats:

- MSTAR was fitted to measurements up to 2002, which likely include some of these (Bimbot et al. 1989 for Ar
  in H₂/D₂), so this is a published fit judged against overlapping data, not an independent test.
- The port reproduces MSTAR's own numbers for water and air to 0.03%. Its H₂/He coefficients have no
  reference value; they were checked line by line against the Fortran and libdEdx's C port.
- At 1 MeV/u, C, O and Ar in CD₂ are now 15–17% below SRIM (before 12–15%); the polymer data there side with
  SRIM (MSTAR runs 5–7% low for C and O in (CH₂)ₙ at 1–10 MeV/u). Near 0.1 MeV/u they side with MSTAR.
- Kr and Xe are outside MSTAR (Z₁ ≤ 18) and unchanged: in H₂/D₂ and Mylar they are still 15–19% low, against
  SRIM's 10%. They need another source.

**Kr and Xe in hydrogen: candidates (2026-11-15, nothing changed in the code).** Median / RMS of model ÷ measured − 1:

| Data | Points | physim | DPASS 21.06 | SRIM-96 (ZBL) | SRIM 2013 |
|---|---|---|---|---|---|
| Kr in H₂/D₂, 1–3 MeV/u (Bimbot et al. 1989, as all Kr points in H₂/D₂) | 12 | −27% / 27% | −28% / 28% | −22% / 22% | −8% / 9% |
| Kr in H₂/D₂, 3–15 MeV/u (data stop at 5.2 MeV/u) | 22 | −16% / 17% | −21% / 21% | −18% / 19% | −11% / 12% |
| Xe in H₂, 25 MeV/u (He1991, errors 4.7–8.6%) | 2 | +4% | −7% | −7% | −4% |
| Ar in H₂/D₂, 1–15 MeV/u (for comparison; physim uses MSTAR) | 39 | −2% / 4% | −13% / 14% | −15% / 16% | −6% / 6% |

DPASS (Sigmund and Schinner's binary theory) and SRIM-96 (the published ZBL algorithm) come from the tables in
ir2-lab/libdedx. MSTAR's H₂/He fit cannot be extrapolated beyond Ar (it gives factors of 50–700 too high). Only
SRIM 2013, whose treatment of hydrogen targets is not published, matches the Kr data. Xe in H₂ has data only at
25 MeV/u, where physim is already within the errors.

The Kr and Xe failures in Mylar are not a hydrogen problem: replacing only physim's hydrogen term changes them by
under 1%, while physim is low in carbon too (Xe in carbon: −20% at 0.1–1 MeV/u, −7% at 1–10 MeV/u, against 212
measured points; DPASS +5% and 0%, SRIM-96 −9% and −2%).

**Why physim is low for Kr in hydrogen: the facts (2026-11-15).** physim's effective-charge model alone (Kr is
outside MSTAR) against the IAEA database, median of model ÷ measured − 1 per data set:

1. *One experiment.* All 34 Kr points in H₂ and D₂ are Bimbot et al. 1989 (Bi1989b), 1.5–5.2 MeV/u. There are none
   between 5.2 and 21 MeV/u.
2. *The same experiment measured Kr in other gases* at the same energies, with the same beam and method:

   | Gas | H₂ | D₂ | He | N₂ | Ar |
   |---|---|---|---|---|---|
   | physim ÷ measured − 1 | −20% | −16% | −3% | +6% | +7% |

   Errors common to the experiment (beam energy, gas thickness) cannot make this spread. H₂ and D₂ should agree per
   electron; their 4% difference is the experiment's scatter, against 23–27% from H₂ to N₂ and Ar.
3. *Solids show the same order, across many groups* (0.5–10 MeV/u, range of the per-source medians):

   | Ion | C | Al | Ni | Au |
   |---|---|---|---|---|
   | Kr | −12% to +1% | −9% to +2% | −8% to 0% | −8% to +11% |
   | Xe | −17% to −8% (5 of 6 sources near −17%) | −11% to −3% | −2% to +1% | −10% to +11% |

4. *Near full stripping the order remains* (22–26 MeV/u, He et al. 1991): Kr H₂ +4%, He +7%, N₂ +6%, Ar +9%;
   Xe H₂ +4%, He +7%, N₂ +14%, Ar +23%.
5. *It grows with the ion's charge:* in H₂ at 2–9 MeV/u, Ar is −6% (Bi1989a) and Kr −20%.
6. *As an effective charge,* the data imply 2–3 more units of charge on Kr in H₂ than physim gives:

   | E (MeV/u) | 1.5 | 1.9 | 3.1 | 4.3 | 22 |
   |---|---|---|---|---|---|
   | ζ, physim | 0.48 | 0.52 | 0.61 | 0.69 | 0.94 |
   | ζ, from data | 0.57 | 0.60 | 0.69 | 0.75 | 0.92 |

7. *physim cannot make a target dependence.* ζ is ZBL's fit, made mostly to solids, with the target Fermi velocity
   set to 1 for every target; only the proton stopping changes with the target. MSTAR has a separate fit for H₂ and
   He gas for this reason, which is why ions from Li to Ar now pass.
8. *Ruled out:* the proton table (PSTAR's H₂ is used as published); the gas–solid effect (ions carry less charge in
   gases, which would make physim too high, not too low); the target Fermi velocity (tested for CD₂); an isotope
   effect (H₂ and D₂ are both low); scatter in the data (item 2).
9. *DPASS is as low* (−28% and −21%, table above) although it has target-dependent stopping physics (Bloch,
   Barkas and shell corrections). This points at its charge state for H₂ rather than at the stopping theory; it
   does not prove it.

Two explanations fit these facts. (a) *Charge state:* Kr may stay more stripped in H₂ than ZBL's solid fit says,
since H₂ has two slow electrons per molecule and captures poorly onto a fast ion. This can be checked against
measured equilibrium charge states in H₂ and N₂ (gas-stripper studies), which are independent of the stopping
data. (b) *Stopping per unit charge:* for Kr here Z₁v₀/v = 2.5–4.6, so close collisions follow Bohr's classical
regime, which lowers the stopping by a roughly fixed amount: a larger fraction for heavy targets than for hydrogen.
A rough estimate gives H₂ about 15% above N₂ relative to the proton scaling, the right direction and size, but it
also predicts Ar about 15% below N₂, which the data do not show. So (b) alone is not the whole story. Using SRIM
2013's tables for these ions was considered and rejected: SRIM 2013 is fitted to data that include Bi1989b, does
not publish its treatment of hydrogen, and would make the SRIM checks circular. (a) is to be tested first.

Validated 2026-10-04; SRIM comparison 2026-11-14; MSTAR 2026-11-15.

## Known limitations

- Heavy ions in CD₂ differ from SRIM by up to 20% at 1 MeV/u and up to 27% at 0.1 MeV/u, beyond the stated bands
  (see the validation table). Kr and Xe in hydrogen-rich targets are 15–19% below measured data at 1–10 MeV/u.
- The effective charge of ions heavier than Ar does not depend on the target, but measured stopping says it
  does: below ~10 MeV/u physim is up to 20% low in H₂ and C and a few percent high in N₂ and Ar (see "Why physim
  is low for Kr in hydrogen").
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
