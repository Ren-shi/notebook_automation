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

**Charge states of Kr in H₂, He and N₂ (2026-11-15).** Explanation (a) tested against measured equilibrium charge
states, which are independent of the stopping data. Kuboki et al. (2014) measured ⁸⁶Kr at 2.7 MeV/u, inside
Bimbot's range: mean charge 25.1 in H₂, 23.2 in He, and less than 20 in N₂ (estimated). Gas-stripper work at GSI
also finds the highest equilibrium charge of U at 1.4 MeV/u in H₂ (29.2 ± 1.2; Scharrer et al. 2017). The charge
the stopping needs was found by inverting physim's Brandt–Kitagawa ζ for q, gas by gas, from Bimbot's data (linear
fit over 1.5–5.3 MeV/u, read at 2.7 MeV/u; scatter about the fit 0.3–0.5):

| | H₂ | D₂ | He | N₂ | Ar |
|---|---|---|---|---|---|
| Charge the stopping needs | 23.3 | 23.1 | 22.4 | 20.8 | 19.8 |
| Measured charge (Kuboki et al.) | 25.1 | | 23.2 | < 20 | |
| physim (ZBL, the same for every target) | 20.8 | 20.8 | 20.8 | 20.8 | 20.8 |
| physim's stopping ÷ data − 1 | −20% | | −13% | 0% | |
| … with the measured charge put in its place | +15% | | +7% | ≤ −7% | |

- *Direction and order agree:* the measured charge is highest in H₂, then He, then N₂, the same order as the charge
  the stopping needs and as physim's shortfall. So the charge state is at least part of the cause.
- *The size does not:* the measured differences (H₂ − N₂ > 5.1, He − N₂ > 3.2) are at least twice those the stopping
  needs (2.5 and 1.6), and putting the measured charge into physim's formula overshoots. ZBL's q is a quantity
  fitted to stopping in solids, not a measured mean charge, so it cannot simply be replaced in one target.
- *Not settled:* why only about half of the measured charge difference shows up in the stopping. A charge-state
  model with a target dependence, used for every target with the same screening, has to be tested against all of
  Bimbot's gases and the solids before anything changes in the code. For U at 1.4 MeV/u, He gives a lower charge
  than N₂ (Gerhard et al. 2022), so the order between He and N₂ depends on the ion.

**Charge-state formulas tested (2026-11-15; nothing changed in the code).** Published mean-charge formulas for
⁸⁶Kr at 2.7 MeV/u, against Kuboki's measurement and the charge Bimbot's stopping needs:

| Target | Measured | Needed | Schiwietz–Grande 2001, gas | Schiwietz–Grande 2001, solid | Schiwietz et al. 2004, solid | Winger et al. 1992 | ZBL (physim) |
|---|---|---|---|---|---|---|---|
| H₂ | 25.1 | 23.3 | 22.0 | 25.9 | 26.4 | 33.2 | 20.8 |
| He | 23.2 | 22.4 | 22.0 | 25.7 | 26.1 | 31.2 | 20.8 |
| N₂ | < 20 | 20.8 | 22.1 | 25.2 | 25.5 | 27.8 | 20.8 |
| Ar | | 19.8 | 22.1 | 24.9 | 25.0 | 26.0 | 20.8 |

- *Schiwietz–Grande, gas:* reproduces U in N₂ at 1.4 MeV/u (28.4 against 27), but has almost no target dependence
  and puts H₂ below N₂ for U (27.4 against 28.4; measured 28–29 against 27). Used as physim's q for Kr in Bimbot's
  gases: H₂ −22% → −12%, but N₂ +4% → +16% and Ar +7% → +20%. Rejected.
- *The solid formulas* give the charge of ions leaving a foil (25–26, above ZBL's 20.8, which is fitted to stopping)
  and change by only about 1 unit from H to N. They do not describe gases.
- *Winger et al.* has a target dependence of about the measured size (H − N 5.4, He − N 3.4, against > 5.1 and > 3.2
  measured), at solid-like levels. Used directly, stopping is 52–84% too high. Used only as a factor on physim's
  charge, normalised to N₂, Bimbot's gases go from −22%…+7% to −6%…+11%, but heavy solids break: Kr and Xe in Ni,
  Ag and Au go from within 4% to 13–24% low, and Ar ions likewise. Rejected.

No published formula gives the measured gas charges, and none whose target dependence carries into the
Brandt–Kitagawa formula works for gases and solids together. This agrees with the finding above that the stopping
follows only about half of the measured charge differences: the charge state points the right way, but the step
from charge to stopping is where the model is missing something.

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
- R. Bimbot et al., heavy-ion stopping in gases (1989), via the IAEA Stopping Power Database (ref_id Bi1989b).
- H. Kuboki et al., Phys. Rev. ST Accel. Beams **17**, 123501 (2014) (⁸⁶Kr charge states in H₂ and He).
- P. Scharrer et al., Phys. Rev. Accel. Beams **20**, 043503 (2017); P. Gerhard and M. Maier, HIAT 2022, MO3I3
  (heavy-ion charge states in H₂, He and N₂ gas strippers).
- G. Schiwietz and P. L. Grande, Nucl. Instr. Meth. B **175–177**, 125 (2001); G. Schiwietz et al., Nucl. Instr.
  Meth. B **225**, 4 (2004); J. A. Winger, B. M. Sherrill and D. J. Morrissey, Nucl. Instr. Meth. B **70**, 380
  (1992) (mean-charge formulas).
- J. Lindhard and M. Scharff, Phys. Rev. **124**, 128 (1961); Mat. Fys. Medd. Dan. Vid. Selsk. **27**, no. 15 (1953).
- N. Bohr, Mat. Fys. Medd. Dan. Vid. Selsk. **18**, no. 8 (1948); C. Tschalär, Nucl. Instr. Meth. **61**, 141 (1968).
- V. L. Highland, Nucl. Instr. Meth. **129**, 497 (1975); O. I. Dahl, as given by the Particle Data Group.
