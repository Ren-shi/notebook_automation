# Atomic masses and material data

Backlog item 34 · module `physim.nuclear.data` · tests `tests/python/test_nuclear_data.py`

## Model

- **Masses.** Every nuclide's mass excess Δ is read from the AME2020 table, shipped unchanged. Atomic mass
  M = A·u + Δ, with u = 931 494.102 42 keV. Nuclear mass M_nuc = M − Z·mₑ + B_e(Z), where the total electron binding
  energy is B_e(Z) = 14.4381 Z^2.39 + 1.55468×10⁻⁶ Z^5.35 eV (Lunney, Pearson and Thibault 2003).
- **Q-values.** Q = ΣΔ(entrance) − ΣΔ(exit) − E*, from atomic mass excesses. Electron masses cancel because charge
  is conserved; electron binding energies cancel to within a few eV for light projectiles.
- **Natural elements.** NIST isotopic abundances weight the AME2020 atomic masses, so molar masses and isotope
  compositions are consistent with the masses used in kinematics.
- **Materials.** An element, isotope, formula or NIST compound is turned into nuclides per formula unit. Areal
  densities follow from ρ·t; atoms per cm² from (ρt)·N_A / M_formula. A pure isotope gets the natural element's
  density scaled by atomic mass (same number of atoms per volume).

## Assumptions and range of validity

- Neutral-atom masses are exact to the AME2020 uncertainties. Masses marked `#` are estimates from systematics,
  flagged by `Nuclide.estimated`.
- B_e(Z) is a fit good to a few tens of eV; it matters only for nuclear masses of heavy ions.
- Natural abundances are "representative" values; real materials vary slightly, which matters only for isotope
  ratios, not for planning count rates.
- Densities are bulk values (NIST graphite: 1.70 g/cm³). Thin foils often differ, so give `density` when the
  thickness is a length and the difference matters.

## Validation

| Check | Reference | Tolerance | Result | Test |
|---|---|---|---|---|
| Q-values of 8 reactions and decays (D–T, D–D both branches, D–³He, ⁶Li(n,t), deuteron binding, ²³⁸U α decay, ¹²C(α,γ)) | Tabulated values, quoted to 1 keV (deuteron binding to 1 eV) | 1 keV | all within 0.5 keV | `test_q_values_match_published_to_1_kev` |
| α-particle and proton nuclear masses | CODATA 2018 | 0.2 keV (α), 1 eV (p) | pass | `test_masses_in_other_units` |
| Standard atomic weights of all 84 natural elements with Z ≤ 92 | NIST standard atomic weights | inside the published interval or 1σ | 82 pass; Se 1.5σ, Hg 2.4σ (see limitations) | `test_atomic_weights_from_isotopes_agree_with_nist` |
| Areal densities and atoms/cm² for Au, C, Si, CH₂, CD₂ | Hand calculation | 1e-4 relative | pass | `test_gold`, `test_carbon_and_silicon`, `test_polyethylene_and_deuterated_polyethylene` |
| Every data file has its source and licence recorded | `data/SOURCES.md` | — | pass | `test_every_data_file_has_its_source_recorded` |

Validated 2026-10-04.

## Known limitations

- Se and Hg: NIST's isotopic compositions (2009) give atomic weights 1.5σ and 2.4σ from the standard atomic weights
  revised in 2013. The effect on any planning quantity is below 10⁻⁴.
- Isotopic enrichment is available from Python (`Material.enriched`) but not yet in the setup file.
- Densities stop at uranium (Z = 92), and NIST's densities for At and Fr are placeholders.
- No nuclear level data yet (needed from backlog item 43).

## References

- W. J. Huang et al., Chinese Physics C **45**, 030002 (2021); M. Wang et al., Chinese Physics C **45**, 030003
  (2021): AME2020.
- J. S. Coursey et al., *Atomic Weights and Isotopic Compositions with Relative Atomic Masses*, NIST (2015).
- J. H. Hubbell and S. M. Seltzer, *Tables of X-Ray Mass Attenuation Coefficients and Mass Energy-Absorption
  Coefficients*, NIST Standard Reference Database 126 (2004), Tables 1 and 2.
- D. Lunney, J. M. Pearson and C. Thibault, Rev. Mod. Phys. **75**, 1021 (2003): electron binding energies.
