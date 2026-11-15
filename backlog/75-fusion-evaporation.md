# 75 · Fusion-evaporation: the compound nucleus, its decay and the residues' γ rays

**Priority:** P3 · **Size:** L (three stages, each L; split into items when work starts) · **Area:** Nuclear planner
(third slice)

**Status: Built, validation open** (2026-11-15). All three stages and the detectors are in the app. Left:
- Stage 3 misses its targets for heavy beams at high spin (³²S + ⁹²Mo, ⁴⁸Ca + ¹²⁴Sn are off by up to 10× against
  PACE4) and gives neutrons 10–25% softer than PACE4; the cause found so far is the spin removed along the chain
  (PACE4's statistical E2 γ rays).
- Stage 2 against PACE4's excitation functions (five systems, 8–10 energies each, `fusion_vs_pace4`): within 15%
  above 1.1 V_B only for ¹²C + ¹⁵⁰Nd; ⁴He + ¹⁶⁵Ho 16% high at 1.14 V_B, and the heavier beams up to 43% low at
  2 V_B, where physim's capture stops at the Bass potential's critical ℓ and PACE4's keeps rising. Against measured
  excitation functions (`fusion_vs_data`, ¹⁶O + ¹⁴²,¹⁵⁰Nd and ¹⁵⁴Sm): within 20% from 1.05 V_B (spherical) or
  1.15 V_B (deformed) up to 1.3 V_B, worst 9.6%; above 1.3 V_B ¹⁴²Nd is up to 1.6× high (fusion hindrance, which
  Visakh et al. 2021 fit with a larger diffuseness), below the barrier 3–500× low (coupled channels).
- In the validation suite since 2026-11-15: 11 checks over four register rows, PACE4 runs saved by
  `scripts/pace4_reference.py`. The exact and internal checks pass (✅ in the register's literature column); the
  three PACE4 comparisons fail their item-75 tolerances and are known failures: σ_fus (above), the three strongest
  channels (no system passes; ⁴He and ¹²C fail on the third channel), ⟨ε_n⟩ (10–24% low; ⟨ε_p⟩ passes).
- Wave 2 of the suite (2026-11-15), six more checks (17 in all): GEMINI++ 3.3.0 for three systems (it hangs on ³²S and ⁴⁸Ca) fails
  the same targets as PACE4 (the strongest three residues for no system, shares off by up to 65%; ⟨ε_n⟩ 10–17%
  high, ⟨ε_p⟩ passes) and is a known failure too. GEMINI++'s neutrons are softer than physim's and PACE4's harder:
  physim's ⟨ε_n⟩ lies inside the band of the two codes (they differ by 25–35%), which `evaporation_neutron_energy_band`
  checks (band widened by 10%). Measured ¹⁶⁵Ho(α,xn) residue cross sections (`residues_vs_data`, three data sets
  from EXFOR) pass: 2n, 3n, 4n where each is the strongest channel are within 30% of the band the data span (worst
  21%; physim 1.1–1.6× the Usman 2020 values), while the tails above each peak, fed by pre-equilibrium emission,
  are reported (3n at 50 MeV 0.4–0.5 of measured). The ENDF/B-VIII.0 carbon check is registered (`carbon_vs_endf`,
  literature).
- Wave 3 (2026-11-15): a published experiment reproduced end to end through the app, the check
  `reproduction_usman2020` (18 checks in all). Usman et al., NIM B 469, 42 (2020) (⁴He + ¹⁶⁵Ho stacked-foil
  activation at RIKEN, 26 Ho foils, 1n–4n): each foil as the new example `fusion_tm165` at the foil's energy, the
  channel table from *Fill from physim*, the residues simulated and σ inferred from them
  (`physim.nuclear.reproduction`, `scripts/reproduce_usman2020.py`). Above 1.1 V_B the peaks of 2n, 3n, 4n lie
  within the energy spread + 1 MeV (27.6, 36.9, 47.6 MeV against 25.8, 37.6, 47.1), the strongest channel is right
  at all 20 foils where the data resolve it, and the strongest three channels are within 30% of the band of the
  three measurements up to each peak (worst +15%; physim/Usman 0.6–1.6). PACE4 at the same energies misses two of
  the peaks (2n +2.4, 4n +3.1 MeV) and the strongest channel at 3 foils. Reported as the documented limit: the
  tails above the peaks (pre-equilibrium), σ(3n)/σ(4n) at 49.6 MeV 0.055 against 0.16–0.27 measured. physim's
  stopping through the stack order inferred from the energies gives each foil's energy within its spread.

Done for the neutron detectors (2026-11-15):
- The carbon cross sections are now ENDF/B-VIII.0 (natural carbon from ¹²C and ¹³C, bin-averaged), checked by
  `tests/reference/nuclear/neutron_carbon_endf.csv` (`scripts/endf_neutron_reference.py`). The old typical table
  was 31% low at 3.4 MeV and 24% high at 8.9 MeV, so cell efficiencies move by up to −12% (3.4 MeV) and +6%
  (8.9 MeV). Gammel's n–p formula agrees with ENDF/B-VIII.0 H-1 within 0.8% from 0.5 to 20 MeV.
- Channel selection sorts events into detected classes (0n, 1n, p2n, …) rather than requiring every particle, and
  now applies the particle detectors' thresholds (each detector's `threshold`, else `[neutron_wall]
  particle_threshold`, default 1 MeV, with dead layer and punch-through) and γ rays firing neutron cells
  (`gamma_multiplicity`, default J/2 + 4 from the entry spin; `gamma_efficiency`, `gamma_as_neutron`). It returns
  the analytic matrix of the same events beside the Monte Carlo; the two agree within statistics. Still left out:
  energy loss in the target, the particle detectors' own p/α identification, γ-ray cross-talk.

## Why
Fusion-evaporation is how most nuclei away from stability, and most high-spin states, are made: a heavy-ion beam
above the Coulomb barrier fuses with the target, the compound nucleus boils off neutrons, protons and α particles,
and the residue it leaves emits a cascade of γ rays. Planning such an experiment needs numbers the app cannot give
today: how much of each residue is made, at what beam energy and target thickness, how fast the residues move and
in which direction (the Doppler correction), how many γ rays come per event, and how well the light-particle and
neutron detectors tell one exit channel from another. Students now take these from PACE4 inside LISE++ and stitch
them together by hand.

## What is already there
- Masses and Q-values (34), two-body kinematics (35), stopping and straggling (36): the compound nucleus's
  excitation energy, the recoil velocity and the energy loss through the target.
- The Coulomb barrier at the interaction radius (`Rutherford.coulomb_barrier_cm`, 37), a first estimate of the
  fusion threshold.
- Level schemes from ENSDF (50), γ-detector response (53), orientation and angular distributions (54), Doppler
  correction and coincidences (55), add-back and Compton suppression (62), cascades with γ–γ correlation (71),
  and the lifetimes and summing work in `gamma_events` (72, 73).
- The run model (63) and its tabs, which read events and do not care which reaction made them.

## Scope

### Stage 1: residues from a channel table
The first useful version uses cross sections it does not compute.
- `[reaction] type = "fusion-evaporation"`, added to `REACTIONS`. The setup names the compound nucleus (from beam
  and target) and lists the exit channels, each with a cross section: typed in, or read from a PACE4 output file.
- Exact kinematics of the compound nucleus: E* = E_cm + Q, recoil velocity v_CN = v_beam m_b / (m_b + m_t) along
  the beam, and how both change through the target's thickness (the stopping, 36).
- Each event picks a channel by its cross section and emits the light particles with an isotropic CM distribution
  and an evaporation (Maxwellian) energy spectrum whose temperature the table gives or a default sets. The
  residue's direction and speed follow from momentum conservation, particle by particle; this sets the Doppler
  cone.
- The residue's γ rays: the cascade (71) starts from an entry state of given spin and energy, fed into the known
  levels by spin (yrast and side feeding), with the compound nucleus's spin aligned perpendicular to the beam (54).
- Rates and beam time per channel (39), and a Plan table of residues: cross section, rate, the strongest lines.

### Stage 2: the app's own fusion cross section
- The fusion barrier from the Bass potential and the capture cross section from Wong's formula (a parabolic barrier
  summed over partial waves): σ_fus(E) and its partial-wave distribution σ_ℓ = πƛ²(2ℓ+1)T_ℓ.
- The excitation function averaged over the target's thickness, shown against beam energy so the user can pick
  the energy.
- Warnings: below the barrier (fusion falls by orders of magnitude and coupled-channels effects matter), above
  about 10 MeV/u (incomplete fusion), and heavy systems where quasi-fission competes.

### Stage 3: the app's own evaporation Monte Carlo
- Statistical decay of the compound nucleus in the Hauser–Feshbach picture, one random chain per event as PACE4
  does: at each step neutron, proton, α, γ (and fission for heavy nuclei) compete by their widths.
- Level densities from the back-shifted Fermi gas (or Gilbert–Cameron), transmission coefficients from a global
  optical model or a fitted parametrisation, a rigid-rotor yrast line, giant-dipole E1 for statistical γ rays,
  Bohr–Wheeler for fission.
- Outputs per event: the channel, each particle's energy and angle, the residue's momentum, and the entry state
  (E*, J) that Stage 1 then feeds into the cascade. The channel table of Stage 1 becomes a result instead of an
  input, and either can be used.
- On CPU cores: chains differ in length and each step is a choice, which suits threads better than a GPU.

### Detectors for channel selection
- Neutron detectors as a new kind: liquid or plastic scintillator, an efficiency against energy, neutron–γ
  discrimination, scattering between neighbouring detectors (cross-talk).
- A charged-particle ball (silicon or CsI) is already possible with the existing detectors; what is new is
  counting how many protons and α particles each event shows.
- The channel-selection matrix: for each true channel, the share of events assigned to each channel from the
  detected multiplicities. This is the number an experimentalist designs around.
- A recoil separator only as an acceptance (a transmission number per residue); its ion optics are left out.

## Left out
- Incomplete fusion, pre-equilibrium emission and deep-inelastic reactions.
- Quasi-fission and fission fragments beyond counting fission as a lost share.
- Coupled-channels fusion below the barrier (CCFULL would be the outside check if it is added).
- The separator's optics and focal-plane detectors.
- Prompt γ rays from the compound nucleus before the first particle.

## Verification
- PACE4 (in LISE++) for the channel cross sections and the particle spectra; GEMINI++ as a second opinion. The app
  runs PACE4 if the installed version has a silent mode, otherwise the comparison is by hand.
- Published: one measured fusion excitation function across the barrier, and measured residue cross sections for
  one well-studied system (to choose).
- References: C. Y. Wong, Phys. Rev. Lett. 31, 766 (1973); R. Bass, *Nuclear Reactions with Heavy Ions* (1980);
  W. Hauser and H. Feshbach, Phys. Rev. 87, 366 (1952); A. Gilbert and A. G. W. Cameron, Can. J. Phys. 43, 1446
  (1965); A. Gavron, Phys. Rev. C 21, 230 (1980) (PACE).
- Each capability gets its row in the physics register.

## Depends on
34–39, 50, 54, 55, 63, 71. Stage 1 also benefits from 72 and 73 (high γ multiplicity makes summing and pile-up
matter; fast residues decay in flight).

## Done when
- **Stage 1:** E* and the recoil velocity equal the hand calculation within 1 keV and 1e-6 c; with a channel table,
  the residues counted per channel equal table × luminosity within statistics; the Doppler-corrected width of a
  residue line equals what the opening angle and the detector size give, within statistics; seeded runs are
  bit-identical.
- **Stage 2:** above the barrier, σ_fus for the chosen system agrees with PACE4 and the measured excitation function
  within 15% (the suite's measured check allows 20%: this 15% combined with the data's 10–15% errors); far above, it approaches the geometric limit; the thick-target average equals the integral over the
  energy loss within 1%.
- **Stage 3:** for the chosen system, the strongest three channels agree with PACE4 in order and within 30%; the
  mean neutron and proton energies within 10%; turning off charged-particle emission gives Weisskopf's neutron
  spectrum.
- **Detectors:** the channel-selection matrix from the event generator agrees with the binomial expectation from the
  efficiencies within statistics.
