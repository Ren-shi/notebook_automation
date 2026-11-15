# Fusion-evaporation: residues from a channel table

`physim.nuclear.residues` plans a fusion-evaporation experiment from a table of exit channels (backlog item 75,
Stage 1). A heavy-ion beam above the Coulomb barrier fuses with a target nucleus into a compound nucleus. The
compound nucleus boils off neutrons, protons and α particles and leaves a residue, whose γ-ray cascade a germanium
array records. Stage 1 takes each channel's cross section from the setup (typed in, from PACE4 in LISE++ or from a
measurement) and works out the rest: the compound nucleus's excitation and speed through the target, each event's
particles and the residue's recoil, the residue's γ rays in the crystals with their Doppler shift, and rates and
beam times per channel. The app's own fusion cross section (Stage 2) and evaporation code (Stage 3) will fill the
same channel table.

## The setup

```toml
[reaction]
type = "fusion-evaporation"
alignment = 0.3           # σ/J of the residues' substates about the beam (0.3 if left out)
side_feeding = "25 %"     # share of the cascades entering below the top yrast level (25 % if left out)
temperature = "1.5 MeV"   # optional: the evaporation temperature of every channel that gives none

[[reaction.channels]]
residue = "158Dy"         # optional: worked out from the compound nucleus, checked if given
particles = "4n"          # "4n", "p2n", "a n", "αn" ...
cross_section = "520 mb"
entry_spin = 12           # optional, ħ
temperature = "1.4 MeV"   # optional
source = "PACE4"          # optional, in words; "placeholder" makes the planner warn
```

The compound nucleus is the beam plus the target's most abundant nuclide. Each residue may have its level scheme
as `[levels.158Dy]`, in the same form as `[levels.beam]`; without one the local ENSDF copy is read (levels up to
6 MeV). The example `fusion_dy158` is ¹²C at 60 MeV on 1 mg/cm² of ¹⁵⁰Nd, making ¹⁶²Dy at 50 MeV of excitation.
Its cross sections are illustrative values of the size PACE4 gives, not a measurement.

## The compound nucleus

The excitation energy is the invariant mass of beam and target less the compound nucleus's mass, with the nuclear
masses of AME2020:

$$
\sqrt{s} = \sqrt{(m_b + m_t)^2 + 2 m_t T_b}, \qquad E^* = \sqrt{s} - M_{CN} = E_{cm} + Q,
\qquad E_{cm} = \sqrt{s} - m_b - m_t .
$$

The compound nucleus moves along the beam at

$$
\beta_{CN} = \frac{p_b}{T_b + m_b + m_t}, \qquad p_b = \sqrt{T_b^2 + 2 m_b T_b} .
$$

This is exact. The non-relativistic forms E_cm = T_b m_t / (m_b + m_t) and v_CN = v_b m_b / (m_b + m_t) are 10 keV
and 4 × 10⁻⁵ c off for the example, more than the 1 keV and 10⁻⁶ c the checks ask for.

**Nuclear, not atomic, Q.** Q from the atomic mass excesses, as `data.q_value` gives it, also contains the binding
energy the electrons gain when the compound atom forms: 68 keV for ¹⁶²Dy from ¹²C and ¹⁵⁰Nd. That energy belongs to
the electron cloud, not the nucleus, so E* uses the nuclear masses (`Nuclide.nuclear_mass_mev`). For the example
Q = −5.567 MeV and E* = 49.977 MeV at the front of the target.

**Through the target.** Each event happens at a depth drawn evenly through the target. The beam's mean energy
there comes from the stopping powers (see {doc}`stopping`); its spread is the beam's own spread carried through the
layer plus straggling. E* falls from 49.98 MeV at the front to 48.96 MeV at the back of the example's target, and
β_CN from 0.007687 to 0.007616.

## The evaporation

Each event picks a channel with probability σ_i / Σσ and so stands for the rate Σσ × L / N, with L the luminosity
(beam particles per second × target nuclei per cm²). The channel's particles leave one after the other. In the
frame of the emitting nucleus a particle leaves isotropically, with relative kinetic energy ε drawn from the
evaporation spectrum

$$
P(\varepsilon) \propto (\varepsilon - V)\, e^{-(\varepsilon - V)/T}, \qquad \varepsilon > V ,
$$

V being the Coulomb barrier Z_x Z_d e² / (1.5 fm (A_x^⅓ + A_d^⅓)) for a charged particle (0 for a neutron). T is the
channel's temperature, the reaction's, or √(E*/a) of the emitting nucleus with a = A/8 MeV⁻¹. The mean neutron energy
in the example is 3.0 MeV, that of the protons 12.7 MeV and that of the α particles 22.7 MeV.

Energy is conserved exactly. The daughter keeps E* − S − ε, S being the separation energy from the masses. Each
particle leaves at least what the ones after it need (their separation energies and barriers), so a channel that
the compound nucleus can afford always completes. The residue must end unable to emit another particle: its E*
below its lowest separation energy plus barrier, for neutrons, protons and α particles. A chain that ends hotter is
drawn again, up to 20 times; after that its last particle is drawn from the tail of the spectrum above what it must
carry away, and the notes say how many events needed it. A channel that needs this often is one the compound
nucleus has too much energy for: in the example, 3n at 50 MeV.

Momentum is conserved particle by particle: the two-body decay is solved in the emitter's frame with the masses
including their excitation, and both products are boosted to the laboratory. The residue's direction and speed are
the sum of the recoils. Its mean squared angle to the beam is about Σ (2/3) 2μ⟨ε⟩ / P²: 6.4° r.m.s. for 4n in the
example, 17.6° for α2n.

## The residue's γ rays

**Entry.** A residue enters its known levels by spin. The top entry is the lowest level of the highest spin not above
the entry spin (the yrast entry). A share `side_feeding` of the cascades enters instead at one of the lower yrast
levels, evenly. The entry spin is the table's, or two thirds of the grazing angular momentum ℓ_max of a sharp
cut-off, ℓ_max(ℓ_max + 1)ħ² = 2μR²(E_cm − V_B), less 1 ħ per nucleon and 3 ħ per α particle emitted.

**Cascade.** The cascade follows the level scheme's branching ratios. A transition with internal conversion α emits
a γ ray in 1/(1 + α) of its decays.

**Alignment.** The compound nucleus's spin is perpendicular to the beam, so its substates m about the beam are mostly
small. The entry level gets the distribution exp(−m²/2σ²) with σ = `alignment` × J, and its statistical tensors
B_k = t_k0 / t_00. Each γ ray has the angular distribution

$$
W(\theta) = 1 + B_2 A_2 P_2(\cos\theta) + B_4 A_4 P_4(\cos\theta)
$$

about the beam, in the residue's frame, with A_k of the transition's multipolarities and mixing. The level below
keeps B_k U_k. Down a cascade of stretched E2 transitions B_2 A_2 stays the same, as the yrast lines of real
residues show: a₂ = +0.31 for J = 10 and σ/J = 0.3.

**Detection.** The γ rays leave where the residue was made, at the speed it was made with, and are boosted to the
laboratory (energy and direction). The crystal each meets and what it records come from the response of the γ-ray
detectors (see {doc}`detectors`).

**Doppler correction.** As an experiment without a recoil detector does it, the correction assumes every residue
moves along the beam at β_CN for the nominal beam energy at the middle of the target:
E₀ = E γ (1 − β cos θ_crystal). What is left of a line's width at 90° is

$$
\sigma_E = E_0\, \beta \sqrt{\langle\psi_x^2\rangle + \alpha^2/4} ,
$$

ψ_x being the residue's angle towards the crystal (the residues' opening angle) and α the crystal's half-angle
(the spread of a disc in one direction is α/2). For the 320.5 keV 6⁺ → 4⁺ line of ¹⁵⁸Dy in a 70 mm crystal at
200 mm, the simulated width is 0.287 keV against 0.290 keV from this formula.

## Rates and the Plan tab

A channel makes σ × L residues per second; for the example's 4n channel that is 6.5 × 10⁴ /s. The Plan tab's table of
residues gives, per channel, its cross section and share, the residues per second and in the run, the threshold E*
it needs, and its three strongest lines: their γ rays per residue (from the feeding and branching) and their
full-energy counts per second in all γ-ray detectors together. The beam time is for the counts wanted in the
channel's strongest line. A channel the compound nucleus cannot afford at the middle of the target is shown in
red. The "Predicted spectra" check shows the residues' simulated γ-ray spectra in each detector, measured and
corrected.

## Limits

- The cross sections are the table's, the same at every depth of the target. Stage 2 adds σ_fus(E) and its average
  over the target.
- The evaporation spectrum is a fixed shape with one temperature (or √(E*/a)), not the Hauser–Feshbach competition
  of Stage 3. The particles are not correlated with the residue's spin.
- No γ rays above the yrast line (the statistical cascade) and none from the compound nucleus before the first
  particle.
- Lifetimes, slowing down of the residues, summing and add-back of their γ rays are left out: every γ ray is emitted
  where the residue was made.
- The particles are not followed into the particle detectors; those see the elastic scattering as for any setup.
- The runs of the Run tab hold the particle detectors' events; the residues' γ rays are simulated by the planner,
  not stored with a run.
- No PACE4 reader in the planner: a channel table from PACE4 is typed in. `scripts/pace4_reference.py` runs PACE4
  in silent mode and reads its HTML output, but only to write the validation suite's reference files.

## Checks

| Check | Reference | Tolerance | Test |
|---|---|---|---|
| E* and β_CN for ¹²C + ¹⁵⁰Nd, ¹⁶O + ¹⁴⁸Nd, ⁴⁸Ca + ¹²⁴Sn | Invariant mass worked out by hand from AME2020; E_cm + Q with Q from the atomic mass excesses and the electron binding | 1 keV, 1e-6 c | `test_excitation_and_speed_equal_the_hand_calculation` |
| E* through the target | Compound nucleus at the beam energy after the target, from `Stopping.energy_after` | 2 keV | `test_excitation_through_the_target_follows_the_stopping` |
| Residues per channel | σ × luminosity | 4σ (total exact) | `test_counts_per_channel_equal_the_table_times_the_luminosity` |
| Momentum and energy, event by event | The compound nucleus's four-momentum | 1e-8 relative | `test_momentum_and_energy_are_conserved_event_by_event` |
| The residues' opening angle | Σ (2/3) 2μ⟨ε⟩ / P² | 3% | `test_the_residues_opening_angle_follows_from_the_particles_momenta` |
| Doppler-corrected width at 90° | E₀ β √(⟨ψ_x²⟩ + α²/4) | 5% (achieved 1%) | `test_doppler_corrected_width_equals_the_geometry` |
| Cascade yields | Feeding and branching, 1/(1 + α) | 1e-9 (analytic), 4σ (Monte Carlo) | `test_cascade_yields_and_feeding` |
| Alignment | B₂ < 0; a₂ = B₂ A₂ the same down a stretched E2 cascade | 1e-9 | `test_alignment_perpendicular_to_the_beam` |
| Seeded runs | The same seed twice | bit-identical | `test_seeded_runs_are_bit_identical` |
| Against PACE4 | Channel cross sections and particle spectra (Stage 3), five systems | three strongest in order and within 30%; ⟨ε_n⟩, ⟨ε_p⟩ 10% | `evaporation_channels_vs_pace4`, `evaporation_spectra_vs_pace4` (known failures, see {doc}`evaporation`) |
| Against GEMINI++ | The same quantities, three systems (GEMINI++ hangs on ³²S and ⁴⁸Ca), physim from GEMINI++'s E* and spin distribution | the same; residues tied within 2 standard errors of GEMINI++'s statistics may swap places | `evaporation_channels_vs_gemini`, `evaporation_spectra_vs_gemini` (known failures) |
| ⟨ε_n⟩ against both codes | PACE4 and GEMINI++ disagree by 25–35%; physim from PACE4's E* and σ_J | inside their band, widened by 10% at each edge | `evaporation_neutron_energy_band` |
| Measured residue cross sections | ¹⁶⁵Ho(α,xn) 1n–4n, 6n (Usman 2020, Kormazeva 2022, Gadkari 1997); physim σ_fus × share at each energy | 2n, 3n, 4n where each is the strongest measured channel: within 30% of the band of the three data sets ± their errors; tails reported | `residues_vs_data` |
| Measured fusion | ¹⁶O + ¹⁴²,¹⁵⁰Nd (Visakh 2021), ¹⁶O + ¹⁵⁴Sm (Stokstad 1980) | 20% from 1.05 / 1.15 V_B to 1.3 V_B | `fusion_vs_data` (see {doc}`fusion`) |
| A published experiment, end to end | Usman et al. 2020 (⁴He + ¹⁶⁵Ho stacked-foil activation): each Ho foil as the `fusion_tm165` example, channels from *Fill from physim*, residues simulated, σ inferred from the residues made | above 1.1 V_B: peaks of 2n, 3n, 4n within the energy spread + 1 MeV; strongest channel the same where the data resolve it; strongest three within 30% of the band of the three measurements up to each peak; tails and σ(3n)/σ(4n) reported | `reproduction_usman2020` (see {doc}`evaporation`, "Reproduced experiment") |

The first and third rows are also checks of the validation suite (`compound_nucleus_hand_calculation`,
`residue_counts_vs_luminosity` in `physim.nuclear.validation`), which the physics register reads, as are the last
six rows.

**Measured residue cross sections.** ⁴He + ¹⁶⁵Ho is the well-studied system of item 75: three independent
stacked-foil activation measurements of the Tm residues, from EXFOR (`evaporation_4He_165Ho_literature.csv`). They
disagree by up to 1.5× at the (α,4n) peak and in its position (45–50 MeV), so `residues_vs_data` judges physim
against the band they span, each set interpolated in energy and widened by its own error. It judges only what
physim claims to model, the compound nucleus's evaporation: each xn channel where it is the strongest measured
channel (2n at 19.6–31 MeV, 3n at 31.7–41.8 MeV, 4n at 42.6–49.7 MeV; 4n only up to its measured peak, since
(α,5n) was not measured and takes over above it). Above each peak pre-equilibrium emission, which physim and PACE4
leave out, feeds the channels with fewer neutrons: physim's 3n at 50 MeV is 0.4–0.5 of the measured value. At the
35 judged points physim is 1.1–1.6× the Usman values (up to 2.3× Gadkari's, which are the lowest); the worst
point is 21% above the band (4n at 46.2 MeV), inside the 30% tolerance.
