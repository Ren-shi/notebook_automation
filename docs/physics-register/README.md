# Physics register

What physics Physim can simulate today, how far each capability has got, and the evidence that it is right. Written
for anyone auditing the project: students deciding whether to trust a number, supervisors, reviewers.

Physim brings together calculations usually done with LISE++, SRIM and similar tools. It does not replace them.
Every capability is checked against those tools first, then against published measurements, and this register
records where each one stands.

## Status levels

| Column | Meaning |
|---|---|
| **Engine** | Implemented and usable from Python. |
| **App** | Available in the experiment planner web app ({doc}`../planner-app`), without writing code; every tab is rendered by `tests/python/test_nuclear_app.py`. |
| **vs tools** | Agrees with SRIM, LISE++ or another established code within the stated tolerance. |
| **vs literature** | Agrees with published data or analytic results within the stated tolerance. |

✅ done · 🟡 in progress · ⬜ planned · — not applicable

## Rules

1. **Every ✅ in a validation column points to an automated test** that checks the stated tolerance and runs in CI.
   If the test fails, the entry is no longer ✅. A claim without a test is 🟡 at best.
2. **Tolerances are stated, not implied.** "Agrees with SRIM" is not an entry; "agrees with SRIM within 3% for
   α in Au from 1 to 50 MeV" is.
3. **Limits are listed as clearly as capabilities.** Each page says where its model stops being valid.
4. **Reference data has provenance**: tool and version (or paper and table), exact inputs, date (see backlog item 40).
5. **The table is checked by a test.** `physim.nuclear.validation` holds every comparison behind this table, and
   `tests/python/test_nuclear_validation.py` fails if a ✅ here is not backed by passing checks, or if a deliberately
   broken formula slips through. `python -c "from physim.nuclear import validation; print(validation.report())"`
   prints the current state, and the {doc}`validation notebooks <../validation>` plot every comparison.

## Nuclear experiment planner

First slice: elastic scattering (backlog items 33–42).

| Capability | Engine | App | vs tools | vs literature | Page | Backlog |
|---|---|---|---|---|---|---|
| Experiment setup file | ✅ | ✅ | — | — | [Guide](../nuclear-setup.md) | 33 |
| Atomic masses and material data | ✅ | ✅ | ✅ LISE++ masses, 46 nuclides, within 2 keV or the AME2020 uncertainty | ✅ Q-values within 1 keV; atomic weights vs NIST | [Data](nuclear-data.md) | 34 |
| Two-body reaction kinematics | ✅ | ✅ | 🟡 LISE++: 1 case passes, 8 requested | ✅ analytic limit, Lorentz boost, threshold | [Kinematics](kinematics.md) | 35 |
| Stopping power and range | ✅ | ✅ | 🟡 LISE++ (1680 points, tolerances by regime) and SRIM for Au, Si (50 points, same bands) pass; SRIM compounds: Mylar and CsI pass, CD₂ fails (up to 27% at 0.1 MeV/u, 20% at 1 MeV/u) | 🟡 NIST CSDA ranges within 1% pass; measured slow heavy ions (Ar–U) in elemental solids (IAEA, 3620 points) pass, medians within 10% per ion class and band; measured heavy-ion stopping in H₂/D₂, (CH₂)ₙ and Mylar (IAEA, 744 points) fails 3 of 21 cells, all Kr and Xe (outside MSTAR): Kr in H₂/D₂ −21%, Kr and Xe in Mylar −12% and −16% at 1–10 MeV/u | [Stopping](stopping.md) | 36 |
| Energy and angular straggling | ✅ | ✅ | ✅ LISE++ ≥ 30 MeV/u (25%); TRIM widths under Monte Carlo spectra | ✅ Bohr and Highland formulas | [Stopping](stopping.md) | 36 |
| Rutherford cross section | ✅ | ✅ | ✅ LISE++ (1 case, 0.03%) | ✅ analytic, Geiger–Marsden 1913 | [Rutherford](rutherford.md) | 37 |
| Distance of closest approach and validity checks | ✅ | ✅ | ✅ LISE++ grazing angle (1 case) | 🟡 above-barrier data requested | [Rutherford](rutherford.md) | 37 |
| Coulomb trajectories | ✅ | ✅ | — exact closed form instead | ✅ analytic b(θ) to 1e-6° | [Rutherford](rutherford.md) | 37 |
| Detector solid angles and response | ✅ | ✅ | — exact closed forms instead | ✅ closed forms to 1e-10, Monte Carlo | [Detectors](detectors.md) | 38 |
| Count rates and beam time | ✅ | ✅ | 🟡 LISE++ cross sections at the detector angles requested | ✅ I n (dσ/dΩ) Ω within 1%; Monte Carlo within 4σ | [Rates and events](rates-and-events.md) | 39 |
| Monte Carlo spectra | ✅ | ✅ | 🟡 TRIM, 4 layers: 2 within mean 1% and σ 20%; α in Mylar mean +1.05%, slow ¹³²Xe σ −63%. All 6 exported TRIM.IN files run unchanged | ✅ peak positions and widths vs analytic (5%); seeded, thread-independent | [Rates and events](rates-and-events.md) | 39 |
| Beam-time report, CSV and ROOT export | ✅ | ✅ | 🟡 opening the ROOT file in ROOT requested (`pending/root_check.md`); ✅ uproot read-back | ✅ key numbers vs hand calculations (`test_nuclear_report.py`) | [Guide](../nuclear-setup.md) | 42 |

Second slice: Coulomb excitation (backlog item 43).

| Capability | Engine | App | vs tools | vs literature | Page | Backlog |
|---|---|---|---|---|---|---|
| Coulomb excitation, first-order semiclassical | ✅ | ✅ | ✅ GOSIA 20110524.14, first order (weak-coupling limit), ¹⁶O + ⁵⁸Ni at 30 MeV, five angles, E1 polarisation in both (2%; worst 0.23%); the E1 polarisation correction, blind, eight cases at `DIP,` 3.8 and 20: 47 probabilities (1% above 10⁻⁴, 3% above 10⁻⁶; worst 0.48%) and 112 changes by it (2% of the change; worst 0.5%); the coupled equations, blind, with feeding and reorientation: ⁶²Ni (Q_s = ±0.3 e b, 0⁺₂, 2⁺₂, 4⁺₁) on ¹²C and ¹⁹⁴Pt, ¹⁵²Sm on ¹²C, 48 probabilities (1% above 10⁻⁴, 3% above 10⁻⁶; worst 0.46%) | 🟡 closed forms at 180° (1e-6), engine-orbit integration (5e-8) and adiabatic cutoff pass; the reproduced experiment (Allmond 2014, ⁶²Ni on carbon, below) passes; measured data requested | [Theory](../theory/coulex.md) | 43 |
| γ-ray Doppler shift and broadening | ✅ | ✅ | — exact formula instead | ✅ Lorentz boost (1e-12) | [Theory](../theory/coulex.md) | 43 |
| Coulomb-excitation rates and spectra | ✅ | ✅ | — | ✅ Monte Carlo against analytic rates per channel within 4σ (`coulex_ni58`) | [Theory](../theory/coulex.md) | 43 |
| γ-ray yields with deorientation in vacuum | ✅ | ✅ | 🟡 GOSIA 20110524.14, OP,POIN at a point, blind: every line of ⁶²Ni and ¹⁹⁴Pt (⁶²Ni at 237 MeV on ¹⁹⁴Pt) and ⁹⁸Ru (on ²⁴Mg, the recoil detected), three angles each, ten point γ detectors, GOSIA's parameters and a given g: 1320 on/off ratios (2% of the change where it is above 1%, else 1%) and 651 yields without deorientation relative to each other (1%), 6 missed, all in the 2⁺₂ → 0⁺ lines (the two codes' Coulomb excitation of the 2⁺₂ level, not the deorientation; worst 3.2 bands, where that line's yield is 3% of its mean); `test_nuclear_gosia_yields` | — | [Theory](../theory/coulex.md) | 43 |

Third slice: fusion-evaporation (backlog item 75). The checks are in `physim.nuclear.validation`, with the PACE4
runs saved in `tests/reference/nuclear/pace4_*.csv` by `scripts/pace4_reference.py`, the GEMINI++ runs in
`gemini_*.csv` by `scripts/gemini_reference.py`, measured data from EXFOR in `fusion_*_literature.csv` and
`evaporation_4He_165Ho_literature.csv`, and ENDF/B-VIII.0 carbon in `neutron_carbon_endf.csv`. The reproduced experiment is `reproduction_usman2020.csv` (`scripts/reproduce_usman2020.py`), with PACE4 at the same energies in `pace4_excitation_4He_165Ho.csv`. The three PACE4
comparisons and the two GEMINI++ ones fail their tolerances and are listed as known failures in
`test_nuclear_validation.py`.

| Capability | Engine | App | vs tools | vs literature | Page | Backlog |
|---|---|---|---|---|---|---|
| Fusion-evaporation residues from a channel table | ✅ | ✅ | — exact kinematics instead | ✅ E* and β vs the hand calculation from AME2020 masses (1 keV, 1e-6 c; exact), residues per channel vs σ × luminosity (4 standard errors; worst 0.8); Doppler-corrected width vs geometry in `test_nuclear_fusion_evaporation.py` | [Theory](../theory/fusion-evaporation.md) | 75 |
| Fusion cross section (Bass barrier, Wong) | ✅ | ✅ | 🟡 PACE4 4.36.0 excitation functions, 5 systems × 8–10 energies, 15% above 1.1 V_B: ¹²C + ¹⁵⁰Nd passes (7%); ⁴He + ¹⁶⁵Ho 16% at 1.14 V_B; ¹⁶O, ³²S, ⁴⁸Ca beams 24%, 42%, 43% low at 1.75–2 V_B (physim's capture stops at the Bass potential's critical ℓ) | ✅ Bass barriers vs Dutt and Puri 2010 (0.03 MeV, 0.03 fm), Wong's closed form vs the partial-wave sum (1%), classical limit (1e-6), thick-target average vs the integral (1%); measured ¹⁶O + ¹⁴²,¹⁵⁰Nd (Visakh 2021) and ¹⁶O + ¹⁵⁴Sm (Stokstad 1980) within 20% from 1.05 V_B (¹⁴²Nd) or 1.15 V_B (deformed) to 1.3 V_B (worst 9.6%; ¹⁵⁴Sm has no point there). Reported: 3–500× low below the barrier (coupled channels); ¹⁴²Nd 1.1–1.6× high at 1.33–1.57 V_B (fusion hindrance) | [Theory](../theory/fusion.md) | 75 |
| Evaporation Monte Carlo (Hauser–Feshbach) | ✅ | ✅ | 🟡 PACE4 4.36.0, 5 systems, physim from PACE4's E* and σ_J: the strongest three residues in order and within 30% for none (⁴He, ¹²C: the two strongest within 20%, the third is 5n or 6n instead; ³²S off by up to 40×, ⁴⁸Ca by 2×); ⟨ε_p⟩ within 10% for all five (worst −9%), ⟨ε_n⟩ 10–24% low. GEMINI++ 3.3.0, 3 systems (it hangs on ³²S, ⁴⁸Ca): the strongest three in order and within 30% for none (shares off by up to 65%); ⟨ε_p⟩ within 10% (worst +9%), ⟨ε_n⟩ 10–17% high. ✅ ⟨ε_n⟩ inside the PACE4–GEMINI++ band widened by 10% (the codes differ by 25–35%), all three | ✅ Weisskopf's neutron spectrum (χ²/dof 1.1 < 2), energy and momentum conserved in every event; measured ¹⁶⁵Ho(α,xn) (Usman 2020, Kormazeva 2022, Gadkari 1997): 2n, 3n, 4n where each is the strongest channel within 30% of the band of the three data sets (35 points, worst 21%); tails above the peaks (pre-equilibrium) reported: 3n at 50 MeV 0.4–0.5 of measured; the reproduced experiment (Usman 2020 through the planner, below) passes | [Theory](../theory/evaporation.md) | 75 |
| Neutron detectors and channel selection | ✅ | ✅ | — closed forms instead (GEANT4 would be the tool) | ✅ efficiency vs Knoll's closed form (1e-12), Monte Carlo selection matrix vs the binomial expectation (4.5 standard errors; worst 1.9), carbon total cross section vs ENDF/B-VIII.0 (0.2% below 1.8 MeV, 7% / 5% of the ±2.5% / ±10% averages above; worst 91% of the tolerance) | [Theory](../theory/neutrons.md) | 75 |

**Reproduced experiments** (the referee's first question: does the app reproduce a published experiment end to
end?). The check `reproduction_usman2020` (literature, under the evaporation row) reruns it; the result is saved in
`tests/reference/nuclear/reproduction_usman2020.csv` by `scripts/reproduce_usman2020.py`, and drawn in the
validation notebook *Reproduced experiment*. The check `reproduction_allmond2014` (literature, under the Coulomb-excitation
row) reruns the Coulomb-excitation one from the `coulex_ni62_c_allmond` example (`scripts/reproduce_allmond2014.py`);
it was a known failure, for one band it missed, until physim had the E1 polarisation (2026-11-25).

| Experiment | Through the app | Rows it exercises | Status |
|---|---|---|---|
| A. R. Usman et al., Nucl. Instrum. Methods B 469, 42 (2020), EXFOR E2650: ⁴He + ¹⁶⁵Ho stacked-foil activation at RIKEN, ¹⁶⁵Ho(α,xn)Tm, 1n–4n at 26 Ho foils from 49.6 to 13.2 MeV (87 pnA, 2.0 h, 10.81 mg/cm²) | Each foil as the `fusion_tm165` example at the foil's mean energy; channel table from *Fill from physim*; residues simulated; σ inferred from the residues made with the paper's beam, foil and time | Residues from a channel table, fusion, evaporation (not the neutron detectors: the residues were counted off line) | ✅ above 1.1 V_B: peaks of 2n, 3n, 4n within the energy spread + 1 MeV (27.6, 36.9, 47.6 MeV against 25.8, 37.6, 47.1); strongest channel right at all 20 foils where the data resolve it; the strongest three channels within 30% of the band of the three measurements up to each peak (20 points, worst +15%; physim/Usman 0.6–1.6). Reported, the documented limit: the tails above each peak (pre-equilibrium emission, which physim and PACE4 leave out), σ(3n)/σ(4n) at 49.6 MeV 0.055 against 0.16–0.27 measured; 1n (peak at the barrier) 2.4× high at its peak, none above 31 MeV. The band is wide: a common factor 0.39–1.13 on physim's cross sections would still pass. physim's stopping through the inferred stack gives every foil's energy within its spread (worst −0.8 MeV) |
| J. M. Allmond et al., Phys. Rev. C 90, 034309 (2014): ⁶²Ni at 1.8 MeV/u (111.6 MeV) on ~1 mg/cm² natural carbon at the HRIBF, carbon recoils in two BareBall rings (14–28°, 28–44°), nine clovers; <0⁺‖E2‖2⁺> 0.3038(37) and 0.3002(24) e b per ring, B(E2↑) 0.0906(37) e²b² | The `coulex_ni62_c_allmond` example (thickness tuned to the measured ⁵⁸Ni/⁶⁰Ni losses); the check `reproduction_allmond2014` (`physim.nuclear.coulex_reproduction`, `scripts/reproduce_allmond2014.py`, references and bands in `coulex/allmond2014.toml`) | Stopping, Coulomb excitation through a thick target to all orders, the rates, the events, stored runs and the analysis | ✅ all orders judged, first order in brackets. σ_C/σ_R per ring +0.79/+0.69% [+1.49/+1.09%] against GOSIA as run (integer masses, DIP 3.8; the E1 polarisation in physim too), band 1% as the verification's plan fixed it beforehand (this register said 2% for a time, written after about 1.3% had been found). Until physim had the E1 polarisation (2026-11-25) this was a known miss, +1.32/+1.16%; GOSIA's integer masses (0.6%) and its non-relativistic kinematics (0.24%) are what is left. GOSIA with nuclear masses and DIP 3.8, not judged: +0.20/+0.12%. ⁵⁸Ni, ⁶⁰Ni, ⁶²Ni losses 42.40 MeV each (42.7(8), 42.1(8), assumed 42.4(8)); +0.20/+0.12% [+0.90/+0.53%] against GOSIA with nuclear masses and DIP off, physim without it too (0.5%, not a blind band either: 0.2% was known); closure, the mean of four seeds (1001–1004) −0.36/−0.16/−0.19% for ring 2/ring 3/both (both −0.28% until 2026-11-25, with an angular-correlation factor averaged over the particles seen instead of the excitations), MC errors 0.73/0.44/0.38%, scatter errors 0.84/0.19/0.36% (a quarter of the paper's statistical error + 3 MC errors: 2.8/1.7/1.5%; 160 other seeds, off line, before the E1 polarisation and that: +0.15(11)/−0.06(7)/−0.09(6)%, and without the half-life 60: +0.04(17)/+0.33(10)/+0.14(8)%); that sample as a stored run, excitations per particle +0.04/+0.04% from the rates (0.1% + 3 errors) and B(E2) −0.2/+0.6/+0.4% (errors 1.8/1.1/0.9%); the paper's per-ring matrix element through GOSIA's ratio(M) and back through physim's: B(E2) −0.8/−0.7% (4.1%; pulls −0.3/−0.4 against the statistical errors; [−1.5/−1.1%]). The paper gives no yields, so this tests that physim extracts the same B(E2) from the same data as GOSIA did, not B(E2) itself. Q_s = ±0.1 e b moves the ratio by ±2.7/±1.9% (GOSIA ±2.7/±1.9%; 0.1% of the change), ENSDF's higher states feed +1.17/+0.77% (GOSIA +1.17/+0.77%; physim +0.12/+0.27% off it, +0.21/+0.18% with twice the check's angles and beam energies; both schemes without internal conversion, as GOSIA's input: with ENSDF's, α = 1.72·10⁻⁴ of 2⁺ → 0⁺, in the scheme with the higher states alone, until 2026-11-25, it was −1.4/−1.98%) (0.5 percentage points and 2% of GOSIA's change; with one orbit for all pairs of levels, until 2026-11-24, the changes were 6% and 19–20% too large and passed the percentage points). Worst of its band 0.79 (GOSIA as run in ring 2) |

Fourth slice: shell-model level schemes. physim does not solve the shell model itself: it drives KSHELL (the
2019 release of Comput. Phys. Commun. 244, 372, GPLv3, bundled as an outside program), writes its inputs and reads
its logs. Limits: spaces a laptop can diagonalise (a run may use at most half the free RAM, unless overridden), no
truncation yet. The benchmarks are in `tests/reference/nuclear/shellmodel/benchmarks.csv` with KSHELL's logs as
fixtures. The comparison with ENSDF is information, not a test: the shell model is a model, and its rms deviation
from experiment is reported as the papers report it.

| Capability | Engine | App | vs tools | vs literature | Page | Backlog |
|---|---|---|---|---|---|---|
| Shell-model calculations (KSHELL) | ✅ | ⬜ | ✅ KSHELL as physim runs and reads it against published calculations with the same Hamiltonian: ²⁸Si/USDB vs ANTOINE (Frycz et al. 2024: E(0⁺₁), B(E2; 2⁺→0⁺), B(E2; 4⁺→2⁺), Q_s(2⁺₁) with e_p 1.31, e_n 0.46) and ⁴⁸Cr/KB3 vs the exact energies of five states (Pittel et al. 2009), within one unit of the digit each source prints (0.1 MeV, 1 keV, 0.1 e²fm⁴, ...; worst 99%: B(E2; 2⁺→0⁺) 78.20 vs 78.3); M-scheme dimensions exact | ✅ physim's dimension count vs 8 published dimensions of ²⁸Si, ⁴⁸Cr, ⁵⁶Ni (exact). Reported, not judged: ²⁸Si/USDB vs ENSDF, excitation energies rms 320 keV over 7 states matched by J^π and order (largest 0⁺₃ +636 keV, at the prolate band head that Frycz et al. find USDB does not reproduce, and 2⁺₃ +485 keV), B(E2) 1.2–3.1× ENSDF's | [Shell model](../shell-model.md) | — |

Fifth slice: lifetime measurements, the plunger (RDDS) and the Doppler-shift attenuation method (DSAM), backlog item
73, built in waves. The checks are in `physim.nuclear.validation` with their references in
`tests/reference/nuclear/lifetimes/` (closed forms, TRIM, APCAD's and AhKin's examples, NaPyTau, the published
experiments); the page is [Lifetimes](lifetimes.md). Since wave 3C the DDCM is checked against NaPyTau on
physim's own simulated plunger spectra, τ comes back from simulated runs of known τ (closure), and the three
published experiments (⁹⁸Ru, ⁵⁰Cr, ¹⁴³Sm) are simulated from examples made after the papers and analysed as the
papers did (`physim.nuclear.lifetime_reproduction`). Three DSAM checks fail their tolerances and are listed as known
failures in `test_nuclear_validation.py`: the two line-shape comparisons and the ¹⁴³Sm reproduction (a stopping
difference, below).

| Capability | Engine | App | vs tools | vs literature | Page | Backlog |
|---|---|---|---|---|---|---|
| Lifetimes: recoil distance (RDDS) | ✅ physics, recoil transport, Coulex and fusion plunger runs, decay-curve fit and DDCM (singles, gated, three-foil) | ⬜ | ✅ TRIM (SRIM 2013), ⁹⁸Ru at β = 4.5% into the Cu stopper with SRIM's stopping: range 3%, straggling 15%, lateral spread 6%, energy 2% and angles 6–10% after a quarter and half of the range pass (worst 72% of a tolerance); NaPyTau on physim's simulated ⁹⁸Ru plunger areas, τ(x) at 31 distance×degree points within 1% (worst 1.3e-6) | ✅ closed forms: one level with and without time dilation and Bateman's two-level cascade (1e-9; worst 4e-10); ⁹⁸Ru deorientation G_k and G̃_k from Radeck et al.'s Table I (1e-8); closure within 2σ (one level 8.48 ± 0.17 and 8.42 ± 0.11 ps for 8.36, fed from above 12.26 ± 0.40 for 12.20); reproduced ⁹⁸Ru 2⁺₁ 8.35 ± 0.27 ps against 8.36(29) with the deorientation correction (2σ combined) and ⁵⁰Cr 2⁺₁ 12.27 ± 0.38, 4⁺₁ 2.25 ± 0.17 ps inside the published values' 2σ spread; the planner's recommended distances inside the papers' ranges to within a factor 1.5 at each end (worst 1.24); the recoils' stopping against measured data is `stopping_heavy_vs_data` | [Lifetimes](lifetimes.md) | 73 |
| Lifetimes: Doppler-shift attenuation (DSAM) | ✅ physics, recoil transport, DSAM kernel and line shapes, line-shape fit | ⬜ | 🟡 TRIM, ¹⁴³Sm at β = 1.6% into the Au backing, same tolerances, pass (worst 84%). Line shapes fail: APCAD's example (τ = 1 ps, 29 angles), widths agree (mean ratio 1.00) but the first-order shift is 17% larger (10% allowed; physim's ¹²⁰Sn leave the Cu backing where Geant4's stopping keeps them); AhKin's ⁸³Y line, shift inside AhKin's 1σ band but 25% wider (10% allowed; physim's scattering turns the residues, AhKin follows the velocity along the beam only) | 🟡 closed form F = α/(α + τ) for stopping ∝ speed and the distribution of v/v₀ at emission (1e-6; bin shares 2e-5, worst 97%); closure on two ¹⁴³Sm lines within 2σ (−0.8σ, −0.0σ). The ¹⁴³Sm reproduction fails: fitted with SRIM's stopping as the paper did, spectra simulated with physim's stopping give lifetimes 3–27% above the paper's (physim slows ¹⁴³Sm in Au faster, range 19% shorter); seven levels within the published uncertainty plus ±5% (2σ), the two upper limits above their limits; the planner's backing advice puts the paper's Au within 80% of the best backing's sensitivity for every ¹⁴³Sm level (worst 0.19) | [Lifetimes](lifetimes.md) | 73 |

Next slices: inelastic scattering through the nuclear interaction, then transfer reactions.

## Foundations: the classical and quantum engine

The general-purpose engine underneath the planner. Its validation is against analytic solutions and conservation
laws rather than nuclear-physics tools; the results are listed in the "Done" note at the top of each backlog item,
and the tests are in `tests/` (Rust) and `tests/python/`. None of it is in the planner app yet.

| Capability | Engine | Validated against | Backlog |
|---|---|---|---|
| Integrators (symplectic, high-order, adaptive) | ✅ | Convergence orders, energy conservation | 08, 17 |
| Gravity (direct and Barnes–Hut), springs, built-in forces | ✅ | Tree error scaling vs direct sum; resonance curve, 1PN periapsis advance | 10, 19 |
| Charged particles: Lorentz force, Boris pusher, Coulomb | ✅ | Analytic gyration, E×B drift, Rutherford deflection angle | 09 |
| Constraints, collisions, rigid bodies | ✅ | Pendulum periods, textbook collision velocities, Euler-top flip period | 07, 11, 13 |
| Molecular dynamics | ✅ | Thermostat statistics, conservation | 12 |
| Chaos indicators | ✅ | Kepler λ → 0, symplectic ± pairing of Lyapunov spectra | 18 |
| Fields on grids: waves, diffusion, Poisson | ✅ | Analytic solutions, convergence orders | 20 |
| Single-particle quantum mechanics | ✅ | Analytic levels, transmission, norm conservation | 32 |
| Seeded random numbers (thermostats, initial conditions) | ✅ | χ² uniformity, correlations; bit-identical across thread counts and restarts | 31 |

## Adding an entry

When a capability reaches **Engine**, give it a page in this folder (`kinematics.md`, `stopping-power.md`, …) with
these sections, and link it from the table:

- **Model**: the equations used, and any data (with version).
- **Assumptions and range of validity**: where the model holds and where it does not.
- **Validation**: for each comparison, the reference (tool and version, or paper), the cases compared, the
  tolerance, a plot or a link to the validation notebook, the test that enforces it, and the date and commit.
- **Known limitations**: what is missing or approximate.
- **References**: papers and data sources to cite.
