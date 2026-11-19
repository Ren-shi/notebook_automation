# Lifetime measurements: recoil distance (RDDS) and Doppler-shift attenuation (DSAM)

Backlog item 73 · modules `physim.nuclear.lifetimes`, `physim.nuclear.transport` · tests
`tests/python/test_nuclear_lifetimes_physics.py`, `test_nuclear_transport.py`, `test_nuclear_validation.py` · theory
{doc}`../theory/lifetimes` · validation notebook {doc}`../validation/lifetimes`

Built in waves (November 2026). What exists: the lifetime physics (B ↔ τ, Bateman decay curves with feeding, the
plunger's unshifted share against distance, the DSAM attenuation and line shape from a slowing-down kernel,
deorientation in vacuum) and the recoil transport (slowing through a stack with vacuum gaps and a stopper, time
dilation, the ZBL scattering Monte Carlo, the DSAM kernel). Simulated plunger and DSAM spectra, the lifetime analysis
(decay curve, DDCM, line-shape fit) and the planner come in the next waves; their checks are listed below as
waiting.

## Model

- **Plunger (RDDS).** A level of mean life τ formed at the target decays after the flight time t = d/(βc cos θ)
  to the stopper with probability e^(−t/(γτ)) (the unshifted share). With feeding, Bateman's equations for the
  cascade (a branching matrix, row = from), solved as a matrix exponential that stays stable when two lifetimes are
  nearly equal.
- **DSAM.** The kernel is the recoils' velocities against time, independent of τ (AhKin's idea): the line shape of a
  level is the kernel re-weighted by its decay-time distribution (with feeding and side feeding), each γ ray at
  E₀ √(1 − β²)/(1 − β cos ψ), ψ between the recoil's velocity and the detector. The attenuation F(τ) is the mean
  β·cos at emission over its value at formation.
- **Transport.** Slowing from physim's stopping powers ({doc}`stopping`); vacuum gaps at constant speed; the ZBL
  universal potential's classical scattering for the recoil's nuclear collisions above a cut, the rest as continuous
  nuclear stopping and a Gaussian turn (as LINESHAPE and TRIM do in spirit).
- **Deorientation in vacuum.** The static model Radeck et al. use: G_k(t) = α_k + (1 − α_k) e^(−Γ_k t); the shifted
  component sees G_k averaged over the decays in flight (their Eq. 7), the unshifted one G_k at the stopper.

## Assumptions and range of validity

- Recoil speeds up to a few per cent of c (time dilation included; DSAM leaves it out, below 0.2% at β < 6%).
- The stopping powers set the DSAM time scale: physim's are checked against measured slow heavy ions in elemental
  solids (`stopping_heavy_vs_data`, 3620 points: medians within 25% below 0.1 MeV/u, 12% at 0.1–1 MeV/u, 8% above),
  and they differ from SRIM's for these recoils (⁹⁸Ru in Cu: ranges 7% shorter; ¹⁴³Sm in Au: 19% shorter) and from
  Geant4's and AhKin's tables (below). A ±10% scaling of the stopping is the uncertainty to carry into a fitted τ.
- Deorientation: the static model only, with α_k and Γ_k given (from a measurement such as Radeck's), not computed.

## Validation

The checks are in `physim.nuclear.validation`, with their references in `tests/reference/nuclear/lifetimes/` (its
README gives each file's provenance). Run them with
`python -c "from physim.nuclear import validation as v; print(v.report(v.run(capability=v.RDDS) + v.run(capability=v.DSAM)))"`.
The recoils' stopping against measured data is the stopping page's check `stopping_heavy_vs_data`, not repeated
here.

### Recoil distance (RDDS)

| Check | Reference | Tolerance | Result |
|---|---|---|---|
| `rdds_one_level_analytic` | Closed form exp(−d/(γβcτ)), with and without time dilation; τ = 1, 10, 100 ps, β = 2–10%, d = 1 µm–5 mm (`analytic/rdds_one_level.csv`) | 1e-9 relative | pass (worst 3e-10) |
| `rdds_two_levels_bateman` | Bateman's equations, a level fed from one above, three lifetime pairs (one with τ₁ ≈ τ₂) and two population splits (`analytic/rdds_two_levels.csv`) | 1e-9 relative | pass (worst 4e-10) |
| `deorientation_ru98` | ⁹⁸Ru 2⁺₁, Radeck et al. 2012 Table I (α₂ 0.18, Γ₂ 0.142/ps, α₄ 0.24, Γ₄ 0.314/ps; v = 13.49 µm/ps, τ = 8.36 ps): G_k(d) and the shifted G̃_k(d), Eq. 7 by Gauss–Legendre quadrature, 0–450 µm | 1e-8 | pass (worst 1e-15); at 450 µm G₂ = 0.187 unshifted, G̃₂ = 0.562 shifted |
| `rdds_transport_vs_trim` | TRIM (SRIM 2013), 5000 ions per run: ⁹⁸Ru at β = 4.5% into Cu (the stopper), physim slowing with SRIM's stopping so that only the scattering differs | range 3%, range straggling 15%, lateral spread 6%, transmitted 0.02, mean energy 2%, rms and 90% angle 6%, median angle 10% | pass (worst 72% of a tolerance: the median angle after a quarter of the range +7.2%) |
| `rdds_ddcm_vs_napytau` | NaPyTau's DDCM on the same simulated spectra (`napytau/`) | τ within 1% | waits for wave 3C |
| `rdds_closure` | τ back from simulated plunger runs of known τ, decay curve and DDCM | 2σ of the fit | waits for wave 3C |
| `rdds_ru98_reproduction` | ⁹⁸Ru after Coulex (Radeck et al. 2012): 2⁺₁ 8.36(29) ps with deorientation, 4⁺₁ 2.31(16), 2⁺₂ 1.7(2) ps | 2σ | waits for wave 3C |
| `rdds_cr50_reproduction` | ⁵⁰Cr after fusion, three-foil plunger (Beckers et al. 2022): 2⁺₁ 12.2(4), 4⁺₁ 2.4(3) ps, against Coulex, RDDS and DSAM values | 2σ spread of the published values | waits for wave 3C |

### Doppler-shift attenuation (DSAM)

| Check | Reference | Tolerance | Result |
|---|---|---|---|
| `dsam_exponential_analytic` | Stopping ∝ speed, v = v₀ e^(−t/α): F = α/(α + τ) exactly, and the distribution of v/v₀ at emission as the 0° line shape in 20 bins; α = 0.2–1 ps, τ = 0.01–5 ps (`analytic/dsam_exponential.csv`) | F 1e-6 relative, bin shares 2e-5 | pass (F 8e-8; shares 1.9e-5, 97% of the tolerance: the line shape spreads each recoil's energy linearly between the kernel's times) |
| `dsam_transport_vs_trim` | TRIM (SRIM 2013): ¹⁴³Sm at β = 1.6% into Au (the backing), SRIM's stopping | as for ⁹⁸Ru | pass (worst 84%: the range +2.5%) |
| `dsam_lineshape_vs_apcad` | APCAD's OneShape example (StopSim with Geant4): 1000 keV, τ = 1 ps, the excited ¹²⁰Sn after 400 MeV ¹²⁰Sn on 1.5 mg/cm² C on 10 mg/cm² Cu, 29 angles 20°–160°. physim: its kinematics with APCAD's cross-section table, its stopping and scattering, `lifetimes.lineshape`, 3 keV FWHM; compared with the spectra less their stated background, inside each angle's kinematic window | first-order centroid shift (half the difference between θ and 180° − θ, 20°–85°) and the mean rms width within 10%; χ²/bin reported | **fails**: widths agree (mean ratio 1.00), shift 17% larger (+33.0 against +28.9 keV at 20°, −31.6 against −26.5 at 160°); χ²/bin 4.7 (1.0 off the line) |
| `dsam_lineshape_vs_ahkin` | AhKin's example: ⁸³Y from ³²S + ⁵⁸Ni at 135 MeV slowing in Ta, the 1193 keV line at 52.8°, τ = 0.1313 ps. physim with AhKin's inputs (its stopping table, ten-level cascade with side feeding, detector cosines, 10% energy spread, resolution, background) | shift inside AhKin's 1σ band (its lines at 0.0986 and 0.1652 ps: +7.72 to +8.80 keV); rms width within 10%; χ²/bin against the measured line no larger than AhKin's at the band's edges (6.6) | **fails**: shift +8.12 keV against +8.25 (inside), width 6.30 against 5.06 keV (+25%), χ²/bin 11.0 (AhKin 5.6) |
| `dsam_closure` | τ back from simulated DSAM spectra of known τ, with side feeding | 2σ of the fit | waits for wave 3C |
| `dsam_sm143_reproduction` | ¹⁴³Sm after fusion on a Au backing (Rajbanshi et al. 2018): seven lifetimes 0.41–1.02 ps, two upper limits | published uncertainty plus SRIM's ±5% stopping | waits for wave 3C |

**Why the two line-shape comparisons fail.** Both are judged on physim's full chain, and both disagreements are
understood; neither tolerance was moved.

- *APCAD.* In physim the range of the ¹²⁰Sn (about 300 MeV) in Cu is 11.4 mg/cm², so 39% of them leave the
  10 mg/cm² backing and fly on at their exit speed, which makes the line's fully shifted edge stronger and its
  stopped peak weaker. APCAD's spectra have more of the stopped peak and the flying component lower in energy:
  StopSim slows the ions more. Its stopping is Geant4's `G4IonParametrisedLossModel`, whose ICRU 73 tables stop at
  Z = 26 and are scaled for Sn; physim's is checked against measured heavy-ion data. With physim's electronic
  stopping raised by 15% every ion stops and the shifts come within 4–12% of APCAD's (a test run, not a check). The
  widths, which depend on the kinematics and scattering more than on the stopping, agree.
- *AhKin.* With AhKin's own stopping table the centroid lands inside AhKin's 1σ band, but physim's line is 25% wider.
  physim's Monte Carlo turns the residues (its angular spread agrees with TRIM's within 6%), and their sideways
  motion broadens the line at 52.8°. AhKin follows only the velocity along the beam, reduced by Blaugrund's mean
  cosine; with straight paths physim's line is 9% narrower than AhKin's and its shift 12% larger, so AhKin lies
  between physim's straight and scattered lines. The measured line prefers AhKin's width at AhKin's τ; whether
  physim's line fits it at another τ is the line-shape fit's question (wave 2C, closure in 3C). With physim's own
  stopping (5–17% above AhKin's table at 1–40 MeV) the shift is +7.36 keV, just below the band.

`test_dsam_line_shape_disagreements_are_the_documented_ones` pins both disagreements, so a change in the transport,
the stopping or the line shape shows up. Broken formulas (the unshifted share, the deorientation in flight, the
attenuation, the recoil's scattering angles) make their checks fail (`test_a_broken_formula_fails_its_check`).

## Known limitations

- No simulated plunger or DSAM spectra yet, no lifetime fit: the DDCM, closure and published reproductions wait
  for waves 2 and 3.
- The DSAM line shape treats the detector as a point at (θ, φ = 0); a finite detector is an average over angles (as
  the AhKin comparison does with five).
- Deorientation is the static model with given parameters; dynamic deorientation and α_k, Γ_k from atomic physics
  are not modelled.
- The stopping powers are the largest systematic uncertainty of a DSAM lifetime (above).

## References

- D. Radeck et al., Phys. Rev. C 85, 014301 (2012): ⁹⁸Ru, RDDS after Coulex, deorientation (Table I).
- C. Beckers et al., Nucl. Instrum. Methods A 1042, 167418 (2022): ⁵⁰Cr, three-foil plunger.
- S. Rajbanshi et al., Phys. Lett. B 782, 143 (2018): ¹⁴³Sm, DSAM.
- A. Garzón, W. Rodríguez, F. Cristancho, M. Tao, Comput. Phys. Commun. (2020): AhKin (MIT licence).
- C. Stahl, J. Leske, M. Lettmann, N. Pietralla, Comput. Phys. Commun. (2017): APCAD and StopSim (CC BY 4.0).
- NaPyTau, github.com/h-mayr/NaPyTau_Public (EUPL-1.2); J. F. Ziegler, J. P. Biersack, U. Littmark, *The Stopping
  and Range of Ions in Solids* (1985), the ZBL potential; SRIM 2013.
