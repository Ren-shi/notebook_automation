# Coulomb excitation

`physim.nuclear.coulex` computes the probability that the electric field of one nucleus excites the other as they
pass. `physim.nuclear.gamma` computes the Doppler-shifted energy of the γ ray emitted afterwards. Below the Coulomb
barrier the nuclei never touch, so the excitation is purely electromagnetic and calculable. This makes Coulomb
excitation a standard way of measuring B(E2) values and other collective properties.

## First-order semiclassical theory

The relative motion is a classical Rutherford orbit (see {doc}`rutherford`), valid when the Sommerfeld parameter
η = Z₁Z₂e²/ħv ≫ 1. Along that orbit, the time-dependent multipole field of one nucleus (charge Z) excites the other
from a 0⁺ ground state to a state of energy E* and spin λ. To first order in perturbation theory (Alder and Winther)
the excitation amplitude is

$$
a_\mu = \frac{4\pi Z e}{i\hbar(2\lambda+1)} \langle \lambda\mu | \mathcal{M}(E\lambda,\mu) | 00 \rangle^*
\int_{-\infty}^{\infty} e^{i\omega t}\, \frac{Y_{\lambda\mu}(\hat r(t))}{r(t)^{\lambda+1}}\, dt ,
\qquad \omega = E^*/\hbar ,
$$

and since $|\langle \lambda\mu | \mathcal{M} | 00\rangle|^2 = B(E\lambda\uparrow)/(2\lambda+1)$ for each μ, the
probability is

$$
P(\theta) = \left(\frac{4\pi Z e^2}{\hbar c\,(2\lambda+1)}\right)^2 \frac{B(E\lambda\uparrow)}{2\lambda+1}
\sum_\mu \left|\frac{Y_{\lambda\mu}(\tfrac{\pi}{2}, 0)\, I_\mu(\theta, \xi)}{\beta\, a^\lambda}\right|^2 .
$$

**The orbit.** It is parametrised in closed form: r = a(ε cosh w + 1), t = (a/v)(ε sinh w + w), with
ε = 1/sin(θ/2). The position is x = a(cosh w + ε), y = a√(ε² − 1) sinh w from the scattering centre, which gives
the azimuth φ(w) in the orbit plane. The orbit integral is then

$$
I_\mu(\theta, \xi) = \int_{-\infty}^{\infty} \frac{e^{i\xi(\varepsilon\sinh w + w)}\, e^{i\mu\varphi(w)}}
{(\varepsilon\cosh w + 1)^\lambda}\, dw .
$$

**Symmetrisation.** The state takes energy from the motion, so the orbit is symmetrised between the initial and
final velocities: a = Z₁Z₂e²/(μ v_i v_f), β = √(v_i v_f)/c. The adiabaticity is

$$
\xi = \frac{Z_1Z_2e^2}{\hbar}\left(\frac{1}{v_f} - \frac{1}{v_i}\right).
$$

When ξ ≳ 1 the collision is slower than the nuclear motion, and the excitation dies away (the adiabatic cutoff).

**The cross section.** It is dσ/dΩ = P(θ) dσ_R/dΩ, with the symmetrised Rutherford cross section
(a²/4)/sin⁴(θ/2), in the CM frame. The lab cross section and the total follow as for elastic scattering.

**Numerics.**
- The integrand at −w is the conjugate of that at w, so I_μ = 2 Re ∫₀^∞. That half is done in Gauss–Legendre
  panels (16 points each) spanning about four radians of the phase, out to a cut; beyond the cut the oscillating
  tail is taken by its endpoint expansion (integration by parts, twice), and the cut is placed where the
  expansion's remainder is below 10⁻⁹. P(θ) is tabulated every 1° (181 orbit integrals, about 50 ms) and
  interpolated with four-point cubics, to ~10⁻⁵.
- `probability(θ, exact=True)` integrates each angle directly.

## When it applies

- **Cline's safe distance.** The closest approach a(1 + 1/sin(θ/2)) must be at least 1.25 (A₁^⅓ + A₂^⅓) + 5 fm.
  Closer in, nuclear forces interfere. `Coulex.max_safe_angle()` gives the largest safe CM angle, and the warnings
  name it.
- **η ≫ 1** for the semiclassical orbit; a warning is given below 5.
- **First order needs P ≪ 1.** At P ≳ 0.1, multi-step excitation and reorientation matter. GOSIA is the standard
  code for that, and a warning says so. With a level scheme of the excited nucleus the simulated events go beyond
  first order (next section).
- **Not included in first order:**
  - nuclear–Coulomb interference;
  - excitation from states other than a 0⁺ ground state;
  - feeding from higher states.

## Beyond first order: the events of a setup with a level scheme

First order excites each level in one step from the ground state. A 4⁺ state has no E2 matrix element to a 0⁺
ground state, so first order never reaches it, and the same holds for an excited 0⁺ state. They are reached in two
steps through the 2⁺ state. `physim.nuclear.coupled.CoupledChannels` solves the coupled equations for the
amplitudes of every magnetic substate of every level along the orbit, which holds all orders: the two-step
excitations, the depletion of the ground state and of the 2⁺ state, and the reorientation through the diagonal
matrix elements (the quadrupole moments).

When the setup has a level scheme of the excited nucleus (`[levels]`) whose decays follow more than one
transition, the simulated events take that solution (`gamma.cascade_excitation`):

- **Which level an excitation starts in** is drawn from the populations the collision leaves at the event's CM
  angle, P_n(θ), before any decay. The decay then follows the scheme's branches as before.
- **The excitation probability** of the channel, in the event generator's table and in the analytic rates alike
  (`rates.excitation_probability`), is the sum of those populations, Σ P_n(θ). `Rates` and the events therefore
  agree; "excitations" counts every collision that leaves the nucleus excited, in whichever level.
- **The angular distribution of each γ ray** comes from the same solution: the amplitudes give the density
  matrix of every level, in two steps or one, its statistical tensors follow, and the decay hands them down the
  cascade. The 4⁺ → 2⁺ γ rays of a 4⁺ state reached in two steps are emitted with that state's own alignment;
  nothing is taken as isotropic.
- **The γ rays of one cascade are correlated** (`physim.nuclear.correlation`). Each excitation carries the
  density matrix ρ of its level down the cascade. A γ ray from spin I to I_f by multipoles L, L′ (mixing ratio
  δ) is drawn along n with the probability W(n) ∝ Σ_σ Tr[A_σ ρ_n A_σ†], where ρ_n is ρ quantised along n, σ = ±1
  the photon's helicity and (A_σ)_{mM} = Σ_L c_L √(2L + 1) ⟨I_f m L σ | I M⟩ with c_L = 1, c_L′ = −σδ. The level
  it feeds is left in Σ_σ A_σ ρ_n A_σ† rotated back and normalised, and the next γ ray starts from that. Summed
  over directions this is the unobserved feeding above, so the singles are unchanged and only the γ–γ
  coincidences between crystals change: for an unoriented 4⁺ → 2⁺ → 0⁺ cascade the two γ rays follow
  1 + 0.1020 P₂ + 0.0091 P₄ of the angle between them. A converted transition hands its orientation on as an
  unobserved one, and `emission = "isotropic"` draws every γ ray on its own, without the correlation.
- **First order stays available:** `[reaction] order = "first"` gives the events of first-order theory (bit for
  bit those of earlier versions). A setup without a level scheme is first order in any case, since its one matrix
  element allows nothing else.

The solution is made once per setup, at the beam energy in the middle of the target, on a grid of 90 CM angles
(2° apart) with a relative accuracy of 10⁻⁴ in the probabilities, and kept. The equations of all angles are
integrated together (`CoupledChannels.solve_many`), and the forward angles at which first order excites less than
10⁻⁶ of the collisions are left in first order (below 14° in the ⁶²Ni example, where no population changes by
more than 3 × 10⁻⁹). The table then takes about 2 s for the six levels of ⁶²Ni, where one angle after the other
takes 16 s; first order takes 0.2 s.

For `coulex_ni62_pt194` (⁶²Ni at 237 MeV on ¹⁹⁴Pt, the beam excited), at θ = 150° in the CM frame, with the E1
polarisation correction of the default setup (d = 0.0038, below) and without it:

| Level of ⁶²Ni | First order | All orders | First order, no E1 polarisation | All orders, no E1 polarisation |
|---|---|---|---|---|
| 2⁺, 1172.9 keV | 0.1827 | 0.1756 | 0.1932 | 0.1844 |
| 0⁺, 2048.6 keV | 0 | 0.0074 | 0 | 0.0081 |
| 2⁺, 2301.9 keV | 0.0005 | 0.0016 | 0.0006 | 0.0017 |
| 4⁺, 2336.5 keV | 0 | 0.0053 | 0 | 0.0058 |
| any level | 0.1833 | 0.1899 | 0.1938 | 0.2001 |

The 2⁺ state loses 3.9% to the levels above it (4.5% without the correction) and gets it back, with more, when
they decay through it: 3% more 2⁺ → 0⁺ γ rays than first order gives, 4⁺ → 2⁺ γ rays (1163.5 keV) at 3% of them
and 0⁺ → 2⁺ γ rays (875.6 keV) at 4%.

**One particle per collision.** A collision either excites the nucleus or leaves it as it was, so elastic
scattering is the Rutherford cross section times 1 − P(θ): 20% fewer elastic particles at backward angles in the
⁶²Ni example, and elastic and excited particles together are Rutherford scattering.

**The E1 polarisation correction** (`[reaction] e1_polarisation`, on by default with d = 0.0038). During the
collision the field of the other nucleus also excites the giant dipole resonance, virtually, and this second-order
E1 process acts on the low-lying levels like an extra E2 interaction: proportional to the E2 matrix elements,
falling off as 1/r⁴ instead of 1/r³, and opposing the direct one. physim uses GOSIA's approximation of it
(D. Cline, T. Czosnyka et al., *GOSIA user manual* (2012), eq. 3.27–3.28, after K. Alder and A. Winther,
*Electromagnetic Excitation* (North-Holland, 1975), appendix J): every E2 collision function is multiplied along
the orbit by

$$
1 - z\,\frac{a}{r} = 1 - \frac{z}{\varepsilon\cosh w + 1},
\qquad
z = d\,\frac{E_p A_2}{Z_2^2\,(1 + A_1/A_2)} = d\,\frac{E_\mathrm{cm} A}{Z^2},
$$

with E_p the beam energy in MeV, 1 the projectile and 2 the target, and A, Z those of the excited nucleus
(E_cm from the mass numbers). The manual gives the formula for target excitation; for projectile excitation
GOSIA's code (subroutine `LOAD` of `gosia_20110524.14.f`) takes the same d E_cm A/Z² with A and Z of the
projectile. The parameter d is the setup's `e1_polarisation`. A setup that does not set it has d = 0.0038, what
GOSIA uses since version 20110524.11 (April 2023) when its input does not set it, so that a GOSIA run with default
settings compares like with like; GOSIA's manual has 0.005, its default before then. `e1_polarisation = 0` leaves
the correction out. (GOSIA's `DIP,` option is 1000 d, and there `DIP,0` counts as not set.) The classes themselves
(`CoupledChannels`, `Excitation`, `Coulex`) leave it out unless given a d. The factor is largest at the closest approach (a/r = 1/(ε + 1), a half at 180°). It applies to every
E2 coupling, the transitions and the diagonal elements alike, in the coupled equations and in first order
(GOSIA puts it into the collision functions every calculation of its starts from); E1 and E3 are unchanged. At
180° in the sudden limit (ξ → 0) the first-order amplitude is ∫ dw/(1 + cosh w)² = 2/3 and the correction takes
z ∫ dw/(1 + cosh w)³ = 4z/15 from it, so the probability falls by (1 − 0.4 z)², which the three solvers give.

What the default d = 0.0038 changes on the templates, against `e1_polarisation = 0` (populations and γ rays
integrated over all angles, the plan's rates in each particle detector, and B(E2) read from one simulated run made
with the correction, analysed with and without it):

| | z | Populations | γ rays | Plan's excitations | B(E2) read without the correction |
|---|---|---|---|---|---|
| `coulex_ni58`: ⁵⁸Ni, one state, first order | 0.0064 | 2⁺ −0.6% | −0.6% | −0.7% CD, −0.5% strips | −0.6% |
| `coulex_ni62_pt194`: ⁶²Ni (the beam) | 0.053 | 2⁺₁ −3.8%, 0⁺₂ −8.7%, 2⁺₂ −6.4%, 4⁺ −7.9% | 1172.9 keV −4.0%, 1163.5 keV −7.9% | −5.1% CD, −4.9% strips | −7.2% (clovers: −7.1%) |
| the same: ¹⁹⁴Pt (the target) | 0.022 | 2⁺₁ −0.6%, 2⁺₂ −1.9%, 4⁺₁ −2.1%, 6⁺ −4.0% | 328.5 keV −0.9% | — | — |

z is at the beam energy in the middle of the target (¹⁶O at 28.9 MeV, ⁶²Ni at 233.8 MeV). The levels reached in two
or three steps lose about twice what the first 2⁺ state loses, since each step is weakened. B(E2) moves more than
the yield because the excitation saturates (step 7 of the analysis), and analysed with the correction the same run
gives back the B(E2) it was made with. With the manual's d = 0.005 every change is about a third larger (the ⁶²Ni
2⁺₁ state −5.0%, B(E2) read without the correction −9.4%).

The reorientation of the ⁶²Ni 2⁺ state (⟨2‖E2‖2⟩ = +6.6 e fm², Q = +5 e fm², raising its population by 8.1% at
170° against Q = 0) keeps its sign and shrinks by only 1.2% of itself, since the diagonal coupling is weakened with
the transitional one. But the correction falls off faster than the excitation and so weighs the backward angles
more: the ratio of the 2⁺ populations at 160° and at 60° falls by 1.4%, which is what a diagonal element 28%
smaller (4.8 e fm², Q 1.4 e fm² lower) does without the correction (with d = 0.005: 1.9%, 37%, Q 1.8 e fm² lower).
That is why quadrupole moments are where it matters.

**What this leaves out.**

- Every excitation has the kinematics of the `[reaction]` state (Q = −E* of that state), whichever level it
  reached: with ⁶²Ni left in its 4⁺ state the two nuclei share 1.2 MeV less kinetic energy than the events give
  them.
- The other nucleus of the collision (¹⁹⁴Pt here) is excited as well when the setup has its level scheme, with
  the summed populations of its own coupled equations (96% at backward angles here, 0.85 γ rays of 328.5 keV per
  backward particle). The two excitations are treated as independent: a collision can excite both, and the
  excitation of one does not change the probability of the other. The particle keeps the energy of the
  collision without the other nucleus's excitation.
- First-order theory reads B(E2) too high from events with all orders: the 2⁺ → 0⁺ peak holds the feeding from
  above, and a line the Doppler-corrected peak does not resolve from it (the 4⁺ → 2⁺ line of ⁶²Ni lies 9.4 keV
  below a peak 15 keV wide, and is 3% of its window). The analysis therefore solves the coupled equations of the
  level scheme for the B(E2) that gives the measured yield, 4 to 9% below the first-order reading on the ⁶²Ni
  examples, with every other matrix element taken from the scheme. The excitation saturates: the yield grows as
  B(E2) to the power 0.67 there, so the uncertainties of B(E2) are 1.5 times those of the yield. With the ENSDF scheme of ⁵⁸Ni added to `coulex_ni58` (¹⁶O at 30 MeV), all orders give 5% fewer
  2⁺ → 0⁺ γ rays in the CD and 2% fewer in the strip detectors than first order, mostly through the reorientation.
- The coupled equations have no excitation of both nuclei in one collision and no nuclear interference above the
  safe energy. The E1 polarisation correction is GOSIA's approximation, with one parameter d for every nucleus,
  not a calculation from the giant dipole resonance of the nucleus at hand.
- The plan's coincidence rate takes the angular correlation from the same solution as the events: with a level
  scheme the 2⁺ state as the whole scheme leaves it, to all orders. The state alone in first order is more
  sharply aligned: on the ⁶²Ni example it would send 14% fewer γ rays to the crystal at 90° with the CD, and 3%
  more to those at 45° and 135°. Without a level scheme the plan and the events both take the state alone.

## γ rays and the Doppler shift

The excited nucleus decays in flight. A photon emitted at angle α to its velocity β has

$$
E_\gamma = E_0\, \frac{\sqrt{1-\beta^2}}{1 - \beta\cos\alpha} ,
$$

the Lorentz boost of the rest-frame photon (tested against the explicit four-vector transformation to 10⁻¹²).
`gamma.doppler_table(experiment)` pairs each particle detector, which fixes the direction of the excited nucleus
through the kinematics, with each `[[gamma_detectors]]` entry. It samples:

- the particle detector's face (10 × 10 points);
- the γ detector's opening cone (144 directions);
- 8 depths through the target, with exact kinematics at each and the emitter slowed to the target exit.

Each sample is weighted by the excitation cross section. The weighted mean is the Doppler-shifted energy, and the
weighted spread is the Doppler broadening, to which the γ detector's resolution is added in quadrature. The state is
assumed to decay after leaving the target, isotropically; lifetimes and γ angular distributions are not modelled.

## Rates, spectra, app and report

A Coulomb-excitation setup adds excitation channels next to the elastic ones (`rates.channels`):

- target excitation: the main nuclide of the target;
- projectile excitation: the beam on every nuclide of the target and backing.

Each excitation channel works like an elastic one:

- **Kinematics** use Q = −E*.
- **Cross section:** Rutherford's at each depth times P(θ). P is computed once, at the beam energy in the middle of
  the layer, because it changes slowly through a thin target. With a level scheme it is the probability of
  exciting any level (previous section).
- **Event generator:** the excitation channel is picked as often as the elastic one, the masses include E*, and
  each event's weight is multiplied by P(θ*) from a table every 0.25°. The rare inelastic events therefore get as
  many samples as the elastic ones.
- **Expected peaks:** Coulomb-excitation peaks are always listed, however weak.

The app has an **Excitation and γ rays** tab (P(θ), events per detector, the Doppler table). The report adds a
Coulomb-excitation section and `gamma.csv`.

For `coulex_ni58` (0.1 pnA of 30 MeV ¹⁶O on 0.5 mg/cm² ⁵⁸Ni): the CD records 1.28 excitations per second next to
660 elastic counts. The simulated rates agree with the analytic ones within 1.5σ for every channel and detector,
at 0.7% statistical error from 2 × 10⁶ events.

## Checks

| Check | Reference | Tolerance | Test |
|---|---|---|---|
| Orbit integrals at θ = 180°, ξ = 0, for E1, E2, E3 | ∫ dw/(1 + cosh w)^λ = 2, 2/3, 4/15 | 1e-8 | `test_backscattering_in_the_sudden_limit_is_a_closed_form` |
| P(180°) as ξ → 0 | Closed form above | 1e-6 | `test_excitation_probability_closed_form_at_180_degrees` |
| P at 60°, 120°, 170°, for ξ = 0 and 0.4 | The amplitude integrated in time along orbits from physim's engine (no hyperbolic parametrisation) | 1e-6 (achieved 5e-8) | `test_agrees_with_integration_along_engine_orbits` |
| Adiabatic cutoff | Monotonic fall with ξ; below 1% by ξ ≈ 6 | — | `test_adiabatic_suppression` |
| Scaling | P ∝ B exactly; projectile vs target excitation in the same orbit ∝ (Z₂/Z₁)² | 1e-12 | `test_scaling_with_b_and_with_the_exciting_charge` |
| Safe distance | Closest approach equals the safe distance at the largest safe angle | 1e-9 | `test_safe_distance` |
| Doppler formula | Four-vector Lorentz boost | 1e-12 | `test_doppler_formula_is_a_lorentz_boost` |
| Angles solved together | Each angle solved alone | 1e-12 | `test_many_angles_solved_together_are_the_single_solutions` |
| E1 polarisation strength z, target and projectile excitation | GOSIA manual eq. 3.28; GOSIA's subroutine `LOAD` | 1e-6 | `test_z_is_gosias_for_target_and_projectile_excitation` |
| E1 polarisation at 180°, ξ → 0 | (1 − 0.4 z)² from the integrals above | 1e-6 (first order), 1e-5 (coupled equations) | `test_its_size_in_sudden_backscattering_is_the_hand_worked_one` |
| E1 polarisation with a diagonal element | Second-order perturbation theory by nested quadrature | 2e-3 | `test_with_a_diagonal_element_the_solver_and_second_order_theory_agree_on_the_correction` |
| B(E2) from a run made with the E1 polarisation, analysed with it | The value put in | 8% (the run's statistics) | `test_the_analysis_gives_back_the_b_e2_of_a_run_made_with_the_correction` |
| The 4⁺ → 2⁺ and 0⁺ → 2⁺ lines in a 40 min run of ⁶²Ni on ¹⁹⁴Pt | Multi-step yields (`Multistep.yields`) times what the crystals record in the line's window | 3σ and 3% (absolute); 3σ and 1% (relative to the 2⁺ line) | `test_a_run_of_62ni_on_194pt_shows_the_4plus_and_0plus_lines_of_the_multistep_solution` |
| Excitation rate in each particle detector, events against `Rates` | The same excitation probability in both | 3σ and 0.5% | `test_the_rates_and_the_generator_share_the_excitation_probability` |
| Monte Carlo against analytic rates, every channel and detector of `coulex_ni58` | `Rates` | 4σ | `test_monte_carlo_rates_agree_with_analytic_rates[coulex_ni58]` |
| Inelastic peak below the elastic one by the kinematic difference | `TwoBody` with and without E* | 10% | `test_coulomb_excitation_in_the_report` |
| Against GOSIA and a published measurement | Requested from the user (`tests/reference/nuclear/pending/coulex_gosia.md`) | — | `validation.coulex_vs_gosia` (pending) |
