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
- **First order stays available:** `[reaction] order = "first"` gives the events of first-order theory (bit for
  bit those of earlier versions). A setup without a level scheme is first order in any case, since its one matrix
  element allows nothing else.

The solution is made once per setup, at the beam energy in the middle of the target, on a grid of 90 CM angles
(2° apart) with a relative accuracy of 10⁻⁴ in the probabilities, and kept. The equations of all angles are
integrated together (`CoupledChannels.solve_many`), and the forward angles at which first order excites less than
10⁻⁶ of the collisions are left in first order (below 14° in the ⁶²Ni example, where no population changes by
more than 3 × 10⁻⁹). The table then takes about 2 s for the six levels of ⁶²Ni, where one angle after the other
takes 16 s; first order takes 0.2 s.

For `coulex_ni62_pt194` (⁶²Ni at 237 MeV on ¹⁹⁴Pt, the beam excited), at θ = 150° in the CM frame:

| Level of ⁶²Ni | First order | All orders |
|---|---|---|
| 2⁺, 1172.9 keV | 0.1932 | 0.1845 |
| 0⁺, 2048.6 keV | 0 | 0.0081 |
| 2⁺, 2301.9 keV | 0.0006 | 0.0017 |
| 4⁺, 2336.5 keV | 0 | 0.0058 |
| any level | 0.1938 | 0.2001 |

The 2⁺ state loses 4.5% to the levels above it and gets it back, with more, when they decay through it: 3% more
2⁺ → 0⁺ γ rays than first order gives, 4⁺ → 2⁺ γ rays (1163.5 keV) at 3% of them and 0⁺ → 2⁺ γ rays (875.6 keV)
at 4%.

**What this leaves out.**

- Every excitation has the kinematics of the `[reaction]` state (Q = −E* of that state), whichever level it
  reached: with ⁶²Ni left in its 4⁺ state the two nuclei share 1.2 MeV less kinetic energy than the events give
  them.
- Elastic scattering is not reduced by the excitation probability (20% at backward angles in the ⁶²Ni example).
- The analysis (`physim.nuclear.analysis`) still extracts B(E2) from the 2⁺ → 0⁺ peak with first-order theory. On
  events with all orders the peak holds the feeding from above, and a line the Doppler-corrected peak does not
  resolve from it (the 4⁺ → 2⁺ line of ⁶²Ni lies 9.4 keV below a peak 15 keV wide): for the ⁶²Ni example it
  returns a B(E2) 5% above what it returns on first-order events (3% from the feeding, the rest from the
  unresolved line). With the ENSDF scheme of ⁵⁸Ni added to `coulex_ni58` (¹⁶O at 30 MeV), all orders give 5% fewer
  2⁺ → 0⁺ γ rays in the CD and 2% fewer in the strip detectors than first order, mostly through the reorientation.
- The coupled equations have no E1 polarisation (virtual excitation of the giant dipole resonance), no excitation
  of both nuclei in one collision, and no nuclear interference above the safe energy.
- The plan's coincidence rate takes the angular correlation of the `[reaction]` state in first order; the events
  use the solution's. In a 40 min run of the ⁶²Ni example the two agree within the run's statistics (1% in the
  CD).

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
| The 4⁺ → 2⁺ and 0⁺ → 2⁺ lines in a 40 min run of ⁶²Ni on ¹⁹⁴Pt | Multi-step yields (`Multistep.yields`) times what the crystals record in the line's window | 3σ and 3% (absolute); 3σ and 1% (relative to the 2⁺ line) | `test_a_run_of_62ni_on_194pt_shows_the_4plus_and_0plus_lines_of_the_multistep_solution` |
| Excitation rate in each particle detector, events against `Rates` | The same excitation probability in both | 3σ and 0.5% | `test_the_rates_and_the_generator_share_the_excitation_probability` |
| Monte Carlo against analytic rates, every channel and detector of `coulex_ni58` | `Rates` | 4σ | `test_monte_carlo_rates_agree_with_analytic_rates[coulex_ni58]` |
| Inelastic peak below the elastic one by the kinematic difference | `TwoBody` with and without E* | 10% | `test_coulomb_excitation_in_the_report` |
| Against GOSIA and a published measurement | Requested from the user (`tests/reference/nuclear/pending/coulex_gosia.md`) | — | `validation.coulex_vs_gosia` (pending) |
