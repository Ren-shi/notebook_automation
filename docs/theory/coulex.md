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
  code for that, and a warning says so.
- **Not included:**
  - nuclear–Coulomb interference;
  - excitation from states other than a 0⁺ ground state;
  - feeding from higher states.

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
  the layer, because it changes slowly through a thin target.
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
| Monte Carlo against analytic rates, every channel and detector of `coulex_ni58` | `Rates` | 4σ | `test_monte_carlo_rates_agree_with_analytic_rates[coulex_ni58]` |
| Inelastic peak below the elastic one by the kinematic difference | `TwoBody` with and without E* | 10% | `test_coulomb_excitation_in_the_report` |
| Against GOSIA and a published measurement | Requested from the user (`tests/reference/nuclear/pending/coulex_gosia.md`) | — | `validation.coulex_vs_gosia` (pending) |
