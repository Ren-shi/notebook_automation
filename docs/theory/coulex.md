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

**The cross section.** In the CM frame (Alder, Bohr, Huus, Mottelson and Winther, Rev. Mod. Phys. 28, 432 (1956),
section II A; Alder and Winther, *Electromagnetic Excitation*, 1975, ch. II):

$$
\frac{d\sigma}{d\Omega} = P(\theta)\,\frac{v_f}{v_i}\,\frac{a^2}{4}\,\frac{1}{\sin^4(\theta/2)}
= P(\theta)\,\frac{v_i}{v_f}\,\frac{d\sigma_R^{\rm el}}{d\Omega}(E_{cm}),
$$

with the symmetrised a above. The second form is the same thing: the elastic Rutherford cross section at the
initial CM energy, (Z₁Z₂e²/4E_cm)²/sin⁴(θ/2), times v_i/v_f = √(E_cm/(E_cm − E*)). It is the convention of GOSIA
too (its CMLAB: 250 r₃ √(E/(E − ΔE)) d² ε⁴, d the head-on half-distance at the initial energy); integrated over
the ⁶²Ni-on-carbon rings and target below, the app's Rutherford factor and GOSIA's printed integrated Rutherford
agree to 0.03 %. The factor v_i/v_f is 1.04 for ⁶²Ni at 90 MeV on carbon (E_cm 14.7 MeV, E* 1.17 MeV), 1.03 in
`coulex_ni58`, and under 1 % in the other examples. Every place the app forms a Coulomb-excitation cross section
uses it: `Coulex.cross_section_cm` and `rutherford_cm`, the rates and the event generator (`rates.
rutherford_factor`), the Doppler table, and the solutions of a level scheme. (The other nucleus's γ rays are the
exception: they come per particle of the `[reaction]`, with no factor of their own.) Which level's v_f it is:

- `cross_sections` of `Excitation` and of `CoupledChannels` (mb, over all angles): **each level its own**. The
  coupled equations give a cross section to every excited level the collision populates, in one step or several.
- `Multistep.yields` (γ rays in the detectors): that of the **reference level**, the first that a matrix element
  joins to the ground state, for every level, with its kinematics.
- The rates and the events: that of the **`[reaction]` state**, whose kinematics the channel follows, for the
  populations of all levels.

The last two do as GOSIA does with its NCM level. A level at E_n then gets √((E_cm − E_n)/(E_cm − E_ref)) of what it
has with its own v_f: for the 2.3 MeV levels of ⁶²Ni 0.3 % less at 237 MeV on ¹⁹⁴Pt and 4 % less at 100 MeV on carbon.

**Which Rutherford a ratio divides by.** An experiment divides its excited counts by the elastic ones, so the
app's ratios σ_Coulex/σ_Rutherford divide by the **elastic** Rutherford cross section: `Coulex.ratio_to_elastic`
(P × v_i/v_f at one angle) and `DepthProbability.ratio` (the same integrated through a thick layer). GOSIA's
ratios divide by its own integrated Rutherford, which is the inelastic one above, with the inelastic kinematics
(for the ⁶²Ni rings, 11 % and 16 % above the elastic integral); P alone is the ratio to that one at one angle.

The lab cross section and the total follow as for elastic scattering.

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

## The E1 polarisation correction

During the collision the field of the other nucleus also excites the giant dipole resonance, virtually. This
second-order E1 process acts on the low-lying levels like an extra E2 interaction: proportional to the E2 matrix elements,
falling off as 1/r⁴ instead of 1/r³, and opposing the direct one. physim has GOSIA's approximation of it (D. Cline, T.
Czosnyka et al., *GOSIA user manual* (2012), eq. 3.27–3.28, after K. Alder and A. Winther, *Electromagnetic Excitation*
(North-Holland, 1975), appendix J; subroutines `LOAD`, `SNAKE` and `QE` of `gosia_20110524.14.f`). Every E2 coupling, of
the transitions and of the diagonal elements alike, in first order and in the coupled equations, is multiplied along the
orbit by

$$
1 - z\,\frac{a}{r} = 1 - \frac{z}{\varepsilon\cosh w + 1},
\qquad
z = d\,\frac{E_\mathrm{cm} A}{Z^2} ,
$$

with A and Z of the excited nucleus and E_cm in MeV. E1 and E3 are unchanged. In the first-order orbit integrals the factor
makes I_μ(2) into I_μ(2) − z I_μ(3).

- **The parameter d** is the setup's `[reaction] e1_polarisation`, 0.0038 when it sets none: what GOSIA takes when its
  input sets none (since version 20110524.11; Phys. Rev. C 91, 064602 (2015); GOSIA's manual and earlier versions have
  0.005). So a GOSIA run with default settings compares like with like. 0 leaves the correction out. GOSIA's input `DIP,`
  is 1000 d, and there `DIP,0` means the default. The classes themselves (`Coulex`, `Excitation`, `CoupledChannels`) leave
  it out unless given a d.
- **z.** GOSIA's two branches (target or projectile excited) both reduce to d E_cm A/Z² of the excited nucleus. E_cm is the
  one the collision has everywhere else in physim, from the nuclear masses and relativistic, and A is the mass number
  (the formula's A comes from the A^(5/3) of the dipole polarisability). GOSIA takes the masses that its input gives and
  non-relativistic kinematics; the difference is 0.1 to 0.3 % of z, that is 0.1 to 0.3 % of the correction.
- **One z for every pair of levels**, as in GOSIA, with each pair's own a (next section). Strictly it is z a that does not
  depend on the pair: the polarisation term has one size, and the orbit's a grows as the energy falls. A pair at the
  symmetrised energy √(E_n E_m) has z √(E_n E_m)/E_cm then, a few percent of z less in a light system. physim follows
  GOSIA, so that the two can be compared.
- **Its size.** The factor is smallest at the closest approach, 1 − z/(ε + 1), and so weakens the excitation more at
  backward angles than at forward ones, as a change of the quadrupole moment does. At 180° in the sudden limit (ξ → 0) the
  first-order amplitude is ∫ dw/(1 + cosh w)² = 2/3 and the correction takes z ∫ dw/(1 + cosh w)³ = 4z/15 from it, so the
  probability falls by (1 − 0.4 z)². A level reached in two steps loses about twice what the first 2⁺ level loses.

Against GOSIA, in blind runs with `DIP,` off, at 3.8 and at 20 (checks below), what the correction changes agrees to 0.5 %
of the change (0.2 % with GOSIA's E_cm), from 0.25 % to 46 % of a probability. What the default changes on the templates,
against `e1_polarisation = 0` (populations and γ rays integrated over all angles, excitations of the plan's rates, and B(E2)
read from one sample made with the correction and analysed with and without it):

| | z | Populations | γ rays | Plan's excitations | B(E2) read without it |
|---|---|---|---|---|---|
| `coulex_ni58`: ⁵⁸Ni, one state, first order | 0.0064 | 2⁺ −0.57% | −0.57% | −0.75% CD, −0.53% strips | −0.60% |
| `coulex_ni62_pt194` (and `_clovers`): ⁶²Ni, the beam | 0.053 | 2⁺₁ −3.8%, 0⁺₂ −8.7%, 2⁺₂ −6.4%, 4⁺ −7.9% | 1172.9 keV −4.0%, 1163.5 keV −7.9% | −5.1% CD, −4.8% strips | −7.2% (clovers −7.1%) |
| the same: ¹⁹⁴Pt, the target | 0.021 | 2⁺₁ −0.6%, 2⁺₂ −1.9%, 4⁺₁ −2.1%, 6⁺ −4.0% | 328.5 keV −0.9% | — | — |
| `coulex_ni62_c_allmond`: ⁶²Ni on ¹²C, one state | 0.0044 | 2⁺ −0.40% | −0.40% | −0.53% ring 2, −0.47% ring 3 | −0.49% |
| `coulex_ru98_plunger`: ⁹⁸Ru, one state | 0.011 | 2⁺ −0.70% | −0.70% | −0.72% | — |

B(E2) moves more than the yield for ⁶²Ni because the excitation saturates. With d = 0.005 every change is about a third
larger. The reorientation of the ⁶²Ni 2⁺ state (Q = +5 e fm²; it raises the population at 170° by 8%) keeps its sign and
loses one to two percent of itself, since the diagonal coupling is weakened with the transitional one; but the ratio of the
2⁺ populations at 160° and at 60° falls by 1.5%, as a smaller diagonal element would make it. That is why quadrupole moments
are where it matters.

**Stored runs.** A run's `setup.toml` says which d its events have, written out when the run is taken even when the setup
sets none. A run stored before physim had the correction says nothing; it is loaded with d = 0, which is what its events
have, so that its analysis, its simulated peak and its continuation stay without it. Compared with a current setup that
has the default, the change is reported ("reaction.e1_polarisation: 0.0 → 0.0038").

## Beyond first order: the events of a setup with a level scheme

First order excites each level in one step from the ground state. A 4⁺ state has no E2 matrix element to a 0⁺
ground state, so first order never reaches it, and the same holds for an excited 0⁺ state. They are reached in two
steps through the 2⁺ state. `physim.nuclear.coupled.CoupledChannels` solves the coupled equations for the
amplitudes of every magnetic substate of every level along the orbit, which holds all orders: the two-step
excitations, the depletion of the ground state and of the 2⁺ state, and the reorientation through the diagonal
matrix elements (the quadrupole moments).

**The orbit of each pair of levels.** The coupling of levels n and m takes the orbit symmetrised between the speeds
the nucleus has in the two (Alder and Winther), a_nm = Z₁Z₂e²/(μ v_n v_m) and v_nm = √(v_n v_m), with
½μv_n² = E_cm − E_n:

$$
K_{\lambda,nm} = \frac{4\pi Z e^2}{\hbar\, v_{nm}\,(2\lambda+1)\, a_{nm}^\lambda}
\;\propto\; \big[(E_{cm}-E_n)(E_{cm}-E_m)\big]^{(2\lambda-1)/4} ,
$$

and ξ_nm = (Z₁Z₂e²/ħ)(1/v_n − 1/v_m). It is GOSIA's ψ (its LOAD: (pp₁ pp₂)^((2λ−1)/4), pp = E − (1 + A₁/A₂)E_n,
which is (1 + A₁/A₂)(E_cm − E_n) with A₁ the projectile's mass, whichever nucleus is excited), and the orbit first
order gives each level. A diagonal element acts at the speed the nucleus has in its own level, v_n. The
eccentricity 1/sin(θ/2) is the same for every pair.

Until 2026-11-24 every pair had the orbit of the ground state and the reference level (the first 2⁺), which is
too strong by [E_cm(E_cm − E_ref) / ((E_cm − E_n)(E_cm − E_m))]^((2λ−1)/4) for the others. The excitation of a level
in one step in a scheme of two levels without a quadrupole moment was right. A reorientation effect was too large by
(E_cm/(E_cm − E*))^¾ (6% for ⁶²Ni at 111.6 MeV on carbon, 5% for ¹⁶O at 30 MeV on ⁵⁸Ni, 0.5% for ⁶²Ni at 237 MeV on
¹⁹⁴Pt), and the populations of the levels above by more (12 to 26% for ⁶²Ni at 100 MeV on carbon).

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
  density matrix ρ of its level down the cascade. A γ ray from spin I to I_f by multipoles L < L′ (mixing ratio
  δ, ENSDF's: the amplitude of L′ over L with Krane and Steffen's sign, whatever order the multipolarity is written
  in; E0 has no γ ray and is left out) is drawn along n with the probability W(n) ∝ Σ_σ Tr[A_σ ρ_n A_σ†], where
  ρ_n is ρ quantised along n, σ = ±1 the photon's helicity and (A_σ)_{mM} = Σ_L c_L √(2L + 1) ⟨I_f m L σ | I M⟩
  with c_L = 1, c_L′ = −σδ: for either sign of δ this is Krane and Steffen's A_k = [F_k(LL) + 2δ F_k(LL′) +
  δ² F_k(L′L′)]/(1 + δ²) (checked against their F coefficients as exact surds for 2 → 2 and 4 → 4 lines, and against
  GOSIA's yields of four mixed lines, which miss by up to 500 bands with δ's sign flipped). The level
  it feeds is left in Σ_σ A_σ ρ_n A_σ† rotated back and normalised, and the next γ ray starts from that. Summed
  over directions this is the unobserved feeding above, so the singles are unchanged and only the γ–γ
  coincidences between crystals change: for an unoriented 4⁺ → 2⁺ → 0⁺ cascade the two γ rays follow
  1 + 0.1020 P₂ + 0.0091 P₄ of the angle between them. A converted transition hands its orientation on as an
  unobserved one, and `emission = "isotropic"` draws every γ ray on its own, without the correlation.
- **First order stays available:** `[reaction] order = "first"` gives the events of first-order theory (bit for
  bit those of earlier versions). A setup without a level scheme is first order in any case, since its one matrix
  element allows nothing else.

**Through a thick target** (one the generator cuts into slices) the excitation probability and the orientation are
those of each event's own beam energy: the equations are solved at 3 to 8 beam energies through the target, the
probability is interpolated between them in ln P, and the density matrix of the level the excitation starts in by the
polynomial through them (`gamma.DepthOrientation`), as the correlation table's coefficients are. A matrix that the
polynomial leaves with an eigenvalue below −10⁻⁴ (close to 0°, or at the back of a very thick target; none of the
⁶²Ni-on-carbon events) has its negative eigenvalues set to zero. Which level the excitation starts in is still drawn from
the populations in the middle of the target. A target of one slice, a lifetime measurement and the other nucleus have the
middle's solution alone.

The solution is made once per setup and beam energy (in the middle of the target, and the others of a thick one), on a
grid of 90 CM angles (2° apart) with a relative accuracy of 10⁻⁴ in the probabilities, and kept. The equations of all
angles are integrated together (`CoupledChannels.solve_many`), and the forward angles at which first order excites less than
10⁻⁶ of the collisions are left in first order (below 14° in the ⁶²Ni example, where no population changes by
more than 3 × 10⁻⁹). The table then takes about 2 s for the six levels of ⁶²Ni, where one angle after the other
takes 16 s; first order takes 0.2 s.

For `coulex_ni62_pt194` (⁶²Ni at 237 MeV on ¹⁹⁴Pt, the beam excited), at θ = 150° in the CM frame:

| Level of ⁶²Ni | First order | All orders | First order, no E1 polarisation | All orders, no E1 polarisation |
|---|---|---|---|---|
| 2⁺, 1172.9 keV | 0.1826 | 0.1748 | 0.1931 | 0.1836 |
| 0⁺, 2048.6 keV | 0 | 0.0073 | 0 | 0.0080 |
| 2⁺, 2301.9 keV | 0.0005 | 0.0015 | 0.0006 | 0.0017 |
| 4⁺, 2336.5 keV | 0 | 0.0052 | 0 | 0.0057 |
| any level | 0.1831 | 0.1889 | 0.1936 | 0.1990 |

The 2⁺ state loses 4.3% to the levels above it and gets it back, with more, when they decay through it: 2.9% more
2⁺ → 0⁺ γ rays than first order gives, 4⁺ → 2⁺ γ rays (1163.5 keV) at 3% of them and 0⁺ → 2⁺ γ rays (875.6 keV)
at 4%.

**One particle per collision.** A collision either excites the nucleus or leaves it as it was, so elastic
scattering is the Rutherford cross section times 1 − P(θ): 19% fewer elastic particles at backward angles in the
⁶²Ni example, and elastic and excited particles together are Rutherford scattering (to the inelastic orbit's
(v_i/v_f − 1)·P, a few 10⁻⁴: semiclassical theory conserves the flux only on the symmetrised orbit).

**What this leaves out.**

- Every excitation has the kinematics of the `[reaction]` state (Q = −E* of that state), whichever level it
  reached: with ⁶²Ni left in its 4⁺ state the two nuclei share 1.2 MeV less kinetic energy than the events give
  them.
- The other nucleus of the collision (¹⁹⁴Pt here) is excited as well when the setup has its level scheme, with
  the summed populations of its own coupled equations (95% at backward angles here, 0.85 γ rays of 328.5 keV per
  backward particle). The two excitations are treated as independent: a collision can excite both, and the
  excitation of one does not change the probability of the other. The particle keeps the energy of the
  collision without the other nucleus's excitation.
- First-order theory reads B(E2) too high from events with all orders: the 2⁺ → 0⁺ peak holds the feeding from
  above, and a line the Doppler-corrected peak does not resolve from it (the 4⁺ → 2⁺ line of ⁶²Ni lies 9.4 keV
  below a peak 15 keV wide, and is 3% of its window). The analysis therefore solves the coupled equations of the
  level scheme for the B(E2) that gives the measured yield, 4 to 8% below the first-order reading on the ⁶²Ni
  example, with every other matrix element taken from the scheme. The excitation saturates: the yield grows as
  B(E2) to the power 0.67 there, so the uncertainties of B(E2) are 1.5 times those of the yield. With the ENSDF scheme of ⁵⁸Ni added to `coulex_ni58` (¹⁶O at 30 MeV), all orders give 5% fewer
  2⁺ → 0⁺ γ rays in the CD and 2% fewer in the strip detectors than first order, mostly through the reorientation.
- The coupled equations have no excitation of both nuclei in one collision and no nuclear interference above the safe
  energy. The E1 polarisation correction is GOSIA's approximation, with one parameter d for every nucleus, not a
  calculation from the giant dipole resonance of the nucleus at hand.
- The plan's coincidence rate takes the angular correlation from the same solution as the events: with a level
  scheme the 2⁺ state as the whole scheme leaves it, to all orders. The state alone in first order is more
  sharply aligned: on the ⁶²Ni example it would send 14% fewer γ rays to the crystal at 90° with the CD, and 3%
  more to those at 45° and 135°. Without a level scheme the plan and the events both take the state alone.

## Deorientation in vacuum

A nucleus that leaves the target is an ion stripped of many electrons, and its electron shells are still decaying. The
hyperfine field they make at the nucleus turns the nuclear spin, so the orientation the collision gave fades while the
level lives in vacuum, and the γ rays come out more evenly. physim has GOSIA's two-state model of it
(`physim.nuclear.deorientation`; F. Bosch and H. Spehl, Z. Phys. A 280, 329 (1977); R. Brenn et al., Z. Phys. A 281,
219 (1977); the GOSIA manual (2012), section 4.4; subroutines `GKVAC`, `GKK`, `XSTATIC` and function `ATS` of
`gosia_20110524.14.f`, whose arithmetic it follows). The statistical tensors of rank k of a level that has been in vacuum
for a time t are G_k(t) times what the collision, or the decay that fed it, left. The field first *fluctuates*, turned at
random every τ_c, until at the rate Λ* the ion reaches its ground state and the field is *static*:

$$
G_k(t) = e^{-\Lambda^* t} G_k^\mathrm{fl}(t) + \int_0^t \Lambda^* e^{-\Lambda^* t'} G_k^\mathrm{fl}(t')\,
G_k^\mathrm{st}(t - t')\,dt' ,
\qquad
G_k^\mathrm{st}(t) = \langle\alpha_k\rangle + (1 - \langle\alpha_k\rangle)\,e^{-\Gamma t} ,
$$

with the hard core α_k averaged over the ion's charge states and the field H = K Z β^x. Averaged over the decay of a level
of mean life τ, (1/τ)∫ G_k(t) e^{−t/τ} dt, this is GOSIA's closed form, the G_k GOSIA prints.

- **The parameters** are GOSIA's defaults (its input `VAC,`): atomic spin J = 3, Γ = 0.02 /ps, Λ* = 0.0345 /ps,
  τ_c = 3.5 ps, K = 600 T, x = 0.6. The g factor is Z/A, as in GOSIA, unless the setup gives one:
  `[reaction] deorientation_g`, for every level of the nucleus the reaction names. The other nucleus of the collision
  has Z/A. `[reaction] deorientation = "none"` leaves it out; left unset it is `"vacuum"`.
- **Only in vacuum.** A level is deoriented while it lives outside the target and what is behind it: from where its
  nucleus leaves the stack (in the events each nucleus with its own depth, direction and speed, so its own time of exit
  and speed in vacuum), or from its feeding if it is fed outside. A nucleus that decays in the target, or stops there, is
  not deoriented, and a level without a half-life in the level scheme decays at once and keeps G_k = 1. Times are the
  nucleus's own, the speed that at the exit.
- **Down a cascade** each level is deoriented before it emits and before it feeds the next one (all ranks of its density
  matrix), as GOSIA feeds with the deoriented tensors; a converted step takes the average over its level's life. The plan's
  angular correlation is what the events average to, the feeding in and out of the target in closed form (`gamma.Faded`),
  and so the coincidence rates, the efficiency and correlation of the analysis and the simulated peaks have it too.
  The GOSIA input of `Multistep` writes `VAC,` with the setup's parameters (Λ* = 0 for none).
- **Its size** on `coulex_ni62_pt194`: G₂ of the ⁶²Ni 2⁺₁ (1.44 ps) is 0.92 at β = 0.05, that of the ¹⁹⁴Pt 2⁺₁ (42 ps)
  0.16. The coincidences fall by 1.5% (CD) and 2.9% (strips) for ⁶²Ni and by 8.3% and 10.6% for ¹⁹⁴Pt; with clovers by 5.3%
  and 6.7%, and by 20.3% and 23.0%. B(E2) from samples made with it comes back (12 seeds: +0.39 ± 0.49%; clovers
  −0.22 ± 0.27%); analysed without it the same samples read 2.65% low (clovers 7.7%). The other templates have no half-life
  in a scheme their events follow and do not change. The Physics page of the app gives G₂ and G₄ of every level with a
  half-life at v/c 0.02, 0.05 and 0.08.
- **Against GOSIA, blind** (`tests/reference/nuclear/coulex/gosia_gk.csv`, the app's numbers stored before GOSIA was run):
  174 coefficients G₂, G₄, G₆ (k up to 2I) of the ⁶²Ni, ¹⁹⁴Pt, ¹⁵²Sm, ⁴⁸Ti and ⁹⁸Ru levels at v/c from 0.02 to 0.08, with
  GOSIA's defaults and with a set in which every parameter was changed: none outside 0.002 + half the last printed digit.
  The largest difference, 0.00058, is GOSIA's g = Z/mass (the mass typed in its input) against physim's Z/A; where `VAC,`
  gives g the two agree to the last digit.
- **Where physim departs from GOSIA:** g = Z/A with the mass number; every speed has its own ⟨α_k⟩ (GOSIA works it out
  once in a run, for the first experiment); g = 0 is no interaction (GOSIA's formula gives 0/0).
- **Not included:** the fields in solids (a nucleus in the target, its backing or a stopper), the transient field of a
  nucleus crossing a magnetised foil, and the deorientation of fusion-evaporation residues. A lifetime measurement has its
  own: a plunger only the static model of its `lifetime.deorientation` ({doc}`lifetimes`), or none; DSAM none, since the
  nuclei stop. With a `[lifetime]` section the reaction's setting does nothing, and the app and the summary say so.
- **Stored runs.** A run's `setup.toml` says which deorientation its events have, written out when it is taken. A run
  stored before physim had it says nothing and is loaded with `"none"`, also when it is continued, taken again, compared
  and for its simulated peak; against the current setup the change is reported ("deorientation in vacuum: none →
  vacuum").

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
- **Cross section:** Rutherford's at each depth times v_i/v_f times P(θ) (the cross section above). P follows the
  beam energy through the layer (`rates.depth_probability`; a thin layer keeps the one P at its middle). With a
  level scheme it is the probability of exciting any level (previous section).
- **Event generator:** the excitation channel is picked as often as the elastic one, the masses include E*, and
  each event's weight is multiplied by P(θ*) × v_i/v_f, averaged over its slice of the target, from a table every
  0.25°. The rare inelastic events therefore get as many samples as the elastic ones.
- **Expected peaks:** Coulomb-excitation peaks are always listed, however weak.

The app has an **Excitation and γ rays** tab (P(θ), events per detector, the Doppler table). The report adds a
Coulomb-excitation section and `gamma.csv`.

For `coulex_ni58` (0.1 pnA of 30 MeV ¹⁶O on 0.5 mg/cm² ⁵⁸Ni): the CD records 1.32 excitations per second next to
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
| E1 polarisation strength z, target and projectile excitation | GOSIA's subroutine `LOAD` | 0.2% (E_cm) | `test_z_is_gosias_for_target_and_projectile_excitation` |
| E1 polarisation at 180°, ξ → 0 | (1 − 0.4 z)² from the integrals above | 1e-6 (first order), 1e-5 (coupled equations) | `test_its_size_in_sudden_backscattering_is_the_hand_worked_one` |
| E1 polarisation with a diagonal element | Second-order perturbation theory by nested quadrature | 2e-3 | `test_the_diagonal_coupling_has_it_as_second_order_theory_does` |
| E1 polarisation against GOSIA, blind: ⁶²Ni (Q_s = ±0.3 e b) on ¹²C (either nucleus the beam) and ¹⁹⁴Pt, ¹⁵²Sm on ¹²C, ¹⁹⁴Pt by ⁶²Ni, ⁴⁸Ti on ²⁰⁸Pb | GOSIA 20110524.14 (nuclear masses; `DIP,` 0.001, 3.8, 20), the app's numbers stored before the runs, `coulex/gosia_e1.csv` | Each level at 3.8: 1% above 10⁻⁴, 3% above 10⁻⁶ (worst 0.48%); the change at 3.8 and 20 (above 0.2%): 2% of it (worst 0.5%) | `test_it_is_gosias_in_blind_runs` |
| B(E2) from a run made with the E1 polarisation | The value put in, and 7% less when the analysis leaves it out | 4σ; 1% | `test_the_events_are_made_with_it_and_the_analysis_gives_the_b_e2_back` |
| A run stored before the correction | Read without it, B(E2) of its own events back | 4σ | `test_a_run_stored_before_the_correction_is_read_without_it` |
| Against GOSIA, first order, ¹⁶O + ⁵⁸Ni at 30 MeV, five angles, E1 polarisation in both | GOSIA 20110524.14 with the matrix element ÷10 (weak coupling), `tests/reference/nuclear/coulex_gosia.csv` | 2% (worst 0.23%) | `validation.coulex_vs_gosia` |
| Each pair of levels on its own orbit | K_λ ∝ [(E_cm − E_n)(E_cm − E_m)]^((2λ−1)/4) for E1, E2, E3, either nucleus excited; first order for every level where the coupling is weak | 1e-12; 1e-4 | `test_each_pair_of_levels_is_coupled_on_its_own_symmetrised_orbit`, `test_every_level_is_first_order_where_the_coupling_is_weak` |
| Against GOSIA, all orders, the same case | GOSIA 20110524.14, `coulex_gosia_all_orders.csv` | 1% (worst 0.24%) | `test_all_orders_agree_with_gosia` |
| Against GOSIA, all orders, feeding and reorientation: ⁶²Ni with Q_s = ±0.3 e b, 0⁺₂, 2⁺₂, 4⁺₁ on ¹²C (100 MeV, either nucleus the beam) and ¹⁹⁴Pt (200 MeV), a rotor ¹⁵²Sm on ¹²C (400 MeV) | GOSIA 20110524.14 (nuclear masses, E1 polarisation off), the app's numbers stored before the runs, `coulex/gosia_pairs.csv` | 1% above 10⁻⁴, 3% above 10⁻⁶ (worst 0.5%; 0.01 to 0.2% with GOSIA's non-relativistic E_cm) | `test_feeding_and_reorientation_agree_with_gosia` |
| Deorientation in vacuum: GOSIA's defaults, atomic spins and charge states | GOSIA's `VAC,` defaults, function `ATS`, subroutine `XSTATIC` | exact | `test_the_defaults_are_gosias`, `test_the_atomic_spins_are_gosias_function_ats`, `test_the_charge_states_are_gosias` |
| Deorientation: the average over the decay | GOSIA's closed form, and the time-dependent G_k(t) averaged by quadrature | 1e-6 | `test_the_average_over_the_decay_is_gosias_closed_form` |
| Deorientation against GOSIA, blind: 174 G_k of ⁶²Ni, ¹⁹⁴Pt, ¹⁵²Sm, ⁴⁸Ti, ⁹⁸Ru levels, v/c 0.02 to 0.08, default and changed parameters | GOSIA 20110524.14, `tests/reference/nuclear/coulex/gosia_gk.csv` | 0.002 + half the last digit (worst 0.00058) | `test_it_is_gosias_in_blind_runs` |
| Deorientation in the events against the plan | `correlation_table` with the feeding in and out of the target; decays drawn one by one | 4σ | `test_the_events_of_one_state_follow_the_plan`, `test_the_feeding_is_what_drawn_decays_average_to` |
| B(E2) from deoriented events | The value put in; 2.65% low when the analysis leaves it out | 4σ | `test_the_analysis_gives_b_e2_back_from_deoriented_events` |
| A run stored before deorientation | Read with `"none"`, and stays so when continued | exact | `test_a_run_stored_before_is_read_without_it_and_goes_on_without_it` |
| Against a published measurement | Requested from the user (`tests/reference/nuclear/pending/coulex_literature.md`) | — | not yet written |
