# Evaporation of the compound nucleus

`physim.nuclear.evaporation.evaporate` follows the statistical decay of a hot compound nucleus, one random chain per
event, as PACE4 does (Gavron 1980). The compound nucleus boils off neutrons, protons and α particles, emits
statistical γ rays and, if it is heavy, may fission. The input is the compound nucleus, its excitation energy E* in
MeV and its spin distribution, usually the partial-wave cross sections σ_ℓ in mb from the fusion model. The output,
event by event, is the exit channel, every emitted particle's energy and direction, the residue's recoil, and the
entry state (E*, J) from which the residue's discrete γ cascade starts.

## One step of the chain

A nucleus at excitation energy $E$ and spin $J$ can emit a particle $x$ (n, p, α), a γ ray, or fission. In the
Hauser–Feshbach picture (Hauser and Feshbach 1952) the width for emitting $x$ with kinetic energy $\varepsilon$ and
orbital angular momentum $\ell$, leaving the daughter at spin $J'$, is

$$
\Gamma_x(E, J) = \frac{1}{2\pi\rho(E, J)} \sum_{\ell}\sum_{J'=|J-\ell|}^{J+\ell} \int g_x\,T_\ell^x(\varepsilon)\,
\rho_d(E - S_x - \varepsilon, J')\,d\varepsilon ,
$$

with $S_x$ the separation energy from AME2020 nuclear masses, $T_\ell$ the transmission coefficient, $\rho_d$ the
daughter's level density and $g_x = 2s+1$ (2 for n and p, 1 for α). The particle's own spin is not coupled; its
spin states are counted by $g_x$. The γ-ray width has the same form with $S = 0$, $\ell = 1$ and $T_{E1}$; fission
is in the next sections. All widths share the factor $1/2\pi\rho(E, J)$, so the branching never needs the parent's
own level density.

At each step the code draws the decay from these widths, then the energy $\varepsilon$, then $\ell$, then $J'$, each
from its exact conditional weight. The widths are tabulated once per nucleus on a grid of $E$ (0.5 MeV) and $J$
(1 ħ), for the states the chains actually reach. Inside a 0.5 MeV bin, $\varepsilon$ is drawn by rejection against
the exact integrand, so the spectra have no grid structure. Energy is conserved exactly: the daughter keeps
$E - S_x - \varepsilon$, whatever grid point was used to choose the decay.

The emission is isotropic in the emitter's frame. The two-body break-up is relativistic, and each particle is boosted
to the compound nucleus's frame with the emitter's velocity. So the residue's recoil is the exact sum of all the
kicks, and in every event the particles plus the residue carry $M_\text{CN} + E^*$ and zero momentum.

The chain ends when particles have less than $10^{-3}$ of the total width (PACE4 has a similar discriminator,
0.002), or when the state sits on the yrast line where no statistical γ ray is possible. The state after the last
particle is the entry state. Statistical γ rays emitted between particles are part of the event; those after the
last particle belong to the discrete cascade that the event generator starts from the entry state.

## Level densities

The density of states with angular-momentum projection $M$ is a back-shifted Fermi gas measured from the
rigid-rotor yrast energy of $M$:

$$
\rho_M(E, M) = \frac{(\hbar^2/2I)^{1/2}}{12}\,\frac{e^{2\sqrt{aU}}}{(U + t_0)^{3/2}},\qquad
U = E - \delta - \frac{\hbar^2 M^2}{2I}.
$$

The density of levels of spin $J$ is $\rho(E, J) = \rho_M(E, J) - \rho_M(E, J+1)$. For $U \gg \hbar^2J^2/2I$ this is
the familiar $\rho(E,J) \propto (2J+1)\exp(-J^2/2\sigma^2)\,\rho(E)$ with $\sigma^2 = IT/\hbar^2$ (Bohr and Mottelson
1969, appendix 2B). The same form gives the yrast line: no level of spin $J$ lies below $\hbar^2J^2/2I$. The sum over
$J'$ in the width telescopes to two terms, $\rho_M(E', |J-\ell|) - \rho_M(E', J+\ell+1)$, which keeps the tables cheap.

- Level-density parameter $a = A/10$ MeV⁻¹, PACE4's default. It suits excitation energies of tens of MeV. The
  low-energy value is nearer $A/8$, and no shell correction is applied.
- Back shift $\delta$: $2\Delta$ for even-even nuclei, $\Delta$ for odd-A and 0 for odd-odd, with $\Delta = 12/\sqrt A$
  MeV (the smoothed pairing energies of Gilbert and Cameron 1965).
- Moment of inertia: a rigid sphere, $I = \tfrac25 m_u A (1.2\,A^{1/3}\ \text{fm})^2$. For ¹⁷²Yb this gives
  $\hbar^2/2I = 6.8$ keV, about PACE4's yrast line (yrast spin 105 ħ at 71 MeV, that is 6.4 keV).
- Below the matching energy $U_x = 2.5 + 150/A$ MeV, Gilbert and Cameron's constant temperature,
  $\rho_M \propto \exp(U/T)$, takes over, matched in value and slope, down to the yrast line. The pure Fermi gas has far
  too few levels a few MeV above the ground state. That made γ rays win over the last neutron: Γ_n/Γ_γ about 2 at
  3 MeV above the neutron threshold of ¹⁶⁸Yb, where the constant temperature gives about 5.
- $t_0 = 1$ MeV keeps the prefactor finite next to the yrast line.

## Transmission coefficients

**Neutrons: the black nucleus.** Blatt and Weisskopf (1952, ch. VIII). An ingoing wave only, with wave number
$K = \sqrt{k^2 + k_0^2}$ inside the radius $R = 1.3\,A^{1/3}$ fm and $k_0 = 1$ fm⁻¹, gives

$$
T_\ell = \frac{4 P_\ell\, KR}{(KR + P_\ell)^2 + S_\ell^2},
$$

with the penetration $P_\ell$ and shift $S_\ell$ of free spherical waves at $\rho = kR$, from the usual recursion
($P_0 = \rho$, $S_0 = 0$). At low energy $T_0 \to 4k/K \propto \sqrt\varepsilon$, the 1/v law. The inverse cross
section $\pi\bar\lambda^2 \sum (2\ell+1)T_\ell$ is about 2.1 b for $A = 168$, nearly flat from 0.5 to 10 MeV. A Numerov
solution of the Wilmore–Hodgson optical model, the one PACE4 uses, gives 3.4 b at 0.5 MeV falling to 2.6 b at 5 MeV.

**Protons and α particles: the barrier of the real optical potential.** The potential is Coulomb (a uniform sphere
of radius $R_c$) + a real Woods–Saxon + $\hbar^2\ell(\ell+1)/2\mu r^2$:

- protons: Becchetti and Greenlees (1969) at zero energy, $V = 54.0 + 0.4Z/A^{1/3} + 24(N-Z)/A$ MeV,
  $r_0 = 1.17$ fm, $a = 0.75$ fm, $r_c = 1.3$ fm;
- α particles: Huizenga and Igo (1962), $V = 50$ MeV, $R = 1.17A^{1/3} + 1.77$ fm, $a = 0.576$ fm,
  $R_c = 1.17A^{1/3}$ fm (the same parameters as PACE4).

For each $\ell$ the barrier top $V_\ell$ and its curvature $\hbar\omega_\ell$ are found outside $R$. Above the top
the parabolic barrier of Hill and Wheeler (1953) gives $T_\ell = [1 + \exp(2\pi(V_\ell - \varepsilon)/\hbar\omega_\ell)]^{-1}$.
Below it Kemble's WKB form $T_\ell = [1 + \exp 2\theta]^{-1}$, with
$\theta = \int\sqrt{2\mu(V - \varepsilon)}/\hbar\,dr$ across the barrier, keeps the long Coulomb tail. A parabola
there overestimates the tunnelling by many orders of magnitude: cold heavy nuclei then shed strings of 2 MeV α
particles. Both forms give ½ at the top. For α on ¹⁵⁶Gd the s-wave barrier is 17 MeV at 10.4 fm.

## Statistical γ rays: giant-dipole E1

$T_{E1}(\varepsilon) = 2\pi\varepsilon^3 f_{E1}(\varepsilon)$, with the Brink–Axel Lorentzian

$$
f_{E1}(\varepsilon) = 8.674\times10^{-8}\ \frac{\sigma_0\,\Gamma_G^2\,\varepsilon}{(\varepsilon^2 - E_G^2)^2 +
\varepsilon^2\Gamma_G^2}\ \ \text{MeV}^{-3}\quad(\sigma_0\ \text{in mb, energies in MeV}),
$$

$E_G = 31.2A^{-1/3} + 20.6A^{-1/6}$ MeV (Berman and Fultz 1975), $\Gamma_G = 5$ MeV, and $\sigma_0$ from the
Thomas–Reiche–Kuhn sum rule, $\int\sigma\,dE = 60NZ/A$ mb MeV. The resulting radiative width at the neutron
separation energy of ¹⁶⁸Yb is 0.1–0.2 eV; measured values for rare-earth nuclei are 50–100 meV. There is no E2 or
M1: collective E2 strength near the yrast line is left to the discrete cascade.

## Fission: Bohr–Wheeler with a liquid-drop barrier

$$
\Gamma_f(E, J) = \frac{1}{2\pi\rho(E, J)}\int_0^{E - B_f(J)} \rho_s(E - B_f(J) - \varepsilon, J)\,d\varepsilon
$$

(Bohr and Wheeler 1939), with the saddle-point level density $\rho_s$ of the same form, $a_f = a$. The barrier at
zero spin is the liquid drop's: fissility $x = (Z^2/A)/[50.88(1 - 1.7826I^2)]$, surface energy
$E_s = 17.9439(1 - 1.7826I^2)A^{2/3}$ MeV with $I = (N-Z)/A$ (Myers and Swiatecki 1966), and

$$
B_f(0) = 0.38(0.75 - x)E_s\ \ (x \le 2/3),\qquad B_f(0) = 0.83(1 - x)^3 E_s\ \ (2/3 < x < 1),
$$

a fit to the saddle energies of Cohen and Swiatecki (1963). For ¹⁷²Yb this gives 30.3 MeV. Sierk's barrier, which
PACE4 prints, is 30.04 MeV. The barrier falls with spin because the saddle shape has the larger moment of inertia:
the saddle sits at $B_f(0) + \hbar^2J^2/2I_s$ with $I_s = 2I$. That is a simple, documented stand-in for Sierk's
rotating-liquid-drop barriers. No shell correction is applied. For ²²⁴Th at 50 MeV, 44% of the chains fission at
spin 0 and 79% at spin 20 ħ.

## Checks

The tests (`tests/python/test_nuclear_evaporation.py`) check:

- **Conservation.** In every event, energy is conserved to 0.1 keV of 1.5 × 10⁵ MeV and momentum to 10⁻⁶ MeV/c. Each
  residue matches the particles emitted.
- **Weisskopf's spectrum.** With only neutrons allowed, the first neutron from ¹⁶⁸Yb at 30 MeV and spin 0 follows
  $\varepsilon\,\sigma_\text{inv}(\varepsilon)\,\rho(U - \varepsilon)$ from 0.5 to 8 MeV (χ²/dof < 2). From 0.5 to
  4 MeV it follows $\varepsilon\,\sigma_\text{inv}\exp(-\varepsilon/T)$, with $T$ within 5% of the slope of ln ρ and
  within 15% of $\sqrt{U/a}$.
- **Reproducibility.** Seeded runs are bit-identical. The channel shares sum to one, and the cross sections to σ_fus.
- **Fission** takes over in ²²⁴Th, grows with spin, and is absent in ¹⁷²Yb.
- **PACE4.** 12C + 150Nd at 70 MeV (E* = 59.3 MeV in ¹⁶²Dy, PACE4's own σ_J): see below.

The validation suite (`physim.nuclear.validation`) holds the comparisons for all five systems below, and two exact
checks: Weisskopf's spectrum as above (χ²/dof 1.1 against the full form, 0.5–8 MeV) and conservation in every
event.

### Comparison with PACE4

PACE4 4.36.0 (LISE++) was run in silent mode with its defaults, 20 000 cascades each, by
`scripts/pace4_reference.py`, which saves each run as `tests/reference/nuclear/pace4_<beam>_<target>_<E>MeV.csv`
(rerun on 15 November 2026; the residue shares moved within PACE4's statistics). Physim used PACE4's E* and σ_J,
20 000 events. Shares of the three strongest PACE4 residues, then mean CM energies (MeV):

| System (E*) | PACE4 | physim | ⟨ε_n⟩ PACE4 / physim | ⟨ε_p⟩ | ⟨ε_α⟩ |
|---|---|---|---|---|---|
| ⁴He + ¹⁶⁵Ho, 50 MeV (47.6) | ¹⁶⁵Tm 0.90, ¹⁶⁶Tm 0.051, ¹⁶⁵Er 0.015 | 0.85, 0.041, 0.016 | 2.93 / 2.22 | 9.54 / 9.31 | 17.0 / 17.9 |
| ¹²C + ¹⁵⁰Nd, 70 MeV (59.3) | ¹⁵⁷Dy 0.76, ¹⁵⁸Dy 0.147, ¹⁵⁶Dy 0.053 | 0.65, 0.142, 0.153 | 2.80 / 2.43 | 10.2 / 9.65 | 17.6 / 17.9 |
| ¹⁶O + ¹⁴⁴Sm, 90 MeV (52.5) | ¹⁵⁶Yb 0.31, ¹⁵⁷Yb 0.20, ¹⁵⁴Er 0.148 | 0.20, 0.22, 0.12 | 2.77 / 2.27 | 10.1 / 9.2 | 17.9 / 18.0 |
| ³²S + ⁹²Mo, 150 MeV (63.4) | ¹²⁰Ba 0.23, ¹²⁰Cs 0.117, ¹¹⁷Cs 0.088 | 0.023, 0.13, 0.002 | 3.16 / 2.84 | 9.34 / 8.49 | 16.8 / 16.0 |
| ⁴⁸Ca + ¹²⁴Sn, 200 MeV (71.0) | ¹⁶⁶Yb 0.49, ¹⁶⁷Yb 0.40, ¹⁶⁴Er 0.033 | 0.24, 0.28, 0.017 | 3.08 / 2.56 | 10.6 / 10.0 | 18.8 / 18.9 |

The checks `evaporation_channels_vs_pace4` and `evaporation_spectra_vs_pace4` hold item 75's targets, fixed before
these runs: the strongest three residues the same as PACE4's, in the same order, each within 30%; ⟨ε_n⟩ and ⟨ε_p⟩
within 10%. Both fail and are known failures of the suite:

- **Channels: no system passes.** For ⁴He and ¹²C the two strongest residues agree within 20%, but physim's third is
  the next xn channel (¹⁶⁴Tm, 5n, share 0.081; ¹⁵⁶Dy, 6n, 0.153) instead of PACE4's. For ¹⁶O, physim puts ¹⁵⁷Yb
  (3n) above ¹⁵⁶Yb (4n) and ¹⁵⁷Tm (p2n) third. For ³²S and ⁴⁸Ca see below.
- **Energies: protons pass for all five** (worst −9%, ¹⁶O and ³²S). **Neutrons are 10–24% low** in all five (worst
  ⁴He + ¹⁶⁵Ho, −24%). α energies (not judged) agree within 5%.

What agrees and what does not:

- **Charged particles.** α energies agree within 5% everywhere, protons within 10%.
- **Neutrons** come out 10–25% softer than PACE4's. The first neutron of a chain does follow $2T$ of the model's own
  level density (3.06 MeV against the direct integral 3.04 MeV for ¹⁶⁹Tm at 47.6 MeV). The difference is in the
  later, colder steps, and so in the models, not the sampling. Softer neutrons leave energy for one more neutron in
  some chains, which is why ¹⁵⁶Dy (6n) is three times PACE4's share for ¹²C + ¹⁵⁰Nd.
- **Light and medium systems at moderate spin** (⁴He, ¹²C): the strongest two channels agree within 15% and in order.
- **High spin.** For heavy projectiles physim leaves more energy in rotation at the end of the chain. Entry spins are
  about 45 ħ after ⁴⁸Ca + ¹²⁴Sn, so fewer neutrons are emitted (5n and 4n instead of 6n and 5n). Event by event the
  model is consistent: chains that start below 40 ħ give 6n as in PACE4. PACE4 removes much more spin on the way,
  with strong statistical E2 emission (12 E2 γ rays per event) that physim does not have. Adding a constant-strength
  statistical E2 here made the neutron deficit worse, because it also competes with the neutrons near threshold.
- **³²S + ⁹²Mo** (proton-rich ¹²⁴Ce at up to 69 ħ): physim emits too many protons and α particles and too few
  neutrons (multiplicities 2.5, 1.0 and 0.5 against PACE4's 1.9, 0.6 and 1.4).

The unit test (`test_nuclear_evaporation.py`) asserts the light-system agreement only: for ¹²C + ¹⁵⁰Nd, the
strongest residue the same and within 20%, the second within 30% and in order, ⟨ε_p⟩ and ⟨ε_α⟩ within 10%, ⟨ε_n⟩
within 25%. Item 75's target, held by the validation suite, is all three strongest channels within 30% and in
order, and mean n and p energies within 10%. Today physim meets it for protons, and for the first two channels of
light systems; not for the third channel, neutron energies or heavy, high-spin systems.

## Speed

About 7 000 events per second on one CPU core for 10⁵ events of ¹²C + ¹⁵⁰Nd, ⁴⁸Ca + ¹²⁴Sn or ³²S + ⁹²Mo, including
building the tables. The chains of all events advance together, one step at a time, as NumPy arrays.

## Limits

- No coupling of the particles' intrinsic spins (counted by 2s + 1 instead) or of the beam's and target's
  ground-state spins. Parity is not followed.
- Emission is isotropic in the emitter's frame. The compound nucleus's spin alignment, which makes evaporation
  anisotropic, is not followed.
- No shell corrections in the level densities or the fission barrier, and no collective enhancement.
- Neutron transmission from the black nucleus, not an optical model. Charged-particle barriers from the real
  potential only: no absorption above the barrier, no thermal or deformation lowering.
- Only E1 γ rays (giant dipole). No statistical E2 or M1.
- No deuterons, tritons, ³He or heavier fragments; no pre-equilibrium emission (item 75 leaves it out); fission
  only as a lost share, with no fragments.
- Nuclides missing from AME2020 close the decays that would produce them.

## References

- W. Hauser and H. Feshbach, Phys. Rev. 87, 366 (1952): the statistical (Hauser–Feshbach) widths.
- A. Gilbert and A. G. W. Cameron, Can. J. Phys. 43, 1446 (1965): constant temperature + Fermi gas, pairing energies.
- A. Gavron, Phys. Rev. C 21, 230 (1980): PACE, the Monte Carlo chain this follows.
- A. Bohr and B. R. Mottelson, *Nuclear Structure* vol. 1 (1969), appendix 2B: the Fermi-gas level density with spin.
- J. M. Blatt and V. F. Weisskopf, *Theoretical Nuclear Physics* (1952), ch. VIII: the black nucleus; V. F.
  Weisskopf, Phys. Rev. 52, 295 (1937): the evaporation spectrum.
- F. D. Becchetti and G. W. Greenlees, Phys. Rev. 182, 1190 (1969): the proton potential.
- J. R. Huizenga and G. Igo, Nucl. Phys. 29, 462 (1962): the α potential.
- D. L. Hill and J. A. Wheeler, Phys. Rev. 89, 1102 (1953): the parabolic barrier.
- B. L. Berman and S. C. Fultz, Rev. Mod. Phys. 47, 713 (1975): giant-dipole systematics.
- N. Bohr and J. A. Wheeler, Phys. Rev. 56, 426 (1939): the fission width.
- W. D. Myers and W. J. Swiatecki, Nucl. Phys. 81, 1 (1966): liquid-drop constants; S. Cohen and W. J. Swiatecki,
  Ann. Phys. 22, 406 (1963): saddle-point energies.
