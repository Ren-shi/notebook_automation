# Fusion

`physim.nuclear.fusion.Fusion` gives the fusion barrier of two nuclei, the capture cross section and its
partial waves, the compound nucleus they make, and the excitation function averaged over a target's thickness.
It is the first step of a fusion-evaporation plan: the evaporation code takes the compound nucleus's excitation
energy and spin distribution from here.

## The potential

The nuclei are spheres that interact through the Coulomb force of point charges and the Bass nuclear potential,
plus the centrifugal term of partial wave $\ell$:

$$
V_\ell(r) = \frac{Z_1Z_2e^2}{r} - \frac{R_1R_2}{R_1+R_2}\, g(r - R_1 - R_2) + \frac{\hbar^2\ell(\ell+1)}{2\mu r^2},
\qquad g(s) = \left[A\,e^{s/d_1} + B\,e^{s/d_2}\right]^{-1}.
$$

$\mu$ is the reduced mass from the nuclear AME2020 masses and $e^2 = 1.43996$ MeV fm. Two parameter sets are
available:

| Set | $A$ (fm/MeV) | $B$ (fm/MeV) | $d_1$ (fm) | $d_2$ (fm) | Radii $R_i$ (fm) |
|---|---|---|---|---|---|
| `bass80` (default) | 0.033 | 0.007 | 3.5 | 0.65 | $R_s(1 - 0.98/R_s^2)$, $R_s = 1.28A_i^{1/3} - 0.76 + 0.8A_i^{-1/3}$ |
| `bass77` | 0.0300 | 0.0061 | 3.30 | 0.65 | $1.16A_i^{1/3} - 1.39A_i^{-1/3}$ |

Both are as written out by Dutt and Puri (2010), whose Table II also lists Bass 1980 barriers for 60 systems.
Physim reproduces those within 0.01 MeV in height and 0.02 fm in position.

## Barriers

For each $\ell$ the barrier is the outer maximum of $V_\ell(r)$: its height $V_\ell$, position $R_\ell$ and
curvature

$$
\hbar\omega_\ell = \hbar\sqrt{-V_\ell''(R_\ell)/\mu}.
$$

The s-wave barrier ($\ell = 0$) is *the* fusion barrier $V_B$, $R_B$, $\hbar\omega$. As $\ell$ grows, the barrier
rises and moves inward, and the pocket inside it gets shallower. Above the **critical angular momentum**
$\ell_\text{crit}$ the pocket is gone: the potential falls monotonically with $r$, the nuclei cannot be captured,
and those partial waves do not fuse (Bass's picture of the high-energy limit).

## Cross section

Each barrier is replaced by a parabola of the same height and curvature. Its transmission is the Hill–Wheeler
formula, and the capture cross section is the sum over partial waves:

$$
T_\ell(E) = \frac{1}{1 + \exp[2\pi(V_\ell - E)/\hbar\omega_\ell]}, \qquad
\sigma_\ell = \frac{\pi}{k^2} (2\ell+1)\, T_\ell, \qquad
\sigma_\text{fus} = \sum_{\ell=0}^{\ell_\text{crit}} \sigma_\ell, \qquad k^2 = \frac{2\mu E}{\hbar^2}.
$$

$E$ is the CM kinetic energy, $\sqrt s - m_1 - m_2$. Wong (1973) took the s-wave parabola for every $\ell$, raised by
$\hbar^2\ell(\ell+1)/(2\mu R_B^2)$, and replaced the sum by an integral:

$$
\sigma_\text{Wong}(E) = \frac{\hbar\omega R_B^2}{2E} \ln\left\{1 + \exp\left[\frac{2\pi(E - V_B)}{\hbar\omega}\right]\right\}.
$$

`Fusion.wong_cm` gives this closed form, and `Fusion.partial_waves_cm(E, l_dependent=False)` the explicit sum with
Wong's fixed barrier. The two agree within 1% at all energies. The default, `l_dependent=True`, uses each
$\ell$'s own barrier. It agrees with Wong near the barrier (2% at $V_B + 3$ MeV for ¹⁶O + ²⁰⁸Pb), falls below it
higher up, and saturates at $\pi(\ell_\text{crit}+1)^2/k^2$ far above.

Limits that the tests check:

- **Far above the barrier** ($E - V_B \gg \hbar\omega$), Wong's formula becomes the sharp cut-off
  $\sigma = \pi R_B^2 (1 - V_B/E)$ (`Fusion.classical_cm`).
- **Below the barrier**, $\sigma \approx (\hbar\omega R_B^2/2E)\exp[2\pi(E - V_B)/\hbar\omega]$: each MeV lower costs a
  factor $e^{2\pi/\hbar\omega}$, about 4 for ¹⁶O + ²⁰⁸Pb ($\hbar\omega = 4.3$ MeV).

## The compound nucleus

`Fusion.compound(E_lab)` returns the compound nucleus $(Z_1+Z_2, A_1+A_2)$ with

- its excitation energy $E^* = \sqrt s - M_\text{CN} = E_\text{cm} + Q$, where $Q = m_1 + m_2 - M_\text{CN}$, with
  exact kinematics and nuclear AME2020 masses;
- its velocity $\beta = p_\text{beam}/(E_\text{beam} + m_1 + m_2)$, which is that of the CM frame;
- the array `sigma_l_mb` of $\sigma_\ell$ in mb indexed by $\ell$, and `spin_distribution()`, which is
  $\sigma_\ell/\sigma_\text{fus}$.

With spin-0 ground states the compound nucleus's spin is $J = \ell$. Ground-state spins of beam or target are not
coupled in.

## Thick targets

The beam slows down in the target, so each depth $x$ sees its own energy $E(x)$. Weighting each depth equally
(every target atom counts the same), the yield per target atom is the depth average:

$$
\bar\sigma = \frac{1}{t}\int_0^t \sigma\big(E(x)\big)\,dx = \frac{1}{t}\int_{E_\text{out}}^{E_\text{in}}
\frac{\sigma(E)}{S(E)}\,dE,
$$

with $S$ the stopping power in MeV/(mg/cm²) and $t$ the thickness in mg/cm² (see {doc}`stopping`).
`Fusion.thick_target` evaluates the first form with Simpson's rule over 201 depths, taking $E(x)$ from the range
table. The test evaluates the second form on a fine energy grid; the two agree within 1%. If the beam stops
inside the target, the depth beyond its range counts as giving no fusion. The same average gives the mean
$\sigma_\ell$ over the target, and $E^*$ at each depth.

## Where the model stops being valid

`Fusion.warnings(E_lab)` reports these:

- **Below or near the barrier.** Coupling to the nuclei's low-lying states and to transfer channels spreads the
  single barrier into a distribution of barriers. This raises sub-barrier fusion by up to orders of magnitude over
  a one-barrier estimate (coupled channels; CCFULL would be the outside check). The parabola is also too thin deep
  below the barrier. Physim warns below $V_B$, and near it (up to $V_B + \hbar\omega$).
- **Above about 10 MeV/u.** Incomplete fusion (break-up of the projectile, with only part of it fusing) and
  pre-equilibrium emission take a growing share of the cross section, so complete fusion is overestimated.
- **Heavy systems.** Above $Z_1Z_2 \approx 1600$ the system often re-separates before forming a compound nucleus
  (quasi-fission), and capture overestimates fusion, by orders of magnitude for the heaviest systems. The threshold
  is empirical (Schmidt and Morawek 1991; Swiatecki's extra-push model, 1982). Quasi-fission is also seen in some
  lighter systems with deformed or heavy targets, such as ³²S + ²⁰⁸Pb ($Z_1Z_2 = 1312$), so the warning is a
  guide, not a sharp boundary.
- **At the critical angular momentum.** When partial waves up to $\ell_\text{crit}$ fuse, the cross section is set
  by the pocket of the potential, not by the barrier. The highest spins may fission or end in deep-inelastic
  collisions instead.

Left out: deformation and orientation of the nuclei, coupled channels, incomplete fusion, quasi-fission and the
fusion hindrance deep below the barrier.

## References

- C. Y. Wong, Phys. Rev. Lett. 31, 766 (1973): the closed-form cross section.
- R. Bass, Phys. Rev. Lett. 39, 265 (1977): the Bass 1977 potential.
- R. Bass, *Nuclear Reactions with Heavy Ions* (Springer, 1980): the Bass 1980 potential.
- I. Dutt and R. K. Puri, Phys. Rev. C 81, 064609 (2010), arXiv:1005.5213: both Bass parameter sets written out,
  and Bass 1980 barriers for 60 systems (Table II), used in the tests.
- D. L. Hill and J. A. Wheeler, Phys. Rev. 89, 1102 (1953): transmission through a parabolic barrier.
- K.-H. Schmidt and W. Morawek, Rep. Prog. Phys. 54, 949 (1991); W. J. Swiatecki, Nucl. Phys. A 376, 275
  (1982): fusion hindrance and quasi-fission in heavy systems.
