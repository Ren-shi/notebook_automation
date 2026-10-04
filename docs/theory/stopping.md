# Stopping powers and energy loss

`physim.nuclear.stopping.Stopping` gives the stopping power, range, energy left after a layer and energy straggling
of any ion in any material. Protons and α particles use the NIST tables. Heavier ions, and elements or compounds NIST
does not tabulate, are built from those tables with the models below.

## Protons and α particles: the NIST tables

PSTAR and ASTAR (NIST Standard Reference Database 124, following ICRU Reports 49 and 90) tabulate the electronic
and nuclear stopping power, and the CSDA range, of protons and α particles. They cover 74 materials, from 1 keV to
10 GeV (protons) and 1 GeV (α). Physim reads them as published and interpolates log-log.

- **Isotopes** are taken at the same velocity, since electronic stopping depends on the ion's speed, not its mass: a
  deuteron at energy $E$ stops like a proton at $E\,m_p/m_d$, and ³He like ⁴He.
- **Below 1 keV** the stopping is extrapolated as $S \propto \sqrt E$.
- **Ranges** are $R(E) = \int_0^E dE'/S(E')$, integrated on a fine grid in $\ln E$. For protons and α in a NIST
  material, the range is matched to NIST's own CSDA range at 1 keV, so that the extrapolation below the table
  does not shift ranges at low energies.

## Elements NIST does not tabulate

NIST lists 26 elements. For any other element, the stopping per electron is interpolated between the tabulated
neighbours, linearly in $\ln I$ (the mean excitation energy). Bethe's formula,

$$
-\frac{dE}{dx} = K z^2 \frac{Z}{A}\frac{1}{\beta^2}\left[\ln\frac{2m_ec^2\beta^2\gamma^2}{I} - \beta^2 - \ldots\right],
$$

makes the stopping per electron exactly linear in $\ln I$ at high energy. Near and below the Bragg peak each
element's shell structure matters, and no interpolation captures it (see the register for the measured accuracy).

## Compounds

A compound NIST tabulates uses its table, which includes chemical binding effects. Any other compound or
formula, including isotopic ones such as CD₂, uses **Bragg additivity**: the stopping per atom of each element,
weighted by the number of atoms per formula unit.

## Heavy ions: effective charge

An ion heavier than helium is scaled from the proton at the same velocity:

$$
S_\text{ion}(v) = (\zeta Z_1)^2\, S_p(v),
$$

where $\zeta$ is the fractional effective charge in the Brandt–Kitagawa form with Ziegler's fit for the ionised
fraction $q$ (Ziegler, Biersack and Littmark 1985):

$$
y_r = \max\!\left(0.13,\ \frac{v_r}{v_0 Z_1^{2/3}}\right), \quad
q = 1 - \exp\!\left(0.803y_r^{0.3} - 1.3167y_r^{0.6} - 0.38157y_r - 0.008983y_r^2\right),
$$

$$
\zeta = q + \frac{1-q}{2v_F^2}\ln\!\left[1 + \left(\frac{4\Lambda v_F}{1.919}\right)^2\right],
\qquad \Lambda = \frac{2 \times 0.24005\,(1-q)^{2/3}}{Z_1^{1/3}\left[1 - (1-q)/7\right]},
$$

with velocities in units of the Bohr velocity $v_0$ and the constants as in ZBL's code.

Here $v_r$ is the ion's speed relative to the target electrons. The target Fermi velocity $v_F$ is set to 1 for
all targets, since per-element values are not tabulated here.

At low velocity this charge collapses. There the larger of it and **Lindhard–Scharff** stopping is used, joined to
the full-charge stopping as $1/S = 1/S_\text{LS} + 1/(Z_1^2 S_p)$, with

$$
S_\text{LS} = 8\pi e^2 a_0\,\frac{Z_1^{7/6} Z_2}{\left(Z_1^{2/3} + Z_2^{2/3}\right)^{3/2}}\,\frac{v}{v_0}
\quad\text{per atom.}
$$

## Nuclear stopping

Elastic collisions with whole atoms use the **ZBL universal** stopping. In the reduced energy

$$
\varepsilon = \frac{32.53\,M_2 E}{Z_1 Z_2 (M_1 + M_2)\left(Z_1^{0.23} + Z_2^{0.23}\right)} \quad (E\text{ in keV}),
$$

the stopping is $S_n(\varepsilon) = \ln(1 + 1.1383\varepsilon) / \left[2(\varepsilon + 0.01321\varepsilon^{0.21226}
+ 0.19593\varepsilon^{1/2})\right]$ for $\varepsilon \le 30$, and $\ln\varepsilon / (2\varepsilon)$ above. Protons and
α in NIST materials use the tabulated nuclear stopping instead.

## Energy after a layer

The range–energy relation gives the energy left after a path $x$ exactly:

$$
R(E_\text{out}) = R(E_\text{in}) - x .
$$

A layer tilted by $\theta$ from its normal is crossed along $x = t/\cos\theta$. Near 90° this path grows without
limit, which is what makes LISE++'s "energy at detector entrance" curve dip sharply at 90°.

## Energy straggling

In a thin layer, Bohr's variance is

$$
\Omega_B^2 = 4\pi r_e^2 (m_ec^2)^2 N_A\, z_\text{eff}^2 \frac{Z_2}{A_2}\, t\;\frac{1 - \beta^2/2}{1 - \beta^2}
\approx 0.157\, z_\text{eff}^2 \frac{Z_2}{A_2}\, t \ \text{MeV}^2 \quad (t\text{ in g/cm}^2).
$$

It is reduced for slow ions by Lindhard and Scharff's factor $L(\chi)/2$, with $L = 1.36\chi^{1/2} -
0.016\chi^{3/2}$ for $\chi = v^2/(Z_2 v_0^2) < 3$. Here $z_\text{eff}$ is the effective charge implied by the
stopping.

Through a thick layer, two ions whose energies differ by $\delta$ drift apart as $d\delta/dx = -S'(E)\,\delta$. So
the variance follows

$$
\frac{d\sigma^2}{dx} = -2S'(E)\,\sigma^2 + \frac{d\Omega^2}{dx} .
$$

Above the Bragg peak ($S' < 0$) the spread grows faster than Bohr's estimate; below it, it shrinks. The FWHM of a
Gaussian is $2.355\sigma$.

## Multiple scattering

The angular spread from multiple scattering uses Highland's formula,

$$
\theta_0 = \frac{13.6\ \text{MeV}}{\beta c p}\, z\sqrt{\frac{x}{X_0}}\left[1 + 0.038\ln\frac{x z^2}{X_0\beta^2}\right],
$$

with the radiation length from Dahl's fit, $X_0 = 716.4\,A / [Z(Z+1)\ln(287/\sqrt Z)]$ g/cm². It is accurate to
about 11% for fast particles; for slow heavy ions it is only indicative.
