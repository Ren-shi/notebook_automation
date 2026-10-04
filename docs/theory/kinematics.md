# Two-body reaction kinematics

`physim.nuclear.kinematics.TwoBody` gives the lab energies and angles of the two particles leaving a reaction
$A(a, b)B$: beam $a$ on target $A$ at rest, ejectile $b$ and recoil $B$. Elastic scattering is $b = a$, $B = A$;
inelastic scattering leaves an excitation energy $E^*$ in one of them. Everything is exact special relativity, in
units with $c = 1$ (masses and momenta in MeV).

## Masses and Q-value

Masses are nuclear (bare-nucleus) masses from AME2020 (see the physics register entry for the data), with $E^*$
added to the excited particle:

$$
m_1 = m_a,\quad m_2 = m_A,\quad m_3 = m_b\,(+E^*),\quad m_4 = m_B\,(+E^*), \qquad Q = m_1 + m_2 - m_3 - m_4 .
$$

## Into the centre-of-mass frame

A beam of kinetic energy $T$ has momentum $p_1 = \sqrt{T(T + 2m_1)}$. The invariant mass squared and the CM
velocity are

$$
s = (m_1 + m_2)^2 + 2 m_2 T, \qquad \beta = \frac{p_1}{T + m_1 + m_2}, \qquad \gamma = \frac{T + m_1 + m_2}{\sqrt s}.
$$

The reaction can happen only if $\sqrt s \ge m_3 + m_4$, which gives the threshold

$$
T_\text{th} = \frac{(m_3 + m_4)^2 - (m_1 + m_2)^2}{2 m_2} \approx -Q\left(1 + \frac{m_1}{m_2}\right) .
$$

The outgoing particles share the CM momentum

$$
p^* = \frac{\sqrt{\left[s - (m_3 + m_4)^2\right]\left[s - (m_3 - m_4)^2\right]}}{2\sqrt s},
\qquad E_k^* = \sqrt{p^{*2} + m_k^2}.
$$

The two brackets are evaluated as $Q(m_1 + m_2 + m_3 + m_4) + 2m_2T$ and
$(m_1 + m_2 - m_3 + m_4)(m_1 + m_2 + m_3 - m_4) + 2m_2T$. At low energy $s$ and $(m_3+m_4)^2$ agree to many digits,
and subtracting them directly would lose precision.

## Back to the lab

A particle leaving at CM angle $\theta^*$ is boosted along the beam:

$$
p_\parallel = \gamma p^*(\cos\theta^* + g), \qquad p_\perp = p^*\sin\theta^*, \qquad
E = \gamma(E^* + \beta p^*\cos\theta^*), \qquad g = \frac{\beta E^*}{p^*} = \frac{\beta}{\beta^*},
$$

so $\tan\theta_\text{lab} = \sin\theta^* / [\gamma(\cos\theta^* + g)]$ and the kinetic energy is $E - m$. The
parameter $g$ compares the speed of the CM frame with the particle's speed in it:

- $g < 1$: the particle can go anywhere, and each lab angle has one solution.
- $g > 1$: it is swept forward. There is a **maximum lab angle**,
  $\tan\theta_\text{max} = 1 / (\gamma\sqrt{g^2 - 1})$, and every smaller angle has **two** solutions with different
  energies. This happens to the beam in inverse kinematics (heavy beam on a light target) and near threshold.
- $g = 1$ exactly for every elastic recoil and for elastic scattering of equal masses: the maximum angle is 90°.

The lab angle is turned back into CM angles by writing $\sin\theta^*\cos\theta_\text{lab} -
\gamma\sin\theta_\text{lab}\cos\theta^* = \gamma g \sin\theta_\text{lab}$ as $R\sin(\theta^* - \delta) = \gamma g
\sin\theta_\text{lab}$, which gives both solutions at once. Candidates whose momentum would point backwards are
discarded.

In the non-relativistic limit the elastic energy ratio is the familiar kinematic factor

$$
K = \frac{E'}{E} = \left[\frac{m_1\cos\theta + \sqrt{m_2^2 - m_1^2\sin^2\theta}}{m_1 + m_2}\right]^2 .
$$

## Solid angle and kinematic broadening

Cross sections are calculated in the CM frame and measured in the lab. The conversion factor is

$$
\frac{d\Omega^*}{d\Omega_\text{lab}} = \frac{\sin\theta^*}{\sin\theta_\text{lab}}\frac{d\theta^*}{d\theta_\text{lab}},
\qquad
\frac{d\theta_\text{lab}}{d\theta^*} = \frac{\gamma(1 + g\cos\theta^*)}{\gamma^2(\cos\theta^* + g)^2 + \sin^2\theta^*},
$$

so $\sigma_\text{lab} = \sigma_\text{cm}\, d\Omega^*/d\Omega_\text{lab}$. At $\gamma = 1$ this becomes the textbook
$(1 + g^2 + 2g\cos\theta^*)^{3/2} / |1 + g\cos\theta^*|$. At $\theta^* = 0$ the ratio of sines is replaced by its
limit, $(d\theta^*/d\theta_\text{lab})^2$.

The energy changes across a detector at the rate

$$
\frac{dE}{d\theta_\text{lab}} = \frac{-\gamma\beta p^*\sin\theta^*}{d\theta_\text{lab}/d\theta^*},
$$

reported per degree. Multiplied by a detector's angular size, it gives the kinematic contribution to a peak's
width.

## Checks

The tests (`tests/python/test_nuclear_kinematics.py`) check:

- the lab energies and angles against an explicit Lorentz boost of the CM four-momenta;
- energy and momentum conservation between ejectile and recoil;
- the non-relativistic limit against $K$ to $10^{-6}$;
- the Jacobian and $dE/d\theta$ against numerical derivatives;
- that the lab → CM conversion inverts CM → lab on both branches;
- the threshold of ³H(p, n)³He against its tabulated 1.019 MeV.

Comparisons with LISE++ are listed in the physics register.
