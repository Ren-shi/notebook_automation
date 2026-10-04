# Rutherford scattering

`physim.nuclear.rutherford.Rutherford` gives the elastic Coulomb cross section, the classical orbit geometry
(closest approach, impact parameter), the checks that tell you when the formula no longer applies, and orbits
integrated by physim's own engine.

## The cross section

Two point charges $Z_1e$ and $Z_2e$ scatter along hyperbolas. In the centre-of-mass frame, with kinetic energy
$E_\text{cm}$,

$$
\frac{d\sigma}{d\Omega_\text{cm}} = \left(\frac{Z_1Z_2e^2}{4E_\text{cm}}\right)^2 \frac{1}{\sin^4(\theta/2)},
\qquad e^2 = 1.43996\ \text{MeV fm}.
$$

The lab cross section is the CM one times the solid-angle Jacobian from the kinematics (see {doc}`kinematics`):
$d\sigma/d\Omega_\text{lab} = (d\sigma/d\Omega_\text{cm})\, d\Omega_\text{cm}/d\Omega_\text{lab}$. For the recoil,
the CM angle is that of its partner, $180° - \theta^*$. Where the kinematics are double-valued, each lab angle has
two cross sections. Between two CM angles, over an azimuthal range $\Delta\phi$,

$$
\int \frac{d\sigma}{d\Omega}\,d\Omega = \left(\frac{Z_1Z_2e^2}{4E_\text{cm}}\right)^2 2\,\Delta\phi
\left[\frac{1}{\sin^2(\theta_1/2)} - \frac{1}{\sin^2(\theta_2/2)}\right],
$$

which diverges as $\theta_1 \to 0$: the total Rutherford cross section is infinite.

## Orbit geometry

With $d_0 = Z_1Z_2e^2/E_\text{cm}$ (the head-on distance of closest approach):

$$
b = \frac{d_0}{2}\cot\frac{\theta}{2}, \qquad
r_\text{min} = \frac{d_0}{2}\left(1 + \frac{1}{\sin(\theta/2)}\right) = \frac{d_0}{2} + \sqrt{\frac{d_0^2}{4} + b^2}.
$$

## When Rutherford stops being valid

- **Nuclear forces.** Physim takes an interaction radius $R = r_0(A_1^{1/3} + A_2^{1/3}) + \Delta$, by default
  $r_0 = 1.2$ fm and $\Delta = 2$ fm, just outside the touching nuclear surfaces. Orbits with $r_\text{min} < R$ feel
  the nuclear force. That happens beyond the **grazing angle** $\sin(\theta_\text{gr}/2) = 1/(2R/d_0 - 1)$, and
  never if $E_\text{cm}$ is below the Coulomb barrier $Z_1Z_2e^2/R$; the grazing angle is then reported as 180°.
- **Classical orbits** need a large Sommerfeld parameter $\eta = Z_1Z_2\alpha/\beta$ (β the relative velocity).
  Physim warns below $\eta = 5$.
- **Electron screening** lowers the cross section by $1 - 0.049\,Z_1Z_2^{4/3}/E_\text{cm}[\text{keV}]$ at large
  angles (L'Ecuyer, Davies and Matsunami 1979), and more at small angles. Physim warns above 1%.
- **Identical particles** (Mott scattering) interfere; Rutherford's formula does not apply.

## Integrated orbits

`Rutherford.trajectories(b)` integrates the relative motion with physim's engine. One particle of reduced mass
$\mu$ moves in the field of a fixed charge, with lengths in fm, energies in MeV and time in fm/c, using the adaptive
Dormand–Prince integrator. The orbit is picked up at a large distance $r_0$, with speed and sideways offset chosen so
that the energy is exactly $E_\text{cm}$ and the angular momentum $\mu v_\infty b$. It is therefore exactly the orbit
of impact parameter $b$, not a slightly different one.

The deflection is read from the orbit, not from where the integration stops. Through any point, the hyperbola is
fixed by

$$
e\cos\psi = 1 + \frac{L^2}{\mu k r}, \qquad e\sin\psi = \frac{L\dot r}{k}, \qquad \psi = \phi - \phi_p,
$$

with $k = Z_1Z_2e^2$. Its asymptotes lie at $\phi_p \mp \arccos(1/e)$. Taking the incoming asymptote from the first
point and the outgoing one from the last makes the measured deflection a test of the integration as a whole,
including any drift of the orbit's orientation.
