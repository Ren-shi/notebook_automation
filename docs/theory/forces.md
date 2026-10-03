# Forces

A force adds accelerations $\mathbf a_i$ and, when it is conservative, reports a potential energy $U$ so that
`total_energy()` and recorded energies include it. Built-in forces are exact negative gradients of their potentials
(the test suite checks every one against finite differences), and pairwise forces obey Newton's third law
pair by pair, so total momentum is conserved to round-off. Below, $\mathbf d = \mathbf x_j - \mathbf x_i$,
$r = |\mathbf d|$ and $\hat{\mathbf n} = \mathbf d / r$.

## Gravity

**`NewtonianGravity(G, softening=ε)`**:
$U = -\sum_{i<j} G m_i m_j / \sqrt{r^2 + \varepsilon^2}$. Plummer softening $\varepsilon$ removes the singularity of
close encounters (used for collisionless systems such as clusters; leave it at 0 for planets). A particle of mass 0
feels gravity but does not source it, which makes test particles free. Direct summation costs $O(N^2)$, split into
deterministic parallel blocks.

**`TreeGravity(G, softening, theta, quadrupole)`** (Barnes & Hut 1986): particles are sorted into an octree, each
node storing its mass and centre of mass (and quadrupole tensor if asked). A node of size $s$ at distance $d$ is used
as a single body when $s/d < \theta$; otherwise it is opened. The cost is $O(N\log N)$ and the error is controlled by
$\theta$: at $\theta = 0.5$ the rms force error is about $10^{-3}$. Particles in the same leaf group share one
interaction list (a group walk), which keeps the inner loop vectorisable. Tree forces are not exactly antisymmetric,
so momentum is conserved only to the force error.

**Orbital perturbations** from a central particle, added to Newtonian gravity:
- `PostNewtonian(central, c)`: the first post-Newtonian correction in the test-particle limit (harmonic gauge),
  $\mathbf a = \frac{GM}{c^2 r^3}\left[\left(\frac{4GM}{r} - v^2\right)\mathbf r + 4(\mathbf r\cdot\mathbf v)\mathbf v\right]$.
  It gives the relativistic periapsis advance $6\pi GM/(c^2 a(1-e^2))$ per orbit.
- `J2Oblateness(central, J2, radius, axis)`: $\Phi = \frac{GM J_2 R^2}{2r^3}\left(\frac{3z^2}{r^2} - 1\right)$ with
  $z$ along the spin axis; the central body takes the reaction, so momentum is conserved.

## External potentials (per unit mass)

Acting on every particle about a `center`:

| force | potential $\Phi(r)$ | notes |
|---|---|---|
| `PowerLaw(k, n)` | $k r^n$ ($k\ln r$ for $n=0$) | $n=2$: harmonic; $n=-1$, $k=-GM$: Kepler |
| `Yukawa(k, length)` | $-k e^{-r/\lambda}/r$ | screened Coulomb or nuclear |
| `PlummerPotential(GM, a)` | $-GM/\sqrt{r^2+a^2}$ | cored cluster or halo |
| `HernquistPotential(GM, a)` | $-GM/(r+a)$ | galactic bulge |
| `HarmonicTrap(omega)` | $\tfrac12\sum_k \omega_k^2 x_k^2$ | per-axis frequencies |
| `HenonHeiles(lam)` | $\tfrac12(x^2+y^2) + \lambda(x^2y - y^3/3)$ | chaos test bed; escape at $E = 1/6$ |
| `UniformField(g)` | $-\mathbf g\cdot\mathbf x$ | constant acceleration |

## Springs

- `Spring(i, j, k, L)`: $U = \tfrac k2 (r - L)^2$.
- `AnchorSpring(i, anchor, k, L=0)`: the same to a fixed point.
- `DampedSpring(i, j, k, L, c)`: adds a dashpot along the bond, $c\,(\mathbf v_j - \mathbf v_i)\cdot\hat{\mathbf n}$,
  which damps stretching only.
- `ModulatedSpring(i, to, k, depth, omega)`: $k(t) = k(1 + \text{depth}\cos(\omega t + \varphi))$, the Mathieu
  equation; parametric resonance near $\omega = 2\omega_0/n$.
- `SpringNetwork(i, j, k, L)`: many springs from index arrays, bit-for-bit identical to separate `Spring`s and
  about 1.75× faster.

## Dissipation and driving

- `LinearDrag(gamma)`: $\mathbf a = -\gamma\mathbf v$.
- `QuadraticDrag(c)`: $\mathbf F = -c|\mathbf v|\mathbf v$.
- `PeriodicForce(i, A, omega, phase)`: $\mathbf F = \mathbf A\cos(\omega t + \varphi)$ on one particle.

These have no potential. Pair them with `rk4` or `dopri5`, which handle velocity dependence consistently.

## Electromagnetism

Particles carry charges $q_i$. `ElectricField(E)` gives $\mathbf a = (q/m)\mathbf E$ (potential $-q\mathbf E\cdot\mathbf x$);
`MagneticField(B)` gives $\mathbf a = (q/m)\,\mathbf v\times\mathbf B$, which does no work and has no potential;
`FieldForce(E=f, B=g)` takes both from Python functions of $(t, \text{positions})$. `Coulomb(k, softening)` is the
pairwise $U = k\sum_{i<j} q_iq_j/\sqrt{r^2+\varepsilon^2}$ (like charges repel), sharing gravity's parallel pair loop.
Use the `boris` integrator for magnetic forces.

## Pair potentials for molecular dynamics

`LennardJones(epsilon, sigma, cutoff, shift, box)`: $V = 4\epsilon[(\sigma/r)^{12} - (\sigma/r)^6]$;
`Morse(depth, a, r0)`: $V = D[(1 - e^{-a(r-r_0)})^2 - 1]$; `TabulatedPair`: any $V(r)$ from a table, interpolated by
cubic Hermite polynomials in $V$ and $V'$ so forces are continuous.

- **Cutoff and shift**: pairs beyond $r_c$ are skipped. `shift=True` subtracts $V(r_c)$ so the energy is continuous;
  the force still jumps at $r_c$, which adds a small energy noise.
- **Periodic box**: with `box=L`, separations use the minimum image
  $\mathbf d \leftarrow \mathbf d - L\,\mathrm{round}(\mathbf d/L)$ (requires $r_c \le L/2$). Coordinates are not
  wrapped, so displacements and diffusion are measured directly.
- **Neighbour lists**: a cell list builds Verlet lists of all pairs within $r_c + \text{skin}$ ($\text{skin} =
  0.12\,r_c$), rebuilt only when some particle has moved $\text{skin}/2$. Pairs are re-filtered by $r_c$ at every
  evaluation, so results never depend on when the list was built and restarts stay bit-exact.
- **Virial pressure**: $P = (2K + W)/(3V)$ with the virial $W = \sum_{i<j}\mathbf d\cdot\mathbf F_{ij}$, from
  `w.pressure(V)`.

## Contact

`SoftContact(k, damping, law)` pushes overlapping spheres apart: with overlap $\delta = r_i + r_j - r > 0$ the elastic
force is $k\delta$ (`"linear"`, potential $k\delta^2/2$) or $k\delta^{3/2}$ (`"hertz"`, potential $\tfrac25 k\delta^{5/2}$),
plus a dashpot on the normal approach speed (never attractive). A linear contact lasts $\pi\sqrt{\mu/k}$ with
$\mu$ the reduced mass, so the step must be a small fraction of that. Hard, instantaneous collisions are a different
mechanism; see {doc}`methods`.

## Writing your own

Any `acceleration(t, pos, vel, mass)` in Python becomes a force with `CustomForce`; a Rust closure becomes one with
`ClosureForce`; and a struct implementing the `Force` trait is a first-class built-in. See {doc}`../new-ideas`.
