# Integrators

All integrators advance $\dot{\mathbf x} = \mathbf v$, $\dot{\mathbf v} = \mathbf a(t, \mathbf x, \mathbf v)$ over a
step $h$. They differ in **order** (global error $\propto h^p$), **cost** (force evaluations per step), and
**structure**: whether they preserve the geometry of Hamiltonian flow. `w.integrator_info` reports the order and
whether a scheme is symplectic; the orders below are measured by the test suite and in the first example notebook.

## Why structure matters

For a Hamiltonian system $H = T(\mathbf p) + V(\mathbf q)$, a **symplectic** integrator is the exact flow of a nearby
"shadow" Hamiltonian $\tilde H = H + O(h^p)$. Its energy error therefore stays bounded, oscillating at $O(h^p)$
for exponentially long times, instead of drifting. A non-symplectic method of the same order (RK4) has a smaller
error per step but it accumulates: the energy drifts linearly in time. For long runs of conservative systems
(orbits, molecular dynamics, chaos studies) prefer a symplectic method; for short, dissipative or strongly
velocity-dependent problems, Runge–Kutta methods are simpler and fine.

Explicit symplectic schemes are built from **kicks** $\mathbf v \leftarrow \mathbf v + h\,\mathbf a(\mathbf x)$ and
**drifts** $\mathbf x \leftarrow \mathbf x + h\,\mathbf v$. Each is the exact flow of part of $H$, so any product
of them is symplectic, *provided the force does not depend on velocity*. With drag or magnetic forces they remain
consistent but lose the property; use `gauss*` (any Hamiltonian force) or `boris` (magnetic fields) instead.

## The schemes

| name | order | evals/step | symplectic | use it for |
|---|---|---|---|---|
| `explicit_euler` | 1 | 1 | no | teaching: energy grows |
| `symplectic_euler` | 1 | 1 | yes | the simplest structure-preserving scheme |
| `verlet` | 2 | 1* | yes | the default: MD, springs, constraints (RATTLE) |
| `yoshida4` | 4 | 3* | yes | smooth long runs needing more accuracy |
| `pefrl` | 4 | 4 | yes | like `yoshida4`, ~60× smaller error at equal cost |
| `blanes_moan4` | 4 | 6 | yes | highest accuracy per cost at 4th order |
| `yoshida6`, `yoshida8` | 6, 8 | 7, 15 | yes | very high accuracy on smooth problems |
| `gauss2/4/6` | 2, 4, 6 | iterated | yes | velocity-dependent Hamiltonian forces; stiff-ish problems |
| `wisdom_holman` | 2 | 1 | yes | planetary systems around a dominant mass |
| `boris` | 2 | 1 | volume-preserving | charged particles in magnetic fields |
| `rk4` | 4 | 4 | no | short runs, drag, any velocity dependence |
| `dopri5` | 5 | 6 | no | fixed-step high order; the adaptive method |

\* `verlet` and `yoshida4` reuse the end-of-step acceleration when no force depends on velocity.

### Euler methods

Explicit Euler, $\mathbf x_{n+1} = \mathbf x_n + h\mathbf v_n$, $\mathbf v_{n+1} = \mathbf v_n + h\mathbf a_n$, is first
order and on an orbit spirals outwards: its energy grows every step. Symplectic Euler kicks first and drifts with the
*new* velocity; the same cost and order, but its energy error is bounded.

### Velocity Verlet (leapfrog)

$$
\mathbf v_{1/2} = \mathbf v_n + \tfrac h2 \mathbf a(\mathbf x_n),\quad
\mathbf x_{n+1} = \mathbf x_n + h\,\mathbf v_{1/2},\quad
\mathbf v_{n+1} = \mathbf v_{1/2} + \tfrac h2 \mathbf a(\mathbf x_{n+1}).
$$

Second order, symmetric (time reversible) and symplectic, with one force evaluation per step since the final
acceleration is the next step's first. It conserves linear and angular momentum exactly for pairwise central forces.
With rods it becomes RATTLE: Lagrange multipliers project positions back onto the constraints after the drift and
velocities onto their tangent space after the final kick, keeping the scheme symplectic on the constraint manifold.

### Compositions: Yoshida, PEFRL, Blanes–Moan

A symmetric second-order step $S(h)$ composed as $S(w_1 h) S(w_0 h) S(w_1 h)$ with
$w_1 = 1/(2 - 2^{1/3})$, $w_0 = 1 - 2w_1$ cancels the $h^3$ error term: **Yoshida's triple jump**, fourth order.
Repeating the trick gives orders 6 and 8 (7 and 15 substeps). The negative middle weight means a backward substep,
which is why the error constant is large. Optimised splittings place kicks and drifts with free coefficients chosen to
minimise the leading error: **PEFRL** (Omelyan, Mryglod & Folk 2002) and the 6-stage scheme of **Blanes & Moan**
(2002), whose errors are about 60× and 190× smaller than Yoshida's at equal cost. `w.use_composition(weights, order)`
and `w.use_splitting([("kick", b1), ("drift", a1), ...], order)` build your own.

### Gauss–Legendre (implicit Runge–Kutta)

The $s$-stage Gauss collocation method has order $2s$ and is symplectic for *every* Hamiltonian, including
velocity-dependent ones, and conserves all quadratic invariants (angular momentum, the energy of linear systems)
exactly. Its stage equations are implicit and solved by fixed-point iteration, so a step costs several force
evaluations and needs $h$ small enough for the iteration to contract. `gauss2` is the implicit midpoint rule.

### Wisdom–Holman

For planets around a star, $H = H_\text{Kepler} + H_\text{interaction}$ with the second part smaller by the
planet/star mass ratio $\mu \sim 10^{-3}$. Wisdom–Holman alternates exact Kepler drifts (solved with universal
variables in democratic heliocentric coordinates, Duncan, Levison & Lee 1998) with interaction kicks and a small
recoil ("jump") term, so its error is $O(\mu h^2)$ rather than
$O(h^2)$: the outer solar system keeps $|\Delta E/E| \approx 4\times10^{-6}$ over a million years at 1/20 of Jupiter's
period, against $2\times10^{-3}$ for Verlet. Particle 0 must be the central body and gravity must be a single
`NewtonianGravity`.

### Boris

For $\mathbf a = (q/m)(\mathbf E + \mathbf v \times \mathbf B)$: half an electric kick, an exact rotation of
$\mathbf v$ about $\mathbf B$ (by $2\arctan(\omega h/2)$ with $\omega = qB/m$, a relative phase error of
$(\omega h)^2/12$), another half kick, all between two half drifts. It is second order and volume preserving, and in a pure magnetic field it conserves $|\mathbf v|$ to round-off for any step size,
so gyration never spirals. Other forces must not depend on velocity.

### Runge–Kutta: RK4 and Dormand–Prince

Classical RK4 (4 evaluations, order 4) and Dormand–Prince 5(4) (6 evaluations with FSAL reuse, order 5) treat any
$\mathbf a(t, \mathbf x, \mathbf v)$ the same way, which makes them the safe choice with drag or other dissipation.
They are not symplectic: on a conservative system the energy drifts linearly in time.

**Adaptive stepping** (`w.run_adaptive`) uses Dormand–Prince's embedded 4th-order solution as an error estimate.
A step is accepted when $\|\mathbf e\|_\text{rms} \le 1$ with $e_i = \text{err}_i / (\text{atol} + \text{rtol}\,|y_i|)$,
and the next step size comes from Hairer's PI controller: $h_\text{new} = 0.9\,h\,
\|\mathbf e_n\|^{-0.17}\,\|\mathbf e_{n-1}\|^{0.04}$, kept within $[0.2h, 10h]$ (and at most $h$ after a rejection). Outputs at requested times and event
locations use the method's 4th-order continuous extension. Adaptive steps win when time scales vary (close
encounters, eccentric orbits: 350 to 16 000 times fewer evaluations at $e = 0.99$); for very long smooth runs a
fixed-step symplectic method's bounded error wins.

### Thermostats

`use_langevin(T, friction, seed)` integrates the Langevin equation
$\dot{\mathbf v} = \mathbf a - \gamma\mathbf v + \sqrt{2\gamma k_BT/m}\,\boldsymbol\xi$ with the BAOAB splitting
(Leimkuhler & Matthews), which samples configurations with small bias even at large steps; the noise is counter-based,
so runs are reproducible and restart exactly. `use_nose_hoover(T, tau)` adds a thermostat variable $\xi$ with
$\dot\xi = (2K - g k_BT)/Q$; the extended energy $E + $ `thermostat_energy()` is conserved. A single Nosé–Hoover
thermostat is not ergodic for small or stiff systems (its temperature fluctuations come out too large there); prefer
Langevin when canonical fluctuations matter.

## Choosing

- Conservative, long, smooth: `verlet` (cheap) or `yoshida4`/`pefrl`/`blanes_moan4` (accurate).
- Planets around a star: `wisdom_holman`.
- Close encounters, eccentric orbits, transients: `run_adaptive`.
- Drag or other dissipation: `rk4` or `dopri5`.
- Magnetic fields: `boris`; general velocity-dependent Hamiltonians: `gauss4`.
- Constant temperature: `use_langevin`.
- Rods: `verlet` or `yoshida4` (RATTLE); rigid bodies: `RigidSystem` (its own splitting, below in
  {doc}`methods`).
