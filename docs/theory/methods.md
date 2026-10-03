# Numerical methods

How the engine does the things that are not a plain force evaluation: events, constraints, collisions, chaos
indicators, rigid bodies, parallelism and restarts.

## Events

An event is a function $g(t, \mathbf x, \mathbf v)$ whose zeros matter (a periapsis, a crossing, a collision). After
every step the engine compares the sign of $g$ at the two ends. On a sign change in the requested `direction` it finds
the root by **re-stepping the integrator from the start of the step** with a fraction $\theta h$, using the Illinois
variant of regula falsi on $\theta$, to about $10^{-12}$ of a step. The located state therefore lies on the
integrator's own trajectory and is as accurate as the integrator. A terminal event leaves the world exactly at the
event. Two zeros of the same event inside one step are invisible, so the step must be short compared with the spacing
of events. With adaptive stepping the root is found on the continuous extension instead, at no extra force cost.

## Constraints: RATTLE

A rod is a holonomic constraint $\sigma(\mathbf x) = |\mathbf x_i - \mathbf x_j|^2 - L^2 = 0$. RATTLE (Andersen 1983)
adds forces $-G^T\boldsymbol\lambda$ with $G = \partial\sigma/\partial\mathbf x$ to velocity Verlet:

1. kick and drift with the multipliers $\boldsymbol\lambda$ chosen so that every rod has its length after the drift
   (SHAKE; a nonlinear system);
2. after the final kick, project the velocities onto the constraints' tangent space, $G\,\mathbf v = 0$ (a linear
   system, $G M^{-1} G^T \boldsymbol\mu = G\mathbf v$).

Both systems couple all rods that share particles, and all rods are solved together (M-SHAKE): each outer iteration
solves the linearised system $G_0 M^{-1} G_0^T \boldsymbol\lambda = \text{residual}$ by conjugate gradients, which
converges on long chains and closed loops where rod-by-rod iteration is slow, down to the relative
`constraint_tolerance` (default $10^{-10}$). The result is
symplectic on the constraint manifold and second order; `yoshida4` composes it to fourth order. The multipliers
give the tension in each rod.

## Hard collisions

Spheres with radii collide instantaneously. After each step the engine looks for the **earliest contact** in it:
candidate pairs come from a uniform grid over swept spheres (radius plus the distance moved), and contact times from
the straight path between the start and end positions, which is exact without forces. With forces, the contact time is
refined on the integrator's trajectory by Illinois root finding on the gap. The world is re-stepped to the contact,
the impulse $J = -(1 + e)\,u_n/(m_i^{-1} + m_j^{-1})$ is applied along the line of centres (or the wall normal), and
the rest of the step is taken, until no contact remains. Pinned particles act as infinitely massive. With $e = 1$
energy is conserved to round-off (a 500-sphere gas keeps it to $3\times10^{-15}$).

## Lyapunov exponents and MEGNO

The variational equations $\dot{\delta\mathbf x} = \delta\mathbf v$,
$\dot{\delta\mathbf v} = \frac{\partial\mathbf a}{\partial\mathbf x}\delta\mathbf x + \frac{\partial\mathbf a}{\partial\mathbf v}\delta\mathbf v$
are integrated alongside the trajectory with the world's own integrator: tangent vectors ride along as extra
pseudo-particles whose acceleration is a Jacobian–vector product (analytic for most built-in forces, central finite
differences otherwise). **Benettin's method** re-orthonormalises the $n$ tangent vectors by QR decomposition at
intervals and averages the logarithms of the diagonal: the $n$ largest exponents. For a Hamiltonian system the
spectrum comes in $\pm$ pairs summing to zero. **MEGNO** (Cincotta & Simó 2000),
$Y(t) = \frac{2}{t}\int_0^t \frac{\dot\delta}{\delta}s\,ds$, has a time average $\langle Y\rangle \to 2$ for regular
(quasi-periodic) motion and $\langle Y\rangle \approx \lambda t/2$ for chaos, which separates the two much sooner than
the exponent itself, whose finite-time estimate decays only like $\ln t/t$ on regular orbits.

## Rigid bodies

`RigidSystem` stores each body's orientation as a unit quaternion $q$ and its angular momentum $\mathbf L$ in the space
frame. The Hamiltonian $H = \sum \frac{p^2}{2m} + \sum_k \frac{L_k^2}{2I_k} + V$ (with $L_k$ the body-frame components)
is split:

- **kick**: $\mathbf p \mathrel{+}= h\mathbf F$, $\mathbf L \mathrel{+}= h\boldsymbol\tau$ (torques about the reference
  point);
- **drift**: positions move with their velocities, and the free rotation is split into the three exactly solvable
  rotations about the body axes, $R_1(h/2)R_2(h/2)R_3(h)R_2(h/2)R_1(h/2)$ (McLachlan; Dullweber, Leimkuhler &
  McLachlan 1997). Under $L_k^2/2I_k$ alone the body turns about its own axis $k$ at the constant rate $L_k/I_k$, and
  $\mathbf L$ does not change.

Composed as kick–drift–kick it is second order, time reversible and symplectic, and a torque-free body's angular
momentum is exactly conserved; `order=4` applies Yoshida's composition. A pivoted body (a top) uses the moments of
inertia about the pivot (parallel-axis theorem) and its centre-of-mass offset for the gravity torque.

## Deterministic parallelism

Force loops are split into blocks whose boundaries depend only on the number of particles, not on the number of
threads, and the partial sums are combined in a fixed order. Floating-point addition is not associative, so this is
what makes results **bit-identical for any thread count**, and identical with the `parallel` feature off.

## Exact restarts

A checkpoint stores the state, every force's parameters, the integrator (including composition weights, thermostat
variables and the random-number counter), constraints and collision settings. Nothing else influences a step (forces
are pure functions; neighbour lists only cache, never change results; the random numbers are a pure function of
seed, step counter and particle index), so continuing from a checkpoint reproduces the uninterrupted run bit for bit.
