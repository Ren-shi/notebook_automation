# Quantum mechanics

`physim.fields::quantum` (Rust) and the Python class `Schrodinger` solve the time-dependent Schrödinger equation for
one particle on a 1D, 2D or 3D grid, find bound states, and measure the usual observables. `physim.quantum` adds
helpers (momentum distributions, probability currents, analytic transmission coefficients and energy levels) and
`physim.plot.wavefunction` / `animate_wavefunction` draw the results.

## The equation

$$
i\hbar\,\frac{\partial\psi}{\partial t} = H\psi, \qquad
H = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf x, t) - iW(\mathbf x),
$$

with $\hbar$ and $m$ as parameters (both 1 by default: atomic-style units). $V$ is a static array or a function of
time, and $W \ge 0$ is an optional *absorbing potential* that removes outgoing waves near the edges so they do not
reflect or wrap around. The wavefunction lives on the same grids as the other field solvers (see {doc}`fields`),
normalised so that $\sum|\psi|^2 h^d = 1$; by default the grid is centred on zero.

Why only one particle: the state of $N$ particles is a function on a $3N$-dimensional space, so the grid grows as
$n^{3N}$. One particle in 3D at $64^3$ is 262 144 complex numbers; two particles would need $64^6 \approx 7\times10^{10}$.
Many-body quantum mechanics needs different methods (Hartree–Fock, density functional theory, tensor networks,
quantum Monte Carlo) and is out of scope.

## Split-step Fourier

On periodic grids whose sizes are powers of two (the default there), `method="split_step"` advances
$\psi \leftarrow e^{-iH\tau/\hbar}\psi$ with the Strang splitting

$$
e^{-iH\tau/\hbar} \approx e^{-iV\tau/2\hbar}\; e^{-iT\tau/\hbar}\; e^{-iV\tau/2\hbar},
$$

applying the potential factor pointwise and the kinetic factor $e^{-i\hbar k^2\tau/2m}$ exactly in Fourier space (two
FFTs per step). Space is therefore resolved *spectrally*: a free packet disperses exactly, and smooth bound states
converge exponentially with the grid spacing. Every factor has modulus one, so the norm is conserved to round-off
and the energy error is bounded ($O(\tau^2)$) with no drift. With `order=4`, three Strang steps of $w_1\tau$,
$w_0\tau$, $w_1\tau$ (Yoshida's weights, $w_1 = 1/(2 - 2^{1/3})$, $w_0 = 1 - 2w_1$) give fourth order at three times
the cost. A time-dependent potential is evaluated at the ends of each half kick, which keeps the composition
symmetric and the orders intact.

**Step size with sharp potentials.** For smooth potentials split-step is accurate at any step the dynamics allow. A
discontinuous potential (a step or barrier) kicks the wave into the shortest wavelengths the grid holds,
$k_\text{max} = \pi/h$, and their phases $\hbar k_\text{max}^2\tau/2m$ must be resolved: keep
$\tau \lesssim 2mh^2/\pi\hbar$, or grid-scale ripple appears (it is a phase error, not an instability).

## Crank–Nicolson

On any grid (and the default for Dirichlet and Neumann boundaries), `method="crank_nicolson"` solves

$$
\Big(1 + \frac{i\tau}{2\hbar}H\Big)\psi^{n+1} = \Big(1 - \frac{i\tau}{2\hbar}H\Big)\psi^n
$$

with the second-order Laplacian and $H$ at the midpoint time. This is the Cayley transform of $H$: unitary, so the
norm is conserved, and it commutes with $H$, so for a static potential the *discrete* energy is conserved exactly
too. It is second order in space and time and stable for any step. On Dirichlet grids $\psi$ vanishes on the edge
nodes (infinite walls); Neumann walls reflect with zero slope; periodic grids may have any size.

In 1D the system is tridiagonal and is solved directly (the Thomas algorithm, with a Sherman–Morrison correction on
periodic grids). In 2D and 3D it is solved by BiCGSTAB with a Jacobi preconditioner, to a relative residual of
$10^{-12}$, starting from the previous $\psi$; a few iterations per step are typical. The finite-difference
Laplacian makes waves slightly too slow (group velocity $\sin(kh)/h$ instead of $k$), an $O(h^2)$ error that shows
up as a lag of fast packets.

## Bound states

`eigenstates(count)` returns the lowest eigenvalues and (real, normalised) eigenfunctions of the *same* discrete
$H$ the time stepping uses: spectral for split-step, finite differences for Crank–Nicolson. So an eigenstate placed in
$\psi$ is stationary, picking up only the phase $e^{-iEt/\hbar}$. The solver is LOBPCG (locally optimal block
preconditioned conjugate gradients) on a block of `count` vectors plus guard vectors, shifted by
$\sigma = V_\text{min} - \hbar^2/(2mL^2)$ so that $H - \sigma$ is positive definite. Preconditioners: $(\hbar^2k^2/2m +
c)^{-1}$ in Fourier space for the spectral operator ($c$ the mean of $V - \sigma$); an exact tridiagonal solve of
$H - \sigma$ in 1D with finite differences; a few Jacobi-preconditioned conjugate-gradient iterations in 2D and 3D. It
stops when every residual $\lVert H\psi - E\psi\rVert$ is below `tol` $(E - \sigma)$; eigenvalue errors are then of
order `tol`², far below the discretisation error.

## Observables

All expectation values are normalised by the current norm (which an absorber reduces):
$\langle x\rangle$ and $\langle x^2\rangle - \langle x\rangle^2$ from $|\psi|^2$; $\langle p\rangle$ and
$\langle T\rangle$ spectrally for split-step (from $|\hat\psi(k)|^2$) and by central differences otherwise;
$\langle V\rangle$ at the current time; `energy()` $= \langle T\rangle + \langle V\rangle$ (the absorber is not part
of it). `physim.quantum.momentum_density` returns $|\phi(p)|^2$ on the momentum grid $p = \hbar k$,
`probability_current` the current $\mathbf j = (\hbar/m)\,\mathrm{Im}(\psi^*\nabla\psi)$, and
`probability(s, region)` the probability inside a region.

## Validation

Measured by `tests/quantum.rs` and `tests/python/test_quantum.py` ($\hbar = m = 1$):

| Check | Result |
| --- | --- |
| Norm over $10^5$ steps (1D, both methods) | conserved to $< 10^{-11}$ |
| Energy over $10^5$ steps, static potential | Crank–Nicolson $< 10^{-10}$; split-step $< 10^{-6}$, no drift |
| Free Gaussian packet, centre and width $\sigma(t) = \sigma_0\sqrt{1 + (\hbar t/2m\sigma_0^2)^2}$ | split-step to $10^{-9}$; Crank–Nicolson second order in $h$ |
| Coherent state in a harmonic trap | $\langle x\rangle = x_0\cos\omega t$ and constant width to $10^{-7}$ (`order=4`) |
| Time-step convergence (anharmonic well) | Strang 2.0, Yoshida 4.0, Crank–Nicolson 2.0 |
| Oscillator levels $(n + \tfrac12)\hbar\omega$, spectral | to $10^{-9}$ (1D, 5 levels); 2D degeneracies 1, 2, 2, 3, 3, 3 |
| Oscillator and infinite-well levels, finite differences | converge at order 2.0 in $h$ |
| Rectangular-barrier transmission (packet-averaged analytic $T$) | within 0.2% above and below the barrier |
| Driven trap $V = \tfrac12(x - a\sin\Omega t)^2$ | $\langle x\rangle$ matches the classical solution |
| Ehrenfest: double well, small $\hbar$ | $\langle x\rangle$ follows a classical `World` run to 0.03 for $t \le 1$, then departs |
| Absorbing layer | an outgoing packet leaves $< 10^{-3}$ of its probability |

## Cost

On the development machine (12 logical cores; the FFTs run in parallel):

| Case | Time per step |
| --- | --- |
| Split-step, 1D, 4096 points | 0.6 ms |
| Split-step, 2D, $256^2$ | 10 ms (29 ms with `order=4`) |
| Split-step, 3D, $64^3$ | 16 ms |
| Crank–Nicolson, 1D, 4096 points | 0.17 ms |
| Crank–Nicolson, 2D, $256^2$ | 20–35 ms (2–4 BiCGSTAB iterations) |

Ten oscillator eigenstates take 0.3 s on a 1024-point spectral grid and 0.13 s with 2001-point finite differences;
2D and 3D problems take seconds. Python potential functions are called once or twice per step (vectorised), which
is negligible next to the FFTs for grids of more than a few thousand points.
