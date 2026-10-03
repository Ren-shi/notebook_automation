# Fields and continua

`physim.fields` (Rust) and the Python classes `WaveEquation`, `HeatEquation`, `solve_poisson` and `ParticleMesh`
solve partial differential equations on regular grids and couple them to particles.

## Grids and boundaries

A grid has 1 to 3 axes with $n_a$ points each, a common spacing $h$, and one boundary condition on every face:

- **Dirichlet**: the first and last point of each axis are boundary nodes whose values are held fixed (clamped
  strings and membranes, grounded conductors).
- **Neumann**: zero normal derivative. Points are cell centres $(i + \tfrac12)h$ and the boundary lies half a spacing
  outside the edge points; the ghost value equals its neighbour. This keeps the discrete Laplacian symmetric and
  conserves $\sum u\,h^d$ exactly in diffusion.
- **Periodic**: the grid wraps around.

The Laplacian is the second-order stencil $\nabla^2 u \approx \sum_a (u_{i+e_a} - 2u_i + u_{i-e_a})/h^2$.

## Wave equation

$u_{tt} = c^2\nabla^2 u$ is integrated with leapfrog (velocity Verlet on the semi-discrete system):
$v \mathrel{+}= \tfrac{\Delta t}2 c^2\nabla^2u$, $u \mathrel{+}= \Delta t\,v$, $v \mathrel{+}= \tfrac{\Delta t}2 c^2\nabla^2u$.
It is second order in space and time, symplectic (the discrete energy
$\tfrac12\sum v^2 + \tfrac12 c^2\sum_\text{edges}(\Delta u/h)^2$ has a bounded $O(\Delta t^2)$ error), and stable for
the CFL condition $c\,\Delta t/h \le 1/\sqrt d$; larger steps are refused.

## Heat equation

$u_t = D\nabla^2u$:

- **explicit** (forward Euler): $u \mathrel{+}= \Delta t\,D\nabla^2u$, first order in time, stable for
  $D\Delta t/h^2 \le 1/(2d)$ (refused otherwise). With $\Delta t \propto h^2$ the total error is $O(h^2)$.
- **Crank–Nicolson**: $(I - \tfrac{\Delta t}2 D\nabla^2)u^{n+1} = (I + \tfrac{\Delta t}2 D\nabla^2)u^n$, second order in
  time and unconditionally stable. The symmetric positive definite system is solved by conjugate gradients to a
  relative residual of $10^{-12}$ (with Dirichlet nodes eliminated, so the operator stays symmetric).

## Poisson's equation

$\nabla^2\varphi = f$:

- **periodic**: by FFT (sizes must be powers of two). Each Fourier mode is divided by the *discrete* Laplacian's
  eigenvalue $\sum_a (2\cos(2\pi k_a/n_a) - 2)/h^2$, so the result solves the finite-difference equations exactly
  (second order against the continuum). The mean of $f$ has no periodic solution and is removed; $\varphi$ has zero
  mean.
- **Dirichlet**: conjugate gradients on $-\nabla^2$ with the boundary values taken from the initial $\varphi$.
- **Neumann**: conjugate gradients on the singular but consistent system after removing the mean of $f$.

## Particle-mesh gravity

`ParticleMesh(box_size, cells, G, center, periodic, softening)` is a force that couples particles to a grid:

1. **Deposit**: each particle's mass is shared among the 8 surrounding grid nodes with cloud-in-cell (trilinear)
   weights.
2. **Solve**:
   - *isolated* (default): the mass grid is convolved with the softened Newtonian kernel
     $-G\,\mathbf r/(r^2 + \varepsilon^2)^{3/2}$ by FFT on a zero-padded grid of twice the size (Hockney & Eastwood),
     so there are no periodic images. The kernel transforms are computed once and cached; the zero padding is
     skipped in the transforms.
   - *periodic*: $\nabla^2\varphi = 4\pi G(\rho - \bar\rho)$ by FFT with the discrete Laplacian's eigenvalues, and
     $\mathbf a = -\nabla\varphi$ by central differences.
3. **Interpolate** the node accelerations back to the particles with the same weights.

Using identical weights for deposit and interpolation with an antisymmetric force kernel means a particle exerts no
force on itself and the total momentum is conserved to round-off. The cost is $O(N + M\log M)$ for $M$ grid points,
independent of how the particles are arranged.

**Accuracy.** Cloud-in-cell smoothing weakens forces at short range: the two-body force is 11% low at 2 cells, 4% at
4 cells, 1% at 8 cells and 0.1% at 32 cells. On a Plummer sphere of $10^5$ stars the median force error against direct
summation (with the same softening) is 4.5% on a $32^3$ grid and 1.1% on $64^3$, against 0.5% and 0.2% for
Barnes–Hut at $\theta = 0.5$; particle mesh is 9× and 1.6× faster than the tree there (`cargo run --release --example
particle_mesh`). Use it for large, smooth systems where the grid scale is an acceptable resolution limit; use
`TreeGravity` or direct summation when close encounters matter. In periodic mode, keep the grid no finer than the
mean particle spacing: a finer grid resolves the discreteness of the particles and aliases it into the forces.
