# physim

A classical-mechanics engine: the time stepping and forces run in **Rust**, and **Python** sets up systems and
analyses results. Point particles under gravity, springs, fields and pair potentials; constraints, collisions and
events; rigid bodies; symplectic, adaptive and stochastic integrators; chaos indicators; and fields on grids
(waves, diffusion, Poisson's equation, particle-mesh gravity).

```python
import physim as ps

w = ps.scenarios.two_body(m1=1.0, m2=1e-3, a=1.0, e=0.5)
traj = w.run(dt=1e-3, steps=10_000, record_every=10)
ps.plot.orbits(traj)
ps.relative_energy_error(traj).max()
```

```{toctree}
:maxdepth: 2
:caption: Using physim

install
quickstart
guide
new-ideas
```

```{toctree}
:maxdepth: 2
:caption: Theory

theory/integrators
theory/forces
theory/methods
theory/fields
```

```{toctree}
:maxdepth: 1
:caption: Examples

examples/01_getting_started
examples/02_nbody_cluster
examples/03_three_body_chaos
examples/04_driven_oscillator
examples/05_charged_particles
examples/06_fields
examples/07_accelerators
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/python
api/rust
changelog
```
