# Quickstart

## A world, a force, a run

```python
import numpy as np
import physim as ps

w = ps.World(integrator="yoshida4")            # pick any name from ps.INTEGRATORS
w.add_particle([0.5, 0, 0], vel=[0, 0.4, 0], mass=0.5)
w.add_particle([-0.5, 0, 0], vel=[0, -0.4, 0], mass=0.5)
w.add_force(ps.NewtonianGravity(G=1.0))

traj = w.run(dt=1e-3, steps=10_000, record_every=10)
traj.t.shape, traj.pos.shape                   # (1001,), (1001, 2, 3): frames x particles x xyz
ps.relative_energy_error(traj).max()           # bounded, ~1e-12 for this symplectic integrator
```

A `World` holds the state (positions, velocities, masses, charges, radii, pinned flags), a list of forces and an
integrator. `run` takes fixed steps and records every `record_every`-th one; the world keeps its final state, so
calling `run` again continues the simulation.

## Ready-made systems and units

```python
w = ps.scenarios.solar_system(["Jupiter", "Saturn"])   # AU, solar masses, years (G = 4π²)
traj = w.run(0.5, 2000, record_every=10)                # Wisdom-Holman, 1000 years
ps.units.ASTRONOMICAL.to_si(traj.vel[-1, 1], "velocity")  # Jupiter's velocity in m/s
```

## Finding moments exactly

```python
periapsis = ps.Event.radial_velocity(1, 0, direction=+1)
traj = w.run(0.01, 100_000, record_every=100, events=[periapsis])
traj.event_t                                            # located to ~1e-12 of a step
```

## Plotting

```python
ps.plot.orbits(traj, labels=["Sun", "Jupiter", "Saturn"])
ps.plot.energy_error(traj)
anim = ps.plot.animate(traj, trail=30); anim.save("orbits.gif", writer="pillow")
```

## Your own force, in Python first

```python
def accel(t, pos, vel, mass):                  # (N, 3) arrays in, (N, 3) accelerations out
    return -pos * np.linalg.norm(pos, axis=1, keepdims=True)

w.add_force(ps.CustomForce(accel, velocity_dependent=False))
```

When it works, port it to Rust for speed: see {doc}`new-ideas`.

## Where next

- {doc}`guide`: everything the engine does, feature by feature.
- {doc}`theory/integrators`: which integrator to choose and why.
- The example notebooks: integrator comparisons, an N-body cluster, chaos in the three-body problem, a driven
  oscillator's route to chaos, charged particles in a magnetic bottle, fields on grids, particle accelerators
  (linac, cyclotrons, Rutherford scattering), and quantum wave packets, tunnelling and bound states.
