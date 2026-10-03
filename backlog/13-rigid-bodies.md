# 13 · Rigid bodies

**Priority:** P3 · **Size:** L · **Area:** Physics

## Why
Spinning tops, tumbling bodies and the intermediate-axis (Dzhanibekov) instability need rotational
degrees of freedom, not just point particles.

## Scope
- Bodies with orientation (unit quaternion), angular momentum, and an inertia tensor.
- Torques from forces applied at body points.
- Symplectic free-rigid-body integrator (e.g. splitting into exactly solvable rotations).

## Notes
- Alternative route: rigid bodies as particle clusters held by constraints (item 07). Simpler, but slower and less accurate.

## Done when
- A torque-free asymmetric top shows the intermediate-axis flip and conserves energy and |L| to tolerance.
- A heavy symmetric top precesses at the rate predicted in the fast-top limit.
