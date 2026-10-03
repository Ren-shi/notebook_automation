# 13 · Rigid bodies

**Priority:** P3 · **Size:** L · **Area:** Physics · **Status: Done**

> **Done** (`src/rigid.rs`, `src/python/rigid.rs`; Python `RigidSystem`, `RigidTrajectory`, `BodyGravity`,
> `BodySpring`, `ps.inertia_*`, `ps.quaternion_from_axis_angle`).
> - A separate `RigidSystem` (not part of `World`): bodies with mass, principal moments, position, velocity, unit
>   quaternion and space-frame angular momentum. Free bodies are referenced to their centre of mass; pivoted bodies
>   turn about a fixed point with the centre of mass on a principal axis (parallel-axis theorem applied for you).
> - Forces are `BodyForce`s returning a force and a torque per body: `BodyGravity`, `BodySpring` between body points
>   or fixed anchors, closures (Rust) or callables (Python).
> - Step: kick(h/2) – drift(h) – kick(h/2) with the free rotation split as R1(h/2) R2(h/2) R3(h) R2(h/2) R1(h/2), each
>   an exact rotation about a body axis. Symplectic, time reversible, second order; `order=4` is the Yoshida
>   composition. The angular momentum of a torque-free body never changes.
>
> Results (`tests/rigid.rs`, `tests/python/test_rigid.py`, notebook §14, `cargo run --release --example spinning_tops`):
> - Torque-free box I = (1, 2, 3) spun about axis 2: flips every 16.60 time units, matching the elliptic-integral
>   period (Landau & Lifshitz §37) to 1e-4; energy error 7e-9 over 2e5 steps at dt = 1e-3; L exact.
> - Heavy symmetric top (spin 50, tilt 0.5): precession 0.04009 vs the fast-top rate mgl/(I3 ω3) = 0.04 (+0.23%, the
>   size of the expected next-order correction); energy to 1e-6 and L_z to 1e-10.
> - Torques equal minus the angular gradient of the potential; two spring-coupled tumbling bodies conserve linear
>   and total angular momentum to round-off; observed convergence orders 2 and 4.
>
> Left out: contact and collisions between bodies, rigid bodies inside `World` (mixing with particles), checkpoints
> for `RigidSystem`, general (non-principal) pivot offsets.

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
