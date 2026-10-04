# 35 · Two-body reaction kinematics

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner · **Status: Done** (LISE++ comparison pending)

> **Done** (`python/physim/nuclear/kinematics.py`, `tests/python/test_nuclear_kinematics.py`,
> `docs/theory/kinematics.md`, `docs/physics-register/kinematics.md`, `tests/reference/nuclear/`).
> - `TwoBody(beam, target, beam_energy, ejectile=None, recoil=None, excitation_mev=0, excite="recoil")`: exact
>   relativistic kinematics with AME2020 nuclear masses; `at_cm`, `theta_cm`/`at_lab` (both branches),
>   `recoil_for`, `max_angle`, `double_valued`, `threshold_mev`, `q_value_mev`; NumPy-vectorised.
>   `kinematic_factor` (non-relativistic K) and `elastic(experiment)` for each target isotope.
> - CM momentum written as Q(Σm) + 2m₂T so it stays exact at low energy (the direct form lost about 6 digits at 1 keV);
>   g = 1 (elastic recoils, equal masses) handled exactly, giving the 90° maximum.
>
> Results (38 tests): lab values agree with an explicit four-vector boost to 1e-9° / 1e-10 for seven reactions;
> energy and momentum conserved to 1e-9; the non-relativistic limit matches K to 1e-6 (4 mass ratios); Jacobian and
> dE/dθ match numerical derivatives to 1e-6 / 1e-5; lab → CM inverts CM → lab to 1e-7° on both branches;
> ³H(p, n)³He threshold 1.0190 MeV (tabulated 1.019).
>
> Not done here: the **LISE++ comparison** needs LISE++ runs by hand. The eight cases (inverse kinematics and two
> double-valued cases among them), with physim's values, are in `tests/reference/nuclear/pending/kinematics_lisepp.md`;
> the test that reads the results is written and skipped until they exist. The **Rust port** moves to item 39,
> where the event generator needs it.

## Why
Every output of the planner starts from kinematics: which energy a scattered particle has at a given angle, where
the recoil goes, how lab and centre-of-mass quantities convert. Writing it for a general reaction a(b, c)d with a
Q-value now means inelastic scattering (item 43) and transfer reactions only add a term, not a rewrite.

## Scope
- Fully relativistic two-body kinematics for a(b, c)d with Q = (m_a + m_b − m_c − m_d)c² − E*; elastic scattering is
  the case c = b, d = a, E* = 0.
- Lab energy and angle of the ejectile and the recoil, as functions of lab or CM angle.
- Both solutions where the kinematics are double-valued (beam heavier than the target, or endothermic reactions near
  threshold), and the maximum lab angle.
- Lab ↔ CM angle conversion and the solid-angle Jacobian dΩ_cm/dΩ_lab (needed by item 37).
- Kinematic broadening dE/dθ, used to estimate peak widths from a detector's angular size.
- Threshold energy for reactions with Q < 0.
- Prototype in Python, then port to Rust (the `docs/new-ideas.md` workflow) once validated, because item 39 calls it
  per event.

## Depends on
34.

## Done when
- The non-relativistic limit matches the analytic elastic formulas (kinematic factor K) to 1e-6 at low energy.
- Energies and angles match LISE++ (or CATKIN) for at least five cases, including inverse kinematics and a
  double-valued case, within 1 keV and 0.01°. Reference files stored as described in item 40.
- Jacobian checked against a numerical derivative of the angle conversion.
- Theory note in `docs/theory/`, and an entry in the physics register.
