# Two-body reaction kinematics

Backlog item 35 · module `physim.nuclear.kinematics` · tests `tests/python/test_nuclear_kinematics.py` · theory
{doc}`../theory/kinematics`

## Model

Exact relativistic two-body kinematics for target(beam, ejectile)recoil with nuclear masses from AME2020 and an
optional excitation energy on either outgoing particle. It gives:

- lab angle and kinetic energy from the CM angle, and both CM solutions from a lab angle;
- maximum lab angle;
- threshold energy;
- solid-angle Jacobian dΩ_cm/dΩ_lab;
- kinematic broadening dE/dθ.

## Assumptions and range of validity

- Two bodies out. Breakup and three-body final states are not covered.
- The target is at rest, and the beam has a single energy. Energy loss in the target is applied separately
  (backlog item 36).
- Valid at any energy: there is no non-relativistic approximation.

## Validation

| Check | Reference | Tolerance | Result | Test |
|---|---|---|---|---|
| Lab angle and energy for 7 reactions (elastic, inelastic, inverse kinematics, (d,p), near threshold, 167 MeV/u) | Independent explicit Lorentz boost of the CM four-momenta | 1e-9° and 1e-10 relative | pass | `test_matches_explicit_lorentz_boost` |
| Energy and momentum conservation | Exact | 1e-9 relative | pass | `test_energy_and_momentum_are_conserved` |
| Non-relativistic limit, 4 mass ratios incl. m₁ > m₂ | Analytic kinematic factor K | 1e-6 | pass | `test_non_relativistic_limit_is_the_kinematic_factor` |
| Jacobian and dE/dθ, ejectile and recoil | Numerical derivatives; non-relativistic Jacobian formula | 1e-6 / 1e-5 relative | pass | `test_jacobian_and_broadening_match_numerical_derivatives`, `test_jacobian_non_relativistic_formula` |
| Both branches, maximum angle | Round trip CM → lab → CM; brute-force maximum; arcsin(m₂/m₁) limit | 1e-7° / 1e-6° / 0.2% | pass | `test_lab_to_cm_inverts_cm_to_lab`, `test_double_valued_cases` |
| ³H(p, n)³He threshold | Tabulated 1.019 MeV | 0.5 keV | 1.0190 MeV | `test_q_value_and_threshold` |
| ⁴He + ¹⁹⁷Au at 20 MeV, 50° CM: ejectile and recoil | LISE++ kinematic calculator (screenshots) | 1 keV or the digits shown (2 keV) | pass: every displayed digit agrees | `test_against_lise_reference_files` |
| 8 more cases incl. inverse kinematics and two double-valued cases | LISE++ | 1 keV, 0.01° | 🟡 request written (`tests/reference/nuclear/pending/kinematics_lisepp.md`) | same test, once the results are added |

Validated 2026-10-04.

## Known limitations

- Python (NumPy, vectorised) only for now; the Rust port comes with the event generator (backlog item 39), which
  calls it per event.
- An ejectile at rest in the CM frame (exactly at threshold) has no direction; the code reports its lab angle as the
  beam direction.

## References

- R. Hagedorn, *Relativistic Kinematics* (Benjamin, 1963), ch. 7.
- E. Byckling and K. Kajantie, *Particle Kinematics* (Wiley, 1973), ch. 4.
