# Rutherford scattering, closest approach and Coulomb trajectories

Backlog item 37 · module `physim.nuclear.rutherford` · tests `tests/python/test_nuclear_rutherford.py` · theory
{doc}`../theory/rutherford`

## Model

- Rutherford's cross section in the CM frame from the CM kinetic energy. Lab cross sections via the kinematic
  Jacobian, for the ejectile or the recoil, on both branches where double-valued.
- The integrated cross section over an angular range.
- Closest approach and impact parameter.
- Validity checks: interaction radius and grazing angle, Coulomb barrier, Sommerfeld parameter, electron screening,
  identical particles.
- Classical orbits integrated by physim's engine.

## Assumptions and range of validity

- Point charges, classical orbits, no nuclear force: valid below the grazing angle. The default interaction radius
  R = 1.2 (A₁^⅓ + A₂^⅓) + 2 fm is a convention; `r0` and `margin` can be changed.
- The cross section is unscreened. The screening correction is reported as a warning, not applied.
- Mott scattering of identical particles is not implemented (warned).

## Validation

| Check | Reference | Tolerance | Result | Test |
|---|---|---|---|---|
| CM cross section; 180° value (d₀/4)² | Analytic formula | 1e-12 | pass | `test_cross_section_formula` |
| Lab cross section, 4 angles, ¹⁶O + ²⁰⁸Pb | Independent count of solid angle (numerical cos θ differences) | 1e-6 | pass | `test_lab_cross_section_is_the_cm_one_times_the_jacobian` |
| Recoil at 0° lab = 4 σ_cm(180°) | Kinematics | 1e-3 | pass (also matches LISE++'s plot, ~340 mb/sr) | same |
| Integrated cross section | Numerical integration | 1e-8 | pass | `test_integrated_cross_section` |
| Orbit geometry: b(θ), r_min(θ) = r_min(b) | Analytic | 1e-12 | pass | `test_orbit_geometry` |
| Deflection of integrated orbits, 3 systems × 5 impact parameters | Analytic 2 arctan(d₀/2b) | 1e-6° | pass (worst ≈ 2e-7°) | `test_integrated_orbits_match_the_analytic_deflection` |
| ⁴He + ¹⁹⁷Au, 20 MeV, 50° CM: σ_cm, σ_lab, grazing angle | LISE++ kinematic calculator (user's screenshot) | 0.1% | pass (0.03%) | `test_lise_value_for_alpha_on_gold` |
| Angular distribution 15°–150°, gold | Geiger and Marsden, Phil. Mag. **25**, 604 (1913), Table II | data/theory within ±25% of its mean; spread < 1.5 (a 1/sin² law: > 16) | pass (worst 23%, at 22.5°; the paper's own column scatters ±18%) | `test_geiger_marsden_angular_distribution` |
| Departure from Rutherford above the barrier where the warning says | Farwell and Wegner, Phys. Rev. **95**, 1212 (1954) | — | 🟡 data requested (`tests/reference/nuclear/pending/rutherford_above_barrier.md`); physim predicts the onset for α + Pb at 60° at ≈ 32.5 MeV | — |

The LISE++ value (2640.611 mb/sr) lies between physim's values at the full beam energy (−0.03%) and at the
mid-target energy LISE++ quoted (+0.07%). No combination of energy, masses or e² reproduces it exactly, so the
difference is in LISE++'s internal energy convention and is far below any measurable effect: screening alone changes
this case by 0.2%.

Validated 2026-10-04.

## Known limitations

- No optical-model (nuclear) elastic scattering above the grazing angle; the planner warns there instead.
- No Mott scattering for identical particles, and no screening correction applied to the numbers.

## References

- E. Rutherford, Phil. Mag. **21**, 669 (1911); H. Geiger and E. Marsden, Phil. Mag. **25**, 604 (1913).
- J. L'Ecuyer, J. A. Davies and N. Matsunami, Nucl. Instr. Meth. **160**, 337 (1979): electron screening.
- G. W. Farwell and H. E. Wegner, Phys. Rev. **95**, 1212 (1954).
