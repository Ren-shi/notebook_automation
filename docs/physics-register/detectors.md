# Detector geometry and response

Backlog item 38 · modules `physim.nuclear.detectors`, `physim.nuclear.plot` · tests
`tests/python/test_nuclear_detectors.py` · theory {doc}`../theory/detectors`

## Model

- Flat rectangular (double-sided strips), disc and annular (rings × sectors) faces, placed and oriented as in the
  setup file.
- Ray hits with segment indices and incidence angles.
- Solid angles by Gauss–Legendre quadrature over each face or segment.
- Angular coverage, with exact edge extremes for whole faces.
- Shadowing and beam-blocking warnings.
- Response: dead layer, punch-through, Gaussian resolution, threshold.
- The exit path from a reaction inside the target.
- Pictures: the setup in 3D, and detector coverage in (θ, φ).

## Assumptions and range of validity

- Point source at the target unless a source point is given. The event generator (item 39, {doc}`rates-and-events`)
  samples the beam spot.
- Flat faces, no gaps between strips, no inter-strip charge sharing. Resolution is constant in energy.
- Pulse-height defects in silicon for heavy ions, and ΔE–E telescopes, are not modelled.

## Validation

| Check | Reference | Tolerance | Result | Test |
|---|---|---|---|---|
| Disc on axis, 3 sizes | 2π(1 − cos α) | 1e-10 | pass | `test_disc_on_axis` |
| Centred rectangle facing the target, 3 sizes, off-axis direction | 4 arctan(ab / (d√(a² + b² + d²))) | 1e-10 | pass | `test_rectangle_facing_the_target` |
| Off-centre rectangle | Corner-rectangle decomposition (exact) | 1e-10 | pass | `test_off_centre_rectangle` |
| Annulus on axis, and each of its 384 segments | 2π(cos θ_in − cos θ_out) | 1e-9 | pass | `test_annulus_on_axis_and_its_segments` |
| Tilted and rotated faces (no closed form) | Independent Monte Carlo, 4 × 10⁵ directions | 4σ (≈ 0.3%) | pass | `test_monte_carlo_agrees_for_tilted_faces` |
| Strip indices, hit positions, θ and φ ranges | Construction; dense grid; atan(25/150) | 1e-9° | pass | `test_strips_add_up_and_are_hit_where_expected`, `test_angular_coverage` |
| Shadowing fraction of a disc wholly in front of another | Ω_front / Ω_back | 0.1% | pass | `test_shadowing_and_beam_warnings` |
| Punch-through energy | Range from {doc}`stopping` | 1e-6 | pass | `test_response_dead_layer_punch_through_threshold_resolution` |
| Resolution | σ = FWHM/2.355 from 40 000 samples | 2% | pass | same |
| Example setups drawn in 3D and in (θ, φ) | Checked by eye (2026-10-04); outlines match the computed θ ranges | 0.05° | pass | `test_setup_pictures` |

Validated 2026-10-04.

## Known limitations

- Segment θ ranges are taken over quadrature points, to a fraction of the segment's size. Only whole faces get the
  exact edge search.
- No inter-strip effects, dead strips or pulse-height defects.
