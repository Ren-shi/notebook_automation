# Request: Coulomb-excitation cross sections from GOSIA (or CLX), and a published measurement

For backlog item 43. physim's first-order semiclassical Coulomb excitation (`physim.nuclear.coulex`) passes its
closed-form and engine-orbit checks. It still needs comparing with the standard code and with data.

## 1. GOSIA (or Winther's CLX), first order

Use the setup of the example `coulex_ni58`:
- beam ¹⁶O at 30 MeV on ⁵⁸Ni, a thin target (no energy loss);
- the 2⁺₁ state at 1.454 MeV, with B(E2; 0⁺ → 2⁺) = 0.0695 e²b²;
- no other states, no reorientation (quadrupole moment 0).

Give the **excitation probability** P(θ_CM) and dσ/dΩ (mb/sr, CM) at CM angles 60°, 90°, 120°, 150° and 170°. In
GOSIA, OP,POIN with the beam energy fixed and the matrix element ⟨2⁺‖E2‖0⁺⟩ = √(B(E2)) = 0.2636 eb. Note the GOSIA
version and whether quantum corrections are on.

Save as `tests/reference/nuclear/coulex_gosia.csv` with the usual header (tool, produced by, settings) and the
columns `beam,target,beam_energy_mev,state_mev,b_e2_e2b2,theta_cm_deg,probability,dsigma_cm_mb_sr`.

physim's values for comparison:

```python
from physim.nuclear.coulex import Coulex
c = Coulex("16O", "58Ni", 30.0, energy=1.454, b_up="0.0695 e2b2")
[(t, c.probability(t, exact=True), c.cross_section_cm(t)) for t in (60, 90, 120, 150, 170)]
```

## 2. A published Coulomb-excitation measurement

Any first-order case with a quoted cross section or excitation probability against angle, with its beam, energy,
target, state and B(E2). Classic choices are ⁵⁸,⁶⁰Ni or ¹¹⁴Cd with α particles or ¹⁶O. Transcribe the table into
`tests/reference/nuclear/coulex_literature.csv` with the same columns plus `uncertainty`, and record the reference
(authors, journal, table) in the header.
