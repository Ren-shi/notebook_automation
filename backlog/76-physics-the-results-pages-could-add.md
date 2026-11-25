# 76 · Physics the results pages could add

**Priority:** P3 · **Size:** M (three parts, each S–M with its own validation) · **Area:** Nuclear planner (results)

## Why
The results ledger (`physim.nuclear.results`), the registry of observables (`physim.nuclear.observables`) and the
comparison engine (`physim.nuclear.compare`) put a value beside the evaluation, other measurements and theory, and
convert between observables. By the owner's decision of 2026-11-24 they exist without any formula of their own:
they call the conversions the code already had (B ↔ τ, B↑ ↔ B↓, Weisskopf units, the matrix element, Q_s from the
diagonal matrix element) and change none of them. Three things a results page would naturally show were left out on
purpose, because each needs a new formula and its own validation before a number from it stands beside a
measurement.

## Scope
- **The effective charge a shell-model value would need.** A shell-model B(E2) below the measured one is usually
  read as "which effective charges would match": the calculated proton and neutron parts of the matrix element,
  rescaled until B(E2) meets the measured value (one condition; with two charges, along the isoscalar or a fixed
  ratio). Needs the proton and neutron matrix elements separately from the shell-model result, the convention for
  the pair of charges, and the uncertainty carried from the measurement.
- **B(E2) expected from E(2⁺).** Grodzins/Raman-type systematics: the B(E2; 0⁺ → 2⁺₁) that the energy of the first
  2⁺ state suggests for a given Z and A, with the scatter of the systematics as its band. A third kind of
  expectation beside the evaluation and theory, useful most where nothing is measured. Needs the choice of the
  parametrisation, its region of validity (not near closed shells) and a check against the evaluated values.
- **Sum rules.** How much of a sum rule a set of transitions exhausts (the energy-weighted E2 sum rule first), as a
  check on a level scheme's B(E2) values or a shell-model calculation's. Needs the sum rule's form and constants
  written down and checked, and a statement of which levels the sum runs over.

Each becomes a route or a check of the registry (`observables.py`) once it is validated: a function with a
reference, a check in the validation suite, and a line in the physics register.

Also left out of the registry for the same reason (the formula exists but only inside another class or inline, and
the first waves did not refactor physics code): half-life ↔ mean life, and the branching of a level with several γ
rays (so B(E2) ↔ τ is converted for a level with one γ ray only). The rotor reading of a B(E2) (β₂, Q₀ and the rotor's
Q_s) is the Coulex analysis's `Analysis._shape`, called with a stand-in for the analysis
(`observables.rotor_reading`); making it a function of its own, with no number changed, would be cleaner.

Left out when the routes were tested for each nucleus and the checks between observables were added (the checks are
plain ratios of two rows):

- **The uncertainty of a ratio** (E(4⁺)/E(2⁺), B(E2; 4→2)/B(E2; 2→0), a calculation over the reference), with the
  correlation of two rows that share a source; of a matrix element from a B value (a square root); and of β₂, Q₀ and
  the rotor's Q_s. Without them a check has no verdict.
- **Comparing the rotor's Q_s with the measured one as a pull**, which needs that uncertainty.
- **B(M1) from a lifetime**: a conversion τ → B(M1) in the comparison engine and a B(M1) row of the lifetime analysis
  (it works the B(M1) of a pure M1 line out but writes τ alone), with the mixing ratio; and the E2/M1 mixing ratio from
  B(E2) and B(M1).
- **What a magnetic moment is read against**: the g factor from μ and J, the Schmidt values, Z/A.
- **Whether a target or a beam can be had.** The code knows the natural abundance and ENSDF's half-life, and no limit
  on how long-lived a target must be nor which radioactive beams exist: for a nucleus that is neither stable nor in the
  natural element the route says "not checked".
- **Whether a valence space suits the nucleus** (⁶²Ni in a space with a ⁵⁶Ni core has no excitations across Z = 28):
  the estimate tests the size of the calculation only.
- **Whether Coulomb excitation reaches the state at all** (its yield in a beam time): the planner works it out for a
  setup, and a route is not a setup.

## Depends on
The results ledger, registry and comparison engine (2026-11-24); the shell model (`physim.nuclear.shellmodel`) for
the effective charges.

## Done when
- Each of the three has a function with its reference, a validation check against published numbers, and an entry
  in the physics register.
- The registry lists it as a conversion, a route or a check, and the nucleus page shows it labelled as what it is
  (an expectation or a reading, not a measurement).
