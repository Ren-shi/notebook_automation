# 57 · The run record: every number explains itself

**Priority:** P1 · **Size:** M · **Area:** Experiment workbench

## Why
A result is of little use to someone who cannot see how it was obtained. Not every user knows how a particle–γ
coincidence is counted, or how a peak area becomes a B(E2) and what that says about the nucleus. The report of
item 42 lists results; the link from each result to the theory behind it is missing.

## Scope
- **Reader:** an MSc student who knows nuclear physics and is new to Coulomb excitation.
- **An explanation for every calculated quantity,** in four parts:
  1. the formula;
  2. the same formula with this run's numbers substituted;
  3. what the quantity means physically;
  4. the assumptions made and a reference.
- **In the app:** each number in the side panel opens its explanation.
- **The run record,** viewable in the app and saved as a PDF, extends the report of item 42 with:
  - the setup, including the scene;
  - every piece of nuclear data and its provenance (ENSDF, derived, assumed, supplied by the user);
  - the method, step by step, from scattering to the Doppler-corrected spectrum;
  - the chain from peak area to yield, to excitation probability, to matrix element, to B(E2);
  - what B(E2) and the quadrupole moment mean: Weisskopf units, collectivity, β₂, prolate or oblate, lifetime;
  - the uncertainty budget;
  - what the simulation leaves out.
- **Left out:** full derivations, which belong in the theory pages of the documentation, linked from each
  explanation.

## Design notes
- The numbers in an explanation come from the calculation that produced the result, never from a second
  calculation, so the two cannot disagree.
- The report is HTML printed to PDF today. Keep one source for the app's view and the PDF.

## Depends on
42, 56.

## Done when
- Every number in the results panels and in the record has an explanation.
- A test recomputes each substituted formula from the numbers shown and obtains the reported result.
- A reader from the intended audience, other than the author, follows the B(E2) chain without help.
- The PDF and the view in the app show the same content.
