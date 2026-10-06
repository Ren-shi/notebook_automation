# 57 · The run record: every number explains itself

**Priority:** P1 · **Size:** M · **Area:** Experiment workbench

**Status: Done, except the reading by a reader from the intended audience (the user).**

> **Done** (`python/physim/nuclear/record.py`, `Planner.explanations` and `Planner.record_html`, `record.html`
> in the report's zip, the **?** buttons and the record view in the app, `tests/python/test_nuclear_record.py`,
> guide section "The run record" in `docs/nuclear-setup.md`).
> - **Explanations:** eight for the plan (solid angle, rate, counts, beam time, excitation probability, γ-ray
>   efficiency, coincidences, Doppler shift) and nine for the analysis chain (area, yield, ⟨P⟩, B(E2), Weisskopf
>   units, β₂, Q₀ and Q_s, lifetime, uncertainty), each with the formula, the numbers substituted, the meaning,
>   the assumptions and a link to the theory page. The numbers are passed in from the calculation that produced
>   the result.
> - **One source for the app and the page:** `record_html` builds the page; the app shows its body with the
>   same style scoped to its box, and the report's zip holds it as `record.html` to print to PDF.
> - **Results:**
>   - A test recomputes every substituted formula from the numbers shown and gets the reported value to 10⁻⁹,
>     for the plan and for the analysis chain.
>   - The values are those of the results panels: the solid angle, counts, beam time, coincidence rate,
>     efficiency and Doppler shift equal the planner's own numbers; B(E2), W.u., β₂ and the uncertainty equal
>     the analysis result's.
>   - The record and the app's view are the same page (the app shows its body).
>   - Checked in the running app: the **?** by the γ-ray efficiency opens its explanation; the record shows on
>     the Report tab.
> - **Not done here, or to confirm:**
>   - A reader from the intended audience has not read it (the user, or a student they choose).
>   - Numbers in the result tables (Rates, Doppler table) have no **?** of their own; the side panel's and the
>     analysis's do. The explanations are for one detector at a time (the selected one).
>   - The theory links point at the documentation's pages, which hold the derivations only as far as they
>     exist today.

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
