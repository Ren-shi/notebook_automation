# 46 · Guided workflow, field help and "how to read this"

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner

**Status: Done, except a trial with a new user.**

## Why
From the first look at the app (2026-10-05): it does everything it needs to, but it reads like LISE++ or SRIM.
You enter values and everything is shown at once, so a new user needs the manual open beside them. Users don't
know:
- what each value means;
- what changing it does;
- what to look for in each plot.

## Scope
**Guided mode, the default.** A stepper in the left panel works through the plan in order, with the result that
step affects shown on the right:

| Step | Inputs | Result shown |
|---|---|---|
| 1. What do you want to measure? | Rutherford scattering or Coulomb excitation; starting example | — |
| 2. Beam | nuclide, energy, current | Kinematics |
| 3. Target | material, thickness, backing | Energy loss |
| 4. Particle detectors | angles, sizes, thresholds | Geometry, then Kinematics |
| 5. γ-ray detectors (Coulomb excitation only) | the crystals | Excitation and γ rays |
| 6. Rates and beam time | beam time, counts wanted | Rates and beam time |
| 7. Spectra | — | Spectra |
| 8. Report | — | Report |

- Each step opens with one line on what you are deciding, and starts filled with values that work.
- **Next** is blocked while the step has a setup error; the error is shown inside the step.
- Steps can be revisited in any order once reached.
- An **Expert view** switch in the header gives the current layout. The browser remembers the choice.

**Field help.** Each input has a help icon with three short lines:
- what the value is;
- its typical range;
- what happens to the results if you raise it.

The help is practical, not physics teaching; the derivations stay in the Explain panels. All of it lives in one
table in `planner.py`, which a test checks covers every field.

**"How to read this."** One to three sentences above each result, written for the current setup with its numbers.
For example: "The excited ⁵⁸Ni line sits 0.64 MeV below the elastic one in DSSD-L; with 50 keV resolution they
separate easily." These come from `Planner.reading(tab)`, so they can be tested without the GUI.

## Depends on
41, 45.

## Done when
- From a fresh start, a new user can reach a report by pressing **Next** through the steps without opening the
  docs.
- Every input has help (tested), and every tab has a reading for every example (tested).
- The expert view is unchanged.

> **Done** (`python/physim/nuclear/guide.py`: `HELP`, `STEPS`, `GOALS`, `steps_for`, `reading`; `app.py`: the guided
> stepper, help icons, "How to read this" boxes, the Guided/Expert switch remembered in the browser's local storage
> and `?mode=` in the address).
> - **Tests** (`tests/python/test_nuclear_guide.py`):
>   - every input of the setup panel has help;
>   - every tab has a reading for every example, with no exponent-notation numbers;
>   - the readings quote this setup's numbers and change their verdict with it;
>   - the steps differ for elastic and Coulomb excitation, and each kind of setup problem lands in its own step.
>
>   `test_the_app_serves_every_example` loads every example in both modes.
> - **Checked by hand in a browser:**
>   - the eight Coulomb-excitation steps through to the report;
>   - the help tooltips;
>   - switching modes keeps the setup and survives a reload.
> - **Readings found real points about the examples:**
>   - in `coulex_ni58` the elastic and excited peaks overlap over the whole CD (2.9 MeV wide, 0.88 MeV apart), but
>     each ring sees a narrow range, so the reading says to separate them ring by ring;
>   - `alpha_on_gold`'s A20 is above the pile-up rate.
> - **Also fixed:** setup errors suggested "40 mm" as the example unit for every field. They now suggest a unit for
>   that field ("'12' has no unit; write it with one, e.g. '12 h'").
> - **Waiting on:** a trial by someone new to the planner, as for item 41.
