# 56 · Automatic analysis: from peak areas to B(E2) and the shape

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

**Status: Done.**

> **Done** (`python/physim/nuclear/analysis.py`, `Planner.analysis` and `Planner.analyse`, the Analyse block of
> the app's Spectra tab, `tests/python/test_nuclear_analysis.py`, guide section "Automatic analysis" in
> `docs/nuclear-setup.md`).
> - **Steps:** particle gate (inelastic group, chosen detectors and rings), Doppler correction, peak fit
>   (Gaussian on a line, then the counts within ±3σ over the line, so a flat-topped corrected peak is counted
>   whole), random coincidences subtracted, yield over efficiency × correlation × branch × gate share,
>   normalisation, B(E2), budget, shape, truth, beam time. Each run is kept with its settings; a change is
>   named with the previous B(E2).
> - **The inelastic gate** predicts each particle's measured energy as the generator does (two-body kinematics,
>   the way out of the target, the dead layer), and is as wide as the spread the reaction depth gives. In the
>   ⁵⁸Ni example the elastic and inelastic groups overlap for backward ¹⁶O, so the gate keeps about half of
>   the excitations; the share is known from the simulation and divided out, with a note.
> - **Results:**
>   - ⁵⁸Ni example, 400 000 reactions: B(E2↑) = 737 ± 2.8% (stat.) ± 7.8% (syst.) e²fm⁴ against 695 put in,
>     a pull of +0.7 with the gate and +0.7 without; 670 (−0.4) with the backward ring alone.
>   - Target normalisation, checked on ¹⁶O excited at 70 MeV with the ⁵⁸Ni 2⁺ state as reference: agrees with
>     the Rutherford normalisation to 1% and with the input within 2σ.
>   - Over six seeds the scatter of B(E2) equals the Monte Carlo uncertainty the result quotes (6.5% for
>     150 000 reactions), which is larger than the experiment's statistical uncertainty (2% for the 24 h run);
>     the two meet when the sample is as large as the run.
>   - Each systematic is found by changing its input by one standard deviation: beam energy 4.8% (the orbit
>     changes), detector positions 1.6%, efficiency and correlation as assumed.
>   - A changed fit window, ring selection or detector selection changes the result and is recorded.
> - **Not done here, or to confirm:**
>   - The reference transition of the target normalisation is analytic with Poisson noise; the simulation
>     excites one state.
>   - A high-energy γ ray on a large Compton continuum (the 6.9 MeV ¹⁶O line without a gate) is not always
>     found by the peak finder; with the gate it is.
>   - The gate's share of excitations comes from the simulation's truth; a real analysis would estimate it the
>     same way, from a simulation.
>   - "Matrix elements assumed" is zero in first order; item 58 fills it.

## Why
The purpose of the simulated experiment is the number it would measure. The app should analyse its own events as an
experimentalist would, and state how precisely the planned beam time determines the result.

## Scope
- **Automatic steps,** each shown and adjustable by the user:
  1. particle gates on the S3 (energy against ring), and the choice of rings;
  2. the Doppler correction;
  3. peak finding and fits, with a background under each peak;
  4. subtraction of random coincidences.
- **Yields:** peak areas corrected for efficiency, internal conversion and the angular correlation.
- **Normalisation,** suggested by the app and chosen by the user:
  - to a known transition in the target (for example the 2⁺ state of ¹⁹⁴Pt);
  - to the scattered particles in the same rings (Rutherford), for a target that is hardly excited, such as
    ²⁰⁸Pb.
- **B(E2):** in first order the excitation probability is proportional to B(E2), which gives the matrix element
  directly. Reported in e²b², e²fm⁴ and Weisskopf units. (Multi-step fitting is item 58.)
- **Uncertainties:**
  - statistical, from the peak areas and the background;
  - systematic: efficiency, target thickness, beam energy, the normalising B(E2), detector positions, and the
    matrix elements that were assumed.
  - A budget table lists each contribution.
- **Shape:**
  - the size of β₂ from B(E2);
  - the B(E2) in Weisskopf units and the ratio of the 4⁺ and 2⁺ energies, as signs of collectivity;
  - the quadrupole moment a rigid rotor of that B(E2) would have;
  - a statement that these readings depend on the rotor model.
- **Truth comparison:** the extracted B(E2) next to the value that was put in, with the difference in units of the
  uncertainty.
- **Beam time:** the counts in each peak per shift, and the time needed for a statistical uncertainty the user
  chooses.
- **Left out:** fitting more than one matrix element; the sign of the quadrupole moment (item 58).

## Depends on
55.

## Done when
- For a simulated experiment in the first-order regime, the extracted B(E2) agrees with the input within its
  uncertainty, with both normalisations.
- Over many seeds, the scatter of the extracted values matches the quoted statistical uncertainty.
- Each systematic contribution is checked by changing its input by one standard deviation and running again.
- A gate or fit range changed by the user changes the result, and the change is recorded.
