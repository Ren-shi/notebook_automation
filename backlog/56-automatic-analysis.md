# 56 · Automatic analysis: from peak areas to B(E2) and the shape

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

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
