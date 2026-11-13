# 71 · Cascades: several γ rays per excitation

**Priority:** P2 · **Size:** L · **Area:** Experiment workbench (physics)

## Why
One γ ray per excitation is the first version (55). A populated state above the first one decays by a cascade,
and the level scheme (50) and the multistep solver (58) already know the populations and the branches. Without
cascades there are no true γ–γ coincidences and the spectra miss the lines an experimentalist expects.

## Scope
- The γ chain follows the decay path of each populated state through the level scheme's branches, with the
  angular correlation of each step (54) in the frame of the previous one; each γ ray goes to the crystals.
- True γ–γ coincidences appear in the matrix (55); gates on one line select another (67).
- Left out: internal conversion beyond the scheme's coefficients; lifetimes (73).

## Depends on
50, 54, 55, 63.

## Done when
- For a scheme with two levels populated, the γ–γ coincidence counts between the two lines agree with the
  populations, branches and efficiencies within statistics; with one level the results equal item 55's.

## Status
- Cascades through the level scheme, true γ–γ coincidences and gates on one line: done.
- The angular correlation between successive γ rays (each in the frame the previous one leaves): done, in
  `physim.nuclear.correlation`, tested in `tests/python/test_nuclear_correlation.py` (the textbook 4⁺ → 2⁺ → 0⁺
  correlation, mixed multipoles, singles unchanged).
- Still left out: the conversion electrons' own correlation, and a level's own kinematics (every excitation keeps
  the kinematics of the `[reaction]` state).
