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
