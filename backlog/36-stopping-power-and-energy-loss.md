# 36 · Stopping power and energy loss

**Priority:** P1 · **Size:** L · **Area:** Nuclear planner

## Why
Beam particles lose energy on the way into the target and scattered particles lose energy on the way out, through
backings and through detector dead layers. Energy loss shifts and broadens every peak in a spectrum and decides
whether a particle stops in a detector or punches through. This is the part students currently do in SRIM by hand.

## Scope
- Electronic stopping power dE/dx for any ion in any material from item 34:
  - Bethe–Bloch with shell, Barkas and density-effect corrections at high energy (above roughly 1 MeV/u).
  - A semi-empirical low-energy model for protons and α particles, and effective-charge scaling for heavier ions.
  - Nuclear stopping at low energy.
  - Bragg additivity for compounds.
- Range, and energy after a given thickness by integrating dE/dx (not the thin-layer linear approximation), with
  tilted layers handled through the path length.
- Energy straggling (Bohr, with a correction for thick layers) and angular straggling (multiple scattering).
- **Decision to make first:** implement the models ourselves, or adopt an existing open-source stopping library and
  wrap it. Check the licence of each candidate before choosing. The choice and the reasons go in the physics register.
- Left out: channelling, charge-state distributions (a fixed or equilibrium mean charge is used), nuclear reactions
  in the target.

## Depends on
34.

## Done when
- Proton and α stopping powers and ranges agree with NIST PSTAR/ASTAR across their tabulated energy range, within
  a tolerance fixed after the first comparison (target: 2%) and stated in the register.
- Heavy-ion stopping (at least C, O, Ar, Kr, Xe in Au, C, Si, Mylar) agrees with SRIM within a stated tolerance
  (target: 5%), with any systematic difference explained.
- Energy loss through a thick layer agrees with SRIM/TRIM for the same layer.
- Straggling widths are compared with TRIM output and the comparison is shown in a validation notebook.
- Theory note in `docs/theory/`, and an entry in the physics register.
