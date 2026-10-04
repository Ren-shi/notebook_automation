# 36 · Stopping power and energy loss

**Priority:** P1 · **Size:** L · **Area:** Nuclear planner · **Status: Done** (SRIM/TRIM comparison pending)

> **Done** (`python/physim/nuclear/stopping.py`, `python/physim/nuclear/data/nist_pstar.csv` and `nist_astar.csv`,
> `scripts/fetch_star_tables.py`, `scripts/lise_reference.py`, `tests/python/test_nuclear_stopping.py`,
> `tests/reference/nuclear/stopping_lise.csv`, `docs/theory/stopping.md`, `docs/physics-register/stopping.md`).
> - **Decision (with the user):** own models plus the NIST PSTAR/ASTAR (ICRU 49/90) tables, downloaded with
>   approval, rather than formulas only or wrapping a library.
> - `Stopping(ion, material)`: stopping power (electronic/nuclear/total), range, `energy_after` / `energy_loss`
>   through tilted layers, `straggling` (Bohr + relativistic factor + Lindhard–Scharff slow-ion factor, propagated
>   through thick layers with dσ²/dx = −2S′σ² + dΩ²/dx); `angular_straggling_mrad` (Highland).
> - p/α from the tables (isotopes at equal velocity; ranges anchored to NIST's CSDA range at 1 keV); untabulated
>   elements interpolated in ln I; Bragg additivity; heavy ions by Brandt–Kitagawa/Ziegler effective charge with a
>   Lindhard–Scharff low-velocity limit; ZBL nuclear stopping.
> - LISE++ comparison automated: its Excel library is called from Python (`scripts/lise_reference.py`, local
>   only); 1680 values (3 models × 7 ions × 8 targets × 10 energies) committed. Found on the way: until `set_loss`
>   is called the library's residual energies are ~8× too low, and `straggling_energy` is a FWHM.
>
> Results (31 tests): ranges match NIST CSDA ranges within 1% from 2 keV to 100 MeV (typically < 0.2%). Against
> LISE++ ATIMA 1.4, worst cases above 1 MeV/u: stopping 7% (p/α), 11–12% (heavy ions); energy after a fifth of the
> range 1.2% (p/α), 4.8–6% (heavy ions); straggling FWHM within 25% at ≥ 30 MeV/u. Interpolated elements within 7%
> above 0.5 MeV (leave-one-out). Building a table takes ~0.1 s.
>
> Not met / left out: the **SRIM/TRIM** comparison (compounds, slow heavy ions, transmitted-energy spread) needs
> runs by hand, requested in `tests/reference/nuclear/pending/stopping_srim.md`. The planned 2% (p/α vs NIST) holds for
> ranges; the 5% heavy-ion target holds only above 10 MeV/u: below 1 MeV/u heavy ions are within 25–55% (Kr/Xe ranges
> within a factor 2), the limit of a model without Ziegler's per-ion coefficient tables, which could be added
> later. Straggling below ~10 MeV/u differs from LISE++ by up to 2–4×. Channelling and charge-state distributions
> are ignored.

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
