# 36 · Stopping power and energy loss

**Priority:** P1 · **Size:** L · **Area:** Nuclear planner · **Status: Done** (SRIM/TRIM compared 2026-11-14: compounds outside tolerance for heavy ions)

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
> **SRIM/TRIM comparison (2026-11-14,** `scripts/srim_reference.py`, `stopping_srim.csv`, `straggling_trim.csv`**):**
> SRIM 2013's SR Module and TRIM both run unattended from Python (SR Module has no window; TRIM with `TRIMAUTO` = 1).
> - Stopping, 5 ions × 5 targets × 0.1–10 MeV/u, judged by the LISE++ bands: Au and Si pass (worst 90% of a band,
>   Xe in Au at 1 MeV/u −10.7%). CsI passes (Xe at 10 MeV/u exactly at its 10%). **CD₂ and Mylar fail for heavy
>   ions:** CD₂ +31% (¹⁶O) and +39% (⁴⁰Ar) at 0.1 MeV/u, −14% (Ar) and −20% (Xe) at 1 MeV/u, −4.2% (¹⁶O) at
>   10 MeV/u; Mylar ⁴⁰Ar at 0.1 MeV/u +25.4% against 25%. α particles agree within 8% in every target. The gap is in
>   hydrogen: physim's heavy-ion stopping in pure H is +80–97% (0.1 MeV/u) and −25 to −35% (1 MeV/u) off SRIM's,
>   while carbon agrees within ~6%. SRIM's implied effective charge in H differs from its own in C (above 1 at
>   10 MeV/u). Ruled out: the target Fermi velocity (SRIM's v_F moves CD₂ by ≤ 0.3%). Against measured data (IAEA
>   database): in H₂/D₂ gas at 1–3 MeV/u physim is 17–27% low (RMS 23%, SRIM 8%); in (CH₂)ₙ polymers both codes are
>   about 10% RMS. The data are now a reference file (`stopping_measured_hydrogen.csv`) and a literature check,
>   `stopping_hydrogen_vs_data`, which fails for Kr in H₂/D₂ (−19%) and Kr, Xe in Mylar (−15%, −18%) at 1–10 MeV/u.
>   Carbon is within 7% of measured data; the hydrogen term is 6–60% low. A fitted factor on the H term fixes the
>   gas data but fails a held-out test, so the open fix is a hydrogen-specific effective charge from an independent
>   source, judged by that check unchanged.
>   **Done for Li–Ar (2026-11-15):** MSTAR 3.12 (ICRU 73's method) is ported (`physim.nuclear.mstar`) and used for
>   ions Li–Ar in hydrogen, ICRU 49 compounds and their deuterated versions (CD₂). RMS against the data for C–Ar:
>   H₂/D₂ 7.9% → 4.1%, (CH₂)ₙ 10.3% → 8.5%, Mylar 6.9% → 5.9%; Mylar now passes against SRIM, CD₂ still fails.
>   **Open:** Kr and Xe in hydrogen (still 15–19% low; outside MSTAR). Tried 2026-11-15: DPASS 21.06 and
>   SRIM-96 are as low as physim for Kr in H₂ (−18 to −28%); only SRIM 2013 matches (−8 to −11%); Xe in H₂ is
>   fine at 25 MeV/u, the only data. Kr/Xe in Mylar fail because of carbon and oxygen (Xe in carbon −20% at
>   0.1–1 MeV/u), not hydrogen. Using SRIM 2013's tables was rejected (fitted to the same data, unpublished,
>   circular). The facts (register, "Why physim is low for Kr in hydrogen"): in one experiment physim goes from
>   −20% in H₂ to +7% in Ar, so the effective charge depends on the target and physim's does not. **Next:** test
>   whether Kr's equilibrium charge state in H₂ is higher than ZBL's solid fit, against measured charge states.
>   SRIM was run with compound
>   correction 1 (Bragg additivity) throughout; SR Module cannot compute SRIM's Core-and-Bond factor.
> - SRIM's ranges are projected, physim's are path lengths, so they are reported, not judged.
> - TRIM: every exported TRIM.IN (4 cases and the 2 example setups) ran unchanged. Mean energies within 1% except α
>   in 10 µm Mylar (+1.05%; physim loses 4% less energy, likely TRIM's Bragg-additive Mylar against ICRU's); widths within 20%
>   except 660 MeV ¹³²Xe in 2 µm C (−63%: charge-exchange straggling is not modelled, see below).
> - Both comparisons are known failures in `tests/python/test_nuclear_validation.py`; the tolerances were not changed.
>
> Not met / left out: the 5% heavy-ion target against SRIM (see above). The planned 2% (p/α vs NIST) holds for
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
