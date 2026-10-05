# 50 · Level schemes from ENSDF

**Priority:** P1 · **Size:** L · **Area:** Experiment workbench

**Status: in progress.** Built (2026-10-05), and tested against hand-typed data in ENSDF's layout:
- **Reader** (`python/physim/nuclear/ensdf.py`): the adopted dataset of a nuclide from a local copy, zipped or not.
- **Level schemes** (`python/physim/nuclear/levels.py`): sources on every value, derived matrix elements, the
  user's values kept, truncation, and the `[levels.beam]` and `[levels.target]` sections of the setup file.
- **Planner and app:** look-up, level diagram, matrix-element table, and "Plan" for a chosen state.
- **Download:** `scripts/fetch_ensdf.py`, into `~/.physim/ensdf` (outside the repository) or `$PHYSIM_ENSDF`.

What remains:
- **The real ENSDF copy has not been downloaded,** so nothing is checked against real data yet. Six tests in
  `tests/python/test_nuclear_levels.py` wait for it, and so do the five nuclides to check by hand.
- **A table of moments.** Quadrupole moments are read only where ENSDF's adopted levels give them.
- **Editing in the app.** A matrix element can be changed in Python and in the setup file, not yet in the app.
- **The request to the NNDC** (the user).
- The rates and spectra still use the single state of `[reaction]`; items 54 and 58 use the whole scheme.

## Why
Typing a level scheme by hand is the most tedious part of setting up a Coulomb-excitation calculation, and the data
already exist. Today the setup file holds one excited state whose energy and B(E2) the user types (item 43). Every
later item (50–58) needs levels, transitions and matrix elements for both nuclei.

## Scope
- **Bundled data:** a processed copy of the ENSDF adopted levels and γ rays ships with the app, so the lookup works
  offline. The app may check for a newer copy when online, and never requires it.
- **What is read for each nucleus:**
  - levels: energy, spin and parity, half-life;
  - γ transitions: energy, branching, multipolarity, mixing ratio, conversion coefficient;
  - B(Eλ) values where evaluated.
- **Matrix elements:** the size of each reduced matrix element is derived from B(Eλ), or from the half-life and
  branching. ENSDF does not give signs, and gives few diagonal elements (quadrupole moments); those come from a
  separate table of moments where one exists, and are otherwise assumed and marked as such.
- **Provenance:** every value is marked as from ENSDF, derived, assumed or supplied by the user, with its
  uncertainty. The user's value always takes precedence.
- **Editing:** the user can change any value and add or remove levels and transitions.
- **Truncation:** by default, levels below an energy limit that connect to the ground state by E2 or E3, directly
  or in steps. The user can change the limit and the selection.
- **Data model:** any ground-state spin, any number of levels, matrix elements with signs. The setup file gains a
  level-scheme section; a setup with the single state of item 43 still loads.
- **In the app:** choosing the beam and the target shows both level schemes as diagrams, with the provenance of
  each value.
- **Left out:** fetching from the web on demand; decay data of radioactive nuclei; evaluating conflicting
  measurements ourselves.

## Design notes
- **Permission to bundle is not yet confirmed** (checked 2026-10-05).
  - Found: the data are free to download from the NNDC; there are citation guidelines; the work is funded by the US
    Department of Energy; open-source packages under MIT and BSD licences (radioactivedecay, PyNE, paceENSDF)
    already bundle data derived from ENSDF.
  - Not found: a written statement from the NNDC or the IAEA that permits redistribution.
  - **To do (the user, not yet sent):** ask the NNDC (nndc@bnl.gov) to confirm in writing, and keep the reply with
    the data sources. Nothing is bundled in the installer until the reply arrives.
  - **Until then (decided 2026-10-05):** the ENSDF copy is downloaded to this machine and kept in a folder that
    `.gitignore` excludes, so it is used for building and testing and never enters the public repository. The
    repository holds the script that downloads and processes it. Tests that need the data are skipped where the
    copy is absent, as on CI; small hand-typed level schemes cover the logic there.
  - **If permission is refused or does not arrive:** the app downloads the data from the NNDC on first use and
    keeps its own copy. It then works offline from the second start, and nothing is redistributed.
- **Citation:** follow the NNDC's guideline: the Nuclear Data Sheets evaluation for up to ten nuclides, otherwise
  the ENSDF database with the date of the copy. The run record (item 57) prints it.
- Record the terms of use of the moments table in the same way.
- Record the download size the bundle adds to the installer (item 44).

## Depends on
34.

## Done when
- Choosing ¹⁹⁴Pt gives its 2⁺ state at 328 keV, with a B(E2) that matches the evaluated value within its
  uncertainty.
- Five nuclei are checked by hand against the ENSDF web pages: even-even, odd-mass, light and heavy.
- A value changed by the user is used in every calculation and is marked as the user's in the report.
- The lookup works with the network disconnected.
- Setup files written before this item load unchanged.
