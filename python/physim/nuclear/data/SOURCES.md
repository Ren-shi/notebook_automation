# Nuclear and material data: sources and licences

Every data file in this folder, where it came from, and the terms it is redistributed under. Retrieved
2026-10-04. The NIST CSV files were produced from the downloads by `scripts/convert_nuclear_data.py`, which
changes the layout only (HTML or text listing → CSV) and copies the values as published.

## `mass_1.mas20.txt`: atomic masses (AME2020)

- **Source:** Atomic Mass Data Center, IAEA Nuclear Data Section,
  <https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt>. Shipped **unchanged**.
- **Contents:** mass excess, binding energy per nucleon, β-decay energy and atomic mass, with uncertainties, for
  every nuclide in the evaluation. `#` in place of a decimal point marks an estimated (non-experimental) value.
- **Cite:**
  - W. J. Huang, M. Wang, F. G. Kondev, G. Audi and S. Naimi, "The AME 2020 atomic mass evaluation (I)",
    Chinese Physics C **45**, 030002 (2021).
  - M. Wang, W. J. Huang, F. G. Kondev, G. Audi and S. Naimi, "The AME 2020 atomic mass evaluation (II). Tables,
    graphs and references", Chinese Physics C **45**, 030003 (2021), doi:10.1088/1674-1137/abddaf.
- **Licence:** the evaluation is published open access under the Creative Commons Attribution 3.0 licence; further
  distribution must keep the attribution above (authors, title, journal citation and DOI).

## `nist_isotopes.csv`: natural isotopic compositions

- **Source:** J. S. Coursey, D. J. Schwab, J. J. Tsai and R. A. Dragoset, "Atomic Weights and Isotopic Compositions
  with Relative Atomic Masses", NIST Physical Measurement Laboratory (last updated January 2015),
  <https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=all>.
  Isotopic compositions there are from M. Berglund and M. E. Wieser, "Isotopic compositions of the elements 2009".
- **Modified:** reduced to the isotopes that have a natural isotopic composition (288 rows), as CSV.
  The relative atomic masses are kept for reference only; physim takes masses from AME2020.

## `nist_elements.csv`: element densities and mean excitation energies (Z = 1–92)

- **Source:** J. H. Hubbell and S. M. Seltzer, "Tables of X-Ray Mass Attenuation Coefficients and Mass
  Energy-Absorption Coefficients from 1 keV to 20 MeV for Elements Z = 1 to 92 and 48 Additional Substances of
  Dosimetric Interest", NIST Standard Reference Database 126 (last updated July 2004), doi:10.18434/T4D01F,
  Table 1, <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html>.
- **Modified:** converted to CSV, with the standard atomic weight of each element added from the isotopic
  composition listing above. NIST notes that some densities are nominal, and those for Z = 85 and 87 were set to
  10 g/cm³ arbitrarily.

## `nist_compounds.csv`: compounds and mixtures

- **Source:** the same database, Table 2, <https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html>
  (48 materials, with composition by mass fraction).
- **Modified:** converted to CSV.

## `nist_pstar.csv` and `nist_astar.csv`: proton and α stopping powers and ranges

- **Source:** M. J. Berger, J. S. Coursey, M. A. Zucker and J. Chang, "Stopping-Power & Range Tables for Electrons,
  Protons, and Helium Ions", NIST Standard Reference Database 124 (last updated July 2017), doi:10.18434/T4NC7P,
  PSTAR and ASTAR, <https://physics.nist.gov/cgi-bin/Star/ap_table.pl>. The methods are those of ICRU Reports 37 and
  49; graphite, air and water were re-evaluated following ICRU Report 90.
- **Contents:** for 74 materials (26 elements, 48 compounds and mixtures), on NIST's default energy grid (protons
  1 keV to 10 GeV, 132 energies; α particles 1 keV to 1 GeV, 121 energies): electronic, nuclear and total stopping
  power (MeV cm²/g), CSDA range and projected range (g/cm²), detour factor.
- **Modified:** retrieved 2026-10-04 by `scripts/fetch_star_tables.py`, one request per material, and converted from
  the HTML answers to CSV; values copied as published.

## NIST terms

Data from NIST Standard Reference Databases: ©Copyright by the U.S. Secretary of Commerce on behalf of the
United States of America. All rights reserved. NIST data that are not Standard Reference Data are works of the
U.S. Government and not subject to copyright in the United States. Both are redistributed here with
attribution to NIST and with the modifications noted above. NIST provides the data "AS IS" and makes no warranty
of any kind, express or implied; NIST does not warrant that the data are accurate or fit for any purpose and is not
liable for any damages arising from their use. Contact: NIST Physical Measurement Laboratory, <https://www.nist.gov/pml>.

## ENSDF: level schemes (not in this folder)

- **Not shipped.** Level schemes (`physim.nuclear.levels`) are read from a copy of the Evaluated Nuclear Structure
  Data File that each user downloads with `scripts/fetch_ensdf.py` into `~/.physim/ensdf` (or `$PHYSIM_ENSDF`).
  The data are free to use, but no statement that permits redistribution was found (checked 2026-10-05), so the
  copy is kept out of the repository and out of the installer until the NNDC confirms in writing.
- **Source:** National Nuclear Data Center, Brookhaven National Laboratory, on behalf of the international Nuclear
  Structure and Decay Data network, <https://www.nndc.bnl.gov/ensdfarchivals/>, doi:10.18139/nndc.ensdf/1845010.
- **Used:** the "ADOPTED LEVELS, GAMMAS" dataset of a nuclide: level energies, spins, parities, half-lives and
  quadrupole moments; γ-ray energies, intensities, multipolarities, mixing ratios, conversion coefficients and
  B(Eλ) in Weisskopf units. Values are read as published; what physim derives from them is marked "derived".
- **Cite:** for up to ten nuclides, the Nuclear Data Sheets evaluation named in each dataset (physim keeps it as
  the scheme's reference); otherwise the ENSDF database with the date of the copy.
