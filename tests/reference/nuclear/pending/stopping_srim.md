# Request: stopping powers, ranges and straggling from SRIM/TRIM

For backlog item 36 (see `../README.md` for the file format). LISE++'s library has already been compared
automatically (`../stopping_lise.csv`), but it handles **elemental targets only**. These cases cover what it cannot:
compounds, and slow heavy ions, where the models differ most. Note the SRIM version used.

## 1. SRIM "Stopping/Range Tables" (`SR Module`)

For each ion and target below, ion energies **0.1, 0.3, 1, 3, 10 MeV/u** (in keV, as SRIM asks: multiply by the mass
number), stopping units **MeV/(mg/cm²)**. Save as `tests/reference/nuclear/stopping_srim.csv` with columns
`ion,target,density_g_cm3,e_mev_u,electronic_mev_mg_cm2,nuclear_mev_mg_cm2,projected_range_um,longitudinal_straggling_um`.

| Ions | Targets (density to enter in SRIM) |
|---|---|
| ⁴He, ¹²C, ¹⁶O, ⁴⁰Ar, ¹³²Xe | Au (19.32), Si (2.33), Mylar (SRIM compound "Mylar", 1.38), CD₂ (C₁D₂, 1.06), CsI (4.51) |

## 2. TRIM: energy spread after a layer

Ready-made inputs are in `trim/` (written by `physim.nuclear.export.trim_in`): copy `TRIM_n.IN` to `TRIM.IN` in the
SRIM folder and start TRIM, which reads it at start-up. Each runs 2000 ions, "Ion distribution and quick calculation
of damage", saving transmitted ions to `TRANSMIT.txt`. Report the mean and the standard deviation of the
transmitted energy. **Also note whether each file opened and ran without changes**: that is how the export itself is
checked (backlog item 40).

| # | File | Ion and energy | Layer | physim mean (MeV) | physim σ (MeV) |
|---|---|---|---|---|---|
| 1 | `trim/TRIM_1.IN` | ⁴He, 5.5 MeV | Mylar, 10 µm | 4.3653 | 0.0231 |
| 2 | `trim/TRIM_2.IN` | ¹⁶O, 64 MeV | Au, 0.5 µm | 62.2827 | 0.0571 |
| 3 | `trim/TRIM_3.IN` | ¹H, 10 MeV | Si, 300 µm | 7.2970 | 0.0850 |
| 4 | `trim/TRIM_4.IN` | ¹³²Xe, 660 MeV (5 MeV/u) | C, 2 µm (density 2.0) | 630.377 | 0.180 |

Save as `tests/reference/nuclear/straggling_trim.csv` with the usual header and the columns
`ion,material,density,energy_mev,thickness,mean_mev,sigma_mev,opened_unchanged` (thickness as written above, e.g.
`10 um`; density empty unless given; `opened_unchanged` yes/no). The validation check
`physim.nuclear.validation.spectra_vs_trim` reads it: mean within 1%, σ within 20%.

`trim/TRIM_alpha_on_gold.IN` and `trim/TRIM_oxygen_on_lead_array.IN` are the beams of the two example setups through
their targets (`export.trim_in_for`). They are there to check that multi-layer files open too; their results are
not needed.
