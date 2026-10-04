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

Full-cascade off ("Ion distribution and quick calculation of damage"), 2000 ions, with "transmitted ions" saved
(`TRANSMIT.txt`). Report the mean and the standard deviation of the transmitted energy.

| # | Ion and energy | Layer |
|---|---|---|
| 1 | ⁴He, 5.5 MeV | Mylar, 10 µm |
| 2 | ¹⁶O, 64 MeV | Au, 0.5 µm |
| 3 | ¹H, 10 MeV | Si, 300 µm |
| 4 | ¹³²Xe, 660 MeV (5 MeV/u) | C, 2 µm (density 2.0) |

physim's values for the TRIM cases are listed when the request is run, by `Stopping(...).energy_after` and
`.straggling` on the same layers.
