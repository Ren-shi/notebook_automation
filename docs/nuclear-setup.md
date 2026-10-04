# Nuclear experiment setup files

The nuclear experiment planner starts from one description of an experiment: the beam, the target, the detectors
and how long the beam runs. It is a plain-text [TOML](https://toml.io) file that the planner app, Python scripts
and notebooks all read and write, so a setup can be saved, shared with a supervisor and attached to a beam-time
proposal.

```{note}
The planner is being built in stages (see the physics register for what is available). This page describes the
setup file, which is complete; kinematics, cross sections, energy loss and count rates are added on top of it.
```

## A first setup

```toml
schema = "physim.experiment/1"
title = "Alpha scattering on gold"

[beam]
nuclide = "4He"
energy = "5.5 MeV"
current = "1 pnA"

[target]
material = "Au"
thickness = "0.5 mg/cm2"

[run]
beam_time = "2 h"

[[detectors]]
name = "A30"
shape = "circle"
radius = "2.5 mm"
theta = "30 deg"
distance = "100 mm"
thickness = "300 um"
```

Each `[[detectors]]` block adds one detector; repeat it for as many as you need. Two complete examples are in
physim (`physim.nuclear.example_names()` lists them).

From Python:

```python
from physim.nuclear import Experiment

exp = Experiment.example("alpha_on_gold")    # or Experiment.load("my_setup.toml")
exp.beam.energy_mev                 # 5.5
exp.detectors[0].theta = "25 deg"   # edit in Python...
exp.save("my_setup.toml")           # ...checked again before saving
```

## Units

Every physical value is text with its unit, for example `"40 mm"`. A bare number such as `40` is rejected, because
a setup file should never leave the unit to guesswork.

| Quantity | Units accepted |
|---|---|
| Energy | `eV`, `keV`, `MeV`, `GeV` |
| Energy per nucleon | `keV/u`, `MeV/u`, `GeV/u`, `AMeV` |
| Length | `nm`, `um` (or `µm`), `mm`, `cm`, `m` |
| Areal density | `ug/cm2` (or `µg/cm2`), `mg/cm2`, `g/cm2` |
| Angle | `deg`, `rad`, `mrad` |
| Beam current (particles) | `ppA`, `pnA`, `puA` (or `pµA`) |
| Beam current (electrical) | `epA`, `enA`, `euA` (or `eµA`) |
| Time | `s`, `min`, `h`, `d` |
| Density | `g/cm3`, `mg/cm3`, `kg/m3` |
| Fraction | `%` |

Particle current counts beam particles (1 pnA = 6.24 × 10⁹ particles/s). Electrical current is what a Faraday
cup reads; the planner divides it by the charge state. Plain `nA` is rejected because it is ambiguous.

## Coordinates

The beam travels along **+z** through the target at the origin. **x** is horizontal and **y** is vertical. A
detector at polar angle **θ** (measured from the beam direction) and azimuth **φ** (from +x towards +y) sits at
distance **d** from the target. θ = 0° is straight downstream and θ = 180° straight upstream.

## Fields

### Top level

| Field | Required | Meaning |
|---|---|---|
| `schema` | yes | Always `"physim.experiment/1"` for this version. |
| `title` | yes | A short name for the setup. |
| `description` | no | Free text: what the experiment is for. |

### `[reaction]`

| Field | Required | Meaning |
|---|---|---|
| `type` | yes, if the section is given | `"elastic"` (the default when the section is left out). More reactions are added in later stages. |

### `[beam]`

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `nuclide` | yes | — | `"4He"`, `"He-4"`, `"16O"`, or `p`, `d`, `t`, `alpha`. |
| `energy` | yes | energy or energy per nucleon | Kinetic energy, total (`"64 MeV"`) or per nucleon (`"4 MeV/u"`). |
| `current` | yes | particle or electrical current | Beam intensity. |
| `charge_state` | no | — | Charge of the beam ions, 1 to Z. Fully stripped (Z) if left out. |
| `energy_spread` | no | energy or `%` | Energy spread of the beam (FWHM). |
| `spot_size` | no | length | Beam spot diameter on target (FWHM). |

### `[target]` and `[target.backing]`

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `material` | yes | — | See [Materials](#materials). |
| `thickness` | yes | areal density or length | Thickness along the target normal. A length needs a known density. |
| `density` | no | density | Overrides the tabulated density, or supplies one for a formula. |
| `tilt` | no | angle | Rotation of the target about the vertical (y) axis, between −90° and 90°. |

`[target.backing]` takes `material`, `thickness` and `density` for a backing layer on the downstream side of the target.

### `[[detectors]]`

Placement: give **either** `theta`, `distance` (and optionally `phi`), **or** `position`.

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `name` | no | — | A label; `D1`, `D2`, ... if left out. Names must be unique. |
| `shape` | yes | — | `"rectangle"`, `"annular"` or `"circle"`. |
| `theta` | see above | angle | Polar angle of the detector centre, 0° to 180°. |
| `phi` | no | angle | Azimuth of the detector centre; 0° if left out. |
| `distance` | see above | length | Distance from the target to the centre of the detector face. |
| `position` | see above | 3 lengths | `["x mm", "y mm", "z mm"]` of the centre of the detector face. |
| `facing` | no | — | Direction the sensitive face looks along, as three numbers `[x, y, z]`. Towards the target if left out. |
| `rotation` | no | angle | Rotation of the detector about its own normal (turns the strips). |
| `thickness` | yes | length or areal density | Active thickness; thinner than a particle's range means it punches through. |
| `material` | no | — | Detector material; silicon (`"Si"`) if left out. |
| `density` | no | density | Detector density, if not tabulated. |
| `dead_layer` | no | length or areal density | Inactive entrance layer. |
| `resolution` | no | energy | Energy resolution (FWHM). |
| `threshold` | no | energy | Lowest energy recorded. |

Size, depending on the shape:

| Shape | Fields |
|---|---|
| `rectangle` | `width`, `height` (required); `strips_x`, `strips_y`: number of strips on each side (1 if left out). |
| `annular` | `inner_radius`, `outer_radius` (required); `rings`, `sectors` (1 if left out). |
| `circle` | `radius` (required). |

### `[run]`

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `beam_time` | yes | time | How long the beam runs. |
| `counts_wanted` | no | — | Counts needed per detector; the planner reports the beam time this takes. |

## Materials

A material can be written as:

- **an element**, `"Au"`, with natural isotopic abundances;
- **an isotope**, `"208Pb"`, meaning 100% of that isotope. Its density is the natural element's, scaled by
  atomic mass, because the number of atoms per volume is the same;
- **a formula**, `"CD2"`, `"C10H8O4"`, `"13CH4"`. `D` and `T` mean deuterium and tritium; unmarked elements have
  natural abundances. Formulas have no tabulated density, so give `density` or write the thickness in `mg/cm2`;
- **a named compound** from NIST's table of 48 materials, by its short name (`"Mylar"`, `"Polyethylene"`,
  `"Polystyrene"`, `"PMMA"`, `"Teflon"`, `"PVC"`, `"CsI"`, `"LiF"`, `"GaAs"`, `"CdTe"`, `"Water"`, `"Air"`) or
  its full NIST name. `"Kapton"` is accepted as the formula C₂₂H₁₀N₂O₅, without a density.

Masses come from the AME2020 atomic mass evaluation, and abundances, densities and mean excitation energies from
NIST. The source of every number is listed in `physim/nuclear/data/SOURCES.md`. Tabulated densities are those of the
bulk material; NIST's carbon is graphite at 1.70 g/cm³. Thin foils can differ, so give `density` when you quote a
thickness in µm and the difference matters. From Python, `physim.nuclear.data` gives the same data directly:

```python
from physim.nuclear import data

data.q_value(["2H", "3H"], ["4He", "n"])          # 17.589 MeV
data.material("CD2", density="1.06 g/cm3").atoms_per_cm2("200 ug/cm2")
data.Material.enriched("C", {13: 0.99, 12: 0.01})   # isotopically enriched (Python only, for now)
```

## Rates, beam time and simulated spectra

Once a setup reads cleanly, the planner works out what it will measure:

```python
from physim.nuclear import Experiment, plot
from physim.nuclear.rates import Rates
from physim.nuclear.events import simulate

exp = Experiment.example("oxygen_on_lead_array")
r = Rates(exp)
r.per_detector()              # counts per second in each detector
r.beam_time_for("DSSD3")      # seconds for run.counts_wanted counts
r.peaks("DSSD1", (8, 8))      # peaks in one strip pair: mean energy, width and what makes it up
r.warnings()                  # too fast, too few counts, Rutherford's limits, detectors in the beam

ev = simulate(exp, events=1_000_000, seed=1)   # under a second
plot.spectra(ev)              # measured-energy spectra, counts in the planned beam time
plot.theta_energy(ev)         # energy against angle: the kinematic line of each channel
```

Rates count the scattered beam and the recoiling target nuclei, from the target and from its backing. The beam
current sets the absolute scale. `energy_spread` and `spot_size` (both FWHM) enter the simulation and the peak
widths. The same `seed` always gives the same events. {doc}`theory/events` explains the method and
{doc}`physics-register/rates-and-events` how it is checked.

## The beam-time report and data files

One command writes the report and every table:

```bash
python -m physim.nuclear.report my_setup.toml -o my-report --seed 1
```

From Python, `physim.nuclear.report.build(exp, seed=1).write("my-report")` does the same. The folder then holds:

- **`report.html`**: the whole report on one page, with figures embedded. It covers the setup, every warning, the
  detector table, kinematics, expected peaks, energy loss, simulated spectra, the validation status of each model,
  data sources and references to cite. Print it from the browser ("Save as PDF") for the PDF version; the page has
  print styles.
- **`figures/`**: the geometry, coverage, kinematics and spectra figures, as PNG (200 dpi) and PDF.
- **`setup.toml`**: the setup the report was made from. With the same seed it reproduces every number.
- **CSV tables**, with units in the column names:

| File | One row per | Columns |
|---|---|---|
| `detectors.csv` | detector | `theta_min_deg`, `theta_max_deg`, `theta_mean_deg`, `phi_min_deg`, `phi_max_deg`, `solid_angle_msr`, `dsigma_domega_lab_mb_sr` (scattered beam on the main target nuclide at the mean angle, mid-target energy), `rate_per_s` (counted above threshold), `rate_all_per_s` (everything reaching the face), `mc_rate_per_s` and `mc_rate_error_per_s` (Monte Carlo), `counts_in_run`, `beam_time_s` (for `counts_wanted`), `relative_error` |
| `strips.csv` | strip or ring–sector | `detector`, `segment_i`, `segment_j`, `rate_per_s` |
| `peaks.csv` | expected peak | `detector`, `channel` (nuclide and layer), `particle` (ejectile or recoil), `branch`, `mean_MeV`, `fwhm_keV`, `rate_per_s`, `above_threshold` |
| `kinematics.csv` | lab angle, 0–180° in 1° steps | `theta_lab_deg`, then `E_<particle>_MeV` for the scattered beam and the recoil of every target nuclide (higher-energy solution; empty past the maximum angle) |
| `energy_loss.csv` | layer | `layer`, `material`, `thickness_mg_cm2`, `energy_in_MeV`, `energy_out_MeV`, `loss_MeV`, `straggling_fwhm_MeV` (beam) |

## When something is wrong

The planner reads the whole file and lists every problem at once, each naming the field:

```text
the setup has 2 problems:
  - beam: energy 5.5 has no unit; write it as text with a unit, e.g. "5.5 MeV"
  - detector 3 (D3): distance must be positive, got '-40 mm'
```

Misspelt field names are reported with a suggestion (`unknown field 'enrgy'; did you mean 'energy'?`).
