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
`examples/nuclear/` in the repository.

From Python:

```python
from physim.nuclear import Experiment

exp = Experiment.load("examples/nuclear/alpha_on_gold.toml")
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
| `material` | yes | — | An element (`"Au"`), an isotope (`"208Pb"`), a formula (`"CD2"`, `"C10H8O4"`; `D` and `T` mean deuterium and tritium) or a named material (`"Mylar"`, `"Kapton"`, `"Polyethylene"`). |
| `thickness` | yes | areal density or length | Thickness along the target normal. |
| `tilt` | no | angle | Rotation of the target about the vertical (y) axis, between −90° and 90°. |

`[target.backing]` takes `material` and `thickness` for a backing layer on the downstream side of the target.

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

## When something is wrong

The planner reads the whole file and lists every problem at once, each naming the field:

```text
the setup has 2 problems:
  - beam: energy 5.5 has no unit; write it as text with a unit, e.g. "5.5 MeV"
  - detector 3 (D3): distance must be positive, got '-40 mm'
```

Misspelt field names are reported with a suggestion (`unknown field 'enrgy'; did you mean 'energy'?`).
