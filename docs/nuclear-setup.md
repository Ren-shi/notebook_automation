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
| `type` | yes, if the section is given | `"elastic"` (the default when the section is left out) or `"coulex"` (Coulomb excitation of one state). |
| `excite` | no | `"target"` (default) or `"projectile"`: which nucleus is excited (`coulex` only). |
| `energy` | for `coulex` | Energy of the excited state, e.g. `"1.454 MeV"`. |
| `multipolarity` | no | `"E2"` (default), `"E1"` or `"E3"`: the transition from the 0⁺ ground state. |
| `emission` | no | `"correlated"` (default): the γ rays follow the particle–γ angular correlation. `"isotropic"`: they leave evenly in all directions. |
| `b_up` | for `coulex` | B(Eλ↑), the reduced transition probability up, with its unit: `"0.0695 e2b2"` or `"695 e2fm4"` for E2 (`e2b`/`e2fm2` for E1, `e2b3`/`e2fm6` for E3). Take it from ENSDF; physim does not look it up. |

A Coulomb-excitation setup:

```toml
[reaction]
type = "coulex"
excite = "target"
energy = "1.454 MeV"
multipolarity = "E2"
b_up = "0.0695 e2b2"
```

See {doc}`theory/coulex` for the physics, and the example `coulex_ni58`.

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
| `position` | no | length | Where the target sits along the beam, from the chamber's centre (0 if left out). Detectors are placed from the chamber's centre; every angle, distance and Doppler correction in the results is from the target. |
| `ladder` | no | pairs | A target ladder: `[["58Ni", "0.5 mg/cm2"], ["208Pb", "1 mg/cm2"]]`; `selected` (from 1) is the one in the beam and gives `material` and `thickness`. |
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

### `[[gamma_detectors]]` (optional)

γ-ray detectors, for Doppler shifts and broadening of the γ rays from an excited nucleus (`physim.nuclear.gamma`)
and for their efficiency and spectrum (`physim.nuclear.response`, see [γ-ray response](#γ-ray-response)). Each
faces the target: a single disc, or the crystals of a [model](#detector-models-and-the-chamber).

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `name` | no | — | A label. |
| `theta`, `phi` | `theta` yes | angle | Direction of the detector centre (φ = 0 if left out). |
| `distance` | yes | length | From the target to the detector face. |
| `radius` | yes | length | Radius of the face, which sets the opening angle. |
| `resolution` | no | energy | Resolution (FWHM) at 1332 keV; it is scaled to other energies and added to the Doppler broadening. |
| `material` | no | — | The crystal: `"Ge"` (if left out) or `"LaBr3"`. |
| `absorbers` | no | pairs | Material between the target and the detector: `[["Pb", "1 mm"], ["Cu", "0.5 mm"]]`. A thickness is a length or an areal density. |
| `efficiency_curve` | no | pairs | A measured full-energy-peak efficiency of the detector where it stands: `[["122 keV", "1.2 %"], ["1408 keV", "0.2 %"]]`, at least two points. |
| `efficiency` | no | fraction (`%`) | One full-energy-peak efficiency, used at every energy. |
| `threshold` | no | energy | Energy below which the crystal records nothing. |

Without `efficiency` or `efficiency_curve`, the efficiency comes from the typical response of such a crystal.

### `[run]`

| Field | Required | Unit | Meaning |
|---|---|---|---|
| `beam_time` | yes | time | How long the beam runs. |
| `coincidence_window` | no | time | Full width of the particle–γ coincidence window (`"100 ns"` if left out). |
| `dead_time` | no | time | A non-paralysable dead time per count; every count is scaled by the live fraction 1/(1 + τ × total rate). |
| `room_background` | no | rate (`/s`) | Room background in each γ-ray crystal: the lines of ⁴⁰K and the thorium and uranium series, this many counts per second together. |
| `extra_lines` | no | pairs | Lines added by hand to every crystal: `[["1274.5 keV", "0.5 /s"]]`. |
| `counts_wanted` | no | — | Counts needed per detector; the planner reports the beam time this takes. For Coulomb excitation these are particle–γ coincidences (excitation events seen with their γ ray), or excitation events if there are no γ detectors. |

## Detector models and the chamber

A detector can name a model of the catalogue. The model fills in every field the setup does not give, and any
field can be overridden:

```toml
[[detectors]]
name = "S3"
model = "S3"
theta = "180 deg"
distance = "30 mm"
thickness = "300 um"          # overrides the model's 1000 um

[[gamma_detectors]]
name = "Clover1"
model = "clover"
theta = "135 deg"
phi = "90 deg"
distance = "200 mm"           # from the target to the front face of the crystals
```

| Model | For | What it is |
|---|---|---|
| `S3` | `[[detectors]]` | Micron S3: an annular silicon detector with 24 rings and 32 sectors, active from 22 to 70 mm in diameter (rings 1 mm apart), on its circuit board |
| `clover` | `[[gamma_detectors]]` | a germanium clover of four crystals, each 50 mm in diameter and 70 mm long |
| `LaBr3_2x2` | `[[gamma_detectors]]` | a LaBr₃(Ce) scintillator, one crystal of 2 inches by 2 inches |

The chamber's `wall_material` and `wall_thickness` attenuate the γ rays of every detector outside it.

**Typical values.** Not every dimension is in a data sheet. Each model lists which of its values are typical and
should be replaced by those of the detector you have: `catalogue.model("S3").typical` in Python, and the `typical`
and `notes` entries in `physim/nuclear/data/detector_models.toml`, which also names the source of every number.
For the S3 these are the thickness, the dead layer and the size of the circuit board; for the clover, the distance
between crystal centres and the housing.

**Crystals.** A γ-ray detector with `crystals = 4` is a clover: four crystal faces in a square, `crystal_pitch`
apart, labelled A to D. Each crystal has its own row in the Doppler table, since each is corrected for the Doppler
shift on its own. A γ-ray detector without a model is a single disc of the `radius` given, as before.

**Dead material.** A model's circuit board, and the front of a γ-ray detector's housing (`housing_side`), stop
particles. A detector behind them is hidden in the rates and in the simulated events alike, and the warnings name
what hides it.

**The chamber.** An optional `[chamber]` section describes the scattering chamber:

```toml
[chamber]
radius = "120 mm"
wall_thickness = "3 mm"
wall_material = "Al"
beam_pipe_radius = "20 mm"
```

With a chamber, a particle detector that reaches beyond its radius, or a γ-ray detector whose front lies inside
its wall, is reported as a problem.

## γ-ray response

A γ-ray detector's efficiency and spectrum are built from a few numbers per crystal material; no photon is followed
inside the crystal. **These are typical responses, not the calibration of a real detector.** Give your own
`efficiency_curve` when you have measured one.

For a γ ray of energy E emitted at the target:

| Step | Rule | Where the numbers come from |
|---|---|---|
| Solid angle | the fraction of all directions the crystal faces cover | the setup |
| Attenuation | exp(−Σ μᵢ xᵢ) through the chamber wall (for a detector outside it), the `absorbers` and the housing window | NIST mass attenuation coefficients |
| Interaction | k (1 − exp(−μ L)) for a crystal of length L | NIST; k from the clover's data sheet (germanium), 1 for LaBr₃ (typical) |
| Peak-to-total | P/T(E) = min(0.95, P/T(1332 keV) (E / 1332 keV)^−q) | typical: 0.18 and q = 0.68 for germanium, 0.20 and 0.75 for LaBr₃ |
| The rest | a Compton continuum (Klein–Nishina, one scattering), a flat part up to the peak, escape peaks above 1.022 MeV | typical shares |
| Resolution | FWHM² = noise² + (FWHM(1332 keV)² − noise²) E / 1332 keV | the data sheets |

For germanium, k is fixed so that a 50 mm × 70 mm crystal has the 21.5% relative efficiency of Mirion's clover
sheet. The default clover then has a full-energy-peak efficiency of 0.099% at 1332 keV and 25 cm, all four
crystals together and without add-back. A detector given only a `radius` is taken as long as it is wide.

```python
from physim.nuclear import Experiment, response

exp = Experiment.example("coulex_ni58")
r = response.Response(exp, exp.gamma_detectors[0])
r.peak_efficiency(1.332)          # fraction of all γ rays emitted that end in the full-energy peak
r.total_efficiency(1.332)         # ... that leave any energy
r.fwhm(1.332)                     # MeV
response.transmission("Pb", "1 mm", 0.122)
```

**Calibration sources.** A source at the target position, in place of the beam, lets you measure the efficiency
of your arrangement as you would in the experimental hall:

```python
run = response.source_run(exp, "152Eu", activity="37 kBq", time="1 h")
run.efficiency_points("Ge90")     # per strong line: the efficiency from its peak area, and the one put in
response.source_run(exp, "60Co").peak_to_total("Ge90")
```

- The sources are ²²Na, ⁶⁰Co, ⁸⁸Y, ¹³³Ba, ¹³⁷Cs and ¹⁵²Eu, with the lines and intensities of the DDEP evaluations.
- The counts in each bin are drawn from a Poisson distribution, and the source decays during the run.
- A peak area is the counts within ±3σ less the level of the bands beside it. Lines that overlap a neighbour are
  left out of the efficiency points.
- Not included: two γ rays of one decay summing in a crystal, the room background, dead time, add-back.

Units of activity are `Bq`, `kBq`, `MBq`, `uCi` and `mCi`.

## Orientation and the particle–γ correlation

Coulomb excitation leaves the nucleus oriented: which magnetic substates are populated depends on the orbit, and
so on the angle the particle scattered to. The γ rays that follow are therefore not emitted evenly, and a crystal
sees more or fewer of them than its efficiency alone says. `physim.nuclear.orientation` keeps the amplitude of
every substate, in first order:

```python
from physim.nuclear.orientation import Excitation, simple_scheme

scheme = simple_scheme("58Ni", "1.454 MeV", "E2", "0.0695 e2b2")     # 0⁺ → 2⁺ → 0⁺
ex = Excitation("16O", "58Ni", "30 MeV", scheme, excite="target")
state = ex.at(150.0)                 # the projectile scattered to 150° in the CM frame
state.population(1)                  # probability that the 2⁺ state is populated
state.gamma_yield(1, 0)              # γ rays of 2⁺ → 0⁺ per collision
state.w_lab(1, 0, directions, phi_particle_deg=0.0, velocity=(0, 0, 0.03))   # per steradian, in the laboratory
ex.cross_sections()                  # mb, over all angles: each level directly and with feeding, each γ ray
```

- **Any level scheme.** `Excitation` takes a `LevelScheme` (see [Level schemes](#level-schemes)) and excites every
  level that an E1, E2 or E3 matrix element connects to the ground state, whatever the ground state's spin.
- **Decay.** Each level decays by its branching ratios. Internal conversion takes α/(1 + α) of a transition. A
  level fed from above takes over its parent's orientation, reduced as the theory of angular correlations gives.
  Mixed transitions use the mixing ratio (ENSDF's sign convention, that of Krane and Steffen).
- **In the laboratory.** The γ rays are thrown forward by the motion of the emitting nucleus, which also changes
  the solid angle.
- **In the rates.** `correlation_table(experiment)` in `physim.nuclear.gamma` gives, for every particle detector
  and γ-ray crystal, the γ rays seen in coincidence relative to isotropic emission. The coincidence rates and beam
  times use these factors. For the ⁵⁸Ni example a detector at 90° sees 0.4 times the isotropic number when the
  particle is in the backward ring, and those at 45° and 135° see 1.5 times.
- **Not included:** excitation in two or more steps and reorientation (first order only, so the matrix elements
  enter through their size and not their sign), deorientation in vacuum, and lifetimes: every state decays in
  flight with the orientation it was given.

## Particle–γ events

`physim.nuclear.gamma_events` follows the γ ray of every simulated Coulomb-excitation event whose particle
reached a detector:

```python
from physim.nuclear import Experiment
from physim.nuclear.gamma_events import simulate_gammas

g = simulate_gammas(Experiment.example("coulex_ni58"), events=400_000, seed=1)
g.counts("Ge90", "CD")                        # particle–γ coincidences in the run
g.spectrum("Ge90", "CD", corrected="recoil")  # counts per bin, random coincidences included
g.coincidences()                              # particle × γ matrix: true, random, and random γ–γ
g["corrected_recoil"], g["weight"]            # per γ ray (see gamma_events.GAMMA_COLUMNS)
```

- **The chain.** The excited nucleus leaves the target with its velocity (its own simulated path, or the
  two-body partner's of the particle that was detected), emits its γ ray in its rest frame with the angular
  correlation for that orbit, and the γ ray is transformed to the laboratory. A crystal it meets records it with
  the response: full energy, Compton deposit or escape peak, then the resolution and the threshold.
- **Doppler correction** from what the detectors know: the centre of the segment that fired, the centre of the
  crystal, and two-body kinematics at the nominal beam energy. It is made for the scattered beam and for the
  target recoil (`corrected_projectile`, `corrected_recoil`): the right one puts the peak at the transition
  energy, the wrong one smears it. The width that remains comes from the sizes of segments and crystals.
- **Coincidences.** Every γ ray is a true coincidence with its particle, weighted as the particle events are.
  Random coincidences follow from the singles rates and the window: 2τ × (particle rate) × (γ singles), the γ
  singles being all excitations (detected particle or not) times the crystal's total efficiency, plus the room
  background and the extra lines. True γ–γ coincidences need a cascade, which one excited state does not give;
  random γ–γ rates are given.
- **Dead time** scales every count by the live fraction.
- **Statistics.** Each excited event emits its γ ray ten times over (`gammas_per_event`), each with a share of
  its weight, since coincidences are rare. On the ⁵⁸Ni example 400 000 reactions take about 5 s and give about
  6 000 γ rays in coincidence.
- **Not included:** cascades (one γ ray per excitation), summing, pile-up, lifetimes, contaminant reactions.

In the app, **Simulate the γ rays** on the Spectra tab shows the particle × γ matrix and each detector's raw and
corrected spectrum. `write_root(..., gammas=g)` adds a `gammas` tree and `gamma_<detector>` histograms to the
ROOT file.

## Automatic analysis

`physim.nuclear.analysis` analyses the simulated events as an experimentalist would, and says how precisely the
planned beam time determines B(E2):

```python
from physim.nuclear.analysis import Analysis, Settings

a = Analysis(exp, simulate_gammas(exp, 400_000, seed=1))
r = a.run()                                # the default settings
r.b_e2fm4, r.statistical, r.systematic     # B(E2↑) and its relative uncertainties
r.budget                                   # each systematic contribution
r.steps                                    # what was done, step by step, in words and numbers
a.run(Settings(rings=range(4, 16), normalisation="target",
               reference={"energy": "328.5 keV", "b_up": "1.65 e2b2", "unc": 0.03}))
```

The steps, each adjustable through `Settings`:

1. **Particle gate:** the particle must sit in the inelastic group of its detector (its measured energy nearer
   the kinematic value for a particle that excited the state than for an elastic one, within `gate_width`
   resolutions), in the chosen `detectors` and `rings`. The share of excitations the gate keeps is known from the
   simulation and divided out, as an experimentalist takes it from one.
2. **Doppler correction** for the excited nucleus, or for the other one to see what that does.
3. **Peak fit:** a Gaussian on a straight line within `fit_half_width` resolutions, then the counts within ±3σ
   over the line, with the random coincidences expected in the window subtracted.
4. **Yield:** the area over the full-energy-peak efficiency, the γ-ray branch and the angular-correlation factor.
5. **Normalisation:** `"rutherford"`, to the elastic particles counted in the same rings (beam current and target
   thickness cancel); or `"target"`, to a `reference` transition of known B(E2↑). Since the simulation excites one
   state, the reference peak is its analytic expectation with Poisson noise.
6. **B(E2↑)** from the first-order proportionality of the excitation probability to B(E2), in e²fm⁴, e²b² and
   Weisskopf units; the statistical uncertainty from the fit and the normalisation; a budget of systematics
   (efficiency, angular correlation, beam energy, detector positions, the reference's B(E2)), each found by
   changing its input by one standard deviation and running again; and the Monte Carlo sample's own uncertainty,
   quoted separately.
7. **Shape:** β₂ from B(E2), the strength in Weisskopf units, E(4⁺)/E(2⁺) when a level scheme is there, and the
   quadrupole moment a rigid rotor would have, with the note that these readings depend on the rotor model.

The result also compares the extracted value with the one put in (the pull), and gives the counts in the peak per
8-hour shift and the beam time for the `wanted_precision`. Every run is kept in `Analysis.history`, and a changed
setting is named in the result's notes with the previous value of B(E2).

In the app, **Analyse** on the Spectra tab runs it with the choices made in the row above the button.

## The run record

Every number of the plan and of the analysis comes with an explanation in four parts: the formula, the formula
with this run's numbers in it, what the quantity means physically, and the assumptions with a link to the theory
page. The numbers are those of the calculation that produced the result, and each explanation carries a rule that
recomputes its value from the numbers shown, which the tests check.

```python
from physim.nuclear.record import explanations, analysis_explanations, record_html

for x in explanations(planner, "CD", "Ge45"):   # solid angle, rate, counts, beam time, P(θ), ε, coincidences, Doppler
    print(x.title, x.value, x.unit)
    print(x.formula); print(x.substituted); print(x.meaning)
analysis_explanations(planner, result)          # area → yield → ⟨P⟩ → B(E2) → W.u., β₂, Q₀, lifetime, uncertainty
page = record_html(planner, result)             # the record as one page, printable to PDF
```

The record holds the setup (with the scene), every piece of nuclear data with its provenance (ENSDF, derived,
assumed, user), the method step by step from scattering to the Doppler-corrected spectrum, every number explained,
the chain from the peak area to B(E2) and the shape, the uncertainty budget, and what the simulation leaves out.
It is written for a reader who knows nuclear physics and is new to Coulomb excitation.

In the app, the **?** beside a number in the scene's side panel and in the analysis opens its explanation, and
**Show the run record** on the Report tab shows the page; the report's zip holds the same page as `record.html`.

## Multi-step excitation and reorientation

`physim.nuclear.coupled` integrates the coupled equations for the amplitudes of every magnetic substate of every
level along the orbit, so excitation in several steps, the reorientation effect of the quadrupole moments and
interference between paths are all included. It has the interface of first order, so the populations, tensors,
decay and correlation work unchanged:

```python
from physim.nuclear.coupled import CoupledChannels

cc = CoupledChannels("16O", "194Pt", "60 MeV", scheme, excite="target")
cc.probabilities(150.0)           # every level, summing to one
cc.at(150.0).gamma_yield(1, 0)    # with the decay through the scheme
```

- **The orbit** is the symmetrised one of the ground state and a reference level (the lowest with a matrix element
  to the ground state); with weak coupling the solution is first order exactly.
- **Matrix elements** enter with their signs, from the level scheme; a diagonal E2 element is the quadrupole
  moment's, through `quadrupole_factor`.
- **Speed:** a Runge–Kutta integration whose steps follow the fastest phase. 25 levels (265 substates) take
  about 0.4 s per angle; a grid of 45 angles and 2 energies for an experiment about 40 s.

`physim.nuclear.multistep` applies it to the planned experiment:

```python
from physim.nuclear.multistep import Multistep

ms = Multistep(exp, scheme, role="target")
y = ms.yields()                     # γ rays per second of each transition in each detector and ring
ms.shapes(level=1)                  # Q(2⁺) prolate, zero, oblate: can the run tell them apart?
ms.fit({"CD": {(1, 0): (counts, unc)}}, free=[(0, 1, "E2"), (1, 1, "E2")])
ms.gosia_input(y)                   # a GOSIA input file for this setup, with the yields
```

- The probabilities are solved on a grid of CM angles at a few beam energies through the target and interpolated
  to every quadrature direction of every segment, weighted by the Rutherford cross section of the symmetrised
  orbit (the analytic rates use the elastic one, 6% apart at the ⁵⁸Ni example's energy).
- **Prolate or oblate:** the same experiment with the 2⁺ state's quadrupole moment at the rigid-rotor value,
  at zero and at its opposite; the yields per ring are compared with the counting uncertainties of the run.
- **The fit** adjusts up to three matrix elements by Gauss–Newton least squares to measured counts; its
  uncertainties come from the curvature at the minimum.
- **GOSIA** is not part of physim. The input file follows the published format so that GOSIA, run elsewhere, can
  check these numbers or fit the matrix elements; it has not been checked against a run, since no GOSIA
  installation exists on the development machine.

In the app, **Solve with all orders** in the "Excitation and γ rays" tab (once a level scheme is looked up) shows
the yields with all orders and with first order, the prolate–zero–oblate comparison per ring, and offers the
GOSIA input file.

Not included: the E1 polarisation correction, deorientation, lifetimes, mutual excitation of both nuclei, and
nuclear interference above the safe energy.

## Tracks of simulated events

`physim.nuclear.tracks` turns a sample of simulated events into straight paths for the scene:

```python
from physim.nuclear.tracks import sample_tracks, describe

tracks = sample_tracks(exp, gammas, n=40, select="coincidences", weighted=True)
tracks[0].paths          # beam, scattered beam, recoil and γ ray, each with the hit it ends in
tracks[0].particles      # the event's rows of the particle events; .gammas its γ rays
describe(tracks, "coincidences", True, gammas.events.n_events)
```

Each particle's path runs from the target to where its track crosses the detector face (inside the segment the
event record names); an undetected partner runs to the edge of the scene; a γ ray runs to the crystal it hit.
With `weighted=True` events are picked with probability proportional to the rate they stand for, so the sample
looks like a run; with `weighted=False` they are picked as generated, so the rare large-angle scatterings show.
The animation itself runs in the browser.

## A misplaced target

The analysis may assume the target somewhere else than it is. `physim.nuclear.alignment` shows what that does and
finds the offset back, as one does with real data:

```python
from physim.nuclear.alignment import diagnostic, fit_offset, overlay, with_offset

assumed = with_offset(exp, 2.0)                  # the target assumed 2 mm further along the beam
overlay(gammas, exp, assumed)                    # the corrected peak with each geometry: shift and broadening
diagnostic(gammas, assumed)                      # centroid against ring, per crystal: flat when right, sloped when not
fit_offset(gammas, assumed)                      # the offset that flattens it, with its uncertainty
```

- `gamma_events.recorrect(gammas, experiment)` makes the Doppler correction again with another geometry.
- The diagnostic's centroids are judged against the pattern of a simulation made with the assumed geometry and
  corrected with it, because large crystals and the energy loss leave the corrected centroids a little
  ring-dependent even when the geometry is right.
- `Settings.target_offset_mm` in the analysis makes it assume the offset; the budget's "detector positions" entry
  includes the change of the peak area when the target moves by `position_unc_mm`.
- In the app, **Check the alignment** on the Spectra tab shows the overlay, the diagnostic plot and the fitted
  offset, and the scene shows the assumed target as a faint outline.

A calibration source can be placed anywhere: `source_run(exp, "152Eu", position=(0, 60, 0))` puts it 60 mm along
y from the target, and each crystal's rate follows its solid angle from there.

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

## Level schemes

A setup can carry the level scheme of the beam nucleus and of a target nucleus: levels, γ-ray transitions and
reduced matrix elements. Schemes are read from ENSDF, so they need not be typed.

**Getting ENSDF.** physim does not ship ENSDF. Download a copy once (about 40 MB), into `~/.physim/ensdf` or the
folder named by the `PHYSIM_ENSDF` environment variable:

```
python scripts/fetch_ensdf.py
```

**Looking up a scheme.** In the app, the *Excitation and γ rays* tab has a *Look up in ENSDF* button for the target
and for the beam. Each matrix element is shown in a field: type a new value with its unit to use your own. In
Python:

```python
from physim.nuclear.levels import LevelScheme

pt = LevelScheme.from_ensdf("194Pt", max_energy_kev=1500)
pt.levels[1].energy              # Value(328.464, 0.012, "ensdf")
pt.b(0, 1, "E2")                 # B(E2; 0+ → 2+) in e² fm⁴
pt.b_weisskopf(1, 0, "E2")       # B(E2; 2+ → 0+) in Weisskopf units
exp.levels["target"] = pt        # saved with the setup
```

**Where each value comes from.** Every value has a source:

| Source | Meaning |
|---|---|
| `ensdf` | as evaluated in ENSDF |
| `derived` | computed by physim from ENSDF values; the note says how |
| `assumed` | chosen because the data do not say |
| `user` | set by you; it takes precedence, and a new look-up keeps it |

**Matrix elements.** ENSDF gives transition strengths, not matrix elements, so physim derives the size of each
E1, E2 and E3 matrix element:

- from B(Eλ) in Weisskopf units where ENSDF gives it;
- otherwise from the level's half-life and its γ-ray branching, with conversion coefficients and mixing ratios;
- a quadrupole moment in ENSDF becomes the diagonal E2 matrix element of its level. Where ENSDF writes the moment
  without a sign, the note says that the sign is assumed.

ENSDF gives no signs. Derived matrix elements are positive, and the note says the sign is assumed.

**Which levels are kept.** A look-up keeps the levels below the energy limit that are joined to the ground state by
a chain of matrix elements, and every level their γ rays feed. Levels whose energy or spin ENSDF does not pin down
cannot carry a matrix element.

**In the setup file.** A scheme is written as `[levels.target]` or `[levels.beam]`:

```toml
[levels.target]
nuclide = "58Ni"

[[levels.target.level]]
energy = "0 keV"
jpi = "0+"

[[levels.target.level]]
energy = "1454.21 keV"
jpi = "2+"
half_life = "6.52e-13 s"

[[levels.target.transition]]
from = 1
to = 0
energy = "1454.2 keV"
intensity = 100
multipolarity = "E2"

[[levels.target.matrix_element]]
from = 0
to = 1
multipolarity = "E2"
value = "0.26 eb"
```

Levels are numbered from 0 in order of energy, the ground state first. A matrix element is written in `efm2` or
`eb` for E2, `efm` or `eb0.5` for E1, and `efm3` or `eb1.5` for E3. A value without a `source` is taken as yours.

The rates and spectra still use the single state of `[reaction]`; *Plan* in the app (or `Planner.use_state`) fills
that state in from a scheme.

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
- **`figures/`**: the geometry, coverage, kinematics and spectra figures, as PNG (600 dpi) and PDF. They are
  drawn in a journal's style at its column width: Physical Review, one column, unless you say otherwise with
  `--journal nature --width double` (or `build(exp, journal="nature", width="double")`). The styles are
  `physical_review`, `nature`, `science`, `elsevier` and `springer`; see {doc}`planner-app` for their widths.
- **`setup.toml`**: the setup the report was made from. With the same seed it reproduces every number.
- **`events.root`** (when `uproot` is installed: `pip install physim-engine[root]`), for analysis in ROOT. It holds:
  - the simulated particles as a TTree `events`, one branch per column (detector, strip, energies, angles, depth,
    weight);
  - a TH1D `spectrum_<detector>` per detector, in counts for the planned run, with Monte Carlo errors;
  - the setup file, as `setup`.

  In ROOT: `events->Draw("measured", "weight*(detector==2 && counted)")` or `spectrum_A45->Fit("gaus")`. See
  `physim.nuclear.rootio`.
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
