# Lifetime measurements

physim plans, simulates and analyses measurements of the lifetimes of excited states by the two Doppler-shift
methods, after Coulomb excitation or fusion-evaporation:

- **Recoil distance (RDDS, a plunger)**, for mean lives of about 1 ps to 1 ns. The excited nucleus flies out of the
  target across a gap of a few µm to a few mm into a stopper (or a degrader). γ rays emitted in flight are
  Doppler-shifted; those emitted after the stopper are not. The share of each, against the distance, gives τ.
- **Doppler-shift attenuation (DSAM)**, for about 10 fs to 1 ps. The nucleus slows down and stops in the target's
  backing; the γ-ray line shape records how fast it was moving when it decayed.

A lifetime is the independent check on a Coulex B(E2): the same matrix element from a different observable. The
physics is in [the theory note](theory/lifetimes.md); the checks against outside tools and published experiments are
in the [physics register](physics-register/lifetimes.md).

## Setting up the measurement

A `[lifetime]` section goes with a `coulex` or `fusion-evaporation` reaction:

```toml
[lifetime]
method = "rdds"                        # or "dsam": the target's [target.backing] slows the recoil
distances = ["10 um", "20 um", "50 um", "100 um", "200 um", "500 um"]
offset = "8 um"                        # the distance at electrical contact, added to every distance
deorientation = { alpha2 = 0.18, gamma2 = "0.142 /ps", alpha4 = 0.24, gamma4 = "0.314 /ps" }   # optional

[lifetime.stopper]
kind = "stopper"                       # "stopper", "degrader" or "differential" (the first of two foils)
material = "Cu"
thickness = "15.7 mg/cm2"

# [lifetime.degrader2]                 # a differential plunger's second foil (CoCoDiff)
# material = "Ta"
# thickness = "6.5 mg/cm2"
# gap = "100 um"
```

From Python, the planner's draft methods make the same edits (each is checked; problems go to `p.problems`):

```python
from physim.nuclear.planner import Planner

p = Planner.example("coulex_ru98_plunger")
p.set_lifetime(distances=[3, 10, 25, 60, 160, 450], offset="8 um")   # numbers in µm, or text with units
p.set_stopper("stopper", "Cu", "15.7 mg/cm2")
p.set_lifetime(deorientation={})                                      # {} removes it
p.set_lifetime(method="dsam")      # the plunger's fields are kept aside, and come back with method="rdds"
p.add_ancillary("plunger")         # a fusion setup gets a default plunger the same way
```

The setup's problems say when the recoil stops in the target, when the stopper does not stop it, or when a degrader
does. In the app, the plunger card in Setup edits the same section.

## Planning: which method, which distances, how much beam

```python
plan = p.lifetimes()                # every level with a strong line, the studied one first
plan = p.lifetimes(level="2+")      # one level: its name, its number in the scheme, or an energy ("652 keV")
```

The nucleus is Coulomb excitation's excited one (the projectile or the target) or fusion-evaporation's studied
residue (`[run] channel`, else the largest cross section). For each level the plan gives:

| Key | What it is |
|---|---|
| `tau_ps`, `tau_source` | τ from the level's half-life (ENSDF, typed) or else its B(Eλ↓) (the setup's Coulex B(E2), the scheme's matrix elements, the shell model), with branching and conversion |
| `method` | "rdds" for τ of about 1 ps to 1 ns, "dsam" for 10 fs to 1 ps, else "outside" (also when τ is unknown) |
| `beta` | the recoil's speed: leaving the target for a plunger, where it is made for DSAM |
| `recommended_distances_um`, `sensitive_um` | log-spaced over the range where a measured unshifted share gives τ within twice the best precision (0.18 to 5.1 decay lengths γβcτ) |
| `separation` | per γ-ray detector: its angle, the shifted − unshifted separation, its FWHM, and whether the peaks are `resolved` (at least two FWHM apart) |
| `attenuation` | F(τ) on the setup's backing (DSAM) |
| `backing_advice` | DSAM: the usual backings ranked for this level (below); None for a plunger level |
| `hours_per_distance` | the beam time at each distance for the wanted precision on τ (default 10%) |

and for the whole measurement `time_to_stop_ps` (the recoil's time to stop in the backing, or in the plunger's
stopper, where decays make a DSAM-like tail under the unshifted peak), `stopper` (the beam in the plunger's foils,
below), `problems` and `warnings`.

**Which backing for DSAM.** A measured attenuation F with uncertainty σ_F gives σ_τ/τ = σ_F / |dF/d ln τ|, so the
best backing for a level is the one on which F changes fastest with τ: for stopping proportional to speed,
F = α/(α + τ) and |dF/d ln τ| = ατ/(α + τ)² is largest, ¼, when the slowing-down time α matches τ and F = ½. For
each level DSAM would measure (every level with a τ in a DSAM setup, and the levels whose τ suits DSAM in a plunger
setup) `backing_advice` lists Pt, Au, W, Ta, Pb, Ag, Cu, Ni, Fe and C (and the setup's own backing material if it is
none of them), best first, each behind the setup's target and with its slowing-down kernel (ZBL scattering, 250
recoils, the same random numbers for every backing): `material`, `density_g_cm3`, `attenuation` F(τ),
`sensitivity` |dF/d ln τ| (`lifetimes.dsam_sensitivity`), `relative_sensitivity` (its share of the best's), `good`
(at least 80% of the best's: σ_τ at most 25% larger), `slowing_time_ps` α (`lifetimes.slowing_time`: ∫β cos θ dt
over its value at the start), `time_to_stop_ps`, `thickness_to_stop_mg_cm2` and `_um` (the typical recoil's range in
it) and `thickness_mg_cm2` and `_um` (1.5 times that, rounded up: room for the recoils' spread and the stopping's
5–10%), and `setup` (the setup's backing material). For decays mostly inside the target (a few tens of fs after
fusion) every backing does about as well; heavy and medium backings often lie within a few per cent of each other,
so read the ranking with the `good` mark. The setup's backing below 80% of the best is warned about; a backing that
does not stop the recoil is already a problem of the setup.

**The beam in the plunger's foils.** The beam that passes the target flies on into the plunger's stopper (for a
differential plunger the degrader and the second foil), in normal and inverse kinematics alike. Per foil `stopper`
gives the beam's energy entering and leaving it (`beam_in_mev`, `beam_out_mev`, 0 when it stops there),
`deposited_mev`, `power_w` = particle current × deposited energy against `power_limit_w` = 0.5 W (a stretched foil
of a few mg/cm², cooled only by radiation, sheds about half a watt from a 1 cm² beam spot at about 900 K; above it
the foil heats, sags and the distance drifts), and the foil's own Coulomb excitation by the beam: `lines` of its
isotopes (those with at least 10% of its atoms), each with `yield_per_beam` (first-order Coulex cross sections of the
level scheme, integrated over the foil from the entry energy down to where the beam leaves or half the entry
energy) and `rate_per_s` (full-energy, all γ-ray detectors, singles); `missing` lists isotopes without level data
(the setup's scheme for that nuclide, else the local ENSDF copy): physim says so rather than guess. A foil over the
power limit and a foil line within 3 FWHM of a studied line's unshifted peak or its shifted peak at any detector are
warnings, in the plan and the planner's banner; no level data is a note. For ⁹⁸Ru on Cu, for instance, ⁶³Cu's 670 keV
line sits next to the 652 keV line shifted at 41.5°: gate on the particles or the feeding γ rays.

**How the numbers are worked out.** The recoil is one typical nucleus
(`Experiment.lifetime_recoil`): made in the middle of the target along the beam (Coulex in a
head-on collision; fusion at the compound nucleus's speed), slowed through the rest of the target and the backing.
The events and runs have the full spread of speeds and angles; the plan has this one recoil. F(τ) is
`lifetimes.dsam_attenuation` on the slowing-down kernel of that recoil, with ZBL scattering.
The beam time per distance comes from the line's full-energy rate in all γ-ray detectors (Coulex: the planner's
particle–γ coincidences; fusion: the residue rate times the line's γ rays per residue times the efficiency) and the
binomial precision of the unshifted share at each distance, √(e^x − 1)/(x√N) for x = t/τ, combined over the
setup's distances (or the recommended ones). It is statistics only, for one level: background and the systematic
entries of the analysis's budget come on top.

**Warnings** the plan gives: a Coulex plunger without deorientation for a level that lives longer than 1 ps (the
recoil's alignment fades in vacuum, which biases the shifted/unshifted ratio at each angle); a level whose τ suits
the other method; fewer than three distances in the sensitive range; peaks that no detector resolves; levels with no
known τ; a DSAM backing far from the advice; a plunger foil taking more than 0.5 W of beam, or with Coulex lines near
the studied ones.

## Predicted spectra

```python
sp = p.lifetime_spectra()                  # every distance (DSAM: the line shape) and every ring of detectors
sp = p.lifetime_spectra(distance_um=25.0)  # one distance
```

`sp["spectra"]` holds one entry per distance and ring of γ-ray detectors (their distinct polar angles): the counts
per bin on `sp["edges_kev"]` (expected full-energy counts in the beam time per distance), where the shifted and
unshifted peaks are, and the unshifted share. `sp["model"]` says how they were made:

- **fusion-evaporation**: `residues.plunger_spectra` or `residues.dsam_lineshapes` (residues
  drawn from the evaporation, the cascade with its side-feeding bands), scaled by the residue rate and each ring's
  efficiency and folded with its resolution;
- **Coulomb excitation, a plunger**: analytic, fast enough for the app. The unshifted share is
  `lifetimes.unshifted` of the level's cascade at the typical recoil's flight time; the shifted peak is a Gaussian at
  its Doppler energy as wide as the resolution and the spread of angles over a crystal. Real recoils spread in speed
  and angle (a 0° particle detector also catches forward-scattered beam), which broadens and moves the shifted peak:
  take a run to see that;
- **Coulomb excitation, DSAM**: `lifetimes.lineshape` on the typical recoil's slowing-down kernel.

## Simulating the runs

A plunger measurement is one run per distance. Runs remember their distance, and the experiment folder groups them
into a distance series:

```python
p.create_experiment(root="~/.physim/experiments")
for d in p.experiment.lifetime.distances_um:
    p.start_run("beam", duration="2 h", distance_um=d)         # "fusion" for a fusion-evaporation setup
series = p.folder.distance_series(p.experiment)[0]             # runs.DistanceSeries: {distance: [runs]}
```

A DSAM measurement needs only runs of the setup (no distance). In the app, the Run step takes the series and the
Data step shows the spectra per distance.

## Analysing them

```python
a = p.analyse_lifetimes()                         # the studied level, this setup's series, the default method
a = p.analyse_lifetimes(level="4+", method="decay_curve")
a["tau_ps"], a["stat_ps"], a["budget"], a["b_down"], a["b_unit"], a["compare"]
```

What it does, as an experimenter would:

1. **Spectra per distance and ring** (`dataviews.distance_spectra`), with the Doppler shifts kept;
   the rings closer than about 75° to 90° are left out (too small a shift).
2. **The recoil's speed** from where the shifted peaks are, not from the plan.
3. **Shifted and unshifted intensities** per ring (`lifetime_analysis.peak_areas`): counted in windows over a
   straight background when the shift is at least 3 FWHM (the shifted line is broad and flat-topped, and the
   stopper's slowing tail fills the gap), else Gaussian fits. Coulomb excitation: the shifted counts are weighted
   by the speed each bin shows (the emitters' speeds spread widely), each ring's shifted/unshifted efficiency
   ratio is taken from the setup, and, with deorientation in the setup, the components are corrected for it (the
   alignment of the particles the detectors see, iterated with τ). A differential plunger: the fully shifted,
   reduced-shifted and unshifted peaks.
4. **τ**: `method="ddcm"` (default; the differential decay-curve method, with the strongest feeding line's
   intensities for the feeding correction), `"ddcm_three_foil"` (the default for a differential plunger), or
   `"decay_curve"` (one exponential). DSAM: `"lineshape"`, the line shape at each ring fitted on the transport
   kernel, τ the weighted mean of the rings.
5. **B(λ↓)** from τ with the scheme's branching and conversion (`lifetime_analysis.b_from_result`), and `compare`:
   τ and B from the level scheme (ENSDF, typed), the local ENSDF copy, the setup's Coulex B(E2) and shell-model
   matrix elements in the scheme.

The result has `points` per distance (the unshifted and shifted shares, their error and the DDCM's τ(x)), `curve`
(the fitted curve, or the fitted line shape), a `budget` (statistics, distance, velocity, feeding, stopping) and
`notes` that say what was assumed.

## The examples

Three example setups reproduce published measurements (their numbers and validation status are in the
[physics register](physics-register/lifetimes.md)):

- `coulex_ru98_plunger`: ⁹⁸Ru on ²⁴Mg at Cologne (Radeck et al., PRC 85, 014301 (2012)), a plunger after Coulomb
  excitation with deorientation in vacuum;
- `fusion_cr50_plunger`: ⁵⁰Cr from ³²S + ²⁴Mg with the CoCoDiff differential (three-foil) plunger (Beckers et al.,
  NIM A 1042, 167418 (2022));
- `fusion_sm143_dsam`: ¹⁴³Sm from ²⁴Mg + ¹²⁴Sn on a gold backing (Rajbanshi et al., PLB 782, 143 (2018)), DSAM.

```python
p = Planner.example("fusion_cr50_plunger")
p.lifetimes()["levels"][:3]       # the 2+, 4+ and 6+ of the yrast band, each with its method and distances
```

The fusion examples take their level schemes from a local ENSDF copy; without one, add `[levels.50Cr]` (or the
residue you study) to the setup.

## From an LLM, through the MCP server

The MCP server (`physim mcp`; see [Shell-model calculations](shell-model.md) for adding it to Claude Code) has four
lifetime tools. Quantities are text with units, and the spectra's bins are left out.

| Tool | Does |
|---|---|
| `plan_lifetime_measurement(level?, method?, precision?)` | the plan above, per level (with the backing advice: the first five backings and the setup's); the beam in the plunger's foils; `method` plans as if the setup used that method, changing nothing |
| `set_plunger(distances?, offset?, kind?, material?, thickness?, density?, second_foil?, deorientation?, method?, remove?)` | sets the `[lifetime]` section through the draft methods; returns the summary with problems and warnings |
| `lifetime_spectra(distance?, level?)` | where the peaks are, the line's counts and the unshifted share per distance and angle |
| `analyse_lifetimes(level?, method?)` | τ with its budget, the points, B(λ↓) and the comparison |

## In the report

The beam-time report (`python -m physim.nuclear.report`) has a **Lifetime measurement** section for a setup with
`[lifetime]`: the measurement as set up, the plan per level, its warnings, the analysis of the runs when the
experiment has some, and the caveats: the limits this setup reaches and the lifetime checks of the validation suite
that fail. The models table lists every check of the lifetime capabilities, picked by capability, so a check added
to the suite appears there without a change to the report.

## Limits

- **Deorientation in vacuum** follows the static model (G_k(t) = α_k + (1 − α_k)e^(−Γ_k t), α_k and Γ_k given);
  dynamic deorientation and an a-priori calculation of α_k and Γ_k are left out. The planner's analysis corrects
  Coulomb excitation's first excited state for it; the correction is good to about 2% of τ (the alignment is
  averaged over a ring, with one flight time per distance at the mean speed), which is not in the budget. Other
  levels are not corrected, and the result's notes say so.
- **Stopping powers** set the DSAM line shapes and F(τ): known to about 5–10% for heavy ions; the line-shape fit's
  budget scales the kernel's time axis by ±7%. Two DSAM checks are known failures, with their tolerances kept: the
  APCAD comparison (physim's first-order shift larger) and the AhKin comparison (physim's line wider, from the ZBL
  scattering's sideways velocity). The report quotes the suite's own numbers.
- **One typical recoil in the plan.** The plan's speeds and separations are those of one recoil; the analysis takes
  the speed from the shifted peaks.
- **The planner's choices against the papers'.** The recommended distances fall inside the papers' ranges to within
  a factor 1.5 at each end (⁹⁸Ru 19–560 µm against 3–450 µm, ⁵⁰Cr 27–780 µm against 17–750 µm), and the paper's Au
  backing for ¹⁴³Sm is within 80% of the best backing's sensitivity for every level (validation checks
  `rdds_planner_distances`, `dsam_planner_backing`).
- **The beam in the plunger's foils** is a warning, not part of the simulation: the runs do not add the foils' own
  γ rays, the power limit is a rule of thumb for radiation-cooled foils, and the foils' Coulex is first order, below
  half the entry energy left out.
- **Line shapes near 90°** leave out the recoil's sideways velocity; time dilation is left out of DSAM (below 0.2%).
- Not covered: fast electronic timing (LaBr₃), intermediate-energy plungers beyond the relativistic Doppler formula,
  lifetimes from Coulex reorientation.
