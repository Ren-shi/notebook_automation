# Results and comparison

physim keeps every result as a row: one value of one observable of one nucleus, by one method, with its
uncertainty and where it came from. The rows of a nucleus are read together, across its experiments, so that a
value is always seen beside what else is known about it. This page says what a row holds, where the rows are, how
they are compared, and what was left out on purpose. The screens that show them are described in
{doc}`planner-app`: the nuclei, the page of a nucleus, and *Beside what is known* on an experiment's Analysis step.

No number physim computes was changed to make this: the rows store the numbers the analyses already gave, and the
conversions call the functions the code already had.

## A row

A row holds:

- **the nucleus** (`62Ni`);
- **the observable** (`B(E2)`, `tau`, `Qs`, `B(M1)`, `mu`, `E_level`; see [the registry](#the-registry-of-observables));
- **the level or transition** it belongs to (see [Names](#names-of-levels-and-transitions));
- **the value** with its unit and, for a B value, its direction (↑ or ↓ is always stated, never implied);
- **the uncertainty**, absolute, in one of three forms: one total, statistical and systematic, or plus and minus;
  and the uncertainty budget where the analysis has one;
- **the method**, from a fixed list, with a line of detail:

  | Method | |
  |---|---|
  | Coulex from the yield (first order) | physim's Coulomb-excitation analysis of a run |
  | Coulex fitted to all orders | the matrix elements fitted to the yields |
  | RDDS (decay curve, DDCM) | a plunger |
  | DSAM (line shape) | line shapes in a backed target |
  | Shell model | a calculation |
  | Evaluation | ENSDF's adopted values |
  | Converted from another row | a B(E2) from a τ, or the other way |
  | NRF, (n,n′γ), Electron scattering, Inelastic scattering, Transfer reaction, Moment measurement, Other | for values from the literature or your own data points: physim has no analysis that gives them |

- **the evidence** (below);
- **where it came from:** a run of an experiment with the analysis settings; a reference (authors, year, DOI, where
  in the paper); a shell-model job (interaction, space, effective charges); your own measurement (a date and a
  note); or the row it was converted from;
- **physim's version and commit**, and the time it was written;
- for a projected row, **the value the simulation was given**;
- **the rows it depends on** (a value converted from another row depends on it) and **the row it supersedes**.

### The two kinds of evidence

Rows are in two groups, which are never mixed.

- **Independent:** *measured* (a published value, or your own), *evaluated* (ENSDF's adopted values), *theory* (the
  shell model) and *outside program* (GOSIA and the like on the same inputs).
- **Projected:** a simulated run analysed back.

A projected value comes from a run that physim simulated with a value it was given, then analysed as a real run
would be. Getting that value back shows that the chain from the spectra to the result works (the closure), and the
uncertainty shows what precision the setup would reach. It does not measure the nucleus. So a projected value is
never set against a measured or evaluated value as if it tested it. It is closed against the value its simulation
was given, and its precision is set beside how well the value is already known:

> This setup reaches ±7.1 %; the value is known to ±1.5 % (the reference), and to ±4.1 % by the best single
> measurement.

### Where the rows are

Rows are kept in JSON-lines files, one row per line:

- **`results.jsonl` in an experiment's folder:** the rows the experiment's analyses produced. The folder stays
  complete on its own, to zip or move.
- **`points.jsonl` in the nucleus's folder under the nuclei folder** (`~/.physim/nuclei/62Ni/points.jsonl`; the
  environment variable `PHYSIM_NUCLEI` moves the nuclei folder, which is otherwise `nuclei` beside the experiments
  folder): the rows that belong to no experiment. These are ENSDF's adopted values, your data points, and the
  shell-model results.

There is no central database: the page of a nucleus reads its own file and the files of every experiment.

What writes rows:

| What you do | Rows | Where |
|---|---|---|
| Analyse a Coulomb-excitation run | B(E2)↑ of the excited transition, projected | the experiment's `results.jsonl` |
| Analyse a lifetime measurement | τ of the level, and the B(E2)↓ derived from it, projected | the experiment's `results.jsonl` |
| Fit the matrix elements to a run's yields (all orders) | B(E2)↑ and Q_s, projected | the experiment's `results.jsonl` |
| Read or apply a shell-model result | level energies, B(E2)↑, B(M1)↑, Q_s and μ, theory | the nucleus's `points.jsonl` |
| Open the page of a nucleus the first time | ENSDF's adopted level energies, τ, Q_s, μ, B(E2) and B(M1), evaluated | the nucleus's `points.jsonl` |
| Add, correct or retract a data point | measured, outside program or theory | the nucleus's `points.jsonl` |

An experiment from before the results existed gets its rows from its kept analyses when it is opened, or when the
nuclei are listed. A result that cannot be written is logged and never stops the analysis.

**Rows are appended, never rewritten.** A new analysis of the same thing supersedes the row that stood; an analysis
that gives exactly what the standing row says adds nothing, so pressing Analyse twice does not grow the file. A
data point is corrected by a new row that supersedes it, and withdrawn by a row that says why. Taking the evaluation
in again adds nothing, unless a newer ENSDF copy changed a value: then the new value supersedes the old one. The
superseded and withdrawn rows stay in the history.

ENSDF values without an uncertainty (limits, approximate values) are left out of the evaluated rows. A B(E2) is
derived from an ENSDF half-life only for a level whose γ rays the reader could all place.

## Names of levels and transitions

A level is named by its nucleus, J^π and ordinal: `62Ni:2+1` is the first 2⁺ state of ⁶²Ni. The name is anchored
to the level of the local ENSDF copy within 2 keV, so the ordinal does not depend on how far a level scheme was cut.
A level that ENSDF does not have (or any level, without a local ENSDF copy) keeps its energy and J^π instead,
`62Ni:2+@1173`, and is marked unanchored.

A transition is its upper level, its lower level and the multipolarity: `62Ni:2+1>0+1:E2`. The direction of a B
value is stated in the row, not in the name.

You write a level or transition as you would say it: `2+1`, `2+ at 1173 keV`, `2+1 > 0+1`. The pages write them as
B(E2; 0⁺₁ → 2⁺₁)↑ and τ(2⁺₁).

## The registry of observables

Each observable has one convention in which all its rows are compared. A value entered in another unit or
direction is converted for the comparison and keeps what it was given as ("given as 12.13 W.u.↓").

| Observable | Belongs to | Convention | Units accepted |
|---|---|---|---|
| B(E2) | a transition | B(E2)↑ in e²fm⁴, from the lower level to the upper one | e²fm⁴, e²b², W.u. |
| τ | a level | the level's mean life in ps (not the half-life) | fs, ps, ns, s |
| Q_s | a level | the spectroscopic quadrupole moment in e fm², with its sign | e fm², e b |
| B(M1) | a transition | B(M1)↑ in μ_N², from the lower level to the upper one | μ_N², W.u. |
| μ | a level | the magnetic dipole moment in nuclear magnetons μ_N, with its sign | μ_N |
| E | a level | the level's excitation energy in keV | keV, MeV |

**A value in Weisskopf units is always the ↓ value** (the decay's), and is marked so. B↑ and B↓ differ by the ratio
of the spin factors (2J + 1) of the two levels: for 0⁺ → 2⁺, B↑ = 5 × B↓. The page of a nucleus can write B(E2)↑ in
e²b² instead of e²fm⁴; that changes the display only. (One table does not follow this yet: the matrix elements of a
level scheme, on the Plan step and in Physics, give B↑ with B↑ over the Weisskopf unit beside it, unmarked.)

**What a value gives.** From the reference value of an observable the page works out B↓, the value in Weisskopf
units, the size of the matrix element, and τ from B(E2) (or B(E2) from τ) through the γ-ray energy, the branching and
the conversion coefficient. B(E2) and τ are converted into each other only for a level with one γ ray, where the
branching is 1. A converted value says what it was made from, depends on it, and is made again on reading: it is
never written to a file. The evaluation is not converted, nor is a value whose partner is already a row (a lifetime
analysis writes its B(E2) beside its τ).

**What it is used for.** Each observable lists what its value tells, with a reference: for B(E2), for example,
collectivity in Weisskopf units, the lifetime it fixes, the size of the deformation read as a rotor, a test of a
shell-model interaction and its effective charges, the next step up the band, and its use as an input to other
experiments.

## Comparison

The rows of one observable of one level or transition are compared with a **reference**: the evaluated value,
unless you pick another row. A projected row cannot be the reference.

Each independent row gets:

- the **difference** from the reference, in percent;
- the **pull**: the difference over the two uncertainties added in quadrature. Of an asymmetric uncertainty, the
  side that faces the other value counts;
- the **verdict**, by a rule fixed beforehand:

  | Pull (size) | Verdict |
  |---|---|
  | below 2 | agree |
  | 2 to 3 | tension |
  | above 3 | differ |

- **no uncertainty, no verdict:** a row without an uncertainty (a calculation, for example) gets its ratio to the
  reference and is not judged;
- a mark (◦) when it is **not independent of the reference**, because the two share a source: their agreement is
  no test.

Each projected row gets its **closure**: its value minus the value its simulation was given, over its own total
uncertainty, judged by the same rule (*closes*, *closes poorly*, *does not close*). Beside it stand the precision
the run reaches, the reference's precision, and that of the best single measurement.

Rows that hold now are compared. A superseded or withdrawn row is shown only in the history.

## Routes

A **route** is a way to reach an observable: a method, and the kind of experiment that plans it where physim has
one. For B(E2) the routes are Coulex from the yield, Coulex fitted to all orders, RDDS, DSAM (these four can be
planned in physim), inelastic scattering (a kind physim does not have yet), and the shell model, the evaluation,
NRF, (n,n′γ) and electron scattering (not experiments of the planner: their values come in as data points). The
other observables have their own lists.

Each route has a state for the level or transition:

| State | Meaning |
|---|---|
| done | an independent value by this method holds now |
| projected | only simulated runs by this method hold; nothing is measured by it yet |
| in physim | physim has this kind of experiment and it can be planned (**Plan it**) |
| later | a kind of experiment physim does not have yet |
| no experiment | the evaluation, a calculation, or a method known from the literature |

A route is also **still open** when physim can plan it and no run or measurement of your own has taken it, so a
published value marks a route done and still leaves it open to plan.

Each route says whether it **applies** to this nucleus, with the reason, or that this is not tested:

- whether the nucleus can be the target or the beam (stable, with its share of the natural element; or ENSDF's
  half-life);
- for RDDS and DSAM, whether the level's τ (evaluated, or the best measured) is in the range physim plans a plunger
  for (1 ps to 1 ns) or DSAM (10 fs to 1 ps);
- for the shell model, whether a calculation with an interaction physim has fits this computer now (half the RAM
  free); this tests the size only, not whether the valence space suits the nucleus;
- whether ENSDF has the nucleus;
- for Q_s, the spin (no Q_s for J = 0 or ½);
- for the yield analysis and electron scattering, whether the transition starts from the ground state.

A route with a value by it is never hidden as "does not apply": a published value by a method outside physim's
range shows as done, with "not tested" and the reason.

## Checks between observables

The page of a nucleus sets a value beside others of the same nucleus. Only independent rows take part; where a check
needs one value it takes the evaluated one, else the measured one with the smallest uncertainty.

- **E(4⁺₁)/E(2⁺₁):** 2.0 is a harmonic vibrator, 2.5 a γ-soft rotor, 3.33 a rigid axial rotor; below 2 the nucleus
  is near a closed shell. A ratio of energies alone: it shows neither the sign of the deformation nor how
  collective the transitions are.
- **B(E2; 4⁺₁ → 2⁺₁)/B(E2; 2⁺₁ → 0⁺₁), both ↓:** 10/7 = 1.43 is the rigid axial rotor and 2 the harmonic vibrator;
  below 1 the 4⁺₁ state does not continue the collective motion of the 2⁺₁ state. Read it with E(4⁺₁)/E(2⁺₁).
- **The rotor's Q_s(2⁺₁) from B(E2), beside the measured Q_s:** read as a rigid axial rotor, B(E2; 0⁺₁ → 2⁺₁)
  gives β₂, Q₀ and the size of Q_s(2⁺₁); the sign of the measured Q_s says prolate (negative) or oblate
  (positive). A measured Q_s of another size says the state is not that rotor. Agreement in size alone does not
  prove a rotor.
- **A calculation over the reference** (a shell-model value over the evaluated one), for each observable: 1 means
  the interaction, its valence space and its effective charges (or g factors) together give the value; it does not
  say which of them is off.

These are plain ratios. No uncertainty is worked out for them, so a check has no verdict: each value in it is shown
with its own uncertainty and source.

## What was left out on purpose

The results pages add no physics of their own. What would need a new formula, with its own validation, is backlog
item 76 (*Physics the results pages could add*):

- the effective charge a shell-model value would need to meet the measured one;
- the B(E2) expected from E(2⁺) (Grodzins/Raman-type systematics);
- sum rules;
- the uncertainty of a ratio, of a matrix element taken from a B value, and of the rotor's numbers (so the checks
  have no verdict, and the rotor's Q_s is not pulled against the measured one);
- B(M1) from a lifetime, and the E2/M1 mixing ratio;
- the g factor, the Schmidt values and Z/A beside a magnetic moment;
- half-life to mean life as a function of its own, and the branching of a level with several γ rays (so B(E2) and τ
  are converted only for a level with one γ ray);
- whether a target or beam of an unstable nucleus can be had, whether a valence space suits the nucleus, and
  whether Coulomb excitation's yield reaches the state.

Also left out, by decision: fetching literature values from NNDC or XUNDL (they come in as data points); averaging
measured points into an evaluation of physim's own (the page shows the spread, and ENSDF stays the evaluation); and
raw measured spectra.

## From Python

```python
from physim.nuclear import nucleus, results, compare

nucleus.snapshot_evaluation("62Ni")                     # ENSDF's adopted values as evaluated rows
nucleus.add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 906, "e2fm4", direction="up", uncertainty=37,
                       method="coulex-fit", reference={"authors": "Allmond et al.", "year": 2014,
                                                       "doi": "10.1103/PhysRevC.90.034309"})
view = nucleus.observable_view("62Ni", "B(E2)", "2+1 > 0+1")
[(e.label, e.difference_percent, e.pull, e.verdict) for e in view["comparison"].entries]
[(p.closure, p.reach_percent, p.known_percent) for p in view["comparison"].projected]
nucleus.import_points("points.csv")                     # many points from a file
```

`results.of_experiment(folder)` and `results.of_nuclide("62Ni")` read the rows; `results.current(rows)` drops the
superseded ones; `compare.compare(rows)` compares them.

## From an LLM, through the MCP server

The MCP server (`physim mcp`; see {doc}`shell-model` for adding it to Claude Code) has seven tools for the results.
They act on a nucleus and need no open experiment.

| Tool | Does |
|---|---|
| `list_nuclei()` | the nuclei physim holds something about, each with its experiments (studied or partner), runs and analysed runs, what is known per observable, and what is still open; and the experiments that belong to no nucleus, with the reason |
| `list_results(nucleus, observable?, key?, history?)` | the rows of a nucleus that hold now, each with its id, evidence, method, value in the observable's convention and source; `history` lists every row ever written with its status (holds, superseded, retracted) |
| `snapshot_evaluation(nucleus)` | takes ENSDF's adopted values in as evaluated rows; says how many were added, superseded and unchanged |
| `add_data_point(nucleus, observable, key, value, method, unit?, direction?, uncertainty? or stat and sys or plus and minus, evidence?, reference? or own_measurement?, detail?, corrects?)` | adds a published value or your own measurement; everything wrong is said at once; `corrects` names the point it replaces |
| `retract_data_point(nucleus, row, reason)` | withdraws a data point; nothing is deleted |
| `compare_results(nucleus, observable, key, reference?)` | every result for one observable of one level or transition: the independent ones against the reference with difference, pull and verdict; the projected ones apart with their closure and precision; the routes with their state and whether they apply; what the value gives; the checks; what it is used for |
| `import_data_points(path)` | reads data points from a JSON or CSV file, or from the literature file of a nucleus; says what was added, what was already there and what was refused |

Quantities follow the same conventions as on the page: B(E2) in e2fm4 (or e2b2, or W.u., which is always down), τ
in ps, Q_s in efm2, B(M1) in mu_N2, μ in mu_N, E in keV.

A short session ("how does the Allmond 2014 B(E2) of ⁶²Ni compare, and what would my setup add?"):

```text
snapshot_evaluation("62Ni")
  → reference "ENSDF adopted, NDS 204, 1 (2025)", added 456, superseded 0, unchanged 0
add_data_point("62Ni", "B(E2)", "2+1 > 0+1", 906, method="coulex-fit", unit="e2fm4", direction="up",
               uncertainty=37, reference={"authors": "Allmond et al.", "year": 2014,
                                          "doi": "10.1103/PhysRevC.90.034309"})
  → added "B(E2; 0+₁ → 2+₁)↑ = 906 ± 37 e2fm4↑", key "62Ni:2+1>0+1:E2", evidence "measured",
    source "Allmond et al. 2014"
compare_results("62Ni", "B(E2)", "2+1 > 0+1")
  → reference: "ENSDF adopted, NDS 204, 1 (2025): 884 ± 13 e2fm4↑"
    independent: Allmond et al. 2014, Coulex fitted to all orders, 906 ± 37 e2fm4↑,
                 "+2.5 %, pull +0.56: agree"
    projected:   62Ni-on-194Pt · 03-beam-24h, Coulex from the yield (first order), 919 ± 65 e2fm4↑,
                 simulation was given 884.105 e2fm4,
                 "+0.53 standard deviations from the value the simulation was given: agree",
                 setup reaches ±7.1 %, known to ±1.5 %, best single measurement ±4.1 %
    routes:      RDDS (decay curve, DDCM): in physim, applies: "τ = 2.08 ps (ENSDF adopted, NDS 204, 1 (2025))
                 is in the range physim plans a plunger for (1 ps to 1 ns)"; ...
```

The projected run is reported with its closure and the precision it reaches, never with a pull against the
reference. (The experiment's rows are there once its runs have been analysed, or once the app has opened it.)
