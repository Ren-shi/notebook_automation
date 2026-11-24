# The experiment planner app

The planner app is for planning beam time and keeping what is known about the nuclei you study. You describe an
experiment by pointing and clicking, and it shows what you will measure: geometry, kinematics, count rates and beam
time, energy loss, simulated spectra and Coulomb orbits. It analyses the simulated runs as you would analyse real
ones, sets every result beside what is already known about the nucleus, and writes the report. You don't need to
write any code.

It runs on your own computer, in your web browser. There are no accounts, and nothing is uploaded.

## Install and start it

**On Windows, without Python:** download `physim-planner-…-windows-x64-setup.exe` from the
[releases page](https://github.com/Ren-shi/notebook_automation/releases) and run it.

1. Windows may say *"Windows protected your PC"*, because the installer is not code-signed yet. Click **More info**,
   then **Run anyway**.
2. The installer needs no administrator rights. It puts everything, including its own copy of Python, in
   `%LOCALAPPDATA%\Programs\physim planner`, and adds **physim planner** to the Start menu (and, if you tick the
   box, to the desktop).
3. Click **physim planner**. The planner opens in its own window after a few seconds (on Edge WebView2, which
   Windows 10 and 11 bring). Closing the window stops it. Clicking the shortcut again while it runs opens the
   running planner.
4. To remove it: *Settings → Apps → physim planner → Uninstall*. To update, run a newer installer over the old one.

If the window does not open, the log is in `%LOCALAPPDATA%\physim\planner.log`. The browser is still there:
`physim app` (or the shortcut's command with `--browser`) opens the planner in a browser tab, for a planner on a
server that others reach from their own machines.

**If you have never used Python (macOS, Linux, or Windows without the installer):**

1. Install Python 3.11 or newer from <https://www.python.org/downloads/>. On Windows, tick **"Add python.exe to
   PATH"** on the first screen of the installer.
2. Open a terminal: on Windows, *Terminal* or *Command Prompt*; on macOS, *Terminal*.
3. Install physim with the app:

   ```bash
   python -m pip install "physim-engine[app,root]"
   ```

4. Start the planner:

   ```bash
   physim app
   ```

   Your browser opens at <http://localhost:8080>. Leave the terminal open while you use the app, and close it (or
   press Ctrl+C in it) to stop.

**If you already use Python:** `pip install "physim-engine[app,root]"`, then `physim app` (or `python -m physim app`).
`physim app --example oxygen_on_lead_array` opens an example instead of the nuclei; `--port` picks another port.

`physim app --window` opens the planner in its own window instead of a browser tab (`pip install
"physim-engine[app,window]"` adds pywebview, 5 MB; without it the browser opens and the command says so).

`physim app` listens only on this computer (127.0.0.1). `--host 0.0.0.0` serves the local network too, for
example to show the planner on a lab PC; anyone on that network can then open it.

## The nuclei: the first screen

The planner opens on the nuclei, under the title **What are you studying?** Each nucleus holds its experiments and
everything known about it, by method. The most recent experiment is open behind this screen: **Back to the
experiment** returns to it.

Each nucleus that an experiment studies, or that has results of its own, has a card:

- **The head line:** how many experiments, how many analysed runs, how many data points, and whether ENSDF's
  adopted values have been taken in (*evaluation: ENSDF adopted, NDS …*, or *evaluation not taken yet*). Click the
  nucleus's name to open its page. **New experiment for ⁶²Ni** starts an experiment for it (see
  [Starting an experiment](#starting-an-experiment)).
- **What is known and projected:** the first few observables, each with the evaluated value, how many measured
  values there are, and what the simulated runs reach (*projected: reaches ±7.1 % (2 simulated runs · …)*). Click
  one to open the nucleus's page on it.
- **Still open:** the routes physim can plan that no run of your own has taken yet, the runs not analysed, and
  whether the evaluation is still to be taken in.
- **Experiments:** a table with each experiment's kind, runs and latest result. A row opens the experiment.
  Experiments with the same name show their folder too.

Below the cards:

- *Only as the partner of an experiment:* the nuclei that appear only as the other nucleus of a collision (the
  target of a Coulomb excitation of the beam, for example), each a link to its page.
- **Experiments tied to no nucleus**, with the reason: elastic scattering, and fusion-evaporation without a studied
  channel.

**Another nucleus** and **Open its page** open the page of any nucleus, also one with no experiment yet (it starts
with ENSDF's adopted values). Write it as `58Ni` or `Ni-58`. The page of a nucleus also has an address of its own:
`/?nucleus=62Ni`.

**Start from a kind of experiment** opens the kinds (below). The very first time, with no experiment and no
nucleus yet, the planner opens on the kinds.

An experiment belongs to the nucleus it studies. In Coulomb excitation that is the nucleus that is excited; the
other one is its partner. A fusion-evaporation experiment belongs to the residue of the channel its setup marks as
studied (`[run] channel`); without that mark it belongs to no nucleus. Elastic scattering belongs to none.

## The page of a nucleus

The page of a nucleus (from its card, or `/?nucleus=62Ni`) holds everything known and projected about it, by
method. **Nuclei** at the top goes back to the nuclei. The page reads ENSDF's adopted values in the first time it
is opened (from the local ENSDF copy; the line under the head says which data set, how many values and when).

From the top down:

- **The head:** Z, N, how many experiments, analysed runs and data points. **Add a data point** scrolls to the form;
  **New experiment for ⁶²Ni** starts an experiment.
- **The read-outs:** the headline observables (for ⁶²Ni, B(E2; 0⁺₁ → 2⁺₁)↑ and τ(2⁺₁)) with their evaluated value
  and what the simulated runs reach, and how many things are still open.
- **Observables:** one row per level or transition, with the columns *Evaluated*, *Measured*, *Theory*, *Projected
  (simulated runs)* and *Still open* (and *Outside program* when there is such a value). Click a row to see every
  value of it below. The rows that hold only the evaluation are folded under *… more hold the evaluation only*.
  The switch **B(E2)↑ in e²fm⁴ / B(E2)↑ in e²b²** changes how a B(E2) is written on the page; nothing is changed in
  the stored values.
- **Every value, side by side**, for the chosen row, against a reference (the evaluated value unless you pick
  another row in **against**). Two tables, never mixed:
  - *Independent:* the evaluated, measured, theory and outside-program values. Each has its source, its method in
    words, its value, how precisely it is known, and its difference from the reference, its pull and its verdict
    (✓ agree, △ tension, ✗ differ). A value converted from another quantity says so (*DSAM (line shape), from τ =
    2.01 ± 0.12 ps*). A value that shares a source with the reference is marked ◦: their agreement is no test.
  - *Projected:* the simulated runs analysed back. Each has its statistical and systematic uncertainty, the value
    its simulation was given, the closure (✓ closes, △ closes poorly, ✗ does not close) and the precision it
    reaches. A projected value is never set against the reference.
- **The values:** a plot of every value with its uncertainty, the reference's uncertainty as a band, and the value
  the simulations were given.
- **What the projected runs say:** in one sentence, the precision the setup reaches (and how much of it is
  systematic), how well the value is already known, and how well the best single measurement knows it.
- **Routes to …:** every way to reach the observable, each with its state: *done* (an independent value by this
  method is here), *projected only* (only simulated runs), *in physim* (an experiment of that kind can be planned),
  *later* (a kind of experiment physim does not have yet) or *does not apply*. Each route says whether it applies to
  this nucleus and why: whether the nucleus can be a target or a beam, whether the level's lifetime is in the range of
  a plunger or of DSAM, whether a shell-model calculation fits this computer. **Plan it for ⁶²Ni** (or **Plan another
  for ⁶²Ni** where a run of your own exists) opens the form for that kind of experiment. The routes that are not
  experiments of the planner (the evaluation, the shell model, methods known from the literature) are listed in one
  line under the cards.
- **What this value gives:** what follows from the reference value (B↓, Weisskopf units, the matrix element, τ from
  B(E2) or B(E2) from τ), the **checks against other observables** (for example E(4⁺₁)/E(2⁺₁) and the rotor's Q_s
  from B(E2) beside the measured Q_s), each with how to read it, and **what it is used for**, with references.
- **Experiments:** every experiment with this nucleus; a row opens it.
- **Add a data point**, and the **History** of every row ever written for the observable, newest first, with what
  superseded or withdrew it.

{doc}`results` explains the rows, the comparison, the routes and the checks in full.

### Data points

A data point is a value you bring in: a published value, or a measurement of your own. It is entered as a value
with its uncertainty; physim keeps no raw spectra. The form **Add a data point** asks for:

- the **Observable** and the **Transition** or **Level**, as you would say it: `2+1 > 0+1`, `2+1`, or `2+ at 1173
  keV`;
- the **Value**, its **Unit** and, for a B value, the **Direction** (*↑ excitation (lower → upper)* or *↓ decay*; a
  value in W.u. is always ↓);
- the **Uncertainty**, as **One total**, **Stat and sys**, or **Plus and minus**;
- the **Evidence** (*Measured*, *Outside program* or *Theory*) and the **Method**, with an optional line of detail;
- **Where it is from:** **A reference** (authors, year, DOI, where in the paper) or **Own measurement** (date and a
  note on where and how).

**Add the point** adds it. If something is wrong, the form lists everything at once under *Not a data point yet*.
A point is never changed: the ⋮ menu beside it offers **Correct…**, which adds a new row that supersedes it, and
**Retract…**, which withdraws it with a reason. Both stay in the history.

**Import points from a file** reads many points at once: from a JSON or CSV file, or from the literature file of a
nucleus. It says what was added, what was already there and what was refused.

### The report of a nucleus

**Report for this nucleus** builds a report of the nucleus across its experiments: what is known about it by every
route, and its experiments. The report of one experiment is on that experiment's Report step.

## Starting an experiment

There are two ways in.

**From a nucleus.** **New experiment for ⁶²Ni** on its card or page opens the kinds of experiment with a note
*An experiment for ⁶²Ni: pick its kind*, and **Back to ⁶²Ni**. **Plan it for ⁶²Ni** beside a route goes straight
to the form of that route's kind. Where the form has a field for the nucleus that is excited (Coulomb excitation
and lifetimes), the nucleus is filled in: as the beam, or as the target where it is the kind's own usual target.
The other nucleus, the energy and the thickness stay the kind's usual ones, and the form says so: set them for
your nucleus (the Setup step says whether the energy is safe). Nothing else is computed for you.

**From a kind.** **Start from a kind of experiment** shows the kinds under **What experiment are you planning?**:

| Kind | What the planner gives |
|---|---|
| Coulomb excitation | B(E2)↑, Qs, excitation probability, γ yields per angle, Doppler correction |
| Elastic (Rutherford) scattering | dσ/dΩ, rates per strip, energy loss, safe distance, beam time |
| Fusion-evaporation | channel cross sections, residue yields, γ-line rates, channel selection |
| Lifetimes (RDDS, DSAM) | τ per level, plunger distances, peak separation, F(τ), beam time per distance, decay curves, B(E2) |

*Transfer reactions* and *Inelastic scattering* are shown greyed: they come later.

Pick a kind and its form opens:

- **Beam and target:** beam, energy, target, thickness. For Coulomb excitation and lifetimes, **Which nucleus is
  excited** (the target or the beam). For fusion-evaporation, a line on the compound nucleus (its excitation, the
  highest spin that fuses, the beam against the Bass barrier, the likely xn channel) and the **Exit channels**:
  **Fill from physim** (physim's own fusion and evaporation, an estimate) or **Enter by hand**.
- **Detectors from:** a template. For Coulomb excitation: *Coulomb excitation of the beam (62Ni on 194Pt) with
  silicon at backward angles*, the same *and four clovers*, *Coulomb excitation of a target nucleus with a CD and
  clovers*, or *Blank: one detector, chosen from the kinematics*. Elastic scattering has *Elastic scattering with
  silicon pads* and *with an array of strip detectors*; fusion-evaporation a γ-ray array with or without
  ancillaries; lifetimes *Plunger after Coulomb excitation*, *Plunger after fusion-evaporation* and *Line shapes
  (DSAM) after fusion-evaporation, backed target*. A template brings its detectors; the beam and target are the
  ones you give.
- **Or start from a worked example:** the examples of that kind, opened as they are.
- **Name** (made from the beam and target if left empty).

**Make the experiment** makes it and opens its Setup step. For Coulomb excitation the first 2⁺ state is taken from
the local ENSDF copy; without one the setup says what to enter. **Change the kind** goes back to the kinds, and
**The nuclei** to the first screen. Below the kinds, **Or open an experiment** lists every experiment, and **Import
a setup file…** starts from a `.toml` setup (see {doc}`nuclear-setup`).

## Experiments

An experiment has a name, its current setup and the runs taken with it. It lives in a folder under
`~/.physim/experiments` (the environment variable `PHYSIM_EXPERIMENTS` moves it), where the setup is saved every
time you press **Apply**. Its results go to `results.jsonl` in the same folder (see {doc}`results`).

At the top of the rail on the left:

- the nucleus the experiment studies, a link to its page (nothing for an experiment that studies none);
- the experiment's name, which is a menu:
  - **Open an experiment…** and **New experiment…** both open the nuclei screen: open an experiment from its table,
    or start one with **New experiment for …** or **Start from a kind of experiment**;
  - **Import a setup file…** replaces the current setup with a `.toml` file; **Export the setup…** downloads it
    (see {doc}`nuclear-setup`);
- the beam, its energy, the target and its thickness, and the kind of experiment (*Coulex · beam*, *Lifetimes ·
  plunger after Coulex*, *Fusion-evaporation*, *Elastic*).

`physim app --example coulex_ni58` opens an example as an unsaved setup instead.

## The steps

The rail on the left is the experiment in the order it is done, from the top down: **1 Setup → 2 Plan → 3 Run →
4 Data → 5 Analysis → 6 Report**. The step on screen has its number on the accent colour.

- **Each step has a line of status under its name:** the detectors, the beam time the counts need, the runs
  taken, the current run, the result of the analysed run (B(E2) for Coulomb excitation, τ for a lifetime
  measurement), what the report holds.
- **A green tick in place of the number marks a step that has its result:** a setup without errors, a beam time
  that is enough, a run, a run's data, an analysis of the current run. A step that needs a run first says *needs a
  run* and is greyed.
- **A status in the warning colour needs attention:** warnings in the setup, a beam time too short, a run taken
  with a setup that has changed since.
- **Each step ends with one button that leads on:** *Next: check the plan*, *Next: take a run*, *Next: look at
  run 2's data*, *Next: analyse run 2*, *Next: the report*. Where there is nothing to go on to yet, the line
  says what is missing. For fusion-evaporation and elastic scattering, whose runs physim does not analyse back yet,
  the Data step leads to the report.
- The rail stays in place as you scroll, and you can choose any step at any time.

The steps:

- **Setup:** the experiment to scale; click a detector for its numbers and drag it to move it (see
  [The scene](#the-scene)). **Edit the setup** opens the setup as decisions (see
  [The setup as a list of decisions](#the-setup-as-a-list-of-decisions)).
- **Plan:** what you check before asking for beam, answers first. It opens with four numbers: the highest rate,
  the longest beam time for the counts you asked for, the head-on closest approach, and how many warnings the setup
  has. Then the beam time your counts need against the beam time planned; each particle detector's rate, busy
  fraction, share of the excitations and safe rings; each γ-ray detector's efficiency and coincidence rate; the
  particle × γ coincidence rates. The **?** beside a number shows how it is worked out. Below the answers, **The
  checks** open one at a time: *Kinematics: where the particles go*, *Energy loss in the target and the dead
  layers*, *The orbit and the safe distance*, *The excitation, the particle–γ correlation and the Doppler shifts*
  (with the level schemes), *γ-ray efficiency curves*, *Rates per strip, and the beam-time sweep*, and *Predicted
  spectra: what each detector will record*. The next run keeps these numbers, so the Analysis can put them beside
  what the run measured. A fusion-evaporation experiment shows the channel studied, the fusion, the residues (click
  a row to study that channel), the studied residue's band, the channel selection and what the compound nucleus
  boils off; a lifetime measurement adds its own panel (see [Lifetime measurements](#lifetime-measurements)).
- **Run:** **Run** takes the setup for a duration (its beam time by default). Kinds of run:
  - a beam run;
  - a source run: beam off, a calibration source at a position;
  - an alignment check: the target assumed off by an amount you give.

  Before you start, the step says how much is simulated event by event (10 min by default, up to an hour) and what
  that costs. While the run goes:
  - a clock and a bar in beam time;
  - the counts, rate and busy fraction of each detector, the singles and coincidences of each crystal, and the
    particle × γ matrix filling in;
  - a bar per detector for the counts you want;
  - the last particles counted.

  **Stop** keeps what is accumulated and **Extend** adds beam time to the same run. **Watch events** draws a
  sample of the run's events as tracks in the scene. The run list sets which run is current for the Data and
  Analysis steps. A run taken before the other nucleus of the collision was simulated is marked *one
  nucleus only*; **Take again with both nuclei** makes a new run with its setup, seed and beam time that holds
  both, and keeps the old one. **Not simulated yet** lists what the runs leave out. See also
  [Experiments and runs](#experiments-and-runs-from-python).
- **Data:** every spectrum of the current run at once, one panel per detector (per crystal on request). Without
  a run it says so and points to the Run step and to the Plan's predicted spectra. Click a
  panel's corner to enlarge it and download it as CSV. Above the panels you choose what is done to them:
  - the Doppler correction (off, for the projectile, for the recoil), with a line on why it is right or wrong for
    this setup;
  - a particle gate;
  - singles or coincidences, with the random coincidences shown, subtracted, or left out;
  - add-back and suppression on or off;
  - the binning, and a view window for the γ-ray panels (the spectra always hold the whole range of the run,
    up to the highest line recorded; the window only zooms);
  - the real part of the run, or the whole run scaled.

  **Gates** are named and kept with the experiment, and the analysis uses them by name. Make one from a detector,
  a ring range, the elastic or excited group and an energy window, or drag a box on the **energy against ring**
  view, where each group's kinematic line is drawn. The **γ-ray energy against crystal** view, raw or corrected,
  shows a misplaced target as lines out of step. **Compare runs** overlays the same spectrum from two runs and
  names what differs between their setups. **Export** writes a ROOT file (events, spectra, the gates as cuts) or
  the run's `.npz`.
- **Analysis:** a row of buttons at the top (*On this step*) jumps to its parts: *Lifetimes* (for a lifetime
  measurement), *Predicted against measured*, *B(E2) from the peak*, *Target alignment*, *All orders*. Without a
  run, the analysis and the alignment check work on a fresh simulated sample and the step says so. First,
  **predicted against measured**: what the Plan said before the run (kept with it) beside what the run measured.
  That covers the rates, the coincidences in the peak, the beam time for your counts, and the γ-ray efficiency
  (measured if the experiment has a source run). Each row says whether the two agree within the run's statistics.
  Then the analysis of the current run, from the γ-ray peaks to B(E2), with a gate from the Data step if you pick
  one. **Analyse** shows it in three parts, and the analysis kept with a run is shown as soon as you open the step:
  1. *Find the peaks:* the spectrum corrected for the projectile and the one corrected for the recoil, each over
     its whole range, with the peaks the search found marked and listed (energy, width, area, significance).
     A peak on a transition of the level scheme is named. Each nucleus's lines are sharp in its own spectrum; in
     the other they are smeared, and the table says so.
  2. *Integrate the peaks:* the peak of each nucleus with its fit (a Gaussian on a line, with tails where the
     peak has them) and the ±3σ window shaded, and under it the numbers: the counts over the background, the
     random coincidences subtracted, the area, the efficiency, the angular correlation, the γ rays emitted, and the
     yield per particle. For the other nucleus (¹⁹⁴Pt in the opening template) the yield its level scheme gives is
     beside the measured one.
  3. *From the area to B(E2):* the chain in numbers and the result. **Normalisation** chooses between the
     particles in the same rings and the other nucleus's peak. When the events follow a level scheme, the last
     step solves its coupled equations for B(E2) (the feeding of the state and a line the peak blends with are
     counted), and says what first order would have read and how much the excitation saturates.

  The result shows B(E2)↑, its uncertainty (statistical and systematic), the value given to the simulation with
  the closure against it, and the beam time for the precision you want. A value in Weisskopf units is the decay's,
  marked ↓. It says the statistical uncertainty comes from the run's own events, and what the whole run would give.

  Under the result, **Beside what is known** sets it beside the nucleus's other results, in short:
  - *PROJECTED*: the result in the comparison's convention, from a simulated run, so it is the precision the setup
    would reach and not a measurement;
  - *Reference*: the evaluated value with its source (or none, if no evaluated value is entered);
  - *Closure*: whether the analysis gives back the value the simulation was given;
  - *Reach*: *This setup reaches ±x %; the value is known to ±y %*, and how well the best single measurement knows
    it;
  - *Other routes*: how many routes hold a value, and which are still open;
  - **Compare on the ⁶²Ni page →**, a link to the nucleus's page on this transition.

  A projected result is never pulled against the reference. A result with no row in the results ledger says why.
  The alignment check reads an alignment run's offset. The all-orders solution can fit the matrix element to the
  run's yields, and writes the GOSIA input file.

  For fusion-evaporation and elastic scattering the step says *No analysis for … yet*: what an analysis of a run
  would need that physim does not have, with **The plan** and **The report**, and the predicted-against-measured
  table.
- **Report:** **the report** of the experiment, one document that leads with the results:
  1. *The result*: the value with its statistical and systematic uncertainties, its method and its runs, and plainly
     whether it is projected from a simulated run. A lifetime measurement leads with τ and the B(E2)↓ that follows.
  2. *Beside the other routes*: the independent values against the reference, with difference, pull and verdict;
     apart from them, the projected values with their closure and "this setup reaches ±x %; the value is known to
     ±y %".
  3. *Comparisons still to run*: the routes with no value yet for this nucleus and transition, the runs of the
     experiment not analysed, and the independent evidence that is missing.
  4. *What limits it*: the uncertainty budget, and the limits that bear on this value.
  5. *What it means.*

  Then the appendices A to H: the setup; the plan; the runs and their spectra; the gates; the Doppler correction;
  the analysis; every number, explained; the data, the method and what is left out. An elastic or
  fusion-evaporation experiment has no result of this kind: its report opens with what the setup gives (rates;
  the fusion cross section and the residues), then the appendices. Every B(E2) in the report says ↑ or ↓, and a
  value in Weisskopf units is always ↓. **Show the report** puts it on the page; **Save as PDF** and **Save as
  HTML** save it; **Save everything (zip)** adds the setup file, the plan's tables and figures and the experiment
  folder with every run; **Save the experiment folder** zips the folder alone. The PDF is printed by Edge, Chrome
  or Chromium; if none is installed the app says so, and the HTML page prints to the same PDF from any browser.
  Building a report writes nothing into the experiment's folder.

**Physics**, at the foot of the rail, is the reference behind the numbers; it is not a step. It follows the order
of the physics register: nuclear data and level schemes, kinematics, energy loss, the Rutherford orbit, the
detectors and the γ-ray response, rates and events, and Coulomb excitation. Each section links to its register page
and names the tests behind it. **See the physics** in a number's explanation, or the flask on a setup card, opens
the section it belongs to.

Each result opens with **How to read this**, using your numbers, and ends with an **Explain** panel giving the
formula, the assumptions and where they stop being valid.

## Lifetime measurements

A plunger (RDDS, for about 1 ps to 1 ns) or line shapes in a backed target (DSAM, for about 10 fs to 1 ps), after
Coulomb excitation or fusion-evaporation. {doc}`lifetimes` has the physics and the Python interface; in the app:

- **Setup:** start from the kind *Lifetimes (RDDS, DSAM)*, or add the measurement to a Coulomb-excitation or
  fusion-evaporation setup: **Add** under *Lifetime measurement* in the setup window offers **Plunger (RDDS)** and
  **Line shapes (DSAM)**. In a fusion-evaporation experiment the panel **Measure the lifetimes** on the Plan step
  does the same. The *Lifetime measurement* card holds the distances, the stopper or degrader foils and the
  deorientation (or the DSAM backing).
- **Plan:** the panel *Lifetimes · Plunger (RDDS)* (or *Line shapes (DSAM)*): per level its method, τ and where τ
  comes from, the recommended distances, the separation of the shifted and unshifted peaks at each detector, F(τ)
  for DSAM, the beam time per distance, the beam in the plunger's foils, and the predicted spectra per distance.
- **Run:** a run per plunger distance; **Run each distance not taken** takes the ones still missing.
- **Data:** the runs grouped by distance, with the spectra per distance and ring.
- **Analysis:** the part *Lifetimes*, with the **Level** and the **Analysis** (for a plunger: *Differential decay
  curve (DDCM)*, *Decay curve fit*, or *DDCM, three foils* for a differential plunger). It shows τ, its total
  uncertainty, B(Eλ)↓ from τ, the decay curve, τ at each distance, the uncertainty budget, and **Beside what is
  known** for τ and for the B(E2) that follows from it, as under a Coulomb-excitation result. The analysis is kept
  with the experiment: when the app is opened again it is shown, marked *Kept from an earlier session* with its
  date, and is not run again until you press **Analyse again**. The rail says when runs were taken that the kept
  analysis did not read.

## The shell model

Beside the level scheme of the target and of the beam (in the check *The excitation, the particle–γ correlation and
the Doppler shifts* on the Plan step, and in Physics under nuclear data) is the panel **Shell model: calculate or
import the levels of …**. It is closed until you open it.

- **Calculate:** choose the **Interaction** (with its valence protons and neutrons and the M-scheme dimension), the
  **States** (`0+ ×3, 2+ ×3, 4+ ×2`), the **Transitions** (E2, M1) and the **Effective charges**. *Before running*
  says what the calculation needs (dimension, memory, disk, a rough time) against the limit of half the RAM free
  now; over that limit **Run KSHELL** stays disabled unless you tick the override. The run goes on in the
  background, with its progress.
- **Or import results run elsewhere:** a KSHELL output folder, one log, KSHELL's summary file, or a physim
  shell-model table.
- **The results:** the states beside the known levels (paired by J^π and order). **Apply to the level scheme** puts
  the values you tick into it. They fill gaps and never replace ENSDF, derived or your own values unless you tick
  **Replace existing values**.

Every result read or applied is also copied into the nucleus's own results as theory (level energies, B(E2), B(M1),
Q_s and μ), so it shows on the nucleus's page and survives the cleaning of the job cache. {doc}`shell-model` has the
interactions, the rules and the import formats.

## The scene

The **Setup** step shows the experiment as it would stand in the chamber, to scale: the beam (red) along the
axis, the target at the origin, and every detector with its rings, strips or crystals. The numbers along the beam
are millimetres from the target. A detector's circuit board is green and a γ-ray detector's housing is the
transparent box around its crystals.

**Look around.** Drag the background to turn the view and scroll to zoom. **3D**, **Side**, **Top** and **Along
beam** set the camera.

**Full screen.** The ⛶ button at the right of the scene's toolbar gives the scene the whole window: the setup
panel, the rail with the steps and the numbers beside and below the scene are hidden, and the scene fills the
height. Selecting, dragging, the camera buttons and the tracks work as before. Esc, or the same button, brings the
page back.

**Select.** Click a detector: it turns orange, its distance, angle and size appear beside it, and the panel on the
right shows its numbers.
- A particle detector shows the angles it covers, its solid angle, how much of it other detectors hide, its rate
  and its counts in the run. **Simulate its spectrum** runs the Monte Carlo for it.
- Clicking a ring, a strip or a crystal also shows that one piece: its angles, solid angle and rate.
- A γ-ray detector shows its half-angle, its coverage and efficiency, and its particle–γ coincidence rate.
  - Below them is its efficiency against energy, with what the γ rays pass through on the way.
  - **Run the source** puts a calibration source (¹⁵²Eu, ⁶⁰Co and others) at the target in place of the beam. The
    spectrum appears, and the efficiency from each strong peak is drawn on the curve.
- Click the background for the whole experiment: the totals, and advice on where the kinematics send the particles.

**Move.** Drag a selected detector.
- With **Drag changes the angle**, it moves over a sphere around the target, keeps its distance and keeps facing
  the target. With **the distance**, it moves along its own direction.
- A detector around the beam (an S3, a CD) always slides along the beam.
- The panel follows the move: angles, solid angle, shadowing and a quick estimate of the rate. When you let go,
  the detector stays where you put it and the move joins the changes that wait for **Apply** (see
  [Change, then Apply](#change-then-apply)); the panel keeps the quick numbers until then.
- A drag cannot put a detector in the beam, inside another detector or the target, or through the chamber wall. It
  stops at the limit, or turns red and returns to the last allowed place when you let go, and a message says why.
- For an exact value, type it under **Exact values** in the panel, or in the setup window (**Edit the setup**).
  Typed values are not held back; once applied, the warnings report a detector in the beam and detectors that
  overlap.

**Tracks.** **Show tracks** draws a sample of simulated events and animates them: the beam particle to the
target, then the scattered beam, the recoil and the γ ray to what they hit, which lights up. Choose how many
tracks (up to about 60 stay smooth on a laptop), which events (all, particle–γ coincidences, or one channel) and
how the sample is picked: *as in a run* picks events by the rate each stands for, *as generated* shows the rare
large-angle scatterings the generator over-represents. The scene states how the sample was chosen. Pause, the
speed slider and Clear control the animation; clicking a track shows that event's numbers in the side panel.

**A moved target.** **Position along the beam** in the Target section of the setup panel moves the target; every
angle, distance and correction follows it, and the scene draws it where it is. A target ladder goes in the
setup file (`ladder` and `selected`).

**Rings that are not safe.** For Coulomb excitation, a ring or strip that can see collisions closer than Cline's
safe distance is drawn in violet, and the panel says how many there are.

**What is drawn with typical sizes.** A setup file does not give everything a drawing needs. The target foil is
drawn 10 mm across, a circuit board 1.6 mm thick, and a γ-ray detector without crystal dimensions as long as it is
wide. These affect only the picture and the overlap check, not the rates.

## The setup as a list of decisions

**Edit the setup**, at the top right of the Setup step, opens the setup in a window of its own at the right of
the page; the scene stays in view beside it and follows what you change. **Done editing** applies what you
changed and closes the window; **Cancel** closes it without the changes. The window is the setup, as the
decisions you make in order: **Beam**, **Target and reaction**, **Particle detectors**, **γ-ray detectors**
(and, where they apply, **Ancillary detectors** and **Lifetime measurement**), **Run conditions**. Each is a card
with two lines:

- what you chose, for example *CD: annular at 30 mm, behind the target (θ 180 deg)*;
- what it leads to, for example *126°–163°, 2.31 sr, 662/s*. This line is worked out when you press **Apply**.

Open a card to change its fields. The rarely touched ones (dead layers, crystal pitch, housing, absorbers, the
backing, the coincidence window and dead time) are behind **More**. The **?** in each field explains what the
value is, a typical value, and what raising it does. **Add** above the detectors offers ready-made ones.

### Change, then Apply

Change as many fields as you like, then press **Apply**: the numbers are worked out once, for all the changes
together.

- **A change waits.** When you leave a field (Tab, Enter or a click elsewhere), its value goes into the setup you
  are editing and nothing is computed. The field gets an outline, its card a *not applied* mark and the line
  *Not applied yet*, and a bar at the top of the panel lists the changes: *Beam: energy 30 MeV → 31 MeV*. Adding,
  copying, removing or dragging a detector waits in the same way.
- **What follows at once**, because it needs no rates: the first line of each card, a problem with a value (in
  red at the top; Apply waits until it is fixed), the scene, the tables of angles beside it, and for a selected
  detector its place, what it covers and a quick estimate of its rate.
- **Everything else shows the setup as it is applied**: the steps in the rail, every step's content, the
  warnings. A note at the top of the step says so while changes wait.
- **Apply** (in the bar, in the note, or Ctrl+Enter from anywhere, a field you are typing in included) works out
  the numbers and brings every step up to date. If you pause for a second or two first, they are worked out in the
  background meanwhile and Apply shows them at once. In an experiment, Apply also saves the setup; changes that
  are not applied are not saved, and are dropped if you open another experiment or close the page.
- **Discard** puts the fields, the cards and the scene back to the applied setup.
- **Run** takes the applied setup: it asks you to apply or discard first.
- Some changes have a button of their own beside a result and apply at once, leaving waiting changes as they are:
  **Look up in ENSDF** and **Remove** for a level scheme, a matrix element typed into the level table, **Plan**
  beside a state, and the *Correlated / Isotropic* choice.
- **Apply each change as I make it**, the switch at the top of the panel, brings the old behaviour back: every
  change is applied when you leave its field. The planner remembers the switch in this browser.

**The checks are always at the top of the panel:** problems with the setup in red, then warnings (detectors
counting too fast, too few counts, angles where Rutherford's formula fails, detectors in the beam or behind
another), and *The setup changed since run 3* once you change the physics after a run.

The page itself shows the rail and the step only. The experiment's name, beam and target stay at the top of the
rail, and the line under **1 Setup** there says how many detectors there are and whether the setup has
warnings. Above the scene, a strip gives the setup in a few numbers: the particles counted per second, how many
of them come with a γ ray, the solid angle of the particle detectors and the γ-ray efficiency. Under it,
errors and warnings of the setup are written out; notes are a count that opens the window, where every check
is listed.

## Figures for a paper

Every figure has a **Paper figure** button. It opens the figure as it would be printed and lets you choose:

- **Journal style**, which sets the column widths, the typeface and the label size:

  | Style | Widths | Lettering |
  |---|---|---|
  | Physical Review (APS) | 86 mm, 178 mm | 8 pt serif |
  | Nature | 89 mm, 183 mm | 7 pt sans-serif |
  | Science | 57 mm, 121 mm, 184 mm | 7 pt sans-serif |
  | Elsevier (Nucl. Phys. A, Phys. Lett. B, NIM) | 90 mm, 140 mm, 190 mm | 8 pt sans-serif |
  | Springer (Eur. Phys. J. A) | 84 mm, 174 mm | 8 pt sans-serif |

- **Width:** one column, two columns, or the journal's middle width.
- **Format:** PDF or SVG (vector; the text stays text, so you can still edit it), or PNG at 600 dpi.

**Download figure** saves the file at exactly that width. **Download the plotted data (CSV)** saves the numbers
behind the curves, for replotting in your own tool.

The paper figures are white with black, boxed axes and inward ticks, whichever theme the page uses. Every curve has
its own line style as well as its own colour, so the figure still reads in greyscale and for colour-blind readers.
Journals change their guidelines: check the current instructions for authors before you submit.

The Report step has the same choice under **The style of the plan's figures** (**Journal style** and **Figure
width**): the report's figures (in `report.html` and in the `figures` folder of the zip) are drawn in that style.

From Python, the same figures come from `physim.nuclear.paper`:

```python
from physim.nuclear import paper
from physim.nuclear.planner import Planner

p = Planner.example("alpha_on_gold")
paper.figure(p, "kinematics", journal="nature", width="double").savefig("kinematics.pdf")
```

## Light and dark

**Light** and **Dark** at the foot of the rail on the left switch the whole page, figures included. The first
time, the planner follows your system's setting; after that it remembers your choice in this browser.
`?theme=dark` (or `light`) in the address sets it for one visit.

## While the planner computes

The page keeps answering while numbers are worked out; nothing you click has to wait for them.

- **Opening the page:** it is there at once, with a spinner in the setup panel and in the step for the second or so
  the count rates take. The nuclei and a nucleus's page are read in the background too, with a spinner meanwhile;
  the first look at a nucleus also takes ENSDF's adopted values in.
- **After Apply:** the dots beside *physim* at the top of the rail turn while the rates are worked out (about a
  second for the examples; nothing if they were worked out while you paused; a fraction of a second if you changed
  only the beam current or the run conditions, which reuse the rates). Then the checks, the cards, the rail and the
  step on screen follow, and after them the other steps, one by one, so the step you choose next is already built.
  If you choose one sooner, it shows a spinner and then its content.
- **The checks on the Plan step and the cards of the setup panel** are filled when you open them. A check you left
  open stays open after an edit.
- **A run's events** are loaded after the steps are built: seven to nine seconds, and about 1.5 GB of memory, for
  an hour of beam simulated event by event. Open the Data step before that and it shows a spinner until they are
  there. A choice on the Data step recomputes only the views it changes.

## A Coulomb-excitation plan

The example *coulex_ni58* is a complete one. To build your own (in the setup window, **Edit the setup**):

1. In the card **Target and reaction**, set *What happens in the target* to **Coulomb excitation**.
2. Choose the **Excited nucleus**: the target, or the beam.
3. Set the **Multipolarity**, usually E2.
4. Enter the **State energy** (`1.454 MeV`) and **B(Eλ↑)** (`0.0695 e2b2` or `695 e2fm4`). Take both from ENSDF
   for your nucleus. Until both are filled in, a red line names what is missing.
5. **Particle detectors** are silicon detectors. They measure the scattered beam particles and recoils, so you know
   each excited nucleus's direction and speed. Close collisions excite the state, so put detectors where those
   particles go: **Add → Backward ring around the beam (CD, 180°)** for the backscattered beam, and forward strip
   detectors for the target recoils.
6. Under **γ-ray detectors** (germanium), **Add** a detector for each crystal, with its angle, distance, radius and
   resolution.
   - Left to itself, the planner uses the typical response of such a crystal at that distance: about 0.1% for
     a clover at 25 cm and 1.3 MeV.
   - If you have measured your detector, set its **Efficiency** for your γ ray, or give an `efficiency_curve` in
     the setup file.
   - **Absorbers** takes sheets between the target and the detector, such as `Pb 1 mm, Cu 0.5 mm`.
   - A clover has two switches, **Add-back** and **Compton suppression (BGO shield)**, with the add-back and
     suppression factors beside them; the scene draws the shield around the housing.
7. On the **Plan** step, the beam time counts what the measurement uses: excitation events seen in a particle
   detector together with their γ ray. *Counts wanted* is the number of these particle–γ coincidences, so the beam
   time is usually hours. The check *The excitation, the particle–γ correlation and the Doppler shifts* holds:
   - the level schemes of the target and the beam, with **Look up in ENSDF**, their matrix elements and the Shell
     model panel;
   - the excitation probability against angle, with the excitation events per particle detector;
   - **Particle energies: elastic and after exciting the state**, at each particle detector's edges and centre,
     with the speed β of the excited nucleus that the Doppler correction needs;
   - **γ rays in coincidence: the angular correlation**, for each particle detector and crystal, relative to γ rays
     sent evenly in all directions, with the choice *Correlated with the particle* or *Isotropic*;
   - **γ rays: Doppler-shifted energy and width** for every particle-detector and γ-detector pair.
8. *Predicted spectra* shows what each detector will record. Random coincidences and the room background come from
   the `[run]` section of the setup file (`coincidence_window`, `room_background`, `extra_lines`, `dead_time`).
9. Take a run, look at its spectra on the **Data** step, and **Analyse** it on the **Analysis** step: particle gate,
   Doppler correction, peak fit, yield, normalisation (to the elastic particles, or to a known transition),
   uncertainties with a budget of systematics, the shape readings, and the beam time for the precision you want.
   Change a choice and analyse again: the result says what changed. **Check the alignment** shows what a target
   assumed off its true place does: the corrected peak with each geometry, the centroid against ring (flat when
   the geometry is right), and the offset the fit gives back. **All orders: multi-step excitation and
   reorientation**, with a level scheme looked up, solves the coupled equations (**Solve with all orders**): the γ
   yields of every transition next to first order, the 2⁺ state run prolate, spherical and oblate to see whether
   the planned run can tell them apart, **Fit the matrix element to this run's yields**, and **GOSIA input file for
   this setup**.
10. The **?** beside a number opens its explanation: the formula, the formula with this run's numbers, what it
    means and what it assumes. The report holds all of them in its appendix "Every number, explained", with the
    data's provenance, the method and what is left out.

## Experiments and runs (from Python)

What the app shows of an experiment and its runs is in `physim.nuclear.planner.Planner`, and can be used from a
notebook:

```python
from physim.nuclear.planner import Planner

p = Planner.example("coulex_ni58")
p.create_experiment()                          # ~/.physim/experiments/16O-on-58Ni, every edit saved there
p.start_run("source", duration="1 h", source="152Eu")
run = p.start_run("beam", duration="24 h")     # the first 10 min event by event, the rest scaled
run.describe()                                 # 'Beam run, 24h; the first 10min simulated event by event, ...'
p.gamma_spectra()                              # read from the run; nothing simulates on its own
p.set("target", "thickness", "1.0 mg/cm2")
p.run_status()["changes"]                      # ['target.thickness: 0.5 mg/cm2 → 1.0 mg/cm2']
```

Each run folder holds `setup.toml` (the setup it was taken with), `summary.json` (counters, rates, the Plan's
predictions at the time) and `events.npz` (named NumPy arrays: one row per counted particle and per γ ray). The
results of the experiment's analyses are rows of `results.jsonl` in the experiment's folder ({doc}`results`).

## What it does not do yet

- Elastic (Rutherford) scattering, Coulomb excitation, fusion-evaporation and lifetime measurements are planned.
  Transfer reactions and inelastic scattering come later.
- Runs of Coulomb excitation and lifetime measurements are analysed back. Runs of fusion-evaporation and of elastic
  scattering are not yet: their Analysis step says what that would need.
- Fusion-evaporation's own cross sections (physim's fusion and evaporation) are an estimate, checked against PACE4
  for light beams only; for heavy beams at high spin, take the channel table from PACE4 or a measurement.
- Literature values are not fetched from NNDC or XUNDL: they come in as data points, typed in or imported from a
  file.
- One person at a time, on their own computer; it is not meant to be hosted for a group.
