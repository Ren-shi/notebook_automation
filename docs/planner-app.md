# The experiment planner app

The planner app lets you describe an experiment by pointing and clicking. It then shows what you will measure:
geometry, kinematics, count rates and beam time, energy loss, simulated spectra and Coulomb orbits. It also writes
the report. You don't need to write any code.

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
`physim app --example oxygen_on_lead_array` opens an example; `--port` picks another port.

`physim app --window` opens the planner in its own window instead of a browser tab (`pip install
"physim-engine[app,window]"` adds pywebview, 5 MB; without it the browser opens and the command says so).

`physim app` listens only on this computer (127.0.0.1). `--host 0.0.0.0` serves the local network too, for
example to show the planner on a lab PC; anyone on that network can then open it.

## Light and dark

**Light** and **Dark** at the top switch the whole page, figures included. The first time, the planner follows your
system's setting; after that it remembers your choice in this browser. `?theme=dark` (or `light`) in the address
sets it for one visit.

The Plan step opens with four numbers: the highest rate, the longest beam time for the counts you asked for, the
head-on closest approach, and how many warnings the setup has.

## The scene

The **Setup** step shows the experiment as it would stand in the chamber, to scale: the beam (red) along the
axis, the target at the origin, and every detector with its rings, strips or crystals. The numbers along the beam
are millimetres from the target. A detector's circuit board is green and a γ-ray detector's housing is the
transparent box around its crystals.

**Look around.** Drag the background to turn the view and scroll to zoom. **3D**, **Side**, **Top** and **Along
beam** set the camera.

**Full screen.** The ⛶ button at the right of the scene's toolbar gives the scene the whole window: the setup
panel, the step bar and the numbers beside and below the scene are hidden, and the scene fills the
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
- For an exact value, type it under **Exact values** in the panel, or in the setup panel on the left. Typed values
  are not held back; once applied, the warnings report a detector in the beam and detectors that overlap.

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

The **Report** tab has the same journal and width choice: the report's figures (in `report.html` and in the
`figures` folder of the zip) are drawn in that style.

From Python, the same figures come from `physim.nuclear.paper`:

```python
from physim.nuclear import paper
from physim.nuclear.planner import Planner

p = Planner.example("alpha_on_gold")
paper.figure(p, "kinematics", journal="nature", width="double").savefig("kinematics.pdf")
```

## Experiments

The planner works on **experiments**. An experiment has a name, its current setup and the runs taken with it, and
lives in a folder under `~/.physim/experiments`, where the setup is saved every time you press **Apply**.

- The first time, the planner offers **New experiment**. It asks for the beam and its energy, the target, and what
  to measure: elastic scattering, or Coulomb excitation of the target or the beam. It names the experiment from
  them (for example *16O on 208Pb*) and puts one particle detector where the kinematics send the particles. For
  Coulomb excitation it takes the first 2⁺ state from your local ENSDF copy; without one it says what to enter.
- **Detectors from** starts from a template instead: *Coulomb excitation of the beam (62Ni on 194Pt) with
  silicon at backward angles*, *Coulomb excitation of a target nucleus with a CD and clovers*, *Elastic
  scattering with silicon pads*, *Elastic scattering with an array of strip detectors*, or *Blank*. A template
  brings its own beam and target. The first of them is what the planner shows before there is any experiment.
- The experiment's name at the top left is a menu:
  - **Open an experiment** lists your experiments, newest first; **Open** one to carry on where you left it.
  - **New experiment** makes another.
  - **Import a setup file** replaces the current setup with a `.toml` file; **Export the setup** downloads it
    (see {doc}`nuclear-setup`).

`physim app --example coulex_ni58` opens an example as an unsaved setup instead.

## The setup as a list of decisions

On the **Setup** step, the panel on the left is the setup, as the decisions you make in order: **Beam**, **Target and reaction**, **Particle
detectors**, **γ-ray detectors**, **Run conditions**. Each is a card with two lines:

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
  red at the top; Apply waits until it is fixed), the scene, the tables of angles under it, and for a selected
  detector its place, what it covers and a quick estimate of its rate.
- **Everything else shows the setup as it is applied**: the step bar, every step, the warnings. A note above
  the step bar says so while changes wait.
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

The setup panel stands beside the Setup step only. The other steps have the whole width of the window; the
experiment's name, beam and target stay in the header, and the line under **1 Setup** in the step bar says how
many detectors there are and whether the setup has warnings.

## The steps

The step bar across the top is the experiment in the order it is done: **1 Setup → 2 Plan → 3 Run → 4 Data →
5 Analysis → 6 Report**.

- **Each step has a line of status under its name:** the detectors, the beam time the counts need, the runs
  taken, the current run and its gates, the B(E2) of the analysed run, what the report holds.
- **A tick marks a step that has its result:** a setup without errors, a beam time that is enough, a run, a run's
  data, an analysis of the current run. A step that needs a run first says *needs a run* and is greyed.
- **A status in the warning colour needs attention:** warnings in the setup, a beam time too short, a run taken
  with a setup that has changed since.
- **Each step ends with one button that leads on:** *Next: check the plan*, *Next: take a run*, *Next: look at
  run 2's data*, *Next: analyse run 2*, *Next: the report*. Where there is nothing to go on to yet, the line
  says what is missing.
- The bar stays at the top of the window as you scroll, and you can choose any step at any time.

**Physics**, at the top right, is the reference behind the numbers; it is not a step.

The steps:

- **Setup:** the experiment to scale; click a detector for its numbers and drag it to move it (see
  [The scene](#the-scene)).
- **Plan:** what you check before asking for beam, answers first: the beam time your counts need against the
  beam time planned; each particle detector's rate, busy fraction, share of the excitations and safe rings; each
  γ-ray detector's efficiency and coincidence rate; the particle × γ coincidence rates. The **?** beside a number
  shows how it is worked out. Below the answers, the checks open one at a time: kinematics, energy loss, the orbit
  and safe distance, the correlation and Doppler shifts, the γ efficiency curves, the rates per strip and the
  sweep, and the **predicted spectra** (what each detector will record, from a simulated sample of the setup).
  The next run keeps these numbers, so the Analysis can put them beside what the run measured.
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
  Analysis steps. **Not simulated yet** lists what the runs leave out. See also
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
- **Analysis:** a list at the top jumps to its parts (predicted against measured, B(E2) from the peak, target
  alignment, all orders). Without a run, the analysis and the alignment check work on a fresh simulated sample
  and the step says so. First, **predicted against measured**: what the Plan said before the run (kept with it) beside
  what the run measured. That covers the rates, the coincidences in the peak, the beam time for your counts, and
  the γ-ray efficiency (measured if the experiment has a source run). Each row says whether the two agree within
  the run's statistics. Then the analysis of the current run: from the γ-ray peak to B(E2), with a gate from the
  Data step if you pick one. It says the statistical uncertainty comes from the run's own events, and what the
  whole run would give. The alignment check reads an alignment run's offset. The multi-step solution can fit the
  matrix element to the run's yields, and writes the GOSIA input file.
- **Report:** **the report**, one document. The setup with its level schemes and warnings; the plan (rates, beam
  time, kinematics, expected peaks, predicted spectra); every run with its spectra (each particle detector, the
  energy against ring with your gates drawn on it, each γ-ray detector as measured and Doppler-corrected to the
  top of the level scheme); what each gate does to the γ-ray spectrum; the Doppler correction and its
  parameters; the analysis of every analysed run; every number explained; and the method. **Show the report**
  puts it on the page; **Save as PDF** and **Save as HTML** save it; **Save everything** adds the setup file, the
  plan's tables and figures and every run's data in a zip. The PDF is printed by Edge, Chrome or Chromium; if
  none is installed the app says so, and the HTML page prints to the same PDF from any browser.
**Physics** (the button at the top right) is the physics behind the numbers, as a reference, in the order of
the physics register: nuclear
data and level schemes, kinematics, energy loss, the Rutherford orbit, the detectors and the γ-ray response,
rates and events, and Coulomb excitation. Each section links to its register page and names the tests behind it.
**See the physics** in a number's explanation, or the flask on a setup card, opens the section it belongs to.

Each result opens with **How to read this**, using your numbers, and ends with an **Explain** panel giving the
formula, the assumptions and where they stop being valid.

## While the planner computes

The page keeps answering while numbers are worked out; nothing you click has to wait for them.

- **Opening the page:** it is there at once, with a spinner in the setup panel and in the step for the second or so
  the count rates take.
- **After Apply:** the dots in the header turn while the rates are worked out (about a second for the examples;
  nothing if they were worked out while you paused; a fraction of a second if you changed only the beam current
  or the run conditions, which reuse the rates). Then the checks, the cards, the step bar and the step on
  screen follow, and after them the other steps, one by one, so the step you choose next is already built. If you
  choose one sooner, it shows a spinner and then its content.
- **The checks on the Plan step and the cards of the setup panel** are filled when you open them. A check you left
  open stays open after an edit.
- **A run's events** are loaded after the steps are built: seven to nine seconds, and about 1.5 GB of memory, for
  an hour of beam simulated event by event. Open the Data step before that and it shows a spinner until they are there. A
  choice on the Data step recomputes only the views it changes.

## A Coulomb-excitation plan

The example *coulex_ni58* is a complete one. To build your own:

1. Under **Reaction**, set *What happens in the target* to **Coulomb excitation**.
2. Choose the **Excited nucleus**: the target, or the beam.
3. Set the **Multipolarity**, usually E2.
4. Enter the **State energy** (`1.454 MeV`) and **B(Eλ↑)** (`0.0695 e2b2` or `695 e2fm4`). Take both from ENSDF
   for your nucleus. Until both are filled in, a red line names what is missing.
5. **Particle detectors** are silicon detectors. They measure the scattered beam particles and recoils, so you know
   each excited nucleus's direction and speed.
   - Close collisions excite the state, so put detectors where those particles go: **Add → Backward ring around the
     beam** for the backscattered beam, and forward strip detectors for the target recoils.
   - *Geometry* names the detector with the largest share of excitation events.
6. Under **γ-ray detectors** (germanium), **Add** a detector for each crystal, with its angle, distance, radius and
   resolution. They appear in *Geometry* as bronze crystals.
   - Left to itself, the planner uses the typical response of such a crystal at that distance: about 0.1% for
     a clover at 25 cm and 1.3 MeV.
   - If you have measured your detector, set its **Efficiency** for your γ ray, or give an `efficiency_curve` in
     the setup file.
   - **Absorbers** takes sheets between the target and the detector, such as `Pb 1 mm, Cu 0.5 mm`.
7. **Rates and beam time** counts what the measurement uses: excitation events seen in a particle detector
   together with their γ ray. *Counts wanted* is the number of these particle–γ coincidences, so the beam time is
   usually hours. The table shows, per detector:
   - all particles;
   - the excitation events;
   - the excitation events with a γ ray.
8. Open **Excitation and γ rays**:
   - **Excitation probability:** against angle, with the excitation events per particle detector.
   - **Particle energies:** the energy of each particle at each particle detector's edges and centre, elastic and
     after exciting the state. It also gives the speed β of the excited nucleus, which is what you need to correct
     the γ-ray energies for the Doppler shift.
   - **All levels excited from the ground state:** with a level scheme looked up, the first-order cross section
     of every level, directly and with feeding from above, and of every γ ray.
   - **Solve with all orders:** with a level scheme looked up, the coupled equations give the γ yields of
     every transition with multi-step excitation and reorientation, next to first order; the 2⁺ state is run
     prolate, spherical and oblate to see whether the planned run can tell them apart; and a GOSIA input file
     for the setup can be saved.
   - **The angular correlation:** for each particle detector and crystal, the γ rays seen in coincidence relative
     to γ rays sent evenly in all directions. A switch turns the correlation off (isotropic emission).
   - **Doppler table:** the shifted γ-ray energy and peak width for every particle-detector and γ-detector pair.
9. On **Spectra**, **Simulate the γ rays** follows the γ ray of every excited event to the crystals: the
   particle × γ matrix of true and random coincidences, and each γ-ray detector's spectrum as measured and
   Doppler-corrected from the segment and crystal that fired. Random coincidences and the room background come
   from the `[run]` section of the setup file (`coincidence_window`, `room_background`, `extra_lines`,
   `dead_time`). A clover has two switches in the setup panel, **Add-back** and **Compton suppression (BGO
   shield)**, with the add-back and suppression factors beside them; its spectrum is then shown with and
   without them, and the scene draws the shield around the housing.
10. **Analyse**, below it, extracts B(E2) from the simulated γ-ray peak as an experimentalist would: particle gate,
    Doppler correction, peak fit, yield, normalisation (to the elastic particles, or to a known transition),
    uncertainties with a budget of systematics, the shape readings, and the beam time for the precision you want.
    Change a choice and analyse again: the result says what changed.
11. **Check the alignment**, below the analysis, shows what a target assumed off its true place does: the
    corrected peak with each geometry, the diagnostic plot of centroid against ring (flat when the geometry is
    right), and the offset the fit gives back. The scene shows the assumed target as a faint outline.
12. The **?** beside a number in the scene's side panel, or in the analysis, opens its explanation: the formula,
    the formula with this run's numbers, what it means and what it assumes. The report on the Report tab holds all
    of them in its section "Every number, explained", with the data's provenance, the method and what is left
    out.

## Experiments and runs (from Python)

The app's model is being rebuilt around experiments and runs (backlog 63–70). The model is in place and can be
used from a notebook; the tabs that show it follow.

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
predictions at the time) and `events.npz` (named NumPy arrays: one row per counted particle and per γ ray).

## What it does not do yet

- Elastic (Rutherford) scattering and Coulomb excitation of one state. Transfer and fusion-evaporation come
  later.
- One person at a time, on their own computer; it is not meant to be hosted for a group.
