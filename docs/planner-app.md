# The experiment planner app

The planner app lets you describe an experiment by pointing and clicking. It then shows what you will measure:
geometry, kinematics, count rates and beam time, energy loss, simulated spectra and Coulomb orbits. It also writes
the beam-time report. You don't need to write any code.

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
`physim app --example oxygen_on_lead_array` starts from another example; `--port` picks another port.

`physim app --window` opens the planner in its own window instead of a browser tab (`pip install
"physim-engine[app,window]"` adds pywebview, 5 MB; without it the browser opens and the command says so).

`physim app` listens only on this computer (127.0.0.1). `--host 0.0.0.0` serves the local network too, for
example to show the planner on a lab PC; anyone on that network can then open it.

## Light and dark

**Light** and **Dark** at the top switch the whole page, figures included. The first time, the planner follows your
system's setting; after that it remembers your choice in this browser. `?theme=dark` (or `light`) in the address
sets it for one visit.

The strip above the results always shows four numbers: the highest rate, the longest beam time for the counts you
asked for, the head-on closest approach, and how many warnings the setup has.

## The scene

The **Geometry** tab shows the experiment as it would stand in the chamber, to scale: the beam (red) along the
axis, the target at the origin, and every detector with its rings, strips or crystals. The numbers along the beam
are millimetres from the target. A detector's circuit board is green and a γ-ray detector's housing is the
transparent box around its crystals.

**Look around.** Drag the background to turn the view and scroll to zoom. **3D**, **Side**, **Top** and **Along
beam** set the camera.

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
  the scene and the panel settle first and the other tabs follow, the Monte Carlo included.
- A drag cannot put a detector in the beam, inside another detector or the target, or through the chamber wall. It
  stops at the limit, or turns red and returns to the last allowed place when you let go, and a message says why.
- For an exact value, type it under **Exact values** in the panel, or in the setup panel on the left. Typed values
  are not held back; the warnings report a detector in the beam and detectors that overlap.

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

## The guided mode

The planner opens in **Guided** mode: the left panel is a list of steps, in the order you would plan an experiment.

1. **What do you want to measure?** Rutherford (elastic) scattering, or Coulomb excitation. **Start here** loads a
   complete example of that kind, so every result is filled in from the start. You then change it one step at a
   time.
2. **Beam**, with the kinematics beside it.
3. **Target**, with the energy loss.
4. **Particle detectors**, with the geometry and kinematics.
5. **γ-ray detectors** (Coulomb excitation only), with the excitation and Doppler results.
6. **Rates and beam time**.
7. **Spectra**.
8. **Report**.

In every step:
- **The ? in each field** explains what the value is, a typical value, and what raising it does to the results.
  Hover over it, or tap it on a touch screen.
- **How to read this** opens every result. It says what to look for, using your setup's numbers: which detector
  counts fastest, whether two peaks separate, how much beam time you need.
- **Problems show in red inside the step.** **Next** stays greyed out until they are fixed. You can go back to any
  step you have already reached by clicking its title.

**Expert view** (top right) shows every input on the left and every result as tabs on the right, as in the
walkthrough below. The browser remembers which mode you used last.

## A four-detector setup in the expert view

1. **Start from an example.** At the top right, pick *alpha_on_gold* in **Start from example**. The setup appears on
   the left and the results on the right.
2. **Set the beam.** Under *Beam*, type an **Energy** such as `6 MeV`, or `1.5 MeV/u`, and press Enter. Every
   value needs its unit; the field shows an example when it is empty. The results update at once.
3. **Set the target.** Under *Target*, type the **Material** (`Au`, `208Pb`, `CD2`, `Mylar`) and **Thickness**
   (`0.5 mg/cm2` or `250 nm`). For a carbon backing, fill in *Backing material* and *Backing thickness*.
4. **Arrange the detectors.**
   - Open a detector in the *Detectors* list to change its angle (θ), distance, size, thickness, resolution or
     threshold.
   - **Duplicate** copies a detector, and **Remove** deletes it.
   - **Add** offers a ready-made detector for each region, which you then edit:
     - a forward strip detector (45°);
     - a side pad (90°);
     - a backward pad (150°);
     - a backward ring around the beam (a CD at 180°, covering about 125–165°).

     θ below 90° is forward and above 90° is backward. Backward detectors count slowly but see the closest
     collisions. Make four, for example at 30°, 60°, 90° and 150°.
5. **Read the warnings.** The yellow box above the results lists everything to check: detectors counting too fast,
   too few counts in the planned beam time, angles where Rutherford's formula fails, detectors in the beam.
   - A red line means the setup has a problem, and the message names the field.
   - While a red line is showing, the results stay on the last valid setup.
6. **Look at each tab.**
   - **Geometry:** the experiment to scale; click a detector for its numbers and drag it to move it (see
     [The scene](#the-scene)).
   - **Kinematics:** energy against angle, with each detector's range shaded.
   - **Rates and beam time:** the rate per detector and per strip, the counts in the run, and the beam time for the
     counts you want.
   - **Energy loss:** what the target, backing and dead layers take from the beam.
   - **Spectra:** press *Simulate* for more events or another seed.
   - **Trajectories:** Coulomb orbits.
   - **Excitation and γ rays:** for a Coulomb-excitation setup (see the example *coulex_ni58*), the excitation
     probability against angle, excitation events per detector, and the Doppler-shifted γ-ray energy and width
     for every pair of particle and γ detector.

   Each tab opens with **How to read this** and ends with an **Explain** panel giving the formula, the assumptions
   and where they stop being valid.
7. **Try a sweep.** In *Rates and beam time*, under **Sweep one parameter**, choose *beam energy*, enter
   `4 MeV, 5 MeV, 6 MeV, 7 MeV`, and press *Run sweep* to see how the rate changes.
8. **Save your setup.** **Save setup** downloads `setup.toml`, and **Load setup** opens it again later. The file is
   plain text (see {doc}`nuclear-setup`).
9. **Export the report.** In **Report**, press *Build and download the report*. You get a zip archive with
   `report.html` (open it in the browser, and print it to PDF if you need one), CSV tables, figures, the setup file
   and, if `root` was included in the install, `events.root` for ROOT.

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
    the formula with this run's numbers, what it means and what it assumes. **Show the run record** on the
    Report tab puts all of them on one page with the setup, the data's provenance, the method and what is left
    out; the report's zip holds the same page as `record.html`, to print to PDF.

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
