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
3. Click **physim planner**. Your browser opens the planner after a few seconds. There is no window to keep open:
   the planner stops by itself a minute after you close its last browser tab. Clicking the shortcut again while it
   runs just opens another tab.
4. To remove it: *Settings → Apps → physim planner → Uninstall*. To update, run a newer installer over the old one.

If the browser does not open, the log is in `%LOCALAPPDATA%\physim\planner.log`.

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

`physim app` listens only on this computer (127.0.0.1). `--host 0.0.0.0` serves the local network too, for
example to show the planner on a lab PC; anyone on that network can then open it.

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
   - **Geometry:** turn the 3D view with the mouse.
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
   resolution. They appear in *Geometry* as dashed circles.
   - Set each crystal's **Efficiency**: its full-energy-peak efficiency for your γ ray, typically 0.5–3% per crystal
     at 1.3 MeV.
   - Left empty, the planner uses the crystal's geometric coverage, an upper limit that makes the beam time look
     shorter than it will be.
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
   - **Doppler table:** the shifted γ-ray energy and peak width for every particle-detector and γ-detector pair.

## What it does not do yet

- Elastic (Rutherford) scattering and Coulomb excitation of one state. Transfer and fusion-evaporation come
  later.
- One person at a time, on their own computer; it is not meant to be hosted for a group.
