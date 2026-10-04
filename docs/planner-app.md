# The experiment planner app

The planner app lets you describe an experiment by pointing and clicking. It then shows what you will measure:
geometry, kinematics, count rates and beam time, energy loss, simulated spectra and Coulomb orbits. It also writes
the beam-time report. You don't need to write any code.

It runs on your own computer, in your web browser. There are no accounts, and nothing is uploaded.

## Install and start it

**If you have never used Python:**

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

A one-click installer that needs no Python is planned (backlog item 44).

## A four-detector setup, step by step

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
   - **Add** makes a new one, a 5 mm disc at 45° and 100 mm, which you then edit. Make four, for example at 30°,
     60°, 90° and 150°.
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

   Each tab ends with an **Explain** panel giving the formula, the assumptions and where they stop being valid.
7. **Try a sweep.** In *Rates and beam time*, under **Sweep one parameter**, choose *beam energy*, enter
   `4 MeV, 5 MeV, 6 MeV, 7 MeV`, and press *Run sweep* to see how the rate changes.
8. **Save your setup.** **Save setup** downloads `setup.toml`, and **Load setup** opens it again later. The file is
   plain text (see {doc}`nuclear-setup`).
9. **Export the report.** In **Report**, press *Build and download the report*. You get a zip archive with
   `report.html` (open it in the browser, and print it to PDF if you need one), CSV tables, figures, the setup file
   and, if `root` was included in the install, `events.root` for ROOT.

## What it does not do yet

- Only elastic (Rutherford) scattering. Inelastic scattering and Coulomb excitation are next (backlog item 43).
- One person at a time, on their own computer; it is not meant to be hosted for a group.
