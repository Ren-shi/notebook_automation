# 44 · One-click installer for the planner app

**Priority:** P2 · **Size:** M · **Area:** Nuclear planner

**Status: in progress.** The Windows installer is built and test-installed in CI. What remains:
- a test on a clean machine by someone other than the author;
- a first tagged release, so the installer appears on the releases page.

**Size (first CI run, 2026-10-04):** the installer is 57 MB to download and 289 MB once installed.

## Why
Many students who would use the planner have never installed Python. Today they follow the steps in
`docs/planner-app.md`: install Python, then run `pip install` in a terminal. A download that starts the app with a
double-click removes that step.

## Scope
- Windows first, then macOS: a single download that installs everything needed and adds a "physim planner" shortcut,
  which starts `physim app` and opens the browser.
- **Options to compare:**
  - a bundled Python ("python-build-standalone" plus the physim wheel);
  - PyInstaller or Briefcase;
  - NiceGUI's native-window mode (`ui.run(native=True)`, via pywebview) for a desktop window instead of a browser
    tab.

  Compare on download size, antivirus false positives, code-signing needs and update path, and record the decision
  here.
- Built in CI for each release, and attached to the GitHub release.
- Left out: Linux packages (users there can use pip), automatic updates.

## Depends on
41.

## Done when
- On a clean Windows machine without Python, downloading and double-clicking starts the planner. Tested by someone
  other than the author.
- The installer is built automatically for each release, and its size is recorded here.

## Decision (2026-10-04)

**Windows: the official embeddable Python plus wheels, packed with Inno Setup.**

| | Embeddable Python + Inno Setup (chosen) | PyInstaller / `nicegui-pack` | Briefcase | Native window (pywebview) |
|---|---|---|---|---|
| What runs | python.org's `pythonw.exe`, signed by the PSF | PyInstaller's own unsigned bootloader | its own stub app | adds pywebview and pythonnet, on Edge WebView2 |
| Antivirus false positives | rare: no unsigned executable of ours apart from the setup program | common, a known PyInstaller problem, worst in one-file mode | moderate | as the packaging underneath |
| Size | about the installed packages (NumPy, Matplotlib, Plotly, NiceGUI, uproot) | about the same | about the same | +20–30 MB |
| Code signing | only the setup program would need it; unsigned it shows SmartScreen's "More info → Run anyway" | the bootloader and the setup | the MSI | — |
| Update path | run a newer installer over the old one; it replaces the bundled packages | the same | MSI upgrade | — |
| Build | a 150-line script and an `.iss` file; reproducible in CI | a spec file per hidden import (NiceGUI, uvicorn) | a Toga-oriented project layout | — |

- The planner stays a **browser tab**. A native window would add a dependency and a WebView2 requirement for little
  gain. `ui.run(native=True)` can be revisited if students ask for a separate window.
- The shortcut runs `pythonw.exe -m physim app --desktop`. Desktop mode:
  - logs to `%LOCALAPPDATA%\physim\planner.log`, as there is no console;
  - opens the running planner if one is already up;
  - takes a free port if 8080 is busy;
  - stops 60 s after the last browser tab closes, so no hidden process lingers.
- The app now listens on 127.0.0.1 only, so Windows Firewall does not ask on first start.
- Per-user install (no administrator rights) to `%LOCALAPPDATA%\Programs\physim planner`.
- Python 3.12.10, the last 3.12 release with Windows binaries; `build.py` pins it, and the workflow's Python must
  match (a test checks this).

**Deferred by the user (2026-10-04): code signing and macOS.**
- Unsigned, the installer shows SmartScreen's "Windows protected your PC". The install guide tells students to
  click **More info → Run anyway**.
- On macOS, `pip install` works meanwhile. A later bundle would use python-build-standalone in a `.app` inside a
  `.dmg`. Without an Apple Developer ID ($99 per year) and notarisation, Gatekeeper would block it until the user
  right-clicks → Open.

## Progress (2026-10-04)

- **Bundle:** `installer/windows/build.py` downloads the embeddable Python and enables `site-packages` in its
  `._pth`. It then pip-installs the physim wheel with `[app,root]` for win_amd64, prunes test folders, precompiles,
  draws the icon (a Coulomb orbit) and checks that the bundled Python imports everything from inside the bundle.
- **Installer:** `installer/windows/planner.iss` does a per-user install with Start-menu and optional desktop
  shortcuts, an uninstaller, a clean replacement of the old packages on upgrade, and "start now" at the end.
- **CI:** `.github/workflows/installer.yml` runs on version tags, PRs touching the installer or the app, and manual
  runs.
  - It builds the wheel, the bundle and the installer, and installs it silently.
  - With Python removed from the PATH, it starts the planner as the shortcut does, loads every example page,
    checks that a second launch reuses the running planner, and uninstalls.
  - It writes the installer and installed sizes to the job summary, and on a tag attaches the installer to the
    GitHub release.
- **Tests:** `tests/python/test_nuclear_app.py`
  - `test_desktop_launch`: the log file, reuse of a running planner, and idle exit;
  - `test_ports`;
  - `test_installer_bundle_helpers`: the `._pth` edit, the icon, and the Python version matching the workflow.
- **First CI run (#40):** the bundled Python imported everything from inside the bundle. With Python off the PATH,
  the installed planner started and served all three examples, a second launch reused it, and the uninstaller
  removed it.
- **Checked by hand:** in a browser, an open tab keeps the planner running past `--idle-exit`; closing it stops
  the planner.
