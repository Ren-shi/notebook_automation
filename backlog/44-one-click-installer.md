# 44 · One-click installer for the planner app

**Priority:** P2 · **Size:** M · **Area:** Nuclear planner

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
