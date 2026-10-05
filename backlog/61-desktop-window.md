# 61 · The planner in its own desktop window

**Priority:** P2 · **Size:** S · **Area:** Experiment workbench

## Why
The planner is an application that calculates on the user's machine, but it opens in a browser tab, which makes it
look like a website. It should open as its own window.

## Scope
- **A desktop window** with no browser tab or address bar, running the same app. The calculations already run
  locally in Python and Rust; only the window changes.
- **The browser remains available** (`physim app`), for a server with more hardware that is reached from another
  machine.
- **Installer:** the shortcut of item 44 opens the window.
- **Left out:** rewriting the interface in a native toolkit.

## Design notes
- Item 44 compared NiceGUI's native-window mode (pywebview on Edge WebView2) and did not choose it for the
  installer. Revisit that comparison: download size, antivirus false positives, and whether the scene of item 52
  draws at full speed inside the window.

## Depends on
44, 52.

## Done when
- The installed shortcut opens the planner in its own window on a clean Windows machine.
- The scene of item 52 performs in the window as it does in a browser.
- The added download size is recorded here.
