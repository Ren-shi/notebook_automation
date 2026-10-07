# 61 · The planner in its own desktop window

**Priority:** P2 · **Size:** S · **Area:** Experiment workbench

**Status: Done, except the test of the installed shortcut on a clean Windows machine (the user, or someone
other than the author).**

> **Done** (`--window`, `--browser` and `wants_window` in `app.py`, the `window` extra in `pyproject.toml`, the
> installer bundle with `[app,root,window]`, the shortcut's comment and the guide's install steps).
> - **Decision: NiceGUI's native mode (pywebview on Edge WebView2),** revisited as the design note asked:
>   - **Download size:** pywebview 6.2.1 with pythonnet, clr-loader, bottle and proxy-tools adds 4.6 MB
>     installed (item 44 had guessed 20–30 MB). WebView2 itself comes with Windows 10 and 11.
>   - **Antivirus:** nothing new is executed but the bundled Python; pywebview is pure Python plus pythonnet's
>     signed .NET loader.
>   - **The scene:** draws in the window as in a browser (the same Edge engine): checked with a screenshot of
>     the window on the Geometry tab, WebGL scene included.
> - **Behaviour:** `--desktop` (the shortcut) and `--window` open the window; `--browser` or `--no-browser`
>   keep the browser (the CI test and the idle-exit test use `--desktop --no-browser`, unchanged). Without
>   pywebview the browser opens and the command says so. Closing the window stops the planner.
> - **Results:** the window opens here with the full page (1400 × 900) and no address bar; the server behind it
>   answers `/physim-planner` as before; the scene renders inside it.
> - **Not done here, or to confirm:**
>   - The installed shortcut on a clean Windows machine without Python: the installer workflow builds the
>     bundle with the `window` extra, but the window itself is not opened in CI (no display); the user's
>     test stands, as in item 44.
>   - The scene's frame rate in the window was not measured, only that it draws.
>   - macOS and Linux: pywebview uses the system's web view there too, untested here.

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
