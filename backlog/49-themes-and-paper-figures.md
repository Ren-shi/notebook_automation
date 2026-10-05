# 49 · Light and dark themes, and paper figures in a journal's style

**Priority:** P1 · **Size:** M · **Area:** Nuclear planner

**Status: Done.**

## Why
From a look at the app (2026-10-05):
- **No light or dark choice:** the page was forced white, even in a dark-mode browser.
- **Figures not paper-ready:** the on-screen Plotly figures could only be saved as screenshots. Figures taken out
  of the planner should go into an article as they are.
- **Generic look:** the page was a stock NiceGUI layout.

A new design was drawn first as an HTML mockup and agreed, then applied here.

## Done
- **Themes:** **Light** / **Dark** in the header, and `?theme=` in the address.
  - The first visit follows the system's setting; the choice is remembered in the browser.
  - The page's colours are CSS variables, set once per theme. Figures are redrawn in the theme's colours
    (`app.themed`).
  - No fonts or scripts are fetched, so the app still works offline.
- **Look:** a quiet header, a strip of key results above the tabs (highest rate, longest beam time, closest
  approach, warnings), numbers in a monospace face, boxed axes with inward ticks, and the Okabe–Ito palette.
- **Paper figures:** `physim.nuclear.paper` draws each of the app's figures with Matplotlib in a journal's style.
  - Styles: Physical Review, Nature, Science, Elsevier, Springer (`paper.JOURNALS`), each with its column widths,
    typeface and label size.
  - Formats: PDF and SVG with editable text, PNG at 600 dpi, and the plotted data as CSV.
  - Every curve has its own line style as well as its own colour.
  - In the app, each figure's **Paper figure** button opens a preview at the chosen width and downloads the file.
- **Report figures:** the beam-time report draws its geometry, coverage, kinematics and spectra figures with
  `paper` too (`--journal`, `--width`; a journal choice in the app's **Report** tab). The default is Physical
  Review, one column.
- **Tests:** `test_nuclear_paper.py` checks every figure for every example, the width and type size of every
  journal style, the three file formats, and that curves differ without colour. `test_nuclear_app.py` checks the
  themed figures.

## Not done
- `physim.nuclear.plot` (the notebook helpers `setup_3d`, `coverage`, `spectra`, `kinematics`) keeps the house
  style; the report no longer uses it.
- The journal widths and type sizes were taken from the publishers' guidelines as remembered, not re-read for this
  change. They should be checked against the current instructions for authors.
