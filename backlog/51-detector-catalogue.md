# 51 · Detector catalogue: S3, clover and LaBr₃ as solids

**Priority:** P1 · **Size:** M · **Area:** Experiment workbench

**Status: Done, except the values marked typical, which need the drawings of the detectors used (the user).**

> **Done** (`python/physim/nuclear/catalogue.py`, `python/physim/nuclear/data/detector_models.toml`, models,
> crystals and the chamber in `experiment.py`, dead material in `detectors.py` and `events.py`, crystals in
> `gamma.py`, drawing in `planner.py` and `app.py`, presets in `guide.py`,
> `tests/python/test_nuclear_catalogue.py`, guide section "Detector models and the chamber" in
> `docs/nuclear-setup.md`).
> - **Models:** `S3`, `clover` and `LaBr3_2x2`, each with its sources. A setup names a model and overrides what
>   differs.
> - **From the manufacturers' documents (read 2026-10-05):**
>   - S3 (Micron's product page): 24 rings, 32 sectors, active area 22 to 70 mm in diameter, chip 20 to 76 mm.
>   - Clover (Mirion's specification sheet): four crystals of 50 × 70 mm, gaps of at most 0.7 mm, 21–22%
>     relative efficiency per crystal, below 2.1 keV at 1.33 MeV.
>   - LaBr₃ (Saint-Gobain's technical note): 2.9% at 662 keV, 2.1% at 1332 keV, 1.6% at 2615 keV.
> - **Ring pitch:** recorded as 1 mm. Micron's "junction pitch" of 886 µm does not fit 24 rings in 24 mm; it is
>   taken as the strip width, and the rings are treated as 1 mm wide with no gap, because charge between strips
>   is shared by the neighbours. This reading should be confirmed with Micron's drawing.
> - **Dead material:** a model's board sits 0.01 mm behind its active face; a γ-ray detector's housing is a
>   square at its front window. Both stop particles in the analytic rates and in the Monte Carlo, with no change
>   to the Rust generator (a face whose threshold cannot be reached).
> - **Crystals:** a clover is four crystal faces; the Doppler table has a row per crystal.
> - **Chamber:** a `[chamber]` section, with detectors checked against it.
> - **Results:**
>   - Each S3 ring's solid angle matches the hand calculation to 0.1% (the whole detector to 1e-12).
>   - A pad behind the S3's board, and a pad behind a clover, count nothing in the rates and in the events.
>   - The Monte Carlo agrees with the analytic rate of an S3 within 4σ.
>   - Setups without a model give unchanged results (all earlier tests pass).
> - **Typical values, not from a document** (marked in the catalogue):
>   - S3: thickness (1000 µm), dead layer (0.5 µm), and the board's size (50 mm from the axis).
>   - Clover: the distance between crystal centres (45 mm), the housing's side (101 mm) and the window gap.
>   - LaBr₃: the housing and the single resolution value (item 53 makes it depend on energy).
> - **Not done here:** models added by the user in a file of their own (a setup can still give every dimension
>   itself); solids for drawing to scale (the scene of item 52 builds them from these numbers); the
>   chamber wall's effect on γ rays (item 53); tapered crystal fronts.

## Why
A detector today is an infinitely thin face (item 38), and a γ detector is a circle with an efficiency number. There
are no named models. The workbench needs real hardware: its dimensions, its segmentation and the volume it occupies.

## Scope
- **Catalogue:** a data file of named models. Each entry gives the dimensions, the segmentation, the material and
  the source of the numbers. Every dimension can be edited, and users can add their own models.
- **First models:**
  - **S3 silicon detector:** an annulus with 24 rings and 32 sectors, with its real inner and outer active radii,
    ring pitch and thickness, on a circuit board with the central hole.
  - **Germanium clover:** four crystals in a common housing, each with its own position, diameter and length.
  - **LaBr₃:** a single cylindrical crystal in a housing (2″ × 2″ by default).
- **Ring pitch is stored,** where today it is implied by equal division between the inner and outer radius.
- **Solids:** each model has active volumes and blocking volumes (housing, circuit board). A particle or γ ray that
  meets a blocking volume first does not reach what lies behind it.
- **γ detectors join the array,** so they take part in the shadowing check of item 48 and can be oriented.
- **Chamber:** a generic chamber with a radius, a wall thickness and material, and a beam pipe. No facility is
  built in.
- **Setup file:** a detector may name a model (`model = "S3"`) and override any dimension. The existing shapes
  (rectangle, circle, annular) keep working.
- **Left out:** curved or tapered detail of real housings; mesh import; BGO shields (item 62); the γ response
  (item 53).

## Design notes
- One description serves both the drawing (item 52) and the physics, in Python and in the Rust core.
- **First task of this item:** look up each model on its manufacturer's website (Micron Semiconductor for the S3;
  the makers of the clover and of the LaBr₃ crystal) and take the default dimensions from the data sheets. Cite
  each sheet in the catalogue.
  - Where a data sheet does not give a dimension (crystal lengths and housings are often missing), use a published
    description of such a detector, and mark the value as typical.
  - No dimensions have been looked up yet; the numbers in this file are placeholders until then.

## Depends on
38, 48.

## Done when
- The solid angle of each S3 ring matches a hand calculation from its radii and distance to 0.1%.
- A clover placed in front of another detector hides it in the rates and in the Monte Carlo alike.
- Each model's dimensions are checked against its data sheet, and the source is recorded.
- Setup files written before this item give unchanged results.
