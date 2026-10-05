# 51 · Detector catalogue: S3, clover and LaBr₃ as solids

**Priority:** P1 · **Size:** M · **Area:** Experiment workbench

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
