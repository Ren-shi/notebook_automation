"""Inputs for the trusted tools, so a physim result can be checked in SRIM/TRIM or LISE++ in a few clicks.

::

    from physim.nuclear import Experiment, export

    exp = Experiment.example("oxygen_on_lead_array")
    open("TRIM.IN", "w").write(export.trim_in_for(exp))      # the beam through the target and backing
    print(export.lise_settings(exp))                         # what to enter in LISE++, and physim's own numbers

``trim_in`` writes the ``TRIM.IN`` file that TRIM reads at start-up ("TRIM-96 / TRIM Calculation from a
setup file"): put it in the SRIM folder and start TRIM. Only the ion, energy, layers and the transmitted-ions file
matter for energy loss; damage energies are SRIM's usual defaults and do not affect the transmitted energies.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

from . import data
from .quantity import Quantity

#: Displacement energies (eV) SRIM proposes for a few light elements; 25 eV otherwise.
_DISPLACEMENT = {1: 10.0, 6: 28.0, 7: 28.0, 8: 28.0}


def _element_masses(material: data.Material) -> dict:
    """Mean atomic mass (u) of each element as it occurs in the material (an enriched isotope keeps its own)."""
    num: dict = {}
    den: dict = {}
    for z, a, n in material.atoms:
        m = data.nuclide((z, a)).atomic_mass_u if a is not None else data.element(z).molar_mass
        num[z] = num.get(z, 0.0) + n * m
        den[z] = den.get(z, 0.0) + n
    return {z: num[z] / den[z] for z in num}


def _layer_atoms(material: data.Material) -> dict:
    """Atom fractions by element Z (isotopes of one element merged, as TRIM treats elements)."""
    out: dict = {}
    for z, _, n in material.atoms:
        out[z] = out.get(z, 0.0) + n
    total = sum(out.values())
    return {z: n / total for z, n in out.items()}


def trim_in(ion: str, energy_mev: float, layers: Sequence, ions: int = 2000, title: Optional[str] = None,
            angle_deg: float = 0.0, seed: int = 0) -> str:
    """The text of a ``TRIM.IN`` file: ``ion`` at ``energy_mev`` into ``layers``, a list of
    (material, thickness) or (material, thickness, density) with thickness a :class:`Quantity` or text such as
    ``"10 um"`` or ``"0.5 mg/cm2"``. Calculation type "Ion distribution and quick calculation of damage", with
    transmitted ions written to ``TRANSMIT.txt``."""
    n = data.nuclide(ion)
    rows = []
    for layer in layers:
        mat_name, thickness = layer[0], layer[1]
        density = layer[2] if len(layer) > 2 else None
        mat = mat_name if isinstance(mat_name, data.Material) else data.material(mat_name, density)
        if mat.density_g_cm3 is None:
            raise ValueError(f"TRIM needs a density for {mat.name}: give one")
        q = Quantity.parse(thickness) if isinstance(thickness, str) else thickness
        mg = mat.areal_density_mg_cm2(q)
        width_ang = mg * 1e-3 / mat.density_g_cm3 * 1e8
        rows.append((mat, width_ang, _layer_atoms(mat)))
    elements: list = []
    masses: dict = {}
    for mat, _, atoms in rows:
        for z, m in _element_masses(mat).items():
            masses.setdefault(z, m)
        for z in atoms:
            if z not in elements:
                elements.append(z)
    title = title or f"{n.name} ({energy_mev * 1e3:g} keV) into " + " + ".join(m.name for m, _, _ in rows)
    title = title.replace('"', "'")[:60]
    depth = sum(w for _, w, _ in rows)
    lines = [
        "==> SRIM-2013.00 This file controls TRIM Calculations.",
        "Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number.",
        f"{n.Z:6d} {n.atomic_mass_u:9.4f} {energy_mev * 1e3:11.3f} {angle_deg:6.2f} {ions:6d} 1 {ions:8d}",
        "Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders",
        f"{1:23d} {seed:35d} {0:7d}",
        "Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), "
        "Special EXYZ.txt file",
        f"{0:27d} {0:7d} {1:11d} {0:7d} {0:15d} {0:31d}",
        "Target material : Number of Elements & Layers",
        f'"{title}" {len(elements):16d} {len(rows):15d}',
        "PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]",
        f"{5:8d} {0:25d} {depth:14.0f}",
        "Target Elements:    Z   Mass(amu)",
    ]
    for i, z in enumerate(elements, start=1):
        lines.append(f"Atom {i} = {data.element(z).symbol} = {z:10d} {masses[z]:8.3f}")
    head = "".join(f"{data.element(z).symbol}({z})".rjust(9) for z in elements)
    lines.append(f"Layer   Layer Name /               Width Density {head}")
    lines.append("Numb.   Description                (Ang) (g/cm3) " + "   Stoich" * len(elements))
    for k, (mat, width, atoms) in enumerate(rows, start=1):
        name = f'"{mat.name[:20]}"'
        stoich = "".join(f"{atoms.get(z, 0.0):9.6f}" for z in elements)
        lines.append(f"{k:2d}      {name:<24s} {width:11.1f} {mat.density_g_cm3:7.4f}{stoich}")
    lines += [
        "0  Target layer phases (0=Solid, 1=Gas)",
        " ".join("0" for _ in rows),
        "Target Compound Corrections (Bragg)",
        " ".join("1" for _ in rows),
        "Individual target atom displacement energies (eV)",
        "".join(f"{_DISPLACEMENT.get(z, 25.0):8.0f}" for z in elements),
        "Individual target atom lattice binding energies (eV)",
        "".join(f"{3.0:8.0f}" for _ in elements),
        "Individual target atom surface binding energies (eV)",
        "".join(f"{2.0:8.2f}" for _ in elements),
        "Stopping Power Version (1=2011, 0=2011)",
        " 0",
    ]
    return "\r\n".join(lines) + "\r\n"


def _stack_layers(experiment) -> list:
    t = experiment.target
    layers = [(t.material_data(), t.thickness)]
    if t.backing is not None:
        layers.append((t.backing.material_data(), t.backing.thickness))
    return layers


def trim_in_for(experiment, ions: int = 2000) -> str:
    """``TRIM.IN`` for an experiment's beam going through its target (and backing), at the target tilt."""
    from .rates import tilt_deg

    beam = experiment.beam
    ion = data.nuclide((beam.Z, beam.A)).name
    return trim_in(ion, beam.energy_mev, _stack_layers(experiment), ions=ions, angle_deg=tilt_deg(experiment),
                   title=f"{experiment.title}"[:60])


def lise_settings(experiment) -> str:
    """What to enter in LISE++ to reproduce physim's numbers for an experiment, with physim's values to compare."""
    from .kinematics import TwoBody
    from .rates import Rates, beam_ion, stack, tilt_deg

    beam = experiment.beam
    ion = beam_ion(experiment)
    layers = stack(experiment)
    r = Rates(experiment)
    lines = [
        f"LISE++ settings for: {experiment.title}",
        "",
        f"Projectile: {ion}, {beam.energy_mev / beam.A:.4f} MeV/u ({beam.energy_mev:g} MeV), charge state "
        f"{beam.charge}.",
    ]
    for lay in layers:
        lines.append(f"{lay.name.capitalize()}: {lay.material.name}, {lay.thickness:.4g} mg/cm2"
                     + (f", tilted {tilt_deg(experiment):g} deg" if lay.name == "target" and tilt_deg(experiment)
                        else "") + ".")
    lines += [
        "Kinematics calculator: reaction 'Scattering'; set the target very thin (or the reaction point to the target",
        "entrance) to compare at the full beam energy, as below.",
        "",
        "Detector  lab angle (deg)  physim E_ejectile (MeV)  physim dsigma/dOmega lab (mb/sr)  physim counts/s",
    ]
    for g in r.array:
        theta = g.mean_theta()
        for ch in r.channels:
            if ch.layer != 0:
                continue
            tb = TwoBody(ion, ch.nuclide.name, beam.energy_mev)
            first, _ = tb.at_lab(theta)
            if not math.isfinite(first.energy):
                continue
            from .rutherford import Rutherford

            sigma = Rutherford(ion, ch.nuclide.name, beam.energy_mev).cross_section_lab(theta)[0]
            lines.append(f"{g.name:9s} {theta:15.3f} {first.energy:24.4f} {sigma:33.5g} "
                         f"{r.by_channel(g.name).get((ch.label, 'ejectile'), 0.0):16.4g}   ({ch.nuclide.name})")
    return "\n".join(lines) + "\n"


__all__ = ["lise_settings", "trim_in", "trim_in_for"]
