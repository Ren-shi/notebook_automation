"""The reproduced fusion-evaporation experiment (backlog 75, D10): Usman et al., NIM B 469, 42 (2020), through the
planner. The check itself (``reproduction_usman2020``) runs in test_nuclear_validation.py; these test its parts."""

from pathlib import Path

import numpy as np
import pytest

from physim.nuclear import Experiment, reproduction
from physim.nuclear.exit_channels import SOURCE

REFERENCE = Path(__file__).resolve().parents[2] / "tests" / "reference" / "nuclear"
needs_reference = pytest.mark.skipif(not REFERENCE.is_dir(), reason="reference files are not shipped with the wheel")


def test_the_example_is_the_first_foil_with_the_paper_beam():
    exp = Experiment.example(reproduction.EXAMPLE)
    assert exp.beam.particles_per_second == pytest.approx(reproduction.BEAM_PER_S, rel=1e-3)
    assert exp.run.beam_time.value == 2 and exp.run.beam_time.unit == "h"
    assert all(ch.source == SOURCE for ch in exp.fusion.channels)
    # The beam enters the first foil so that its mean energy there is the paper's 49.6 MeV.
    assert reproduction.mean_energy(exp.beam.energy_mev) == pytest.approx(49.6, abs=0.01)


def test_entry_energy_inverts_the_mean_energy_in_the_foil():
    for e in (14.9, 30.0, 49.6):
        assert reproduction.mean_energy(reproduction.entry_energy(e)) == pytest.approx(e, abs=1e-6)


def test_the_peak_of_a_parabola_is_its_vertex():
    e = np.linspace(40, 50, 11)
    s = 1000 - 30 * (e - 46.3) ** 2
    pe, ps, edge = reproduction.peak(e, s)
    assert (pe, ps, edge) == (pytest.approx(46.3), pytest.approx(1000), False)
    assert reproduction.peak(e, e)[2]  # still rising at the end: flagged


def test_the_stack_has_26_ho_foils_between_ti_and_cu():
    assert reproduction.STACK.count("Ho") == 26
    assert reproduction.STACK[:5] == ("Ti", "Ho", "Cu", "Ho", "Ti")


@needs_reference
def test_physim_stopping_through_the_stack_gives_the_paper_energies():
    """From 50.4 MeV through the inferred stack, physim's mean energy in each Ho foil lies within the paper's energy
    spread (EN-RSL-HW) of the tabulated one."""
    rows = [r for r in reproduction.measured(REFERENCE) if r["ref"] == reproduction.REF]
    spread = {float(r["e_lab_mev"]): float(r["e_lab_unc_mev"]) for r in rows}
    paper = reproduction.foil_energies(rows)
    ours = reproduction.stack_energies()
    assert len(ours) == len(paper) == 26
    for s, e in zip(ours, paper):
        assert abs(s["e_mean_mev"] - e) <= spread[e], (s["foil"], s["e_mean_mev"], e)


@needs_reference
def test_the_saved_reproduction_passes_its_comparisons():
    saved = reproduction.read_result(REFERENCE / reproduction.RESULT_FILE)
    c = reproduction.compare(saved, reproduction.measured(REFERENCE))
    assert c["passed"], reproduction.summary(c)
    assert {x["channel"] for x in c["peaks"] if x["judged"]} == {"2n", "3n", "4n"}
    # The documented limit: physim's σ(3n)/σ(4n) at 49.6 MeV is far below every measurement (pre-equilibrium).
    at50 = [x for x in c["ratio_3n_4n"] if x["e_lab_mev"] >= 49]
    assert at50 and all(x["physim"] < 0.5 * x["measured"] for x in at50)
