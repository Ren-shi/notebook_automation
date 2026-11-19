"""The recoil's path for a lifetime measurement (rates.recoil_stack; lifetimes wave 2 base)."""

import numpy as np
import pytest

from physim.nuclear import Experiment, transport
from physim.nuclear.rates import recoil_stack, stack


def names(segments):
    return [s.name for s in segments]


def test_a_plunger_adds_the_gap_with_its_offset_and_the_stopper():
    exp = Experiment.example("coulex_ru98_plunger")
    seg = recoil_stack(exp, 40.0)
    assert names(seg) == ["target", "gap", "stopper"]
    assert seg[0].thickness_mg_cm2 == pytest.approx(0.7)
    assert seg[1].material is None and seg[1].gap_um == pytest.approx(48.0)  # 40 um plus the 8 um offset
    assert seg[2].thickness_mg_cm2 == pytest.approx(15.7)
    assert seg[2].material is not None
    # The beam's path does not change.
    assert [lay.name for lay in stack(exp)] == ["target"]


def test_a_differential_plunger_has_two_gaps_and_two_foils():
    exp = Experiment.example("fusion_cr50_plunger")
    seg = recoil_stack(exp, 42.0)
    assert names(seg) == ["target", "gap", "degrader", "gap2", "degrader2"]
    assert seg[1].gap_um == pytest.approx(42.0)
    assert seg[2].thickness_mg_cm2 == pytest.approx(4.4)
    assert seg[3].gap_um == pytest.approx(100.0)
    assert seg[4].thickness_mg_cm2 == pytest.approx(6.5)


def test_dsam_and_no_lifetime_are_the_target_stack():
    exp = Experiment.example("fusion_sm143_dsam")
    assert names(recoil_stack(exp)) == ["target", "backing"]
    assert names(recoil_stack(Experiment.example("coulex_ni58"))) == [lay.name for lay in stack(
        Experiment.example("coulex_ni58"))]


def test_a_plunger_needs_a_distance():
    exp = Experiment.example("coulex_ru98_plunger")
    with pytest.raises(ValueError, match="distance"):
        recoil_stack(exp)
    with pytest.raises(ValueError, match="negative"):
        recoil_stack(exp, -1.0)


def test_the_recoil_flies_the_gap_and_stops_in_the_stopper():
    # 98Ru at 4.5 % of c (about 93 MeV) made in the middle of the target: it crosses 48 um of vacuum in 3.6 ps and
    # stops in the copper.
    exp = Experiment.example("coulex_ru98_plunger")
    seg = recoil_stack(exp, 40.0)
    mass = 98 * 931.494
    energy = mass * (1 / np.sqrt(1 - 0.045**2) - 1)
    late = transport.slow(seg, "98Ru", mass, energy, [0.0, 0.0, 1.0], 0.35, 1e4)
    assert bool(late.stopped[0]) and int(late.segment[0]) == 2
    mid = transport.slow(seg, "98Ru", mass, energy, [0.0, 0.0, 1.0], 0.35, 2.0)
    assert int(mid.segment[0]) == 1 and not bool(mid.stopped[0])
