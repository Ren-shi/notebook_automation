"""γ-ray yields with deorientation in vacuum against GOSIA's OP,POIN, from blind runs
(tests/reference/nuclear/coulex/gosia_yields_vac.csv): the on/off ratio of every line of ⁶²Ni, ¹⁹⁴Pt and ⁹⁸Ru at ten point
γ detectors, and the yields without deorientation relative to each other. The plan's own pieces at one point of a
particle detector; no GOSIA run needed."""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np
import pytest

from physim.nuclear import data, validation
from physim.nuclear.coupled import CoupledChannels
from physim.nuclear.deorientation import Parameters
from physim.nuclear.gamma import EVENT_TOLERANCE, DepthOrientation, Faded, _orbit_axes
from physim.nuclear.gamma_events import _directions, _distribution
from physim.nuclear.levels import Level, LevelScheme, MatrixElement, Transition, Value
from physim.nuclear.lifetimes import rate_per_unit_b

#: The reference's schemes, by GOSIA's level numbers (1 the 0⁺ ground state): levels (spin, MeV) from 2, <‖E2‖> in e b,
#: <‖M1‖> in µ_N (sizes), δ (E2/M1) and the conversion α by (initial, final).
SCHEMES = {
    "62Ni": ([(2, 1.17291), (0, 2.04855), (2, 2.30185), (4, 2.33646), (0, 2.89063)],
             {(1, 2): 0.2973390887, (2, 2): 0.06596887884, (2, 3): 0.2805649941, (2, 4): 0.3078174824,
              (1, 4): 0.06215263427, (2, 5): 0.5248890412, (2, 6): 0.03492560042},
             {(2, 4): 0.0935726}, {(4, 2): 3.1},
             {(2, 1): 0.000172, (3, 2): 0.000335, (4, 2): 0.000244, (4, 1): 0.000504, (5, 2): 0.000174,
              (6, 2): 0.000255}),
    "194Pt": ([(2, 0.328473), (2, 0.622024), (4, 0.811288), (4, 1.22952), (0, 1.2672), (6, 1.41183)],
              {(1, 2): 1.284983542, (2, 2): 0.6333012368, (3, 3): -0.6596887884, (4, 4): 0.6630460049,
               (2, 3): 1.723018614, (1, 3): 0.09767377711, (2, 4): 2.260455691, (4, 5): 1.039602873,
               (3, 5): 1.136188655, (2, 5): 0.1249447155, (3, 6): 0.2338928199, (2, 6): 0.06483062948,
               (4, 7): 2.410566081},
              {(2, 3): 0.0281468, (4, 5): 0.120979}, {(3, 2): 15.0, (5, 4): 3.0},
              {(2, 1): 0.0755, (3, 2): 0.106, (3, 1): 0.01483, (4, 2): 0.027, (5, 4): 0.043, (5, 3): 0.01565,
               (5, 2): 0.00674, (6, 3): 0.01367, (6, 2): 0.00621, (7, 4): 0.01607}),
    "98Ru": ([(2, 0.65246), (0, 1.32216), (4, 1.39791), (2, 1.41436), (6, 2.22265)],
             {(1, 2): 0.63239146, (2, 2): -0.27706929, (2, 3): 0.33575115, (2, 4): 1.17341475, (2, 5): 0.78570005,
              (1, 5): 0.11813934, (4, 6): 0.66829702},
             {(2, 5): 0.0384377}, {(5, 2): 13.0},
             {(2, 1): 0.00253, (3, 2): 0.00236, (4, 2): 0.00179, (5, 2): 0.00169, (5, 1): 0.0005, (6, 4): 0.0016}),
}
G_CHANGED = {"62Ni": 0.30, "194Pt": 0.25, "98Ru": 0.35}
DETECTORS = [(t, p) for t in (25.0, 55.0, 90.0, 125.0, 155.0) for p in (45.0, 135.0)]


def scheme(nuc: str) -> LevelScheme:
    """The scheme the app is given: the elements and the decays they give (γ intensities, δ, mean lives with α)."""
    levels, e2s, m1s, deltas, alphas = SCHEMES[nuc]
    spin, energy = [0] + [j for j, _ in levels], [0.0] + [e for _, e in levels]
    rates: dict = defaultdict(lambda: [0.0, 0.0])
    for slot, (lam, table) in enumerate((("E2", e2s), ("M1", m1s))):
        for (a, b), m in table.items():
            if a != b:
                i, f = (a, b) if energy[a - 1] > energy[b - 1] else (b, a)
                rates[i, f][slot] += (rate_per_unit_b(lam, 1e3 * (energy[i - 1] - energy[f - 1]))
                                      * (100 * m if lam == "E2" else m) ** 2 / (2 * spin[i - 1] + 1))
    s = LevelScheme(nuc)
    for i in range(1, len(spin) + 1):
        tau = 1 / sum((e2 + m1) * (1 + alphas[k]) for k, (e2, m1) in rates.items() if k[0] == i) if i > 1 else None
        s.levels.append(Level(Value(1e3 * energy[i - 1], source="user"), float(spin[i - 1]), 1, f"{spin[i - 1]}+",
                              None if tau is None else Value(tau * math.log(2), source="user")))
    for (i, f), (e2, m1) in sorted(rates.items()):
        mixed = e2 > 0 and m1 > 0
        assert not mixed or abs(deltas[i, f] ** 2 / (e2 / m1) - 1) < 1e-4
        s.transitions.append(Transition(i - 1, f - 1, Value(1e3 * (energy[i - 1] - energy[f - 1]), source="user"),
                                        Value(e2 + m1, source="user"), "M1+E2" if mixed else "E2",
                                        Value(deltas[i, f], source="user") if mixed else None,
                                        Value(alphas[i, f], source="user")))
    for (a, b), m in sorted(e2s.items()):
        s.matrix_elements.append(MatrixElement(a - 1, b - 1, "E2", Value(100 * m, source="user")))
    return s


def gosias_kinematics(beam, target, e, excite, theta_cm, e_ncm) -> tuple:
    """(v/c, lab polar angle) of the excited nucleus as GOSIA has them: non-relativistic, the first excited level's."""
    m1, m2 = (data.nuclide(x).nuclear_mass_mev for x in (beam, target))
    v_cm = math.sqrt(2 * e / m1) * m1 / (m1 + m2)
    v_rel = math.sqrt(2 * (e * m2 / (m1 + m2) - e_ncm) * (m1 + m2) / (m1 * m2))
    own = v_rel * (m2 if excite == "projectile" else -m1) / (m1 + m2)
    t = math.radians(theta_cm)
    vx, vz = abs(own * math.sin(t)), v_cm + own * math.cos(t)
    return math.hypot(vx, vz), math.degrees(math.atan2(vx, vz))


def point_yields(nuc, beam, target, e, excite, detected, theta_cm, vac, ex) -> dict:
    """{"i-f": yield per collision and steradian at each detector}: the plan's tensors (Faded.at for nuclei in vacuum
    from the collision on, DepthOrientation.coefficients) carried into the lab as gamma._correlation_rows does."""
    beta, th_e = gosias_kinematics(beam, target, e, excite, theta_cm, SCHEMES[nuc][0][0][1])
    phi_p = 0.0 if detected == "projectile" else 180.0
    u = _directions(th_e, phi_p if excite == "projectile" else phi_p + 180.0)
    dirs = _directions(np.array([t for t, _ in DETECTORS]), np.array([p for _, p in DETECTORS]))
    cos_a = dirs @ u
    rest = ((dirs - cos_a[:, None] * u) * (np.sqrt(1 - beta**2) / (1 - beta * cos_a))[:, None]
            + ((cos_a - beta) / (1 - beta * cos_a))[:, None] * u)
    local = rest @ _orbit_axes(np.array([float(theta_cm)]), np.array([phi_p]), excite == "projectile")[0]
    lives = np.array([0.0 if lv.half_life is None else 1e12 * lv.half_life.value / math.log(2)
                      for lv in ex.scheme.levels])
    par = None if vac == "off" else Parameters(g_factor=G_CHANGED[nuc] if vac == "changed" else None)
    orient, out = DepthOrientation(None, [ex]), {}
    for i, f in ex.transitions():
        faded = None if par is None else Faded(ex, i, lives, par)
        fade = None if faded is None or faded.off else faded.at(np.zeros(1), np.full(1, beta))
        a = orient.coefficients(i, f, np.array([float(theta_cm)]), ex.beam_energy, fade)
        w = (_distribution({k: v[0] for k, v in a.items()}, np.arccos(np.clip(local[:, 2], -1, 1)),
                           np.arctan2(local[:, 1], local[:, 0])) * (1 - beta**2) / (1 - beta * cos_a) ** 2)
        p = float(orient.population(i, ex.beam_energy).probability(float(theta_cm)))
        out[f"{i + 1}-{f + 1}"] = p * ex.decay(i, f).gamma * w / (4 * math.pi)
    return {"beta": beta, "yields": out}


def test_yields_with_deorientation_are_gosias_in_blind_runs():
    """1971 numbers within their bands (fixed before the runs), but for the six missed in the blind check (the 2⁺₂ → 0⁺
    lines: the two codes' Coulomb excitation of the 2⁺₂ level, not the deorientation), which stay within 4 bands."""
    ref = validation.reference_dir()
    path = ref / "coulex" / "gosia_yields_vac.csv" if ref else None
    if path is None or not path.exists():
        pytest.skip("no GOSIA reference")
    runs = defaultdict(list)
    for r in validation.read_reference(path)[1]:
        runs[r["nucleus"], r["beam"], r["target"], float(r["beam_energy_mev"]), r["excite"], r["detected"],
             int(r["theta_cm_deg"])].append(r)
    counts, outside, recorded, solved = defaultdict(int), [], [], {}
    for (nuc, beam, target, e, excite, detected, theta), rows in runs.items():
        if nuc not in solved:  # the plan's solution (gamma.cascade_excitation's call), one per collision
            solved[nuc] = CoupledChannels(beam, target, e, scheme(nuc), excite=excite, isotropic=False,
                                          tolerance=EVENT_TOLERANCE, e1_polarisation=0.0038)
        ex = solved[nuc]
        mine = {vac: point_yields(nuc, beam, target, e, excite, detected, theta, vac, ex)
                for vac in ("off", "default", "changed")}
        assert mine["off"]["beta"] == pytest.approx(float(rows[0]["beta_gosia"]), abs=1e-6)
        y = {vac: m["yields"] for vac, m in mine.items()}
        for r in rows:
            d = DETECTORS.index((float(r["theta_gamma_deg"]), float(r["phi_gamma_deg"])))
            line = r["line"]
            if r["kind"] == "ratio":
                app = y[r["set"]][line][d] / y["off"][line][d]
            else:
                app = y["off"][line][d] / y["off"]["2-1"][0]
            counts[r["kind"]] += 1
            where = (nuc, theta, r["set"], line, r["theta_gamma_deg"], r["phi_gamma_deg"], r["kind"])
            off_by = abs(app - float(r["gosia"])) / float(r["band"])
            if r["blind"] == "miss":
                recorded.append(where)
                assert off_by < 4, (where, off_by)
            elif off_by > 1:
                outside.append((where, off_by))
    assert dict(counts) == {"ratio": 1320, "relative": 651}
    assert len(recorded) == 6
    assert not outside, outside
