"""Multi-step Coulomb excitation of the planned experiment: the yields each ring would see with all orders,
the prolate–zero–oblate comparison that the reorientation effect allows, a fit of a few matrix elements to
measured yields, and a GOSIA input file for the setup.

::

    from physim.nuclear.multistep import Multistep

    ms = Multistep(experiment, scheme, role="target")
    y = ms.yields()                       # γ rays per second of each transition in each particle detector
    ms.shapes()                           # the same with Q(2⁺) prolate, zero and oblate, against the uncertainties
    ms.fit({"CD": {(1, 0): (area, unc)}}, free=[(1, 1, "E2")])
    ms.gosia_input()                      # text of a GOSIA input file for this setup

The excitation probabilities come from :class:`~physim.nuclear.coupled.CoupledChannels` on a grid of CM angles
(``angle_step``) at a few beam energies through the target (``energies``, Gauss–Legendre), interpolated to each
quadrature direction of each detector segment, weighted by the Rutherford cross section of the orbit, and the
decay of :mod:`physim.nuclear.orientation` shares the populations out into γ rays.

GOSIA itself is not here: the input file is written from the published format so that GOSIA, run elsewhere,
can check these numbers or fit the matrix elements from measured yields.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from . import data
from .coupled import CoupledChannels
from .detectors import Array
from .levels import LevelScheme, quadrupole_factor
from .orientation import Excitation
from .rates import Rates, _q, beam_energy_at, beam_ion, stack


def _target_name(experiment) -> str:
    layers = stack(experiment)
    return data.nuclide(max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]).name


@dataclass
class Yields:
    """γ rays per second of each transition in coincidence with each particle detector (and each of its rings
    or strips), from a solution of the coupled equations."""

    #: Rate per second by detector name, then by transition (a pair of level indices).
    detectors: dict
    #: The same by detector, then ring or strip number, then transition.
    rings: dict
    #: Excitations per second by detector, then level (populated, with feeding).
    levels: dict
    #: The transitions followed and the levels, for labels.
    transitions: list
    notes: list = field(default_factory=list)

    def total(self, transition: tuple) -> float:
        return sum(d.get(transition, 0.0) for d in self.detectors.values())


class Multistep:
    """The planned experiment with the coupled equations. ``role`` says which nucleus the scheme belongs to
    ("target" or "beam"); ``first_order`` solves with first-order theory instead, for comparison."""

    def __init__(self, experiment, scheme: LevelScheme, role: str = "target", angle_step: float = 4.0,
                 energies: int = 2, tolerance: float = 1e-4, first_order: bool = False, order: int = 4):
        self.experiment = experiment
        self.scheme = scheme
        self.role = role
        self.angle_step = angle_step
        self.n_energies = energies
        self.tolerance = tolerance
        self.first_order = first_order
        self.order = order
        self.layers = stack(experiment)
        self.beam = beam_ion(experiment)
        self.target = _target_name(experiment)
        self.excite = "target" if role == "target" else "projectile"
        self.array = Array.from_experiment(experiment)
        self.notes: list = []

    # -- the solution on a grid -------------------------------------------------------------------------------------

    def _energies(self) -> tuple:
        """Beam energies through the target and their weights (Gauss–Legendre over the thickness)."""
        x, w = np.polynomial.legendre.leggauss(self.n_energies)
        th = self.layers[0].thickness
        z = (x + 1) / 2 * th
        e = beam_energy_at(self.experiment, 0, z, self.layers)
        return np.asarray(e, dtype=float), w / 2

    def solver(self, energy: float, scheme: Optional[LevelScheme] = None):
        scheme = scheme or self.scheme
        if self.first_order:
            return Excitation(self.beam, self.target, energy, scheme, excite=self.excite)
        return CoupledChannels(self.beam, self.target, energy, scheme, excite=self.excite, tolerance=self.tolerance)

    def yields(self, scheme: Optional[LevelScheme] = None) -> Yields:
        """γ rays per second of every transition in every particle detector and ring."""
        scheme = scheme or self.scheme
        exp = self.experiment
        energies, weights = self._energies()
        pps = exp.beam.particles_per_second
        atoms = sum(ch.atoms_per_cm2 for ch in Rates(exp).channels if ch.excitation is None
                    and ch.nuclide.name == self.target) or self.layers[0].nuclides.get(self.target, 0.0)
        from .kinematics import TwoBody

        dets, rings, levels = {}, {}, {}
        transitions = None
        for e, we in zip(energies, weights):
            ex = self.solver(float(e), scheme)
            if transitions is None:
                transitions = ex.transitions()
                self.notes += ex.notes
            grid, states = ex.table(self.angle_step)
            cx = ex.cx if hasattr(ex, "cx") else ex.paths[min(ex.paths)][0][2]
            ruth = np.asarray(cx.rutherford_cm(grid)) * 1e-27  # cm²/sr
            yields_grid = {t: np.array([s.gamma_yield(*t) for s in states]) for t in transitions}
            pops_grid = {n: np.array([s.population(n) for s in states]) for n in range(len(ex.scheme.levels))}
            # The kinematics of the reference level's excitation map each lab direction to the orbit's angle.
            ref = ex.reference if hasattr(ex, "reference") else min(ex.paths)
            e_ref = (ex.scheme.levels[ref].energy.value - ex.scheme.levels[ex.ground].energy.value) * 1e-3
            tb = TwoBody(self.beam, self.target, float(e), excitation_mev=e_ref,
                         excite="recoil" if self.excite == "target" else "ejectile")
            for g in self.array:
                dets.setdefault(g.name, {})
                rings.setdefault(g.name, {})
                levels.setdefault(g.name, {})
                for seg in g.segments:
                    dirs, dom = g.directions(seg, self.order)
                    theta = np.degrees(np.arccos(np.clip(dirs[:, 2], -1, 1)))
                    for particle in ("ejectile", "recoil"):
                      for point in tb.at_lab(theta, particle):
                        if point is None:
                            continue
                        th_cm = np.asarray(point.theta_cm, dtype=float)
                        jac = np.nan_to_num(np.asarray(point.jacobian, dtype=float))
                        ok = ~np.isnan(th_cm) & (jac > 0)
                        if not ok.any():
                            continue
                        th_ej = th_cm[ok] if particle == "ejectile" else 180.0 - th_cm[ok]
                        weight = pps * atoms * np.interp(th_ej, grid, ruth) * jac[ok] * dom[ok] * we
                        ring = rings[g.name].setdefault(seg[0], {})
                        for t, yg in yields_grid.items():
                            v = float(np.sum(weight * np.interp(th_ej, grid, yg)))
                            dets[g.name][t] = dets[g.name].get(t, 0.0) + v
                            ring[t] = ring.get(t, 0.0) + v
                        for n, pg in pops_grid.items():
                            v = float(np.sum(weight * np.interp(th_ej, grid, pg)))
                            levels[g.name][n] = levels[g.name].get(n, 0.0) + v
        return Yields(dets, rings, levels, transitions or [], list(self.notes))

    # -- the shape ------------------------------------------------------------------------------------------------

    def shapes(self, level: int = 1, beam_time_s: Optional[float] = None) -> dict:
        """The 2⁺ state's γ yields in each ring with its quadrupole moment at the prolate rotor value, zero and
        the oblate rotor value, and whether the three differ by more than the counting uncertainties in the
        planned beam time: {"cases": {"prolate": Yields, ...}, "q_efm2": {...}, "rings": [...], "separable": bool}."""
        scheme = self.scheme
        b_up = scheme.b(0, level, "E2")
        if b_up is None:
            raise ValueError("the scheme needs an E2 matrix element from the ground state to the level")
        q0 = math.sqrt(16 * math.pi * b_up / 5)
        qs = -2 / 7 * q0
        j = float(scheme.levels[level].spin)
        factor = quadrupole_factor(j)
        t = beam_time_s if beam_time_s is not None else _q(self.experiment.run.beam_time).to("s")
        cases, q_values = {}, {}
        import copy

        for name, q in (("prolate", qs), ("spherical", 0.0), ("oblate", -qs)):
            s = copy.deepcopy(scheme)
            if q:
                s.set_matrix_element(level, level, "E2", q / factor, source="assumed",
                                     note=f"rigid rotor, {name}")
            else:
                s.matrix_elements = [me for me in s.matrix_elements if not (me.a == level and me.b == level)]
            cases[name] = self.yields(s)
            q_values[name] = q
        transition = next((tr for tr in cases["spherical"].transitions if tr[0] == level), None)
        rows = []
        separable = False
        for det in cases["spherical"].detectors:
            for ring in sorted(cases["spherical"].rings[det]):
                counts = {n: cases[n].rings[det][ring].get(transition, 0.0) * t for n in cases}
                unc = math.sqrt(max(counts["spherical"], 1.0))
                diff = abs(counts["prolate"] - counts["oblate"])
                rows.append({"detector": det, "ring": ring, "counts": counts, "uncertainty": unc,
                             "difference_sigma": diff / unc if unc > 0 else 0.0})
                separable = separable or diff > 3 * unc
        total = {n: sum(r["counts"][n] for r in rows) for n in cases}
        total_unc = math.sqrt(max(total["spherical"], 1.0))
        return {"cases": cases, "q_efm2": q_values, "transition": transition, "rings": rows, "totals": total,
                "separable": separable or abs(total["prolate"] - total["oblate"]) > 3 * total_unc,
                "total_difference_sigma": abs(total["prolate"] - total["oblate"]) / total_unc}

    # -- the fit ----------------------------------------------------------------------------------------------------

    def fit(self, measured: dict, free: list, beam_time_s: Optional[float] = None, iterations: int = 8) -> dict:
        """Adjust the ``free`` matrix elements [(a, b, "E2"), ...] (at most three) until the calculated counts
        match ``measured``: {detector: {(initial, final): (counts, uncertainty)}} in the planned beam time.
        Gauss–Newton least squares with numerical derivatives; the uncertainties come from the curvature at
        the minimum. Returns {"values": {...}, "uncertainties": {...}, "chi2": ..., "iterations": n}."""
        if not 1 <= len(free) <= 3:
            raise ValueError("fit one to three matrix elements")
        import copy

        t = beam_time_s if beam_time_s is not None else _q(self.experiment.run.beam_time).to("s")
        keys = [(d, tr) for d, trs in measured.items() for tr in trs]
        y = np.array([measured[d][tr][0] for d, tr in keys], dtype=float)
        sig = np.array([max(measured[d][tr][1], 1e-9) for d, tr in keys], dtype=float)
        scheme = copy.deepcopy(self.scheme)
        x = np.array([scheme.matrix_element(a, b, m).value for a, b, m in free], dtype=float)

        def model(params: np.ndarray) -> np.ndarray:
            s = copy.deepcopy(scheme)
            for (a, b, m), v in zip(free, params):
                s.set_matrix_element(a, b, m, float(v), source="user")
            yl = self.yields(s)
            return np.array([yl.detectors.get(d, {}).get(tr, 0.0) * t for d, tr in keys])

        n_iter = 0
        f = model(x)
        chi2 = float(np.sum(((y - f) / sig) ** 2))
        for n_iter in range(1, iterations + 1):
            jac = np.empty((len(keys), len(x)))
            for k in range(len(x)):
                step = 0.02 * abs(x[k]) if x[k] else 1.0
                xp = x.copy()
                xp[k] += step
                jac[:, k] = (model(xp) - f) / step
            a = (jac / sig[:, None]).T @ (jac / sig[:, None])
            b = (jac / sig[:, None]).T @ ((y - f) / sig)
            try:
                dx = np.linalg.solve(a + 1e-12 * np.eye(len(x)), b)
            except np.linalg.LinAlgError:
                break
            x_new = x + dx
            f_new = model(x_new)
            chi2_new = float(np.sum(((y - f_new) / sig) ** 2))
            if chi2_new > chi2:  # halve the step until it helps
                for _ in range(4):
                    dx /= 2
                    x_new = x + dx
                    f_new = model(x_new)
                    chi2_new = float(np.sum(((y - f_new) / sig) ** 2))
                    if chi2_new <= chi2:
                        break
            x, f, done = x_new, f_new, abs(chi2 - chi2_new) < 1e-3 * max(chi2, 1.0)
            chi2 = chi2_new
            if done:
                break
        cov = np.linalg.inv((jac / sig[:, None]).T @ (jac / sig[:, None]) + 1e-12 * np.eye(len(x)))
        return {"values": {free[k]: float(x[k]) for k in range(len(x))},
                "uncertainties": {free[k]: float(math.sqrt(max(cov[k, k], 0.0))) for k in range(len(x))},
                "chi2": chi2, "degrees_of_freedom": len(keys) - len(x), "iterations": n_iter,
                "calculated": {key: float(v) for key, v in zip(keys, f)}}

    # -- GOSIA ------------------------------------------------------------------------------------------------------

    def gosia_input(self, yields: Optional[Yields] = None, title: Optional[str] = None) -> str:
        """A GOSIA input file for this setup, in the format of the GOSIA manual (OP.TITL, OP.GOSI with LEVE,
        ME, EXPT and CONT, then OP.YIEL and OP.INTI). The levels and matrix elements are the scheme's (with
        GOSIA's level numbering from 1 and its units, e b^λ/2 ... given as e fm^λ converted to eb and e b²), each
        experiment is one particle detector at its mean angle, and ``yields`` (if given) are written after
        OP.YIEL as the measured yields to fit. It is written from the manual, not checked against a run, since
        GOSIA is not on this machine."""
        exp = self.experiment
        s = self.scheme
        lines = ["OP,TITL", title or exp.title, "OP,GOSI", "LEVE"]
        for i, lev in enumerate(s.levels, start=1):
            parity = "+1" if (lev.parity or 1) > 0 else "-1"
            lines.append(f"{i} {parity} {lev.spin if lev.spin is not None else 0:g} {lev.energy.value * 1e-3:.6f}")
        lines.append("0 0 0 0")
        lines.append("ME")
        for mult in ("E1", "E2", "E3"):
            mes = [me for me in s.matrix_elements if me.multipolarity == mult]
            if not mes:
                continue
            lam = int(mult[1])
            lines.append(f"{lam} 0 0 0 0")
            unit = 100.0 ** lam / 1e4 if lam == 2 else (10.0 if lam == 1 else 1e6)  # e fm^λ → e b^(λ/2)... see note
            for me in sorted(mes, key=lambda m: (m.a, m.b)):
                value = me.value.value / {1: 10.0, 2: 100.0, 3: 1000.0}[lam]
                # GOSIA takes ⟨a‖M‖b⟩ in e b^(λ/2): e fm → e b^½ is /10, e fm² → e b is /100, e fm³ → /1000.
                lines.append(f"{me.a + 1} {me.b + 1} {value:.5f} {-5 * abs(value) - 0.1:.5f} {5 * abs(value) + 0.1:.5f}")
            del unit
        lines.append("0 0 0 0 0")
        lines.append("EXPT")
        dets = list(self.array)
        z_p, a_p = exp.beam.Z, exp.beam.A
        nuc_t = data.nuclide(self.target)
        lines.append(f"{len(dets)} {nuc_t.Z if self.role == 'target' else z_p} {nuc_t.A if self.role == 'target' else a_p}")
        for g in dets:
            theta = g.mean_theta()
            sign = -1 if self.role == "target" else 1  # GOSIA: negative Z_p for target excitation (projectile detected)
            lines.append(f"{sign * (z_p if self.role == 'target' else nuc_t.Z)} {a_p if self.role == 'target' else nuc_t.A} "
                         f"{exp.beam.energy_mev:.3f} {theta:.2f} 1 0 0 0 360 0 1")
        lines.append("CONT")
        lines.append("END,")
        lines.append("0 0 0")
        lines.append("OP,YIEL")
        lines.append("0")
        lines.append("1 1")
        lines.append(f"{len(exp.gamma_detectors) or 1}")
        for gd in exp.gamma_detectors or [None]:
            theta = _q(gd.theta).to("deg") if gd is not None else 90.0
            phi = _q(gd.phi).to("deg") if gd is not None and gd.phi is not None else 0.0
            lines.append(f"{theta:.2f} {phi:.2f}")
        if yields is not None:
            lines.append("OP,INTI")
            for i, g in enumerate(dets, start=1):
                for (a, b), v in sorted(yields.detectors.get(g.name, {}).items()):
                    counts = v * _q(exp.run.beam_time).to("s")
                    lines.append(f"{i} {a + 1} {b + 1} {counts:.4g} {math.sqrt(max(counts, 1.0)):.4g}")
        lines.append("OP,EXIT")
        return "\n".join(lines) + "\n"


__all__ = ["Multistep", "Yields"]
