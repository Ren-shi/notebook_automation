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

GOSIA itself is not here: :meth:`Multistep.gosia_input` writes an input (checked against runs of GOSIA 20110524)
that integrates the γ yields over the target and each particle detector, so that GOSIA can check these numbers;
:func:`physim.external.gosia.integrated` reads its integrated yields and Rutherford cross sections back.
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
from .rates import Rates, _q, beam_energy_at, beam_ion, stack, stopping


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

    def gosia_experiments(self, mesh_theta: int = 9) -> list:
        """The GOSIA experiments :meth:`gosia_input` writes, in order: one for each particle detector and each
        nucleus that can reach it (the scattered beam, the target recoil, or both, as :meth:`yields` counts them),
        and for a beam detected below its largest angle in inverse kinematics one for each branch of the
        kinematics. Each is a dict: ``detector``, ``particle`` ("ejectile" for the scattered beam, "recoil" for the
        target nucleus), ``ikin``, ``theta`` (the (lowest, highest) lab angle of the detected nucleus, degrees,
        clipped below the largest angle it reaches), ``theta_expt`` (the EXPT line's angle: negative for a detected
        recoil, and then the beam's lab angle, as GOSIA wants), ``phi`` (None for the full azimuth, else the
        azimuth ranges (degrees) at each angle of the mesh) and ``notes``."""
        a_p, a_t = data.nuclide(self.beam).A, data.nuclide(self.target).A
        e_in, e_out = self._gosia_energy_range()
        e_x = self._gosia_ncm_energy()
        out = []
        for g in self.array:
            lo, hi = g.theta_range()
            for particle in ("ejectile", "recoil"):
                top = _gosia_max_angle(e_out, e_x, a_p, a_t, particle)
                if lo >= top:
                    continue
                notes = []
                t_lo, t_hi = lo, (hi if hi < top else top - 0.01 * (top - lo))
                if hi >= top:
                    notes.append(f"{g.name}: the {'beam' if particle == 'ejectile' else 'recoil'} reaches only "
                                 f"{top:.2f}° at the target's exit energy; its range is cut at {t_hi:.2f}°.")
                if particle == "ejectile" and t_lo < 0.5:
                    notes.append(f"{g.name}: the beam's Rutherford cross section diverges at 0°; the range starts "
                                 "at 0.5°.")
                    t_lo = 0.5
                thetas = np.linspace(t_lo, t_hi, mesh_theta)
                phi = _phi_intervals(g, thetas)
                if all(len(p) == 1 and p[0][1] - p[0][0] >= 359.9 for p in phi):
                    phi = None
                mean = float(np.clip(g.mean_theta(), t_lo, t_hi))
                if particle == "recoil":
                    tp, ikin = _gosia_invkin(e_in, e_x, a_p, a_t, mean)
                    out.append(dict(detector=g.name, particle="recoil", ikin=ikin, theta=(t_lo, t_hi),
                                    theta_expt=-tp, phi=phi, notes=notes))
                    continue
                tau = _gosia_tau(e_out, e_x, a_p, a_t)
                branches = (1, 0) if tau > 1.0 else (1,)
                for ikin in branches:
                    out.append(dict(detector=g.name, particle="ejectile", ikin=ikin, theta=(t_lo, t_hi),
                                    theta_expt=mean, phi=phi, notes=list(notes)))
        return out

    def _gosia_energy_range(self) -> tuple:
        """Beam energy (MeV) entering and leaving the target layer."""
        th = self.layers[0].thickness
        e_in, e_out = (float(x) for x in beam_energy_at(self.experiment, 0, np.array([0.0, th]), self.layers))
        if e_in - e_out < 1e-3:
            e_out = e_in - 1e-3
        return e_in, e_out

    def _gosia_ncm_energy(self) -> float:
        """Energy (MeV) of GOSIA's NCM level, whose excitation sets the kinematics: the first excited level."""
        s = self.scheme
        return (s.levels[1].energy.value - s.levels[0].energy.value) * 1e-3 if len(s.levels) > 1 else 0.0

    def gosia_input(self, title: Optional[str] = None, gamma: str = "array", mesh: tuple = (9, 9),
                    integration: tuple = (40, 40), stopping_points: int = 20, accuracy: int = 7,
                    substates: tuple = (5, 5), dipole: float = 3.8) -> str:
        """A GOSIA (20110524) input that integrates the γ yields of this setup over the target's thickness and
        each particle detector, to check :meth:`yields` against GOSIA.

        - OP,GOSI with LEVE and ME from the scheme (levels numbered from 1; matrix elements in e b^(λ/2)).
        - EXPT: the investigated nucleus is the scheme's. GOSIA's sign rule (its CMLAB and INVKIN): a *negative*
          partner Z means the investigated nucleus is the *projectile*, a positive one the target. So the beam's
          scheme (``role="beam"``) writes the target's Z negative, the target's scheme the beam's Z positive. A
          negative angle means the target recoil is detected, with the beam's lab angle given (the mesh of OP,INTI
          then holds recoil angles). One experiment per detector and per nucleus that reaches it
          (:meth:`gosia_experiments`); IAX = 0 for a ring around the beam, else the azimuths covered at each angle.
        - CONT: ``ACC``; ``DIP`` (the E1 polarisation; 3.8 is GOSIA's own default, written out so that it does not
          warn); ``SPL,1.`` (splines: GOSIA's default, one polynomial through the whole stopping table, puts the
          integrated Rutherford of a thick target well off, 17 % low on ⁶²Ni + C); ``END,`` and the blank line
          GOSIA needs after it.
        - OP,YIEL: no internal conversion. ``gamma="array"`` places one GOSIA detector on the axis of each of the
          setup's γ detectors, sized by OP,GDET (written first, so the file runs as it is) as a cylinder of the
          crystals' front area and length; ``gamma="4pi"`` uses one detector with all attenuation coefficients
          zero, so its yield ×4π is the angle-integrated γ yield. GOSIA reads that detector from ``fort.9``:
          write :data:`GOSIA_4PI_FORT9` there beside the input.
        - OP,INTI: the beam energy from the target's exit to its entry (``mesh[0]`` points), the detected nucleus's
          lab angles (``mesh[1]`` points), the beam's stopping power in the target from the app's tables
          (``stopping_points``), and ``integration`` steps in energy and angle.

        GOSIA's "integrated Rutherford cross section" uses the kinematics of the inelastic scattering (its orbit
        of the NCM level), not the elastic: a ratio γ yield / Rutherford built from its output divides by that
        one (11–16 % above the elastic integral for ⁶²Ni on carbon at 14–44° recoil angles)."""
        if gamma not in ("array", "4pi"):
            raise ValueError("gamma must be 'array' or '4pi'")
        exp = self.experiment
        s = self.scheme
        lines = ["OP,TITL", title or exp.title]
        gdets = list(exp.gamma_detectors) if gamma == "array" else []
        if gamma == "array" and not gdets:
            raise ValueError("the setup has no γ detectors; use gamma='4pi'")
        if len(gdets) > 32:
            raise ValueError("GOSIA takes at most 32 γ detectors per experiment")
        if gdets:
            lines += ["OP,GDET", str(len(gdets))]
            for gd in gdets:
                n = gd.crystals or 1
                radius_cm = gd.crystal_radius_mm() * math.sqrt(n) / 10.0
                length_cm = _q(gd.crystal_length).to("cm") if gd.crystal_length is not None else 7.0
                lines.append(f"0.01 {radius_cm:.4f} {length_cm:.4f} {gd.distance_mm() / 10.0:.4f}")
                lines.append("0 0 0 0 0 0 0")
        lines += ["OP,GOSI", "LEVE"]
        for i, lev in enumerate(s.levels, start=1):
            parity = "+1" if (lev.parity or 1) > 0 else "-1"
            lines.append(f"{i} {parity} {lev.spin if lev.spin is not None else 0:g} {lev.energy.value * 1e-3:.6f}")
        lines.append("0 0 0 0")
        lines.append("ME")
        lams = []
        for mult in ("E1", "E2", "E3"):
            mes = [me for me in s.matrix_elements if me.multipolarity == mult]
            if not mes:
                continue
            lam = int(mult[1])
            lams.append(lam)
            lines.append(f"{lam} 0 0 0 0")
            for me in sorted(mes, key=lambda m: (m.a, m.b)):
                # ⟨a‖Eλ‖b⟩ in e fm^λ → e b^(λ/2): /10, /100, /1000 for λ = 1, 2, 3.
                value = me.value.value / 10.0 ** lam
                lim = 5 * abs(value) + 0.1
                lines.append(f"{me.a + 1} {me.b + 1} {value:.6f} {-lim:.5f} {lim:.5f}")
        lines.append("0 0 0 0 0")

        expts = self.gosia_experiments(mesh[1])
        if not expts:
            raise ValueError("no particle detector sees the beam or the target recoil")
        a_p, a_t = data.nuclide(self.beam).A, data.nuclide(self.target).A
        z_p, z_t = data.nuclide(self.beam).Z, data.nuclide(self.target).Z
        if self.excite == "projectile":
            inv, partner = (z_p, a_p), (-z_t, a_t)        # negative partner Z: the projectile is investigated
        else:
            inv, partner = (z_t, a_t), (z_p, a_p)         # positive partner Z: the target is investigated
        e_in, e_out = self._gosia_energy_range()
        lines += ["EXPT", f"{len(expts)} {inv[0]} {inv[1]}"]
        dets = [g.name for g in self.array]
        for x in expts:
            iax, phi1, phi2 = (0, 0.0, 360.0) if x["phi"] is None else (1, *_phi_span(x["phi"]))
            lines.append(f"{partner[0]} {partner[1]} {e_in:.4f} {x['theta_expt']:.4f} {substates[0]} "
                         f"{substates[1]} {iax} {phi1:.2f} {phi2:.2f} {x['ikin']} {dets.index(x['detector']) + 1}")
        # GOSIA needs the blank line after END,
        lines += ["CONT", f"ACC,{accuracy}.", f"DIP,{dipole:g}", "SPL,1.", "END,", ""]

        # OP,YIEL: no internal conversion; the γ detectors at each experiment; normalised to the first transition.
        offd = sorted((me.b, me.a) for me in s.matrix_elements if me.a != me.b)
        trans = offd or [(1, 0)]
        e_top = max([lev.energy.value * 1e-3 for lev in s.levels] + [0.1])
        energies = [0.05, e_top / 2, e_top + 0.5]
        if gdets:
            ids = list(range(1, len(gdets) + 1))
            angles = [(_q(gd.theta).to("deg"), _q(gd.phi).to("deg") if gd.phi is not None else 0.0) for gd in gdets]
        else:
            ids, angles = [1], [(90.0, 0.0)]
        lines += ["OP,YIEL", "0", f"{len(energies)} {len(lams)}", " ".join(f"{e:.4f}" for e in energies)]
        for lam in lams:
            lines += [str(lam), " ".join("0.0" for _ in energies)]
        lines.append(" ".join(str(len(ids)) for _ in expts))
        for _ in expts:
            lines += [" ".join(map(str, ids)), " ".join(f"{a[0]:.2f}" for a in angles),
                      " ".join(f"{a[1]:.2f}" for a in angles)]
        lines.append(f"{trans[0][0] + 1} {trans[0][1] + 1}")
        for _ in expts:
            lines += [str(len(ids)), " ".join("1.0" for _ in ids), " ".join("1.0" for _ in ids)]
        lines += ["0", "0 0", "0 0", "0 0", "0 0"]

        # OP,INTI: the meshes of all experiments, then the stopping powers and the integration steps.
        lines.append("OP,INTI")
        e_mesh = np.linspace(e_out, e_in, mesh[0])
        for x in expts:
            th = np.linspace(*x["theta"], mesh[1])
            nt = mesh[1] if x["phi"] is None else -mesh[1]
            lines.append(f"{mesh[0]} {nt} {e_out:.4f} {e_in:.4f} {x['theta'][0]:.4f} {x['theta'][1]:.4f}")
            lines.append(" ".join(f"{e:.4f}" for e in e_mesh))
            lines.append(" ".join(f"{t:.4f}" for t in th))
            if x["phi"] is not None:
                # GOSIA divides each mesh point's yield by the azimuth it covers; where the face only touches an
                # angle (a corner), a sliver at the face's middle azimuth stands in. The widths that weigh the
                # integration (below) stay exact.
                mid = sum(_phi_span(x["phi"])) / 2
                for ivs in x["phi"]:
                    ivs = [(a, b) for a, b in ivs if b - a > 0.05] or [(mid - 0.05, mid + 0.05)]
                    lines += [str(len(ivs)), " ".join(f"{a:.3f} {b:.3f}" for a, b in ivs)]
        st = stopping(self.beam, self.layers[0].material)
        pad = 0.05 * (e_in - e_out) + 0.1
        e_sp = np.linspace(max(e_out - pad, 0.01), e_in + pad, stopping_points)
        npce, npct = (int(v) + int(v) % 2 for v in integration)
        for k, x in enumerate(expts):
            if k == 0:
                lines += [str(stopping_points), " ".join(f"{e:.4f}" for e in e_sp),
                          " ".join(f"{float(st.stopping_power(e)):.5f}" for e in e_sp)]
            else:
                lines.append("0")
            if x["phi"] is None:
                lines.append(f"{npce} {npct}")
            else:
                lines.append(f"{npce} {-npct}")
                t_int = np.linspace(*x["theta"], npct + 1)
                widths = [sum(b - a for a, b in ivs) for ivs in _phi_intervals(self.array[x["detector"]], t_int)]
                lines.append(" ".join(f"{w:.3f}" for w in widths))
        lines.append("OP,EXIT")
        return "\n".join(lines) + "\n"


#: GOSIA's γ-detector file (``fort.9``, normally written by OP,GDET) for one detector whose attenuation
#: coefficients Q_k are zero for every k ≥ 1: the yield it gives, times 4π, is the angle-integrated γ yield. Put it
#: beside an input written with ``gosia_input(gamma="4pi")``.
GOSIA_4PI_FORT9 = "1\n1.0\n0.0\n" + "0.0 1.0 0.0\n" * 8


def _gosia_tau(e: float, e_x: float, a_p: float, a_t: float) -> float:
    """GOSIA's τ (CMLAB): (A_p/A_t)·sqrt(E/(E − E_x(1 + A_p/A_t))), non-relativistic with integer masses."""
    return a_p / a_t * math.sqrt(e / (e - e_x * (1 + a_p / a_t)))


def _gosia_max_angle(e: float, e_x: float, a_p: float, a_t: float, particle: str) -> float:
    """Largest lab angle (degrees) of the scattered beam or of the target recoil in GOSIA's kinematics."""
    taup = math.sqrt(e / (e - e_x * (1 + a_p / a_t)))
    if particle == "ejectile":
        tau = taup * a_p / a_t
        return math.degrees(math.asin(1.0 / tau)) if tau > 1.0 else 180.0
    # Recoil: tan θ_r = sin z / (cos z + τ'), largest where cos z = −1/τ'.
    return math.degrees(math.atan(1.0 / math.sqrt(taup * taup - 1.0))) if taup > 1.0 else 90.0


def _gosia_invkin(e: float, e_x: float, a_p: float, a_t: float, theta_r: float) -> tuple:
    """GOSIA's INVKIN: the beam's lab angle (degrees) for a target recoil at ``theta_r`` (the larger of the two
    roots, as OP,INTI takes it) and the IKIN branch it lies on."""
    taup = math.sqrt(e / (e - e_x * (1 + a_p / a_t)))
    tau = taup * a_p / a_t
    y = math.tan(math.radians(theta_r))
    t = math.sqrt(taup**2 * y**4 - (1 + y * y) * (taup**2 * y * y - 1))
    xs = [(taup * y * y + sg * t) / (1 + y * y) for sg in (1, -1)]
    tp = max(math.atan2(math.sqrt(max(1 - x * x, 0.0)), tau + x) for x in xs)
    if tau * math.sin(tp) > 1.0:
        tp = math.asin(1.0 / tau)
    if tau < 1.0:
        return math.degrees(tp), 1
    tc = math.acos(-1.0 / tau)
    thres = math.sin(tc) / (taup - math.cos(tc))
    return math.degrees(tp), 1 if abs(y) > abs(thres) else 0


def _phi_intervals(g, thetas, step: float = 0.25) -> list:
    """For each lab angle, the azimuth ranges (degrees, within [0, 360], merged across 0°) the detector face
    covers there."""
    phis = np.arange(0.0, 360.0, step) + step / 2
    lo, hi = g.theta_range()
    eps = 1e-4 * max(hi - lo, 1e-3)
    out = []
    for th in np.atleast_1d(thetas):
        # The ends of the mesh lie on the face's edges: look just inside them.
        t = math.radians(float(np.clip(th, lo + eps, hi - eps)))
        dirs = np.column_stack([np.sin(t) * np.cos(np.radians(phis)), np.sin(t) * np.sin(np.radians(phis)),
                                np.full(phis.shape, np.cos(t))])
        hit = np.array([h is not None for h in g.hits(dirs)])
        if hit.all():
            out.append([(0.0, 360.0)])
            continue
        ivs, start = [], None
        for k, h in enumerate(hit):
            if h and start is None:
                start = phis[k] - step / 2
            if not h and start is not None:
                ivs.append((start, phis[k] - step / 2))
                start = None
        if start is not None:
            if ivs and ivs[0][0] <= 0.0:          # the run through 360° continues at 0°
                ivs[0] = (start - 360.0, ivs[0][1])
            else:
                ivs.append((start, 360.0))
        out.append(ivs)
    return out


def _phi_span(phi: list) -> tuple:
    """The overall azimuth range (degrees) of a detector's mesh, for the EXPT line."""
    a = [iv[0] for ivs in phi for iv in ivs]
    b = [iv[1] for ivs in phi for iv in ivs]
    return (min(a), max(b)) if a else (0.0, 360.0)


__all__ = ["GOSIA_4PI_FORT9", "Multistep", "Yields"]
