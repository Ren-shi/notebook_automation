"""Particle–γ events: γ rays from the simulated Coulomb-excitation events, their detection, the Doppler correction
an experimentalist can make, coincidences and background.

Built on the particle events of :mod:`physim.nuclear.events`: every reaction that put a particle into a detector
and excited the state emits its γ ray, which is followed to the crystals::

    from physim.nuclear import Experiment
    from physim.nuclear.gamma_events import simulate_gammas

    g = simulate_gammas(Experiment.example("coulex_ni58"), events=200_000, seed=1)
    g.counts("Ge90")                            # particle–γ coincidences in the run, any particle detector
    g.spectrum("Ge90", corrected="recoil")      # the Doppler-corrected spectrum, counts per bin in the run
    g.coincidences()                            # particle × γ-detector matrix of true and random coincidences

**The chain, per excited event:**

1. the particle has reached a ring and sector (or strip) of a particle detector, as the generator decided;
2. the excited nucleus leaves the target with its velocity (its own path, or, if it was not the detected
   particle, that of the two-body partner of the one that was);
3. it emits its γ ray in its rest frame with the particle–γ correlation of :mod:`physim.nuclear.orientation`
   for that orbit; the direction and energy are transformed to the laboratory;
4. the γ ray that meets a crystal face interacts with the probability the response of
   :mod:`physim.nuclear.response` gives (after the absorbers), leaves its full energy or a Compton or escape
   deposit, and is recorded with the crystal's resolution; a threshold applies.

**The Doppler correction** uses only what the detectors know: the centre of the segment the particle hit, the
centre of the crystal, and two-body kinematics at the nominal beam energy. It is made twice: as if the emitter
were the scattered beam, and as if it were the target recoil. The width that remains comes from the sizes of the
segments and crystals.

**Background** in the γ spectra: the Compton continua of the response; random coincidences, from the singles
rates and the coincidence window (``[run] coincidence_window``); room-background lines (⁴⁰K and the thorium and
uranium series) at ``[run] room_background`` counts per second per crystal; and ``[run] extra_lines`` added by
hand. A simple non-paralysable dead time (``[run] dead_time``) scales every count by the live fraction.

**Add-back and Compton suppression** (clovers, :mod:`physim.nuclear.response`): a Compton deposit in a clover with
``addback`` is summed with what the scattered photon left in the neighbouring crystal, with the share the add-back
factor asks for, and the γ ray is then attributed to the crystal with the larger deposit, which the Doppler
correction uses. In a detector with a ``shield``, a deposit that is not the full energy is rejected with
probability 1 − 1/S(E). ``simulate_gammas(..., plain=True)`` leaves both out: the crystals as they are.

Everything here runs in NumPy on the particle events; the per-event loop of the generator stays in Rust.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from . import angular, data
from .detectors import Array
from .events import Events, simulate
from .gamma import cascade_excitation, excitation_of, excitation_order, other_excitation, start_levels
from .kinematics import TwoBody
from .quantity import Quantity
from .rates import (Rates, _after, _exit_paths, _q, beam_energy_at, beam_ion, channels, exit_energy, stack, stopping,
                    tilt_deg)
from .response import FWHM_PER_SIGMA, Response, compton_edge

#: Lines of the room background and their typical relative strengths (counts, not intensities): ⁴⁰K, the
#: thorium series (²²⁸Ac, ²¹²Pb, ²⁰⁸Tl), the uranium series (²¹⁴Pb, ²¹⁴Bi) and annihilation. The shares are
#: typical of a germanium detector in a hall without shielding; ``room_background`` scales them together.
ROOM_LINES = (
    (1460.820, "40K", 1.00), (2614.511, "208Tl", 0.35), (583.187, "208Tl", 0.25), (911.204, "228Ac", 0.25),
    (968.971, "228Ac", 0.15), (238.632, "212Pb", 0.40), (609.312, "214Bi", 0.45), (1120.287, "214Bi", 0.15),
    (1764.494, "214Bi", 0.15), (351.932, "214Pb", 0.35), (295.224, "214Pb", 0.18), (511.0, "e+e-", 0.30),
)

#: The columns of :class:`GammaEvents`, with units.
GAMMA_COLUMNS = {
    "event": "event number, as in the particle events",
    "particle": "row of the particle events used for the Doppler correction",
    "cascade": "which decay of the excited nucleus the γ ray belongs to (γ rays of one cascade share it)",
    "initial": "level the γ ray leaves (index in the level scheme; 1 without a scheme)",
    "final": "level the γ ray goes to",
    "nucleus": "which nucleus emitted it: 0 the nucleus of the [reaction], 1 the other one of the collision",
    "crystal": "index into GammaEvents.crystals",
    "energy0": "γ-ray energy in the nucleus's rest frame, MeV",
    "energy_lab": "γ-ray energy in the laboratory (Doppler shifted), MeV",
    "deposited": "energy left in the crystal, MeV",
    "measured": "measured energy (deposited, with resolution), MeV",
    "counted": "above the crystal's threshold",
    "theta": "lab polar angle of the γ ray, degrees",
    "phi": "lab azimuth of the γ ray, degrees",
    "beta": "speed of the emitting nucleus at emission, units of c",
    "corrected_projectile": "measured energy corrected as if the scattered beam had emitted it, MeV",
    "corrected_recoil": "measured energy corrected as if the target recoil had emitted it, MeV",
    "weight": "rate this coincidence stands for, per second",
}


def _unit(v: np.ndarray) -> np.ndarray:
    """Unit vectors; a zero vector is taken along the beam."""
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    safe = np.where(norm > 0, norm, 1.0)
    return np.where(norm > 0, v / safe, np.array([0.0, 0.0, 1.0]))


def _directions(theta_deg, phi_deg) -> np.ndarray:
    th, ph = np.radians(theta_deg), np.radians(phi_deg)
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=-1)


@dataclass
class Crystal:
    name: str
    detector: int
    element: Optional[int]
    centre: np.ndarray
    radius: float
    response: Response
    threshold: float


@dataclass
class GammaEvents:
    """γ rays that reached a crystal in coincidence with a detected particle (see :data:`GAMMA_COLUMNS`), with
    the singles spectra and rates that set the random coincidences."""

    columns: dict
    events: Events
    crystals: list
    seed: int
    #: The transition's energy, MeV, and the emitting nucleus.
    energy_mev: float
    emitter: str
    #: Coincidence window (full width), s; dead time, s; live fraction of the run.
    window_s: float
    dead_time_s: float
    live_fraction: float
    #: Per crystal: the singles rate (1/s) of γ rays of the reaction, of the room background and of the extra
    #: lines, and the expected singles spectrum (rate per bin) on ``singles_edges``.
    singles_rate: dict = field(default_factory=dict)
    singles_spectrum: dict = field(default_factory=dict)
    singles_edges: np.ndarray = field(default_factory=lambda: np.zeros(0))
    #: Particle singles rates, 1/s, by detector name (counted particles).
    particle_rate: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    #: The γ-ray amplifiers' shaping time (s), for pile-up; 0 for none.
    shaping_time_s: float = 0.0

    def __len__(self) -> int:
        return len(self.columns["event"])

    def __getitem__(self, name: str) -> np.ndarray:
        return self.columns[name]

    @property
    def crystal_names(self) -> list:
        return [c.name for c in self.crystals]

    def detector_names(self) -> list:
        """The γ-ray detectors (a clover once, not per crystal)."""
        return list(dict.fromkeys(c.name.rsplit(" ", 1)[0] if c.element is not None else c.name
                                  for c in self.crystals))

    def _crystal_mask(self, name: str) -> np.ndarray:
        idx = [i for i, c in enumerate(self.crystals)
               if c.name == name or (c.element is not None and c.name.rsplit(" ", 1)[0] == name)]
        if not idx:
            raise KeyError(f"no γ-ray detector or crystal named {name!r}; the crystals are "
                           f"{', '.join(self.crystal_names)}")
        return np.isin(self.columns["crystal"], idx)

    def select(self, crystal: Optional[str] = None, detector: Optional[str] = None,
               counted: bool = True) -> np.ndarray:
        """Boolean mask over the γ rays: of a crystal (or a whole γ-ray detector), in coincidence with a
        particle detector, above threshold."""
        m = np.ones(len(self), dtype=bool)
        if crystal is not None:
            m &= self._crystal_mask(crystal)
        if detector is not None:
            rows = self.columns["particle"]
            m &= self.events.columns["detector"][rows] == self.events.detectors.index(detector)
        if counted:
            m &= self.columns["counted"]
        return m

    def rate(self, crystal: Optional[str] = None, detector: Optional[str] = None, counted: bool = True) -> tuple:
        """(true particle–γ coincidences per second, statistical error), before dead time."""
        w = self.columns["weight"][self.select(crystal, detector, counted)]
        return float(w.sum()), float(math.sqrt(np.sum(w**2)))

    def counts(self, crystal: Optional[str] = None, detector: Optional[str] = None) -> float:
        """True coincidences in the planned beam time, after dead time."""
        return self.rate(crystal, detector)[0] * self.events.beam_time_s * self.live_fraction

    def random_rate(self, crystal: str, detector: str) -> float:
        """Random particle–γ coincidences per second: 2τ × (particle singles) × (γ singles), τ being the
        window's full width."""
        singles = sum(sum(self.singles_rate[c.name].values()) for c in self.crystals if self._in(c, crystal))
        return 2 * self.window_s * self.particle_rate.get(detector, 0.0) * singles

    def gamma_random_rate(self, first: str, second: str) -> float:
        """Random γ–γ coincidences per second between two crystals or detectors."""
        a = sum(sum(self.singles_rate[c.name].values()) for c in self.crystals if self._in(c, first))
        b = sum(sum(self.singles_rate[c.name].values()) for c in self.crystals if self._in(c, second))
        return 2 * self.window_s * a * b

    def _in(self, c: Crystal, name: str) -> bool:
        return c.name == name or (c.element is not None and c.name.rsplit(" ", 1)[0] == name)

    def _cascades(self) -> np.ndarray:
        """Which cascade each γ ray belongs to (the particle row, for γ rays made before cascades were kept)."""
        return self.columns["cascade"] if "cascade" in self.columns else self.columns["particle"]

    def gamma_gamma(self, first: Optional[tuple] = None, second: Optional[tuple] = None,
                    first_detector: Optional[str] = None, second_detector: Optional[str] = None,
                    corrected: Optional[str] = None) -> tuple:
        """True γ–γ coincidences: cascades with one counted γ ray in ``first_detector`` (any, if None) with an
        energy in the window ``first`` (MeV; any if None) and another in ``second_detector`` in the window
        ``second``. ``corrected`` is None (measured energies), "projectile" or "recoil". Returns (per second, its
        statistical error); each cascade counts once."""
        key = {None: "measured", "projectile": "corrected_projectile", "recoil": "corrected_recoil"}[corrected]
        x = self.columns[key]
        ok = self.columns["counted"]

        def mask(window, name):
            m = ok.copy()
            if window is not None:
                m &= (x >= window[0]) & (x <= window[1])
            if name is not None:
                m &= self._crystal_mask(name)
            return m

        a, b = mask(first, first_detector), mask(second, second_detector)
        cas = self._cascades()
        ia, ib = np.flatnonzero(a), np.flatnonzero(b)
        if not len(ia) or not len(ib):
            return 0.0, 0.0
        # Pairs of different γ rays of one cascade.
        by_cascade: dict = {}
        for i in ib:
            by_cascade.setdefault(int(cas[i]), []).append(i)
        weights = {}
        for i in ia:
            others = [j for j in by_cascade.get(int(cas[i]), ()) if j != i]
            if others:
                weights[int(cas[i])] = float(self.columns["weight"][i])
        w = np.array(list(weights.values()))
        return float(w.sum()), float(math.sqrt(np.sum(w**2)))

    def coincidences(self) -> dict:
        """The particle × γ-detector matrix: {"true": {(detector, γ detector): counts in the run},
        "random": {...}, "gamma_gamma_random": {(γ detector, γ detector): counts}, "gamma_gamma_true": {...}}
        after dead time. True γ–γ coincidences come from cascades (a level scheme with more than one transition);
        one excited state gives none."""
        t = self.events.beam_time_s * self.live_fraction
        true, random = {}, {}
        for d in self.events.detectors:
            for g in self.detector_names():
                true[(d, g)] = self.counts(g, d)
                random[(d, g)] = self.random_rate(g, d) * t
        gg, gg_true = {}, {}
        names = self.detector_names()
        cascades = len(np.unique(self._cascades())) < len(self)
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                gg[(a, b)] = self.gamma_random_rate(a, b) * t
                gg_true[(a, b)] = (self.gamma_gamma(first_detector=a, second_detector=b)[0] * t if cascades
                                   else 0.0)
        return {"true": true, "random": random, "gamma_gamma_random": gg, "gamma_gamma_true": gg_true}

    def spectrum(self, crystal: Optional[str] = None, detector: Optional[str] = None,
                 corrected: Optional[str] = None, bins: int = 400, range: Optional[tuple] = None,
                 randoms: bool = True, seed: Optional[int] = None) -> tuple:
        """Counts per bin in the planned beam time, after dead time: the true coincidences, with the random ones
        (drawn from the singles spectra with the Poisson statistics of the run) unless ``randoms`` is False.
        ``corrected`` is None (the measured energy), "projectile" or "recoil"."""
        key = {None: "measured", "projectile": "corrected_projectile", "recoil": "corrected_recoil"}[corrected]
        m = self.select(crystal, detector)
        x = self.columns[key][m]
        t = self.events.beam_time_s * self.live_fraction
        if range is None:
            top = max(float(x.max()) if len(x) else 0.0, 1.2 * self.energy_mev,
                      float(self.singles_edges[-1]) if len(self.singles_edges) else 0.0)
            range = (0.0, top)
        h, edges = np.histogram(x, bins=bins, range=range, weights=self.columns["weight"][m] * t)
        if randoms and detector is not None and crystal is not None:
            rng = np.random.default_rng(self.seed + 7919 if seed is None else seed)
            expected = self.random_spectrum(crystal, detector, edges) * t
            h = h + rng.poisson(expected)
        if crystal is not None:
            h = self.pile_up(h, edges, crystal)
        return h, edges

    def pile_rate(self, crystal: str) -> float:
        """The singles rate (1/s) of a crystal, or the mean of a detector's crystals, that sets its pile-up."""
        rates = [sum(self.singles_rate[c.name].values()) for c in self.crystals if self._in(c, crystal)]
        return float(np.mean(rates)) if rates else 0.0

    def pile_up(self, counts: np.ndarray, edges: np.ndarray, crystal: str, singles: bool = False) -> np.ndarray:
        """A spectrum of a crystal (or detector) with pile-up from ``[run] shaping_time`` and its singles rate: a
        loss from every bin and a shoulder above each peak (:func:`physim.nuclear.response.pile_up`); a γ ray in
        coincidence piles with one of the crystal's singles, a singles spectrum with itself (``singles``).
        Unchanged without a shaping time, or on bins that do not start at 0."""
        if not self.shaping_time_s or abs(edges[0]) > 1e-12:
            return counts
        from .response import pile_up

        partner = None
        if not singles:
            centres = (self.singles_edges[:-1] + self.singles_edges[1:]) / 2
            partner = np.zeros(len(edges) - 1)
            for c in self.crystals:
                if self._in(c, crystal):
                    partner += np.histogram(centres, bins=edges, weights=self.singles_spectrum[c.name])[0]
        piled, _ = pile_up(counts, self.pile_rate(crystal), self.shaping_time_s, partner)
        return piled

    def random_spectrum(self, crystal: str, detector: str, edges: np.ndarray) -> np.ndarray:
        """The expected rate per bin (1/s) of random coincidences of a particle detector with a crystal (or a
        whole γ-ray detector): the crystal's singles spectrum scaled to :meth:`random_rate`."""
        total = np.zeros(len(edges) - 1)
        for c in self.crystals:
            if not self._in(c, crystal):
                continue
            s = self.singles_spectrum[c.name]
            centres = (self.singles_edges[:-1] + self.singles_edges[1:]) / 2
            rebinned, _ = np.histogram(centres, bins=edges, weights=s)
            total += rebinned
        return total * 2 * self.window_s * self.particle_rate.get(detector, 0.0)


def _reconstruct(experiment, tb_ej, layers, particle_dirs, theta_lab, recoil_detected, emitter: str,
                 e_mid: float) -> tuple:
    """What two-body kinematics says about the emitting nucleus from the detected particle's direction:
    (unit direction of the emitter, its speed β after leaving the target), for the assumed ``emitter``
    ("ejectile" or "recoil")."""
    n = len(theta_lab)
    # The CM angle of the scattered beam on this event's orbit, from whichever particle was detected.
    th_cm = np.full(n, np.nan)
    ej = ~recoil_detected
    if ej.any():
        th_cm[ej] = np.asarray(tb_ej.theta_cm(theta_lab[ej], "ejectile")[0])
    if recoil_detected.any():
        th_cm[recoil_detected] = 180.0 - np.asarray(tb_ej.theta_cm(theta_lab[recoil_detected], "recoil")[0])
    th_cm = np.where(np.isnan(th_cm), 90.0, th_cm)
    phi_p = np.degrees(np.arctan2(particle_dirs[:, 1], particle_dirs[:, 0]))
    phi_ej = np.where(recoil_detected, phi_p + 180.0, phi_p)
    if emitter == "ejectile":
        point = tb_ej.at_cm(th_cm, "ejectile")
        phi = phi_ej
    else:
        point = tb_ej.recoil_for(th_cm)
        phi = phi_ej + 180.0
    dirs = _directions(np.asarray(point.theta_lab), phi)
    energy = np.nan_to_num(np.asarray(point.energy, dtype=float))
    ion = tb_ej.beam if emitter == "ejectile" else tb_ej.target
    depth = np.full(n, layers[0].thickness / 2)
    e_out = exit_energy(experiment, layers, 0, depth, ion.name, energy, dirs)
    mass = data.nuclide(ion.name).nuclear_mass_mev + (tb_ej.excitation_mev if (emitter == "recoil")
                                                      == (tb_ej.excite == "recoil") else 0.0)
    gam = 1 + np.maximum(e_out, 0.0) / mass
    beta = np.sqrt(np.maximum(1 - 1 / gam**2, 0.0))
    return dirs, beta


def _room_spectrum(experiment, crystal: Crystal, edges: np.ndarray, rate: float, extra: list) -> tuple:
    """(rate per bin of the room background and extra lines in a crystal, their total rate)."""
    s = np.zeros(len(edges) - 1)
    total = 0.0
    if rate > 0:
        share = sum(x[2] for x in ROOM_LINES)
        for e_kev, _, strength in ROOM_LINES:
            r = rate * strength / share
            s += r * crystal.response.shape(e_kev * 1e-3, edges) / crystal.response.peak_to_total(e_kev * 1e-3)
            total += r / float(crystal.response.peak_to_total(e_kev * 1e-3))
    for e, r in extra:
        e_mev, per_s = _q(e).to("MeV"), _q(r).to("/s")
        s += per_s * crystal.response.shape(e_mev, edges) / crystal.response.peak_to_total(e_mev)
        total += per_s / float(crystal.response.peak_to_total(e_mev))
    return s, total


@dataclass
class _Emitter:
    """One nucleus of the collision as a source of γ rays."""

    #: 0 for the nucleus the ``[reaction]`` names, 1 for the other one.
    nucleus: int
    name: str
    #: The target recoil emits (else the scattered beam).
    recoil: bool
    #: The excitation energy in the two-body kinematics and in the emitter's mass, MeV.
    e_star: float
    #: The :class:`~physim.nuclear.orientation.Excitation` of the one ``[reaction]`` state, when no scheme is followed.
    single: object
    #: The solution of the level scheme the γ rays follow, or None.
    cascade: object
    tb: TwoBody
    #: Seeds of the emission, the lifetimes and the add-back.
    seeds: tuple


def _first_crystal(origin: np.ndarray, lab: np.ndarray, centres: np.ndarray, normals: np.ndarray,
                   radii: np.ndarray, block: int = 200_000) -> tuple:
    """The crystal each γ ray meets first (its index) and whether it meets one at all, in blocks of rows so the
    (rows × crystals × 3) intersection points never fill the memory."""
    n = len(lab)
    which = np.zeros(n, dtype=np.int64)
    reached = np.zeros(n, dtype=bool)
    dist = np.linalg.norm(centres, axis=1)
    for a in range(0, n, block):
        o, d = origin[a:a + block], lab[a:a + block]
        along = d @ normals.T                                            # (rows, crystals)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(along > 1e-9, (dist[None, :] - o @ normals.T) / along, np.inf)
        point = o[:, None, :] + d[:, None, :] * t[:, :, None]
        inside = np.linalg.norm(point - centres[None, :, :], axis=2) <= radii[None, :]
        t = np.where(inside, t, np.inf)
        w = np.argmin(t, axis=1)
        which[a:a + block] = w
        reached[a:a + block] = np.isfinite(t[np.arange(len(w)), w])
    return which, reached


#: The shares of :func:`_own_shares`, kept by the crystals' places.
_OWN_SHARES: dict = {}


def _own_shares(crystals: list, rings: int = 60, spokes: int = 96) -> np.ndarray:
    """Of each crystal's face, the share that is its own: the part no other crystal of the same detector covers
    first, as the target sees it (1 for a crystal that stands alone).

    The faces are discs, each turned to the target, and those of a clover overlap: crystals 50 mm across stand
    45 mm apart, as the real ones do, which are cut flat where they meet. A γ ray through the overlap goes to the
    crystal it meets first (:func:`_first_crystal`), so a crystal receives the γ rays of its own part only. Its
    efficiency is the data sheet's, measured on the crystal as cut, and the response counts it over the whole
    disc; the event chain therefore lets the γ rays of the own part interact 1 / share as often, and a run counts
    what the response says. Only crystals of one detector share faces: a detector behind another is hidden by it."""
    key = tuple((cr.detector, tuple(float(x) for x in cr.centre), float(cr.radius)) for cr in crystals)
    if key not in _OWN_SHARES:
        shares = np.ones(len(crystals))
        phi = (np.arange(spokes) + 0.5) / spokes * 2 * math.pi
        for det in sorted({cr.detector for cr in crystals}):
            own = [k for k, cr in enumerate(crystals) if cr.detector == det]
            if len(own) < 2:
                continue
            centres = np.array([crystals[k].centre for k in own], dtype=float)
            radii = np.array([crystals[k].radius for k in own], dtype=float)
            normals = _unit(centres)
            for j, k in enumerate(own):
                n = normals[j]
                a = _unit(np.cross(n, [0.0, 0.0, 1.0] if abs(n[2]) < 0.9 else [1.0, 0.0, 0.0]))
                b = np.cross(n, a)
                rho = radii[j] * np.sqrt((np.arange(rings) + 0.5) / rings)       # cells of equal area
                points = (centres[j] + rho[:, None, None] * (np.cos(phi)[None, :, None] * a
                                                             + np.sin(phi)[None, :, None] * b)).reshape(-1, 3)
                r2 = np.einsum("ij,ij->i", points, points)
                line = points / np.sqrt(r2)[:, None]
                solid = (line @ n) / r2                                          # the solid angle of each cell
                # The cell's own crystal is met at its face; another one counts only if it is met sooner.
                which, _ = _first_crystal(np.zeros_like(line), line, centres, normals, radii * (1 - 1e-9))
                shares[k] = float(solid[which == j].sum() / solid.sum())
        _OWN_SHARES[key] = shares
    return _OWN_SHARES[key]


def _distribution(a_kq: dict, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """4π W / yield of a transition along (θ, φ) of the orbit's frame (radians), from its coefficients a_kq at
    each event (arrays): Σ Re[a_kq √(4π/(2k+1)) Y_kq(θ, φ)]. The Legendre functions and the phases are worked out
    once per (k, |q|) and per |q|, not once per term as :func:`physim.nuclear.angular.spherical_harmonic` would."""
    x, sn = np.cos(theta), np.sin(theta)
    phases: dict = {}
    legendre: dict = {}
    total = np.zeros(len(theta))
    for (k, q), v in a_kq.items():
        m = abs(q)
        if m not in phases:
            phases[m] = np.exp(1j * m * phi)
        if (k, m) not in legendre:
            pmm = (-1) ** m * angular._double_factorial(2 * m - 1) * sn**m
            if k == m:
                plm = pmm
            else:
                a, b = pmm, x * (2 * m + 1) * pmm
                for n in range(m + 2, k + 1):
                    a, b = b, ((2 * n - 1) * x * b - (n + m - 1) * a) / (n - m)
                plm = b
            legendre[(k, m)] = math.sqrt((2 * k + 1) / (4 * math.pi) * math.factorial(k - m)
                                         / math.factorial(k + m)) * plm
        y = legendre[(k, m)] * phases[m]
        if q < 0:
            y = (-1) ** m * np.conj(y)
        total += (v * math.sqrt(4 * math.pi / (2 * k + 1)) * y).real
    return total


def _aimed(rng, tables: dict, initial: np.ndarray, final: np.ndarray, theta_cm: np.ndarray, phi_ej: np.ndarray,
           velocity: np.ndarray, beta: np.ndarray, gam: np.ndarray, origin: np.ndarray, flip: bool,
           centres: np.ndarray, normals: np.ndarray, radii: np.ndarray) -> Optional[tuple]:
    """Which γ rays reach a crystal, and where: (indices of those γ rays, the crystal of each, its unit direction
    in the laboratory). The same γ rays, in distribution, as drawing every direction over the whole sphere and
    keeping those that meet a crystal, at a small part of the work:

    - a γ ray is first kept with a probability q that no γ ray's chance of reaching a crystal exceeds: the
      largest value of its distribution, times the forward throw of the fastest emitter, times the solid angle
      of all the crystal faces from the nearest point;
    - a kept one is aimed at a point drawn evenly over the crystal faces, and counts if that crystal is the
      first its line meets, with probability (its distribution in the laboratory along that line) × (the solid
      angle the faces stand for there) / q.

    The chance that a γ ray counts, and the direction it then has, are those of its distribution; the others
    miss every crystal and were never needed. None when q would reach 1 (crystals all round the target): the
    caller then draws every direction."""
    n = len(theta_cm)
    areas = math.pi * radii**2
    total_area = float(areas.sum())
    dist = np.linalg.norm(centres, axis=1)
    off = float(np.linalg.norm(origin, axis=1).max()) if n else 0.0
    nearest = float(dist.min()) - off
    if n == 0 or nearest <= 0:
        return None
    # No distribution exceeds the sum of the sizes of its coefficients, anywhere on the table of angles.
    largest = max(float(np.max(sum(np.abs(v) for v in coeff.values()))) for _, coeff, _ in tables.values())
    b_max = float(beta.max())
    q = largest * (1 + b_max) / (1 - b_max) * total_area / (4 * math.pi * nearest**2)
    if not q < 1.0:
        return None
    idx = np.flatnonzero(rng.uniform(0.0, 1.0, n) < q)
    m = len(idx)
    if not m:
        return idx, np.zeros(0, dtype=np.int64), np.zeros((0, 3))
    # A point drawn evenly over the faces: the crystal by its share of the area, then the point on its disc.
    k = np.minimum(np.searchsorted(np.cumsum(areas) / total_area, rng.uniform(0.0, 1.0, m)), len(areas) - 1)
    e1 = _unit(np.cross(normals, np.where(np.abs(normals[:, 1:2]) < 0.9, [0.0, 1.0, 0.0], [1.0, 0.0, 0.0])))
    e2 = np.cross(normals, e1)
    rho = radii[k] * np.sqrt(rng.uniform(0.0, 1.0, m))
    psi = rng.uniform(0.0, 2 * math.pi, m)
    point = centres[k] + rho[:, None] * (np.cos(psi)[:, None] * e1[k] + np.sin(psi)[:, None] * e2[k])
    line = point - origin[idx]
    r2 = np.einsum("ij,ij->i", line, line)
    lab = line / np.sqrt(r2)[:, None]
    cos_face = np.einsum("ij,ij->i", lab, normals[k])
    which, reached = _first_crystal(origin[idx], lab, centres, normals, radii)
    ok = reached & (which == k) & (cos_face > 0)   # behind another crystal: that one's own points cover the line
    # The direction in the emitter's rest frame, and how much the motion gathers the γ rays along it.
    be, ga = beta[idx], gam[idx]
    nvel = _unit(velocity[idx])
    cos_a = np.einsum("ij,ij->i", lab, nvel)
    par = (cos_a - be) / (1 - be * cos_a)
    perp = lab - cos_a[:, None] * nvel
    rest = _unit(perp / (ga * (1 - be * cos_a))[:, None] + par[:, None] * nvel)
    throw = (1 - be**2) / (1 - be * cos_a) ** 2
    # The orbit's frame in the laboratory for each event (as Populated.axes, vectorised).
    th, ph = np.radians(theta_cm[idx]), np.radians(phi_ej[idx])
    out = _directions(theta_cm[idx], phi_ej[idx])
    beam = np.array([0.0, 0.0, 1.0])
    s_half, c_half = np.sin(th / 2), np.cos(th / 2)
    side = np.c_[np.cos(ph), np.sin(ph), np.zeros(m)]
    x = np.where(s_half[:, None] > 1e-9, (out - beam) / np.maximum(2 * s_half, 1e-12)[:, None], side)
    y = np.where(c_half[:, None] > 1e-9, (out + beam) / np.maximum(2 * c_half, 1e-12)[:, None], side)
    z = np.cross(x, y)
    if flip:
        x, y = -x, -y
    theta = np.arccos(np.clip(np.einsum("ij,ij->i", rest, z), -1.0, 1.0))
    phi = np.arctan2(np.einsum("ij,ij->i", rest, y), np.einsum("ij,ij->i", rest, x))
    w = np.zeros(m)
    for (i_lev, f_lev), (grid, coeff, _) in tables.items():
        mine = ok & (initial[idx] == i_lev) & (final[idx] == f_lev)
        if not mine.any():
            continue
        a_kq = {key: np.interp(theta_cm[idx][mine], grid, v.real) + 1j * np.interp(theta_cm[idx][mine], grid, v.imag)
                for key, v in coeff.items()}
        w[mine] = _distribution(a_kq, theta[mine], phi[mine])
    chance = np.where(ok, w / (4 * math.pi) * throw * total_area * cos_face / r2, 0.0)
    keep = rng.uniform(0.0, 1.0, m) * q < chance
    return idx[keep], k[keep], lab[keep]


def _chain(experiment, ev: Events, layers: list, crystals: list, bare: list, em: _Emitter, rows: np.ndarray,
           base: np.ndarray, k_rep: int, plain: bool, first_cascade: int = 0, aimed: bool = True) -> dict:
    """The γ rays of one nucleus: from the particle rows ``rows`` of ``ev`` on whose orbits it was excited (each
    standing for the rate ``base``) to what the crystals record. Returns the γ-ray columns before the Doppler
    correction, and "info" for the notes."""
    c = ev.columns
    rng = np.random.default_rng(em.seeds[0])
    excite_recoil = em.recoil
    emitter, e0, tb, cex = em.name, em.e_star, em.tb, em.cascade
    m_emitter = data.nuclide(emitter).nuclear_mass_mev + e0
    centres = np.array([cr.centre for cr in crystals])
    normals = _unit(centres)
    radii = np.array([cr.radius for cr in crystals])

    # -- one row per excited event: the detected particle and the emitter's path -------------------------------
    # One γ-ray chain per detected particle: an event whose two particles both reached a detector is a
    # coincidence with each of them.
    n = len(rows)
    theta_cm = c["theta_cm"][rows]
    recoil_det = c["recoil"][rows]
    phi_det = c["phi"][rows]
    phi_ej = np.where(recoil_det, phi_det + 180.0, phi_det)
    # The emitter's own path: the detected particle if it is the emitter, else its two-body partner.
    is_emitter = recoil_det == excite_recoil
    th_em = np.where(is_emitter, c["theta"][rows], np.nan)
    e_em = np.where(is_emitter, c["energy"][rows], np.nan)
    other = ~is_emitter
    if other.any():
        point = tb.recoil_for(theta_cm[other]) if excite_recoil else tb.at_cm(theta_cm[other], "ejectile")
        th_em[other] = np.asarray(point.theta_lab)
        e_em[other] = np.nan_to_num(np.asarray(point.energy, dtype=float))
    ph_em = np.where(is_emitter, phi_det, phi_det + 180.0)
    dir_em = _directions(th_em, ph_em)
    e_out = exit_energy(experiment, layers, 0, c["depth"][rows], emitter, e_em, dir_em)
    # Kept for the lifetimes (backlog item 73): where and how fast each emitter starts.
    start = {"energy": e_em, "direction": dir_em, "depth": c["depth"][rows]}
    gam = 1 + np.maximum(e_out, 0.0) / m_emitter
    beta = np.sqrt(np.maximum(1 - 1 / gam**2, 0.0))
    velocity = dir_em * beta[:, None]
    # Every excited event emits several times over, each γ ray standing for a share of the event's rate.
    rows, theta_cm, phi_ej = np.repeat(rows, k_rep), np.repeat(theta_cm, k_rep), np.repeat(phi_ej, k_rep)
    base = np.repeat(base, k_rep)
    velocity, beta, gam = np.repeat(velocity, k_rep, axis=0), np.repeat(beta, k_rep), np.repeat(gam, k_rep)
    start = {k: np.repeat(v, k_rep, axis=0) for k, v in start.items()}
    n = len(rows)

    excited_rows = n // k_rep
    # -- the cascade: which levels each excitation passes through, and its γ rays (backlog item 71) -----------------
    cascade = np.arange(n)
    initial = np.ones(n, dtype=int)
    final = np.zeros(n, dtype=int)
    e0s = np.full(n, e0)
    if cex is not None:
        src, initial, final, e0s = _cascade_steps(cex, theta_cm, rng)
        cascade = src
        rows, theta_cm, phi_ej, base = rows[src], theta_cm[src], phi_ej[src], base[src]
        velocity, beta, gam = velocity[src], beta[src], gam[src]
        start = {k: v[src] for k, v in start.items()}
        n = len(rows)

    # -- lifetimes: the emitter slows down in the target and its backing before it decays (backlog item 73) ------
    taus = _mean_lives(experiment, cex, initial)
    stopped = np.zeros(n, dtype=bool)
    if taus is not None and np.any(taus > 0):
        lrng = np.random.default_rng(em.seeds[1])
        draws = lrng.exponential(1.0, n) * taus
        # A level fed from above starts its own clock when the level above decays: the times add up the cascade.
        t_decay = _cumulative_by(cascade, draws)
        e_at, stopped = _slowing(experiment, layers, start["depth"], start["direction"], emitter, start["energy"],
                                 t_decay, m_emitter)
        gam = 1 + np.maximum(e_at, 0.0) / m_emitter
        beta = np.sqrt(np.maximum(1 - 1 / gam**2, 0.0))
        velocity = start["direction"] * beta[:, None]

    # The γ ray leaves from where the reaction happened: the beam spot or halo moves it off the axis (item 74).
    origin = np.zeros((n, 3))
    if "x" in c:
        tilt = math.radians(tilt_deg(experiment))
        origin[:, 0], origin[:, 1] = c["x"][rows], c["y"][rows]
        origin[:, 2] = -origin[:, 0] * math.tan(tilt)

    # -- emission with the correlation, and which crystal ---------------------------------------------------------
    # Each transition with its own orientation (with the feeding from above), at the event's CM angle.
    pairs = [(1, 0)] if cex is None else sorted(set(zip(initial.tolist(), final.tolist())))
    tables = {pair: (em.single if cex is None else cex).coefficient_table(*pair) for pair in pairs}
    hit = _aimed(rng, tables, initial, final, theta_cm, phi_ej, velocity, beta, gam, origin, not excite_recoil,
                 centres, normals, radii) if aimed and n else None
    if hit is not None:
        # Only the γ rays that reach a crystal are drawn, each with its direction (see _aimed).
        reached, which, lab = hit
        nvel = _unit(velocity[reached])
        energy_lab = (e0s[reached] * np.sqrt(1 - beta[reached] ** 2)
                      / (1 - beta[reached] * np.einsum("ij,ij->i", lab, nvel)))
    else:
        # Every direction over the whole sphere, in the rest frame, then the crystal it meets, if any.
        local = np.zeros((n, 3))
        for pair in pairs:
            m = (initial == pair[0]) & (final == pair[1]) if cex is not None else np.ones(n, dtype=bool)
            local[m] = _rest_directions(tables[pair], theta_cm[m], rng)
        # The orbit's frame in the laboratory for each event (as Populated.axes, vectorised).
        th = np.radians(theta_cm)
        ph = np.radians(phi_ej)
        out = _directions(theta_cm, phi_ej)
        beam = np.array([0.0, 0.0, 1.0])
        s_half, c_half = np.sin(th / 2), np.cos(th / 2)
        side = np.c_[np.cos(ph), np.sin(ph), np.zeros(n)]
        x = np.where(s_half[:, None] > 1e-9, (out - beam) / np.maximum(2 * s_half, 1e-12)[:, None], side)
        y = np.where(c_half[:, None] > 1e-9, (out + beam) / np.maximum(2 * c_half, 1e-12)[:, None], side)
        z = np.cross(x, y)
        if not excite_recoil:
            x, y = -x, -y
        rest = local[:, 0:1] * x + local[:, 1:2] * y + local[:, 2:3] * z
        # To the laboratory: the photon's momentum boosted along the emitter's velocity.
        nvel = _unit(velocity)
        cos_r = np.einsum("ij,ij->i", rest, nvel)
        par = (cos_r + beta) / (1 + beta * cos_r)
        perp = rest - cos_r[:, None] * nvel
        scale = 1 / (gam * (1 + beta * cos_r))
        lab = _unit(perp * scale[:, None] + par[:, None] * nvel)
        energy_lab = e0s * np.sqrt(1 - beta**2) / (1 - beta * np.einsum("ij,ij->i", lab, nvel))
        which, reached = _first_crystal(origin, lab, centres, normals, radii)
        which, lab, energy_lab = which[reached], lab[reached], energy_lab[reached]
    rows, beta = rows[reached], beta[reached]
    cascade, initial, final, e0s = cascade[reached], initial[reached], final[reached], e0s[reached]
    stopped = stopped[reached]
    weight = base[reached] / k_rep
    deposited = np.zeros(len(rows))
    measured = np.zeros(len(rows))
    counted = np.zeros(len(rows), dtype=bool)
    keep = np.zeros(len(rows), dtype=bool)
    kind = np.zeros(len(rows), dtype=int)
    own = _own_shares(crystals)
    for k, cr in enumerate(crystals):
        m = which == k
        if not m.any():
            continue
        e = energy_lab[m]
        r = bare[k]
        element = cr.element
        # The crystal's efficiency over the part of its face that is its own (see _own_shares).
        cover = (sum(r.coverage) if element is None else r.coverage[element]) * own[k]
        p_int = r.transmission(e) * np.minimum(r.total_efficiency(e, element) / cover / r.transmission(e), 1.0)
        interact = rng.uniform(0.0, 1.0, m.sum()) <= p_int
        dep, kd = _deposit(r, e, rng)
        meas = dep + rng.normal(0.0, 1.0, len(dep)) * r.fwhm(np.maximum(dep, 1e-3)) / FWHM_PER_SIGMA
        idx = np.flatnonzero(m)
        deposited[idx], measured[idx], kind[idx] = dep, meas, kd
        counted[idx] = interact & (meas >= cr.threshold)
        keep[idx] = interact
    # γ rays of one cascade in one crystal add up: one signal at the sum of their deposits (backlog item 72).
    summed = _sum_cascades(cascade, which, deposited, measured, counted, keep, kind, energy_lab, crystals) \
        if cex is not None else 0
    # Add-back and suppression draw their own numbers, so the chain above is the same with and without them.
    done = {} if plain else _addback_and_suppression(crystals, which, energy_lab, deposited, measured, counted, keep,
                                                     kind, np.random.default_rng(em.seeds[2]))
    rows, which, lab, energy_lab, beta = rows[keep], which[keep], lab[keep], energy_lab[keep], beta[keep]
    weight, deposited, measured, counted = weight[keep], deposited[keep], measured[keep], counted[keep]
    cascade, initial, final, e0s = cascade[keep], initial[keep], final[keep], e0s[keep]
    stopped = stopped[keep]
    return {"particle": rows, "cascade": cascade + first_cascade, "initial": initial, "final": final,
            "nucleus": np.full(len(rows), em.nucleus, dtype=np.int64), "crystal": which, "energy0": e0s,
            "energy_lab": energy_lab, "deposited": deposited, "measured": measured, "counted": counted,
            "theta": np.degrees(np.arccos(np.clip(lab[:, 2], -1, 1))),
            "phi": np.degrees(np.arctan2(lab[:, 1], lab[:, 0])), "beta": beta, "stopped": stopped, "weight": weight,
            "info": {"excited_rows": excited_rows, "summed": summed, "done": done,
                     "lifetimes": taus is not None and bool(np.any(taus > 0)), "chains": excited_rows * k_rep}}


def other_rows(experiment, ev: Events, solution, count_once: bool, seed: int) -> tuple:
    """The particle rows on whose orbits the other nucleus of the collision was excited, and the rate each stands
    for: the particles scattered on the target's main nuclide, whether or not the collision also excited the
    ``[reaction]`` nucleus, each with the probability that it excited the other nucleus at its CM angle (the
    summed populations of ``solution``, :func:`physim.nuclear.gamma.other_excitation`). The two excitations are
    taken as independent: a collision that excites both has the γ rays of both on one particle. The probability
    goes into the weight; or, when every particle stands for one count (``count_once``, a run), the nucleus is
    excited in a row with that probability and the row keeps its weight."""
    from .rates import SchemeProbability

    layers = stack(experiment)
    key = max(layers[0].nuclides.items(), key=lambda kv: kv[1])[0]
    label = f"{data.nuclide(key).name} ({layers[0].name})"
    c = ev.columns
    # Elastic scattering on that nuclide, and the collisions on it that excited the [reaction] nucleus.
    on_it = [i for i, name in enumerate(ev.channels) if name == label or name.startswith(label + ",")]
    if not on_it:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    rows = np.flatnonzero(np.isin(c["channel"], on_it))
    kept = solution.__dict__.get("_scheme_probability")
    if kept is None:
        kept = solution.__dict__["_scheme_probability"] = SchemeProbability(solution)
    p = np.clip(np.asarray(kept.probability(c["theta_cm"][rows]), dtype=float), 0.0, 1.0)
    if count_once:
        pick = np.random.default_rng(seed + 5_000_011).uniform(0.0, 1.0, len(rows)) < p
        rows = rows[pick]
        return rows, np.asarray(c["weight"][rows], dtype=float)
    return rows, np.asarray(c["weight"][rows], dtype=float) * p


def simulate_gammas(experiment, events: int = 200_000, seed: int = 1, particle_events: Optional[Events] = None,
                    bin_kev: float = 1.0, gammas_per_event: int = 10, plain: bool = False,
                    rates: Optional[Rates] = None, count_once: bool = False, aimed: bool = True,
                    other: bool = True) -> GammaEvents:
    """Generate the particle events (or take ``particle_events``) and follow the γ ray of every excited event to
    the crystals. The same setup and seed give the same γ rays.

    Coincidences are rare, so each excited event emits its γ ray ``gammas_per_event`` times over, each with a
    share of the event's weight: the rates are unchanged and the spectra are smoother.

    **Both nuclei.** The ``[reaction]`` names one nucleus; when the setup's ``[levels]`` also hold a level scheme
    of the other one (the target's, for a beam that is excited), that nucleus is excited too, on the orbits of
    the elastically scattered particles, with the probability its coupled equations give
    (:func:`physim.nuclear.gamma.other_excitation`), and its γ rays are in the same events: the ``nucleus``
    column is 0 for the ``[reaction]`` nucleus and 1 for the other, and ``initial`` and ``final`` count levels
    in that nucleus's scheme. The γ rays of the ``[reaction]`` nucleus are the same with and without it.
    ``count_once`` says that every particle stands for one count (a run): the other nucleus is then excited
    in a row with its probability, instead of carrying the probability as weight. ``other=False`` leaves the
    other nucleus out (the extension of a run that was taken without it).

    ``plain`` switches off the add-back and the shields of every γ-ray detector (the geometry stays): the same
    γ rays as the setup's own run up to the add-back and the rejections, for a comparison. ``rates`` is the
    setup's :class:`~physim.nuclear.rates.Rates`, if already computed (it sets the singles and the randoms).

    Only the γ rays that reach a crystal are drawn (a thinning that leaves their number and their directions
    as their distribution has them, at a small part of the work); ``aimed=False`` draws every direction over
    the whole sphere and keeps those that meet a crystal, the slow way to the same γ rays in distribution."""
    exc = experiment.excitation
    if exc is None:
        raise ValueError("the setup has no excited state ([reaction] type = \"coulex\")")
    if not experiment.gamma_detectors:
        raise ValueError("the setup has no γ-ray detectors")
    ev = particle_events if particle_events is not None else simulate(experiment, events, seed)
    layers = stack(experiment)
    chans = channels(layers, experiment)
    excited_channels = [i for i, ch in enumerate(chans) if ch.excitation is not None]
    ex = excitation_of(experiment)
    e0 = exc.energy_mev
    emitter = ex.scheme.nuclide
    c = ev.columns
    k_rep = max(int(gammas_per_event), 1)

    # The chain starts from the bare crystals; add-back and suppression are applied to the deposits afterwards.
    crystals, bare = crystals_of(experiment, plain)

    cex = cascade_excitation(experiment)
    main = _Emitter(0, emitter, exc.excite == "target", e0, ex, cex,
                    TwoBody(ex.beam, ex.target, ex.beam_energy, excitation_mev=e0,
                            excite="recoil" if exc.excite == "target" else "ejectile"),
                    (seed + 1_000_003, seed + 4_000_003, seed + 2_000_003))
    rows = np.flatnonzero(np.isin(c["channel"], excited_channels))
    parts = [_chain(experiment, ev, layers, crystals, bare, main, rows, c["weight"][rows], k_rep, plain,
                    aimed=aimed)]

    # -- the other nucleus of the collision, when the setup has its level scheme -----------------------------------
    other = other_excitation(experiment) if other else None
    if other is not None:
        o_rows, o_base = other_rows(experiment, ev, other, count_once, seed)
        levels = other.scheme.levels
        e_first = 1e-3 * min(lv.energy.value - levels[other.ground].energy.value
                             for n, lv in enumerate(levels) if n != other.ground)
        em = _Emitter(1, other.scheme.nuclide, other.excite == "target", e_first, None, other,
                      TwoBody(ex.beam, ex.target, ex.beam_energy, excitation_mev=e_first,
                              excite="recoil" if other.excite == "target" else "ejectile"),
                      (seed + 6_000_011, seed + 7_000_003, seed + 8_000_009))
        parts.append(_chain(experiment, ev, layers, crystals, bare, em, o_rows, o_base, k_rep, plain,
                            first_cascade=parts[0]["info"]["chains"], aimed=aimed))
    info = [p.pop("info") for p in parts]
    cols = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    rows = cols["particle"]

    # -- the Doppler correction from what the detectors know -------------------------------------------------
    corrected = doppler_correct(experiment, c, rows, cols["crystal"], cols["measured"])

    # -- singles, randoms, dead time ---------------------------------------------------------------------------
    bg = backgrounds(experiment, crystals, rates, bin_kev)
    excitation_rate, all_excitations = bg["excitation_rate"], bg["all_excitations"]
    excited_rows, done, summed = info[0]["excited_rows"], info[0]["done"], info[0]["summed"]
    notes = [f"{len(parts[0]['particle'])} γ rays from {excited_rows} excited events with a detected particle "
             f"({k_rep} emissions "
             f"each); {ev.n_events} reactions "
             f"generated. The particle events stand for {excitation_rate:.3g} excitations/s with a detected "
             f"particle out of {all_excitations:.3g}/s in all."]
    if done.get("added_back") or done.get("suppressed"):
        notes.append(f"Add-back returned {done.get('added_back', 0)} γ rays to the full-energy peak; the shields "
                     f"rejected {done.get('suppressed', 0)}.")
    if summed:
        notes.append(f"{summed} pairs of γ rays of one cascade summed in a crystal.")
    if info[0]["lifetimes"]:
        stopped = parts[0]["stopped"]
        notes.append(f"Lifetimes: each emitter slows down in the target and backing and decays after its level's "
                     f"lifetime; {stopped.mean():.1%} of the γ rays kept come from nuclei that had stopped.")
    how = ("in first order (each level in one step from the ground state)" if excitation_order(experiment)
           == "first" else "to all orders (the coupled equations: levels reached in several steps, and the "
           "reorientation)")
    if cex is not None:
        notes.append(f"Each excitation starts in a level of the scheme, populated {how}, and decays by the level "
                     f"scheme's cascade: {len(cex.transitions())} transitions "
                     "followed, each with its own orientation; the angular correlation between successive γ rays of "
                     "one cascade is not included.")
    if other is not None:
        notes.append(f"{len(parts[1]['particle'])} γ rays of {other.scheme.nuclide}, the other nucleus of the "
                     f"collision, from {info[1]['excited_rows']} scattered particles: it is excited on "
                     f"their orbits with the probability its level scheme gives, {how}, and decays by its "
                     f"{len(other.transitions())} transitions. The two nuclei are excited independently: a "
                     "collision can excite both.")
    cols["event"] = c["event"][rows]
    cols["corrected_projectile"], cols["corrected_recoil"] = corrected["ejectile"], corrected["recoil"]
    shaping = _q(experiment.run.shaping_time).to("s") if experiment.run.shaping_time is not None else 0.0
    return GammaEvents(cols, ev, crystals, int(seed), e0, emitter, bg["window_s"], bg["dead_time_s"],
                       bg["live_fraction"], bg["singles_rate"], bg["singles_spectrum"], bg["singles_edges"],
                       bg["particle_rate"], notes, shaping)


def crystals_of(experiment, plain: bool = False) -> tuple:
    """The crystals of the setup's γ-ray detectors with their responses (``plain``: without add-back and
    shields), and the bare response of each crystal: ([Crystal], [Response])."""
    crystals, bare = [], []
    for i, gd in enumerate(experiment.gamma_detectors):
        r = Response(experiment, gd, bare=plain)
        rb = Response(experiment, gd, bare=True)
        thr = _q(gd.threshold).to("MeV") if gd.threshold is not None else 0.0
        for k, (label, centre, radius) in enumerate(r.elements):
            name = (gd.name or f"G{i + 1}") + (f" {label}" if label else "")
            crystals.append(Crystal(name, i, k if len(r.elements) > 1 else None, np.array(centre, dtype=float),
                                    radius, r, thr))
            bare.append(rb)
    return crystals, bare


def spectrum_top(experiment) -> float:
    """The highest energy (MeV) a γ-ray spectrum of the setup has to show: 1.3 × the energy of the state the
    reaction excites, and, when a nucleus of the collision has a level scheme, 1.1 × its highest γ-ray energy. The
    second matters as soon as the scheme holds more than the first state: 62Ni's second 2⁺ state decays to the
    ground state with 2302 keV, well above 1.3 × 1173 keV, and a spectrum that stopped there showed none of it."""
    exc = experiment.excitation
    if exc is None:
        return 0.0
    top = 1.3 * exc.energy_mev
    # Both nuclei: the other one's lines are in the events when the setup has its scheme.
    for scheme in (experiment.levels.get("target"), experiment.levels.get("beam")):
        if scheme is not None and scheme.transitions:
            top = max(top, 1.1e-3 * max(float(t.energy.value) for t in scheme.transitions))
    return float(top)


def backgrounds(experiment, crystals: list, rates: Optional[Rates] = None, bin_kev: float = 1.0) -> dict:
    """What sets the random coincidences and the dead time, from the rates rather than the events: the
    coincidence window and dead time (s), the particle singles rate per detector, the γ-ray singles rate and
    spectrum per crystal (the reaction's γ rays, the room background and the extra lines), the live fraction,
    and the excitation rates (with a detected particle, and in all)."""
    e0 = experiment.excitation.energy_mev
    ex = excitation_of(experiment)
    chans = channels(stack(experiment), experiment)
    array = Array.from_experiment(experiment)
    run = experiment.run
    window = _q(run.coincidence_window).to("s") if run.coincidence_window is not None else 100e-9
    dead = _q(run.dead_time).to("s") if run.dead_time is not None else 0.0
    room = _q(run.room_background).to("/s") if run.room_background is not None else 0.0
    extra = run.extra_lines or []
    rates = rates if rates is not None else Rates(experiment)
    particle_rate = {g.name: rates.rate(g.name) for g in array.geometries}
    excitation_rate = sum(rates.rate(g.name, what="excitations") for g in array.geometries)
    # All excitations, whether or not a particle was detected: the γ-ray singles.
    cx = ex.paths[1][0][2]
    pps = experiment.beam.particles_per_second
    atoms = sum(ch.atoms_per_cm2 for ch in chans if ch.excitation is not None)
    all_excitations = pps * atoms * cx.total() * 1e-27
    top = max(spectrum_top(experiment), max(x[0] for x in ROOM_LINES) * 1e-3 * 1.05 if room > 0 else 0.0,
              max((_q(e).to("MeV") for e, _ in extra), default=0.0) * 1.1) + 0.05
    edges = np.arange(0.0, top, bin_kev * 1e-3)
    singles_rate, singles_spectrum = {}, {}
    other = other_lines(experiment)
    for cr in crystals:
        r = cr.response
        gamma_rate = all_excitations * float(r.total_efficiency(e0, cr.element))
        s = gamma_rate * r.shape(e0, edges)
        bg, bg_rate = _room_spectrum(experiment, cr, edges, room, extra)
        singles_rate[cr.name] = {"reaction": gamma_rate, "background": bg_rate}
        if other:  # the other nucleus's γ rays, line by line (unshifted, as the reaction's own)
            o_rate = 0.0
            for e_line, per_s in other:
                seen = per_s * float(r.total_efficiency(e_line, cr.element))
                s = s + seen * r.shape(e_line, edges)
                o_rate += seen
            singles_rate[cr.name]["other nucleus"] = o_rate
        singles_spectrum[cr.name] = s + bg
    total_rate = sum(particle_rate.values()) + sum(sum(v.values()) for v in singles_rate.values())
    return {"window_s": window, "dead_time_s": dead, "live_fraction": 1 / (1 + dead * total_rate),
            "particle_rate": particle_rate, "singles_rate": singles_rate, "singles_spectrum": singles_spectrum,
            "singles_edges": edges, "excitation_rate": excitation_rate, "all_excitations": all_excitations}


def other_lines(experiment) -> list:
    """The γ rays the other nucleus of the collision emits in all, whether or not a particle is detected:
    [(energy in MeV, γ rays per second)] for each transition of its scheme, from the Rutherford cross section
    times the transition's yield per collision, integrated over the scattering angle. Empty without a level
    scheme for it (:func:`physim.nuclear.gamma.other_excitation`)."""
    from .rutherford import Rutherford

    sol = other_excitation(experiment)
    if sol is None:
        return []
    kept = sol.__dict__.get("_other_lines")
    layers = stack(experiment)
    key, atoms = max(layers[0].nuclides.items(), key=lambda kv: kv[1])
    pps = experiment.beam.particles_per_second
    if kept is None:
        grid, states = sol.table()
        u = np.linspace(math.cos(math.radians(1.0)), -1.0, 400)
        th = np.degrees(np.arccos(np.clip(u, -1, 1)))
        ruth = np.asarray(Rutherford(sol.beam, sol.target, sol.beam_energy).cross_section_cm(th))  # mb/sr
        levels = sol.scheme.levels
        kept = []
        for i, f in sol.transitions():
            y = np.interp(th, np.r_[0.0, grid], np.r_[0.0, [st.gamma_yield(i, f) for st in states]]) * ruth
            mb = float(2 * math.pi * abs(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(u))))
            if mb > 0:
                kept.append((1e-3 * float(levels[i].energy.value - levels[f].energy.value), mb))
        sol.__dict__["_other_lines"] = kept
    return [(e, pps * atoms * mb * 1e-27) for e, mb in kept]


LN2 = math.log(2)
#: The speed of light, cm/s.
C_CM_S = 2.99792458e10


def _mean_lives(experiment, cex, initial: np.ndarray) -> Optional[np.ndarray]:
    """The mean life (s) of the level each γ ray leaves, from the level scheme's half-lives: the scheme of the
    cascade, or for one state the setup's scheme level at the state's energy. None when no level has a half-life
    (every decay then happens after the target, as before)."""
    exc = experiment.excitation
    if cex is not None:
        scheme = cex.scheme
        lives = np.array([(lv.half_life.value / LN2) if lv.half_life is not None else 0.0 for lv in scheme.levels])
        return lives[initial] if np.any(lives > 0) else None
    scheme = experiment.levels.get("target" if exc.excite == "target" else "beam")
    if scheme is None:
        return None
    for lv in scheme.levels:
        if abs(lv.energy.value - 1e3 * exc.energy_mev) < 1.0 and lv.half_life is not None and lv.half_life.value > 0:
            return np.full(len(initial), lv.half_life.value / LN2)
    return None


def _cumulative_by(group: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Running sums of ``x`` within runs of equal ``group`` (``group`` sorted, as the cascade steps are)."""
    if not len(x):
        return x
    total = np.cumsum(x)
    starts = np.r_[0, np.flatnonzero(np.diff(group)) + 1]
    offset = np.repeat(np.r_[0.0, total[starts[1:] - 1]], np.diff(np.r_[starts, len(x)]))
    return total - offset


_CLOCKS: dict = {}


def _clock(species: str, material, density: float, mass_mev: float) -> tuple:
    """The time (s) a nucleus of ``species`` takes to stop from each energy of its stopping table in ``material``
    of density ``density`` (g/cm³): t(E) = ∫ dR / (ρ v), on the table's energies (MeV). Kept once computed."""
    key = (species, material.name, density, mass_mev)
    if key in _CLOCKS:
        return _CLOCKS[key]
    st = stopping(species, material)
    e, r = np.asarray(st._e, dtype=float), np.asarray(st._range, dtype=float)  # MeV, mg/cm²
    beta = np.sqrt(np.maximum(e * (e + 2 * mass_mev), 0.0)) / (e + mass_mev)
    mid = 0.5 * (beta[1:] + beta[:-1])
    dt = np.diff(r) / (1000.0 * density * C_CM_S * np.maximum(mid, 1e-12))
    t = np.r_[r[0] / (1000.0 * density * C_CM_S * max(beta[0], 1e-12)), 0.0]
    t = np.cumsum(np.r_[t[0], dt])
    _CLOCKS[key] = (np.r_[0.0, e], np.r_[0.0, t])
    return _CLOCKS[key]


def _slowing(experiment, layers: list, depth, dirs, species: str, energy, t_decay, mass_mev: float) -> tuple:
    """The emitter's energy (MeV) at its decay time ``t_decay`` (s after the reaction): it slows along its path
    through the rest of the stack, then flies at its exit energy; one that stops first decays at rest. Returns
    (energy at the decay, whether it had stopped). The direction is kept (no angular straggling); time dilation
    (a few parts in a thousand here) is left out."""
    fwd, steps = _exit_paths(experiment, layers, 0, depth, dirs)
    e = np.asarray(energy, dtype=float).copy()
    elapsed = np.zeros(len(e))
    done = np.zeros(len(e), dtype=bool)
    stopped = np.zeros(len(e), dtype=bool)
    out = np.full(len(e), np.nan)
    for j, path, forward in steps:
        m = (fwd == forward) & ~done & (np.asarray(path) > 0)
        if not m.any():
            continue
        mat = layers[j].material
        density = mat.density_g_cm3
        if not density:
            raise ValueError(f"lifetimes need the density of {mat.name}, which is not known")
        grid_e, grid_t = _clock(species, mat, float(density), float(mass_mev))
        idx = np.flatnonzero(m)
        e_in = e[idx]
        e_out = _after(stopping(species, mat), e_in, np.asarray(path)[idx])
        t_in, t_out = np.interp(e_in, grid_e, grid_t), np.interp(e_out, grid_e, grid_t)
        here = t_decay[idx] < elapsed[idx] + (t_in - t_out)
        out[idx[here]] = np.interp(t_in[here] - (t_decay[idx[here]] - elapsed[idx[here]]), grid_t, grid_e)
        done[idx[here]] = True
        rest = (e_out <= 0) & ~here
        out[idx[rest]] = 0.0
        stopped[idx[rest]] = True
        done[idx[rest]] = True
        go = ~here & ~rest
        elapsed[idx[go]] += (t_in - t_out)[go]
        e[idx[go]] = e_out[go]
    out[~done] = e[~done]  # out of the stack: in flight at the exit energy
    return out, stopped


def _sum_cascades(cascade, which, deposited, measured, counted, keep, kind, energy_lab, crystals) -> int:
    """Two γ rays of one cascade that interact in one crystal give one signal: the second's deposit is added to
    the first's, and the second is dropped. In place; returns how many were summed."""
    done = 0
    for _ in range(8):  # a cascade of more than two γ rays in one crystal sums step by step
        idx = np.flatnonzero(keep)
        if len(idx) < 2:
            break
        order = idx[np.lexsort((which[idx], cascade[idx]))]
        same = (cascade[order[1:]] == cascade[order[:-1]]) & (which[order[1:]] == which[order[:-1]])
        if not same.any():
            break
        first, second = order[:-1][same], order[1:][same]
        # One pair per γ ray in this pass.
        _, unique_first = np.unique(first, return_index=True)
        first, second = first[unique_first], second[unique_first]
        clash = np.isin(second, first)
        first, second = first[~clash], second[~clash]
        deposited[first] += deposited[second]
        measured[first] += measured[second]
        energy_lab[first] += energy_lab[second]
        kind[first] = np.where((kind[first] == 0) & (kind[second] == 0), 0, 3)
        thresholds = np.array([c.threshold for c in crystals])
        counted[first] = measured[first] >= thresholds[which[first]]
        keep[second] = False
        counted[second] = False
        done += len(first)
    return done


def _rest_directions(table: tuple, theta_cm: np.ndarray, rng) -> np.ndarray:
    """Unit vectors in the orbit's frame drawn from a transition's distribution W(θ, φ), its coefficients a_kq
    (``table`` from :meth:`~physim.nuclear.orientation.Excitation.coefficient_table`) interpolated at each event's
    CM angle; by acceptance against Σ |a_kq|."""
    grid, coeff, _ = table
    n = len(theta_cm)
    a_kq = {key: np.interp(theta_cm, grid, v.real) + 1j * np.interp(theta_cm, grid, v.imag)
            for key, v in coeff.items()}
    bound = sum(np.abs(v) for v in a_kq.values())  # W × 4π ≤ Σ |a_kq|, since |Y_kq| ≤ √((2k+1)/4π)
    local = np.zeros((n, 3))
    todo = np.ones(n, dtype=bool)
    for _ in range(200):
        idx = np.flatnonzero(todo)
        if not len(idx):
            break
        u = rng.uniform(-1.0, 1.0, len(idx))
        ph = rng.uniform(0.0, 2 * math.pi, len(idx))
        th = np.arccos(u)
        w = _distribution({key: v[idx] for key, v in a_kq.items()}, th, ph)
        accept = rng.uniform(0.0, 1.0, len(idx)) * bound[idx] <= w
        hit = idx[accept]
        local[hit] = np.c_[np.sin(th[accept]) * np.cos(ph[accept]), np.sin(th[accept]) * np.sin(ph[accept]),
                           u[accept]]
        todo[hit] = False
    return local


def _cascade_steps(cex, theta_cm: np.ndarray, rng) -> tuple:
    """The decay of each excitation through the level scheme: the level it starts in, drawn from the populations
    the collision leaves at its CM angle (to all orders, or the one-step probabilities of first order: whichever
    solution ``cex`` is), then each step down by the scheme's branches (a converted transition gives no γ ray).
    Returns (index of the excitation each γ ray belongs to, initial level, final level, energy in the nucleus's
    frame in MeV), one entry per γ ray emitted."""
    grid, states = cex.table()
    starts = start_levels(cex)
    p = np.array([[s.direct.get(lev, 0.0) for s in states] for lev in starts])
    p_at = np.array([np.interp(theta_cm, grid, row) for row in p])  # (levels, events)
    cum = np.cumsum(p_at, axis=0)
    u = rng.uniform(0.0, 1.0, len(theta_cm)) * cum[-1]
    pick = np.minimum((u[None, :] > cum).sum(axis=0), len(starts) - 1)
    current = np.array(starts)[pick]
    levels = cex.scheme.levels
    decays: dict = {}
    for (i, f) in cex.transitions():
        decays.setdefault(i, []).append(cex.decay(i, f))
    src, ini, fin = [], [], []
    alive = np.ones(len(current), dtype=bool)
    for _ in range(4 * len(levels) + 4):
        idx = np.flatnonzero(alive)
        if not len(idx):
            break
        for lev in np.unique(current[idx]):
            here = idx[current[idx] == lev]
            ways = decays.get(int(lev))
            if not ways:
                alive[here] = False
                continue
            share = np.array([d.branch for d in ways])
            cum_b = np.cumsum(share / share.sum())
            k = np.minimum((rng.uniform(0.0, 1.0, len(here))[:, None] > cum_b[None, :]).sum(axis=1), len(ways) - 1)
            gamma_share = np.array([d.gamma / d.branch if d.branch > 0 else 0.0 for d in ways])[k]
            emits = rng.uniform(0.0, 1.0, len(here)) < gamma_share
            finals = np.array([d.final for d in ways])[k]
            src.append(here[emits])
            ini.append(np.full(emits.sum(), int(lev)))
            fin.append(finals[emits])
            current[here] = finals
    if not src:
        return (np.zeros(0, dtype=int), np.zeros(0, dtype=int), np.zeros(0, dtype=int), np.zeros(0))
    src, ini, fin = np.concatenate(src), np.concatenate(ini), np.concatenate(fin)
    energy = np.array([levels[i].energy.value - levels[f].energy.value for i, f in zip(ini, fin)]) * 1e-3
    order = np.argsort(src, kind="stable")
    return src[order], ini[order], fin[order], energy[order]


def doppler_correct(experiment, columns: dict, rows: np.ndarray, crystal: np.ndarray, measured: np.ndarray) -> dict:
    """The Doppler correction an experimentalist makes with the geometry of ``experiment``: from the centre of
    the segment each particle row hit (``columns`` of the particle events, ``rows`` their indices), the centre of
    the crystal (``crystal`` indices) and two-body kinematics at the nominal beam energy. Returns
    {"ejectile": corrected energies, "recoil": ...}.

    Given an experiment whose target sits elsewhere than the one that made the events, this is the correction with
    a misplaced target: the peaks shift and broaden, ring by ring."""
    layers = stack(experiment)
    exc = experiment.excitation
    ex = excitation_of(experiment)
    tb = TwoBody(ex.beam, ex.target, ex.beam_energy, excitation_mev=exc.energy_mev,
                 excite="recoil" if exc.excite == "target" else "ejectile")
    array = Array.from_experiment(experiment)
    centres = []
    for gd in experiment.gamma_detectors:
        centres += [np.array(c, dtype=float) for _, c, _ in gd.elements()]
    normals = _unit(np.array(centres))
    seg_centre = np.zeros((len(rows), 3))
    cache: dict = {}
    for i, g in enumerate(array.geometries):
        m = columns["detector"][rows] == i
        for j in np.flatnonzero(m):
            key = (i, int(columns["segment_i"][rows[j]]), int(columns["segment_j"][rows[j]]))
            if key not in cache:  # one segment's centre serves every hit in it
                cache[key] = g.segment_centre(key[1:], weighted=True)
            seg_centre[j] = cache[key]
    p_dir = _unit(seg_centre)
    theta_lab_p = np.degrees(np.arccos(np.clip(p_dir[:, 2], -1, 1)))
    crystal_dir = normals[crystal]
    out = {}
    for assumed in ("ejectile", "recoil"):
        d_em, b_em = _reconstruct(experiment, tb, layers, p_dir, theta_lab_p, columns["recoil"][rows], assumed,
                                  ex.beam_energy)
        cos_a = np.einsum("ij,ij->i", crystal_dir, d_em)
        out[assumed] = measured * (1 - b_em * cos_a) / np.sqrt(1 - b_em**2)
    return out


def recorrect(gammas: GammaEvents, experiment) -> dict:
    """The γ rays of ``gammas`` corrected with the geometry of ``experiment`` (the assumed one), as
    {"corrected_projectile": ..., "corrected_recoil": ...}."""
    out = doppler_correct(experiment, gammas.events.columns, gammas["particle"], gammas["crystal"], gammas["measured"])
    return {"corrected_projectile": out["ejectile"], "corrected_recoil": out["recoil"]}


def _addback_and_suppression(crystals: list, which: np.ndarray, energy_lab: np.ndarray, deposited: np.ndarray,
                             measured: np.ndarray, counted: np.ndarray, keep: np.ndarray, kind: np.ndarray,
                             rng) -> dict:
    """Add-back and Compton suppression on the deposits, in place. In a clover with add-back, a Compton deposit
    is summed with what the scattered photon left in a neighbouring crystal (one of the two beside it round the
    square), with the share :meth:`Response.addback_share` gives, and the γ ray goes to the crystal with the
    larger deposit; in a shielded detector a deposit that is not the full energy is kept with probability 1/S(E).
    Returns {"added_back": n, "suppressed": n}."""
    done = {"added_back": 0, "suppressed": 0}
    by_detector: dict = {}
    for k, cr in enumerate(crystals):
        by_detector.setdefault(cr.detector, []).append(k)
    for members in by_detector.values():
        r = crystals[members[0]].response
        mine = np.isin(which, members) & keep
        if r.addback and len(members) == 4:
            idx = np.flatnonzero(mine & (kind == 3))
            if len(idx):
                idx = idx[rng.uniform(0.0, 1.0, len(idx)) <= r.addback_share(energy_lab[idx])]
            if len(idx):
                first = deposited[idx]
                second = energy_lab[idx] - first
                noise = np.sqrt(r.fwhm(np.maximum(first, 1e-3)) ** 2
                                + r.fwhm(np.maximum(second, 1e-3)) ** 2) / FWHM_PER_SIGMA
                deposited[idx] = energy_lab[idx]
                measured[idx] = energy_lab[idx] + rng.normal(0.0, 1.0, len(idx)) * noise
                element = np.array([crystals[w].element for w in which[idx]])
                step = np.where(rng.uniform(0.0, 1.0, len(idx)) < 0.5, 1, -1)
                neighbour = np.array(members)[(element + step) % 4]
                switch = second > first
                which[idx[switch]] = neighbour[switch]
                kind[idx] = 0
                counted[idx] = measured[idx] >= np.array([crystals[w].threshold for w in which[idx]])
                done["added_back"] += int(len(idx))
        if r.suppressed:
            idx = np.flatnonzero(mine & (kind != 0))
            if len(idx):
                reject = rng.uniform(0.0, 1.0, len(idx)) > 1.0 / r.suppression_factor(energy_lab[idx])
                keep[idx[reject]] = False
                done["suppressed"] += int(reject.sum())
    return done


def _deposit(r: Response, energy: np.ndarray, rng) -> tuple:
    """The energy each interacting γ ray leaves: the full energy, an escape peak, or a Compton deposit
    (Klein–Nishina for one scattering, or the flat part up to the full energy for several). Returns (deposit,
    kind) with kind 0 for the full energy, 1 and 2 for the single and double escape peaks, 3 for a Compton
    deposit."""
    e = np.asarray(energy, dtype=float)
    pt = r.peak_to_total(e)
    u = rng.uniform(0.0, 1.0, len(e))
    out = e.copy()
    se = np.zeros_like(e)
    de = np.zeros_like(e)
    above = e > 2 * 0.51099895
    if above.any():
        x = e[above] - 2 * 0.51099895
        d = np.minimum(r.crystal.escape_cap, r.crystal.double_escape * x**1.5)
        de[above] = np.minimum(d, (1 - pt[above]) / 2)
        se[above] = np.minimum(r.crystal.single_escape * d, (1 - pt[above]) / 2)
    single = (u >= pt) & (u < pt + se)
    double = (u >= pt + se) & (u < pt + se + de)
    compton = u >= pt + se + de
    out[single] = e[single] - 0.51099895
    out[double] = e[double] - 2 * 0.51099895
    if compton.any():
        ec = e[compton]
        flat = rng.uniform(0.0, 1.0, len(ec)) < r.crystal.multiple
        edge = compton_edge(ec)
        out_c = np.empty(len(ec))
        out_c[flat] = edge[flat] + rng.uniform(0.0, 1.0, flat.sum()) * (ec[flat] - edge[flat])
        k = ec[~flat] / 0.51099895
        # Klein–Nishina electron spectrum sampled by inversion on a grid of the fraction T/E.
        s = np.linspace(0.0, 1.0, 801)[:, None] * (2 * k / (1 + 2 * k))[None, :]  # up to the edge
        kn = 2 + s**2 / (k**2 * (1 - s) ** 2) + s / (1 - s) * (s - 2 / k)
        cdf = np.cumsum(kn, axis=0)
        cdf /= cdf[-1:]
        v = rng.uniform(0.0, 1.0, (~flat).sum())
        # s at v on each γ ray's own cumulative curve: a straight line between the two points that hold v.
        cols = np.arange(len(v))
        j = np.clip((cdf <= v[None, :]).sum(axis=0) - 1, 0, len(cdf) - 2)
        lo, hi = cdf[j, cols], cdf[j + 1, cols]
        frac = np.where(hi > lo, (v - lo) / np.where(hi > lo, hi - lo, 1.0), 0.0)
        pick = np.where(v <= cdf[0], s[0], s[j, cols] + np.clip(frac, 0.0, 1.0) * (s[j + 1, cols] - s[j, cols]))
        out_c[~flat] = pick * ec[~flat]
        out[compton] = out_c
    kind = np.zeros(len(e), dtype=int)
    kind[single], kind[double], kind[compton] = 1, 2, 3
    return out, kind


__all__ = ["Crystal", "GAMMA_COLUMNS", "GammaEvents", "ROOM_LINES", "backgrounds", "crystals_of", "doppler_correct",
           "other_lines", "other_rows", "spectrum_top",
           "recorrect", "simulate_gammas"]
