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
from .gamma import cascade_excitation, excitation_of
from .kinematics import TwoBody
from .quantity import Quantity
from .rates import Rates, _q, beam_energy_at, beam_ion, channels, exit_energy, stack
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


def simulate_gammas(experiment, events: int = 200_000, seed: int = 1, particle_events: Optional[Events] = None,
                    bin_kev: float = 1.0, gammas_per_event: int = 10, plain: bool = False,
                    rates: Optional[Rates] = None) -> GammaEvents:
    """Generate the particle events (or take ``particle_events``) and follow the γ ray of every excited event to
    the crystals. The same setup and seed give the same γ rays.

    Coincidences are rare, so each excited event emits its γ ray ``gammas_per_event`` times over, each with a
    share of the event's weight: the rates are unchanged and the spectra are smoother.

    ``plain`` switches off the add-back and the shields of every γ-ray detector (the geometry stays): the same
    γ rays as the setup's own run up to the add-back and the rejections, for a comparison. ``rates`` is the
    setup's :class:`~physim.nuclear.rates.Rates`, if already computed (it sets the singles and the randoms)."""
    exc = experiment.excitation
    if exc is None:
        raise ValueError("the setup has no excited state ([reaction] type = \"coulex\")")
    if not experiment.gamma_detectors:
        raise ValueError("the setup has no γ-ray detectors")
    ev = particle_events if particle_events is not None else simulate(experiment, events, seed)
    rng = np.random.default_rng(seed + 1_000_003)
    layers = stack(experiment)
    chans = channels(layers, experiment)
    excited_channels = [i for i, ch in enumerate(chans) if ch.excitation is not None]
    ex = excitation_of(experiment)
    e0 = exc.energy_mev
    excite_recoil = exc.excite == "target"
    emitter = ex.scheme.nuclide
    e_mid = ex.beam_energy
    tb = TwoBody(ex.beam, ex.target, e_mid, excitation_mev=e0, excite="recoil" if excite_recoil else "ejectile")
    m_emitter = data.nuclide(emitter).nuclear_mass_mev + e0

    # -- the crystals and their responses ----------------------------------------------------------------------
    # The chain starts from the bare crystals; add-back and suppression are applied to the deposits afterwards.
    crystals, bare = crystals_of(experiment, plain)
    centres = np.array([c.centre for c in crystals])
    normals = _unit(centres)
    radii = np.array([c.radius for c in crystals])

    # -- one row per excited event: the detected particle and the emitter's path -------------------------------
    c = ev.columns
    excited = np.isin(c["channel"], excited_channels)
    # One γ-ray chain per detected particle: an event whose two particles both reached a detector is a
    # coincidence with each of them.
    rows = np.flatnonzero(excited)
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
    gam = 1 + np.maximum(e_out, 0.0) / m_emitter
    beta = np.sqrt(np.maximum(1 - 1 / gam**2, 0.0))
    velocity = dir_em * beta[:, None]
    # Every excited event emits several times over, each γ ray standing for a share of the event's rate.
    k_rep = max(int(gammas_per_event), 1)
    rows, theta_cm, phi_ej = np.repeat(rows, k_rep), np.repeat(theta_cm, k_rep), np.repeat(phi_ej, k_rep)
    velocity, beta, gam = np.repeat(velocity, k_rep, axis=0), np.repeat(beta, k_rep), np.repeat(gam, k_rep)
    n = len(rows)

    excited_rows = n // k_rep
    # -- the cascade: which levels each excitation passes through, and its γ rays (backlog item 71) -----------------
    cex = cascade_excitation(experiment)
    cascade = np.arange(n)
    initial = np.ones(n, dtype=int)
    final = np.zeros(n, dtype=int)
    e0s = np.full(n, e0)
    if cex is not None:
        src, initial, final, e0s = _cascade_steps(cex, theta_cm, rng)
        cascade = src
        rows, theta_cm, phi_ej = rows[src], theta_cm[src], phi_ej[src]
        velocity, beta, gam = velocity[src], beta[src], gam[src]
        n = len(rows)

    # -- emission in the rest frame, with the correlation ----------------------------------------------------
    if cex is None:
        local = _rest_directions(ex.coefficient_table(1, 0), theta_cm, rng)
    else:
        # Each transition with its own orientation (with the feeding from above), at the event's CM angle.
        local = np.zeros((n, 3))
        for i_lev, f_lev in sorted(set(zip(initial.tolist(), final.tolist()))):
            m = (initial == i_lev) & (final == f_lev)
            local[m] = _rest_directions(cex.coefficient_table(i_lev, f_lev), theta_cm[m], rng)
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

    # -- which crystal, and what it records ------------------------------------------------------------------
    along = lab @ normals.T                                          # (n, crystals)
    dist = np.linalg.norm(centres, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(along > 1e-9, dist[None, :] / along, np.inf)
    point = lab[:, None, :] * t[:, :, None]
    inside = np.linalg.norm(point - centres[None, :, :], axis=2) <= radii[None, :]
    t = np.where(inside, t, np.inf)
    which = np.argmin(t, axis=1)
    reached = np.isfinite(t[np.arange(n), which])
    rows, which, lab, energy_lab, beta = rows[reached], which[reached], lab[reached], energy_lab[reached], beta[reached]
    cascade, initial, final, e0s = cascade[reached], initial[reached], final[reached], e0s[reached]
    weight = c["weight"][rows] / k_rep
    deposited = np.zeros(len(rows))
    measured = np.zeros(len(rows))
    counted = np.zeros(len(rows), dtype=bool)
    keep = np.zeros(len(rows), dtype=bool)
    kind = np.zeros(len(rows), dtype=int)
    for k, cr in enumerate(crystals):
        m = which == k
        if not m.any():
            continue
        e = energy_lab[m]
        r = bare[k]
        element = cr.element
        p_int = (r.transmission(e) * np.minimum(r.total_efficiency(e, element)
                                                / (sum(r.coverage) if element is None else r.coverage[element])
                                                / r.transmission(e), 1.0))
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
                                                     kind, np.random.default_rng(seed + 2_000_003))
    rows, which, lab, energy_lab, beta = rows[keep], which[keep], lab[keep], energy_lab[keep], beta[keep]
    weight, deposited, measured, counted = weight[keep], deposited[keep], measured[keep], counted[keep]
    cascade, initial, final, e0s = cascade[keep], initial[keep], final[keep], e0s[keep]

    # -- the Doppler correction from what the detectors know -------------------------------------------------
    corrected = doppler_correct(experiment, c, rows, which, measured)

    # -- singles, randoms, dead time ---------------------------------------------------------------------------
    bg = backgrounds(experiment, crystals, rates, bin_kev)
    excitation_rate, all_excitations = bg["excitation_rate"], bg["all_excitations"]
    notes = [f"{len(rows)} γ rays from {excited_rows} excited events with a detected particle ({k_rep} emissions "
             f"each); {ev.n_events} reactions "
             f"generated. The particle events stand for {excitation_rate:.3g} excitations/s with a detected "
             f"particle out of {all_excitations:.3g}/s in all."]
    if done.get("added_back") or done.get("suppressed"):
        notes.append(f"Add-back returned {done.get('added_back', 0)} γ rays to the full-energy peak; the shields "
                     f"rejected {done.get('suppressed', 0)}.")
    if summed:
        notes.append(f"{summed} pairs of γ rays of one cascade summed in a crystal.")
    if cex is not None:
        notes.append(f"Each excitation decays by the level scheme's cascade: {len(cex.transitions())} transitions "
                     "followed, each with its own orientation; the angular correlation between successive γ rays of "
                     "one cascade is not included.")
    cols = {"event": c["event"][rows], "particle": rows, "cascade": cascade, "initial": initial, "final": final,
            "crystal": which, "energy0": e0s,
            "energy_lab": energy_lab, "deposited": deposited, "measured": measured, "counted": counted,
            "theta": np.degrees(np.arccos(np.clip(lab[:, 2], -1, 1))),
            "phi": np.degrees(np.arctan2(lab[:, 1], lab[:, 0])), "beta": beta,
            "corrected_projectile": corrected["ejectile"], "corrected_recoil": corrected["recoil"],
            "weight": weight}
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
    top = max(1.3 * e0, max(x[0] for x in ROOM_LINES) * 1e-3 * 1.05 if room > 0 else 0.0,
              max((_q(e).to("MeV") for e, _ in extra), default=0.0) * 1.1) + 0.05
    edges = np.arange(0.0, top, bin_kev * 1e-3)
    singles_rate, singles_spectrum = {}, {}
    for cr in crystals:
        r = cr.response
        gamma_rate = all_excitations * float(r.total_efficiency(e0, cr.element))
        s = gamma_rate * r.shape(e0, edges)
        bg, bg_rate = _room_spectrum(experiment, cr, edges, room, extra)
        singles_rate[cr.name] = {"reaction": gamma_rate, "background": bg_rate}
        singles_spectrum[cr.name] = s + bg
    total_rate = sum(particle_rate.values()) + sum(sum(v.values()) for v in singles_rate.values())
    return {"window_s": window, "dead_time_s": dead, "live_fraction": 1 / (1 + dead * total_rate),
            "particle_rate": particle_rate, "singles_rate": singles_rate, "singles_spectrum": singles_spectrum,
            "singles_edges": edges, "excitation_rate": excitation_rate, "all_excitations": all_excitations}


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
        w = np.zeros(len(idx))
        for (k, q), v in a_kq.items():
            w += (v[idx] * math.sqrt(4 * math.pi / (2 * k + 1)) * angular.spherical_harmonic(k, q, th, ph)).real
        accept = rng.uniform(0.0, 1.0, len(idx)) * bound[idx] <= w
        hit = idx[accept]
        local[hit] = np.c_[np.sin(th[accept]) * np.cos(ph[accept]), np.sin(th[accept]) * np.sin(ph[accept]),
                           u[accept]]
        todo[hit] = False
    return local


def _cascade_steps(cex, theta_cm: np.ndarray, rng) -> tuple:
    """The decay of each excitation through the level scheme: the level it starts in, drawn from the direct
    excitation probabilities at its CM angle, then each step down by the scheme's branches (a converted transition
    gives no γ ray). Returns (index of the excitation each γ ray belongs to, initial level, final level, energy in
    the nucleus's frame in MeV), one entry per γ ray emitted."""
    grid, states = cex.table()
    starts = sorted(cex.paths)
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
        pick = np.array([np.interp(v[i], cdf[:, i], s[:, i]) for i in range(len(v))])
        out_c[~flat] = pick * ec[~flat]
        out[compton] = out_c
    kind = np.zeros(len(e), dtype=int)
    kind[single], kind[double], kind[compton] = 1, 2, 3
    return out, kind


__all__ = ["Crystal", "GAMMA_COLUMNS", "GammaEvents", "ROOM_LINES", "backgrounds", "crystals_of", "doppler_correct",
           "recorrect", "simulate_gammas"]
