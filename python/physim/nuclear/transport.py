"""Recoil transport for lifetime measurements: an excited nucleus slowing down (and scattering) through a stack of
layers and vacuum gaps until it decays.

A plunger (RDDS) is a target, a vacuum gap and a stopper (or degraders); a DSAM target is a target on a backing.
The emitter starts at a depth of the target with an energy and a direction; :func:`slow` says where it is, how fast
it moves and in which direction when it decays at a given (lab) time::

    from physim.nuclear.transport import Segment, slow
    stack = [Segment(mg, 0.7, name="target"), Segment(None, gap_um=50.0, name="gap"), Segment(cu, 15.7, name="stopper")]
    s = slow(stack, "98Ru", m_ru, energy, directions, depth, t_decay_ps, scatter=True, rng=1)
    s.beta, s.direction, s.segment, s.stopped

Units: time ps, matter thickness and depth mg/cm² (along the stack's normal, the z axis), gaps µm, energies MeV,
speeds as β. Arrays have one entry per recoil.

**Slowing down.** Without scattering (``scatter=False``) the recoil keeps its direction and loses energy as the
stopping tables say (:class:`~physim.nuclear.stopping.Stopping`, electronic and nuclear): its time to slow from E is
t(E) = ∫ dR / (ρ v), the clock of backlog item 73, and this is the arithmetic of ``gamma_events._slowing`` before it
moved here (the same results to rounding: the clock is now kept in ps). The starting depth counts the matter from
the front of the stack, so a nucleus made in the backing starts in the backing (``_slowing`` took every depth as
one in the target). In a vacuum gap the speed is constant: time = path / v along the recoil's direction.

**Time dilation (where γ enters).** The clocks above integrate dR / v with the lab speed, so they are lab time, and
``t_decay_ps`` is lab time: :func:`slow` needs no γ. γ enters where a decay time is drawn from a lifetime: the
nucleus decays after a proper time s drawn from exp(−s/τ), and the lab time of that decay is t = ∫ γ ds along the
path (γ falls as the recoil slows; it is 1 once it has stopped). :func:`lab_time` makes that conversion along the
same path (a proper clock, ∫ dR / (ρ v γ), in matter; t = γ s in a gap), and :class:`Kernel` carries the proper time
beside each lab time. At the speeds of Coulomb-barrier recoils (β ≈ 2–5%) γ − 1 is 2·10⁻⁴ to 1.3·10⁻³.

**Angular scattering** (``scatter=True``). A Monte Carlo of the recoil's nuclear collisions as it slows, as
LINESHAPE and TRIM do, vectorised over recoils. Each step crosses a path L(E) over which the recoil loses a fraction
:data:`STEP_FRACTION` of its energy (never shorter than the mean distance between atoms, N^(−1/3)). Over that path
it meets exactly one atom with an impact parameter p ≤ p_max, π p_max² N L = 1, drawn uniformly in p² (TRIM's free-
flight convention); the scattering angle comes from the ZBL universal potential (a table of the classical
scattering integral, Θ(ε, b) in reduced energy and impact parameter), the energy it gives the atom is
T = γ_m E sin²(Θ/2), and the recoil turns by the lab angle tan ψ = sin Θ / (cos Θ + M₁/M₂). The collisions with
p > p_max (many, each small) are not sampled: they turn the recoil by a Gaussian of variance N L ∫ ψ² 2πp dp from
p_max up, and their energy loss is continuous. The collision is taken at the step's mid energy. The energy
bookkeeping follows the stopping tables without double counting: electronic stopping is continuous, and the
continuous nuclear part is Stopping's nuclear stopping minus the mean loss of exactly the collisions the Monte Carlo
draws (N ∫₀^p_max T 2πp dp over the same interpolated table), so at every energy a recoil loses S_e + S_n per unit
path on average, as in :class:`~physim.nuclear.stopping.Stopping` (where the sampled part alone would exceed S_n,
within a few parts in a thousand at the lowest energies, the continuous part is zero). Below :data:`BETA_CUT` the
recoil stops along its last direction by the stopping tables. Electronic energy-loss straggling is left out (the
energy spread comes from the nuclear collisions). The random numbers come from ``rng`` (a ``numpy`` Generator or a
seed; None is seed 0), so a run is reproducible. Step size: the range, its straggling, the lateral spread and the
angles change by about 1% between steps of 3% and 1.5% of the energy.

**The DSAM kernel.** :func:`kernel` follows ``n`` sample recoils and records each one's velocity at a list of lab
times: the distribution of β·cos(angle to a detector) against time, which does not depend on the lifetime (AhKin's
idea), so a line-shape fit only re-weights it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Optional, Sequence

import numpy as np

from . import data
from .rates import MAX_PATH_FACTOR, _after, stopping
from .stopping import A0_CM, E2_MEV_CM

#: The speed of light, cm/ps and µm/ps.
C_CM_PS = 2.99792458e-2
C_UM_PS = 299.792458
#: Energy a recoil loses (as a fraction) over one Monte Carlo step.
STEP_FRACTION = 0.03
#: Below this speed a scattering recoil stops along its last direction by the stopping tables (a Doppler shift of
#: 0.05%, well below a germanium detector's resolution).
BETA_CUT = 5e-4


# ---------------------------------------------------------------------------------------------------------------
# The stack


@dataclass(frozen=True)
class Segment:
    """A layer of matter, or a vacuum gap when ``material`` is None."""

    material: object | None  # a physim.nuclear.data.Material with its density; None for vacuum
    thickness_mg_cm2: float = 0.0  # matter
    gap_um: float = 0.0  # vacuum
    name: str = ""  # "target", "backing", "gap", "stopper", ...

    @property
    def vacuum(self) -> bool:
        return self.material is None


@dataclass
class Slowed:
    """Where each recoil is at its decay."""

    energy: np.ndarray  # MeV at the decay
    direction: np.ndarray  # (n, 3) unit vectors at the decay (unchanged when scatter=False)
    segment: np.ndarray  # index of the segment the decay is in (len(segments) = beyond the last, −1 = out the front)
    stopped: np.ndarray  # bool: at rest before the decay
    beta: np.ndarray  # at the decay
    #: Beyond the contract: the proper time (ps) at the decay, and the position (µm, from the starting point; z
    #: along the stack's normal). For an infinite decay time, the position where the recoil stopped or left.
    proper_ps: Optional[np.ndarray] = None
    position_um: Optional[np.ndarray] = None


#: The most detector angles a :class:`Kernel` keeps :meth:`Kernel.beta_cos` for, and the most bytes it spends on
#: them (the latest angle is always kept): a fit asks for one angle hundreds of times, a ring scan a few angles.
KERNEL_ANGLES_KEPT = 8
KERNEL_ANGLE_BYTES = 128 * 2**20

_KERNEL_ARRAYS = ("times_ps", "beta", "direction", "proper_ps", "stopped", "segment", "beta0", "direction0")


@dataclass(frozen=True, eq=False)
class Kernel:
    """The recoils' velocities against lab time, for DSAM line shapes: ``n`` sample recoils followed through the
    stack, their state recorded at each of ``times_ps``. Arrays are (times, n) (directions (times, n, 3)).

    The distribution of β·cos(angle to a detector) at time t is :meth:`beta_cos`; a level of mean life τ decays at
    lab time t with probability density e^(−s/τ) / (τ γ), s the proper time (``proper_ps``) and γ = 1/√(1 − β²), so
    attenuations and line shapes are weighted sums over this table, with no transport per τ.

    A kernel does not change once made: its fields cannot be reassigned and its arrays are read-only (an array
    given as a view of another is copied first), so what is worked out from it once (:meth:`beta_cos` per angle
    here, the binned kernel in :func:`physim.nuclear.lifetimes.lineshape`) stays right. For different arrays, make
    a new kernel (``dataclasses.replace``)."""

    times_ps: np.ndarray  # (T,) lab times after the reaction
    beta: np.ndarray  # (T, n) speed (0 once stopped)
    direction: np.ndarray  # (T, n, 3) unit vectors
    proper_ps: np.ndarray  # (T, n) proper time at each lab time (time dilation)
    stopped: np.ndarray  # (T, n) bool: at rest
    segment: np.ndarray  # (T, n) the segment each sample is in
    beta0: np.ndarray  # (n,) speed at the start
    direction0: np.ndarray  # (n, 3) direction at the start
    _angles: dict = field(default_factory=dict, init=False, repr=False, compare=False)  # (θ, φ) → β·cos ψ

    def __post_init__(self):
        self._freeze()

    def _freeze(self) -> None:
        for name in _KERNEL_ARRAYS:
            a = np.asarray(getattr(self, name))
            if not a.flags.owndata:  # a view: its base could still be written through
                a = a.copy()
            a.flags.writeable = False
            object.__setattr__(self, name, a)

    def __getstate__(self) -> dict:
        state = dict(self.__dict__)
        state["_angles"] = {}
        return state

    def __setstate__(self, state: dict) -> None:
        for name, value in state.items():
            object.__setattr__(self, name, value)
        object.__setattr__(self, "_angles", {})
        self._freeze()

    @property
    def n(self) -> int:
        return self.beta.shape[1]

    def beta_cos(self, theta_deg: float, phi_deg: float = 0.0) -> np.ndarray:
        """β·cos(angle between the recoil's velocity and a detector at (θ, φ) from the beam axis), (T, n),
        read-only. Kept for the last few angles asked (:data:`KERNEL_ANGLES_KEPT`, :data:`KERNEL_ANGLE_BYTES`)."""
        key = (float(theta_deg), float(phi_deg))
        kept = self._angles
        out = kept.pop(key, None)  # put back last: the most recently used
        if out is None:
            th, ph = math.radians(key[0]), math.radians(key[1])
            axis = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
            out = self.beta * (self.direction @ axis)
            out.flags.writeable = False
            while kept and (len(kept) >= KERNEL_ANGLES_KEPT
                            or sum(a.nbytes for a in list(kept.values())) + out.nbytes > KERNEL_ANGLE_BYTES):
                kept.pop(next(iter(kept)), None)
        kept[key] = out
        return out

    def mean_beta_cos(self, theta_deg: float = 0.0, phi_deg: float = 0.0) -> np.ndarray:
        """The mean of :meth:`beta_cos` over the samples at each time, (T,)."""
        return self.beta_cos(theta_deg, phi_deg).mean(axis=1)


# ---------------------------------------------------------------------------------------------------------------
# Inputs


def _rng(rng) -> np.random.Generator:
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(0 if rng is None else rng)


def _beta_gamma(energy, mass_mev: float) -> tuple:
    e = np.maximum(np.asarray(energy, dtype=float), 0.0)
    gam = 1.0 + e / mass_mev
    return np.sqrt(np.maximum(1.0 - 1.0 / gam**2, 0.0)), gam


def _inputs(energy, direction, depth, *times) -> tuple:
    """Broadcast the per-recoil inputs: energy (n,), unit directions (n, 3), depth (n,) and each of ``times``."""
    e = np.atleast_1d(np.asarray(energy, dtype=float))
    u = np.asarray(direction, dtype=float)
    u = u.reshape(1, 3) if u.ndim == 1 else u
    z = np.atleast_1d(np.asarray(depth, dtype=float))
    ts = [np.atleast_1d(np.asarray(t, dtype=float)) for t in times]
    n = max(len(e), len(u), len(z), *(len(t) for t in ts))
    e = np.broadcast_to(e, (n,)).copy()
    u = np.broadcast_to(u, (n, 3)).copy()
    norm = np.linalg.norm(u, axis=1)
    if np.any(norm == 0):
        raise ValueError("a direction is the zero vector")
    u /= norm[:, None]
    z = np.broadcast_to(z, (n,)).copy()
    return (e, u, z, *(np.broadcast_to(t, (n,)).copy() for t in ts))


def _into_stack(direction, tilt_deg: float) -> np.ndarray:
    """Lab directions in the frame of a stack tilted by ``tilt_deg`` about the vertical (y) axis, whose normal
    (sin t, 0, cos t) becomes z (the frame :func:`slow` works in); an untilted stack leaves them as they are."""
    u = np.asarray(direction, dtype=float)
    if not tilt_deg:
        return u
    t = math.radians(tilt_deg)
    c, s = math.cos(t), math.sin(t)
    u = u.reshape(-1, 3)
    return np.stack([c * u[:, 0] - s * u[:, 2], u[:, 1], s * u[:, 0] + c * u[:, 2]], axis=1)


class _Layout:
    """Thicknesses and where each recoil starts: the first layer of matter whose back face is at or below its
    starting depth (mg/cm² of matter from the front of the stack; gaps hold none), so a recoil made at the back
    face of the target starts there, before the gap."""

    def __init__(self, segments: Sequence[Segment]):
        self.segments = list(segments)
        if not self.segments:
            raise ValueError("the stack is empty")
        self.gap = np.array([s.material is None for s in self.segments])
        self.th = np.array([float(s.gap_um if s.material is None else s.thickness_mg_cm2) for s in self.segments])
        if np.any(self.th < 0):
            raise ValueError("a segment has a negative thickness")
        for s in self.segments:
            if s.material is not None and not getattr(s.material, "density_g_cm3", None):
                raise ValueError(f"lifetimes need the density of {s.material.name}, which is not known")
        self.matter = np.flatnonzero(~self.gap)
        if not len(self.matter):
            raise ValueError("the stack needs a layer of matter for the recoil to start in")
        mat_th = np.where(self.gap, 0.0, self.th)
        self.total = float(mat_th.sum())
        self.cum_after = np.cumsum(mat_th)
        self.cum_before = self.cum_after - mat_th

    def start(self, depth: np.ndarray) -> tuple:
        k = np.searchsorted(self.cum_after[self.matter], depth, side="left")
        j0 = self.matter[np.minimum(k, len(self.matter) - 1)]
        return j0, depth - self.cum_before[j0]


# ---------------------------------------------------------------------------------------------------------------
# Clocks: the time to stop from each energy of the stopping table


_CLOCKS: dict = {}


def _clock(species: str, material, density: float, mass_mev: float, proper: bool = False) -> tuple:
    """The time (ps) a nucleus of ``species`` takes to stop from each energy (MeV) of its stopping table in
    ``material`` of density ``density`` (g/cm³): t(E) = ∫ dR / (ρ v), lab time (the clock of backlog item 73); with
    ``proper`` the nucleus's own time, ∫ dR / (ρ v γ). (energies with 0 first, times with 0 first.) Kept once
    computed."""
    key = (species, material.name, density, mass_mev, proper)
    if key in _CLOCKS:
        return _CLOCKS[key]
    st = stopping(species, material)
    e, r = np.asarray(st._e, dtype=float), np.asarray(st._range, dtype=float)  # MeV, mg/cm²
    beta = np.sqrt(np.maximum(e * (e + 2 * mass_mev), 0.0)) / (e + mass_mev)
    speed = beta * (1.0 + e / mass_mev) if proper else beta
    mid = 0.5 * (speed[1:] + speed[:-1])
    dt = np.diff(r) / (1000.0 * density * C_CM_PS * np.maximum(mid, 1e-12))
    t = np.r_[r[0] / (1000.0 * density * C_CM_PS * max(speed[0], 1e-12)), 0.0]
    t = np.cumsum(np.r_[t[0], dt])
    _CLOCKS[key] = (np.r_[0.0, e], np.r_[0.0, t])
    return _CLOCKS[key]


# ---------------------------------------------------------------------------------------------------------------
# Straight paths (no scattering)


def _straight(segments, species: str, mass_mev: float, energy, direction, depth, t, clock: str = "lab") -> dict:
    """Follow each recoil along its direction until the time ``t`` on ``clock`` ("lab" or "proper"): the slowing of
    backlog item 73 (``gamma_events._slowing``, the same arithmetic) with vacuum gaps. Returns the energy at that
    time, segment, stopped, lab and proper time, position (µm), and the lab time each recoil stopped (inf if it did
    not stop before ``t``)."""
    lay = _Layout(segments)
    e, u, z, t = _inputs(energy, direction, depth, t)
    n = len(e)
    mass_mev = float(mass_mev)
    j0, d_in = lay.start(z)
    cos_n = u[:, 2]
    fwd = cos_n > 0
    with np.errstate(divide="ignore"):
        inv = 1.0 / np.abs(cos_n)
    proper = clock == "proper"
    elapsed, pel = np.zeros(n), np.zeros(n)
    done, stopped = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    out = np.full(n, np.nan)
    seg_out = np.full(n, -2)
    lab_at, prop_at = np.full(n, np.nan), np.full(n, np.nan)
    t_stop = np.full(n, np.inf)
    pos = np.zeros((n, 3))
    S = len(lay.segments)
    for forward, order in ((True, range(S)), (False, range(S - 1, -1, -1))):
        for j in order:
            th = lay.th[j]
            if forward:
                dz = np.where(j0 == j, th - d_in, np.where(j0 < j, th, 0.0))
            else:
                dz = np.where(j0 == j, d_in, np.where(j0 > j, th, 0.0))
            with np.errstate(invalid="ignore"):
                if lay.gap[j]:
                    path = np.where(dz > 0, dz * inv, 0.0)
                else:
                    path = np.where(dz > 0, np.minimum(dz * inv, MAX_PATH_FACTOR * lay.total), 0.0)
            m = (fwd == forward) & ~done & (path > 0)
            if not m.any():
                continue
            idx = np.flatnonzero(m)
            if lay.gap[j]:
                _gap_step(idx, path[idx], e, u, t, elapsed, pel, pos, done, out, seg_out, lab_at, prop_at, j,
                          mass_mev, proper)
                continue
            mat = lay.segments[j].material
            density = float(mat.density_g_cm3)
            st = stopping(species, mat)
            grid_e, grid_t = _clock(species, mat, density, mass_mev)
            _, grid_p = _clock(species, mat, density, mass_mev, proper=True)
            um_per_mg = 10.0 / density
            e_in = e[idx]
            e_out = _after(st, e_in, path[idx])
            t_in, t_out = np.interp(e_in, grid_e, grid_t), np.interp(e_out, grid_e, grid_t)
            p_in, p_out = np.interp(e_in, grid_e, grid_p), np.interp(e_out, grid_e, grid_p)
            if proper:
                here = t[idx] < pel[idx] + (p_in - p_out)
                at = np.interp(p_in[here] - (t[idx[here]] - pel[idx[here]]), grid_p, grid_e)
            else:
                here = t[idx] < elapsed[idx] + (t_in - t_out)
                at = np.interp(t_in[here] - (t[idx[here]] - elapsed[idx[here]]), grid_t, grid_e)
            h = idx[here]
            out[h] = at
            done[h] = True
            seg_out[h] = j
            lab_at[h] = elapsed[h] + t_in[here] - np.interp(at, grid_e, grid_t)
            prop_at[h] = pel[h] + p_in[here] - np.interp(at, grid_e, grid_p)
            pos[h] += u[h] * ((st.range(e_in[here]) - st.range(at)) * um_per_mg)[:, None]
            rest = (e_out <= 0) & ~here
            r = idx[rest]
            out[r] = 0.0
            stopped[r] = True
            done[r] = True
            seg_out[r] = j
            t_stop[r] = elapsed[r] + t_in[rest]
            p_stop = pel[r] + p_in[rest]
            if proper:
                prop_at[r] = t[r]
                lab_at[r] = t_stop[r] + (t[r] - p_stop)
            else:
                lab_at[r] = t[r]
                prop_at[r] = p_stop + (t[r] - t_stop[r])
            pos[r] += u[r] * (st.range(e_in[rest]) * um_per_mg)[:, None]
            go = ~here & ~rest
            g = idx[go]
            elapsed[g] += (t_in - t_out)[go]
            pel[g] += (p_in - p_out)[go]
            pos[g] += u[g] * (path[g] * um_per_mg)[:, None]
            e[g] = e_out[go]
    rest = ~done  # out of the stack: in flight at the exit energy
    out[rest] = e[rest]
    seg_out[rest] = np.where(fwd[rest], S, -1)
    beta, gam = _beta_gamma(e[rest], mass_mev)
    if proper:
        prop_at[rest] = t[rest]
        lab_at[rest] = elapsed[rest] + (t[rest] - pel[rest]) * gam
    else:
        lab_at[rest] = t[rest]
        with np.errstate(invalid="ignore"):
            prop_at[rest] = pel[rest] + (t[rest] - elapsed[rest]) / gam
    fly = np.where(np.isfinite(lab_at[rest]), lab_at[rest] - elapsed[rest], 0.0)  # inf: where it left
    pos[rest] += u[rest] * (beta * C_UM_PS * fly)[:, None]
    return {"energy": out, "direction": u, "segment": seg_out, "stopped": stopped, "lab_ps": lab_at,
            "proper_ps": prop_at, "position_um": pos, "t_stop_ps": t_stop}


def _gap_step(idx, path_um, e, u, t, elapsed, pel, pos, done, out, seg_out, lab_at, prop_at, j, mass_mev, proper):
    """Recoils ``idx`` cross a vacuum gap of ``path_um`` at constant speed (inf: they never leave it)."""
    beta, gam = _beta_gamma(e[idx], mass_mev)
    with np.errstate(divide="ignore", invalid="ignore"):
        dt = np.where(path_um > 0, path_um / (np.maximum(beta, 1e-300) * C_UM_PS), 0.0)
    clock, dclock = (pel, dt / gam) if proper else (elapsed, dt)
    here = t[idx] < clock[idx] + dclock
    h = idx[here]
    out[h] = e[h]
    done[h] = True
    seg_out[h] = j
    if proper:
        prop_at[h] = t[h]
        lab_at[h] = elapsed[h] + (t[h] - pel[h]) * gam[here]
    else:
        lab_at[h] = t[h]
        prop_at[h] = pel[h] + (t[h] - elapsed[h]) / gam[here]
    pos[h] += u[h] * (beta[here] * C_UM_PS * (lab_at[h] - elapsed[h]))[:, None]
    g = idx[~here]
    elapsed[g] += dt[~here]
    pel[g] += (dt / gam)[~here]
    pos[g] += u[g] * path_um[~here][:, None]


# ---------------------------------------------------------------------------------------------------------------
# The public functions


def slow(segments, species: str, mass_mev: float, energy, direction, start_depth_mg_cm2, t_decay_ps,
         *, scatter: bool = False, rng=None) -> Slowed:
    """Where each recoil is when it decays at lab time ``t_decay_ps`` after the reaction.

    ``segments`` is the stack (:class:`Segment`, in the order the beam meets them, normal along z); ``species`` the
    recoil (``"98Ru"``), ``mass_mev`` its mass (with its excitation); ``energy`` (MeV), ``direction`` (unit vectors,
    (n, 3) or one (3,)) and ``start_depth_mg_cm2`` (mg/cm² of matter from the front of the stack) where it starts.
    It slows along its path through the rest of the stack (vacuum gaps at constant speed), then flies at its exit
    energy; one that stops first decays at rest. Without ``scatter`` its direction is kept (exactly item 73's
    slowing); with it, the ZBL Monte Carlo (see the module notes) turns it, with random numbers from ``rng``."""
    if scatter:
        e, u, z, t = _inputs(energy, direction, start_depth_mg_cm2, t_decay_ps)
        r = _monte_carlo(segments, species, float(mass_mev), e, u, z, t[:, None], _rng(rng))
        res = {k: (v[:, 0] if isinstance(v, np.ndarray) and v.ndim >= 2 and k != "ledger" else v)
               for k, v in r.items()}
    else:
        res = _straight(segments, species, mass_mev, energy, direction, start_depth_mg_cm2, t_decay_ps)
    beta, _ = _beta_gamma(res["energy"], float(mass_mev))
    return Slowed(energy=res["energy"], direction=res["direction"], segment=res["segment"],
                  stopped=res["stopped"], beta=beta, proper_ps=res["proper_ps"], position_um=res["position_um"])


def time_to_stop(segments, species: str, mass_mev: float, energy, start_depth_mg_cm2) -> np.ndarray:
    """The lab time (ps) each recoil takes to stop, moving along the stack's normal (no scattering), gaps
    included; inf for one that leaves the stack. For planning (the stopper must stop the recoil; DSAM backings)."""
    e = np.atleast_1d(np.asarray(energy, dtype=float))
    res = _straight(segments, species, mass_mev, e, np.array([0.0, 0.0, 1.0]), start_depth_mg_cm2, np.inf)
    return res["t_stop_ps"]


def lab_time(segments, species: str, mass_mev: float, energy, direction, start_depth_mg_cm2,
             proper_ps) -> np.ndarray:
    """The lab time (ps) at which a recoil decays whose own (proper) decay time is ``proper_ps``: t = ∫ γ ds along
    its straight path (time dilation; γ = 1 once it has stopped). This is where γ enters: draw the decay times from
    the lifetime in proper time, convert them here, and give them to :func:`slow`."""
    return _straight(segments, species, mass_mev, energy, direction, start_depth_mg_cm2, proper_ps,
                     clock="proper")["lab_ps"]


def kernel(segments, species: str, mass_mev: float, energy, direction, start_depth_mg_cm2, times_ps,
           *, scatter: bool = True, rng=None, n: int = 20000) -> Kernel:
    """The DSAM kernel: ``n`` sample recoils, drawn from the given starting states (energy, direction, depth: one,
    or arrays to draw from), followed through the stack, with their velocities recorded at each lab time of
    ``times_ps`` (sorted here). See :class:`Kernel`. Memory grows as n × len(times_ps)."""
    g = _rng(rng)
    times = np.sort(np.atleast_1d(np.asarray(times_ps, dtype=float)))
    e, u, z = _inputs(energy, direction, start_depth_mg_cm2)
    m = len(e)
    pick = np.zeros(n, dtype=int) if m == 1 else g.choice(m, size=n, replace=m < n)
    e, u, z = e[pick], u[pick], z[pick]
    T = len(times)
    if scatter:
        r = _monte_carlo(segments, species, float(mass_mev), e, u, z, np.broadcast_to(times, (n, T)), g,
                         positions=False)
        energy_t, dir_t, prop, stop, seg = (r["energy"].T, np.transpose(r["direction"], (1, 0, 2)),
                                            r["proper_ps"].T, r["stopped"].T, r["segment"].T)
    else:
        r = _straight(segments, species, mass_mev, np.tile(e, T), np.tile(u, (T, 1)), np.tile(z, T),
                      np.repeat(times, n))
        energy_t, dir_t = r["energy"].reshape(T, n), r["direction"].reshape(T, n, 3)
        prop, stop, seg = r["proper_ps"].reshape(T, n), r["stopped"].reshape(T, n), r["segment"].reshape(T, n)
    beta, _ = _beta_gamma(energy_t, float(mass_mev))
    beta0, _ = _beta_gamma(e, float(mass_mev))
    return Kernel(times_ps=times, beta=beta.astype(np.float32), direction=dir_t.astype(np.float32),
                  proper_ps=prop.astype(np.float32), stopped=stop, segment=seg.astype(np.int16), beta0=beta0,
                  direction0=u)


# ---------------------------------------------------------------------------------------------------------------
# ZBL scattering: the classical scattering angle in the universal potential, tabulated once


#: The ZBL universal screening function φ(x) = Σ c e^(−d x), x = r / a_U.
_SCREEN_C = np.array([0.18175, 0.50986, 0.28022, 0.028171])
_SCREEN_D = np.array([3.1998, 0.94229, 0.40290, 0.20162])
#: The table's grid: reduced energy ε and reduced impact parameter b = p / a_U, uniform in their logarithms.
_LN_EPS = np.linspace(math.log(1e-5), math.log(1e5), 161)
_LN_B = np.linspace(math.log(1e-4), math.log(1e2), 241)


def _zbl_phi(x):
    x = np.asarray(x, dtype=float)
    return np.sum(_SCREEN_C * np.exp(-_SCREEN_D * x[..., None]), axis=-1)


def _scattering_angle(eps, b, phi=_zbl_phi, nodes: int = 64) -> np.ndarray:
    """CM scattering angle Θ (rad) in a screened Coulomb potential V = Z₁Z₂e² φ(r/a)/r at reduced CM energy ε and
    impact parameter b: Θ = π − 2b ∫ dx / (x² √g), g = 1 − φ(x)/(x ε) − b²/x², from the distance of closest
    approach x₀ (g(x₀) = 0, by bisection), with u = x₀/x = 1 − t² to take out the root and Gauss–Legendre in t."""
    eps, b = np.broadcast_arrays(np.asarray(eps, dtype=float), np.asarray(b, dtype=float))
    shape = eps.shape
    eps, b = eps.ravel(), b.ravel()
    lo = b.copy()
    hi = np.maximum(0.5 * (1 / eps + np.sqrt(1 / eps**2 + 4 * b**2)), b)  # the bare-Coulomb root: g(hi) ≥ 0
    for _ in range(100):
        mid = np.sqrt(lo * hi)
        neg = 1 - phi(mid) / (mid * eps) - (b / mid) ** 2 < 0
        lo, hi = np.where(neg, mid, lo), np.where(neg, hi, mid)
    x0 = hi
    t, w = np.polynomial.legendre.leggauss(nodes)
    t, w = 0.5 * (t + 1), 0.5 * w
    uu = 1 - t**2
    x = x0[:, None] / uu
    g = 1 - uu * phi(x) / (x0[:, None] * eps[:, None]) - (b[:, None] * uu / x0[:, None]) ** 2
    integral = np.sum(w * t / np.sqrt(np.maximum(g, 1e-300)), axis=1)
    return (math.pi - 4 * (b / x0) * integral).reshape(shape)


@lru_cache(maxsize=1)
def _tables() -> dict:
    """Θ, ln(1 − cos Θ) on (ln ε, ln b) and I(ε, b) = ∫₀^b sin²(Θ/2) b' db' (the nuclear stopping of collisions inside b, reduced:
    the ZBL reduced stopping is s_n = 2 ε I(ε, ∞)) on the universal grid."""
    le, lb = np.meshgrid(_LN_EPS, _LN_B, indexing="ij")
    theta = np.clip(_scattering_angle(np.exp(le), np.exp(lb)), 0.0, math.pi)
    b = np.exp(_LN_B)
    f = np.sin(theta / 2) ** 2 * b**2  # integrand of ∫ ... b db in ln b
    dlb = _LN_B[1] - _LN_B[0]
    inner = np.concatenate([f[:, :1] / 2, f[:, :1] / 2 + np.cumsum(0.5 * (f[:, 1:] + f[:, :-1]) * dlb, axis=1)],
                           axis=1)
    # ln(1 − cos Θ) = ln(2 sin²(Θ/2)), interpolated for the collisions (smooth from Θ ≈ 0, where it is 2 ln Θ − ln 2,
    # to Θ = π), from which cos Θ, sin Θ and the energy given to the atom follow without trigonometry.
    ln_vers = np.log(np.maximum(2 * np.sin(theta / 2) ** 2, 1e-300))
    return {"theta": theta, "ln_vers": ln_vers, "inner": inner}


@lru_cache(maxsize=64)
def _tail(mass_ratio: float) -> np.ndarray:
    """C(ε, b) = ∫_b^∞ ψ² b' db' for an ion of mass M₁ = ``mass_ratio`` × M₂: the mean square lab deflection of
    the collisions outside b (reduced)."""
    theta = _tables()["theta"]
    psi = np.arctan2(np.sin(theta), np.cos(theta) + mass_ratio)
    f = psi**2 * np.exp(_LN_B) ** 2
    dlb = _LN_B[1] - _LN_B[0]
    seg = 0.5 * (f[:, 1:] + f[:, :-1]) * dlb
    return np.concatenate([np.cumsum(seg[:, ::-1], axis=1)[:, ::-1], np.zeros((len(_LN_EPS), 1))], axis=1)


def _bilinear(table: np.ndarray, ln_eps, ln_b) -> np.ndarray:
    """``table`` on the universal grid at (ln ε, ln b), linearly in both (clamped at the edges)."""
    x = np.clip((ln_eps - _LN_EPS[0]) / (_LN_EPS[1] - _LN_EPS[0]), 0.0, len(_LN_EPS) - 1.000001)
    y = np.clip((ln_b - _LN_B[0]) / (_LN_B[1] - _LN_B[0]), 0.0, len(_LN_B) - 1.000001)
    i, j = x.astype(np.intp), y.astype(np.intp)
    fx, fy = x - i, y - j
    return ((table[i, j] * (1 - fy) + table[i, j + 1] * fy) * (1 - fx)
            + (table[i + 1, j] * (1 - fy) + table[i + 1, j + 1] * fy) * fx)


# ---------------------------------------------------------------------------------------------------------------
# The Monte Carlo


class _Medium:
    """One recoil species in one material, tabulated on the stopping table's energy grid (uniform in ln E) for the
    Monte Carlo: the step L(E), p_max, the continuous stopping (electronic plus the unsampled nuclear part), its
    range, lab and proper clocks and electronic share, the variance of the unsampled deflections, and what one whole
    step from each energy does (so an ordinary step is one table lookup)."""

    def __init__(self, species: str, material, mass_mev: float, step_fraction: Optional[float] = None):
        step_fraction = STEP_FRACTION if step_fraction is None else step_fraction
        st = stopping(species, material)
        self.st = st
        self.rho = rho = float(material.density_g_cm3)
        self.um_per_mg = 10.0 / rho
        e = np.asarray(st._e, dtype=float)
        self.e = e
        self.ln_e = np.log(e)
        self.ln0, self.dln = float(self.ln_e[0]), float(self.ln_e[1] - self.ln_e[0])
        se, sn = np.asarray(st._se) * 1e-3, np.asarray(st._sn) * 1e-3  # MeV/(mg/cm²)
        stot = se + sn
        # The material's elements: atoms per gram and mean mass, by Z.
        by_z: dict = {}
        for z, a, k in material.atoms:
            w = by_z.setdefault(z, [0.0, 0.0])
            w[0] += k
            w[1] += k * data.nuclide((z, a)).atomic_mass_u
        z2 = np.array(sorted(by_z), dtype=float)
        count = np.array([by_z[int(z)][0] for z in z2])
        m2 = np.array([by_z[int(z)][1] / by_z[int(z)][0] for z in z2])
        apg = count * data.AVOGADRO / material.molar_mass  # atoms per gram
        z1, m1 = st.ion.Z, st.mass_u
        n_atoms = rho * apg.sum()  # per cm³
        l_min = n_atoms ** (-1 / 3) * rho * 1e3  # mg/cm²: the mean distance between atoms
        self.step = np.maximum(step_fraction * e / stot, l_min)
        pmax = 1.0 / np.sqrt(math.pi * apg.sum() * self.step * 1e-3)  # cm: π p_max² × atoms per cm² = 1
        self.ln_pmax = np.log(pmax)
        self.a = 0.8854 * A0_CM / (z1**0.23 + z2**0.23)
        self.ceps = self.a * m2 / (z1 * z2 * E2_MEV_CM * (m1 + m2))  # ε per MeV (lab)
        self.ln_ceps = np.log(self.ceps)
        self.ln_a = np.log(self.a)
        self.gm = 4 * m1 * m2 / (m1 + m2) ** 2
        self.mratio = m1 / m2
        self.cumprob = np.cumsum(apg) / apg.sum()
        sampled = np.zeros_like(e)
        ktail = np.zeros_like(e)
        # The sampled collisions' mean energy loss is that of exactly what the Monte Carlo draws (the same
        # interpolated table, b² uniform up to b_max), so that the continuous nuclear part makes up the difference
        # to Stopping's S_n precisely.
        # Mean of (1 − cos Θ) over u = b²/b_max² uniform, as ∫ vers u d ln u (in ln b, since at high ε it all sits
        # at u ≪ 1), from 16 e-folds of b below b_max.
        s_b = np.linspace(-16.0, 0.0, 641)
        for i in range(len(z2)):
            le, lb = np.log(self.ceps[i] * e), np.log(pmax / self.a[i])
            area = 2 * math.pi * self.a[i] ** 2 * apg[i] * 1e-3  # per (mg/cm²)
            vers = np.minimum(np.exp(_bilinear(_tables()["ln_vers"], le[:, None], lb[:, None] + s_b)), 2.0)
            f = vers * np.exp(2 * s_b) * 2.0  # vers × u per unit ln b
            mean_vers = np.sum(0.5 * (f[:, 1:] + f[:, :-1]), axis=1) * (s_b[1] - s_b[0]) + f[:, 0] / 2
            sampled += (apg[i] / apg.sum()) * self.gm[i] * e * 0.5 * mean_vers / self.step
            ktail += area * _bilinear(_tail(float(self.mratio[i])), le, lb)
        self.sn, self.sn_sampled = sn, sampled
        self.ktail = ktail  # rad² per (mg/cm²)
        sc = se + np.maximum(sn - sampled, 0.0)
        # Continuous range, electronic share and clocks (as Stopping integrates its range: in ln E, S ∝ √E below).
        dln = np.diff(self.ln_e)
        rc = np.r_[0.0, np.cumsum(0.5 * (e[1:] / sc[1:] + e[:-1] / sc[:-1]) * dln)] + 2 * e[0] / sc[0]
        self.ln_rc = np.log(rc)
        fe = e * se / sc
        self.qe = np.r_[0.0, np.cumsum(0.5 * (fe[1:] + fe[:-1]) * dln)] + e[0] * se[0] / sc[0]
        beta, gam = _beta_gamma(e, mass_mev)
        for name, speed in (("tc", beta), ("pc", beta * gam)):
            mid = 0.5 * (speed[1:] + speed[:-1])
            t0 = rc[0] / (1e3 * rho * C_CM_PS * max(speed[0], 1e-12))
            setattr(self, name, np.cumsum(np.r_[t0, np.diff(rc) / (1e3 * rho * C_CM_PS * mid)]))
        # The full stopping, for the last stretch below the cut.
        self.clock_e, self.clock_t = _clock(species, material, rho, mass_mev)
        _, self.clock_p = _clock(species, material, rho, mass_mev, proper=True)
        # One whole step from each energy of the grid: the energy after it, and its lab time, proper time and
        # electronic energy loss.
        # electronic energy loss, each integrated over the step itself (Simpson in ln E between the two energies):
        # differences of the cumulative tables would carry their interpolation error, about (grid step)² × E, into
        # every step.
        ln_after = self.ln_e_at_ln_rc(np.log(np.maximum(rc - self.step, 1e-300)))
        g_t = e / (sc * 1e3 * rho * C_CM_PS * np.maximum(beta, 1e-300))  # dt/d ln E
        simpson = np.r_[1.0, np.tile([4.0, 2.0], 8)[:-1], 1.0] / 48.0  # 16 intervals on [0, 1]
        pts = ln_after[:, None] + np.linspace(0.0, 1.0, 17) * (self.ln_e - ln_after)[:, None]
        width = (self.ln_e - ln_after)[:, None]

        def over_step(g):
            return np.sum(np.exp(np.interp(pts, self.ln_e, np.log(g))) * simpson * width, axis=1)

        dt, dp, dq = over_step(g_t), over_step(g_t / gam), over_step(e * se / sc)
        # Columns read at the start of every step (one gather each for the lower and upper grid node).
        self.cols = np.stack([self.step, ln_after, dt, dp, dq, self.ln_pmax, self.ktail], 1)
        self.dcols = np.diff(self.cols, axis=0)

    def columns(self, ln_e) -> np.ndarray:
        """:attr:`cols` at each ln E (linearly in ln E), (n, columns)."""
        i, w = self.index(ln_e)
        return self.cols[i] + self.dcols[i] * w[:, None]

    def index(self, ln_e) -> tuple:
        x = np.clip((ln_e - self.ln0) / self.dln, 0.0, len(self.e) - 1.000001)
        i = x.astype(np.intp)
        return i, x - i

    @staticmethod
    def at(table, i, w):
        return table[i] * (1 - w) + table[i + 1] * w

    def ln_e_at_ln_rc(self, ln_r):
        return np.interp(ln_r, self.ln_rc, self.ln_e)

    def ln_e_at_tc(self, t):
        return np.interp(t, self.tc, self.ln_e)


def _rotate(u: np.ndarray, cos_t, sin_t, cos_p, sin_p, normalise: bool = True) -> np.ndarray:
    """Turn unit vectors ``u`` by the polar angle (cos, sin) about themselves, at the azimuth (cos, sin)."""
    ux, uy, uz = u[:, 0], u[:, 1], u[:, 2]
    perp = np.sqrt(np.maximum(1.0 - uz * uz, 0.0))
    near = perp < 1e-10
    q = sin_t / np.where(near, 1.0, perp)
    out = np.empty_like(u)
    out[:, 0] = ux * cos_t + q * (ux * uz * cos_p - uy * sin_p)
    out[:, 1] = uy * cos_t + q * (uy * uz * cos_p + ux * sin_p)
    out[:, 2] = uz * cos_t - sin_t * perp * cos_p
    if near.any():
        k = np.flatnonzero(near)
        sgn = np.where(uz[k] >= 0, 1.0, -1.0)
        st = np.broadcast_to(sin_t, uz.shape)[k]
        ct = np.broadcast_to(cos_t, uz.shape)[k]
        out[k, 0], out[k, 1], out[k, 2] = st * cos_p[k], st * sin_p[k], sgn * ct
    if normalise:
        out /= np.sqrt(np.einsum("ij,ij->i", out, out))[:, None]
    return out


def _azimuths(rng, n: int) -> tuple:
    phi = 2 * math.pi * rng.random(n)
    return np.cos(phi), np.sin(phi)


#: The per-recoil energy ledger the Monte Carlo keeps (MeV, and the path length in mg/cm²).
_LEDGER = ("electronic", "nuclear_continuous", "nuclear_collisions", "below_cut", "path_mg_cm2")


class _Run:
    """One Monte Carlo run: the state of the recoils still followed (arrays in a dict, compacted as recoils
    finish), and what is recorded at each observation time, (n, K)."""

    def __init__(self, segments, species, mass_mev, e, u, z, times, rng, positions):
        self.lay = lay = _Layout(segments)
        self.S = len(lay.segments)
        self.mass = mass_mev
        self.rng = rng
        self.media = {int(j): _Medium(species, lay.segments[j].material, mass_mev) for j in lay.matter}
        first = self.media[int(lay.matter[0])]
        self.e_cut = max(1.1 * float(first.e[0]), 0.5 * mass_mev * BETA_CUT**2)
        n, K = times.shape
        self.times, self.K = times, K
        self.rec_e = np.zeros((n, K))
        self.rec_u = np.zeros((n, K, 3))
        self.rec_seg = np.zeros((n, K), dtype=np.int64)
        self.rec_stop = np.zeros((n, K), dtype=bool)
        self.rec_p = np.zeros((n, K))
        self.rec_pos = np.zeros((n, K, 3)) if positions else None
        self.ledger = {k: np.zeros(n) for k in _LEDGER}
        j0, d0 = lay.start(z)
        self.s = {"ids": np.arange(n), "E": e.copy(), "u": u.copy(), "seg": j0.astype(np.int64), "d": d0,
                  "pos": np.zeros((n, 3)), "t": np.zeros(n), "p": np.zeros(n), "k": np.zeros(n, dtype=np.int64),
                  "next": times[:, 0].copy() if K else np.full(n, np.inf)}
        for key in _LEDGER:
            self.s[key] = np.zeros(n)

    # -- recording ---------------------------------------------------------------------------------------------

    def _record(self, sub, rows, tk, state_at):
        e, u, seg, stopped, p, pos = state_at(rows, tk)
        ids, k = sub["ids"][rows], sub["k"][rows]
        self.rec_e[ids, k] = e
        self.rec_u[ids, k] = u
        self.rec_seg[ids, k] = seg
        self.rec_stop[ids, k] = stopped
        self.rec_p[ids, k] = p
        if self.rec_pos is not None:
            self.rec_pos[ids, k] = pos
        k = k + 1
        sub["k"][rows] = k
        more = k < self.K
        sub["next"][rows] = np.where(more, self.times[ids, np.minimum(k, self.K - 1)], np.inf)

    def _due(self, sub, rows, t_end, state_at):
        """Record every observation of ``rows`` of ``sub`` before ``t_end`` (a value per row of ``sub``), with
        ``state_at(rows, observation times) -> (e, u, seg, stopped, p, pos)``."""
        nxt = sub["next"]
        while len(rows):
            rows = rows[nxt[rows] < t_end[rows]]
            if len(rows):
                self._record(sub, rows, nxt[rows], state_at)

    def _flush(self, sub, rows, state_at):
        """Record all remaining observations of ``rows`` (the recoil is no longer followed)."""
        k = sub["k"]
        rows = rows[k[rows] < self.K]
        while len(rows):
            self._record(sub, rows, self.times[sub["ids"][rows], k[rows]], state_at)
            rows = rows[k[rows] < self.K]

    # -- the end of a track ------------------------------------------------------------------------------------

    def _leave(self, sub, rows):
        """Recoils ``rows`` have left the stack: they fly on at constant velocity."""
        if not len(rows):
            return
        beta, gam = _beta_gamma(sub["E"], self.mass)
        where = np.where(sub["seg"] >= self.S, self.S, -1)

        def state(rr, tk):
            dt = np.where(np.isfinite(tk), tk - sub["t"][rr], 0.0)  # an infinite time: where it left
            return (sub["E"][rr], sub["u"][rr], where[rr], False, sub["p"][rr] + dt / gam[rr],
                    sub["pos"][rr] + sub["u"][rr] * (beta[rr] * C_UM_PS * dt)[:, None])

        self._flush(sub, rows, state)

    def _stop(self, sub, rows):
        """Recoils ``rows`` have fallen below the cut: they stop along their direction by the stopping tables
        (CSDA), with no more scattering."""
        if not len(rows):
            return
        segs = sub["seg"][rows]
        for j in np.unique(segs):
            rr = rows[segs == j]
            med = self.media[int(j)]
            e = sub["E"][rr]
            t_rem, p_rem, r_rem = (np.zeros(len(sub["E"])) for _ in range(3))
            t_rem[rr] = np.interp(e, med.clock_e, med.clock_t)
            p_rem[rr] = np.interp(e, med.clock_e, med.clock_p)
            r_rem[rr] = np.asarray(med.st.range(e))
            sub["below_cut"][rr] += e
            sub["path_mg_cm2"][rr] += r_rem[rr]

            def state(r2, tk, med=med, t_rem=t_rem, p_rem=p_rem, r_rem=r_rem, j=j):
                dt = tk - sub["t"][r2]
                moving = dt < t_rem[r2]
                e_obs = np.where(moving, np.interp(t_rem[r2] - np.where(moving, dt, 0.0), med.clock_t,
                                                   med.clock_e), 0.0)
                p = np.where(moving, sub["p"][r2] + p_rem[r2] - np.interp(e_obs, med.clock_e, med.clock_p),
                             sub["p"][r2] + p_rem[r2] + (dt - t_rem[r2]))
                path = r_rem[r2] - np.where(moving, np.asarray(med.st.range(e_obs)), 0.0)
                pos = sub["pos"][r2] + sub["u"][r2] * (path * med.um_per_mg)[:, None]
                return e_obs, sub["u"][r2], j, ~moving, p, pos

            self._flush(sub, rr, state)

    # -- vacuum gaps -------------------------------------------------------------------------------------------

    def _cross(self, sub, rows, forward):
        """Recoils ``rows`` have reached the far face of their segment (``forward``: moving +z): carry them into the
        next one, through any vacuum gaps, or out of the stack. Returns the rows still followed."""
        lay = self.lay
        sub["seg"][rows] += np.where(forward, 1, -1)
        alive = np.ones(len(rows), dtype=bool)
        while True:
            seg = sub["seg"][rows]
            out = alive & ((seg < 0) | (seg >= self.S))
            if out.any():
                self._leave(sub, rows[out])
                alive &= ~out
            seg_c = np.clip(seg, 0, self.S - 1)
            in_gap = alive & lay.gap[seg_c]
            matter = alive & ~lay.gap[seg_c]
            sub["d"][rows[matter]] = np.where(sub["u"][rows[matter], 2] > 0, 0.0, lay.th[seg_c[matter]])
            if not in_gap.any():
                break
            g = rows[in_gap]
            uz = sub["u"][g, 2]
            beta, gam = _beta_gamma(sub["E"], self.mass)
            with np.errstate(divide="ignore"):
                path = lay.th[sub["seg"][g]] / np.abs(uz)
                dt = path / (beta[g] * C_UM_PS)

            def state(rr, tk, beta=beta, gam=gam):
                dtk = np.where(np.isfinite(tk), tk - sub["t"][rr], 0.0)
                return (sub["E"][rr], sub["u"][rr], sub["seg"][rr], False, sub["p"][rr] + dtk / gam[rr],
                        sub["pos"][rr] + sub["u"][rr] * (beta[rr] * C_UM_PS * dtk)[:, None])

            t_end = np.full(len(sub["E"]), np.inf)
            t_end[g] = sub["t"][g] + dt
            self._due(sub, g, t_end, state)
            forever = ~np.isfinite(dt)  # flying along the gap: every remaining observation is in it
            if forever.any():
                self._flush(sub, g[forever], state)
            fin = np.isfinite(dt)
            gf = g[fin]
            sub["pos"][gf] += sub["u"][gf] * path[fin][:, None]
            sub["t"][gf] += dt[fin]
            sub["p"][gf] += (dt / gam[g])[fin]
            sub["seg"][gf] += np.where(uz[fin] > 0, 1, -1)
            alive &= sub["k"][rows] < self.K
        return rows[alive & (sub["k"][rows] < self.K)]

    # -- one step in matter ------------------------------------------------------------------------------------

    def _step(self, sub, j) -> np.ndarray:
        """One step of every recoil in ``sub`` (all in matter segment ``j``); returns which are still followed."""
        med, th, rng = self.media[j], self.lay.th[j], self.rng
        E, u, d = sub["E"], sub["u"], sub["d"]
        n = len(E)
        lne = np.log(E)
        L, lne_a, dt, dp, dq, ln_pmax, ktail = med.columns(lne).T
        uz = u[:, 2].copy()
        with np.errstate(divide="ignore", invalid="ignore"):
            lb = np.where(uz > 0, (th - d) / uz, np.where(uz < 0, d / -uz, np.inf))
        trunc = lb < L
        ls = L.copy()
        if trunc.any():
            # A step cut short by the segment's face: the continuous loss over the shorter path, from the tables.
            c = np.flatnonzero(trunc)
            ls[c] = np.maximum(lb[c], 0.0)
            i, w = med.index(lne[c])
            ln_a = med.ln_e_at_ln_rc(np.log(np.maximum(np.exp(med.at(med.ln_rc, i, w)) - ls[c], 1e-300)))
            ia, wa = med.index(ln_a)
            lne_a[c] = ln_a
            dt[c] = med.at(med.tc, i, w) - med.at(med.tc, ia, wa)
            dp[c] = med.at(med.pc, i, w) - med.at(med.pc, ia, wa)
            dq[c] = med.at(med.qe, i, w) - med.at(med.qe, ia, wa)
        e_a = np.exp(lne_a)
        t_end = sub["t"] + dt
        due = np.flatnonzero(sub["next"] < t_end)
        if len(due):
            def state(rr, tk):
                i0, w0 = med.index(lne[rr])
                tc0, pc0, rc0 = med.at(med.tc, i0, w0), med.at(med.pc, i0, w0), np.exp(med.at(med.ln_rc, i0, w0))
                lne_o = med.ln_e_at_tc(tc0 - (tk - sub["t"][rr]))
                io, wo = med.index(lne_o)
                path = np.maximum(rc0 - np.exp(med.at(med.ln_rc, io, wo)), 0.0)
                return (np.exp(lne_o), u[rr], j, False, sub["p"][rr] + pc0 - med.at(med.pc, io, wo),
                        sub["pos"][rr] + u[rr] * (path * med.um_per_mg)[:, None])

            self._due(sub, due, t_end, state)
        d += uz * ls
        sub["pos"] += u * (ls * med.um_per_mg)[:, None]
        sub["t"] = t_end
        sub["p"] += dp
        sub["electronic"] += dq
        sub["nuclear_continuous"] += (E - e_a) - dq
        sub["path_mg_cm2"] += ls
        E[:] = e_a
        # The unsampled distant collisions: a Gaussian deflection δ (δ² is exponential with that variance).
        delta = np.sqrt(ktail * ls * rng.standard_exponential(n))
        cos_t, sin_t = np.cos(delta), np.sin(delta)
        # The one sampled collision of the step (with the chance ls/L on a step cut short by the segment's face).
        hit = ~trunc
        if trunc.any():
            c = np.flatnonzero(trunc)
            hit[c] = rng.random(len(c)) * L[c] < ls[c]
        hit = np.flatnonzero(hit)
        if len(hit):
            nc = len(hit)
            el = (0 if len(med.cumprob) == 1  # one element: its constants are scalars
                  else np.minimum(np.searchsorted(med.cumprob, rng.random(nc)), len(med.cumprob) - 1))
            ln_b = ln_pmax[hit] - med.ln_a[el] + 0.5 * np.log1p(-rng.random(nc))
            # The collision happens somewhere along the step: at the step's mid energy (no first-order bias from
            # the energy falling along the step), taking no more than the recoil has left.
            lne_mid = 0.5 * (lne[hit] + lne_a[hit])
            vers = np.exp(_bilinear(_tables()["ln_vers"], med.ln_ceps[el] + lne_mid, ln_b))  # 1 − cos Θ
            vers = np.minimum(vers, 2.0)
            transfer = np.minimum(med.gm[el] * np.exp(lne_mid) * (0.5 * vers), E[hit])
            E[hit] -= transfer
            sub["nuclear_collisions"][hit] += transfer
            # The lab angle: tan ψ = sin Θ / (cos Θ + M₁/M₂).
            x = 1.0 - vers + med.mratio[el]
            y = np.sqrt(vers * (2.0 - vers))
            r = np.sqrt(x * x + y * y)
            cos_psi, sin_psi = x / r, y / r
            # Both turns at once: two turns with independent uniform azimuths are one turn by the angle of the
            # spherical law of cosines, at a uniform azimuth.
            cos_a = np.cos(2 * math.pi * rng.random(nc))
            cos_c = np.clip(cos_psi * cos_t[hit] - sin_psi * sin_t[hit] * cos_a, -1.0, 1.0)
            cos_t[hit], sin_t[hit] = cos_c, np.sqrt(1.0 - cos_c * cos_c)
        u[:] = _rotate(u, cos_t, sin_t, *_azimuths(rng, n), normalise=False)
        alive = sub["k"] < self.K
        cut = alive & (E < self.e_cut)
        if cut.any():
            self._stop(sub, np.flatnonzero(cut))
            alive &= ~cut
        cross = np.flatnonzero(alive & trunc)
        if len(cross):
            keep = self._cross(sub, cross, uz[cross] > 0)
            alive[cross] = False
            alive[keep] = True
        return alive & (sub["k"] < self.K)

    def _retire(self, s, gone) -> None:
        for key in _LEDGER:
            self.ledger[key][s["ids"][gone]] = s[key][gone]

    def run(self) -> None:
        s = self.s
        low = np.flatnonzero(s["E"] < self.e_cut)  # recoils that start (almost) at rest
        if len(low):
            self._stop(s, low)
        keep = s["k"] < self.K
        self._retire(s, ~keep)
        s = {key: v[keep] for key, v in s.items()}
        steps = 0
        while len(s["ids"]):
            segs = s["seg"]
            j_lo, j_hi = int(segs.min()), int(segs.max())
            if j_lo == j_hi:
                alive = self._step(s, j_lo)
            else:
                alive = np.ones(len(s["ids"]), dtype=bool)
                for j in np.unique(segs):
                    m = segs == j
                    sub = {key: v[m] for key, v in s.items()}
                    a = self._step(sub, int(j))
                    for key in s:
                        s[key][m] = sub[key]
                    alive[m] = a
            if not alive.all():
                self._retire(s, ~alive)
                s = {key: v[alive] for key, v in s.items()}
            steps += 1
            if steps % 16 == 0:  # the turns keep |u| = 1 to rounding; put it back now and then
                s["u"] /= np.sqrt(np.einsum("ij,ij->i", s["u"], s["u"]))[:, None]


def _monte_carlo(segments, species: str, mass_mev: float, e, u, z, times, rng, positions: bool = True) -> dict:
    """The scattering Monte Carlo: each recoil's state at each of its observation times ``times`` (n, K, sorted
    along K). Returns energy, direction, segment, stopped, proper_ps, position_um (n, K, ...) and the energy
    ``ledger`` (per recoil: electronic, nuclear_continuous, nuclear_collisions, below_cut, and the path length)."""
    times = np.asarray(times, dtype=float)
    run = _Run(segments, species, mass_mev, e, u, z, times, rng, positions)
    run.run()
    return {"energy": run.rec_e, "direction": run.rec_u, "segment": run.rec_seg, "stopped": run.rec_stop,
            "proper_ps": run.rec_p, "position_um": run.rec_pos, "ledger": run.ledger}
