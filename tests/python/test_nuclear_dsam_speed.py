"""DSAM speed (lifetime wave 4C): what transport.Kernel, lifetimes.lineshape and residues.dsam_lineshapes keep
between calls (β·cos ψ per angle, the binned kernel, the cascade's emission weights) gives the same numbers as
working everything out afresh, and a changed kernel or cascade is never answered from what was kept."""

import dataclasses
import math
import pickle

import numpy as np
import pytest

from physim.nuclear import Experiment, data
from physim.nuclear import lifetime_analysis as la
from physim.nuclear import lifetimes as lt
from physim.nuclear import transport as tr
from physim.nuclear.residues import dsam_lineshapes


@pytest.fixture(scope="module")
def kern():
    """143Sm slowing and scattering in 13 mg/cm² Au (the fusion_sm143_dsam backing), recoiling in a cone."""
    sm = "143Sm"
    mass = data.nuclide(sm).nuclear_mass_mev
    rng = np.random.default_rng(3)
    m = 300
    energy = mass * (1 / math.sqrt(1 - 0.016**2) - 1) * rng.normal(1.0, 0.05, m)
    th, ph = np.radians(rng.uniform(0, 15, m)), rng.uniform(0, 2 * np.pi, m)
    u = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], axis=1)
    times = np.r_[0.0, np.geomspace(1e-3, 3.0, 59)]
    return tr.kernel([tr.Segment(data.material("Au"), 13.0, name="backing")], sm, mass, energy, u, 0.0, times,
                     rng=4, n=500)


class _Plain:
    """The kernel as it was before the cache: β·cos ψ worked out at every call, writable arrays."""

    def __init__(self, k):
        self.times_ps = np.array(k.times_ps)
        self.beta = np.array(k.beta)
        self.direction = np.array(k.direction)

    def beta_cos(self, theta_deg, phi_deg=0.0):
        th, ph = math.radians(theta_deg), math.radians(phi_deg)
        axis = np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
        return self.beta * (self.direction @ axis)


class _KeepNothing(dict):
    def __setitem__(self, key, value):
        pass


@pytest.fixture
def uncached(monkeypatch):
    """lifetimes keeping neither binned kernels nor emission weights."""
    monkeypatch.setattr(lt, "_MATRICES", _KeepNothing())
    monkeypatch.setattr(lt, "_WEIGHTS", _KeepNothing())


def _old_spread(lo, hi, weight, edges, group=None, n_groups=1):
    """lifetimes._spread as it was before wave 4C (the loop over offsets on every segment)."""
    nb = len(edges) - 1
    size = n_groups * nb
    group = np.zeros(len(lo), dtype=np.int64) if group is None else np.asarray(group, dtype=np.int64)
    jlo = np.searchsorted(edges, lo, side="right") - 1
    jhi = np.searchsorted(edges, hi, side="right") - 1
    length = hi - lo
    one = (jlo == jhi) | (length <= 1e-13)
    inside = one & (jlo >= 0) & (jlo < nb)
    out = np.bincount(group[inside] * nb + jlo[inside], weights=weight[inside], minlength=size).astype(float)
    many = np.flatnonzero(~one & (jhi >= 0) & (jlo < nb))
    if len(many):
        lo, length, weight, jlo, span, group = (lo[many], length[many], weight[many], jlo[many], (jhi - jlo)[many],
                                                group[many])
        for o in range(int(span.max()) + 1):
            j = jlo + o
            m = (o <= span) & (j >= 0) & (j < nb)
            if not m.any():
                continue
            jj = j[m]
            left = np.clip((edges[jj] - lo[m]) / length[m], 0.0, 1.0)
            right = np.clip((edges[jj + 1] - lo[m]) / length[m], 0.0, 1.0)
            out += np.bincount(group[m] * nb + jj, weights=weight[m] * (right - left), minlength=size)
    return out.reshape(n_groups, nb)


CASCADE = ([0.15, 0.5, 2.0], np.array([[0, 0.7, 0.3], [0, 0, 1.0], [0, 0, 0]]), [0.4, 0.4, 0.2])


def test_spread_is_the_old_loop_exactly():
    rng = np.random.default_rng(8)
    edges = np.linspace(0.97, 1.03, 400)
    a = rng.uniform(0.95, 1.05, 20000)
    b = a + rng.exponential(0.004, a.size) * (rng.random(a.size) < 0.8)
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    w = rng.random(a.size)
    g = rng.integers(0, 7, a.size)
    np.testing.assert_array_equal(lt._spread(lo, hi, w, edges, g, 7), _old_spread(lo, hi, w, edges, g, 7))
    np.testing.assert_array_equal(lt._spread(lo, hi, w, edges), _old_spread(lo, hi, w, edges))


def test_cached_lineshapes_are_the_uncached_ones(kern, monkeypatch):
    edges = np.linspace(980.0, 1020.0, 321)
    angles = [148.0, 30.0, 148.0, 90.0, 30.0, 148.0]
    cached = [lt.lineshape(1000.0, *CASCADE, kern, th, edges) for th in angles]
    monkeypatch.setattr(lt, "_MATRICES", _KeepNothing())
    monkeypatch.setattr(lt, "_WEIGHTS", _KeepNothing())
    plain = _Plain(kern)
    for th, got in zip(angles, cached):
        want = lt.lineshape(1000.0, *CASCADE, plain, th, edges)
        np.testing.assert_allclose(got, want, rtol=1e-10, atol=0)
        np.testing.assert_array_equal(got, want)  # in fact the same numbers


def _counts(kern, theta, e0=800.0):
    fine = np.arange(770.0, 830.0, 0.05)
    feed = la.Feeding.side(0.1, 0.3)
    shape = np.nan_to_num(lt.lineshape(e0, [0.1, 0.4], feed.feeding, feed.populations, kern, theta, fine)[1])
    edges = np.arange(780.0, 820.0, 0.5)
    mu = np.diff(np.interp(edges, fine, np.r_[0.0, np.cumsum(shape)])) * 8000 + 20.0
    return np.random.default_rng(12).poisson(mu).astype(float), edges, feed


def test_cached_fit_is_the_uncached_one(kern, monkeypatch):
    counts, edges, feed = _counts(kern, 145.0)
    got = la.fit_lineshape(counts, edges, 800.0, kern, 145.0, feed, fwhm_kev=2.0)  # with the stopping budget
    monkeypatch.setattr(lt, "_MATRICES", _KeepNothing())
    monkeypatch.setattr(lt, "_WEIGHTS", _KeepNothing())
    want = la.fit_lineshape(counts, edges, 800.0, _Plain(kern), 145.0, feed, fwhm_kev=2.0)
    assert got.tau_ps == pytest.approx(want.tau_ps, rel=1e-10)
    assert got.stat_ps == pytest.approx(want.stat_ps, rel=1e-10)
    assert got.budget["stopping"] == pytest.approx(want.budget["stopping"], rel=1e-10)
    assert 0.2 < got.tau_ps < 0.8


def test_the_kernel_cannot_change(kern):
    with pytest.raises(ValueError):
        kern.beta[0, 0] = 0.5
    with pytest.raises(ValueError):
        kern.direction[0, 0, 0] = 1.0
    with pytest.raises(dataclasses.FrozenInstanceError):
        kern.beta = np.zeros_like(kern.beta)
    # An array given as a view is copied: writing to its base afterwards does not reach the kernel.
    base = np.array(kern.beta)
    k2 = dataclasses.replace(kern, beta=base[:, :])
    before = k2.beta_cos(30.0).copy()
    base *= 2.0
    np.testing.assert_array_equal(k2.beta_cos(30.0), before)
    np.testing.assert_array_equal(k2.beta, kern.beta)


def test_beta_cos_is_kept_per_angle_and_bounded(kern, monkeypatch):
    k = dataclasses.replace(kern)  # its own, empty cache
    a = k.beta_cos(148.0)
    assert k.beta_cos(148.0) is a and not a.flags.writeable
    np.testing.assert_array_equal(a, _Plain(k).beta_cos(148.0))
    assert k.beta_cos(148.0, 10.0) is not a  # φ is part of the key
    np.testing.assert_array_equal(k.beta_cos(148.0, 10.0), _Plain(k).beta_cos(148.0, 10.0))
    for th in range(tr.KERNEL_ANGLES_KEPT + 3):
        k.beta_cos(float(th))
    assert len(k._angles) == tr.KERNEL_ANGLES_KEPT
    assert (float(tr.KERNEL_ANGLES_KEPT + 2), 0.0) in k._angles and (0.0, 0.0) not in k._angles
    monkeypatch.setattr(tr, "KERNEL_ANGLE_BYTES", int(1.5 * a.nbytes))
    k.beta_cos(1.5)
    k.beta_cos(2.5)
    assert list(k._angles) == [(2.5, 0.0)]
    # Pickled: the arrays come back read-only, the cache empty.
    back = pickle.loads(pickle.dumps(k))
    assert back._angles == {} and not back.beta.flags.writeable
    np.testing.assert_array_equal(back.beta_cos(2.5), k.beta_cos(2.5))


def test_a_changed_kernel_is_binned_again(kern, monkeypatch):
    edges = np.linspace(980.0, 1020.0, 321)
    old = lt.lineshape(1000.0, *CASCADE, kern, 148.0, edges)
    slower = dataclasses.replace(kern, beta=kern.beta * np.float32(0.6))
    new = lt.lineshape(1000.0, *CASCADE, slower, 148.0, edges)
    assert np.nanmax(np.abs(new - old)) > 1e-3
    # An array kernel whose array is replaced (a new object) is binned again too.
    times = np.linspace(0.0, 2.0, 41)
    arr = type("K", (), {})()
    arr.times_ps = times
    arr.beta_cos = 0.02 * np.exp(-times / 0.5)[:, None] * np.ones((1, 3))
    first = lt.lineshape(1000.0, [0.3], None, [1.0], arr, 0.0, edges)
    arr.beta_cos = 0.5 * arr.beta_cos
    second = lt.lineshape(1000.0, [0.3], None, [1.0], arr, 0.0, edges)
    monkeypatch.setattr(lt, "_MATRICES", _KeepNothing())
    monkeypatch.setattr(lt, "_WEIGHTS", _KeepNothing())
    np.testing.assert_array_equal(new, lt.lineshape(1000.0, *CASCADE, _Plain(slower), 148.0, edges))
    np.testing.assert_array_equal(second, lt.lineshape(1000.0, [0.3], None, [1.0], arr, 0.0, edges))
    assert np.nanmax(np.abs(second - first)) > 1e-3


def test_emission_weights_are_kept_per_cascade(kern):
    t = np.asarray(kern.times_ps, dtype=float)
    lam, f, p = lt._cascade(*CASCADE)
    one = lt._emission_weights(lam, f, p, t)
    assert lt._emission_weights(lam.copy(), f.copy(), p.copy(), t.copy())[0] is one[0]
    assert all(not a.flags.writeable for a in one)
    for got, want in zip(one, lt._emission_weights_of(lam, f, p, t)):
        np.testing.assert_array_equal(got, want)
    lam2, f2, p2 = lt._cascade([0.15, 0.5, 2.1], CASCADE[1], CASCADE[2])
    other = lt._emission_weights(lam2, f2, p2, t)
    assert other[0] is not one[0] and np.max(np.abs(other[0] - one[0])) > 0


def test_dsam_lineshapes_works_out_the_emission_weights_once(monkeypatch):
    """The cascade's emission weights do not depend on the ring: one Bateman calculation for three rings."""
    calls = []
    real = lt._emission_weights_of

    def counted(*args):
        calls.append(1)
        return real(*args)

    monkeypatch.setattr(lt, "_emission_weights_of", counted)
    monkeypatch.setattr(lt, "_WEIGHTS", {})
    exp = Experiment.example("fusion_sm143_dsam")  # the paper's level scheme typed in: no ENSDF needed
    out = dsam_lineshapes(exp, theta_deg=[148.0, 123.0, 90.0], n=300, events=2000, times=40, resolution=False)
    assert len(out["theta_deg"]) == 3 and out["lines"]
    assert len(calls) == 1
