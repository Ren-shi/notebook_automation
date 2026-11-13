"""The angular correlation between successive γ rays of one cascade (backlog item 71): the emission of one γ ray
along a direction gives the distribution of :meth:`Populated.w` and leaves the next level in the state whose
average over directions is the unobserved feeding; the γ rays of an unoriented 4⁺ → 2⁺ → 0⁺ cascade follow the
textbook correlation; and the simulated cascades keep the singles while correlating the γ rays."""

import math

import numpy as np
import pytest

from physim.nuclear import angular
from physim.nuclear.correlation import Emission, isotropic, rotation, unobserved
from physim.nuclear.gamma_events import _correlated
from physim.nuclear.levels import Level, LevelScheme, Transition, Value
from physim.nuclear.orientation import Excitation, _Decay, fed_density, tensors_of

CASES = [(2.0, 0.0, 2, 0, 0.0), (4.0, 2.0, 2, 0, 0.0), (2.0, 2.0, 2, 1, 0.5), (2.0, 2.0, 2, 1, -1.3),
         (1.5, 0.5, 1, 2, 0.3), (3.0, 2.0, 1, 2, 0.7), (0.0, 1.0, 1, 0, 0.0)]


def random_rho(j, rng):
    d = int(round(2 * j + 1))
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    r = a @ a.conj().T
    return r / np.trace(r)


def sphere(n_theta=40, n_phi=80):
    """Directions and weights of a quadrature over the sphere (weights sum to 4π)."""
    x, wx = np.polynomial.legendre.leggauss(n_theta)
    ph = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    xx, pp = np.meshgrid(x, ph, indexing="ij")
    s = np.sqrt(1 - xx**2)
    dirs = np.c_[(s * np.cos(pp)).ravel(), (s * np.sin(pp)).ravel(), xx.ravel()]
    return dirs, np.repeat(wx, n_phi) * (2 * np.pi / n_phi)


def w_from_tensors(rho, ji, jf, l1, l2, delta, dirs):
    """4π W from the statistical tensors and A_k, as Populated.w works it out."""
    t = tensors_of(rho, ji)
    th, ph = np.arccos(dirs[:, 2]), np.arctan2(dirs[:, 1], dirs[:, 0])
    total = np.zeros(len(dirs))
    for (k, q), v in t.items():
        if k % 2:
            continue
        a = angular.distribution_coefficient(k, ji, jf, l1, l2, delta) if k else 1.0
        total += (a * v * math.sqrt(4 * math.pi / (2 * k + 1)) * angular.spherical_harmonic(k, q, th, ph)).real
    return total


def test_rotation_takes_z_to_the_direction():
    # D^1 rotates the spherical components of a vector: the m = 0 state along n is n_z in the z basis.
    rng = np.random.default_rng(0)
    d = rng.normal(size=(50, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    r = rotation(1.0, d)
    assert np.allclose(r[:, 1, 1].real, d[:, 2])
    assert np.allclose(np.einsum("nab,nac->nbc", r.conj(), r), np.eye(3)[None])


@pytest.mark.parametrize("ji,jf,l1,l2,delta", CASES)
def test_one_gamma_ray_has_the_distribution_of_its_level(ji, jf, l1, l2, delta):
    rng = np.random.default_rng(1)
    rho = random_rho(ji, rng)
    dirs, _ = sphere(12, 24)
    e = Emission(ji, jf, l1, l2, delta)
    w = e.weights(np.repeat(rho[None], len(dirs), axis=0), dirs)
    assert np.allclose(w, w_from_tensors(rho, ji, jf, l1, l2, delta, dirs), atol=1e-12)
    assert w.max() <= e.bound + 1e-12


@pytest.mark.parametrize("ji,jf,l1,l2,delta", CASES)
def test_averaged_over_directions_the_next_level_is_fed_as_by_an_unseen_gamma_ray(ji, jf, l1, l2, delta):
    rng = np.random.default_rng(2)
    rho = random_rho(ji, rng)
    dirs, wq = sphere()
    e = Emission(ji, jf, l1, l2, delta)
    many = np.repeat(rho[None], len(dirs), axis=0)
    w = e.weights(many, dirs)
    after = e.after(many, dirs)
    assert np.allclose(np.einsum("naa->n", after).real, 1.0)
    average = np.einsum("n,n,nab->ab", wq, w, after) / (4 * math.pi)
    share = [(l1, 1 / (1 + delta**2))] + ([(l2, delta**2 / (1 + delta**2))] if l2 else [])
    fed = sum(s * fed_density(rho, ji, jf, lam) for lam, s in share)
    assert np.allclose(average, fed / np.trace(fed).real, atol=1e-12)
    assert np.allclose(unobserved(rho[None], ji, jf, l1, l2, delta)[0], fed / np.trace(fed).real, atol=1e-12)


def correlation_of(e1: Emission, e2: Emission, angles):
    """W(angle between the γ rays) of an unoriented level: the first along z, the second at each angle."""
    first = e1.after(isotropic(e1.ji)[None], np.array([[0.0, 0.0, 1.0]]))[0]
    d2 = np.c_[np.sin(angles), np.zeros(len(angles)), np.cos(angles)]
    return e2.weights(np.repeat(first[None], len(angles), axis=0), d2)


def legendre_fit(angles, w):
    c = np.cos(angles)
    basis = np.c_[np.ones_like(c), (3 * c**2 - 1) / 2, (35 * c**4 - 30 * c**2 + 3) / 8]
    k = np.linalg.lstsq(basis, w, rcond=None)[0]
    return k[1] / k[0], k[2] / k[0]


def test_an_unoriented_4_2_0_cascade_has_the_textbook_correlation():
    # ⁶⁰Co's cascade: W(θ) = 1 + 0.1020 P₂ + 0.0091 P₄.
    angles = np.linspace(0, math.pi, 13)
    a2, a4 = legendre_fit(angles, correlation_of(Emission(4, 2, 2), Emission(2, 0, 2), angles))
    assert a2 == pytest.approx(0.1020, abs=5e-5) and a4 == pytest.approx(0.0091, abs=5e-5)
    # 0 → 1 → 0 by two dipoles: 1 + cos²θ.
    a2, a4 = legendre_fit(angles, correlation_of(Emission(0, 1, 1), Emission(1, 0, 1), angles))
    assert a2 == pytest.approx(0.5, abs=1e-9) and abs(a4) < 1e-9


class _Cascade:
    """What :func:`_correlated` reads of a solution: the scheme's spins, the decays, and the density matrices
    the collision leaves (``start`` for every angle)."""

    def __init__(self, spins, decays, start):
        self.scheme = LevelScheme("x")
        self.scheme.levels = [Level(Value(100.0 * n, source="user"), s, 1) for n, s in enumerate(spins)]
        self._decays = decays
        self._start = start

    def decay(self, i, f):
        return self._decays[(i, f)]

    def density_table(self, level):
        grid = np.array([0.0, 180.0])
        return grid, np.repeat(self._start[level][None], 2, axis=0)


def _steps(n, chain):
    """``taken`` of :func:`_cascade_steps` for n excitations that all go down ``chain`` (levels, top first)."""
    src, ini, fin = [], [], []
    for a, b in zip(chain[:-1], chain[1:]):
        src.append(np.arange(n))
        ini.append(np.full(n, a))
        fin.append(np.full(n, b))
    src, ini, fin = np.concatenate(src), np.concatenate(ini), np.concatenate(fin)
    return (src, ini, fin, np.ones(len(src), dtype=bool)), np.full(n, chain[0])


def test_simulated_cascades_carry_the_correlation_and_keep_the_singles():
    n = 200_000
    e2 = _Decay(2, 1, 1.0, 1.0, 2, 0, 0.0)
    e1 = _Decay(1, 0, 1.0, 1.0, 2, 0, 0.0)
    cex = _Cascade([0.0, 2.0, 4.0], {(2, 1): e2, (1, 0): e1}, {2: isotropic(4.0)})
    taken, first = _steps(n, [2, 1, 0])
    known, deferred = _correlated(cex, np.full(n, 90.0), taken, first, np.random.default_rng(3), aimed=False)
    assert not deferred and not np.isnan(known).any()
    g1, g2 = known[0::2], known[1::2]  # each cascade's two γ rays, in the order they are emitted
    c = np.einsum("ij,ij->i", g1, g2)
    p2, p4 = np.mean((3 * c**2 - 1) / 2), np.mean((35 * c**4 - 30 * c**2 + 3) / 8)
    # ⟨P_k⟩ = a_k / (2k + 1) for W = 1 + a₂P₂ + a₄P₄.
    sigma = 1 / math.sqrt(5 * n)
    assert p2 == pytest.approx(0.1020 / 5, abs=4 * sigma) and p4 == pytest.approx(0.0091 / 9, abs=4 * sigma)
    # Each γ ray alone is isotropic, as the unoriented 4⁺ makes it.
    for g in (g1, g2):
        assert abs(np.mean((3 * g[:, 2] ** 2 - 1) / 2)) < 4 * sigma
    # When aimed, the last γ ray of each cascade is left to aim, with the state the first one left.
    known, deferred = _correlated(cex, np.full(1000, 90.0), *_steps(1000, [2, 1, 0]), np.random.default_rng(4),
                                  aimed=True)
    assert np.isnan(known[1::2]).all() and not np.isnan(known[0::2]).any()
    assert len(deferred) == 1 and np.array_equal(deferred[0][0], np.arange(1, 2000, 2))


def two_states() -> LevelScheme:
    """⁵⁸Ni-like: two 2⁺ states excited from the ground state; the upper decays to the lower by mixed E2/M1."""
    v = lambda x: Value(x, source="user")  # noqa: E731
    s = LevelScheme("58Ni")
    s.levels = [Level(v(0.0), 0.0, 1), Level(v(1454.2), 2.0, 1), Level(v(2775.4), 2.0, 1)]
    s.transitions = [Transition(1, 0, v(1454.2), v(100.0), "E2"),
                     Transition(2, 1, v(1321.2), v(100.0), "M1+E2", mixing=v(-0.6))]
    s.set_matrix_element(0, 1, "E2", 26.4)
    s.set_matrix_element(0, 2, "E2", 40.0)
    return s


def test_an_oriented_cascade_gives_each_line_the_distribution_of_its_level():
    ex = Excitation("16O", "58Ni", "30 MeV", two_states(), excite="target")
    theta = 150.0
    state = ex.at(theta)
    n = 200_000
    taken, first = _steps(n, [2, 1, 0])
    known, _ = _correlated(ex, np.full(n, theta), taken, first, np.random.default_rng(5), aimed=False)
    g1, g2 = known[0::2], known[1::2]
    dirs, wq = sphere()
    th, ph = np.arccos(dirs[:, 2]), np.arctan2(dirs[:, 1], dirs[:, 0])
    # The 2⁺₂ → 2⁺₁ line: the distribution of the directly excited 2⁺₂.
    # The 2⁺₁ → 0⁺ line of these cascades: the 2⁺₁ fed from the 2⁺₂ alone, by the unobserved transition.
    fed_only = ex._populated(theta, {2: state.density[2]})
    for g, pop, pair in ((g1, state, (2, 1)), (g2, fed_only, (1, 0))):
        w = pop.w(*pair, th, ph)
        w = w / np.sum(w * wq)
        for f in (lambda d: d[:, 2] ** 2, lambda d: d[:, 0] ** 2, lambda d: d[:, 0] * d[:, 1]):
            want = float(np.sum(wq * w * f(dirs)))
            got = f(g)
            assert got.mean() == pytest.approx(want, abs=4 * got.std() / math.sqrt(n)), pair
