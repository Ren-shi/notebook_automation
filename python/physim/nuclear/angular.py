"""Angular-momentum algebra for γ-ray angular distributions: Wigner 3j and 6j symbols, Clebsch–Gordan
coefficients, spherical harmonics, and the F and U coefficients of angular-correlation theory.

Spins may be integers or half-integers (given as numbers such as ``2`` or ``1.5``). The spherical harmonics use the
Condon–Shortley phase. The F coefficients follow Krane and Steffen, whose sign convention for mixing ratios is the
one ENSDF uses::

    from physim.nuclear import angular

    angular.clebsch_gordan(2, 0, 2, 0, 2, 0)      # ⟨2 0 2 0 | 2 0⟩ = −√(2/7)
    angular.f_coefficient(2, 2, 2, 0, 2)          # F₂(2 2 0 2) = −0.5976
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np


def _int(x: float, what: str = "value") -> int:
    r = round(x)
    if abs(x - r) > 1e-9:
        raise ValueError(f"{what} {x} is not an integer")
    return int(r)


def _fact(x: float) -> float:
    return float(math.factorial(_int(x)))


def triangle(a: float, b: float, c: float) -> bool:
    """Whether three angular momenta can couple: |a − b| ≤ c ≤ a + b, with a + b + c an integer."""
    return abs(a - b) - 1e-9 <= c <= a + b + 1e-9 and abs((a + b + c) - round(a + b + c)) < 1e-9


@lru_cache(maxsize=None)
def wigner_3j(j1: float, j2: float, j3: float, m1: float, m2: float, m3: float) -> float:
    """The Wigner 3j symbol (j1 j2 j3; m1 m2 m3), by Racah's formula."""
    if abs(m1 + m2 + m3) > 1e-9 or not triangle(j1, j2, j3):
        return 0.0
    if any(abs(m) > j + 1e-9 for j, m in ((j1, m1), (j2, m2), (j3, m3))):
        return 0.0
    if any(abs((j - m) - round(j - m)) > 1e-9 for j, m in ((j1, m1), (j2, m2), (j3, m3))):
        return 0.0
    delta = (_fact(j1 + j2 - j3) * _fact(j1 - j2 + j3) * _fact(-j1 + j2 + j3) / _fact(j1 + j2 + j3 + 1))
    norm = math.sqrt(delta * _fact(j1 + m1) * _fact(j1 - m1) * _fact(j2 + m2) * _fact(j2 - m2)
                     * _fact(j3 + m3) * _fact(j3 - m3))
    lo = max(0, _int(j2 - j3 - m1), _int(j1 - j3 + m2))
    hi = min(_int(j1 + j2 - j3), _int(j1 - m1), _int(j2 + m2))
    total = 0.0
    for t in range(lo, hi + 1):
        total += (-1) ** t / (_fact(t) * _fact(j3 - j2 + t + m1) * _fact(j3 - j1 + t - m2)
                              * _fact(j1 + j2 - j3 - t) * _fact(j1 - t - m1) * _fact(j2 - t + m2))
    return (-1) ** _int(j1 - j2 - m3) * norm * total


def clebsch_gordan(j1: float, m1: float, j2: float, m2: float, j: float, m: float) -> float:
    """The Clebsch–Gordan coefficient ⟨j1 m1 j2 m2 | j m⟩."""
    if abs(m1 + m2 - m) > 1e-9:
        return 0.0
    return (-1) ** _int(j1 - j2 + m) * math.sqrt(2 * j + 1) * wigner_3j(j1, j2, j, m1, m2, -m)


def _delta(a: float, b: float, c: float) -> float:
    return _fact(a + b - c) * _fact(a - b + c) * _fact(-a + b + c) / _fact(a + b + c + 1)


@lru_cache(maxsize=None)
def wigner_6j(j1: float, j2: float, j3: float, j4: float, j5: float, j6: float) -> float:
    """The Wigner 6j symbol {j1 j2 j3; j4 j5 j6}, by Racah's formula."""
    triads = ((j1, j2, j3), (j1, j5, j6), (j4, j2, j6), (j4, j5, j3))
    if not all(triangle(*t) for t in triads):
        return 0.0
    norm = math.sqrt(math.prod(_delta(*t) for t in triads))
    a = [_int(sum(t)) for t in triads]
    b = [_int(j1 + j2 + j4 + j5), _int(j2 + j3 + j5 + j6), _int(j3 + j1 + j6 + j4)]
    total = 0.0
    for t in range(max(a), min(b) + 1):
        total += ((-1) ** t * math.factorial(t + 1)
                  / (math.prod(math.factorial(t - x) for x in a) * math.prod(math.factorial(x - t) for x in b)))
    return norm * total


@lru_cache(maxsize=None)
def f_coefficient(k: int, l1: int, l2: int, j_final: float, j_initial: float) -> float:
    """The F coefficient F_k(L L′ J_f J_i) of a γ ray of multipolarities L and L′ from J_i to J_f (Krane and
    Steffen)."""
    return ((-1) ** _int(j_final + j_initial - 1)
            * math.sqrt((2 * k + 1) * (2 * l1 + 1) * (2 * l2 + 1) * (2 * j_initial + 1))
            * wigner_3j(l1, l2, k, 1, -1, 0) * wigner_6j(l1, l2, k, j_initial, j_initial, j_final))


@lru_cache(maxsize=None)
def distribution_coefficient(k: int, j_initial: float, j_final: float, l1: int, l2: int = 0,
                             mixing: float = 0.0) -> float:
    """A_k of a γ ray from J_i to J_f: F_k for a pure multipole L, and for a mixture with mixing ratio δ (the
    amplitude of L′ = ``l2`` over that of L = ``l1``)

        A_k = [F_k(L L) + 2δ F_k(L L′) + δ² F_k(L′ L′)] / (1 + δ²)."""
    pure = f_coefficient(k, l1, l1, j_final, j_initial)
    if not l2 or not mixing:
        return pure
    return (pure + 2 * mixing * f_coefficient(k, l1, l2, j_final, j_initial)
            + mixing**2 * f_coefficient(k, l2, l2, j_final, j_initial)) / (1 + mixing**2)


def u_coefficient(k: int, j_initial: float, j_final: float, l: int) -> float:  # noqa: E741
    """U_k(J_i J_f L): how much of the orientation of rank k a level keeps when it is fed by an unobserved
    transition of multipolarity L from J_i to J_f. U_0 = 1."""
    return ((-1) ** _int(j_initial + j_final + l + k) * math.sqrt((2 * j_initial + 1) * (2 * j_final + 1))
            * wigner_6j(j_initial, j_initial, k, j_final, j_final, l))


def spherical_harmonic(l: int, m: int, theta, phi) -> np.ndarray:  # noqa: E741
    """Y_lm(θ, φ) with the Condon–Shortley phase; angles in radians."""
    theta, phi = np.asarray(theta, dtype=float), np.asarray(phi, dtype=float)
    am = abs(m)
    x = np.cos(theta)
    # Associated Legendre P_l^|m|(x), by the standard upward recursion from P_m^m.
    pmm = (-1) ** am * _double_factorial(2 * am - 1) * np.sin(theta) ** am
    if l == am:
        plm = pmm
    else:
        pm1 = x * (2 * am + 1) * pmm
        if l == am + 1:
            plm = pm1
        else:
            a, b = pmm, pm1
            for n in range(am + 2, l + 1):
                a, b = b, ((2 * n - 1) * x * b - (n + am - 1) * a) / (n - am)
            plm = b
    norm = math.sqrt((2 * l + 1) / (4 * math.pi) * math.factorial(l - am) / math.factorial(l + am))
    y = norm * plm * np.exp(1j * am * phi)
    return y if m >= 0 else (-1) ** am * np.conj(y)


def _double_factorial(n: int) -> float:
    out = 1.0
    while n > 1:
        out *= n
        n -= 2
    return out


def multipole_pattern(amplitudes, theta, phi) -> np.ndarray:
    """The angular distribution of radiation from a source with multipole amplitudes a_M of one order L,
    |Σ_M a_M X_LM(θ, φ)|² with the vector spherical harmonics X_LM = L Y_LM / √(L(L+1)) (Jackson, section 9.8).

    ``amplitudes`` holds a_M for M = −L … L. This is the distribution of the γ rays of a state Σ a_M |L M⟩
    decaying to a spin-0 state, and integrates to Σ |a_M|². It is the direct formula against which the tensor
    formalism is checked."""
    a = np.asarray(amplitudes, dtype=complex)
    L = (len(a) - 1) // 2  # noqa: N806
    theta, phi = np.asarray(theta, dtype=float), np.asarray(phi, dtype=float)
    lz = np.zeros(np.broadcast(theta, phi).shape, dtype=complex)
    lp, lm = lz.copy(), lz.copy()
    for i, m in enumerate(range(-L, L + 1)):
        if a[i] == 0:
            continue
        lz = lz + a[i] * m * spherical_harmonic(L, m, theta, phi)
        if m < L:
            lp = lp + a[i] * math.sqrt((L - m) * (L + m + 1)) * spherical_harmonic(L, m + 1, theta, phi)
        if m > -L:
            lm = lm + a[i] * math.sqrt((L + m) * (L - m + 1)) * spherical_harmonic(L, m - 1, theta, phi)
    # |L_x f|² + |L_y f|² = (|L_+ f|² + |L_- f|²) / 2
    return (np.abs(lz) ** 2 + (np.abs(lp) ** 2 + np.abs(lm) ** 2) / 2) / (L * (L + 1))


__all__ = ["clebsch_gordan", "distribution_coefficient", "f_coefficient", "multipole_pattern",
           "spherical_harmonic", "triangle", "u_coefficient", "wigner_3j", "wigner_6j"]
