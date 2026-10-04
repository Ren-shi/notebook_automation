"""Helpers for :class:`physim.Schrodinger`: momentum-space densities, probability currents,
probabilities in a region, and analytic results to check simulations against.

::

    import physim as ps

    s = ps.Schrodinger([2048], 0.05, potential=lambda t, x: 0.5 * (abs(x) < 0.5))
    s.set_gaussian(center=-20, width=4, momentum=1.2)
    s.step(0.05, 600)
    ps.quantum.probability(s, lambda x: x > 0.5)        # transmitted probability
    ps.quantum.barrier_transmission(0.72, 0.5, 1.0)     # plane-wave value at E = p²/2m
"""

from __future__ import annotations

import numpy as np


def _coords(s):
    return np.meshgrid(*s.axes, indexing="ij")


def probability(s, region):
    """Probability ``∫_region |ψ|² dV``: ``region`` is a boolean array of the grid's shape or a
    function of the coordinate arrays (``region(x[, y[, z]])``) returning one."""
    mask = region(*_coords(s)) if callable(region) else region
    mask = np.broadcast_to(np.asarray(mask, dtype=bool), tuple(s.shape))
    return float(np.sum(s.density[mask]) * s.cell_volume)


def momentum_density(s):
    """The momentum distribution ``|φ(p)|²``, normalised so that ``Σ |φ|² Δp^d = Σ |ψ|² h^d``.

    Returns ``(p_axes, density)``: one sorted momentum axis per grid axis (``p = ħk``) and the
    density on that grid. The grid is treated as periodic (any boundary is fine as long as ψ
    vanishes at the edges)."""
    psi = np.asarray(s.psi)
    h = s.spacing
    phi = np.fft.fftshift(np.fft.fftn(psi))
    axes = [s.hbar * 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(n, d=h)) for n in s.shape]
    dp = [a[1] - a[0] for a in axes]
    density = np.abs(phi) ** 2
    density *= np.sum(np.abs(psi) ** 2) * s.cell_volume / (np.sum(density) * np.prod(dp))
    return axes, density


def probability_current(s):
    """The probability current ``j = (ħ/m) Im(ψ* ∇ψ)``, an array ``(d, *shape)``. Gradients are
    spectral on periodic grids (matching the split-step solver) and central differences
    otherwise."""
    psi = np.asarray(s.psi)
    periodic = s.boundary == "periodic"
    out = []
    for axis in range(psi.ndim):
        if periodic:
            n = psi.shape[axis]
            k = 2 * np.pi * np.fft.fftfreq(n, d=s.spacing)
            if n % 2 == 0:
                k[n // 2] = 0.0  # the Nyquist mode has no definite derivative
            shape = [1] * psi.ndim
            shape[axis] = n
            grad = np.fft.ifft(1j * k.reshape(shape) * np.fft.fft(psi, axis=axis), axis=axis)
        else:
            grad = np.gradient(psi, s.spacing, axis=axis)
        out.append(s.hbar / s.mass * np.imag(np.conj(psi) * grad))
    return np.array(out)


def barrier_transmission(energy, height, width, mass=1.0, hbar=1.0):
    """Transmission probability of a plane wave of ``energy`` through a rectangular barrier of
    ``height`` and ``width`` (vectorised over ``energy``)::

        T = 1 / (1 + V₀² sinh²(κa) / (4E(V₀ - E)))     E < V₀,  κ = √(2m(V₀ - E))/ħ
        T = 1 / (1 + V₀² sin²(qa)  / (4E(E - V₀)))     E > V₀,  q = √(2m(E - V₀))/ħ

    A negative ``height`` is a square well (resonant transmission)."""
    e = np.asarray(energy, dtype=float)
    if np.any(e <= 0):
        raise ValueError("energy must be positive")
    v0, a = float(height), float(width)
    diff = e - v0
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        q = np.sqrt(2 * mass * np.abs(diff)) / hbar
        s2 = np.where(diff < 0, np.sinh(q * a) ** 2, np.sin(q * a) ** 2)
        t = 1.0 / (1.0 + v0**2 * s2 / (4 * e * np.abs(diff)))
    # E = V₀ exactly: the limit of both branches.
    t = np.where(diff == 0, 1.0 / (1.0 + mass * a * a * v0 / (2 * hbar**2)), t)
    return t if t.ndim else float(t)


def packet_transmission(momentum, width, height, barrier_width, mass=1.0, hbar=1.0, samples=4001):
    """Transmission of a Gaussian packet with mean ``momentum`` and position spread ``width``
    (σ), the plane-wave result averaged over its momentum distribution
    ``∝ exp(-2σ²(p - p₀)²/ħ²)``. Compare with ``probability(s, region beyond the barrier)``
    once the packet has separated."""
    sigma_p = hbar / (2 * width)
    p = momentum + sigma_p * np.linspace(-8, 8, samples)
    p = p[p > 0]
    w = np.exp(-((p - momentum) ** 2) / (2 * sigma_p**2))
    t = barrier_transmission(p**2 / (2 * mass), height, barrier_width, mass, hbar)
    return float(np.sum(w * t) / np.sum(w))


def harmonic_levels(n, omega=1.0, hbar=1.0, dims=1):
    """The lowest ``n`` energies ``ħω(n₁ + … + n_d + d/2)`` of an isotropic oscillator, with
    degeneracies repeated."""
    from itertools import product

    top = n  # enough quanta per axis to cover the lowest n levels
    levels = sorted(sum(q) for q in product(range(top), repeat=dims))[:n]
    return hbar * omega * (np.array(levels, dtype=float) + dims / 2)


def box_levels(n, length, mass=1.0, hbar=1.0):
    """The lowest ``n`` energies ``k²π²ħ²/(2mL²)`` of a 1D infinite square well of ``length``."""
    k = np.arange(1, n + 1)
    return (k * np.pi * hbar / length) ** 2 / (2 * mass)
