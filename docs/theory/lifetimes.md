---
orphan: true
---

# Lifetimes: plunger (RDDS) and Doppler-shift attenuation (DSAM)

`physim.nuclear.lifetimes` holds the physics of lifetime measurements as pure functions on numbers: B(Eλ) ↔ τ, decay
curves of a cascade, the shifted and unshifted shares in a plunger, the DSAM attenuation and line shape, and the
deorientation of a recoil flying in vacuum. Units: time ps, distance µm, γ-ray energy keV, speed β = v/c. (This note
covers wave 1 of the lifetime work; the planner, analysis and validation parts come later.)

## B(λ) and the mean life

A transition's γ-ray rate per unit B(λ↓) is

$$
T(E\lambda) = \frac{8\pi(\lambda+1)}{\lambda[(2\lambda+1)!!]^2}\,\frac{e^2}{\hbar}
\left(\frac{E_\gamma}{\hbar c}\right)^{2\lambda+1} B(E\lambda\downarrow),
$$

the same for Mλ times (µ_N/e)² = (ħc/2m_pc²)², with B(Mλ) in µ_N² fm^{2λ−2}. A level decays through all its
transitions; one transition carries the share b = I_γ(1 + α)/Σ_k I_γ,k(1 + α_k) of its decays (γ rays and conversion
electrons), and of its γ rays the multipole λ carries the share s (1/(1 + δ²) for M1, δ²/(1 + δ²) for E2, with δ
the E2/M1 mixing ratio). So

$$
\frac{1}{\tau} = \frac{T(\lambda)\,(1+\alpha)}{s\,b},
$$

which is the inverse of `LevelScheme.derive`. `level_taus` uses a level's half-life when it has one (τ = T½/ln 2),
otherwise the B(Eλ↓) of its transitions.

## Decay curves of a cascade

With decay constants λ_i = 1/τ_i and a feeding matrix f (row = from, column = to), the numbers in the levels follow
dN/dt = A N with A = fᵀΛ − Λ (Bateman's equations). The closed form of the solution subtracts nearly equal numbers
when two lifetimes are close; physim computes exp(At) by uniformisation instead: with c = max λ_i, A + c·1 has no
negative entry, so exp(Ah) = e^{−ch} Σ_k ((A + c)h)^k/k! and the squarings that take h back up to t add only
non-negative numbers. Every entry keeps its relative accuracy, far down the tail and for τ₁ = τ₂.

In a plunger the γ rays emitted after the flight time t to the stopper form the unshifted peak. For each level, its
share is

$$
R_i(t) = \frac{[(1 - f^T)^{-1} N(t)]_i}{[(1 - f^T)^{-1} N(0)]_i},
$$

the nuclei in the level or still to pass through it, over all that ever do; for one level, R = e^{−t/τ}.

**Flight time.** t_lab = d/(βc cos θ) for a recoil at angle θ to the axis; the decay law runs in the nucleus's own
time, t_lab/γ. `flight_times` gives lab times (with `proper=True`, proper times); the decay functions take times in
the frame of the τ they are given, so the caller decides whether to include the dilation (below 0.2% for β < 6%).

**Distances.** From an unshifted share R = e^{−x} (x = t/τ) measured with N counts, R has the binomial error
√(R(1 − R)/N) and τ = t/(−ln R), so σ_τ/τ = √(e^x − 1)/(x√N): best at x = 1.59, and within a factor 2 of that
between x = 0.18 and 5.1. `sensitive_range` is that interval in distance (d = x γβcτ) and `recommend_distances`
spaces distances logarithmically over it. The ⁹⁸Ru (3–450 µm, vτ = 113 µm) and ⁵⁰Cr (17–750 µm, vτ = 143 µm)
plunger series cover it.

**Peak separation.** The shifted peak lies at E₀√(1 − β²)/(1 − β cos θ).

## DSAM

The kernel is β·cos θ of a sample of recoils (θ from the beam axis) at a grid of times after the state is formed,
followed through the target and backing (`transport.kernel`). It does not depend on τ. The attenuation is

$$
F(\tau) = \frac{\langle \beta\cos\theta \rangle_\text{emission}}{\langle \beta\cos\theta \rangle_{t=0}}
= \int_0^\infty \frac{e^{-t/\tau}}{\tau}\,\frac{\langle\beta\cos\theta\rangle(t)}{\langle\beta\cos\theta\rangle(0)}\,dt,
$$

with ⟨β cos θ⟩ taken as linear between the kernel's times and integrated exactly against the decay law (for a cascade,
against the level's own decay curve). For stopping proportional to the speed, v = v₀e^{−t/α}, F = α/(α + τ) exactly.
The line shape bins each recoil's E/E₀ = √(1 − (β cos)²)/(1 − β cos θ · cos θ_γ) between two times, weighted by the
decays in that interval. This binning is the same for every τ, so it is done once and a fit over τ only re-weights
it (the method of AhKin, Garzón et al. 2020). The recoil's sideways motion is left out, which matters only for
detectors near 90°.

## Deorientation in vacuum

A recoil flying in vacuum loses its alignment through the hyperfine field of its electrons. In the static model
(Radeck et al., PRC 85, 014301 (2012)) the k-th order of the angular distribution is attenuated by
G_k(t) = α_k + (1 − α_k)e^{−Γ_k t}. The unshifted γ rays see G_k at the flight time T (the alignment is frozen once
the nucleus stops); the shifted ones see the average over the decays in flight,

$$
\tilde G_k(T) = \frac{\int_0^T G_k(t)\,e^{-t/\tau}dt}{\int_0^T e^{-t/\tau}dt}
= \alpha_k + (1-\alpha_k)\,\frac{(1 - e^{-(\Gamma_k + 1/\tau)T})/(\Gamma_k + 1/\tau)}{\tau\,(1 - e^{-T/\tau})},
$$

which is 1 at T = 0.
