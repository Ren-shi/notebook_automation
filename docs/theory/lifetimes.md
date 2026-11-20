# Lifetimes: plunger (RDDS) and Doppler-shift attenuation (DSAM)

`physim.nuclear.lifetimes` holds the physics of lifetime measurements as pure functions on numbers: B(Eλ) ↔ τ, decay
curves of a cascade, the shifted and unshifted shares in a plunger, the DSAM attenuation and line shape, and the
deorientation of a recoil flying in vacuum; `physim.nuclear.lifetime_analysis` gets τ back from spectra; the
planner puts them together ([Lifetime measurements](../lifetimes.md)). Units: time ps, distance µm, γ-ray energy
keV, speed β = v/c.

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

## Planning a measurement

**Method.** A plunger measures τ of about 1 ps to 1 ns (decay lengths γβcτ of tens of µm to tens of mm at a few %
of c); DSAM about 10 fs to 1 ps (comparable to the time a recoil takes to stop in a backing).

**Beam time per distance.** For one level, the unshifted share at flight time t is R = e^{−x}, x = t/τ. With N
counts of the line at that distance, R has the binomial error √(R(1 − R)/N), and τ = t/(−ln R), so

$$
\frac{\sigma_\tau}{\tau} = \frac{p(x)}{\sqrt N},\qquad p(x) = \frac{\sqrt{e^x - 1}}{x},
$$

smallest at x ≈ 1.59 (where `recommend_distances` puts its best single distance). Measuring N counts at each of the
distances x_i and combining, (σ_τ/τ)⁻² = N Σ_i p(x_i)⁻², so the counts per distance for a wanted precision ε are
N = 1/(ε² Σ_i p(x_i)⁻²), and the beam time per distance is N over the line's full-energy rate in the γ-ray
detectors. This is statistics only and one level: feeding, background and the systematic entries of the analysis's
budget (distance, velocity, stopping) come on top.

**Resolved peaks.** The planner calls the shifted and unshifted peaks resolved at a detector when their separation
E₀√(1 − β²)/(1 − β cos θ) − E₀ is at least twice the detector's FWHM: two Gaussians of one width are apart at about
1.7 FWHM, and the shifted peak is broader (the spread of speeds and angles over the crystal).

## Analysis

**Peak areas.** Fusion residues fly at nearly one speed, and the shifted and unshifted peaks are fitted as Gaussians
on a straight background. After Coulomb excitation the emitters' speeds spread widely (a particle detector at 0°
sees both the target recoil and forward-scattered beam, so the excited projectile's speed ranges over almost a factor
of two): one Gaussian catches only part of the shifted component, which biases the unshifted share upwards at every
distance. There the planner sums the counts instead: the unshifted side from E₀ − 3σ to E₀ + half the shift (it
holds the decays while the recoil slows in the stopper, which `lifetimes.unshifted` counts as unshifted), the
shifted side from there to twice the planned shift, both over a straight background through the counts beyond. The
recoil's speed for the flight times comes from the shifted peaks' centroid, E_s/E₀ − 1 ≈ β cos θ (the outermost
peak, for a differential plunger).

**DDCM, joint fit.** For a trial τ, the curve P(t) (linear in its coefficients) minimises
Σ(S − P)²/σ_S² + Σ(N − τP′)²/σ_N². Here N = I_u − α I_u^B (singles) or {B_s, A_u} (gated), and S = I_s or
{B_s, A_s}. χ²(τ) is profiled, and the error comes from Δχ² = 1, scaled by χ²/ndf when that is above 1. τ(x) = N/P′
at each point, from the best curve, for the plot.

**α.** With efficiency-corrected intensities normalised by the transition's own I_u + I_s, α is the share of the
feeder's decays that go through the transition (α = b ε_A/ε_B for raw intensities). Gated data are not normalised by
I_u + I_s, because that sum follows the feeder's decay: they take an explicit normalisation.

**Three foils.** R(t), the share of decays after t, is one function of time. The degrader point sits at t₁ = x₁/v
and the stopper point at t₂ = x₁/v + x₂/v′, that is at x_eff = x₁ + x₂ v/v′ when flown at v. This is derived from
first principles; the comparison with Beckers et al.'s Eqs. 13–14 is one of the validation checks.

**Stopping budget (DSAM).** A uniform scaling of the stopping power by k scales every slowing time by 1/k, so it is
a scaling of the kernel's time axis. The fit refits with ±7% by default (heavy-ion stopping at these speeds is known
to about 5–10%); for one level without feeding this moves τ by exactly that share.

**Deorientation correction.** W(θ) = 1 + Σ a_k G_k P_k(cos θ). The shifted component uses G̃_k(T) (above), the
unshifted one G_k(T). With f = W_U/W_S, the corrected unshifted share is R/(R + f(1 − R)). The aberration of the
shifted γ rays' angle is left out.
