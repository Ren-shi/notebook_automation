# Count rates and Monte Carlo events

`physim.nuclear.rates` works out how many counts each detector collects, how long to run, and where the peaks sit
and how wide they are. `physim.nuclear.events` simulates the experiment particle by particle, giving spectra. Both
use the same physics (items 35–38), so they check each other. Energies are in MeV, angles in degrees, rates per
second.

## Rates

A beam of $I$ particles per second on a layer with $n$ nuclei per cm² gives

$$
R = I\, n \int_\text{detector} \frac{d\sigma}{d\Omega}\, d\Omega
$$

counts per second. Physim integrates the lab cross section over each strip or ring–sector with Gauss–Legendre
quadrature (12 × 12 points per segment by default). It includes both kinematic solutions where the kinematics are
double-valued, and the ejectile (scattered beam) and the recoil (target nucleus) separately. The beam slows down
through the target and Rutherford's cross section grows as $1/E^2$, so the integral is averaged over depth with
4-point Gauss–Legendre quadrature, using the mean beam energy at each depth. Every nuclide in the target and in its
backing is a separate channel, so a carbon backing appears as its own rate and peak.

Two cuts keep the integral finite where Rutherford's formula diverges:

- **Lab angles below 0.5°** are left out, since the cross section diverges at 0°.
- **Particles leaving the reaction with less than the lowest detector threshold** are left out (at least 10 keV).
  Distant collisions send recoils out near 90° with almost no energy; they never leave the target or give a signal.

**Detectors hiding each other.** A particle stops in the first detector face it meets. The integral therefore
leaves out every direction that meets another detector first, the same rule the event generator follows track by
track.
- A detector entirely behind another counts nothing, and its warning names the one in front.
- At the edge of a shadow the integrand jumps, so a partly hidden detector is integrated on a finer grid, about 48 × 48
  points over its face. That is within 0.5% of the exact visible solid angle of a disc behind a smaller one.

**Beam time.** The time needed for $N$ counts is $N/R$. The relative statistical error is then $1/\sqrt N$.

**Counted rate.** `Rates.rate(..., counted=True)` scales each channel by the fraction of its particles whose mean
measured energy is above the threshold. Without smearing this is estimated for the whole detector. The event
generator decides threshold and punch-through particle by particle, so it is the exact reference.

**Warnings.** `Rates.warnings()` flags:

- detectors counting faster than 5000/s, where pile-up, dead time and forward-angle radiation damage begin;
- runs too short for the counts the setup asks for;
- the beam stopping in the target;
- Rutherford's validity limits (item 37) for every channel, over the CM angles the detectors see;
- the layout warnings of item 38.

## Expected peaks

`Rates.peaks(detector, segment=None)` predicts each peak's mean measured energy and its width. It samples 24 depths
through the layer and 16 × 16 directions over the face. At each point it:

1. slows the beam to that depth;
2. applies the two-body kinematics at that direction;
3. slows the particle out of the target along that direction (forward tracks leave by the back face, backward ones
   by the front, with paths $\Delta z/|\cos|$ capped near the target plane);
4. passes it through the dead layer and the active thickness (punch-through).

Each point is weighted by $d\sigma/d\Omega\, d\Omega$. The weighted mean is the peak position. The weighted spread
of these mean values is the **geometric** width, which combines **target thickness** (the spread over depth) and
**kinematic broadening** (the spread over the face); both are also reported alone. Three more widths are added in
quadrature:

- **beam energy spread**, carried to the measured energy;
- **straggling** of the beam on the way in and of the particle on the way out, through the dead layer and, for
  punch-through, the active layer;
- **resolution** (FWHM/2.355).

## Energy loss with straggling, from tables

The event generator carries particles through layers with tables of range $R(E)$, stopping power $S(E)$ and

$$
W(E) = \int^E \frac{1}{S(E')^3} \frac{d\Omega^2}{dx}(E')\, dE',
$$

where $d\Omega^2/dx$ is the Bohr straggling rate of {doc}`stopping` (effective charge, relativistic factor,
Lindhard–Scharff reduction). Crossing a path $L$:

$$
E_\text{out} = R^{-1}\big(R(E_\text{in}) - L\big), \qquad
\sigma^2_\text{added} = S(E_\text{out})^2\,\big[W(E_\text{in}) - W(E_\text{out})\big].
$$

This solves the straggling equation $d\sigma^2/dx = -2S'\sigma^2 + d\Omega^2/dx$ that `Stopping.straggling`
integrates step by step. Write $\sigma^2 = S^2 w$; then $dE/dx = -S$ gives $dw/dE = -(d\Omega^2/dx)/S^3$.

A spread already present is carried by the mapping $E_\text{in} \to E_\text{out}$ itself. Bohr straggling added at
energy $E$ reaches the exit multiplied by $S(E_\text{out})^2/S(E)^2$: it grows above the Bragg peak and shrinks
below it. Each crossing then adds a Gaussian of that width. The tables and `Stopping.straggling` agree within 2%
(tested).

## The event generator

Each event is one reaction. The per-event loop is written in Rust (`physim::nuclear`) and runs on all cores:

1. Pick a channel (nuclide and layer) with probability proportional to its rate.
2. Draw the beam energy (Gaussian, the setup's FWHM spread) and the spot position. The spot is Gaussian in x and y
   with the setup's FWHM, on the tilted target plane.
3. Draw the depth uniformly through the layer. Carry the beam to it with straggling, along paths $\Delta z/\cos t$
   for target tilt $t$.
4. Draw the CM angle $\theta^*$ through $u = \sin^2(\theta^*/2)$, over the range where the ejectile or the recoil
   can reach a detector with at least the cut energy. That range is found for the highest and lowest beam energies
   in the layer, widened by the spot size.
5. Do the two-body kinematics (the Rust port of item 35, equal to the Python version to $10^{-11}$) for the
   ejectile at $(\theta^*, \phi)$ and the recoil at $(180° - \theta^*, \phi + 180°)$.
6. Carry each particle out of the target and backing with straggling. Find the first detector face its straight
   track crosses, so shadowing is automatic. Apply the dead layer, active thickness, resolution and threshold.

**Weights.** Rutherford's $d\sigma/du \propto 1/u^2$ is steep, so drawing from it alone would leave backward
detectors almost empty. Half the draws come from $1/u^2$ and half from $1/u$ (flat in $\ln u$), with probability
density $p(u)$ for the mixture. Each event has the weight

$$
w = \frac{R_k(E)}{N\, P_k} \cdot \frac{f(u)}{p(u)},
$$

where $R_k(E)$ is the rate of channel $k$ into the sampled range at that event's beam energy, $P_k$ the
probability of picking channel $k$, $f$ the normalised Rutherford density and $N$ the number of events. The sum of
weights is an unbiased estimate of the rate, so spectra are absolute: counts per bin in the planned beam time. Its
statistical error is $\sqrt{\sum w^2}$.

**Randomness.** Every number is drawn from the counter-based generator of {doc}`integrators` with key
(seed, $2^{61}$ + event number, draw number). The output therefore depends only on the seed: it is bit-identical for
any thread count, and for a run split into pieces. Work is cut into fixed blocks of 4096 events, and the blocks are
joined in order.

**Speed.** $10^6$ events of the α + gold demo take 0.13 s on a 12-thread laptop (`cargo bench -- nuclear_events`).
`simulate(experiment, 1_000_000)` takes 0.2–0.4 s for the example setups, including building the tables in Python.

## Checks

- **Monte Carlo against analytic rates.** Every detector agrees within 4 statistical standard deviations, per
  channel and in total, with $2 \times 10^6$ events. The setups are the examples, a tilted-target variant, and a pad
  half hidden behind a CD.
- **Hidden detectors.** The visible solid angle of a disc behind a smaller one is checked against 2π(cos α_front −
  cos α_back) within 0.5%. A pad entirely behind the CD counts nothing.
- **Peaks.** Mean energies agree with kinematics plus mean energy loss within statistics. Widths agree with the
  quadrature sum within 5%. The α + gold peak at 30° matches a hand calculation within 1 keV.
- **A small detector's rate.** It equals $I\,n\,(d\sigma/d\Omega)\,\Omega$ within 1%.
- **Reproducibility.** Same seed, same events; 1 thread and 4 threads give bit-identical events; a split run equals
  a whole one.

## Limits

- Only elastic (Rutherford) scattering so far. Above the barrier the rates are Rutherford's, and the warnings say
  where that fails.
- Particles travel in straight lines: multiple scattering in the target and dead layer is not sampled (Highland
  widths are in {doc}`stopping`), and the beam has no angular divergence.
- Straggling is Gaussian. That is good for thin layers but not near the end of a particle's range, where the
  distribution is skewed.
- Detector faces have no inter-strip effects, and silicon pulse-height defects are not modelled (see
  {doc}`detectors`). Identical beam and target nuclei are not given Mott interference.
- No background, random coincidences or beam halo.
