# Neutron detectors and channel selection

`physim.nuclear.neutrons` describes the neutron detectors of a fusion-evaporation experiment and how well the
detected particles tell one exit channel from another (backlog item 75, "Detectors for channel selection").
Energies are in MeV, light in MeVee (the energy an electron would need to give the same light), lengths in mm,
times in ns.

## The detector

A cell is a cylinder of organic scintillator of radius $R$ and depth $D$. Its axis points at the target and its
front face is at distance $d$. Three scintillators are in `data/neutron_scintillators.toml`, with hydrogen and
carbon densities from the Eljen data sheets:

| Material | Kind | Density (g/cm³) | H (atoms/cm³) | C (atoms/cm³) | Same as |
|---|---|---|---|---|---|
| EJ-301 | liquid, PSD | 0.874 | 4.82 × 10²² | 3.98 × 10²² | NE213, BC-501A |
| EJ-309 | liquid, PSD | 0.959 | 5.43 × 10²² | 4.35 × 10²² | |
| EJ-200 | plastic | 1.023 | 5.17 × 10²² | 4.69 × 10²² | NE102, BC-408 |

`NeutronDetector()` with no arguments is a typical EJ-301 cell of a neutron wall: 127 mm across, 150 mm deep, 510 mm
from the target, 100 keVee threshold, 1.5 ns FWHM timing, 2% of neutrons called γ rays and 1% of γ rays called
neutrons, 40% of γ rays above threshold. All of these except the material are typical values, listed in
`TYPICAL_FIELDS`, to be replaced by those of the detector used.

## Intrinsic efficiency

A neutron is seen when it scatters on a proton and the recoil proton makes more light than the threshold. For a
neutron of energy $E$ entering the face along the axis,

$$
\varepsilon(E) = \frac{\Sigma_H}{\Sigma}\left(1 - e^{-\Sigma D}\right)\left(1 - \frac{E_{p,\text{th}}}{E}\right),
\qquad \Sigma = n_H\sigma_H(E) + n_C\sigma_C(E),
$$

the closed form of Knoll (*Radiation Detection and Measurement*, 4th ed., 2010, chapter 15):

- $\frac{\Sigma_H}{\Sigma}(1 - e^{-\Sigma D})$ is the probability that the first collision in the cell is on
  hydrogen;
- below about 20 MeV, n–p scattering is isotropic in the centre of mass, so the recoil energy is uniform from 0 to
  $E$, and the share above the proton threshold is $1 - E_{p,\text{th}}/E$.

**Cross sections.** $\sigma_H$ is Gammel's closed form (in *Fast Neutron Physics*, part II, ed. Marion and Fowler,
1963), good to about 1% up to 40 MeV:

$$
\sigma_H(E) = \frac{3\pi}{1.206E + (-1.86 + 0.09415E + 0.0001306E^2)^2} + \frac{\pi}{1.206E + (0.4223 + 0.13E)^2}
\ \text{b}.
$$

It gives 4.26 b at 1 MeV and 0.689 b at 14.1 MeV, the neutron-standard values, and agrees with the n–p total
cross section of ENDF/B-VIII.0 (H-1, MT=1) within 0.8% from 0.5 to 20 MeV.

$\sigma_C$ is the total cross section of natural carbon from **ENDF/B-VIII.0** (MF=3, MT=1 of ¹²C and ¹³C, weighted
98.93% and 1.07%), interpolated in log E from the data file. The file holds 40 values per decade from 10 keV to
30 MeV, each the average of the pointwise cross section over its own bin (a factor 1.029 either side), so a
resonance narrower than a bin keeps its area but is spread over the bin. `scripts/endf_neutron_reference.py`
writes both the table and the comparison points `tests/reference/nuclear/neutron_carbon_endf.csv` from the ENDF-6
files, and `test_carbon_table_is_the_endf_evaluation` holds the table to the evaluation:

| Energy range | physim vs ENDF/B-VIII.0 (averaged over ±2.5% / ±10% of E) |
|---|---|
| 0.5–1.8 MeV (smooth) | within 0.1% / 0.1% |
| 1.8–11 MeV (resonances) | within 6.4% (4.1 MeV) / 4.1% (2.3 MeV, the narrow 2.08 MeV peak spread over its bins) |
| 11–20 MeV | within 1.3% / 0.3% |

Before 15 November 2026 the table held rounded "typical" values that were 31% low at 3.4 MeV, 18% low at 2.8 MeV
and 24% high at 8.9 MeV. Carbon is half of the attenuation near 3.4 MeV, so the fix lowered the intrinsic efficiency
of the typical 150 mm EJ-301 cell by 12% at 3.4 MeV and 6% at 2.8 MeV, and raised it by 6% at 8.9 MeV; below 2 MeV
and at 5 and 14 MeV the change is under 1%, and cross-talk changed by less than 3%.

**Light output.** The threshold is set in light. The proton energy that gives that light comes from the
light-output curve of Cecil, Anderson and Madey (Nucl. Instrum. Methods 161, 439 (1979)):

$$
L(E_p) = a_1E_p - a_2\left(1 - e^{-a_3E_p^{a_4}}\right),
$$

with $(a_1, a_2, a_3, a_4) = (0.83, 2.82, 0.25, 0.93)$ for NE213 and $(0.95, 8.0, 0.1, 0.90)$ for NE102. EJ-309
uses the NE213 values and EJ-200 the NE102 values, both marked typical. 100 keVee is a proton of 0.63 MeV in NE213.

The efficiency for neutrons from the target is $\varepsilon(E)\,\Omega/4\pi$, with $\Omega = 2\pi(1 - d/\sqrt{d^2
+ R^2})$ for the front face.

## Time of flight

A neutron of kinetic energy $E$ has $\gamma = 1 + E/m_nc^2$ and $\beta = \sqrt{1 - 1/\gamma^2}$; over the flight
path $L = d + D/2$ it takes $t = L/\beta c$. A γ ray takes $L/c$ (the γ flash). The energy measured from the time
has

$$
\frac{\sigma_E}{E} = \gamma(\gamma + 1)\sqrt{\left(\frac{\sigma_t}{t}\right)^2 + \left(\frac{\sigma_L}{L}\right)^2},
\qquad \sigma_L = \frac{D}{\sqrt{12}},
$$

where $\sigma_t$ is the time resolution and $\sigma_L$ the unknown depth of the interaction. Without relativity
$\gamma(\gamma + 1) \to 2$.

## Neutron–γ discrimination

Pulse-shape discrimination is described by two probabilities: a neutron called a γ ray (it is lost), and a γ ray
called a neutron (it adds a false neutron). Fusion-evaporation events carry many γ rays, so the second matters even
when it is small. For a measured figure of merit $\text{FOM} = \Delta / (\text{FWHM}_n + \text{FWHM}_\gamma)$ with
equal peak widths and the cut halfway, each probability is $\tfrac12\operatorname{erfc}(2.3548\,\text{FOM}/\sqrt2)$
(`misidentification_from_fom`).

## Cross-talk

A neutron that scatters in one cell and then in a neighbour fires two cells and looks like two neutrons. The
probability $c$ that a neighbour fires, given that the first cell fired, is found by a seeded Monte Carlo:

1. The first collision is at a point uniform on the face and at a depth drawn from $e^{-\Sigma z}$ within the cell.
2. The neutron scatters on a proton at lab angle $\theta$, with $\cos^2\theta$ uniform (isotropic CM). It keeps
   $E\cos^2\theta$ and gives the proton $E\sin^2\theta$, which must be above $E_{p,\text{th}}$.
3. It leaves the first cell with probability $e^{-\Sigma'\ell_1}$ and fires the neighbour with
   $\frac{\Sigma'_H}{\Sigma'}(1 - e^{-\Sigma'\ell_2})(1 - E_{p,\text{th}}/E')$, where $\ell_1$ and $\ell_2$ are the
   chords to the exit and through the neighbour.

Both the proton and the scattered neutron must be above the threshold, so $c = 0$ for $E \le 2E_{p,\text{th}}$.
`NeutronWall` sums $c$ over the six neighbours of a hexagonal wall (or eight of a square one) and can remove an
average share for cells on the edge. For the typical wall (cells touching at a pitch of 130 mm), $c$ is about 2% at
2 MeV, 6% at 5 MeV and 8% at 10 MeV.

## The channel-selection matrix

Each true channel $i$ emits $n_i$ neutrons, $p_i$ protons, $a_i$ α particles and $g_i$ γ rays ($g_i$ fixed, or a
distribution $P_i(g)$ when it varies from event to event). Each particle is detected independently, with one of a
few outcomes:

| Particle | Outcomes and probabilities |
|---|---|
| neutron | 0: $1 - q$; 1: $q(1 - c)$; 2: $qc$, with $q = \varepsilon_n(1 - P_{n\to\gamma})$ |
| γ ray | 1 false neutron: $\varepsilon_\gamma P_{\gamma\to n}$ |
| proton | none: $1 - \varepsilon_p$; proton: $\varepsilon_p(1 - P_{p\to\alpha})$; α: $\varepsilon_pP_{p\to\alpha}$ |
| α | none: $1 - \varepsilon_\alpha$; proton: $\varepsilon_\alpha P_{\alpha\to p}$; α: $\varepsilon_\alpha(1 - P_{\alpha\to p})$ |

The detected numbers of neutrons, protons and α particles are sums over the particles, so their distribution is the
convolution of the per-particle distributions (binomial and multinomial). An event is assigned to the channel whose
$(n, p, a)$ equals what was detected, or left unassigned when no channel matches. The matrix $M_{ij}$ is the share
of true channel $i$ assigned to channel $j$; each row, with its unassigned share, sums to one. With cross sections
$\sigma_i$, the purity of assigned channel $j$ is $\sigma_jM_{jj} / \sum_i\sigma_iM_{ij}$.

With $f = \varepsilon_\gamma P_{\gamma\to n}$ the chance that one γ ray poses as a neutron, the false neutrons of an
event with $g$ γ rays are binomial, $B(g, f)$; when $g$ varies, their distribution is the mixture
$\sum_g P_i(g)\,B(g, f)$, which is convolved with the true neutrons' counts (`gamma_distributions` gives each channel
its $P_i(g)$ from event lists).

**Monte Carlo.** `simulate_selection` draws the same outcomes particle by particle from event lists: per event,
each particle's kind, lab energy and lab direction. The efficiencies can then be functions of each particle's energy
and direction (`NeutronWall.detection_probability`). `Detection.averaged` gives the mean outcome probabilities per
channel over the same events; with these, the analytic matrix is what the Monte Carlo converges to, as long as the
particles of one event are detected independently. The tests check the agreement within statistics (item 75's
"Detectors" done-when), also with thresholds and a γ-ray number that varies per event. The validation suite holds
this as `channel_selection_binomial` (every cell within 4.5 binomial standard errors, 6 × 10⁴ events; worst 1.9)
and the efficiency as `neutron_efficiency_knoll` (Knoll's closed form, 1e-12).

## Channel selection in the planner

`exit_channels.channel_selection` sorts the residue events of a fusion-evaporation setup with a `[neutron_wall]`
into detected classes ("γ", "n", "p2n", ...), by Monte Carlo over the events, and also returns the analytic matrix of
the same events' mean probabilities (the two agree within statistics; the tests hold them to 4.5σ). Per particle:

- **Neutrons**: the wall's efficiency at the neutron's energy inside the wall's cone, cross-talk at the mean
  neutron energy, and the cells' neutron–γ misidentification.
- **Protons and α particles: particle-detector thresholds.** A particle counts when it meets a particle detector
  (dead material in front stops it) and deposits more than that detector's threshold: its own `threshold`, or
  `[neutron_wall] particle_threshold` for detectors that set none (default 1 MeV, typical of silicon and CsI arrays
  for evaporated particles). The deposit is computed with the detector's response: the dead layer first, then the
  thickness, with path lengths $1/\cos$ of the incidence angle, so a particle that punches through deposits only
  part of its energy.
- **γ rays firing neutron cells.** Each event has $M_\gamma$ γ rays: `[neutron_wall] gamma_multiplicity` if set,
  otherwise $J/2 + 4$ from the event's entry spin $J$ (the stretched E2 cascade plus about four statistical γ rays,
  the usual conversion $J \approx 2(M_\gamma - 4)$), rounded up or down at random to keep the mean. Each γ ray fires
  a cell with the wall's geometric efficiency times the cells' `gamma_efficiency` (typical 40%) and is counted as a
  neutron with `gamma_as_neutron` (typical 1%), which is taken to include the time-of-flight cut.

**What they change**, for the example ¹²C + ¹⁵⁰Nd at 60 MeV with 50 typical cells at 510 mm and two 300 µm silicon
annular detectors (forward and backward, 32% of the protons and α particles reach them), 50 000 events:

- Evaporated protons from ¹⁶²Dy (10–17 MeV in the lab, 5th to 95th percentile) punch through 300 µm of silicon
  and deposit 1.8–3 MeV near normal incidence, more at the steeper angles near the edges of the annuli. The default
  1 MeV threshold keeps all of them (no change); a 3 MeV threshold keeps 38% of them, and the share of p3n events
  seen as "p" falls from 26% to 9% (as "pn" from 5.2% to 2.5%). The α particles (18–28 MeV) deposit 17–25 MeV, so
  no threshold below that changes them.
- With $J/2 + 4$ γ rays (9.4 per event here) and 7.7% of γ rays firing a cell, 1% of those called neutrons, a γ ray
  is a false neutron with 0.077%, so about 0.7% of events gain a false neutron: in each row about 0.5–0.8
  percentage points move from "γ" (nothing seen) to "n" (for 4n, 70.3% → 69.8% and 24.7% → 25.0%). The effect grows
  in proportion to $M_\gamma$ and to `gamma_as_neutron`: a wall with 5% γ-as-neutron and 20 γ rays loses 7% of
  its no-neutron events to false neutrons.

## Limits

- One collision only: double n–p scattering, which raises the efficiency at low thresholds, is left out. The
  efficiency is for neutrons along the axis; edge effects of a finite radius are left out.
- Carbon only attenuates. Above 5 MeV, 20–43% of the carbon total is non-elastic in ENDF/B-VIII.0
  (`neutron_carbon_endf.csv`): ¹²C(n,α) and ¹²C(n,n′3α) above about 7 MeV add light that is left out, and so is
  cross-talk through scattering on carbon. A neutron scattered elastically on carbon keeps at least 72% of its
  energy and could still be seen; counting it as lost underestimates the efficiency somewhat.
- Isotropic CM n–p scattering: above about 20 MeV the recoil spectrum is no longer flat.
- Cross-talk assumes identical parallel cells (fine when the pitch is much smaller than the distance) and no
  material between them; γ-ray cross-talk and γ rays scattered into neutron cells are left out.
- In the matrix, particles are detected independently. In a real array two particles can hit the same cell, and
  time-of-flight or neighbour-rejection cuts against cross-talk are not modelled.
- In the planner's selection: protons and α particles lose no energy in the target on their way out, and the
  particle detectors tell protons from α particles perfectly (a ΔE-only silicon cannot, for punch-through
  particles). γ rays are emitted isotropically, their number comes from the entry spin by the rule above (not from
  the residue's level scheme), and a γ ray that fires a cell is either called a γ ray or one false neutron
  (γ-ray cross-talk is left out).
- Neutron–γ discrimination is two fixed probabilities, not a function of the light.
