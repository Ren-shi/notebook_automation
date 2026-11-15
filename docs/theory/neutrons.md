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

It gives 4.26 b at 1 MeV and 0.689 b at 14.1 MeV, the neutron-standard values. $\sigma_C$ is the total cross
section of natural carbon, interpolated in log E from the data file. **The carbon values are typical:** smooth
averages over the resonances, rounded, following the ENDF/B-VIII.0 evaluation to about 15%, not read point by
point. Because carbon only attenuates, they change the efficiency of a 150 mm cell by a few percent at most.

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
2 MeV, 7% at 5 MeV and 8% at 10 MeV.

## The channel-selection matrix

Each true channel $i$ emits $n_i$ neutrons, $p_i$ protons, $a_i$ α particles and $g_i$ γ rays. Each particle is
detected independently, with one of a few outcomes:

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

The charged-particle efficiencies are given as numbers, since the silicon and CsI detectors of the setup provide
them.

**Monte Carlo.** `simulate_selection` draws the same outcomes particle by particle from event lists: per event,
each particle's kind, lab energy and lab direction. The efficiencies can then be functions of each particle's energy
and direction (`NeutronWall.detection_probability`). `Detection.averaged` gives the mean outcome probabilities per
channel over the same events; with these, the analytic matrix is what the Monte Carlo converges to, as long as the
particles of one event are detected independently. The tests check the agreement within statistics (item 75's
"Detectors" done-when). The validation suite holds this as `channel_selection_binomial` (every cell within 4.5
binomial standard errors, 6 × 10⁴ events; worst 1.9) and the efficiency as `neutron_efficiency_knoll` (Knoll's
closed form, 1e-12).

## Limits

- One collision only: double n–p scattering, which raises the efficiency at low thresholds, is left out. The
  efficiency is for neutrons along the axis; edge effects of a finite radius are left out.
- Carbon only attenuates. ¹²C(n,α) and ¹²C(n,n′3α) above about 7 MeV add light that is left out, and so is
  cross-talk through scattering on carbon. The carbon cross sections are typical (above).
- Isotropic CM n–p scattering: above about 20 MeV the recoil spectrum is no longer flat.
- Cross-talk assumes identical parallel cells (fine when the pitch is much smaller than the distance) and no
  material between them; γ-ray cross-talk and γ rays scattered into neutron cells are left out.
- In the matrix, particles are detected independently. In a real array two particles can hit the same cell, and
  time-of-flight or neighbour-rejection cuts against cross-talk are not modelled.
- Neutron–γ discrimination is two fixed probabilities, not a function of the light.
