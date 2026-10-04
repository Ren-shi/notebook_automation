//! The Monte Carlo event generator.
//!
//! Each event is one reaction of a beam particle in the target (or its backing), followed by
//! the ejectile and the recoil to the detectors:
//!
//! 1. pick a reaction channel (a nuclide in a layer) in proportion to its rate;
//! 2. draw the beam energy (energy spread), the beam spot position and the depth of the
//!    reaction, uniformly through the layer;
//! 3. slow the beam to that depth, with straggling;
//! 4. draw the CM scattering angle, restricted to the range of angles the detectors can see,
//!    half the time from the Rutherford cross section and half the time flat in ln sin²(θ*/2),
//!    and weight the event so rates stay absolute (for Coulomb excitation, also by the excitation
//!    probability P(θ*));
//! 5. two-body kinematics for the ejectile and the recoil;
//! 6. slow each outgoing particle through the rest of the target and backing, find the first
//!    detector face it crosses, and apply that detector's response (dead layer, punch-through,
//!    resolution, threshold).
//!
//! Every random number is drawn from [`crate::rng`] with the key (seed, [`EVENTS`] + event
//! number, draw number), so the output depends only on the seed, never on the thread count.

use std::f64::consts::{PI, TAU};

use super::detector::{Face, Hit};
use super::kinematics::TwoBody;
use super::stopping::Table;
use crate::rng;
use crate::vec3::Vec3;

/// Counter offset of the event generator's draws: event `i` uses counter `EVENTS + i`.
pub const EVENTS: u64 = 1 << 61;
/// Events per block of parallel work (fixed, so results do not depend on the thread count).
const BLOCK: u64 = 4096;
/// Z₁Z₂e² is given in MeV fm, so σ comes out in fm²; 1 fm² = 1e-26 cm².
const FM2_TO_CM2: f64 = 1e-26;
/// FWHM = 2√(2 ln 2) σ.
const FWHM_PER_SIGMA: f64 = 2.354_820_045_030_949;

/// A layer of the target stack, in beam order, perpendicular to the target normal.
#[derive(Clone, Debug)]
pub struct Layer {
    /// Thickness along the normal, mg/cm².
    pub thickness: f64,
    /// Index into the material list.
    pub material: usize,
}

/// Rutherford scattering on one nuclide of one layer.
#[derive(Clone, Debug)]
pub struct Channel {
    pub layer: usize,
    /// Species index of the target nuclide (which is also the recoil).
    pub target: usize,
    /// Nuclei of this kind per cm² in the layer.
    pub atoms_per_cm2: f64,
    /// Z₁Z₂e², MeV fm.
    pub k: f64,
    /// Sampled range of u = sin²(θ*/2), θ* the CM angle of the ejectile.
    pub u_min: f64,
    pub u_max: f64,
    /// Relative probability of picking this channel (normalised internally).
    pub probability: f64,
    /// Excitation energy left in the recoil (`excite_recoil`) or in the ejectile, MeV; 0 for
    /// elastic scattering.
    pub excitation: f64,
    pub excite_recoil: bool,
    /// Excitation probability P(θ*) on equally spaced CM angles from 0° to 180° (Coulomb
    /// excitation); empty for elastic scattering. Event weights are multiplied by it.
    pub p_table: Vec<f64>,
}

impl Channel {
    /// The excitation probability at CM angle `theta_cm` (radians), 1 for elastic scattering.
    fn excitation_probability(&self, theta_cm: f64) -> f64 {
        let n = self.p_table.len();
        if n < 2 {
            return 1.0;
        }
        let x = (theta_cm / PI).clamp(0.0, 1.0) * (n - 1) as f64;
        let i = (x as usize).min(n - 2);
        let f = x - i as f64;
        self.p_table[i] + f * (self.p_table[i + 1] - self.p_table[i])
    }
}

/// Everything the generator needs. Tables are indexed `species × materials + material`.
#[derive(Clone, Debug)]
pub struct Generator {
    /// Nuclear masses of the species, MeV; species `beam` is the beam (and the ejectile).
    pub masses: Vec<f64>,
    pub beam: usize,
    /// Mean beam energy (MeV) and its standard deviation.
    pub beam_energy: f64,
    pub energy_sigma: f64,
    /// Standard deviation of the beam spot in x and in y, mm.
    pub spot_sigma: f64,
    pub particles_per_second: f64,
    /// Target tilt about the vertical (y) axis, radians.
    pub tilt: f64,
    pub layers: Vec<Layer>,
    pub channels: Vec<Channel>,
    pub materials: usize,
    pub tables: Vec<Table>,
    pub faces: Vec<Face>,
    /// Longest path through a layer, as a multiple of the whole stack's thickness (for tracks
    /// almost in the target plane).
    pub max_path_factor: f64,
}

/// One particle that reached a detector face.
#[derive(Clone, Copy, Debug, Default)]
pub struct Record {
    pub event: u64,
    pub detector: u32,
    pub segment: (u32, u32),
    pub recoil: bool,
    pub channel: u32,
    /// Depth of the reaction below the front face, along the normal, mg/cm².
    pub depth: f64,
    /// Beam energy at the reaction, MeV.
    pub beam_energy: f64,
    /// Energy of this particle just after the reaction, at the detector face, deposited in
    /// the active volume, and measured (with resolution), MeV.
    pub energy: f64,
    pub energy_face: f64,
    pub deposited: f64,
    pub measured: f64,
    /// Reached the face with energy left and was measured above threshold.
    pub counted: bool,
    /// Lab angles of the track and the CM angle of the ejectile, radians.
    pub theta: f64,
    pub phi: f64,
    pub theta_cm: f64,
    /// Rate this event stands for, per second.
    pub weight: f64,
}

/// The draws of one event: uniforms number 0, 1, 2, ... of its counter.
struct Draw {
    seed: u64,
    counter: u64,
    next: u64,
}

impl Draw {
    fn uniform(&mut self) -> f64 {
        self.next += 1;
        rng::uniform(self.seed, self.counter, self.next - 1)
    }

    fn normal(&mut self) -> f64 {
        let (u1, u2) = (self.uniform(), self.uniform());
        (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
    }
}

impl Generator {
    /// Checks indices and sizes; [`Generator::run`] calls it first.
    pub fn validate(&self) -> Result<(), String> {
        let ns = self.masses.len();
        if self.beam >= ns || self.tables.len() != ns * self.materials {
            return Err("tables must cover every species in every material".into());
        }
        if self.layers.iter().any(|l| {
            l.material >= self.materials
                || l.thickness.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
        }) {
            return Err("layers need a known material and a positive thickness".into());
        }
        if self.channels.is_empty() {
            return Err("at least one reaction channel is needed".into());
        }
        for c in &self.channels {
            if c.layer >= self.layers.len() || c.target >= ns {
                return Err("a channel refers to an unknown layer or species".into());
            }
            if c.excitation < 0.0 || c.p_table.iter().any(|p| !(p.is_finite() && *p >= 0.0)) {
                return Err(
                    "channel excitation energies and probabilities must not be negative".into(),
                );
            }
            if !(0.0 < c.u_min && c.u_min < c.u_max && c.u_max <= 1.0)
                || c.probability.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater)
            {
                return Err(
                    "channel angle ranges need 0 < u_min < u_max <= 1 and a positive probability"
                        .into(),
                );
            }
        }
        if self.faces.iter().any(|f| f.material >= self.materials) {
            return Err("a detector refers to an unknown material".into());
        }
        Ok(())
    }

    fn table(&self, species: usize, material: usize) -> &Table {
        &self.tables[species * self.materials + material]
    }

    /// Energy after crossing `path` (mg/cm²) of `material`, with straggling.
    fn cross(&self, d: &mut Draw, species: usize, material: usize, e: f64, path: f64) -> f64 {
        if e <= 0.0 || path <= 0.0 {
            return e.max(0.0);
        }
        let (out, sigma) = self.table(species, material).cross(e, path);
        if out <= 0.0 {
            return 0.0;
        }
        (out + sigma * d.normal()).clamp(0.0, e)
    }

    /// Rate of channel `c` (per second) into its sampled angle range at CM energy `e_cm` (MeV):
    /// beam particles/s × nuclei/cm² × 4π (Z₁Z₂e²/4E)² (1/u_min − 1/u_max).
    pub fn rate(&self, c: &Channel, e_cm: f64) -> f64 {
        let a = c.k / (4.0 * e_cm);
        let sigma_fm2 = 4.0 * PI * a * a * (1.0 / c.u_min - 1.0 / c.u_max);
        self.particles_per_second * c.atoms_per_cm2 * sigma_fm2 * FM2_TO_CM2
    }

    /// Generates events `first .. first + n` of a run of `total` events (which sets the
    /// weights), with `seed`. Records come in event order, ejectile before recoil.
    pub fn run(&self, n: u64, total: u64, seed: u64, first: u64) -> Result<Vec<Record>, String> {
        self.validate()?;
        let blocks: Vec<(usize, usize)> = (0..n.div_ceil(BLOCK))
            .map(|b| {
                let lo = first + b * BLOCK;
                (lo as usize, (first + ((b + 1) * BLOCK).min(n)) as usize)
            })
            .collect();
        let norm: f64 = self.channels.iter().map(|c| c.probability).sum();
        let cumulative: Vec<f64> = self
            .channels
            .iter()
            .scan(0.0, |acc, c| {
                *acc += c.probability / norm;
                Some(*acc)
            })
            .collect();
        let parts = crate::parallel::map_blocks(&blocks, |lo, hi| {
            let mut out = Vec::new();
            for i in lo..hi {
                self.event(i as u64, total.max(1), seed, &cumulative, norm, &mut out);
            }
            out
        });
        Ok(parts.into_iter().flatten().collect())
    }

    fn event(
        &self,
        i: u64,
        total: u64,
        seed: u64,
        cumulative: &[f64],
        norm: f64,
        out: &mut Vec<Record>,
    ) {
        let mut d = Draw {
            seed,
            counter: EVENTS.wrapping_add(i),
            next: 0,
        };
        let x = d.uniform();
        let k = cumulative
            .partition_point(|&c| c < x)
            .min(cumulative.len() - 1);
        let ch = &self.channels[k];
        let mut e = self.beam_energy + self.energy_sigma * d.normal();
        let (sx, sy) = (self.spot_sigma * d.normal(), self.spot_sigma * d.normal());
        let (st, ct) = self.tilt.sin_cos();
        // The spot lies on the tilted target plane n · p = 0, with n = (sin t, 0, cos t).
        let source = Vec3::new(sx, sy, -sx * st / ct);
        let z_in = d.uniform() * self.layers[ch.layer].thickness;
        let mut depth = 0.0;
        for (j, layer) in self.layers.iter().enumerate().take(ch.layer + 1) {
            let dz = if j == ch.layer { z_in } else { layer.thickness };
            e = self.cross(&mut d, self.beam, layer.material, e, dz / ct);
            depth += dz;
        }
        let (m1, m2) = (self.masses[self.beam], self.masses[ch.target]);
        let (m3, m4) = if ch.excite_recoil {
            (m1, m2 + ch.excitation)
        } else {
            (m1 + ch.excitation, m2)
        };
        let Some(tb) = TwoBody::new(m1, m2, m3, m4, e) else {
            return; // the beam stopped before this depth, or is below the excitation threshold
        };
        let (u0, u1) = (ch.u_min, ch.u_max);
        // Draw u = sin²(θ*/2) half the time from Rutherford's 1/u² (most events at forward angles,
        // where most of the rate is) and half the time from 1/u (flat in ln u, so backward
        // detectors get events too); the weight f/p corrects for the mixture.
        let (a, b) = (1.0 / u0 - 1.0 / u1, (u1 / u0).ln());
        let u = if d.uniform() < 0.5 {
            1.0 / (1.0 / u0 - d.uniform() * a)
        } else {
            u0 * (d.uniform() * b).exp()
        }
        .clamp(u0, u1);
        let f = 1.0 / (u * u * a);
        let p = 0.5 * f + 0.5 / (u * b);
        let theta_cm = 2.0 * u.sqrt().min(1.0).asin();
        let phi = TAU * d.uniform();
        let weight = self.rate(ch, tb.cm_energy()) * f / p * ch.excitation_probability(theta_cm)
            / (total as f64 * ch.probability / norm);
        let normal = Vec3::new(st, 0.0, ct);
        let stack: f64 = self.layers.iter().map(|l| l.thickness).sum();
        let cap = self.max_path_factor * stack;
        for (recoil, th_star, ph, species) in [
            (false, theta_cm, phi, self.beam),
            (true, PI - theta_cm, phi + PI, ch.target),
        ] {
            let lab = tb.at_cm(th_star, recoil);
            let (s, c) = lab.theta.sin_cos();
            let dir = Vec3::new(s * ph.cos(), s * ph.sin(), c);
            // Out through the rest of the stack: forward tracks leave by the back face,
            // backward ones by the front face. Paths grow as 1/cos, up to the cap.
            let cos_n = dir.dot(normal);
            let path = |dz: f64| {
                if cos_n != 0.0 {
                    (dz / cos_n.abs()).min(cap)
                } else {
                    cap
                }
            };
            let mut ep = lab.energy;
            if cos_n > 0.0 {
                for (j, layer) in self.layers.iter().enumerate().skip(ch.layer) {
                    let dz = if j == ch.layer {
                        layer.thickness - z_in
                    } else {
                        layer.thickness
                    };
                    ep = self.cross(&mut d, species, layer.material, ep, path(dz));
                }
            } else {
                for j in (0..=ch.layer).rev() {
                    let layer = &self.layers[j];
                    let dz = if j == ch.layer { z_in } else { layer.thickness };
                    ep = self.cross(&mut d, species, layer.material, ep, path(dz));
                }
            }
            let mut best: Option<(usize, Hit)> = None;
            for (f, face) in self.faces.iter().enumerate() {
                if let Some(h) = face.hit(source, dir) {
                    if best.is_none_or(|(_, b)| h.distance < b.distance) {
                        best = Some((f, h));
                    }
                }
            }
            let Some((f, hit)) = best else { continue };
            let face = &self.faces[f];
            let (mut deposited, mut measured) = (0.0, 0.0);
            if ep > 0.0 {
                let ci = hit.cos_incidence;
                let e_active = self.cross(&mut d, species, face.material, ep, face.dead_layer / ci);
                let e_back = self.cross(
                    &mut d,
                    species,
                    face.material,
                    e_active,
                    face.thickness / ci,
                );
                deposited = e_active - e_back;
                measured = deposited + face.fwhm / FWHM_PER_SIGMA * d.normal();
            }
            out.push(Record {
                event: i,
                detector: f as u32,
                segment: hit.segment,
                recoil,
                channel: k as u32,
                depth,
                beam_energy: e,
                energy: lab.energy,
                energy_face: ep,
                deposited,
                measured,
                counted: ep > 0.0 && measured > 0.0 && measured >= face.threshold,
                theta: lab.theta,
                phi: ph.rem_euclid(TAU),
                theta_cm,
                weight,
            });
        }
    }
}

impl Generator {
    /// A small self-contained setup for tests and benchmarks: 5.5 MeV α particles on 0.5 mg/cm²
    /// of gold, with power-law stopping powers (S = c/√E) and constant straggling, seen by three
    /// silicon discs (30°, 60°, 135°) and an annular strip detector upstream of the target.
    pub fn demo() -> Generator {
        let power_law = |c: f64, straggle: f64| {
            let e: Vec<f64> = (0..=600)
                .map(|k| 1e-4 * 10f64.powf(k as f64 / 100.0))
                .collect();
            let r: Vec<f64> = e.iter().map(|e| 2.0 / 3.0 * e.powf(1.5) / c).collect();
            let s: Vec<f64> = e.iter().map(|e| c / e.sqrt()).collect();
            let w: Vec<f64> = e
                .iter()
                .map(|e| straggle * e.powf(2.5) / (2.5 * c.powi(3)))
                .collect();
            Table::new(&e, &r, &s, &w).unwrap()
        };
        // Species: α, ¹⁹⁷Au. Materials: gold, silicon.
        let tables = vec![
            power_law(0.23, 2e-4),
            power_law(1.3, 1e-4),
            power_law(4.0, 1e-4),
            power_law(9.0, 1e-4),
        ];
        let disc = |theta: f64, material| {
            let t = f64::to_radians(theta);
            let centre = Vec3::new(100.0 * t.sin(), 0.0, 100.0 * t.cos());
            let n = centre * (-0.01);
            Face {
                centre,
                n,
                u: Vec3::new(0.0, 1.0, 0.0).cross(n),
                v: Vec3::new(0.0, 1.0, 0.0),
                shape: super::Shape::Annular {
                    inner: 0.0,
                    outer: 5.0,
                    rings: 1,
                    sectors: 1,
                },
                material,
                dead_layer: 0.1,
                thickness: 70.0,
                fwhm: 0.02,
                threshold: 0.2,
            }
        };
        let mut cd = disc(180.0, 1);
        cd.centre = Vec3::new(0.0, 0.0, -40.0);
        cd.n = Vec3::new(0.0, 0.0, 1.0);
        cd.u = Vec3::new(1.0, 0.0, 0.0);
        cd.v = Vec3::new(0.0, 1.0, 0.0);
        cd.shape = super::Shape::Annular {
            inner: 9.0,
            outer: 41.0,
            rings: 16,
            sectors: 24,
        };
        let (u_min, u_max) = (
            5f64.to_radians().sin().powi(2),
            87.5f64.to_radians().sin().powi(2),
        );
        Generator {
            masses: vec![3727.379, 183473.2],
            beam: 0,
            beam_energy: 5.5,
            energy_sigma: 0.005,
            spot_sigma: 0.5,
            particles_per_second: 6.24e9,
            tilt: 0.0,
            layers: vec![Layer {
                thickness: 0.5,
                material: 0,
            }],
            channels: vec![Channel {
                layer: 0,
                target: 1,
                atoms_per_cm2: 1.529e18,
                k: 2.0 * 79.0 * 1.439_964_48,
                u_min,
                u_max,
                probability: 1.0,
                excitation: 0.0,
                excite_recoil: true,
                p_table: Vec::new(),
            }],
            materials: 2,
            tables,
            faces: vec![disc(30.0, 1), disc(60.0, 1), disc(135.0, 1), cd],
            max_path_factor: 1e3,
        }
    }
}
