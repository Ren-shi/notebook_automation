//! Hard collisions: instantaneous impulses between particles with a radius, and against
//! fixed planar walls, resolved event by event inside each step.
//!
//! After every integrator step the world looks for the earliest contact along the step
//! (swept spheres on a uniform grid, with contact times from the straight path between the
//! start and end positions, refined on the integrator's own trajectory when forces act). It
//! re-steps the world to that moment, applies the collision impulse with coefficient of
//! restitution `e` along the line of centres (or the wall normal), and continues the rest of
//! the step, until no contact remains. With no forces between collisions this is exact
//! event-driven dynamics. Pinned particles act as infinitely massive.

use crate::broadphase::candidate_pairs;
use crate::error::{invalid, Result};
use crate::state::State;
use crate::vec3::Vec3;
use crate::world::World;

/// A fixed plane `normal · x = offset`. Particles stay on the side the unit `normal` points
/// to, i.e. `normal · x ≥ offset + r`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wall {
    pub normal: Vec3,
    pub offset: f64,
}

impl Wall {
    /// The six walls of the box `[lo, hi]`.
    pub fn box_walls(lo: Vec3, hi: Vec3) -> Vec<Wall> {
        let mut walls = Vec::new();
        for axis in 0..3 {
            let mut n = [0.0; 3];
            n[axis] = 1.0;
            walls.push(Wall {
                normal: Vec3::from(n),
                offset: lo[axis],
            });
            n[axis] = -1.0;
            walls.push(Wall {
                normal: Vec3::from(n),
                offset: -hi[axis],
            });
        }
        walls
    }
}

/// Settings for hard collisions (see the module documentation).
#[derive(Debug, Clone, PartialEq)]
pub struct Collisions {
    /// Coefficient of restitution: 1 elastic, 0 perfectly inelastic.
    pub restitution: f64,
    pub walls: Vec<Wall>,
    /// Collide particles with each other (otherwise only with the walls).
    pub between_particles: bool,
    /// Fail a step that needs more collisions than this (guards against a clamped pile).
    pub max_per_step: usize,
}

impl Default for Collisions {
    fn default() -> Self {
        Self {
            restitution: 1.0,
            walls: Vec::new(),
            between_particles: true,
            max_per_step: 100_000,
        }
    }
}

impl Collisions {
    pub(crate) fn validate(&self) -> Result<()> {
        if !(0.0..=1.0).contains(&self.restitution) {
            return invalid(format!(
                "restitution must be in [0, 1], got {}",
                self.restitution
            ));
        }
        for w in &self.walls {
            let n = w.normal.norm();
            if !((n - 1.0).abs() < 1e-9 && w.offset.is_finite()) {
                return invalid(format!(
                    "wall normals must be unit vectors and offsets finite, got {:?}, {}",
                    w.normal, w.offset
                ));
            }
        }
        Ok(())
    }
}

/// What collides in a contact.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Contact {
    Pair(usize, usize),
    Wall(usize, usize),
}

/// Contact time τ ∈ [0, 1] along straight paths from `start` to `end`, if any.
fn pair_time(s0: &State, s1: &State, i: usize, j: usize) -> Option<f64> {
    let r = s0.radius[i] + s0.radius[j];
    let d0 = s0.pos[j] - s0.pos[i];
    let delta = (s1.pos[j] - s1.pos[i]) - d0;
    let c = d0.norm_squared() - r * r;
    if c <= 0.0 {
        // Touching or overlapping at the start: a contact now only if approaching.
        return ((s0.vel[j] - s0.vel[i]).dot(d0) < 0.0).then_some(0.0);
    }
    let a = delta.norm_squared();
    let b = 2.0 * d0.dot(delta);
    if a == 0.0 || b >= 0.0 {
        return None;
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return None;
    }
    let tau = (-b - disc.sqrt()) / (2.0 * a);
    (0.0..=1.0).contains(&tau).then_some(tau)
}

fn wall_gap(s: &State, wall: &Wall, i: usize) -> f64 {
    wall.normal.dot(s.pos[i]) - wall.offset - s.radius[i]
}

fn wall_time(s0: &State, s1: &State, wall: &Wall, i: usize) -> Option<f64> {
    let (g0, g1) = (wall_gap(s0, wall, i), wall_gap(s1, wall, i));
    if g0 <= 0.0 {
        return (wall.normal.dot(s0.vel[i]) < 0.0).then_some(0.0);
    }
    (g1 < 0.0).then(|| g0 / (g0 - g1))
}

/// Signed gap of a contact (negative when overlapping).
fn gap(s: &State, walls: &[Wall], c: Contact) -> f64 {
    match c {
        Contact::Pair(i, j) => (s.pos[j] - s.pos[i]).norm() - s.radius[i] - s.radius[j],
        Contact::Wall(w, i) => wall_gap(s, &walls[w], i),
    }
}

/// Inverse mass: zero for pinned particles (immovable).
fn inverse_mass(s: &State, i: usize) -> Result<f64> {
    if s.pinned[i] {
        return Ok(0.0);
    }
    if s.mass[i] == 0.0 {
        return invalid(format!("particle {i} collides but is massless"));
    }
    Ok(1.0 / s.mass[i])
}

/// Applies the collision impulse; returns whether anything changed.
fn apply(s: &mut State, walls: &[Wall], e: f64, c: Contact) -> Result<bool> {
    match c {
        Contact::Pair(i, j) => {
            let d = s.pos[j] - s.pos[i];
            let n = d / d.norm();
            let approach = (s.vel[j] - s.vel[i]).dot(n);
            let (wi, wj) = (inverse_mass(s, i)?, inverse_mass(s, j)?);
            if approach >= 0.0 || wi + wj == 0.0 {
                return Ok(false);
            }
            let impulse = -(1.0 + e) * approach / (wi + wj);
            s.vel[i] -= n * (impulse * wi);
            s.vel[j] += n * (impulse * wj);
        }
        Contact::Wall(w, i) => {
            let n = walls[w].normal;
            let vn = n.dot(s.vel[i]);
            if vn >= 0.0 || s.pinned[i] {
                return Ok(false);
            }
            s.vel[i] -= n * ((1.0 + e) * vn);
        }
    }
    Ok(true)
}

impl World {
    /// Enables hard collisions with these settings (`None` disables them).
    pub fn set_collisions(&mut self, collisions: Option<Collisions>) -> Result<()> {
        if let Some(c) = &collisions {
            c.validate()?;
        }
        self.collisions = collisions;
        Ok(())
    }

    pub fn collisions(&self) -> Option<&Collisions> {
        self.collisions.as_ref()
    }

    /// Collisions resolved since the world was created.
    pub fn collision_count(&self) -> u64 {
        self.collision_count
    }

    /// Earliest contact between `start` and `end` (step fraction and contact), if any.
    fn earliest_contact(
        &self,
        c: &Collisions,
        start: &State,
        end: &State,
    ) -> Option<(f64, Contact)> {
        let mut best: Option<(f64, Contact)> = None;
        let mut consider = |tau: f64, contact: Contact| {
            if best.is_none_or(|(t, _)| tau < t) {
                best = Some((tau, contact));
            }
        };
        if c.between_particles {
            // Swept spheres: radius plus displacement over the step.
            let reach: Vec<f64> = (0..start.len())
                .map(|i| {
                    if start.radius[i] > 0.0 {
                        start.radius[i] + (end.pos[i] - start.pos[i]).norm()
                    } else {
                        0.0
                    }
                })
                .collect();
            for (i, j) in candidate_pairs(&start.pos, &reach) {
                if let Some(tau) = pair_time(start, end, i, j) {
                    consider(tau, Contact::Pair(i, j));
                }
            }
        }
        for (w, wall) in c.walls.iter().enumerate() {
            for i in 0..start.len() {
                if let Some(tau) = wall_time(start, end, wall, i) {
                    consider(tau, Contact::Wall(w, i));
                }
            }
        }
        best
    }

    /// Resolves the collisions of the step of size `dt` that went from `start` to the
    /// current state.
    pub(crate) fn resolve_collisions(&mut self, start: &State, dt: f64) -> Result<()> {
        let Some(c) = self.collisions.clone() else {
            return Ok(());
        };
        let mut from = start.clone();
        let mut remaining = dt;
        let free_flight = self.forces.is_empty();
        for _ in 0..c.max_per_step {
            let Some((tau, contact)) = self.earliest_contact(&c, &from, &self.state) else {
                return Ok(());
            };
            // The straight path is exact without forces; otherwise refine the contact
            // time on the integrator's own trajectory (Illinois root finding).
            let theta = if free_flight || tau == 0.0 {
                tau
            } else {
                self.refine_contact(&from, &c.walls, contact, tau, remaining)?
            };
            let (mut at, tension) = self.advance(&from, theta * remaining)?;
            if !apply(&mut at, &c.walls, c.restitution, contact)? {
                // Not actually approaching on the true trajectory: carry on past it.
                from = at;
                remaining *= 1.0 - theta;
                let (end, tension) = self.advance(&from, remaining)?;
                self.state = end;
                self.tension = tension;
                continue;
            }
            self.collision_count += 1;
            from = at;
            self.tension = tension;
            remaining *= 1.0 - theta;
            if remaining <= 1e-15 * dt.abs() {
                self.state = from;
                return Ok(());
            }
            let (end, tension) = self.advance(&from, remaining)?;
            self.state = end;
            self.tension = tension;
        }
        invalid(format!(
            "more than max_per_step = {} collisions in one step at t = {}",
            c.max_per_step, start.t
        ))
    }

    /// Step fraction at which `contact` closes, between 0 and the straight-path estimate's
    /// closest approach, on the integrator's trajectory from `from` over `h`.
    fn refine_contact(
        &mut self,
        from: &State,
        walls: &[Wall],
        contact: Contact,
        tau: f64,
        h: f64,
    ) -> Result<f64> {
        let g0 = gap(from, walls, contact);
        // Bracket: find a fraction where the true trajectory overlaps.
        let mut hi = tau;
        let mut g_hi = gap(&self.advance(from, hi * h)?.0, walls, contact);
        while g_hi > 0.0 && hi < 1.0 {
            hi = (hi * 1.5).min(1.0);
            g_hi = gap(&self.advance(from, hi * h)?.0, walls, contact);
        }
        if g_hi > 0.0 {
            return Ok(tau); // no overlap on the true path; apply() will decide
        }
        let (mut lo, mut g_lo) = (0.0, g0);
        let mut side = 0;
        for _ in 0..100 {
            if hi - lo <= 1e-13 {
                break;
            }
            let mut theta = (lo * g_hi - hi * g_lo) / (g_hi - g_lo);
            if !(theta > lo && theta < hi) {
                theta = 0.5 * (lo + hi);
            }
            let g = gap(&self.advance(from, theta * h)?.0, walls, contact);
            if g > 0.0 {
                (lo, g_lo) = (theta, g);
                if side == -1 {
                    g_hi *= 0.5;
                }
                side = -1;
            } else {
                (hi, g_hi) = (theta, g);
                if side == 1 {
                    g_lo *= 0.5;
                }
                side = 1;
            }
        }
        Ok(hi)
    }
}
