use crate::vec3::Vec3;

/// Phase-space state of a system of point particles.
#[derive(Debug, Clone, Default)]
pub struct State {
    pub t: f64,
    pub pos: Vec<Vec3>,
    pub vel: Vec<Vec3>,
    pub mass: Vec<f64>,
    /// Electric charge of each particle (zero unless set).
    pub charge: Vec<f64>,
    /// Pinned particles never move but still exert forces.
    pub pinned: Vec<bool>,
}

impl State {
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a particle and returns its index.
    pub fn add_particle(&mut self, pos: Vec3, vel: Vec3, mass: f64) -> usize {
        self.pos.push(pos);
        self.vel.push(vel);
        self.mass.push(mass);
        self.charge.push(0.0);
        self.pinned.push(false);
        self.pos.len() - 1
    }

    /// Removes particle `i`; particles above it shift down by one index.
    pub fn remove_particle(&mut self, i: usize) {
        self.pos.remove(i);
        self.vel.remove(i);
        self.mass.remove(i);
        self.charge.remove(i);
        self.pinned.remove(i);
    }

    pub fn len(&self) -> usize {
        self.pos.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.vel
            .iter()
            .zip(&self.mass)
            .map(|(v, m)| 0.5 * m * v.norm_squared())
            .sum()
    }

    pub fn total_charge(&self) -> f64 {
        self.charge.iter().sum()
    }

    pub fn total_mass(&self) -> f64 {
        self.mass.iter().sum()
    }

    pub fn momentum(&self) -> Vec3 {
        self.vel
            .iter()
            .zip(&self.mass)
            .fold(Vec3::ZERO, |p, (v, m)| p + *v * *m)
    }

    /// Total angular momentum about the origin.
    pub fn angular_momentum(&self) -> Vec3 {
        self.pos
            .iter()
            .zip(&self.vel)
            .zip(&self.mass)
            .fold(Vec3::ZERO, |l, ((r, v), m)| l + r.cross(*v * *m))
    }

    pub fn center_of_mass(&self) -> Vec3 {
        let m = self.total_mass();
        if m == 0.0 {
            return Vec3::ZERO;
        }
        self.pos
            .iter()
            .zip(&self.mass)
            .fold(Vec3::ZERO, |c, (r, mi)| c + *r * *mi)
            / m
    }
}
