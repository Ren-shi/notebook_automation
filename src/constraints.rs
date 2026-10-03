//! Holonomic distance constraints (rigid rods), integrated with RATTLE.
//!
//! RATTLE is velocity Verlet with Lagrange multipliers: after the drift, positions are
//! projected back onto the constraints (SHAKE), and after the final kick, velocities are
//! projected onto the constraints' tangent space. It is second order, symmetric and
//! symplectic on the constrained phase space, so `yoshida4`'s triple-jump composition of
//! it is fourth order and symplectic too (Reich 1996).

use crate::error::{invalid, Result};
use crate::state::State;
use crate::vec3::Vec3;

/// What the far end of a rod is attached to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Anchor {
    Particle(usize),
    Point(Vec3),
}

/// Keeps `|pos[i] - anchor| = length`: a massless rigid rod.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rod {
    pub i: usize,
    pub anchor: Anchor,
    pub length: f64,
}

impl Rod {
    fn j(&self) -> Option<usize> {
        match self.anchor {
            Anchor::Particle(j) => Some(j),
            Anchor::Point(_) => None,
        }
    }

    fn other_end(&self, pos: &[Vec3]) -> Vec3 {
        match self.anchor {
            Anchor::Particle(j) => pos[j],
            Anchor::Point(p) => p,
        }
    }

    /// Current `pos[i] - anchor`.
    pub fn separation(&self, pos: &[Vec3]) -> Vec3 {
        pos[self.i] - self.other_end(pos)
    }

    pub fn references_particle(&self, k: usize) -> bool {
        self.i == k || self.j() == Some(k)
    }

    pub fn name(&self) -> String {
        match self.anchor {
            Anchor::Particle(j) => format!("Rod({}, {})", self.i, j),
            Anchor::Point(p) => format!("Rod({}, [{}, {}, {}])", self.i, p.x, p.y, p.z),
        }
    }
}

/// Identifies a constraint within a [`Constraints`] set; stays valid until it is removed.
pub type ConstraintId = usize;

/// The constraints of a world and the solver settings.
#[derive(Debug, Clone)]
pub struct Constraints {
    rods: Vec<(ConstraintId, Rod)>,
    next_id: ConstraintId,
    /// Relative tolerance: rod lengths are held to `length * (1 ± tolerance)`, and the
    /// velocity along each rod to `tolerance` times the largest relative speed of rod ends.
    pub tolerance: f64,
    /// Iteration limit of each projection; reaching it fails the step.
    pub max_iterations: usize,
}

impl Default for Constraints {
    fn default() -> Self {
        Self {
            rods: Vec::new(),
            next_id: 0,
            tolerance: 1e-10,
            max_iterations: 100,
        }
    }
}

impl Constraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add(&mut self, rod: Rod) -> ConstraintId {
        let id = self.next_id;
        self.next_id += 1;
        self.rods.push((id, rod));
        id
    }

    pub fn remove(&mut self, id: ConstraintId) -> Result<Rod> {
        match self.rods.iter().position(|(i, _)| *i == id) {
            Some(p) => Ok(self.rods.remove(p).1),
            None => invalid(format!("no constraint with id {id}")),
        }
    }

    pub fn get(&self, id: ConstraintId) -> Result<&Rod> {
        match self.rods.iter().find(|(i, _)| *i == id) {
            Some((_, r)) => Ok(r),
            None => invalid(format!("no constraint with id {id}")),
        }
    }

    /// Rebuilds a set with the given ids, e.g. from a checkpoint.
    pub fn from_parts(rods: Vec<(ConstraintId, Rod)>, next_id: ConstraintId) -> Result<Self> {
        for (k, (id, _)) in rods.iter().enumerate() {
            if *id >= next_id {
                return invalid(format!("constraint id {id} is not below next_id {next_id}"));
            }
            if rods[..k].iter().any(|(other, _)| other == id) {
                return invalid(format!("duplicate constraint id {id}"));
            }
        }
        Ok(Self {
            rods,
            next_id,
            ..Self::default()
        })
    }

    pub fn next_id(&self) -> ConstraintId {
        self.next_id
    }

    pub fn iter(&self) -> impl Iterator<Item = (ConstraintId, &Rod)> {
        self.rods.iter().map(|(id, r)| (*id, r))
    }

    pub fn len(&self) -> usize {
        self.rods.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rods.is_empty()
    }

    pub fn clear(&mut self) {
        self.rods.clear();
    }

    /// Names of the constraints that use particle `k`.
    pub fn referencing(&self, k: usize) -> Vec<String> {
        self.rods
            .iter()
            .filter(|(_, r)| r.references_particle(k))
            .map(|(_, r)| r.name())
            .collect()
    }

    /// Particle `removed` was deleted: shift indices above it down by one.
    pub fn particle_removed(&mut self, removed: usize) {
        let shift = |k: &mut usize| {
            if *k > removed {
                *k -= 1;
            }
        };
        for (_, r) in &mut self.rods {
            shift(&mut r.i);
            if let Anchor::Particle(j) = &mut r.anchor {
                shift(j);
            }
        }
    }

    /// Largest `| |separation| / length - 1 |` over all rods.
    pub fn max_violation(&self, pos: &[Vec3]) -> f64 {
        self.rods
            .iter()
            .map(|(_, r)| (r.separation(pos).norm() / r.length - 1.0).abs())
            .fold(0.0, f64::max)
    }

    /// Checks indices, lengths and masses; massless particles cannot be constrained.
    pub fn validate(&self, s: &State) -> Result<()> {
        if !(self.tolerance > 0.0 && self.tolerance.is_finite() && self.max_iterations > 0) {
            return invalid(format!(
                "constraint tolerance must be positive and max_iterations at least 1, got {} and {}",
                self.tolerance, self.max_iterations
            ));
        }
        for (_, r) in &self.rods {
            check_rod(r, s)?;
        }
        Ok(())
    }

    /// The rods that can move (at least one end is a free particle), with their ends
    /// mapped to compact slots and their inverse masses.
    fn system(&self, s: &State) -> System {
        let w = |k: usize| if s.pinned[k] { 0.0 } else { 1.0 / s.mass[k] };
        let movable = |r: &Rod| w(r.i) + r.j().map_or(0.0, w) > 0.0;
        // Compact slots for the particles involved: sorted, so lookups are a binary search.
        let mut particles: Vec<usize> = self
            .rods
            .iter()
            .filter(|(_, r)| movable(r))
            .flat_map(|(_, r)| [Some(r.i), r.j()])
            .flatten()
            .collect();
        particles.sort_unstable();
        particles.dedup();
        let slot = |k: usize| particles.binary_search(&k).expect("listed above");
        let rods = self
            .rods
            .iter()
            .enumerate()
            .filter(|(_, (_, r))| movable(r)) // both ends fixed: nothing can move
            .map(|(c, (_, r))| SysRod {
                c,
                i: r.i,
                j: r.j(),
                si: slot(r.i),
                sj: r.j().map(slot),
                wi: w(r.i),
                wj: r.j().map_or(0.0, w),
            })
            .collect();
        System {
            rods,
            slots: particles.len(),
        }
    }

    /// SHAKE: moves `s.pos` back onto the constraints along the rod directions at the start
    /// of the substep (`old`), and applies the matching impulses to `s.vel` (`h` is the
    /// substep).
    ///
    /// All rods are solved together (M-SHAKE): each outer iteration solves the linearised
    /// system `G0 M^-1 G0^T lambda = residual` by conjugate gradients, with `G0` built from
    /// the rod directions at the start of the step. That chord approximation converges
    /// quickly because the directions change little in a step, while solving the coupled
    /// system at once avoids the slow convergence of rod-by-rod iteration on long chains.
    pub(crate) fn project_positions(&self, s: &mut State, old: &[Vec3], h: f64) -> Result<()> {
        let sys = self.system(s);
        if sys.rods.is_empty() {
            return Ok(());
        }
        let dirs: Vec<Vec3> = sys
            .rods
            .iter()
            .map(|r| self.rods[r.c].1.separation(old))
            .collect();
        let mut rhs = vec![0.0; sys.rods.len()];
        // After converging, one more correction squares what is left. That matters for the
        // reported tensions: a position residual `delta` becomes a velocity error `delta / h`
        // in the next step, and a tension error of order `delta / h^2`.
        let mut polishing = false;
        for _ in 0..self.max_iterations {
            let mut converged = true;
            for (k, r) in sys.rods.iter().enumerate() {
                let rod = &self.rods[r.c].1;
                let d = rod.separation(&s.pos);
                let l2 = rod.length * rod.length;
                // Linearised: 2 d . delta_d = L^2 - |d|^2, with d ~ d0.
                rhs[k] = 0.5 * (l2 - d.norm_squared());
                if rhs[k].abs() > self.tolerance * l2 {
                    converged = false;
                }
                let d0 = dirs[k];
                if d.dot(d0) <= 0.1 * d.norm() * d0.norm() {
                    return invalid(format!(
                        "{} could not be satisfied: it turned by more than 80 degrees in one step (time step too large?)",
                        rod.name()
                    ));
                }
            }
            if converged && polishing {
                return Ok(());
            }
            if converged {
                polishing = true;
            }
            let floor = 1e-3 * self.tolerance;
            let target: Vec<f64> = sys
                .rods
                .iter()
                .zip(&rhs)
                .map(|(r, b)| (1e-6 * b.abs()).max(floor * self.rods[r.c].1.length.powi(2)))
                .collect();
            let lambda = sys.solve(&dirs, &rhs, |k, res| res.abs() <= target[k])?;
            for ((r, d0), l) in sys.rods.iter().zip(&dirs).zip(&lambda) {
                s.pos[r.i] += *d0 * (l * r.wi);
                s.vel[r.i] += *d0 * (l * r.wi / h);
                if let Some(j) = r.j {
                    s.pos[j] -= *d0 * (l * r.wj);
                    s.vel[j] -= *d0 * (l * r.wj / h);
                }
            }
            if polishing {
                return Ok(());
            }
        }
        invalid(format!(
            "constraint positions did not converge in {} iterations (max violation {:.3e})",
            self.max_iterations,
            self.max_violation(&s.pos)
        ))
    }

    /// Removes the velocity components along the rods. `tension[c]` receives the tension
    /// of rod `c` implied by the impulse, taken as acting over the half kick `h / 2`
    /// (positive: pulling the ends together).
    pub(crate) fn project_velocities(
        &self,
        s: &mut State,
        h: f64,
        tension: &mut Vec<f64>,
    ) -> Result<()> {
        tension.clear();
        tension.resize(self.rods.len(), 0.0);
        let sys = self.system(s);
        if sys.rods.is_empty() {
            return Ok(());
        }
        let dirs: Vec<Vec3> = sys
            .rods
            .iter()
            .map(|r| self.rods[r.c].1.separation(&s.pos))
            .collect();
        let rel_vel = |s: &State, r: &SysRod| s.vel[r.i] - r.j.map_or(Vec3::ZERO, |j| s.vel[j]);
        // Converged when every rod's radial velocity is below `tolerance` times the largest
        // relative speed of rod ends at entry. (Relative to the current speed instead, a
        // system the projection brings to rest could never converge.)
        let scale = sys
            .rods
            .iter()
            .map(|r| rel_vel(s, r).norm())
            .fold(0.0, f64::max);
        let limit: Vec<f64> = dirs
            .iter()
            .map(|d| self.tolerance * d.norm() * scale)
            .collect();
        let mut radial = vec![0.0; sys.rods.len()];
        // The projection is linear, so one solve does it; repeat only if rounding left a
        // residual above tolerance.
        for _ in 0..self.max_iterations {
            for ((b, r), d) in radial.iter_mut().zip(&sys.rods).zip(&dirs) {
                *b = d.dot(rel_vel(s, r));
            }
            if radial.iter().zip(&limit).all(|(b, l)| b.abs() <= *l) {
                return Ok(());
            }
            let mu = sys.solve(&dirs, &radial, |k, res| res.abs() <= 0.1 * limit[k])?;
            for ((r, d), m) in sys.rods.iter().zip(&dirs).zip(&mu) {
                s.vel[r.i] -= *d * (m * r.wi);
                if let Some(j) = r.j {
                    s.vel[j] += *d * (m * r.wj);
                }
                // Impulse on i is -mu d over h/2: a force of magnitude 2 mu |d| / h toward the anchor.
                tension[r.c] += 2.0 * m * d.norm() / h;
            }
        }
        invalid(format!(
            "constraint velocities did not converge in {} iterations",
            self.max_iterations
        ))
    }
}

/// A movable rod in a [`System`].
struct SysRod {
    /// Index in [`Constraints::rods`].
    c: usize,
    i: usize,
    j: Option<usize>,
    /// Compact slots of the ends, for the matrix-vector product.
    si: usize,
    sj: Option<usize>,
    wi: f64,
    wj: f64,
}

/// The coupled constraint equations: `A = G M^-1 G^T` with `G` built from one direction
/// per rod, so `A[c][e]` is nonzero only for rods sharing a particle.
struct System {
    rods: Vec<SysRod>,
    slots: usize,
}

impl System {
    /// `y = A x`, without forming `A`: spread `x` onto the particles, then gather.
    fn apply(&self, dirs: &[Vec3], x: &[f64], u: &mut [Vec3], y: &mut [f64]) {
        u.fill(Vec3::ZERO);
        for ((r, d), xk) in self.rods.iter().zip(dirs).zip(x) {
            u[r.si] += *d * (xk * r.wi);
            if let Some(sj) = r.sj {
                u[sj] -= *d * (xk * r.wj);
            }
        }
        for ((r, d), yk) in self.rods.iter().zip(dirs).zip(y.iter_mut()) {
            *yk = d.dot(u[r.si] - r.sj.map_or(Vec3::ZERO, |sj| u[sj]));
        }
    }

    /// Solves `A x = b` by Jacobi-preconditioned conjugate gradients, until
    /// `done(k, residual_k)` holds for every rod `k`.
    fn solve(
        &self,
        dirs: &[Vec3],
        b: &[f64],
        done: impl Fn(usize, f64) -> bool,
    ) -> Result<Vec<f64>> {
        let n = self.rods.len();
        let finished = |r: &[f64]| r.iter().enumerate().all(|(k, &x)| done(k, x));
        let mut x = vec![0.0; n];
        let mut r = b.to_vec();
        if finished(&r) {
            return Ok(x);
        }
        let diag: Vec<f64> = self
            .rods
            .iter()
            .zip(dirs)
            .map(|(rod, d)| (rod.wi + rod.wj) * d.norm_squared())
            .collect();
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
        let mut z: Vec<f64> = r.iter().zip(&diag).map(|(r, d)| r / d).collect();
        let mut p = z.clone();
        let mut rz = dot(&r, &z);
        let mut u = vec![Vec3::ZERO; self.slots];
        let mut ap = vec![0.0; n];
        // In exact arithmetic CG finishes in n iterations; allow for rounding.
        for _ in 0..(4 * n + 50) {
            self.apply(dirs, &p, &mut u, &mut ap);
            let pap = dot(&p, &ap);
            if pap.is_nan() || pap <= 0.0 {
                break; // singular direction: redundant or degenerate rods
            }
            let alpha = rz / pap;
            for k in 0..n {
                x[k] += alpha * p[k];
                r[k] -= alpha * ap[k];
            }
            if finished(&r) {
                return Ok(x);
            }
            for k in 0..n {
                z[k] = r[k] / diag[k];
            }
            let rz_new = dot(&r, &z);
            let beta = rz_new / rz;
            rz = rz_new;
            for k in 0..n {
                p[k] = z[k] + beta * p[k];
            }
        }
        invalid(format!(
            "the constraint equations could not be solved ({n} rods); are some rods redundant or degenerate?"
        ))
    }
}

pub(crate) fn check_rod(r: &Rod, s: &State) -> Result<()> {
    let n = s.len();
    let ends = [Some(r.i), r.j()];
    for k in ends.into_iter().flatten() {
        if k >= n {
            return invalid(format!(
                "{}: particle index {k} out of range (have {n} particles)",
                r.name()
            ));
        }
        if s.mass[k] == 0.0 {
            return invalid(format!(
                "{}: particle {k} is massless and cannot be constrained",
                r.name()
            ));
        }
    }
    if r.j() == Some(r.i) {
        return invalid(format!("{}: a rod needs two different particles", r.name()));
    }
    if !(r.length.is_finite() && r.length > 0.0) {
        return invalid(format!(
            "{}: length must be positive, got {}",
            r.name(),
            r.length
        ));
    }
    if let Anchor::Point(p) = r.anchor {
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return invalid(format!("{}: anchor must be finite", r.name()));
        }
    }
    Ok(())
}
